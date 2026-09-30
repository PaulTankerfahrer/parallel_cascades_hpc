/* ============================================================================
 *  cascade_cuda.cu  --  CUDA-GPU-Version der Kollisionskaskaden-Simulation
 * ----------------------------------------------------------------------------
 *
 *  ÜBERBLICK: WARUM GPU?
 *  ----------------------
 *  Die Kraftberechnung (compute_forces in der seriellen Version) ist
 *  embarassingly parallel: Atom i braucht nur die Daten seiner Nachbarn,
 *  nicht aller anderen Atome.  Eine GPU mit tausenden kleinen Threads
 *  kann alle N Atome gleichzeitig berechnen.
 *
 *  GRUNDPRINZIP: EIN THREAD PRO ATOM
 *  -----------------------------------
 *  Thread-Index:  i = blockIdx.x * blockDim.x + threadIdx.x
 *  Jeder Thread berechnet alle Kraefte auf Atom i.
 *
 *  Beispiel fuer N=2e6, blockDim=256:
 *    Benoetigte Bloecke: B = ceil(2e6 / 256) = 7813 Bloecke
 *    Jeder Block: 256 Threads -> 7813*256 = 2.000.128 Threads (etwas mehr als N, ok)
 *
 *  GATHER-FORMULIERUNG (kein atomicAdd im Force-Kernel):
 *  -------------------------------------------------------
 *  In der seriellen Version nutzt compute_forces() Newton's 3. Gesetz:
 *    fx[a] += f*dx;  fx[b] -= f*dx;  // beide Atome gleichzeitig
 *  Das geht NICHT auf der GPU, weil zwei Threads gleichzeitig auf fx[b]
 *  schreiben wuerden (Race Condition).
 *
 *  Loesung -- Gather: Jeder Thread i sammelt NUR die Kraefte auf sich selbst.
 *  Atom j schreibt NICHT in fx[i]; stattdessen liest Thread i selbst von px[j].
 *  Das ist sicher, weil jeder Thread seinen eigenen fx[i]-Wert besitzt.
 *
 *  Konsequenz: Jede Wechselwirkung (i,j) wird ZWEIMAL berechnet
 *  (Thread i rechnet die Kraft von j auf i; Thread j rechnet die Kraft von i auf j).
 *  Das ist doppelt so teuer wie der serielle Halb-Stern, aber ermoeglicht die
 *  volle GPU-Parallelisierung ohne atomicAdd.
 *
 *  CELL LIST AUF DER GPU:
 *  ------------------------
 *  Das serielle Aufbauen der Cell List (Schleife ueber alle Atome) wird auf
 *  der GPU durch zwei Kernel ersetzt:
 *  k_clear_cells: setzt alle Zellen parallel auf -1 (ein Thread pro Zelle)
 *  k_build_cells: jedes Atom haengt sich atomar in seine Zelle ein (atomicExch)
 *
 *  ENERGIE-REDUKTION AUF DER GPU:
 *  ---------------------------------
 *  k_energy() berechnet pro Atom einen Energie-Beitrag -> e_atom[i].
 *  k_reduce() summiert e_atom[] auf der GPU (shared-memory Baumreduktion).
 *  Nur EIN double (die Gesamtsumme) wird zum Host uebertragen.
 *
 *  WICHTIGE MESSWERTE:
 *  --------------------
 *  T_CPU_1CORE=1152,9s, A100 fp64: ~1,58s -> Speedup 731x
 *  A100-Auslastung: 45% Compute, 37% Speicher (latenzlimitiert)
 *  Ursache: niedrige Okkupanz (50%) durch 62 Register/Thread bei fp64
 *
 *  Build:  nvcc -O3 -arch=sm_80 -o cascade_cuda cascade_cuda.cu
 *  Run:    ./cascade_cuda
 *  (Parameter sind fest eingebaut, identisch zu params.ini Scaling-Config)
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>

/* Max. Nachbarn im Dreiecksgitter: 6 pro Innenatom */
#define MAXNBR 6

/* ============================================================================
 *  CUDA_CHECK  --  Fehlerbehandlungs-Makro
 * ---------------------------------------------------------------------------
 *  Jeder CUDA-Aufruf kann fehlschlagen.  CUDA_CHECK wrappt den Aufruf:
 *  Bei Fehler: Fehlermeldung mit Datei, Zeile und Fehlertext ausgeben, dann exit(1).
 *
 *  Ohne dieses Makro wuerde man Fehler still ignorieren (der Code laeuft dann
 *  mit falschen oder undefinierten Daten weiter -- sehr schwer zu debuggen).
 *
 *  Beispiel:  CUDA_CHECK(cudaMalloc(&ptr, size));
 *  Wenn cudaMalloc scheitert (z.B. Out-of-Memory), wird der Fehler sofort gemeldet.
 * ========================================================================== */
#define CUDA_CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){       \
    fprintf(stderr,"CUDA %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e));\
    exit(1);} } while(0)

/* ============================================================================
 *  REFERENZWERTE fuer den CPU-GPU-Speedup-Vergleich
 * ---------------------------------------------------------------------------
 *  Diese Werte wurden ECHT GEMESSEN und sind in dieser Datei hartkodiert,
 *  damit das Programm automatisch die Speedup-Zahlen berechnet.
 *
 *  T_CPU_1CORE:  Serielle Wandzeit fuer N=1414^2, 3000 Schritte auf einem
 *                EPYC-7763-Kern (AMD EPYC Milan, Noctua-2).
 *  SPEEDUP_64:   Gemessener Speedup der MPI-Version mit 64 Kernen = voller Knoten.
 *  PEAK_BW_GBS:  Theoretische HBM2-Speicherbandbreite der NVIDIA A100 40GB.
 *  BYTES_PER_ATOM: Modellierter Speicherbedarf pro Atom pro Force-Kernel-Aufruf
 *                  (px,py,fx,fy,nbr,intact,nbr0,head,next -- "heisser Stream").
 *                  VORSICHT: das ist ein Modell, nicht die gemessene ncu-Zahl.
 * ========================================================================== */
#define T_CPU_1CORE    1152.9    /* s,  1 EPYC-Kern                            */
#define SPEEDUP_64     51.2      /* gemessener 64-Kern-Speedup                 */
#define PEAK_BW_GBS    1555.0    /* GB/s, A100 40GB HBM2 Peak                 */
#define BYTES_PER_ATOM 76.0      /* B/Atom (MODELL, nicht ncu-autoritativ)    */

/* ============================================================================
 *  DEVICE-FUNKTIONEN  (laufen NUR auf der GPU, aufrufbar von anderen Kerneln)
 * ---------------------------------------------------------------------------
 *  __device__: Diese Funktionen werden zu GPU-Maschinencode kompiliert.
 *  Sie koennen NICHT vom Host (CPU) aufgerufen werden.
 *  inline: Der Compiler soll die Funktion direkt einbetten (kein Funktionsaufruf-
 *  Overhead auf der GPU, wo Funktionsaufruf teurer als auf der CPU ist).
 *
 *  Identische Formeln wie in der seriellen Version (cascade_serial.c).
 *  rep_force(r): Abstoss-Kraft bei Abstand r  (r < RCUT vorausgesetzt)
 *  rep_pot(r):   Abstoss-Potenzial (fuer Energiebuchhaltung in k_energy)
 * ========================================================================== */
__device__ inline double rep_force(double r,double RCUT,double KREP,double REPN){
    return KREP*(pow(RCUT/r,REPN) - 1.0);
}
__device__ inline double rep_pot(double r,double RCUT,double KREP,double REPN){
    double n=REPN;
    return KREP*((RCUT - pow(RCUT,n)*pow(r,1.0-n))/(1.0-n) - (RCUT - r));
}

/* ============================================================================
 *  k_clear_cells  --  Cell-List-Zellen leeren (Kernel 1 von 2)
 * ---------------------------------------------------------------------------
 *  Setzt alle Eintraege des head[]-Arrays auf -1 (leer).
 *  Jeder Thread loescht eine Zelle.
 *
 *  Parameter:
 *    head[ncell]: Array der Listenkopf-Zeiger, einer pro Zelle
 *    ncell:       Gesamtanzahl Zellen
 *
 *  Grid-Konfiguration: Bc = ceil(ncell/T) Bloecke, T Threads je Block.
 *  Threads die c >= ncell haben: einfach nichts tun (Guard if(c<ncell)).
 * ========================================================================== */
__global__ void k_clear_cells(int* head, int ncell){
    int c = blockIdx.x*blockDim.x + threadIdx.x;
    if (c < ncell) head[c] = -1;
}

/* ============================================================================
 *  k_build_cells  --  Atome in Cell List eintragen (Kernel 2 von 2)
 * ---------------------------------------------------------------------------
 *  Jeder Thread i traegt Atom i in seine Zelle ein.
 *
 *  ATOMICEXCH -- der entscheidende Trick:
 *  Die serielle Variante ist:
 *    cell_next[i] = cell_head[c];  // i zeigt auf alten Kopf
 *    cell_head[c] = i;             // i wird neuer Kopf
 *  Das DARF NICHT von zwei Threads gleichzeitig ausgefuehrt werden.
 *
 *  Loesung: atomicExch(&head[c], i) ist eine ATOMARE Operation:
 *  Sie setzt head[c]=i und gibt den alten Wert von head[c] zurueck.
 *  In einem einzigen ununterbrechbaren Schritt -- keine Race Condition.
 *  Das Ergebnis wird in next[i] gespeichert (i zeigt auf den alten Kopf).
 *
 *  Resultat: Eine verkettete Liste wird korrekt aufgebaut, auch wenn
 *  tausende Threads gleichzeitig in dieselbe Zelle schreiben.
 * ========================================================================== */
__global__ void k_build_cells(const double* px, const double* py,
                              int* head, int* next,
                              double x0, double y0, double cs,
                              int ncx, int ncy, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    /* Zellenindex aus Position berechnen (identisch zur seriellen cell_of()) */
    int cx=(int)((px[i]-x0)/cs), cy=(int)((py[i]-y0)/cs);
    if(cx<0)cx=0; if(cx>=ncx)cx=ncx-1;
    if(cy<0)cy=0; if(cy>=ncy)cy=ncy-1;
    int c = cy*ncx + cx;
    /* Atomar: i wird neuer Listenkopf, der alte Kopf landet in next[i] */
    next[i] = atomicExch(&head[c], i);
}

/* ============================================================================
 *  k_force  --  Kraftberechnung (Haupt-Kernel, dominiert ~90% der Laufzeit)
 * ---------------------------------------------------------------------------
 *  Jeder Thread i berechnet ALLE Kraefte auf Atom i:
 *   (a) Federkraefte von gebundenen Nachbarn (Gather: Atom i liest px[j], py[j])
 *   (b) Abstoss-Kraefte von allen Atomen j in den 9 Nachbarzellen mit r < RCUT
 *
 *  Unterschied zur seriellen Version:
 *   SERIELL: compute_forces benutzt Newton 3: fx[a]+=f; fx[b]-=f (ein Paar, zwei Updates)
 *   GPU:     Thread i schreibt NUR in fx[i], liest aber px[j] fuer alle j (Gather)
 *            -> Jede Wechselwirkung wird ZWEIMAL berechnet (i sieht j, j sieht i)
 *            -> Kein atomicAdd noetig auf fx[], kein Race Condition
 *
 *  BRUCH auf der GPU:
 *  intact[base+e]=0 ist sicher, weil nur Thread i in intact[i*MAXNBR+e] schreibt.
 *  atomicAdd(ebroken, brk): HIER brauchen wir doch atomicAdd, weil viele Threads
 *  gleichzeitig Bruche zaehlen koennen (ebroken ist ein einzelner globaler Wert).
 *  Das passiert selten (nur beim Bruch) -> minimaler Overhead.
 *
 *  Register-Verbrauch: ~62 Register/Thread bei fp64 (viele lokale Variablen:
 *  xi,yi,fxi,fyi,brk,base,cnt,j,dx,dy,r,f,rc2,cx,cy,...).
 *  Konsequenz: max. Okkupanz = 50% bei A100 (weniger parallele Warps).
 * ========================================================================== */
__global__ void k_force(const double* px, const double* py,
                        double* fx, double* fy,
                        const int* nbr, int* intact, const int* nbr0,
                        const int* head, const int* next,
                        double x0, double y0, double cs, int ncx, int ncy,
                        double K, double L0, double maxstr,
                        double RCUT, double KREP, double REPN,
                        double* ebroken, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    double xi=px[i], yi=py[i];
    double fxi=0.0, fyi=0.0, brk=0.0;
    int base=i*MAXNBR, cnt=nbr0[i];   /* cnt = Anfangszahl Nachbarn (kann sich aendern!) */

    /* (a) Federkraefte -- Gather: Thread i liest px[j], py[j] und schreibt nur fxi,fyi */
    for (int e=0;e<cnt;e++){
        if(!intact[base+e]) continue;          /* gerissene Feder ueberspringen */
        int j=nbr[base+e];                    /* lokaler Index des j-ten Nachbarn */
        double dx=px[j]-xi, dy=py[j]-yi, r=sqrt(dx*dx+dy*dy);
        if (r > L0*maxstr){
            /* Bruch: intact auf 0 setzen, Energie akkumulieren */
            intact[base+e]=0;
            brk += 0.5*(0.5*K*(r-L0)*(r-L0));  /* Halbe Energie: Gather zaehlt jeden Bruch 2x */
            continue;
        }
        double f=K*(r-L0)/r;
        fxi += f*dx;  fyi += f*dy;             /* Anziehung Richtung j */
    }

    /* (b) Abstoss-Kraefte via Cell List (voller 3x3-Stern = 9 Nachbarzellen) */
    double rc2=RCUT*RCUT;
    int cx=(int)((xi-x0)/cs), cy=(int)((yi-y0)/cs);
    if(cx<0)cx=0; if(cx>=ncx)cx=ncx-1;
    if(cy<0)cy=0; if(cy>=ncy)cy=ncy-1;
    for (int oy=-1; oy<=1; oy++) for (int ox=-1; ox<=1; ox++){
        int nx=cx+ox, ny=cy+oy;
        if(nx<0||nx>=ncx||ny<0||ny>=ncy) continue;
        int c=ny*ncx+nx;
        /* Alle Atome j in dieser Zelle durchiterieren */
        for (int j=head[c]; j!=-1; j=next[j]){
            if (j==i) continue;                /* kein Atom mit sich selbst */
            double dx=xi-px[j], dy=yi-py[j];  /* Richtungsvektor: von j zu i (abstossend) */
            double r2=dx*dx+dy*dy;
            if (r2>=rc2 || r2==0.0) continue; /* zu weit weg oder Selbstvergleich */
            /* Gebundene Paare ausschliessen (Feder wirkt statt Abstossung) */
            int bonded=0;
            for(int e=0;e<cnt;e++){
                if(intact[base+e] && nbr[base+e]==j){ bonded=1; break; }
            }
            if (bonded) continue;
            double r=sqrt(r2), f=rep_force(r,RCUT,KREP,REPN)/r;
            fxi += f*dx;  fyi += f*dy;         /* Gather: nur Kraft auf i */
        }
    }
    /* Ergebnisse in globale Arrays schreiben (jeder Thread schreibt eigene Indizes) */
    fx[i]=fxi;  fy[i]=fyi;
    /* Bruch-Energie atomar akkumulieren (selten, aber sicher) */
    if (brk!=0.0) atomicAdd(ebroken, brk);
}

/* ============================================================================
 *  k_kick_drift  --  Velocity-Verlet Halber-Kick + Drift (Schritt 1+2)
 * ---------------------------------------------------------------------------
 *  Jeder Thread i:
 *    v(t+dt/2) = v(t) + (F/m) * dt/2    <- halber Kick
 *    x(t+dt)   = x(t) + v(t+dt/2) * dt  <- Drift
 *  Kein Synchronisierungsproblem: Thread i schreibt nur in vx[i],vy[i],px[i],py[i].
 * ========================================================================== */
__global__ void k_kick_drift(double* px,double* py,double* vx,double* vy,
                             const double* fx,const double* fy,const double* m,
                             double dt,int N){
    int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=N)return;
    vx[i]+=0.5*dt*fx[i]/m[i]; vy[i]+=0.5*dt*fy[i]/m[i];
    px[i]+=dt*vx[i];          py[i]+=dt*vy[i];
}

/* ============================================================================
 *  k_kick  --  Velocity-Verlet zweiter halber Kick (Schritt 4)
 * ---------------------------------------------------------------------------
 *  Jeder Thread i:
 *    v(t+dt) = v(t+dt/2) + (F_neu/m) * dt/2
 *  Wird mit den Kraeften aus dem NEUEN Zeitschritt aufgerufen.
 * ========================================================================== */
__global__ void k_kick(double* vx,double* vy,const double* fx,const double* fy,
                       const double* m,double dt,int N){
    int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=N)return;
    vx[i]+=0.5*dt*fx[i]/m[i]; vy[i]+=0.5*dt*fy[i]/m[i];
}

/* ============================================================================
 *  k_energy  --  Energie-Beitrag pro Atom (Gather, Halbes-Zaehlen)
 * ---------------------------------------------------------------------------
 *  Jeder Thread i berechnet seinen Energie-Beitrag e_atom[i]:
 *    - E_kin:    0.5 * m * v^2                   (komplett, Thread i besitzt Atom i)
 *    - E_spring: 0.5 * (0.5 * K * (r-L0)^2)     (halbe Bindungsenergie)
 *    - E_rep:    0.5 * V_rep(r)                  (halbe Abstossungsenergie)
 *
 *  WARUM 0.5 fuer Potential-Energien?
 *  Die Gather-Formulierung berechnet jeden Beitrag ZWEIMAL:
 *  Thread i sieht Paar (i,j), Thread j sieht Paar (j,i).
 *  Beide zaehlen dieselbe Wechselwirkungsenergie.
 *  Mit Faktor 0.5 stimmt die Summe ueber alle Atome:
 *    Summe_i( 0.5 * E_ij ) = E_ij  (einmal gezaehlt)
 *
 *  Ergebnis: e_atom[] enthaelt pro Atom seinen Anteil.
 *  Die Gesamtenergie = sum(e_atom[]) wird von k_reduce berechnet.
 * ========================================================================== */
__global__ void k_energy(const double* px,const double* py,
                         const double* vx,const double* vy,const double* m,
                         const int* nbr,const int* intact,const int* nbr0,
                         const int* head,const int* next,
                         double x0,double y0,double cs,int ncx,int ncy,
                         double K,double L0,double RCUT,double KREP,double REPN,
                         double* e_atom,int N){
    int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=N)return;
    double xi=px[i], yi=py[i];
    double e = 0.5*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]);   /* kinetische Energie (ganz) */
    int base=i*MAXNBR, cnt=nbr0[i];
    /* Halbe Federpotenzial-Energie jeder intakten Bindung */
    for (int k=0;k<cnt;k++){
        if(!intact[base+k]) continue;
        int j=nbr[base+k];
        double dx=px[j]-xi,dy=py[j]-yi,r=sqrt(dx*dx+dy*dy);
        e += 0.5*(0.5*K*(r-L0)*(r-L0));   /* 0.5: doppeltes Zaehlen */
    }
    /* Halbe Abstoss-Potenzial-Energie */
    double rc2=RCUT*RCUT;
    int cx=(int)((xi-x0)/cs), cy=(int)((yi-y0)/cs);
    if(cx<0)cx=0; if(cx>=ncx)cx=ncx-1; if(cy<0)cy=0; if(cy>=ncy)cy=ncy-1;
    for (int oy=-1;oy<=1;oy++) for (int ox=-1;ox<=1;ox++){
        int nx=cx+ox,ny=cy+oy; if(nx<0||nx>=ncx||ny<0||ny>=ncy)continue;
        int c=ny*ncx+nx;
        for (int j=head[c]; j!=-1; j=next[j]){
            if(j==i) continue;
            double dx=xi-px[j],dy=yi-py[j],r2=dx*dx+dy*dy;
            if(r2>=rc2||r2==0.0) continue;
            int bonded=0;
            for(int k=0;k<cnt;k++){ if(intact[base+k]&&nbr[base+k]==j){bonded=1;break;} }
            if(bonded) continue;
            e += 0.5*rep_pot(sqrt(r2),RCUT,KREP,REPN);  /* 0.5: doppeltes Zaehlen */
        }
    }
    e_atom[i]=e;
}

/* ============================================================================
 *  k_reduce  --  Parallele Summenreduktion auf der GPU
 * ---------------------------------------------------------------------------
 *  Summiert in[0..N) zu einem einzigen double in out[0].
 *
 *  WARUM NICHT EINFACH EINE SCHLEIFE?
 *  Eine sequenzielle Summe auf der GPU (ein Thread) wuerde N Schritte
 *  brauchen und die anderen tausend Threads unbenutzt lassen.
 *  Parallele Reduktion: log2(N) Stufen mit N/2 aktiven Threads pro Stufe.
 *
 *  SCHRITT 1 -- Grid-Stride-Load:
 *  Jeder Thread akkumuliert mehrere Elemente (stride = gridDim*blockDim).
 *  Vorteil: Beliebiges N mit fester Blockzahl; coalesced Memory-Zugriff
 *  (Thread 0 liest in[0], Thread 1 liest in[1], ..., in einem Speicherzug).
 *
 *  SCHRITT 2 -- Shared-Memory-Baumreduktion:
 *  Das Block-Ergebnis (in sdata[]) wird in log2(blockDim) Stufen summiert.
 *  Jede Stufe halbiert die Stride-Weite ("sequential addressing"):
 *    Stufe 1: sdata[0]+=sdata[blockDim/2], sdata[1]+=sdata[blockDim/2+1], ...
 *    Stufe 2: sdata[0]+=sdata[blockDim/4], ...
 *    ...bis Stride=0: sdata[0] enthaelt die Block-Summe.
 *  Vorteil: Keine Bank-Konflikte im Shared Memory (beide Operanden liegen
 *  in verschiedenen Speicher-Baenken), keine divergenten Warps.
 *
 *  SCHRITT 3 -- Block-Ergebnis zu globalem Ergebnis:
 *  atomicAdd(out, sdata[0]) addiert die Block-Summe auf das globale out[0].
 *  Hier ist atomicAdd noetig (viele Bloecke schreiben gleichzeitig).
 *  Da es nur einen atomicAdd pro Block gibt, ist der Overhead klein.
 *
 *  ANFORDERUNG: blockDim MUSS eine Zweierpotenz sein (hier 256).
 *  extern __shared__: Groesse wird beim Kernel-Aufruf als 3. Argument angegeben:
 *  k_reduce<<<gridDim, blockDim, blockDim*sizeof(double)>>>(...)
 * ========================================================================== */
__global__ void k_reduce(const double* __restrict__ in, double* out, int N){
    extern __shared__ double sdata[];
    int tid    = threadIdx.x;
    int gid    = blockIdx.x*blockDim.x + tid;
    int stride = blockDim.x*gridDim.x;

    /* Schritt 1: Grid-Stride-Load -- jeder Thread summiert mehrere Elemente */
    double local = 0.0;
    for (int k=gid; k<N; k+=stride) local += in[k];
    sdata[tid] = local;
    __syncthreads();   /* alle Threads muessen fertig sein vor dem Baum */

    /* Schritt 2: Baumreduktion in Shared Memory (sequential addressing) */
    for (int s=blockDim.x/2; s>0; s>>=1){   /* s = Stride, halbiert sich jede Stufe */
        if (tid < s) sdata[tid] += sdata[tid+s];
        __syncthreads();   /* Synchronisation zwischen Stufen */
    }
    /* Schritt 3: Block-Summe (in sdata[0]) atomar zum globalen Ergebnis addieren */
    if (tid==0) atomicAdd(out, sdata[0]);
}

/* ============================================================================
 *  gpu_reduce_sum  --  Wrapper fuer k_reduce
 * ---------------------------------------------------------------------------
 *  Setzt das Ergebnis-Skalar auf 0, startet k_reduce, liest das Ergebnis
 *  zurueck und gibt es als double zurueck.
 *
 *  Parameter:
 *    d_in:     Device-Array der Eingabewerte (e_atom[])
 *    d_scalar: Device-Skalar fuer das Ergebnis (1 double, wiederverwendbar)
 *    N:        Anzahl Elemente
 *    T:        Threads pro Block (muss Zweierpotenz sein)
 *
 *  RB = min(ceil(N/T), 1024): Blockanzahl.  Maximum 1024, damit die Grid-
 *  Stride-Load im Kernel genuegend Arbeit pro Thread hat (nicht zu viele
 *  leere Threads).
 * ========================================================================== */
static double gpu_reduce_sum(const double* d_in, double* d_scalar, int N, int T){
    CUDA_CHECK(cudaMemset(d_scalar, 0, sizeof(double)));  /* Skalar auf 0 */
    int RB = (N + T - 1)/T; if (RB > 1024) RB = 1024;   /* max. 1024 Bloecke */
    k_reduce<<<RB, T, T*sizeof(double)>>>(d_in, d_scalar, N);
    CUDA_CHECK(cudaGetLastError());   /* Kernel-Startfehler abfangen */
    double s;
    CUDA_CHECK(cudaMemcpy(&s, d_scalar, sizeof(double), cudaMemcpyDeviceToHost));
    return s;
}

/* ============================================================================
 *  add_nb  --  Hilfsfunktion: Bindung in Nachbarliste eintragen (Host)
 * ---------------------------------------------------------------------------
 *  Traegt a in die Nachbarliste von b ein und umgekehrt.
 *  c[a] und c[b] sind die aktuellen Zaehler (werden inkrementiert).
 *  Identisch zur seriellen Version.
 * ========================================================================== */
static void add_nb(int* nbr,int* c,int a,int b){ nbr[a*MAXNBR+c[a]++]=b; nbr[b*MAXNBR+c[b]++]=a; }

/* ============================================================================
 *  MAIN  --  Hauptprogramm
 * ---------------------------------------------------------------------------
 *  Ablauf:
 *  1.  Gitter aufbauen (Host): Positionen, Bindungen, PKA
 *  2.  Cell-List-Domain berechnen (Host)
 *  3.  Device-Speicher allokieren und Daten uebertragen (H2D)
 *  4.  Self-Check: GPU-Reduktion == Host-Summe (Korrektheitsbeweis)
 *  5.  Aufwaermphase (ungetaktet, damit GPU-Clocks eingeschwungen sind)
 *  6.  Getaktete Zeitschleife (NSTEPS Schritte)
 *  7.  Isolierter Force-Kernel-Microbench (fuer Bandbreite)
 *  8.  D2H-Transfer messen
 *  9.  Kennzahlen berechnen und ausgeben
 * ========================================================================== */
int main(void){
    /* ---- Konfiguration: identisch zu params.ini Scaling-Config ---------- */
    const int    NX=1414, NY=1414, N=NX*NY;     /* ~2.0 Mio. Atome            */
    const double L0=1.0, K=100.0, MAXSTR=1.15;
    const double RCUT=0.9, KREP=400.0, REPN=12.0;
    const double dt=2e-4, PKA_E=5000.0;
    const int    NSTEPS=3000, LOGEVERY=500;
    const int    WARMUP=50;     /* Aufwaermschritte (Clocks einschwingen) */
    /* FORCE_REP: Force-Kernel FORCE_REP-mal wiederholen fuer stabiles Zeitmittel */
    const int    FORCE_REP=200;
    const double DY=L0*sqrt(3.0)/2.0;
    size_t db=(size_t)N*sizeof(double);          /* Bytes fuer ein N-double-Array */
    size_t ib=(size_t)N*MAXNBR*sizeof(int);      /* Bytes fuer Nachbarliste */

    /* ---- Gitter aufbauen (Host) ----------------------------------------- */
    double *h_px=(double*)malloc(db),*h_py=(double*)malloc(db);
    double *h_vx=(double*)calloc(N,8),*h_vy=(double*)calloc(N,8),*h_m=(double*)malloc(db);
    double *h_ea=(double*)malloc(db);      /* pro-Atom Energie (fuer Self-Check) */
    int *h_nbr=(int*)malloc(ib),*h_int=(int*)malloc(ib);  /* Bindungsliste + intact[] */
    int *h_cnt=(int*)calloc(N,sizeof(int));                /* temporaerer Zaehler */
    int *h_nbr0=(int*)malloc(N*sizeof(int));               /* Anfangs-Nachbarzahlen */

    /* Positionen: Dreiecksgitter */
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){
        int id=j*NX+i;
        h_px[id]=i*L0+(j&1)*(L0*0.5);   /* gerade: i*L0, ungerade: i*L0+L0/2 */
        h_py[id]=j*DY;
        h_m[id]=1.0;
    }
    /* Bindungen: nur nach rechts und oben (jede Bindung einmal) */
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){
        int a=j*NX+i;
        if(i+1<NX) add_nb(h_nbr,h_cnt,a,j*NX+i+1);
        if(j+1<NY){
            int li,ri;
            if((j&1)==0){li=i-1;ri=i;}else{li=i;ri=i+1;}
            if(li>=0&&li<NX) add_nb(h_nbr,h_cnt,a,(j+1)*NX+li);
            if(ri>=0&&ri<NX) add_nb(h_nbr,h_cnt,a,(j+1)*NX+ri);
        }
    }
    /* nbr0: Anfangs-Nachbarzahl (wird nicht veraendert, dient als Schleifengrenze) */
    for(int i=0;i<N;i++){
        h_nbr0[i]=h_cnt[i];
        for(int e=0;e<MAXNBR;e++) h_int[i*MAXNBR+e]=(e<h_cnt[i])?1:0;
    }
    /* PKA in der Gittermitte platzieren */
    int pka=(NY/2)*NX+NX/2;
    double v=sqrt(2.0*PKA_E/h_m[pka]), th=15.0*M_PI/180.0;
    h_vx[pka]=v*sin(th); h_vy[pka]=-v*cos(th);

    /* ---- Cell-List-Domaine bestimmen (Host) ------------------------------ */
    double x0=1e30,y0=1e30,x1=-1e30,y1=-1e30;
    for(int i=0;i<N;i++){
        if(h_px[i]<x0)x0=h_px[i]; if(h_px[i]>x1)x1=h_px[i];
        if(h_py[i]<y0)y0=h_py[i]; if(h_py[i]>y1)y1=h_py[i];
    }
    /* Grosszuegiger Rand: bei PKA_E=5000 fliegen Atome weit aus dem Gitter heraus */
    double margin=40.0*L0; x0-=margin; y0-=margin;
    double cs=RCUT;
    int ncx=(int)((x1-x0+margin)/cs)+1, ncy=(int)((y1-y0+margin)/cs)+1;
    int ncell=ncx*ncy;

    /* ---- Device-Allokation ----------------------------------------------- */
    double *d_px,*d_py,*d_vx,*d_vy,*d_fx,*d_fy,*d_m,*d_eb,*d_ea,*d_esum;
    int *d_nbr,*d_int,*d_nbr0,*d_head,*d_next;
    CUDA_CHECK(cudaMalloc(&d_px,db));   CUDA_CHECK(cudaMalloc(&d_py,db));
    CUDA_CHECK(cudaMalloc(&d_vx,db));   CUDA_CHECK(cudaMalloc(&d_vy,db));
    CUDA_CHECK(cudaMalloc(&d_fx,db));   CUDA_CHECK(cudaMalloc(&d_fy,db));
    CUDA_CHECK(cudaMalloc(&d_m,db));    CUDA_CHECK(cudaMalloc(&d_ea,db));
    CUDA_CHECK(cudaMalloc(&d_eb,sizeof(double)));    /* E_broken Skalar       */
    CUDA_CHECK(cudaMalloc(&d_esum,sizeof(double)));  /* Reduktions-Ergebnis   */
    CUDA_CHECK(cudaMalloc(&d_nbr,ib));  CUDA_CHECK(cudaMalloc(&d_int,ib));
    CUDA_CHECK(cudaMalloc(&d_nbr0,N*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_head,ncell*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_next,N*sizeof(int)));

    /* ---- H2D-Transfer (Host -> Device), gemessen mit CUDA Events --------- */
    cudaEvent_t ev0,ev1;
    CUDA_CHECK(cudaEventCreate(&ev0)); CUDA_CHECK(cudaEventCreate(&ev1));
    /* H2D-Volumen: px,py,vx,vy,m (5 double-Arrays) + nbr,int (2 int-Arrays) + nbr0 */
    size_t h2d_bytes = 5*db + 2*ib + (size_t)N*sizeof(int);
    CUDA_CHECK(cudaEventRecord(ev0));
    CUDA_CHECK(cudaMemcpy(d_px,h_px,db,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_py,h_py,db,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vx,h_vx,db,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vy,h_vy,db,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_m,h_m,db,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_nbr,h_nbr,ib,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_int,h_int,ib,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_nbr0,h_nbr0,N*sizeof(int),cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaEventRecord(ev1)); CUDA_CHECK(cudaEventSynchronize(ev1));
    float t_h2d_ms=0; CUDA_CHECK(cudaEventElapsedTime(&t_h2d_ms,ev0,ev1));

    /* Grid-Konfiguration: T Threads/Block, B Bloecke fuer N Atome, Bc Bloecke fuer ncell */
    int T=256, B=(N+T-1)/T, Bc=(ncell+T-1)/T;

    /* ---- Initiale Cell List + Anfangskraefte ----------------------------- */
    double zero=0.0;
    CUDA_CHECK(cudaMemcpy(d_eb,&zero,sizeof(double),cudaMemcpyHostToDevice));
    k_clear_cells<<<Bc,T>>>(d_head,ncell);
    k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
    k_force<<<B,T>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                     x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N);
    CUDA_CHECK(cudaGetLastError());
    /* E_broken zuruecksetzen: Aufwaerm-Bruche werden nicht zur Simulation gezaehlt */
    CUDA_CHECK(cudaMemcpy(d_eb,&zero,sizeof(double),cudaMemcpyHostToDevice));

    printf("# CUDA | N=%d | PKA_E=%.0f | dt=%g | Schritte=%d | Zellen=%d\n",
           N,PKA_E,dt,NSTEPS,ncell);

    /* ---- SELF-CHECK: GPU-Reduktion == Host-Summe ------------------------- */
    /* Korrektheitsbeweis: Energie-Summe auf GPU und CPU muessen uebereinstimmen.
     * Erlaubte relative Abweichung: < 1e-9 (numerische Prazision). */
    k_energy<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_m,d_nbr,d_int,d_nbr0,d_head,d_next,
                      x0,y0,cs,ncx,ncy,K,L0,RCUT,KREP,REPN,d_ea,N);
    CUDA_CHECK(cudaGetLastError());
    double gpu_sum = gpu_reduce_sum(d_ea, d_esum, N, T);
    /* Referenzwert: e_atom[] auf den Host kopieren und dort summieren */
    CUDA_CHECK(cudaMemcpy(h_ea,d_ea,db,cudaMemcpyDeviceToHost));
    double host_sum=0.0; for(int i=0;i<N;i++) host_sum+=h_ea[i];
    double rel = (host_sum!=0.0) ? fabs(gpu_sum-host_sum)/fabs(host_sum) : 0.0;
    printf("# [Self-Check] GPU=%.9e  Host=%.9e  rel.Diff=%.2e  -> %s\n",
           gpu_sum, host_sum, rel, (rel<1e-9)?"OK":"ABWEICHUNG!");

    /* Startenergie fuer Drift-Berechnung */
    double eb0; CUDA_CHECK(cudaMemcpy(&eb0,d_eb,sizeof(double),cudaMemcpyDeviceToHost));
    double e0 = gpu_sum + eb0;

    /* ---- Aufwaermphase (ungetaktet) -------------------------------------- */
    /* Die ersten Schritte koennen durch Takt-Einschwingens der GPU-Clocks
     * verfaelscht sein.  WARMUP Schritte werden ausgefuehrt, aber nicht gemessen. */
    for (int s=0; s<WARMUP; s++){
        k_kick_drift<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
        k_clear_cells<<<Bc,T>>>(d_head,ncell);
        k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
        k_force<<<B,T>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                         x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N);
        k_kick<<<B,T>>>(d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
    }
    CUDA_CHECK(cudaDeviceSynchronize());  /* warten bis GPU fertig */

    /* ---- Getaktete Zeitschleife (NSTEPS volle Schritte) ------------------ */
    /* cudaEventRecord startet einen "Timestamp" in der GPU-Command-Queue.
     * cudaEventElapsedTime gibt die Zeit ZWISCHEN zwei Events in ms zurueck.
     * Vorteil: misst nur die GPU-Ausfuehrungszeit, nicht Host-Overhead. */
    cudaEvent_t lt0,lt1;
    CUDA_CHECK(cudaEventCreate(&lt0)); CUDA_CHECK(cudaEventCreate(&lt1));
    CUDA_CHECK(cudaEventRecord(lt0));
    for (int s=0; s<NSTEPS; s++){
        /* Velocity-Verlet: Kick + Drift, Cell List, Force, Kick */
        k_kick_drift<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
        k_clear_cells<<<Bc,T>>>(d_head,ncell);
        k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
        k_force<<<B,T>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                         x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N);
        k_kick<<<B,T>>>(d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
        /* Energie-Log (nur alle LOGEVERY Schritte; Energieberechnung ist teuer!) */
        if ((s+1)%LOGEVERY==0){
            k_energy<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_m,d_nbr,d_int,d_nbr0,d_head,d_next,
                              x0,y0,cs,ncx,ncy,K,L0,RCUT,KREP,REPN,d_ea,N);
            double es = gpu_reduce_sum(d_ea,d_esum,N,T);
            double eb; CUDA_CHECK(cudaMemcpy(&eb,d_eb,sizeof(double),cudaMemcpyDeviceToHost));
            double et=es+eb;
            printf("#   Schritt %6d  E_total=%.4f  drift=%.4f%%\n",
                   s+1, et, 100.0*(et-e0)/e0);
        }
    }
    CUDA_CHECK(cudaEventRecord(lt1)); CUDA_CHECK(cudaEventSynchronize(lt1));
    float t_loop_ms=0; CUDA_CHECK(cudaEventElapsedTime(&t_loop_ms,lt0,lt1));

    /* ---- Isolierter Force-Kernel-Microbench (fuer Bandbreite) ------------ */
    /* FORCE_REP Wiederholungen des Force-Kernels mit unveraenderter Cell List.
     * Gibt ein stabiles Zeitmittel fuer die Kernel-Dauer ohne Cell-List-Overhead.
     * Aus der Durchschnittsdauer + Modell-Bandbreite BYTES_PER_ATOM ergibt
     * sich eine Schaetzung der Speicherbandbreiten-Auslastung (als Sanity-Check;
     * autoritativer Wert kommt von ncu SpeedOfLight). */
    k_clear_cells<<<Bc,T>>>(d_head,ncell);
    k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
    CUDA_CHECK(cudaDeviceSynchronize());
    cudaEvent_t f0,f1;
    CUDA_CHECK(cudaEventCreate(&f0)); CUDA_CHECK(cudaEventCreate(&f1));
    CUDA_CHECK(cudaEventRecord(f0));
    for (int r=0;r<FORCE_REP;r++){
        k_force<<<B,T>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                         x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N);
    }
    CUDA_CHECK(cudaEventRecord(f1)); CUDA_CHECK(cudaEventSynchronize(f1));
    float t_force_ms=0; CUDA_CHECK(cudaEventElapsedTime(&t_force_ms,f0,f1));
    double t_force_call_s = (t_force_ms/1000.0)/FORCE_REP;  /* Sekunden pro Aufruf */

    /* ---- D2H-Transfer (Device -> Host), gemessen ------------------------- */
    CUDA_CHECK(cudaEventRecord(ev0));
    CUDA_CHECK(cudaMemcpy(h_ea,d_ea,db,cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaEventRecord(ev1)); CUDA_CHECK(cudaEventSynchronize(ev1));
    float t_d2h_ms=0; CUDA_CHECK(cudaEventElapsedTime(&t_d2h_ms,ev0,ev1));

    /* ---- Abgeleitete Kennzahlen ------------------------------------------ */
    double t_loop_s    = t_loop_ms/1000.0;
    double per_step_us = (t_loop_ms*1000.0)/NSTEPS;           /* us pro Schritt      */
    double matom_s     = ((double)N*NSTEPS)/t_loop_s/1e6;     /* Matom/s (Durchsatz) */
    /* Modell-Bandbreite: Bytes pro Atom pro Force-Kernel / Ausfuehrungszeit */
    double bw_model    = (BYTES_PER_ATOM*(double)N)/t_force_call_s/1e9;  /* GB/s   */
    double peak_pct    = 100.0*bw_model/PEAK_BW_GBS;          /* % der Peak-BW      */
    double h2d_gbs     = (h2d_bytes/1e9)/(t_h2d_ms/1000.0);  /* H2D-Transferrate   */
    double d2h_gbs     = (db/1e9)/(t_d2h_ms/1000.0);         /* D2H-Transferrate   */
    double sp_1core    = T_CPU_1CORE/t_loop_s;                /* GPU vs. 1 CPU-Kern */
    double t_cpu_node  = T_CPU_1CORE/SPEEDUP_64;              /* ~22,5 s (64 Kerne) */
    double sp_node     = t_cpu_node/t_loop_s;                 /* GPU vs. voller Knoten */

    /* Energie-Drift am Ende der Simulation */
    double eb_end; CUDA_CHECK(cudaMemcpy(&eb_end,d_eb,sizeof(double),cudaMemcpyDeviceToHost));
    k_energy<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_m,d_nbr,d_int,d_nbr0,d_head,d_next,
                      x0,y0,cs,ncx,ncy,K,L0,RCUT,KREP,REPN,d_ea,N);
    double es_end = gpu_reduce_sum(d_ea,d_esum,N,T);
    double drift_end = 100.0*((es_end+eb_end)-e0)/e0;

    /* ---- Ergebnis-Block ausgeben ----------------------------------------- */
    printf("\n");
    printf("# =================== ERGEBNIS-BLOCK ===================\n");
    printf("# Problemgroesse        : N = %d Atome (%dx%d), %d Schritte\n",N,NX,NY,NSTEPS);
    printf("# Energie-Drift (Ende)  : %.4f %%\n", drift_end);
    printf("# --- Laufzeit (GPU, ohne Transfer) ---\n");
    printf("# Zeitschleife gesamt   : %.3f s\n", t_loop_s);
    printf("# pro Schritt           : %.1f us\n", per_step_us);
    printf("# Durchsatz             : %.1f Matom/s\n", matom_s);
    printf("# --- Force-Kernel (isoliert, %d Wdh.) ---\n", FORCE_REP);
    printf("# pro Aufruf            : %.1f us\n", t_force_call_s*1e6);
    printf("# Bandbreite (MODELL)   : %.1f GB/s  (= %.1f %% der A100-Peak %.0f GB/s)\n",
           bw_model, peak_pct, PEAK_BW_GBS);
    printf("#   -> Autoritativ: ncu SpeedOfLight (DRAM-Throughput-Messung).\n");
    printf("# --- Transfer ---\n");
    printf("# H2D (Setup, %.0f MB)   : %.2f ms  (%.0f GB/s)\n", h2d_bytes/1e6, t_h2d_ms, h2d_gbs);
    printf("# D2H (e_atom, %.0f MB)  : %.2f ms  (%.0f GB/s)\n", db/1e6, t_d2h_ms, d2h_gbs);
    printf("# --- CPU vs. GPU (gleiche Konfig) ---\n");
    printf("# CPU 1 Kern   T(1)      : %.1f s   -> Speedup GPU = %.1fx\n", T_CPU_1CORE, sp_1core);
    printf("# CPU 64 Kerne T(1)/51,2 : %.1f s   -> Speedup GPU = %.2fx\n", t_cpu_node, sp_node);
    printf("# =======================================================\n");

    /* ---- Aufraeumen -------------------------------------------------------- */
    cudaFree(d_px);cudaFree(d_py);cudaFree(d_vx);cudaFree(d_vy);
    cudaFree(d_fx);cudaFree(d_fy);cudaFree(d_m);cudaFree(d_ea);
    cudaFree(d_eb);cudaFree(d_esum);
    cudaFree(d_nbr);cudaFree(d_int);cudaFree(d_nbr0);
    cudaFree(d_head);cudaFree(d_next);
    free(h_px);free(h_py);free(h_vx);free(h_vy);free(h_m);free(h_ea);
    free(h_nbr);free(h_int);free(h_cnt);free(h_nbr0);
    return 0;
}
