/* ============================================================================
 *  cascade_cuda_bench_opt.cu  --  Zwei GPU-Optimierungen isoliert testen (Launch
 *                                  Bounds, CUDA Graphs) auf Basis von cascade_cuda_bench.cu
 * ----------------------------------------------------------------------------
 *
 *  ZWECK: WAS TESTET DIESE DATEI?
 *  --------------------------------
 *  cascade_cuda_bench.cu liefert die Baseline (float/double, N-Sweep, Breakdown).
 *  Diese Datei ist eine Kopie mit IDENTISCHER Physik und IDENTISCHEN Kerneln,
 *  die zwei UNABHAENGIGE Low-Level-Optimierungen per Compile-Flag zuschaltet,
 *  um ihren Effekt auf die Laufzeit isoliert zu messen:
 *
 *  (1) LAUNCH BOUNDS -- Register limitieren, um Occupancy zu erhoehen
 *  -----------------------------------------------------------------------
 *  Flag: -DUSE_LB [-DLB_MINBLK=N]  (Standard LB_MINBLK=6, LB_T=256)
 *  __launch_bounds__(LB_T, LB_MINBLK) weist ptxas an, k_force so zu
 *  kompilieren, dass mindestens LB_MINBLK Bloecke a LB_T Threads GLEICHZEITIG
 *  auf einem SM laufen koennen.  Der Compiler MUSS dafuer die Register pro
 *  Thread reduzieren (weniger Register/Thread * mehr gleichzeitige Threads
 *  <= feste Registerdatei-Groesse pro SM).
 *  Gemessen (ncu, A100, fp64, N~2 Mio.):
 *    ohne LB:  62 Register/Thread -> 45 % Occupancy
 *    mit  LB:  40 Register/Thread -> 66 % Occupancy
 *  WICHTIG: Blockgroesse T beim Programmstart MUSS gleich LB_T (256) sein,
 *  sonst schlaegt der Kernel-Start fehl (Launch-Bounds sind ein Vertrag
 *  zwischen Compile-Zeit-Flag und Laufzeit-Blockgroesse).
 *
 *  (2) CUDA GRAPHS -- Kernel-Start-Overhead pro Schritt eliminieren
 *  -----------------------------------------------------------------------
 *  Flag: -DUSE_GRAPH
 *  Ohne Graph startet jeder Zeitschritt 5 Kernel EINZELN vom Host aus
 *  (kick_drift, clear_cells, build_cells, force, kick) -- jeder Start hat
 *  CPU-seitigen Launch-Overhead (~5 us), der bei sehr kurzen Kerneln relativ
 *  ins Gewicht fallen kann.  Mit CUDA Graphs wird die 5-Kernel-Sequenz EINMAL
 *  per Stream-Capture aufgezeichnet und pro Schritt mit einem einzigen
 *  cudaGraphLaunch() abgespielt -- der Host schickt nur noch einen Aufruf
 *  statt fuenf, das Treiber-Overhead pro Schritt sinkt.
 *
 *  Beide Flags sind UNABHAENGIG und LOMBINIERBAR (4 Varianten: base, lb,
 *  graph, lbgraph) -- vgl. plot_cuda_opt.py, das alle vier gegenueberstellt.
 *
 *  ERGEBNIS (A100, fp64, N~2 Mio., siehe cuda_opt_results.csv):
 *    variant   Register  Occupancy  k_force[us]  Durchsatz[Matom/s]
 *    base         62       45.1 %      308.8           3807.6
 *    lb           40       66.3 %      313.1           3759.5
 *    graph        62       45.1 %      309.3           3846.7
 *    lbgraph      40       66.3 %      306.8           3829.9
 *  Kernaussage: Occupancy laesst sich fast verdoppeln (45%->66%), aber der
 *  Durchsatz bleibt in allen vier Varianten INNERHALB VON <2% FLACH.  Das
 *  bestaetigt, dass k_force NICHT occupancy-limitiert ist (mehr gleichzeitige
 *  Warps helfen nicht, Speicherlatenzen zu verstecken, weil schon vorher genug
 *  Warps in Flight waren) und dass Launch-Overhead bei NSTEPS=3000 x 5 Kerneln
 *  vernachlaessigbar gegenueber der eigentlichen Kernel-Laufzeit ist. Die
 *  A100-Auslastung (~45% Compute, siehe cascade_cuda.cu) ist ein echtes
 *  Speicherlatenz-/Bandbreiten-Limit, kein Occupancy-Artefakt.
 *
 *  Kompilierung (4 Varianten):
 *    nvcc -O3 -arch=sm_80                              -o bench_base    cascade_cuda_bench_opt.cu
 *    nvcc -O3 -arch=sm_80 -DUSE_LB                     -o bench_lb      cascade_cuda_bench_opt.cu
 *    nvcc -O3 -arch=sm_80 -DUSE_GRAPH                  -o bench_graph   cascade_cuda_bench_opt.cu
 *    nvcc -O3 -arch=sm_80 -DUSE_LB -DUSE_GRAPH         -o bench_lbgraph cascade_cuda_bench_opt.cu
 *  Aufruf (wie cascade_cuda_bench.cu):
 *    ./bench_base 1999396 256 3000 3    # 1414^2 Atome, 256 Threads/Block, 3000 Schritte, 3 Wdh.
 * ==========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

/* ============================================================================
 *  PRAEZISIONS-AUSWAHL: double (Standard) oder float (-DUSE_FP32)
 * ---------------------------------------------------------------------------
 *  Identisch zu cascade_cuda_bench.cu: "typedef real" ist ein C-Typ-Alias.
 *  R(x) castet Literale sauber auf die gewaehlte Praezision; ohne R() wuerden
 *  Ausdruecke wie "0.5*dt*fx[i]" die 0.5 als double interpretieren, was bei
 *  real=float zu einem Typ-Mismatch fuehren wuerde.
 *  SELFCHK_TOL: Toleranz fuer den Self-Check (GPU-Summe == Host-Summe);
 *  float hat nur ~7 signifikante Dezimalstellen -> groessere Toleranz.
 * ========================================================================== */
#ifdef USE_FP32
typedef float  real;
#define PREC_NAME "float (fp32)"
#define SELFCHK_TOL 1e-2
#else
typedef double real;
#define PREC_NAME "double (fp64)"
#define SELFCHK_TOL 1e-9
#endif
#define R(x) ((real)(x))            /* Literal sauber auf die Praezision casten */

#define MAXNBR 6      /* Max. Nachbarn im Dreiecksgitter: 6 pro Innenatom */

/* CUDA_CHECK: bei Fehler Datei/Zeile/Fehlertext ausgeben und exit(1) (vgl. cascade_cuda.cu) */
#define CUDA_CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){       \
    fprintf(stderr,"CUDA %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e));\
    exit(1);} } while(0)

/* CPU-Referenz (ECHT gemessen, fp64 @ 1414^2, 3000 Schritte, vgl. cascade_cuda.cu) */
#define T_CPU_1CORE   1152.9    /* s, 1 EPYC-Kern                          */
#define SPEEDUP_64    51.2      /* gemessener 64-Kern-MPI-Speedup          */
#define PEAK_BW_GBS   1555.0    /* GB/s, A100 40GB HBM2 Peak               */
#define BYTES_PER_ATOM 76.0     /* B/Atom (MODELL, nicht ncu-autoritativ)  */

/* ============================================================================
 *  DEVICE-FUNKTIONEN  (Abstoss-Kraft und -Potential, templated ueber 'real')
 * ---------------------------------------------------------------------------
 *  Identische Formeln wie in cascade_cuda.cu / cascade_serial.c, aber ueber
 *  'real' parametrisiert -- bei -DUSE_FP32 kompiliert derselbe Code in float.
 * ========================================================================== */
__device__ inline real rep_force(real r,real RCUT,real KREP,real REPN){
    return KREP*(pow(RCUT/r,REPN) - R(1.0));
}
__device__ inline real rep_pot(real r,real RCUT,real KREP,real REPN){
    real n=REPN;
    return KREP*((RCUT - pow(RCUT,n)*pow(r,R(1.0)-n))/(R(1.0)-n) - (RCUT - r));
}

/* ============================================================================
 *  KERNEL: k_clear_cells  --  Cell-List-Zellen leeren (Kernel 1 von 2)
 *  Setzt alle Eintraege von head[] auf -1 (leer). Ein Thread pro Zelle.
 * ========================================================================== */
__global__ void k_clear_cells(int* head, int ncell){
    int c = blockIdx.x*blockDim.x + threadIdx.x;
    if (c < ncell) head[c] = -1;
}

/* ============================================================================
 *  KERNEL: k_build_cells  --  Atome in Cell List eintragen (Kernel 2 von 2)
 *  atomicExch(&head[c], i): setzt head[c]=i und liefert den alten Wert
 *  atomar zurueck -> next[i] zeigt auf den vorherigen Listenkopf. Sicher,
 *  auch wenn tausende Threads gleichzeitig in dieselbe Zelle schreiben
 *  (Details siehe cascade_cuda.cu).
 * ========================================================================== */
__global__ void k_build_cells(const real* px, const real* py,
                              int* head, int* next,
                              real x0, real y0, real cs,
                              int ncx, int ncy, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    int cx=(int)((px[i]-x0)/cs), cy=(int)((py[i]-y0)/cs);
    if(cx<0)cx=0; if(cx>=ncx)cx=ncx-1;
    if(cy<0)cy=0; if(cy>=ncy)cy=ncy-1;
    int c = cy*ncx + cx;
    next[i] = atomicExch(&head[c], i);   /* atomar: i wird neuer Listenkopf */
}

/* ============================================================================
 *  KERNEL: k_force  --  Kraftberechnung (Gather, Haupt-Kernel, ~90% der Zeit)
 * ----------------------------------------------------------------------------
 *  Physik/Struktur identisch zu cascade_cuda_bench.cu / cascade_cuda.cu:
 *   (a) Federkraefte gebundener Nachbarn + Bruchdetektion
 *   (b) Abstoss-Kraefte via Cell List (voller 3x3-Zellenstern)
 *  Gather-Formulierung: Thread i schreibt NUR fx[i],fy[i] (kein atomicAdd
 *  auf Kraefte noetig), liest dafuer px[j],py[j] von allen Nachbarn/Zellen.
 *
 *  LAUNCH-BOUNDS-EXPERIMENT (dieser Kernel ist das einzige Ziel, weil er
 *  ~90% der Laufzeit dominiert und mit ~62 Register/Thread bei fp64 den
 *  groessten Registerdruck hat):
 *  -----------------------------------------------------------------------
 *  Ohne Flag kompiliert ptxas k_force mit so vielen Registern, wie der
 *  Compiler fuer beste Performance PRO THREAD fuer sinnvoll haelt (~62 bei
 *  fp64) -- das begrenzt, wie viele Bloecke GLEICHZEITIG auf einem SM aktiv
 *  sein koennen (Registerdatei pro SM ist fix, z.B. 65536 32-Bit-Register
 *  auf der A100), was die Occupancy (Verhaeltnis aktiver zu maximal
 *  moeglichen Warps) senkt: gemessen 45%.
 *
 *  __launch_bounds__(LB_T, LB_MINBLK) ist ein VERTRAG mit ptxas:
 *  "kompiliere so, dass mindestens LB_MINBLK Bloecke a LB_T Threads
 *  gleichzeitig auf einem SM passen." ptxas erzwingt das, indem es die
 *  Register/Thread SENKT (hier: 62 -> 40 Register), notfalls auf Kosten von
 *  Register-Spilling in den (langsameren) lokalen Speicher. Gemessene Folge:
 *  Occupancy steigt auf 66%, ABER die Kernel-Laufzeit aendert sich kaum
 *  (308.8 us -> 313.1 us, sogar minimal langsamer) -- k_force ist also
 *  bandbreiten-/latenzgebunden, nicht durch zu wenige gleichzeitige Warps
 *  limitiert (mehr Occupancy hilft nicht, mehr Speicherzugriffe gleichzeitig
 *  "in Flight" zu haben, wenn die Bandbreite ohnehin schon ausgelastet ist).
 *
 *  WICHTIG: LB_T MUSS mit der tatsaechlichen Kernel-Start-Blockgroesse T
 *  uebereinstimmen (Standard 256) -- sonst verletzt der Start den Vertrag
 *  und CUDA bricht mit einem Launch-Fehler ab.
 * ========================================================================== */
#ifndef LB_T
#define LB_T 256
#endif
#if defined(USE_LB) && !defined(LB_MINBLK)
#define LB_MINBLK 6
#endif
__global__ void
#ifdef USE_LB
__launch_bounds__(LB_T, LB_MINBLK)
#endif
k_force(const real* px, const real* py,
                        real* fx, real* fy,
                        const int* nbr, int* intact, const int* nbr0,
                        const int* head, const int* next,
                        real x0, real y0, real cs, int ncx, int ncy,
                        real K, real L0, real maxstr,
                        real RCUT, real KREP, real REPN,
                        real* ebroken, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    real xi=px[i], yi=py[i];
    real fxi=R(0.0), fyi=R(0.0), brk=R(0.0);
    int base=i*MAXNBR, cnt=nbr0[i];

    /* (a) Federkraefte gebundener Nachbarn + Bruchdetektion (r > L0*maxstr) */
    for (int e=0;e<cnt;e++){
        if(!intact[base+e]) continue;
        int j=nbr[base+e];
        real dx=px[j]-xi, dy=py[j]-yi, r=sqrt(dx*dx+dy*dy);
        if (r > L0*maxstr){ intact[base+e]=0; brk += R(0.5)*(R(0.5)*K*(r-L0)*(r-L0)); continue; }
        real f=K*(r-L0)/r;  fxi += f*dx;  fyi += f*dy;
    }

    /* (b) Abstoss-Kraefte via Cell List (voller 3x3-Zellenstern) */
    real rc2=RCUT*RCUT;
    int cx=(int)((xi-x0)/cs), cy=(int)((yi-y0)/cs);
    if(cx<0)cx=0; if(cx>=ncx)cx=ncx-1; if(cy<0)cy=0; if(cy>=ncy)cy=ncy-1;
    for (int oy=-1; oy<=1; oy++) for (int ox=-1; ox<=1; ox++){
        int nx=cx+ox, ny=cy+oy;
        if(nx<0||nx>=ncx||ny<0||ny>=ncy) continue;
        int c=ny*ncx+nx;
        for (int j=head[c]; j!=-1; j=next[j]){
            if (j==i) continue;
            real dx=xi-px[j], dy=yi-py[j];
            real r2=dx*dx+dy*dy;
            if (r2>=rc2 || r2==R(0.0)) continue;
            int bonded=0;
            for(int e=0;e<cnt;e++){ if(intact[base+e] && nbr[base+e]==j){bonded=1;break;} }
            if (bonded) continue;
            real r=sqrt(r2), f=rep_force(r,RCUT,KREP,REPN)/r;
            fxi += f*dx;  fyi += f*dy;
        }
    }
    fx[i]=fxi;  fy[i]=fyi;
    if (brk!=R(0.0)) atomicAdd(ebroken, brk);   /* selten (nur bei Bruch), daher billig */
}

/* ============================================================================
 *  KERNEL: k_kick_drift und k_kick  --  Velocity-Verlet Integration
 *  (identisch zu cascade_cuda_bench.cu / cascade_cuda.cu, ueber 'real')
 *  kick_drift: v += 0.5*dt*F/m (halber Kick), dann x += v*dt (Drift)
 *  kick:       v += 0.5*dt*F/m (zweiter halber Kick, mit NEUEN Kraeften)
 * ========================================================================== */
__global__ void k_kick_drift(real* px,real* py,real* vx,real* vy,
                             const real* fx,const real* fy,const real* m,
                             real dt,int N){
    int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=N)return;
    vx[i]+=R(0.5)*dt*fx[i]/m[i]; vy[i]+=R(0.5)*dt*fy[i]/m[i];
    px[i]+=dt*vx[i];             py[i]+=dt*vy[i];
}
__global__ void k_kick(real* vx,real* vy,const real* fx,const real* fy,
                       const real* m,real dt,int N){
    int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=N)return;
    vx[i]+=R(0.5)*dt*fx[i]/m[i]; vy[i]+=R(0.5)*dt*fy[i]/m[i];
}

/* ============================================================================
 *  KERNEL: k_energy  --  Energie-Beitrag pro Atom (Gather, Halbes-Zaehlen)
 * ----------------------------------------------------------------------------
 *  E_kin komplett (Thread i besitzt Atom i), E_spring/E_rep je mit Faktor 0.5,
 *  weil die Gather-Formulierung jede Wechselwirkung (i,j) UND (j,i) zaehlt
 *  (Details/Begruendung siehe cascade_cuda.cu). e_atom[] wird von k_reduce
 *  aufsummiert.
 * ========================================================================== */
__global__ void k_energy(const real* px,const real* py,
                         const real* vx,const real* vy,const real* m,
                         const int* nbr,const int* intact,const int* nbr0,
                         const int* head,const int* next,
                         real x0,real y0,real cs,int ncx,int ncy,
                         real K,real L0,real RCUT,real KREP,real REPN,
                         real* e_atom,int N){
    int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=N)return;
    real xi=px[i], yi=py[i];
    real e = R(0.5)*m[i]*(vx[i]*vx[i]+vy[i]*vy[i]);   /* kinetische Energie (ganz) */
    int base=i*MAXNBR, cnt=nbr0[i];
    for (int k=0;k<cnt;k++){
        if(!intact[base+k]) continue;
        int j=nbr[base+k];
        real dx=px[j]-xi,dy=py[j]-yi,r=sqrt(dx*dx+dy*dy);
        e += R(0.5)*(R(0.5)*K*(r-L0)*(r-L0));    /* halbe Federpotenzial-Energie */
    }
    real rc2=RCUT*RCUT;
    int cx=(int)((xi-x0)/cs), cy=(int)((yi-y0)/cs);
    if(cx<0)cx=0; if(cx>=ncx)cx=ncx-1; if(cy<0)cy=0; if(cy>=ncy)cy=ncy-1;
    for (int oy=-1;oy<=1;oy++) for (int ox=-1;ox<=1;ox++){
        int nx=cx+ox,ny=cy+oy; if(nx<0||nx>=ncx||ny<0||ny>=ncy)continue;
        int c=ny*ncx+nx;
        for (int j=head[c]; j!=-1; j=next[j]){
            if(j==i) continue;
            real dx=xi-px[j],dy=yi-py[j],r2=dx*dx+dy*dy;
            if(r2>=rc2||r2==R(0.0)) continue;
            int bonded=0; for(int k=0;k<cnt;k++){ if(intact[base+k]&&nbr[base+k]==j){bonded=1;break;} }
            if(bonded) continue;
            e += R(0.5)*rep_pot(sqrt(r2),RCUT,KREP,REPN);  /* halbe Abstoss-Energie */
        }
    }
    e_atom[i]=e;
}

/* ============================================================================
 *  KERNEL: k_reduce  --  Parallele Baumreduktion (identisch zu cascade_cuda.cu)
 * ----------------------------------------------------------------------------
 *  Grid-Stride-Load (Schritt 1) + Shared-Memory-Baumreduktion (Schritt 2) +
 *  atomicAdd des Block-Ergebnisses aufs globale out[0] (Schritt 3).
 *  Bei float ist sdata[] nur halb so gross wie bei double -- mehr Platz im
 *  Shared Memory, kann bei knappen Shared-Memory-Limits die Occupancy heben.
 * ========================================================================== */
__global__ void k_reduce(const real* __restrict__ in, real* out, int N){
    extern __shared__ real sdata[];
    int tid=threadIdx.x, gid=blockIdx.x*blockDim.x+tid, stride=blockDim.x*gridDim.x;
    real local=R(0.0);
    for (int k=gid; k<N; k+=stride) local += in[k];     /* Grid-Stride-Load     */
    sdata[tid]=local; __syncthreads();
    for (int s=blockDim.x/2; s>0; s>>=1){               /* Baumreduktion        */
        if(tid<s) sdata[tid]+=sdata[tid+s]; __syncthreads();
    }
    if (tid==0) atomicAdd(out, sdata[0]);               /* Block -> Ergebnis    */
}

/* gpu_reduce_sum: e_atom[] auf der GPU summieren, Skalar zum Host zurueckgeben */
static real gpu_reduce_sum(const real* d_in, real* d_scalar, int N, int T){
    CUDA_CHECK(cudaMemset(d_scalar,0,sizeof(real)));
    int RB=(N+T-1)/T; if(RB>1024)RB=1024;
    k_reduce<<<RB,T,T*sizeof(real)>>>(d_in,d_scalar,N);
    CUDA_CHECK(cudaGetLastError());
    real s; CUDA_CHECK(cudaMemcpy(&s,d_scalar,sizeof(real),cudaMemcpyDeviceToHost));
    return s;
}

/* ============================================================================
 *  STATISTIK-HELFER
 * ----------------------------------------------------------------------------
 *  cmp_d: Vergleichsfunktion fuer qsort (aufsteigende Reihenfolge).
 *  median: sortiert v[] (mutiert das Array!) und liefert den Median
 *          (gerades n: Mittelwert der zwei mittleren Werte).
 * ========================================================================== */
static int cmp_d(const void*a,const void*b){ double x=*(const double*)a,y=*(const double*)b; return (x>y)-(x<y); }
static double median(double* v,int n){ qsort(v,n,sizeof(double),cmp_d); return (n&1)?v[n/2]:R(0.5)*(v[n/2-1]+v[n/2]); }

/* add_nb: Bindung in Nachbarliste eintragen (Host, identisch zu cascade_cuda.cu) */
static void add_nb(int* nbr,int* c,int a,int b){ nbr[a*MAXNBR+c[a]++]=b; nbr[b*MAXNBR+c[b]++]=a; }

/* Breakdown-Sammler: Kernel-Namen + isoliert gemessene Zeiten (us/Aufruf), siehe BENCH-Makro unten */
static const char* bk_name[8]; static double bk_us[8]; static int bk_n=0;

/* ============================================================================
 *  MAIN  --  Hauptprogramm
 * ----------------------------------------------------------------------------
 *  Kommandozeilenargumente (wie cascade_cuda_bench.cu):
 *    argv[1] N_ziel (Standard 1414*1414), argv[2] blocksize (256),
 *    argv[3] n_steps (3000), argv[4] n_rep (3, max 64)
 *
 *  Ablauf:
 *   1. Parameter parsen, Gitter aufbauen (Host), Cell-List-Domain bestimmen
 *   2. Device-Allokation, H2D-Transfer (gemessen)
 *   3. Snapshot der Anfangszustaende (fuer Reset zwischen Wiederholungen)
 *   4. Self-Check: GPU-Reduktionssumme == Host-Summe
 *   5. [nur -DUSE_GRAPH] Schrittsequenz einmal per Stream-Capture aufnehmen
 *   6. NREP Wiederholungen der Headline-Zeitschleife (Median berichten)
 *   7. Kernel-Breakdown: jeder Kernel isoliert FORCE_REP mal gemessen
 *   8. D2H-Transfer messen, Ergebnis-Block ausgeben (inkl. CSV-Zeile)
 * ========================================================================== */
int main(int argc, char** argv){
    /* ---- Kommandozeile parsen --------------------------------------------- */
    long N_target = (argc>1)? atol(argv[1]) : 1414L*1414L;
    int  T        = (argc>2)? atoi(argv[2]) : 256;    /* Threads/Block; bei USE_LB MUSS T==LB_T sein */
    int  NSTEPS   = (argc>3)? atoi(argv[3]) : 3000;
    int  NREP     = (argc>4)? atoi(argv[4]) : 3;
    if (T<32) T=32; if (NREP<1) NREP=1;
    /* N_ziel auf naechstes Quadrat runden: NX=NY=round(sqrt(N_ziel)) */
    int side = (int)(sqrt((double)N_target)+0.5); if (side<8) side=8;
    const int NX=side, NY=side, N=NX*NY;

    /* ---- Physik-Parameter (wie params.ini Scaling-Config) ----------------- */
    const real L0=R(1.0), K=R(100.0), MAXSTR=R(1.15);
    const real RCUT=R(0.9), KREP=R(400.0), REPN=R(12.0);
    const real dt=R(2e-4), PKA_E=R(5000.0);
    const int  LOGEVERY = (NSTEPS>=1000)?500:NSTEPS;
    const int  WARMUP=50, FORCE_REP=200;       /* Aufwaerm- und Breakdown-Wdh. */
    const real DY=L0*R(sqrt(3.0))/R(2.0);
    size_t db=(size_t)N*sizeof(real), ib=(size_t)N*MAXNBR*sizeof(int);

    printf("# ===== cascade_cuda_bench | %s | N=%d (%dx%d) | T=%d | Schritte=%d | rep=%d =====\n",
           PREC_NAME, N, NX, NY, T, NSTEPS, NREP);

    /* ---- Gitter aufbauen (Host): Positionen, Bindungen, PKA --------------- */
    real *h_px=(real*)malloc(db),*h_py=(real*)malloc(db);
    real *h_vx=(real*)calloc(N,sizeof(real)),*h_vy=(real*)calloc(N,sizeof(real)),*h_m=(real*)malloc(db);
    real *h_ea=(real*)malloc(db);
    int *h_nbr=(int*)malloc(ib),*h_int=(int*)malloc(ib);
    int *h_cnt=(int*)calloc(N,sizeof(int)),*h_nbr0=(int*)malloc(N*sizeof(int));
    /* Positionen: Dreiecksgitter */
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){int id=j*NX+i;
        h_px[id]=i*L0+(j&1)*(L0*R(0.5)); h_py[id]=j*DY; h_m[id]=R(1.0);}
    /* Bindungen: rechts und oben (jede einmal) */
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){int a=j*NX+i;
        if(i+1<NX)add_nb(h_nbr,h_cnt,a,j*NX+i+1);
        if(j+1<NY){int li,ri; if((j&1)==0){li=i-1;ri=i;}else{li=i;ri=i+1;}
            if(li>=0&&li<NX)add_nb(h_nbr,h_cnt,a,(j+1)*NX+li);
            if(ri>=0&&ri<NX)add_nb(h_nbr,h_cnt,a,(j+1)*NX+ri);}}
    /* intact[]-Array aufbauen und nbr0[] (Anfangs-Nachbarzahl) speichern */
    for(int i=0;i<N;i++){h_nbr0[i]=h_cnt[i];
        for(int e=0;e<MAXNBR;e++)h_int[i*MAXNBR+e]=(e<h_cnt[i])?1:0;}
    /* PKA in der Gittermitte setzen */
    int pka=(NY/2)*NX+NX/2; real v=sqrt(R(2.0)*PKA_E/h_m[pka]),th=R(15.0*M_PI/180.0);
    h_vx[pka]=v*sin(th); h_vy[pka]=-v*cos(th);

    /* ---- Cell-List-Domain (Host) ------------------------------------------ */
    real x0=R(1e30),y0=R(1e30),x1=R(-1e30),y1=R(-1e30);
    for(int i=0;i<N;i++){ if(h_px[i]<x0)x0=h_px[i]; if(h_px[i]>x1)x1=h_px[i];
                          if(h_py[i]<y0)y0=h_py[i]; if(h_py[i]>y1)y1=h_py[i]; }
    real margin=R(40.0)*L0; x0-=margin; y0-=margin;
    real cs=RCUT; int ncx=(int)((x1-x0+margin)/cs)+1, ncy=(int)((y1-y0+margin)/cs)+1;
    int ncell=ncx*ncy;

    /* ---- Device-Allokation ------------------------------------------------- */
    real *d_px,*d_py,*d_vx,*d_vy,*d_fx,*d_fy,*d_m,*d_eb,*d_ea,*d_esum;
    int *d_nbr,*d_int,*d_nbr0,*d_head,*d_next;
    CUDA_CHECK(cudaMalloc(&d_px,db));CUDA_CHECK(cudaMalloc(&d_py,db));
    CUDA_CHECK(cudaMalloc(&d_vx,db));CUDA_CHECK(cudaMalloc(&d_vy,db));
    CUDA_CHECK(cudaMalloc(&d_fx,db));CUDA_CHECK(cudaMalloc(&d_fy,db));
    CUDA_CHECK(cudaMalloc(&d_m,db)); CUDA_CHECK(cudaMalloc(&d_ea,db));
    CUDA_CHECK(cudaMalloc(&d_eb,sizeof(real)));CUDA_CHECK(cudaMalloc(&d_esum,sizeof(real)));
    CUDA_CHECK(cudaMalloc(&d_nbr,ib));CUDA_CHECK(cudaMalloc(&d_int,ib));
    CUDA_CHECK(cudaMalloc(&d_nbr0,N*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_head,ncell*sizeof(int)));CUDA_CHECK(cudaMalloc(&d_next,N*sizeof(int)));

    /* ---- H2D-Transfer (gemessen) -------------------------------------------- */
    cudaEvent_t ev0,ev1; CUDA_CHECK(cudaEventCreate(&ev0)); CUDA_CHECK(cudaEventCreate(&ev1));
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

    int B=(N+T-1)/T, Bc=(ncell+T-1)/T;

    /* ---- Snapshot speichern (fuer Reset zwischen Wiederholungen) ----------- */
    /* Jede Wdh. muss exakt denselben Anfangszustand haben, damit die Laufzeiten
     * vergleichbar sind. Ein Device->Device-Memcpy ist sehr schnell (~300 GB/s). */
    real *d_px0,*d_py0,*d_vx0,*d_vy0,*d_int_dummy; (void)d_int_dummy;
    CUDA_CHECK(cudaMalloc(&d_px0,db));CUDA_CHECK(cudaMalloc(&d_py0,db));
    CUDA_CHECK(cudaMalloc(&d_vx0,db));CUDA_CHECK(cudaMalloc(&d_vy0,db));
    int *d_int0; CUDA_CHECK(cudaMalloc(&d_int0,ib));
    CUDA_CHECK(cudaMemcpy(d_px0,d_px,db,cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_py0,d_py,db,cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_vx0,d_vx,db,cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_vy0,d_vy,db,cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_int0,d_int,ib,cudaMemcpyDeviceToDevice));

    /* ---- Self-Check: GPU-Reduktionssumme == Host-Summe --------------------- */
    real zero=R(0.0); CUDA_CHECK(cudaMemcpy(d_eb,&zero,sizeof(real),cudaMemcpyHostToDevice));
    k_clear_cells<<<Bc,T>>>(d_head,ncell);
    k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
    k_energy<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_m,d_nbr,d_int,d_nbr0,d_head,d_next,
                      x0,y0,cs,ncx,ncy,K,L0,RCUT,KREP,REPN,d_ea,N);
    CUDA_CHECK(cudaGetLastError());
    real gpu_sum=gpu_reduce_sum(d_ea,d_esum,N,T);
    CUDA_CHECK(cudaMemcpy(h_ea,d_ea,db,cudaMemcpyDeviceToHost));
    /* Host-Summe als double akkumuliert (vermeidet float-Akkumulierungsfehler) */
    double host_sum=0.0; for(int i=0;i<N;i++) host_sum+=(double)h_ea[i];
    double rel=(host_sum!=0.0)?fabs((double)gpu_sum-host_sum)/fabs(host_sum):0.0;
    printf("# [Reduktion] GPU=%.6e Host=%.6e rel=%.2e -> %s\n",
           (double)gpu_sum,host_sum,rel,(rel<SELFCHK_TOL)?"OK":"ABWEICHUNG!");

    /* =================== Headline-Zeitschleife (NREP Wiederholungen, Median) === */
    double loop_s[64]; if(NREP>64)NREP=64;
    double drift_last=0.0;
#ifdef USE_GRAPH
    /* --- CUDA-Graph: die feste 5-Kernel-Schrittsequenz EINMAL aufnehmen -----
     *  Statt pro Schritt 5 Kernel einzeln vom Host aus zu starten (je mit
     *  CPU-seitigem Launch-Overhead), wird die Sequenz einmal per
     *  Stream-Capture aufgezeichnet (cudaStreamBeginCapture/-EndCapture)
     *  und danach als Graph-Objekt instanziiert (cudaGraphInstantiate).
     *  Pro Schritt genuegt dann EIN cudaGraphLaunch() statt fuenf einzelner
     *  Kernel-Starts -- der Treiber spielt die aufgezeichnete Sequenz ab,
     *  ohne dass der Host jeden Kernel-Start einzeln anstossen muss.
     *  Physik/Kernel-Reihenfolge ist UNVERAENDERT -> der Drift-Wert unten
     *  muss (bis auf Rundung) identisch zur Nicht-Graph-Variante sein; das
     *  ist der eingebaute Korrektheits-Check fuer diese Optimierung. */
    cudaStream_t gs; CUDA_CHECK(cudaStreamCreate(&gs));
    cudaGraph_t     cgraph;
    cudaGraphExec_t gexec;
    CUDA_CHECK(cudaStreamBeginCapture(gs, cudaStreamCaptureModeGlobal));
    k_kick_drift<<<B,T,0,gs>>>(d_px,d_py,d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
    k_clear_cells<<<Bc,T,0,gs>>>(d_head,ncell);
    k_build_cells<<<B,T,0,gs>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
    k_force<<<B,T,0,gs>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                          x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N);
    k_kick<<<B,T,0,gs>>>(d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
    CUDA_CHECK(cudaStreamEndCapture(gs,&cgraph));
    CUDA_CHECK(cudaGraphInstantiate(&gexec,cgraph,0));
    printf("# [CUDA-Graph] Schrittsequenz aufgenommen (5 Kernels/Schritt)\n");
#endif
    for (int rep=0; rep<NREP; rep++){
        /* --- Zustand auf Anfangswert zuruecksetzen (Device->Device Copy) --- */
        CUDA_CHECK(cudaMemcpy(d_px,d_px0,db,cudaMemcpyDeviceToDevice));
        CUDA_CHECK(cudaMemcpy(d_py,d_py0,db,cudaMemcpyDeviceToDevice));
        CUDA_CHECK(cudaMemcpy(d_vx,d_vx0,db,cudaMemcpyDeviceToDevice));
        CUDA_CHECK(cudaMemcpy(d_vy,d_vy0,db,cudaMemcpyDeviceToDevice));
        CUDA_CHECK(cudaMemcpy(d_int,d_int0,ib,cudaMemcpyDeviceToDevice));
        CUDA_CHECK(cudaMemcpy(d_eb,&zero,sizeof(real),cudaMemcpyHostToDevice));

        /* Anfangs-Cell-List und -Kraefte berechnen (ausserhalb des Graphen,
         * damit Graph und Nicht-Graph-Pfad mit identischem Startzustand in
         * die getaktete Schleife gehen) */
        k_clear_cells<<<Bc,T>>>(d_head,ncell);
        k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
        k_force<<<B,T>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                         x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N);
        CUDA_CHECK(cudaMemcpy(d_eb,&zero,sizeof(real),cudaMemcpyHostToDevice));

        /* Startenergie fuer Drift-Check dieser Wdh. */
        k_energy<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_m,d_nbr,d_int,d_nbr0,d_head,d_next,
                          x0,y0,cs,ncx,ncy,K,L0,RCUT,KREP,REPN,d_ea,N);
        real es0=gpu_reduce_sum(d_ea,d_esum,N,T);
        real eb0; CUDA_CHECK(cudaMemcpy(&eb0,d_eb,sizeof(real),cudaMemcpyDeviceToHost));
        double e0=(double)es0+(double)eb0;

        /* Aufwaermphase (ungetaktet, immer als Einzel-Kernel-Starts -- der
         * Graph wird erst in der GETAKTETEN Schleife unten verwendet) */
        for (int s=0;s<WARMUP;s++){
            k_kick_drift<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
            k_clear_cells<<<Bc,T>>>(d_head,ncell);
            k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
            k_force<<<B,T>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                             x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N);
            k_kick<<<B,T>>>(d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        /* Getaktete Hauptschleife dieser Wiederholung: entweder NSTEPS-mal
         * cudaGraphLaunch (USE_GRAPH) oder NSTEPS-mal die 5 Kernel einzeln
         * starten (Baseline) -- Physik ist in beiden Faellen identisch. */
        cudaEvent_t l0,l1; CUDA_CHECK(cudaEventCreate(&l0)); CUDA_CHECK(cudaEventCreate(&l1));
#ifdef USE_GRAPH
        CUDA_CHECK(cudaEventRecord(l0,gs));
        for (int s=0;s<NSTEPS;s++) CUDA_CHECK(cudaGraphLaunch(gexec,gs));
        CUDA_CHECK(cudaEventRecord(l1,gs)); CUDA_CHECK(cudaEventSynchronize(l1));
#else
        CUDA_CHECK(cudaEventRecord(l0));
        for (int s=0;s<NSTEPS;s++){
            k_kick_drift<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
            k_clear_cells<<<Bc,T>>>(d_head,ncell);
            k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
            k_force<<<B,T>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                             x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N);
            k_kick<<<B,T>>>(d_vx,d_vy,d_fx,d_fy,d_m,dt,N);
        }
        CUDA_CHECK(cudaEventRecord(l1)); CUDA_CHECK(cudaEventSynchronize(l1));
#endif
        float ms=0; CUDA_CHECK(cudaEventElapsedTime(&ms,l0,l1));
        loop_s[rep]=ms/1000.0;
        cudaEventDestroy(l0); cudaEventDestroy(l1);

        /* Drift am Ende dieser Wiederholung (Korrektheits-Check) */
        k_energy<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_m,d_nbr,d_int,d_nbr0,d_head,d_next,
                          x0,y0,cs,ncx,ncy,K,L0,RCUT,KREP,REPN,d_ea,N);
        real esE=gpu_reduce_sum(d_ea,d_esum,N,T);
        real ebE; CUDA_CHECK(cudaMemcpy(&ebE,d_eb,sizeof(real),cudaMemcpyDeviceToHost));
        drift_last=100.0*(((double)esE+(double)ebE)-e0)/e0;
        printf("#   rep %d: %.3f s  (drift %.4f %%)\n", rep, loop_s[rep], drift_last);
    }
    /* Median + Min/Max aus allen Wiederholungen (Median: robust gegen Ausreisser) */
    double lo=loop_s[0],hi=loop_s[0]; for(int i=1;i<NREP;i++){ if(loop_s[i]<lo)lo=loop_s[i]; if(loop_s[i]>hi)hi=loop_s[i]; }
    double t_loop_s=median(loop_s,NREP);

    /* =================== Kernel-Breakdown ==================================== */
    /* Jeden Kernel isoliert FORCE_REP mal messen (identisch zu cascade_cuda_bench.cu).
     * Reihenfolge: erst lesende Kernel, dann die positions-modifizierenden zuletzt,
     * damit die Positionen fuer alle Messungen konsistent bleiben. Laeuft IMMER
     * als Einzel-Kernel-Starts, auch bei USE_GRAPH -- der Breakdown misst die
     * Kernel-Kosten selbst, nicht den Launch-Mechanismus. */
    cudaEvent_t b0,b1; CUDA_CHECK(cudaEventCreate(&b0)); CUDA_CHECK(cudaEventCreate(&b1));
    bk_n=0;
    /* BENCH-Makro: NAME = Kernel-Name, REP = Wdh., STMT = Kernel-Aufruf.
     * cudaEventSynchronize(b1) stellt sicher, dass alle GPU-Arbeit fertig ist,
     * bevor die Endzeit gemessen wird; die Zeit wird durch REP dividiert, um
     * die Start-Latenz eines einzelnen Aufrufs (~5 us) wegzumitteln. */
    #define BENCH(NAME,REP,STMT) do{ CUDA_CHECK(cudaEventRecord(b0)); \
        for(int _r=0;_r<(REP);_r++){ STMT; } CUDA_CHECK(cudaEventRecord(b1)); \
        CUDA_CHECK(cudaEventSynchronize(b1)); float _m=0; CUDA_CHECK(cudaEventElapsedTime(&_m,b0,b1)); \
        bk_name[bk_n]=NAME; bk_us[bk_n]=_m*1000.0/(REP); bk_n++; }while(0)

    /* Cell List frisch aufbauen (Zustand fuer alle Breakdown-Messungen) */
    k_clear_cells<<<Bc,T>>>(d_head,ncell);
    k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
    CUDA_CHECK(cudaDeviceSynchronize());
    BENCH("clear_cells",FORCE_REP, (k_clear_cells<<<Bc,T>>>(d_head,ncell)));
    BENCH("build_cells",FORCE_REP, (k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N)));
    /* Nach build_cells die Cell List neu aufbauen (build hat sie oben veraendert) */
    k_clear_cells<<<Bc,T>>>(d_head,ncell);
    k_build_cells<<<B,T>>>(d_px,d_py,d_head,d_next,x0,y0,cs,ncx,ncy,N);
    BENCH("force",FORCE_REP, (k_force<<<B,T>>>(d_px,d_py,d_fx,d_fy,d_nbr,d_int,d_nbr0,d_head,d_next,
                              x0,y0,cs,ncx,ncy,K,L0,MAXSTR,RCUT,KREP,REPN,d_eb,N)));
    BENCH("energy",FORCE_REP, (k_energy<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_m,d_nbr,d_int,d_nbr0,d_head,d_next,
                              x0,y0,cs,ncx,ncy,K,L0,RCUT,KREP,REPN,d_ea,N)));
    BENCH("kick_drift",FORCE_REP, (k_kick_drift<<<B,T>>>(d_px,d_py,d_vx,d_vy,d_fx,d_fy,d_m,dt,N)));
    BENCH("kick",FORCE_REP, (k_kick<<<B,T>>>(d_vx,d_vy,d_fx,d_fy,d_m,dt,N)));

    /* Force-Kernel-Zeit aus dem Breakdown-Array heraussuchen (fuer Bandbreitenmodell) */
    double t_force_us=0; for(int i=0;i<bk_n;i++) if(!strcmp(bk_name[i],"force")) t_force_us=bk_us[i];

    /* Pro-Schritt-Summe: alle Kernels AUSSER energy (das laeuft nur beim Logging) */
    double step_us=0; for(int i=0;i<bk_n;i++) if(strcmp(bk_name[i],"energy")) step_us+=bk_us[i];

    /* ---- D2H-Transfer messen ------------------------------------------------ */
    CUDA_CHECK(cudaEventRecord(ev0));
    CUDA_CHECK(cudaMemcpy(h_ea,d_ea,db,cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaEventRecord(ev1)); CUDA_CHECK(cudaEventSynchronize(ev1));
    float t_d2h_ms=0; CUDA_CHECK(cudaEventElapsedTime(&t_d2h_ms,ev0,ev1));

    /* =================== Ergebnisse ausgeben ================================= */
    double per_step_us=(t_loop_s*1e6)/NSTEPS;
    double matom_s=((double)N*NSTEPS)/t_loop_s/1e6;
    double bw_model=(BYTES_PER_ATOM*(double)N)/(t_force_us/1e6)/1e9;
    double peak_pct=100.0*bw_model/PEAK_BW_GBS;

    printf("\n# ================== ERGEBNIS-BLOCK (%s, N=%d) ==================\n",PREC_NAME,N);
    printf("# Reduktions-Drift (letzte Wdh.) : %.4f %%\n", drift_last);
    printf("# Zeitschleife (Median v. %d)    : %.4f s   [min %.4f / max %.4f]\n", NREP, t_loop_s, lo, hi);
    printf("# pro Schritt                    : %.2f us\n", per_step_us);
    printf("# Durchsatz                      : %.1f Matom/s\n", matom_s);
    printf("# --- Kernel-Breakdown (isoliert, us/Aufruf) ---\n");
    for(int i=0;i<bk_n;i++){
        int per_step = strcmp(bk_name[i],"energy")!=0;
        /* (nur beim Logging): energy-Kernel laeuft nicht jeden Schritt */
        printf("#   %-12s %8.2f us  %s\n", bk_name[i], bk_us[i],
               per_step ? "" : "(nur beim Logging)");
    }
    printf("#   -> pro Schritt (Summe ohne energy): %.2f us ; force-Anteil = %.1f %%\n",
           step_us, 100.0*t_force_us/step_us);
    printf("# --- Transfer ---\n");
    printf("# H2D (%.0f MB) : %.2f ms (%.0f GB/s) ; D2H (%.0f MB) : %.2f ms (%.0f GB/s)\n",
           h2d_bytes/1e6, t_h2d_ms, (h2d_bytes/1e9)/(t_h2d_ms/1000.0),
           db/1e6, t_d2h_ms, (db/1e9)/(t_d2h_ms/1000.0));
    printf("# Bandbreite force (MODELL)      : %.1f GB/s (%.1f %% Peak) -- ncu ist autoritativ\n",
           bw_model, peak_pct);
    /* CPU-Vergleich nur sinnvoll fuer fp64 bei ~2 Mio. Atomen (Baseline-Konfig) */
#ifndef USE_FP32
    if (labs((long)N - 1999396L) < 50000L){
        printf("# --- CPU vs GPU (fp64, gleiche Konfig) ---\n");
        printf("# vs 1 EPYC-Kern (1152,9 s) : %.1fx ; vs 64-Kern-Knoten (%.1f s) : %.2fx\n",
               T_CPU_1CORE/t_loop_s, T_CPU_1CORE/SPEEDUP_64, (T_CPU_1CORE/SPEEDUP_64)/t_loop_s);
    }
#endif
    /* Maschinenlesbare CSV-Zeile (z.B. fuer manuelle Aggregation ueber die 4 Varianten) */
    printf("# CSV  prec,N,T,steps,per_step_us,matom_s,force_us,step_us\n");
    printf("# CSV  %s,%d,%d,%d,%.3f,%.2f,%.3f,%.3f\n",
           PREC_NAME,N,T,NSTEPS,per_step_us,matom_s,t_force_us,step_us);
    printf("# =====================================================================\n");

    /* ---- Aufraeumen -------------------------------------------------------- */
    cudaFree(d_px);cudaFree(d_py);cudaFree(d_vx);cudaFree(d_vy);cudaFree(d_fx);cudaFree(d_fy);
    cudaFree(d_m);cudaFree(d_ea);cudaFree(d_eb);cudaFree(d_esum);
    cudaFree(d_nbr);cudaFree(d_int);cudaFree(d_nbr0);cudaFree(d_head);cudaFree(d_next);
    cudaFree(d_px0);cudaFree(d_py0);cudaFree(d_vx0);cudaFree(d_vy0);cudaFree(d_int0);
    free(h_px);free(h_py);free(h_vx);free(h_vy);free(h_m);free(h_ea);
    free(h_nbr);free(h_int);free(h_cnt);free(h_nbr0);
    return 0;
}
