/* ============================================================================
 *  mpi_cascade.c  --  Parallele Kaskaden-Simulation mit MPI (Modell A)
 * ----------------------------------------------------------------------------
 *
 *  ÜBERBLICK: WAS MACHT DIESE VERSION?
 *  ------------------------------------
 *  Diese Datei ist die MPI-parallele Version von cascade_serial.c.
 *  Die Physik ist IDENTISCH: 2D-Dreiecksgitter, harmonische Federn,
 *  Potenzgesetz-Abstossung, Velocity-Verlet-Integrator.
 *
 *  Was sich aendert: Das Gitter wird unter mehreren MPI-Prozessen ("Ranks")
 *  aufgeteilt. Jeder Rank rechnet nur einen STREIFEN des Gitters und tauscht
 *  mit seinen Nachbarn die noetigsten Daten aus.
 *
 *  KONSISTENZTEST: Die Ergebnisse mit 1, 2, 4, 8 Ranks muessen BITGLEICH sein.
 *
 *  WARUM MPI SCHNELLER IST:
 *  -------------------------
 *  Mit P Ranks hat jeder Rank nur N/P Atome zu rechnen.
 *  Die Kraftberechnung skaliert ideal mit P (sie haengt nur von lokalen Atomen
 *  und deren Nachbarn ab). Die Kommunikation (Halo-Austausch) ist nur proportional
 *  zum Rand des Streifens -- das ist O(N/P)-mal kleiner als die Arbeit O(N/P).
 *
 *  1D-DOMÄNENZERLEGUNG (y-Richtung):
 *  -----------------------------------
 *  Das Gitter wird in horizontale STREIFEN geschnitten:
 *
 *    Rank 0: Zeilen y in [0,    Ly/4)   <- "unten"
 *    Rank 1: Zeilen y in [Ly/4, Ly/2)
 *    Rank 2: Zeilen y in [Ly/2, 3Ly/4)
 *    Rank 3: Zeilen y in [3Ly/4, Ly)    <- "oben"
 *
 *  Jeder Rank kennt NUR die Atome in seinem Streifen plus einen schmalen
 *  "Halo"-Rand (Kopien der Randatome des Nachbarn, nur fuer Kraftberechnung).
 *
 *  KOMMUNIKATION PRO ZEITSCHRITT:
 *  --------------------------------
 *  (A) MIGRATION: Atome, die durch hohe Geschwindigkeit die Streifengrenze
 *      ueberquert haben, werden zu ihrem neuen Besitzer uebertragen.
 *      Vollstaendiger Atomzustand (19 doubles/Atom) wird verschickt.
 *
 *  (B) HALO-AUSTAUSCH: Kopien der Randatome (nur Position + Geschwindigkeit,
 *      5 doubles/Atom) werden an den Nachbarrank geschickt.
 *      Diese "Geister" (Halo-Atome) erlauben die Berechnung von Kraften
 *      ueber die Streifengrenze hinweg, ohne weiteren Datenaustausch.
 *
 *  ZWEI KOMMUNIKATIONSVARIANTEN:
 *  --------------------------------
 *  Blocking (Default):     MPI_Sendrecv -- simpel, aber seriell: erst Up, dann Down.
 *  Non-blocking (HALO_NB): MPI_Isend/Irecv -- Up und Down gleichzeitig.
 *
 *  HASH-TABELLE fuer O(1)-Nachbarsuche:
 *  ----------------------------------------
 *  Problem: In der MPI-Version kennt jeder Rank nicht die lokalen Array-Indizes
 *  seiner Nachbarn. Die Bindungsliste speichert globale Atom-IDs (gid).
 *  Um schnell von gid -> lokaler Index zu kommen, gibt es pro Strip eine
 *  Hash-Tabelle (Open-Addressing, Fibonacci-Hashing).
 *
 *  Build:  mpicc -O3 -march=znver3 -o mpi_cascade mpi_cascade.c -lm
 *  Run:    srun -n <ranks> ./mpi_cascade params.ini
 * ==========================================================================*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 *  GLOBALE PARAMETER  (von Rank 0 gelesen, dann an alle gebroadcastet)
 * ========================================================================== */
static int    NX=120, NY=120;          /* Gitterdimensionen (Gesamtgitter)      */
static double L0=1.0, K_SPRING=100.0, MAX_STRETCH=1.15;
static double RCUT=0.9, K_REP=400.0, REP_N=12.0;
static int    HEALING=0; static double HEALING_DIST=1.10, HEALING_VREL=0.5;
static double PKA_X=0.5,PKA_Y=0.5,PKA_ENERGY=1000.0,PKA_ANGLE=15.0,PKA_MASS=1.0;
static int    N_PKA=1;                 /* Anzahl gleichzeitiger PKAs            */
#define MAX_PKA 64
typedef struct { double x,y,energy,angle,mass; int hx,hy,he,ha,hm; } PkaDef;
static PkaDef MANUAL_PKA[MAX_PKA];    /* aus [pka1],[pka2],... gefuellt        */
static int    N_MANUAL=0;             /* hoechster definierter [pkaN]-Index    */
static int    TIMING=1;               /* phasenaufgeloeste Zeitmessung an/aus  */
static int    HALO_NB=0;              /* 0=blocking Sendrecv, 1=non-blocking   */
static int    RESORT_EVERY=0;         /* 0=nie, >0=alle N Schritte cell-sortieren */
static double DT=5e-4; static int NSTEPS=14000, LOG_EVERY=0; static unsigned SEED=1;
static char   OUT_PREFIX[128]="run";
static double DY;                     /* L0*sqrt(3)/2 (Zeilenabstand, abgeleitet) */
#define MAXNBR 6                       /* max. Nachbarn im Dreiecksgitter       */

/* ============================================================================
 *  DETERMINISTISCHER ZUFALLSZAHLENGENERATOR
 * ---------------------------------------------------------------------------
 *  Auf ALLEN Ranks identisch initialisiert (gleicher SEED) -> alle Ranks
 *  berechnen dieselbe PKA-Liste und setzen jeweils nur die PKAs, deren
 *  Atom-gid in den eigenen Streifen faellt.  So ergibt 1 Rank + 4 Ranks
 *  exakt dieselbe Anfangsbedingung.
 * ========================================================================== */
static unsigned long RNG;
static double rng_u(void){
    RNG = RNG*6364136223846793005UL + 1442695040888963407UL;
    return ((RNG>>11) & ((1UL<<52)-1)) / (double)(1UL<<52);
}

/* ============================================================================
 *  PHASEN-TIMER
 * ---------------------------------------------------------------------------
 *  Jede Berechnungsphase wird separat gemessen.  Am Ende gibt Rank 0
 *  MAX (Flaschenhals) und AVG (Mittelwert) aller Ranks aus.
 *  MAX >> AVG in einer Phase bedeutet Lastungleichgewicht.
 * ========================================================================== */
static double T_force=0,T_cells=0,T_halo=0,T_migrate=0,T_energy=0,T_hash=0;

/* ============================================================================
 *  CONFIG-PARSER  (wie in der seriellen Version; nur Rank 0 liest die Datei)
 * ========================================================================== */
static void read_config(const char*fn){
    FILE*f=fopen(fn,"r"); if(!f){fprintf(stderr,"WARN: %s fehlt\n",fn);return;}
    char line[256];
    char section[64]="";
    while(fgets(line,sizeof line,f)){
        char*h=strpbrk(line,";#"); if(h)*h=0;
        char*lb=strchr(line,'[');
        if(lb){ section[0]=0; sscanf(lb,"[%63[^]]]",section); continue; }
        char*eq=strchr(line,'='); if(!eq)continue; *eq=0;
        char k[64],v[128];
        if(sscanf(line,"%63s",k)!=1)continue; if(sscanf(eq+1,"%127s",v)!=1)continue;
        #define KV(s) (!strcmp(k,s))
        /* [pka1],[pka2],...-Sektionen in das Manual-Array routen */
        { int pidx;
          if(sscanf(section,"pka%d",&pidx)==1 && pidx>=1 && pidx<=MAX_PKA){
            PkaDef*p=&MANUAL_PKA[pidx-1]; if(pidx>N_MANUAL)N_MANUAL=pidx;
            if      (KV("pka_x"))     {p->x=atof(v);     p->hx=1;}
            else if (KV("pka_y"))     {p->y=atof(v);     p->hy=1;}
            else if (KV("pka_energy")){p->energy=atof(v);p->he=1;}
            else if (KV("pka_angle")) {p->angle=atof(v); p->ha=1;}
            else if (KV("pka_mass"))  {p->mass=atof(v);  p->hm=1;}
            continue;
          } }
        if(KV("healing"))HEALING=(!strcmp(v,"true")||!strcmp(v,"1"));
        else if(KV("healing_dist"))HEALING_DIST=atof(v);
        else if(KV("healing_vrel"))HEALING_VREL=atof(v);
        else if(KV("NX"))NX=atoi(v); else if(KV("NY"))NY=atoi(v);
        else if(KV("L0"))L0=atof(v); else if(KV("K_SPRING"))K_SPRING=atof(v);
        else if(KV("MAX_STRETCH"))MAX_STRETCH=atof(v); else if(KV("RCUT"))RCUT=atof(v);
        else if(KV("K_REP"))K_REP=atof(v); else if(KV("REP_N"))REP_N=atof(v);
        else if(KV("pka_x"))PKA_X=atof(v); else if(KV("pka_y"))PKA_Y=atof(v);
        else if(KV("pka_energy"))PKA_ENERGY=atof(v); else if(KV("pka_angle"))PKA_ANGLE=atof(v);
        else if(KV("pka_mass"))PKA_MASS=atof(v);
        else if(KV("n_pka"))N_PKA=atoi(v);
        else if(KV("timing"))TIMING=(!strcmp(v,"true")||!strcmp(v,"1"));
        else if(KV("halo_nonblocking"))HALO_NB=(!strcmp(v,"true")||!strcmp(v,"1"));
        else if(KV("resort_every"))RESORT_EVERY=atoi(v);
        else if(KV("dt"))DT=atof(v); else if(KV("n_steps"))NSTEPS=atoi(v);
        else if(KV("log_every"))LOG_EVERY=atoi(v); else if(KV("seed"))SEED=(unsigned)atoi(v);
        else if(KV("out_prefix"))strncpy(OUT_PREFIX,v,127);
        #undef KV
    }
    fclose(f);
}

/* ============================================================================
 *  BROADCAST DER KONFIGURATION
 * ---------------------------------------------------------------------------
 *  Rank 0 hat die Konfiguration gelesen.  Alle anderen Ranks kennen nur die
 *  Default-Werte.  MPI_Bcast kopiert jeden Parameter von Rank 0 zu allen.
 *
 *  Warum nicht einfach alle Ranks lesen?  Das wuerde funktionieren, skaliert
 *  aber schlecht: Bei N_RANKS Prozessen auf N_NODES Knoten und einem
 *  gemeinsamen Dateisystem (wie Lustre/NFS) wuerden N_RANKS gleichzeitige
 *  Lesevorgaenge das Dateisystem ueberlasten.  Bcast ist die richtige Loesung.
 * ========================================================================== */
static void bcast_config(void){
    MPI_Bcast(&NX,1,MPI_INT,0,MPI_COMM_WORLD); MPI_Bcast(&NY,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&L0,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&K_SPRING,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&MAX_STRETCH,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&RCUT,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&K_REP,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&REP_N,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&HEALING,1,MPI_INT,0,MPI_COMM_WORLD); MPI_Bcast(&HEALING_DIST,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&HEALING_VREL,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&PKA_X,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&PKA_Y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&PKA_ENERGY,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&PKA_ANGLE,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&PKA_MASS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&N_PKA,1,MPI_INT,0,MPI_COMM_WORLD); MPI_Bcast(&TIMING,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&N_MANUAL,1,MPI_INT,0,MPI_COMM_WORLD);
    /* MANUAL_PKA als rohes Byte-Array senden (alle Felder auf einmal) */
    MPI_Bcast(MANUAL_PKA,(int)sizeof(MANUAL_PKA),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&HALO_NB,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&RESORT_EVERY,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&DT,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&NSTEPS,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&LOG_EVERY,1,MPI_INT,0,MPI_COMM_WORLD); MPI_Bcast(&SEED,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
}

/* ============================================================================
 *  GEOMETRIE-FUNKTIONEN  (arbeiten auf globaler Atom-ID)
 * ---------------------------------------------------------------------------
 *  Eine globale Atom-ID (gid) kodiert Position im Gesamtgitter:
 *    gid = j * NX + i     (Zeile j, Spalte i)
 *
 *  gx_of(gid): x-Koordinate aus gid berechnen
 *    Gerade Zeilen:   x = i * L0
 *    Ungerade Zeilen: x = i * L0 + L0/2  (Versatz nach rechts)
 *    (j&1) ist genau 1 fuer ungerade j, 0 fuer gerade j.
 *
 *  gy_of(gid): y-Koordinate aus gid berechnen
 *    Alle Zeilen: y = j * DY,  mit DY = L0*sqrt(3)/2
 *
 *  neighbors_of(gid, out[6]): liefert die (bis zu 6) Nachbar-gids
 *    Im Dreiecksgitter hat ein Innenatom genau 6 Nachbarn:
 *    links, rechts, plus 4 diagonale je nach gerader/ungerader Zeile.
 *    Randatome haben weniger (Grenzabfragen: i>=0, i<NX, j>=0, j<NY).
 * ========================================================================== */
static inline double gx_of(long gid){ long j=gid/NX,i=gid%NX; return i*L0+(j&1)*(0.5*L0); }
static inline double gy_of(long gid){ long j=gid/NX; return j*DY; }
static int neighbors_of(long gid,long out[6]){
    long j=gid/NX,i=gid%NX; int n=0;
    if(i+1<NX) out[n++]=j*NX+(i+1);    /* rechts */
    if(i-1>=0) out[n++]=j*NX+(i-1);    /* links  */
    long ul,ur,dl,dr;
    /* Versatz der Diagonalen haengt von gerader/ungerader Zeile ab */
    if((j&1)==0){ ul=i-1; ur=i; dl=i-1; dr=i; }
    else         { ul=i;   ur=i+1; dl=i; dr=i+1; }
    if(j+1<NY){ if(ul>=0&&ul<NX)out[n++]=(j+1)*NX+ul; /* oben-links  */
                if(ur>=0&&ur<NX)out[n++]=(j+1)*NX+ur;} /* oben-rechts */
    if(j-1>=0){ if(dl>=0&&dl<NX)out[n++]=(j-1)*NX+dl; /* unten-links */
                if(dr>=0&&dr<NX)out[n++]=(j-1)*NX+dr;} /* unten-rechts*/
    return n;
}

/* ============================================================================
 *  STRIP-STRUKTUR  --  zentraler Datenspeicher jedes Ranks
 * ---------------------------------------------------------------------------
 *  Ein Strip haelt alle Daten fuer den lokalen Streifen:
 *
 *  RESIDENTE Atome (Indizes 0..n-1):
 *    Atome, die wirklich zu diesem Rank gehoeren.  Ihre Kraefte werden hier
 *    berechnet, ihre Positionen hier aktualisiert.
 *
 *  HALO-Atome (Indizes n..n+nh-1):
 *    Kopien von Randatomen der Nachbar-Ranks.  Nur Position und Geschwindigkeit
 *    sind gueltig (mass=1, nb=0).  Werden nur fuer die Kraftberechnung genutzt
 *    und jeden Schritt neu eingelesen.  Kraefte werden NICHT zurueckgeschickt.
 *
 *  SoA (Structure of Arrays): px[],py[],... -- gleich wie in der seriellen Version.
 *
 *  gid[p]:       globale Atom-ID von Atom p (Schluessel in der Hash-Tabelle)
 *  nb[p]:        Anzahl gespeicherter Bindungspartner von Atom p
 *  pg[p*MAXNBR+e]: globale ID des e-ten Bindungspartners
 *  intact[p*MAXNBR+e]: 1=Bindung intakt, 0=gerissen
 *  E_broken:     akkumulierte Bruch-Energie dieses Strips
 *
 *  HASH-TABELLE (hkey/hval):
 *    hkey[h] = gid eines Atoms (oder -1 = leer)
 *    hval[h] = lokaler Index dieses Atoms
 *    Groesse hcap: naechste Zweierpotenz >= 4*cap (niedrige Auslastung -> wenig Kollisionen)
 * ========================================================================== */
typedef struct {
    double ylo,yhi;           /* y-Grenzen des Streifens (exklusiv bei yhi)   */
    int n,nh,cap;             /* n=Residente, nh=Halo-Atome, cap=Array-Kapazitaet */
    long  *gid;               /* globale Atom-IDs [0..n+nh)                   */
    double *px,*py,*vx,*vy,*fx,*fy,*mass;  /* Zustand                         */
    int   *nb;                /* Anzahl Bindungspartner je Atom               */
    long  *pg;                /* Bindungspartner: pg[p*MAXNBR+e] = gid des e-ten */
    int   *intact;            /* intact[p*MAXNBR+e] = 1/0                    */
    double E_broken;          /* kumulierte Bruchenergie dieses Strips        */
    int    hcap;              /* Hash-Tabellen-Kapazitaet (Zweierpotenz)      */
    long  *hkey;              /* Hash-Tabellen-Schluessel (gid oder -1=leer)  */
    int   *hval;              /* Hash-Tabellen-Werte (lokaler Index)          */
} Strip;

/* Allokiert alle Arrays des Strips.  cap = maximale Atomzahl (inkl. Headroom) */
static void strip_alloc(Strip*s,int cap){
    s->cap=cap; s->n=0; s->nh=0; s->E_broken=0;
    s->gid=malloc(cap*sizeof*s->gid);
    s->px=malloc(cap*sizeof*s->px); s->py=malloc(cap*sizeof*s->py);
    s->vx=malloc(cap*sizeof*s->vx); s->vy=malloc(cap*sizeof*s->vy);
    s->fx=malloc(cap*sizeof*s->fx); s->fy=malloc(cap*sizeof*s->fy);
    s->mass=malloc(cap*sizeof*s->mass);
    s->nb=malloc(cap*sizeof*s->nb);
    /* pg und intact: MAXNBR Eintraege pro Atom */
    s->pg=malloc((size_t)cap*MAXNBR*sizeof*s->pg);
    s->intact=malloc((size_t)cap*MAXNBR*sizeof*s->intact);
    /* Hash-Tabelle: groesste Zweierpotenz >= 4*cap (Load-Factor <= 0.25) */
    s->hcap=1; while(s->hcap<4*cap) s->hcap<<=1;
    s->hkey=malloc(s->hcap*sizeof*s->hkey); s->hval=malloc(s->hcap*sizeof*s->hval);
}

/* ============================================================================
 *  HASH-TABELLE  --  globale ID -> lokaler Index
 * ---------------------------------------------------------------------------
 *  Problem: Die Bindungsliste speichert globale Atom-IDs (pg[]).  Um die
 *  Kraft zwischen Atom p und seinem Bindungspartner zu berechnen, muss man
 *  den lokalen Array-Index des Partners kennen.
 *
 *  Loesung: Hash-Tabelle mit Open-Addressing (linear probing).
 *  Fibonacci-Hashing: Multiplikation mit 2654435761 (goldene Zahl * 2^32)
 *  streut die IDs gut ueber den Tabellen-Indexraum.
 *
 *  hash_build(): Tabellengroesse hcap (Zweierpotenz), alle Eintraege auf -1
 *  (leer) setzen, dann alle n+nh Atome eintragen.  Wird jeden Schritt
 *  nach Migration + Halo-Austausch aufgerufen.
 *
 *  hash_get(): Sucht gid in der Tabelle.  Gibt lokalen Index zurueck
 *  oder -1 wenn nicht gefunden.  -1 bedeutet: Partner ist kein lokales
 *  oder Halo-Atom -> Bindungskraft wird von diesem Rank nicht berechnet.
 *  (Das ist korrekt: der andere Rank berechnet sie.)
 * ========================================================================== */
static void hash_build(Strip*s){
    int m=s->hcap;
    for(int t=0;t<m;t++) s->hkey[t]=-1;   /* alle Eintraege leeren */
    int tot=s->n+s->nh;                     /* Residente + Halo-Atome eintragen */
    for(int p=0;p<tot;p++){
        long g=s->gid[p];
        /* Fibonacci-Hashing: ergibt gleichmaessige Streuung */
        unsigned long h=(unsigned long)g*2654435761u & (m-1);
        /* Linear Probing: naechsten freien Platz suchen */
        while(s->hkey[h]!=-1) h=(h+1)&(m-1);
        s->hkey[h]=g; s->hval[h]=p;
    }
}

/* Sucht gid in der Hash-Tabelle und gibt den lokalen Index zurueck (-1 = nicht gefunden) */
static inline int hash_get(Strip*s,long g){
    int m=s->hcap;
    unsigned long h=(unsigned long)g*2654435761u & (m-1);
    while(s->hkey[h]!=-1){
        if(s->hkey[h]==g) return s->hval[h];
        h=(h+1)&(m-1);
    }
    return -1;
}

/* ============================================================================
 *  ABSTOSS-KRAFT und -POTENTIAL  (identisch zur seriellen Version)
 * ========================================================================== */
static inline double rep_force(double r){ return K_REP*(pow(RCUT/r,REP_N)-1.0); }
static inline double rep_pot(double r){ double n=REP_N;
    return K_REP*((RCUT-pow(RCUT,n)*pow(r,1.0-n))/(1.0-n)-(RCUT-r)); }

/* ============================================================================
 *  CELL LIST  --  lokale Variante (pro Strip)
 * ---------------------------------------------------------------------------
 *  Identisch zur seriellen Version, aber:
 *   - Wird in der Cells-Struct gespeichert (kein globaler Zustand).
 *   - Umfasst ALLE Atome des Strips: Residente + Halo-Atome.
 *     Halo-Atome werden benoetigt, weil die Abstossung ueber die Streifengrenze
 *     wirken kann (wenn r < RCUT < L0 ist die Abstossreichweite kleiner als
 *     der Bindungsabstand, also wirkt sie nur bei sehr nah kommenden Paaren).
 *   - Rand: 2*L0 (grosszuegiger als seriell, weil Halo-Atome etwas ausserhalb
 *     des Streifens liegen koennen).
 * ========================================================================== */
typedef struct { int ncx,ncy; double x0,y0,cs; int*head,*next; } Cells;
static void cells_build(Strip*s,Cells*c){
    int tot=s->n+s->nh;
    double x0=1e30,y0=1e30,x1=-1e30,y1=-1e30;
    for(int p=0;p<tot;p++){
        if(s->px[p]<x0)x0=s->px[p]; if(s->px[p]>x1)x1=s->px[p];
        if(s->py[p]<y0)y0=s->py[p]; if(s->py[p]>y1)y1=s->py[p];
    }
    double m=2.0*L0; c->x0=x0-m; c->y0=y0-m; c->cs=RCUT;
    c->ncx=(int)(((x1-x0)+2*m)/c->cs)+1; c->ncy=(int)(((y1-y0)+2*m)/c->cs)+1;
    c->head=malloc((size_t)c->ncx*c->ncy*sizeof(int));
    c->next=malloc(tot*sizeof(int));
    for(int k=0;k<c->ncx*c->ncy;k++)c->head[k]=-1;
    for(int p=0;p<tot;p++){
        int cx=(int)((s->px[p]-c->x0)/c->cs), cy=(int)((s->py[p]-c->y0)/c->cs);
        if(cx<0)cx=0; if(cx>=c->ncx)cx=c->ncx-1;
        if(cy<0)cy=0; if(cy>=c->ncy)cy=c->ncy-1;
        int cc=cy*c->ncx+cx; c->next[p]=c->head[cc]; c->head[cc]=p;
    }
}
static void cells_free(Cells*c){ free(c->head); free(c->next); }

/* Prueft ob Atom p und das Atom mit gid qgid durch eine INTAKTE Bindung verbunden sind.
 * Benoetigt fuer: Abstoss-Ausschluss gebundener Paare (Feder-Abstoss-Doppelzaehlung vermeiden).
 * Lineare Suche in pg[]-Liste (max. MAXNBR=6 Eintraege). */
static inline int bonded_to(Strip*s,int p,long qgid){
    int base=p*MAXNBR, nb=s->nb[p];
    for(int e=0;e<nb;e++) if(s->intact[base+e] && s->pg[base+e]==qgid) return 1;
    return 0;
}

/* ============================================================================
 *  HEALING / FRENKEL-PAAR-REKOMBINATION
 * ---------------------------------------------------------------------------
 *  Identisch zur seriellen Version, aber:
 *  - Nachbar wird per hash_get() gesucht (nicht per direktem Index).
 *  - Wenn der Nachbar nicht im Strip ist (hash_get=-1), kein Healing von hier;
 *    der andere Rank heilt diese Bindung ggf. von seiner Seite.
 *  - E_broken-Korrektur: nur 0.5*(0.5*K*...) weil die Energie beim Bruch
 *    auch nur halb gebucht wurde (jede Bindung wird von BEIDEN Endpunkten
 *    aus gezaehlt, also 2*(0.5/2)=0.5 insgesamt -- identisches Schema wie
 *    in energies_strip()).
 * ========================================================================== */
static void heal_strip(Strip*s){
    if(!HEALING) return;
    double hd=HEALING_DIST*L0;
    for(int p=0;p<s->n;p++){
        int base=p*MAXNBR;
        for(int e=0;e<s->nb[p];e++){
            if(s->intact[base+e]) continue;            /* schon intakt */
            int q=hash_get(s,s->pg[base+e]); if(q<0) continue; /* Partner unbekannt */
            double dx=s->px[q]-s->px[p], dy=s->py[q]-s->py[p], r=sqrt(dx*dx+dy*dy);
            if(r<RCUT||r>=hd) continue;               /* nicht im Heilfenster */
            double dvx=s->vx[q]-s->vx[p], dvy=s->vy[q]-s->vy[p];
            if(sqrt(dvx*dvx+dvy*dvy)>=HEALING_VREL) continue; /* zu schnell */
            s->intact[base+e]=1;
            s->E_broken -= 0.5*(0.5*K_SPRING*(r-L0)*(r-L0)); /* Energie korrigieren */
        }
    }
}

/* ============================================================================
 *  KRAFTBERECHNUNG  (Herzstueck, dominiert die Laufzeit)
 * ---------------------------------------------------------------------------
 *  Unterschied zur seriellen Version:
 *
 *  (1) FEDERN: Jede Bindung wird NUR VON EINER SEITE berechnet.
 *      Der Gather-Ansatz: Atom p berechnet nur die Kraft, die es SELBST spuert,
 *      nicht die entgegengesetzte auf q.
 *      Aber: wenn q ein Halo-Atom ist (hash_get >= n), erhaelt q nie eine
 *      Kraft zurueck (Halo-Atome kriegen keine Kraefte).
 *      Kompensation: der Nachbar-Rank berechnet dieselbe Bindung von der
 *      anderen Seite (q ist dort resident, p ist Halo-Atom).
 *      -> Jede intakte Bindung wird genau EINMAL als Netto-Kraft gezaehlt.
 *
 *  (2) ABSTOSSUNG: Voller Stern (9 Zellen): jedes Atom p iteriert alle 9
 *      Nachbarzellen (inkl. eigene).  Fuer jedes Paar (p,q) berechnet p
 *      nur seine eigene Kraft (Gather).  q kann ein Halo-Atom sein.
 *      Kein Newton-3-Trick hier -- jeder Thread rechnet fuer sich.
 *
 *      9-Zellen-Stern statt Halb-Stern: Im seriellen Code reichte der Halb-Stern,
 *      weil jedes Paar genau einmal gefunden wurde.  Im MPI-Code hat Atom p nur
 *      Zugriff auf die es umgebenden Atome (nicht auf alle Atome global).
 *      Mit dem Halb-Stern wuerde p die Kraft von q auf sich selbst verpassen,
 *      wenn q in einer "hinteren" Nachbarzelle liegt.
 *      -> Voller Stern notwendig, aber jedes Paar wird dann DOPPELT gezaehlt
 *      (einmal als (p,q) und einmal als (q,p)).
 *      NEIN: da wir nur die Kraft auf P berechnen, nicht auf q, ist es korrekt:
 *      Atom p iteriert alle 9 Zellen und sammelt alle Kraefte, die es selbst
 *      erhaelt.  q-Kraefte werden nicht beruehrt.
 *
 *  (3) KEINE DAEMPFUNG: In dieser Version nicht implementiert.
 *      (Leicht ergaenzbar: nach den Kraefte-Schleifen fx[p] -= DAMPING*vx[p].)
 * ========================================================================== */
static void forces_strip(Strip*s,Cells*c){
    /* Kraefte nur fuer RESIDENTE Atome (0..n-1) zuruecksetzen und berechnen.
     * Halo-Atome haben keine Eintraege in fx/fy (kein Zurueckschicken). */
    for(int p=0;p<s->n;p++){ s->fx[p]=0; s->fy[p]=0; }
    double rc2=RCUT*RCUT;
    /* Voller Stern: alle 9 Nachbarzellen + eigene Zelle */
    static const int off[9][2]={{0,0},{1,0},{-1,0},{0,1},{0,-1},{1,1},{1,-1},{-1,1},{-1,-1}};
    for(int p=0;p<s->n;p++){
        int base=p*MAXNBR;
        /* (a) Federn zu allen Bindungspartnern (auch Halo-Atomen) */
        for(int e=0;e<s->nb[p];e++){
            if(!s->intact[base+e]) continue;
            int q=hash_get(s,s->pg[base+e]);
            /* hash_get -1: Partner nicht im Strip (weder resident noch Halo).
             * Das passiert bei Bindungen zum uebernachsten Nachbarstreifen -- sehr
             * selten (nur wenn Atome sehr weit auswandern). Kraft wird ignoriert;
             * der andere Rank berechnet sie von seiner Seite. */
            if(q<0) continue;
            double dx=s->px[q]-s->px[p], dy=s->py[q]-s->py[p], r=sqrt(dx*dx+dy*dy);
            /* Bruch: auch wenn q Halo ist (wir haben die gid -> intakt-Bit liegt bei uns) */
            if(r > L0*MAX_STRETCH){
                s->intact[base+e]=0;
                s->E_broken += 0.5*(0.5*K_SPRING*(r-L0)*(r-L0)); /* halbe Energie! */
                continue;
            }
            double f=K_SPRING*(r-L0)/r;
            s->fx[p]+=f*dx; s->fy[p]+=f*dy;  /* Gather: nur Kraft auf p */
        }
        /* (b) Abstossung via Cell List (voller Stern) */
        int cx=(int)((s->px[p]-c->x0)/c->cs), cy=(int)((s->py[p]-c->y0)/c->cs);
        if(cx<0)cx=0; if(cx>=c->ncx)cx=c->ncx-1;
        if(cy<0)cy=0; if(cy>=c->ncy)cy=c->ncy-1;
        for(int o=0;o<9;o++){
            int nxc=cx+off[o][0], nyc=cy+off[o][1];
            if(nxc<0||nxc>=c->ncx||nyc<0||nyc>=c->ncy) continue;
            for(int q=c->head[nyc*c->ncx+nxc]; q!=-1; q=c->next[q]){
                if(q==p) continue;           /* kein Atom mit sich selbst */
                double dx=s->px[p]-s->px[q], dy=s->py[p]-s->py[q], r2=dx*dx+dy*dy;
                if(r2>=rc2||r2==0.0) continue;
                if(bonded_to(s,p,s->gid[q])) continue; /* gebunden -> keine Abstossung */
                double r=sqrt(r2), f=rep_force(r)/r;
                s->fx[p]+=f*dx; s->fy[p]+=f*dy;    /* Gather: nur Kraft auf p */
            }
        }
    }
}

/* ============================================================================
 *  ENERGIEBERECHNUNG  (nur fuer Diagnose, alle logevery Schritte)
 * ---------------------------------------------------------------------------
 *  Halbes-Zaehlen-Schema:
 *  Jede Bindung (p,q) wird von BEIDEN Seiten aus betrachtet.  Damit jede Energie
 *  nur einmal gezaehlt wird, multipliziert man mit 0.5:
 *    E_spring = Summe ueber alle (p,e): 0.5 * 0.5 * K * (r-L0)^2
 *               ---  Pro Atom 0.5 * Haelfte per Atom = Viertel, aber zwei Seiten -> korrekt.
 *
 *  Fuer die Abstossung gleiches Schema.
 *
 *  Die lokalen Energien werden per MPI_Reduce summiert (in main()).
 * ========================================================================== */
static void energies_strip(Strip*s,Cells*c,double*ek,double*es,double*er,double*eb){
    double sk=0,ss=0,sr=0;
    double rc2=RCUT*RCUT;
    static const int off[9][2]={{0,0},{1,0},{-1,0},{0,1},{0,-1},{1,1},{1,-1},{-1,1},{-1,-1}};
    for(int p=0;p<s->n;p++){
        sk+=0.5*s->mass[p]*(s->vx[p]*s->vx[p]+s->vy[p]*s->vy[p]);
        int base=p*MAXNBR;
        for(int e=0;e<s->nb[p];e++){
            if(!s->intact[base+e])continue;
            int q=hash_get(s,s->pg[base+e]); if(q<0)continue;
            double dx=s->px[q]-s->px[p],dy=s->py[q]-s->py[p],r=sqrt(dx*dx+dy*dy);
            ss+=0.5*(0.5*K_SPRING*(r-L0)*(r-L0));  /* 0.5 wegen beidseitigem Zaehlen */
        }
        int cx=(int)((s->px[p]-c->x0)/c->cs), cy=(int)((s->py[p]-c->y0)/c->cs);
        if(cx<0)cx=0; if(cx>=c->ncx)cx=c->ncx-1;
        if(cy<0)cy=0; if(cy>=c->ncy)cy=c->ncy-1;
        for(int o=0;o<9;o++){
            int nxc=cx+off[o][0],nyc=cy+off[o][1];
            if(nxc<0||nxc>=c->ncx||nyc<0||nyc>=c->ncy)continue;
            for(int q=c->head[nyc*c->ncx+nxc]; q!=-1; q=c->next[q]){
                if(q==p)continue;
                double dx=s->px[p]-s->px[q],dy=s->py[p]-s->py[q],r2=dx*dx+dy*dy;
                if(r2>=rc2||r2==0.0)continue;
                if(bonded_to(s,p,s->gid[q]))continue;
                sr+=0.5*rep_pot(sqrt(r2));  /* 0.5: Abstossung auch beidseitig gezaehlt */
            }
        }
    }
    *ek=sk; *es=ss; *er=sr; *eb=s->E_broken;
}

/* ============================================================================
 *  PACK / UNPACK / REMOVE  --  Atom-Serialisierung fuer MPI-Puffer
 * ---------------------------------------------------------------------------
 *  ATOM_STRIDE=19: Anzahl doubles pro Atom im Migrations-Puffer:
 *    [0]  gid (als double, Praezision reicht: gid < NX*NY < 2e7 << 2^53)
 *    [1]  px, [2] py, [3] vx, [4] vy, [5] mass, [6] nb (Nachbarzahl)
 *    [7..12]  pg[0..5]   (gids der 6 moeglichen Bindungspartner)
 *    [13..18] intact[0..5]
 *
 *  HALO_STRIDE=5: Nur Position + Geschwindigkeit:
 *    [0] gid, [1] px, [2] py, [3] vx, [4] vy
 *    (Halo-Atome brauchen keine Bindungsinformationen: sie kriegen keine Kraefte)
 *
 *  remove_atom(): "Swap-with-last" -- das letzte Element ueberschreibt den
 *  geloeschten Platz.  O(1), aber veraendert die Reihenfolge der Atome.
 *  Kein Problem: Reihenfolge ist willkuerlich; Hash-Tabelle wird neu gebaut.
 * ========================================================================== */
#define ATOM_STRIDE 19
#define HALO_STRIDE 5

static void pack_atom(double*buf,Strip*s,int p){
    buf[0]=(double)s->gid[p]; buf[1]=s->px[p]; buf[2]=s->py[p];
    buf[3]=s->vx[p]; buf[4]=s->vy[p]; buf[5]=s->mass[p]; buf[6]=(double)s->nb[p];
    for(int e=0;e<MAXNBR;e++){
        buf[7+e] =(double)s->pg    [p*MAXNBR+e];
        buf[13+e]=(double)s->intact[p*MAXNBR+e];
    }
}

static void unpack_atom(double*buf,Strip*s){
    int k=s->n++;    /* neues Atom am Ende anhaengen */
    s->gid[k]=(long)buf[0]; s->px[k]=buf[1]; s->py[k]=buf[2];
    s->vx[k]=buf[3]; s->vy[k]=buf[4]; s->mass[k]=buf[5]; s->nb[k]=(int)buf[6];
    for(int e=0;e<MAXNBR;e++){
        s->pg    [k*MAXNBR+e]=(long)buf[7+e];
        s->intact[k*MAXNBR+e]=(int) buf[13+e];
    }
}

static void remove_atom(Strip*s,int p){
    int last=--s->n;
    if(p!=last){
        s->gid[p]=s->gid[last];
        s->px[p]=s->px[last]; s->py[p]=s->py[last];
        s->vx[p]=s->vx[last]; s->vy[p]=s->vy[last];
        s->mass[p]=s->mass[last];
        s->nb[p]=s->nb[last];
        memcpy(&s->pg    [p*MAXNBR],&s->pg    [last*MAXNBR],MAXNBR*sizeof(long));
        memcpy(&s->intact[p*MAXNBR],&s->intact[last*MAXNBR],MAXNBR*sizeof(int));
    }
}

/* ============================================================================
 *  MPI-KOMMUNIKATION
 * ========================================================================== */
static int RANK,SIZE,PREV,NEXT;
static MPI_Comm CART;   /* kartesischer 1D-Communicator */

/* ============================================================================
 *  sendrecv_var  --  Hilffunktion fuer variable Datenmengen
 * ---------------------------------------------------------------------------
 *  Protokoll in zwei Phasen:
 *  Phase 1: Schicke nsend (int) an dst, empfange nrecv (int) von src.
 *           Beide Seiten wissen jetzt, wie viel Daten kommen.
 *  Phase 2: Schicke sbuf[0..nsend-1] an dst, empfange rbuf[0..nrecv-1] von src.
 *
 *  Warum zwei Phasen?  MPI_Sendrecv benoetigt die Empfangsgroesse im Voraus
 *  (anders als MPI_Probe + MPI_Recv, das ist teurer).
 *
 *  MPI_PROC_NULL an PREV (Rank 0) oder NEXT (letzter Rank):
 *  MPI_Sendrecv an MPI_PROC_NULL ist ein No-Op -> kein if-Sonderfall noetig.
 * ========================================================================== */
static int sendrecv_var(double*sbuf,int nsend,int dst,
                        double*rbuf,int rcap,int src,int tag){
    int nrecv=0;
    MPI_Sendrecv(&nsend,1,MPI_INT,dst,tag, &nrecv,1,MPI_INT,src,tag,
                 CART,MPI_STATUS_IGNORE);
    if(nrecv>rcap) nrecv=rcap;  /* Sicherheitsklampe (sollte nie noetig sein) */
    MPI_Sendrecv(sbuf,nsend,MPI_DOUBLE,dst,tag+1, rbuf,nrecv,MPI_DOUBLE,src,tag+1,
                 CART,MPI_STATUS_IGNORE);
    return nrecv;
}

/* ============================================================================
 *  MIGRATION  --  Atome die Streifengrenzen ueberquert haben, verschieben
 * ---------------------------------------------------------------------------
 *  Jedes Schritt: pruefen ob ein Atom ausserhalb [ylo, yhi) liegt.
 *  - y >= yhi: Atom ist nach oben ausgewandert -> an NEXT schicken.
 *  - y <  ylo: Atom ist nach unten ausgewandert -> an PREV schicken.
 *
 *  WARUM kein Ueberspringen von Ranks?
 *  Mit v*dt << Streifenbreite kann kein Atom in einem Schritt ueber zwei
 *  Streifengrenzen wandern.  Das gilt, solange dt klein genug ist (was
 *  ohnehin fuer Stabilitaet noetig ist).
 *
 *  Puffer-Layout: sUp[nUp*ATOM_STRIDE], sDn[nDn*ATOM_STRIDE].
 *  Erst alle Atome in die Sendepuffer packen, dann zweimal sendrecv_var:
 *    - Zuerst nach oben senden (empfange von unten)
 *    - Dann nach unten senden (empfange von oben)
 * ========================================================================== */
static void migrate(Strip*s,double*sUp,double*sDn,double*rbuf,int bufcap){
    int nUp=0,nDn=0;
    for(int p=0;p<s->n;){
        if(s->py[p] >= s->yhi && NEXT!=MPI_PROC_NULL){
            pack_atom(&sUp[nUp*ATOM_STRIDE],s,p); nUp++; remove_atom(s,p);
            /* remove_atom macht swap-with-last, also p NICHT inkrementieren */
        } else if(s->py[p] < s->ylo && PREV!=MPI_PROC_NULL){
            pack_atom(&sDn[nDn*ATOM_STRIDE],s,p); nDn++; remove_atom(s,p);
        } else p++;
    }
    /* Tag 100: nach oben (NEXT) senden, von unten (PREV) empfangen */
    int rn=sendrecv_var(sUp,nUp*ATOM_STRIDE,NEXT, rbuf,bufcap,PREV,100);
    for(int o=0;o<rn;o+=ATOM_STRIDE) unpack_atom(&rbuf[o],s);
    /* Tag 102: nach unten (PREV) senden, von oben (NEXT) empfangen */
    rn=sendrecv_var(sDn,nDn*ATOM_STRIDE,PREV, rbuf,bufcap,NEXT,102);
    for(int o=0;o<rn;o+=ATOM_STRIDE) unpack_atom(&rbuf[o],s);
}

/* ============================================================================
 *  HALO-AUSTAUSCH  (blocking)
 * ---------------------------------------------------------------------------
 *  Schickt Kopien der Randatome (nur gid, px, py, vx, vy) an den Nachbarn.
 *  Der Nachbar braucht diese "Geister"-Atome fuer:
 *    - Federkraefte ueber die Grenze (wenn ein Bindungspartner druebenliegt)
 *    - Abstossung ueber die Grenze (wenn Abstand < RCUT)
 *
 *  halo_w = Halobreite:  max(RCUT, L0*MAX_STRETCH) + 5% Sicherheitsmarge.
 *  Atome innerhalb halo_w vom Rand werden ins Halo-Array des Nachbarn kopiert.
 *
 *  s->nh wird am Anfang auf 0 gesetzt (alter Halo verworfen).
 *  Neue Halo-Atome werden ab Index s->n angehaengt.
 *
 *  BLOCKING: Zuerst Up-Richtung komplett (senden + empfangen), dann Down.
 *  Die zwei Richtungen laufen seriell -- einfacher, aber etwas langsamer
 *  als non-blocking.
 * ========================================================================== */
static void exchange_halo(Strip*s,double*sUp,double*sDn,double*rbuf,int bufcap,double halo_w){
    s->nh=0;  /* Halo-Atome vom letzten Schritt verwerfen */
    int nUp=0,nDn=0;
    for(int p=0;p<s->n;p++){
        /* Nahe der oberen Grenze: an NEXT schicken */
        if(s->py[p] >= s->yhi-halo_w && NEXT!=MPI_PROC_NULL){
            double*b=&sUp[nUp*HALO_STRIDE];
            b[0]=(double)s->gid[p]; b[1]=s->px[p]; b[2]=s->py[p];
            b[3]=s->vx[p]; b[4]=s->vy[p]; nUp++;
        }
        /* Nahe der unteren Grenze: an PREV schicken */
        if(s->py[p] < s->ylo+halo_w && PREV!=MPI_PROC_NULL){
            double*b=&sDn[nDn*HALO_STRIDE];
            b[0]=(double)s->gid[p]; b[1]=s->px[p]; b[2]=s->py[p];
            b[3]=s->vx[p]; b[4]=s->vy[p]; nDn++;
        }
    }
    /* Nach oben: empfange Halo von unten (PREV) */
    int rn=sendrecv_var(sUp,nUp*HALO_STRIDE,NEXT, rbuf,bufcap,PREV,200);
    for(int o=0;o<rn;o+=HALO_STRIDE){
        int k=s->n+s->nh; s->nh++;
        s->gid[k]=(long)rbuf[o]; s->px[k]=rbuf[o+1]; s->py[k]=rbuf[o+2];
        s->vx[k]=rbuf[o+3]; s->vy[k]=rbuf[o+4];
        s->mass[k]=1.0; s->nb[k]=0;  /* mass/nb braucht gueltigen Wert */
    }
    /* Nach unten: empfange Halo von oben (NEXT) */
    rn=sendrecv_var(sDn,nDn*HALO_STRIDE,PREV, rbuf,bufcap,NEXT,202);
    for(int o=0;o<rn;o+=HALO_STRIDE){
        int k=s->n+s->nh; s->nh++;
        s->gid[k]=(long)rbuf[o]; s->px[k]=rbuf[o+1]; s->py[k]=rbuf[o+2];
        s->vx[k]=rbuf[o+3]; s->vy[k]=rbuf[o+4];
        s->mass[k]=1.0; s->nb[k]=0;
    }
}

/* ============================================================================
 *  HALO-AUSTAUSCH  (non-blocking)
 * ---------------------------------------------------------------------------
 *  Wie exchange_halo(), aber beide Richtungen (Up+Down) laufen GLEICHZEITIG.
 *
 *  Vorteil: Die Wartezeit auf den Nachbarn wird halbiert, weil beide
 *  Richtungen parallel ablaufen.  Bei 2 Ranks spart das gar nichts (nur
 *  ein Nachbar), bei vielen Ranks ist der Gewinn gering aber messbar (~5%).
 *
 *  WICHTIG: Zwei getrennte Empfangspuffer (rUp, rDn) sind noetig, weil
 *  beide Empfaenge gleichzeitig laufen.  Wuerde man denselben rbuf-Puffer
 *  nutzen, koennten sich die Daten ueberschreiben.
 *
 *  Protokoll:
 *   Phase 1: Alle 4 Isend/Irecv fuer die Anzahlen starten, Waitall.
 *   Phase 2: Alle 4 Isend/Irecv fuer die Daten starten (wenn >0), Waitall.
 *  MPI_PROC_NULL-Requests sind automatisch No-Ops (kein if noetig,
 *  aber wir pruefen trotzdem ob n>0 um die Request-Liste sauber zu halten).
 *
 *  Naechster Schritt (nicht implementiert): innere Kraefte waehrend des Halos
 *  berechnen (Communication-Computation-Overlap).  Dafuer muesste man die
 *  Rand-Atome von den inneren Atomen trennen.
 * ========================================================================== */
static void exchange_halo_nb(Strip*s,double*sUp,double*sDn,
                             double*rUp,double*rDn,int bufcap,double halo_w){
    s->nh=0;
    int nUp=0,nDn=0;
    for(int p=0;p<s->n;p++){
        if(s->py[p] >= s->yhi-halo_w && NEXT!=MPI_PROC_NULL){
            double*b=&sUp[nUp*HALO_STRIDE];
            b[0]=(double)s->gid[p]; b[1]=s->px[p]; b[2]=s->py[p];
            b[3]=s->vx[p]; b[4]=s->vy[p]; nUp++;
        }
        if(s->py[p] < s->ylo+halo_w && PREV!=MPI_PROC_NULL){
            double*b=&sDn[nDn*HALO_STRIDE];
            b[0]=(double)s->gid[p]; b[1]=s->px[p]; b[2]=s->py[p];
            b[3]=s->vx[p]; b[4]=s->vy[p]; nDn++;
        }
    }
    int sUpN=nUp*HALO_STRIDE, sDnN=nDn*HALO_STRIDE, rUpN=0, rDnN=0;
    MPI_Request rq[4]; int nq;

    /* Phase 1: Anzahlen austauschen -- alle 4 gleichzeitig */
    nq=0;
    MPI_Irecv(&rUpN,1,MPI_INT,PREV,200,CART,&rq[nq++]);  /* empfange von unten */
    MPI_Irecv(&rDnN,1,MPI_INT,NEXT,202,CART,&rq[nq++]);  /* empfange von oben  */
    MPI_Isend(&sUpN,1,MPI_INT,NEXT,200,CART,&rq[nq++]);  /* sende nach oben    */
    MPI_Isend(&sDnN,1,MPI_INT,PREV,202,CART,&rq[nq++]);  /* sende nach unten   */
    MPI_Waitall(nq,rq,MPI_STATUSES_IGNORE);
    if(rUpN>bufcap) rUpN=bufcap;
    if(rDnN>bufcap) rDnN=bufcap;

    /* Phase 2: Daten austauschen -- alle 4 gleichzeitig */
    nq=0;
    if(rUpN>0) MPI_Irecv(rUp,rUpN,MPI_DOUBLE,PREV,201,CART,&rq[nq++]);
    if(rDnN>0) MPI_Irecv(rDn,rDnN,MPI_DOUBLE,NEXT,203,CART,&rq[nq++]);
    if(sUpN>0) MPI_Isend(sUp,sUpN,MPI_DOUBLE,NEXT,201,CART,&rq[nq++]);
    if(sDnN>0) MPI_Isend(sDn,sDnN,MPI_DOUBLE,PREV,203,CART,&rq[nq++]);
    MPI_Waitall(nq,rq,MPI_STATUSES_IGNORE);

    /* Halo-Atome in das Strip-Array eintragen */
    for(int o=0;o<rUpN;o+=HALO_STRIDE){
        int k=s->n+s->nh; s->nh++;
        s->gid[k]=(long)rUp[o]; s->px[k]=rUp[o+1]; s->py[k]=rUp[o+2];
        s->vx[k]=rUp[o+3]; s->vy[k]=rUp[o+4]; s->mass[k]=1.0; s->nb[k]=0;
    }
    for(int o=0;o<rDnN;o+=HALO_STRIDE){
        int k=s->n+s->nh; s->nh++;
        s->gid[k]=(long)rDn[o]; s->px[k]=rDn[o+1]; s->py[k]=rDn[o+2];
        s->vx[k]=rDn[o+3]; s->vy[k]=rDn[o+4]; s->mass[k]=1.0; s->nb[k]=0;
    }
}

/* ============================================================================
 *  CELL-LINEARE UMSORTIERUNG (resort_strip)
 * ---------------------------------------------------------------------------
 *  PROBLEM: Durch remove_atom (swap-with-last) wird die Reihenfolge der
 *  Atome im Array mit jedem Schritt zufaelliger.  Nach 1000 Schritten ist
 *  die raeumliche Cache-Lokalitaet vollstaendig zerstoert.
 *
 *  LOESUNG: Alle RESIDENTEN Atome nach ihrer Zell-Position sortieren
 *  (cell-linear: Zeile fuer Zeile von unten nach oben, von links nach rechts).
 *  Das ist kein physikalisches Sortieren, nur ein Umbenennen der Array-Indizes.
 *  Die Physik ist identisch, weil Bindungen ueber globale IDs laufen und
 *  die Hash-Tabelle + Cell List jeden Schritt neu aufgebaut werden.
 *
 *  WANN?: Nur zwischen Migration und Halo-Austausch, und nur auf residenten
 *  Atomen (nh=0 zu diesem Zeitpunkt). fx/fy muessen nicht umgeordnet werden
 *  (werden ohnehin vor Gebrauch neu berechnet).
 *
 *  KOSTEN: O(N log N) durch qsort, aber sehr kleines N/P.  Der Gewinn durch
 *  bessere Cache-Nutzung amortisiert die Sortierkosten nach ~10-20 Schritten.
 * ========================================================================== */
static long *g_rskey;
static int cmp_rskey(const void*a,const void*b){
    long ka=g_rskey[*(const int*)a], kb=g_rskey[*(const int*)b];
    return (ka>kb)-(ka<kb);
}
static void resort_strip(Strip*s){
    int n=s->n; if(n<2) return;
    double x0=1e30,y0=1e30,x1=-1e30;
    for(int p=0;p<n;p++){
        double X=s->px[p],Y=s->py[p];
        if(X<x0)x0=X; if(X>x1)x1=X; if(Y<y0)y0=Y;
    }
    double cs=RCUT; int ncx=(int)((x1-x0)/cs)+2;
    long*key=malloc(n*sizeof*key); int*perm=malloc(n*sizeof*perm);
    /* Schluessel = linearisierter Zellenindex */
    for(int p=0;p<n;p++){
        int cx=(int)((s->px[p]-x0)/cs), cy=(int)((s->py[p]-y0)/cs);
        if(cx<0)cx=0; if(cy<0)cy=0;
        key[p]=(long)cy*ncx+cx; perm[p]=p;
    }
    g_rskey=key; qsort(perm,n,sizeof(int),cmp_rskey);
    /* Alle Felder gemaess der Sortier-Permutation umordnen */
    double*td=malloc(n*sizeof(double));
    long  *tg=malloc(n*sizeof(long));
    int   *ti=malloc(n*sizeof(int));
    long  *tpg=malloc((size_t)n*MAXNBR*sizeof(long));
    int   *tin=malloc((size_t)n*MAXNBR*sizeof(int));
    #define RS_D(A) do{ for(int p=0;p<n;p++) td[p]=s->A[perm[p]]; \
                        memcpy(s->A,td,n*sizeof(double)); }while(0)
    RS_D(px); RS_D(py); RS_D(vx); RS_D(vy); RS_D(mass);
    #undef RS_D
    for(int p=0;p<n;p++) tg[p]=s->gid[perm[p]]; memcpy(s->gid,tg,n*sizeof(long));
    for(int p=0;p<n;p++) ti[p]=s->nb [perm[p]]; memcpy(s->nb, ti,n*sizeof(int));
    for(int p=0;p<n;p++){
        memcpy(&tpg[p*MAXNBR],&s->pg    [perm[p]*MAXNBR],MAXNBR*sizeof(long));
        memcpy(&tin[p*MAXNBR],&s->intact[perm[p]*MAXNBR],MAXNBR*sizeof(int));
    }
    memcpy(s->pg,    tpg,(size_t)n*MAXNBR*sizeof(long));
    memcpy(s->intact,tin,(size_t)n*MAXNBR*sizeof(int));
    free(key);free(perm);free(td);free(tg);free(ti);free(tpg);free(tin);
}



/* ============================================================================
 *  MAIN  --  Hauptprogramm
 * ---------------------------------------------------------------------------
 *  Ablauf:
 *  1.  MPI initialisieren, Rank/Size bestimmen
 *  2.  Konfiguration lesen (Rank 0) und broadcasten
 *  3.  1D-kartesische Topologie erstellen (PREV/NEXT Nachbarn)
 *  4.  Lokale Atome initialisieren (nur die im eigenen y-Streifen)
 *  5.  PKA(s) setzen (alle Ranks berechnen gleiche Liste; setzen nur wenn lokal)
 *  6.  Startenergie (MPI_Reduce der lokalen Energien)
 *  7.  Hauptschleife: Velocity-Verlet + Kommunikation
 *  8.  Phasen-Timing ausgeben
 * ========================================================================== */
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
    MPI_Comm_size(MPI_COMM_WORLD,&SIZE);
    const char*cfg=(argc>1)?argv[1]:"params.ini";
    if(RANK==0) read_config(cfg);
    bcast_config();
    DY=L0*sqrt(3.0)/2.0;

    /* ---- 1D-Kartesische Topologie ----------------------------------------
       dims={SIZE}: alle Ranks in einer Reihe (y-Richtung).
       periods={0}: keine periodischen Randbedingungen (an den Raendern ist
                    PREV oder NEXT == MPI_PROC_NULL).
       reorder=1:   erlaubt MPI, Ranks umzuordnen fuer bessere physische Naehe
                    (nur relevant auf Multi-Knoten-Systemen).
       MPI_Cart_shift(cart,0,1,...): gibt die Nachbarn in Dimension 0 mit
                    Verschiebung +1 zurueck.                                    */
    MPI_Comm cart; int dims[1]={SIZE}, periods[1]={0};
    MPI_Cart_create(MPI_COMM_WORLD,1,dims,periods,1,&cart);
    int crank; MPI_Comm_rank(cart,&crank); RANK=crank;
    MPI_Cart_shift(cart,0,1,&PREV,&NEXT);
    CART=cart;

    /* Halobreite: gross genug fuer Federn UND Abstossung + kleiner Sicherheitsrand */
    double halo_w=(RCUT>L0*MAX_STRETCH?RCUT:L0*MAX_STRETCH)+0.05*L0;
    double Ly=(NY-1)*DY;    /* Gesamtausdehnung in y */
    Strip s;
    long Ntot=(long)NX*NY;
    /* Puffer-Kapazitaet: N/P + Headroom fuer Migration (mehrere Reihen) + Halo */
    int capeach=(int)(Ntot/SIZE)+8*NX+128;
    strip_alloc(&s,capeach);

    /* y-Grenzen dieses Streifens.
       Rank 0: ylo=-inf (kein PREV), yhi=Ly/P.
       Letzer Rank: yhi=+inf (kein NEXT), ylo=(P-1)*Ly/P.
       Rand = -1e30 bzw. +1e30 -> alle Atome "darunter/darueber" gehoeren hierher
       (Rand des Gesamtgitters hat keine Nachbarn ausserhalb). */
    s.ylo=(RANK==0)      ? -1e30 : (double)RANK*Ly/SIZE - 0.5*DY;
    s.yhi=(RANK==SIZE-1) ?  1e30 : (double)(RANK+1)*Ly/SIZE - 0.5*DY;

    /* ---- Lokale Atome initialisieren ------------------------------------ */
    for(long g=0; g<Ntot; g++){
        double y=gy_of(g);
        if(y<s.ylo||y>=s.yhi) continue;  /* ausserhalb meines Streifens */
        int p=s.n++; s.gid[p]=g;
        s.px[p]=gx_of(g); s.py[p]=y;
        s.vx[p]=0; s.vy[p]=0; s.mass[p]=1.0;
        long nb[6]; int nn=neighbors_of(g,nb); s.nb[p]=nn;
        /* Bindungsliste: alle 6 Nachbarn (auch die in anderen Streifen) */
        for(int e=0;e<nn;e++){ s.pg[p*MAXNBR+e]=nb[e]; s.intact[p*MAXNBR+e]=1; }
    }

    /* ---- PKA(s) setzen --------------------------------------------------
       Alle Ranks berechnen DIESELBE PKA-Liste (gleicher RNG, gleicher SEED).
       Jeder Rank prueft: liegt das PKA in meinem Streifen?  Wenn ja, setzen. */
    if(N_MANUAL>0){
        for(int k=0,placed=0;k<N_MANUAL && placed<N_PKA;k++){
            PkaDef*p=&MANUAL_PKA[k];
            if(!(p->hx||p->hy||p->he||p->ha||p->hm)) continue;
            double xr=p->hx?p->x:PKA_X, yr=p->hy?p->y:PKA_Y;
            double en=p->he?p->energy:PKA_ENERGY, ang=p->ha?p->angle:PKA_ANGLE;
            double mm=p->hm?p->mass:PKA_MASS;
            long pj=(long)(yr*(NY-1)+0.5), pi=(long)(xr*(NX-1)+0.5);
            long pka_gid=pj*NX+pi;
            for(int q=0;q<s.n;q++) if(s.gid[q]==pka_gid){
                s.mass[q]=mm;
                double v=sqrt(2.0*en/s.mass[q]), th=ang*M_PI/180.0;
                s.vx[q]=v*sin(th); s.vy[q]=-v*cos(th);
            }
            placed++;
        }
    } else {
        RNG = (unsigned long)SEED*2654435761UL + 12345UL;
        for(int kp=0; kp<N_PKA; kp++){
            double fx_,fy_,ang;
            if(kp==0 && N_PKA==1){ fx_=PKA_X; fy_=PKA_Y; ang=PKA_ANGLE; }
            else { fx_=0.15+0.70*rng_u(); fy_=0.10+0.80*rng_u(); ang=360.0*rng_u(); }
            long pj=(long)(fy_*(NY-1)+0.5), pi=(long)(fx_*(NX-1)+0.5);
            long pka_gid=pj*NX+pi;
            for(int p=0;p<s.n;p++) if(s.gid[p]==pka_gid){
                s.mass[p]=PKA_MASS;
                double v=sqrt(2.0*PKA_ENERGY/s.mass[p]), th=ang*M_PI/180.0;
                s.vx[p]=v*sin(th); s.vy[p]=-v*cos(th);
            }
        }
    }

    /* ---- Kommunikationspuffer allokieren --------------------------------
       Puffer gross genug fuer mehrere Gitterzeilen Atome (haufigster Migrationsfall) */
    int bufcap=(8*NX+128)*ATOM_STRIDE;
    double*sUp=malloc(bufcap*sizeof(double));   /* Sendepuffer nach oben   */
    double*sDn=malloc(bufcap*sizeof(double));   /* Sendepuffer nach unten  */
    double*rbuf=malloc(bufcap*sizeof(double));  /* Empfangspuffer (blocking) */
    double*rbuf2=malloc(bufcap*sizeof(double)); /* 2. Empfangspuffer (non-blocking, up-Richtung) */

    /* ---- Initialer Halo + Startenergie ---------------------------------- */
    exchange_halo(&s,sUp,sDn,rbuf,bufcap,halo_w);
    hash_build(&s);
    double ek,es,er,eb,loc[4],glob[4],e0;
    { Cells c; cells_build(&s,&c); energies_strip(&s,&c,&ek,&es,&er,&eb); cells_free(&c); }
    loc[0]=ek;loc[1]=es;loc[2]=er;loc[3]=eb;
    MPI_Reduce(loc,glob,4,MPI_DOUBLE,MPI_SUM,0,CART);
    e0=glob[0]+glob[1]+glob[2]+glob[3];

    if(RANK==0){
        printf("# MPI Modell A | ranks=%d | N=%ld | halo_w=%.3f | healing=%d\n",
               SIZE,Ntot,halo_w,HEALING);
        printf("# E0=%.6f\n# step      t        E_total       drift%%\n",e0);
    }
    int logevery=(LOG_EVERY>0)?LOG_EVERY:NSTEPS/40; if(logevery<1)logevery=1;

    /* Rank 0 schreibt CSV-Energieverlauf (wie die serielle Version) */
    FILE*ecsv=NULL;
    if(RANK==0){
        char fn[160]; snprintf(fn,sizeof fn,"%s_energy.csv",OUT_PREFIX);
        ecsv=fopen(fn,"w");
        if(ecsv) fprintf(ecsv,"step,t,E_kin,E_spring,E_rep,E_broken,E_total\n");
    }

    /* ---- Stoppuhr starten ----------------------------------------------- */
    MPI_Barrier(CART);     /* synchronisieren vor Zeitmessung */
    double t_start=MPI_Wtime();

    /* ============ HAUPTSCHLEIFE =========================================== */
    for(int step=0; step<=NSTEPS; step++){
        /* --- Energie-Diagnose ---
           MPI_Reduce summiert lokale Energien auf Rank 0. */
        if(step%logevery==0){
            double te0=MPI_Wtime();
            { Cells c; cells_build(&s,&c); energies_strip(&s,&c,&ek,&es,&er,&eb); cells_free(&c); }
            loc[0]=ek;loc[1]=es;loc[2]=er;loc[3]=eb;
            MPI_Reduce(loc,glob,4,MPI_DOUBLE,MPI_SUM,0,CART);
            T_energy+=MPI_Wtime()-te0;
            if(RANK==0){
                double et=glob[0]+glob[1]+glob[2]+glob[3];
                printf("%6d %8.3f %13.4f %8.4f\n",step,step*DT,et,
                       (e0!=0)?100.0*(et-e0)/fabs(e0):0.0);
                if(ecsv) fprintf(ecsv,"%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                       step,step*DT,glob[0],glob[1],glob[2],glob[3],et);
            }
        }
        if(step==NSTEPS) break;
        double tt;

        /* --- (A) Halber Kick + Drift ---
           (identisch zur seriellen Version, aber nur fuer RESIDENTE Atome) */
        for(int p=0;p<s.n;p++){
            s.vx[p]+=0.5*DT*s.fx[p]/s.mass[p];
            s.vy[p]+=0.5*DT*s.fy[p]/s.mass[p];
            s.px[p]+=DT*s.vx[p];
            s.py[p]+=DT*s.vy[p];
        }

        /* --- (B) Migration: Atome die ihren Streifen verlassen haben, uebergeben --- */
        tt=MPI_Wtime(); migrate(&s,sUp,sDn,rbuf,bufcap);
        /* Optionale Cache-Neusortierung nach Migration (vor Halo, nur Residente) */
        if(RESORT_EVERY>0 && step>0 && step%RESORT_EVERY==0) resort_strip(&s);
        T_migrate+=MPI_Wtime()-tt;

        /* --- (C) Halo-Austausch: Randkopien an Nachbarn schicken --- */
        tt=MPI_Wtime();
        if(HALO_NB) exchange_halo_nb(&s,sUp,sDn,rbuf,rbuf2,bufcap,halo_w);
        else        exchange_halo   (&s,sUp,sDn,rbuf,    bufcap,halo_w);
        T_halo+=MPI_Wtime()-tt;

        /* --- (D) Hash + Cell List + Healing + Kraefte ---
           Hash muss VOR Kraefte-Berechnung aufgebaut werden (Atom-Migration
           aendert lokale Indizes). Cell List umfasst Residente + Halo. */
        tt=MPI_Wtime(); hash_build(&s);                                T_hash+=MPI_Wtime()-tt;
        { Cells c;
          tt=MPI_Wtime(); cells_build(&s,&c);                          T_cells+=MPI_Wtime()-tt;
          heal_strip(&s);                                              /* Healing zuerst! */
          tt=MPI_Wtime(); forces_strip(&s,&c);                         T_force+=MPI_Wtime()-tt;
          cells_free(&c); }

        /* --- (E) Zweiter halber Kick --- */
        for(int p=0;p<s.n;p++){
            s.vx[p]+=0.5*DT*s.fx[p]/s.mass[p];
            s.vy[p]+=0.5*DT*s.fy[p]/s.mass[p];
        }
    }
    /* ============ ENDE HAUPTSCHLEIFE ====================================== */

    MPI_Barrier(CART);
    double t_loc=MPI_Wtime()-t_start, t_max;
    if(ecsv) fclose(ecsv);
    /* t_max = Wandzeit des langsamsten Ranks (der bestimmt die Gesamtlaufzeit) */
    MPI_Reduce(&t_loc,&t_max,1,MPI_DOUBLE,MPI_MAX,0,CART);

    /* ---- Ergebnis: gerissene Bindungen (Half-Edge-Count / 2) ------------ */
    long bloc=0;
    for(int p=0;p<s.n;p++){
        int base=p*MAXNBR;
        for(int e=0;e<s.nb[p];e++) if(!s.intact[base+e]) bloc++;
    }
    /* bloc zaehlt Half-Edges: jede gerissene Bindung zweimal (von beiden Enden). */
    long bglob=0; int nloc=s.n, nglob=0;
    MPI_Reduce(&bloc,&bglob,1,MPI_LONG,MPI_SUM,0,CART);
    MPI_Reduce(&nloc,&nglob,1,MPI_INT,MPI_SUM,0,CART);
    if(RANK==0){
        printf("# RESULT ranks=%d atoms=%d (soll %ld) broken_bonds~=%ld healing=%d seed=%u\n",
               SIZE,nglob,Ntot,bglob/2,HEALING,SEED);
        printf("# TIME max_walltime=%.4f s  (ranks=%d, N=%ld, steps=%d)\n",
               t_max,SIZE,Ntot,NSTEPS);
    }

    /* ---- Phasenaufgeloestes Timing ---
       max=Flaschenhals, avg=Mittelwert ueber alle Ranks.
       max/avg >> 1 in einer Phase = Lastungleichgewicht dort.
       "comm" = halo+migrate = gesamte Kommunikationszeit. */
    if(TIMING){
        double loc6[6]={T_force,T_cells,T_hash,T_halo,T_migrate,T_energy};
        double mx[6],sm[6];
        MPI_Reduce(loc6,mx,6,MPI_DOUBLE,MPI_MAX,0,CART);
        MPI_Reduce(loc6,sm,6,MPI_DOUBLE,MPI_SUM,0,CART);
        if(RANK==0){
            const char*nm[6]={"force","cells","hash","halo","migrate","energy"};
            double commMx=mx[3]+mx[4];
            printf("# PHASE        max_s     avg_s   max/avg\n");
            for(int i=0;i<6;i++){
                double avg=sm[i]/SIZE;
                printf("# %-9s %9.4f %9.4f %8.2f\n",nm[i],mx[i],avg,
                       avg>0?mx[i]/avg:1.0);
            }
            /* Maschinenlesbare Zeile fuer analyze_scaling.py */
            printf("# CSV ranks,walltime,force,cells,hash,halo,migrate,energy,comm\n");
            printf("# DATA %d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                   SIZE,t_max,mx[0],mx[1],mx[2],mx[3],mx[4],mx[5],commMx);
        }
    }
    MPI_Finalize();
    return 0;
}
