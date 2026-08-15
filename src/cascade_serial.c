/* ============================================================================
 *  cascade_serial.c  --  Serieller Prototyp: Kollisionskaskade im 2D-Gitter
 * ----------------------------------------------------------------------------
 *
 *  WAS MACHT DIESES PROGRAMM?
 *  --------------------------
 *  Dieses Programm simuliert eine KOLLISIONSKASKADE in einem zweidimensionalen
 *  Dreiecksgitter.  Man stelle sich ein gleichmaessiges Netz aus Atomen vor, die
 *  ueber Federn miteinander verbunden sind.  Ein einzelnes Atom (das PKA,
 *  "Primary Knock-on Atom") erhaelt einen starken Anstoss.  Es stoesst andere
 *  Atome an, die wiederum andere anstossen -- aehnlich wie beim Billard, aber in
 *  einem dichten Gitter.  Federn, die dabei zu weit gedehnt werden, reissen ab.
 *  Am Ende entstehen Cluster aus gerissenen Bindungen -- das sind die
 *  simulierten Strahlenschaeden.
 *
 *  WARUM 2D-DREIECKSGITTER?
 *  ------------------------
 *  Ein Dreiecksgitter hat 6 Nachbarn pro Innenatom.  Das macht es:
 *   - Schubstabil: Ein Quadratgitter mit reinen Axialfedern hat KEINE
 *     Schersteifigkeit und wuerde bei seitlicher Kraft kollabieren.
 *   - Elastisch nahezu isotrop: gleiche Steifigkeit in alle Richtungen.
 *  Einfachstes Gitter mit realistischem mechanischen Verhalten.
 *
 *  ZWEI KRAEFTE IM MODELL
 *  ----------------------
 *  (1) HARMONISCHE FEDERN zwischen gebundenen Nachbarn (Distanz L0 im Gleichgewicht):
 *         F_spring = K_SPRING * (r - L0) / r   [Richtungsvektor inklusive]
 *      Wenn die Feder zu weit gedehnt wird (r > L0 * MAX_STRETCH), reisst sie.
 *      Die dabei verlorene Energie wird in E_broken gebucht (Energieerhaltung!).
 *
 *  (2) KURZREICHWEITS-ABSTOSSUNG zwischen NICHT-gebundenen Atompaaren:
 *         F_rep(r) = K_REP * ((RCUT/r)^REP_N - 1)  fuer r < RCUT
 *      Normiert auf F(RCUT) = 0.  Gibt Atomen einen effektiven "Radius",
 *      damit das PKA beim Einschlag tatsaechlich auf andere Atome trifft.
 *      Ohne diese Kraft wuerden sich nicht-gebundene Atome einfach durchdringen.
 *
 *  ENERGIEERHALTUNG  (Korrektheitsbeweis)
 *  ----------------------------------------
 *  Mit DAMPING=0 muss die Gesamtenergie KONSTANT bleiben:
 *       E_total = E_kin + E_spring + E_rep + E_broken = const
 *  Das ist der wichtigste Test!
 *
 *  STRUKTUR-OF-ARRAYS (SoA) DATENLAYOUT
 *  ----------------------------------------
 *  Statt eines Arrays aus Structs (AoS: [{x,y,vx,vy}, {x,y,vx,vy}, ...])
 *  speichert SoA alle x-Koordinaten zusammen, alle y-Koordinaten zusammen usw.:
 *       double px[N], py[N], vx[N], vy[N], fx[N], fy[N];
 *  Vorteil: Schleifen, die nur px[] benoetigen, lesen einen zusammenhaengenden
 *  Speicherblock -> SIMD-Vektorisierung und CPU-Cache werden optimal genutzt.
 *  Dieses Layout ist auch GPU-kompatibel (coalesced memory access).
 *
 *  VELOCITY-VERLET INTEGRATOR  (kick-drift-kick)
 *  -----------------------------------------------
 *  Symplektischer Zeitintegrator mit guter Energieerhaltung.  Jeder Zeitschritt:
 *    1. Halber Kick:   v(t+dt/2) = v(t)     + (F/m) * dt/2
 *    2. Drift:         x(t+dt)   = x(t)     + v(t+dt/2) * dt
 *    3. Kraefte:       F(t+dt)   = F(x(t+dt))
 *    4. Halber Kick:   v(t+dt)   = v(t+dt/2) + (F/m) * dt/2
 *  Besser als Euler: kein kumulativer Energiedrift.
 *
 *  Build:  gcc -O2 -o cascade_serial cascade_serial.c -lm
 *  Run:    ./cascade_serial [params.ini]   (Standard: params.ini im CWD)
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ============================================================================
 *  PHYSIKALISCHE PARAMETER
 * ---------------------------------------------------------------------------
 *  Alle Groessen sind DIMENSIONSLOS:
 *    L0   = 1  (Gitterabstand = Laengeneinheit)
 *    MASS = 1  (Atommasse = Masseneinheit)
 *  Abgeleitet:
 *    Schallgeschwindigkeit:  c = sqrt(K_SPRING) * L0 = 10
 *    Zeiteinheit:            t_u = L0 * sqrt(MASS/K_SPRING) = 0.1
 *    dt-Stabilitaet:         dt < 2/omega_max = 0.2  (Sicherheitsfaktor ~100x)
 *
 *  Diese Werte werden beim Programmstart aus params.ini ueberschrieben.
 * ========================================================================== */
static double L0          = 1.0;     /* Ruhelaenge der Federn / Gitterabstand      */
static double MASS        = 1.0;     /* Masse jedes Atoms (alle gleich schwer)     */
static double K_SPRING    = 100.0;   /* Federkonstante (Haerte der Bindungen)      */
static double MAX_STRETCH = 1.15;    /* Bruchbedingung: Feder reisst bei r > L0*1.15 */
static double RCUT        = 0.9;     /* Abstoss-Reichweite (MUSS kleiner als L0 sein!)*/
static double K_REP       = 400.0;   /* Amplitude der repulsiven Kraft             */
static double REP_N       = 12.0;    /* Steilheitsexponent (n=12: sehr harter Kern) */
static double V_THRESH    = 1e30;    /* Schwelle fuer el. Stopping (1e30 = inaktiv) */

/* ============================================================================
 *  GITTERGEOMETRIE
 *  Das Dreiecksgitter hat NX*NY Atome in einem rechteckigen Ausschnitt.
 *  DY = L0 * sqrt(3)/2 ist der Abstand zwischen den Zeilen (nicht L0!),
 *  weil die Zeilen versetzt sind und damit Dreiecke formen.
 * ========================================================================== */
static int    NX = 120, NY = 120;    /* Spalten (x) und Zeilen (y) des Gitters    */
static double DY;                    /* Zeilenabstand = L0*sqrt(3)/2 (berechnet)  */

/* ============================================================================
 *  TEILCHENFELDER  --  Structure-of-Arrays (SoA) Layout
 * ---------------------------------------------------------------------------
 *  Alle N Werte einer physikalischen Groesse liegen zusammenhaengend im Speicher.
 *  Beispiel fuer N=4 Atome:
 *    px = [x0, x1, x2, x3]   (zusammenhaengend im RAM)
 *    py = [y0, y1, y2, y3]
 *    ...
 *  NICHT:  atoms = [{x0,y0,vx0,...}, {x1,y1,vx1,...}]  (das waere AoS)
 * ========================================================================== */
static int     N;                    /* Gesamtzahl der Atome = NX * NY            */
static double *px, *py;              /* Aktuelle x/y-Positionen                   */
static double *vx, *vy;              /* Geschwindigkeitskomponenten               */
static double *fx, *fy;              /* Kraeftekomponenten (jeden Schritt neu)    */
static double *mass;                 /* Massen (normalerweise alle = MASS)        */
static double *px0, *py0;           /* Anfangspositionen (fuer Verschiebungsanalyse)*/

/* ============================================================================
 *  BINDUNGSLISTE
 * ---------------------------------------------------------------------------
 *  Jede Feder verbindet genau zwei Atome.  Wir speichern sie als Index-Paare:
 *    bond_a[b], bond_b[b]    = die beiden Atome der b-ten Bindung
 *    bond_intact[b]          = 1 (intakt) oder 0 (dauerhaft gerissen)
 *  Jede Bindung wird genau EINMAL gespeichert (nicht doppelt in beide Richtungen).
 * ========================================================================== */
static int    n_bonds;               /* Anzahl der Bindungen gesamt               */
static int   *bond_a, *bond_b;      /* Atom-Index-Paare (bond_a[b] -- bond_b[b]) */
static int   *bond_intact;           /* 1 = Feder intakt, 0 = Feder gerissen      */

/* ============================================================================
 *  NACHBARLISTE  (fuer den Abstoss-Ausschluss)
 * ---------------------------------------------------------------------------
 *  Wenn zwei Atome durch eine intakte Feder verbunden sind, soll die
 *  Kurzreichweits-Abstossung NICHT zusaetzlich wirken (die Feder uebernimmt).
 *  Um das effizient zu pruefen, haelt jedes Atom eine Liste seiner aktuellen
 *  Bindungspartner.
 *
 *  Layout: nbr[p*MAXNBR + k] = k-ter gebundener Nachbar von Atom p
 *          n_nbr[p]           = aktuelle Anzahl Nachbarn von Atom p (0..6)
 *
 *  MAXNBR=6: im Dreiecksgitter hat ein Innenatom genau 6 Nachbarn.
 *  Randata haben weniger.  Die Liste schrumpft wenn Federn reissen.
 * ========================================================================== */
#define MAXNBR 6
static int   *nbr;                   /* [N*MAXNBR] gepackte Nachbar-Indizes       */
static int   *n_nbr;                 /* [N] aktuelle Nachbarzahl je Atom          */
static int   *nbr0;                  /* [N] Anfangs-Nachbarzahl (fuer Bruchdiagnose)*/

/* ============================================================================
 *  CELL LIST  (O(N) Nachbarsuche)
 * ---------------------------------------------------------------------------
 *  Problem: Die Abstossung wirkt nur bei r < RCUT.  Ohne Hilfsdatenstruktur
 *  muesste man alle N*(N-1)/2 Paare pruefen -> O(N^2) pro Schritt.
 *
 *  Loesung: Das Simulationsgebiet wird in Zellen der Groesse RCUT eingeteilt.
 *  Jedes Atom wird in seine Zelle eingetragen.  Fuer die Abstossung von Atom p
 *  reicht es, nur die 9 umliegenden Zellen zu durchsuchen, weil kein Atom
 *  weiter als RCUT entfernt sein kann, wenn es in einer benachbarten Zelle ist.
 *  -> Im Schnitt ~10-20 Nachbaratome statt N -> O(N) pro Schritt.
 *
 *  Implementation als verkettete Liste:
 *    cell_head[c] = erster Atom-Index in Zelle c (-1 wenn leer)
 *    cell_next[p] = naechstes Atom nach p in derselben Zelle (-1 am Ende)
 *  Eintragen: cell_next[p]=cell_head[c]; cell_head[c]=p  (Prepend, O(1))
 * ========================================================================== */
static int    ncx, ncy;              /* Anzahl Zellen in x/y-Richtung            */
static double cell_size;             /* Zellengroesse = RCUT                      */
static double dom_x0, dom_y0;        /* Koordinatenursprung des Zellgitters       */
static int   *cell_head;             /* [ncx*ncy] Kopf der Atom-Kette je Zelle   */
static int   *cell_next;             /* [N] Zeiger auf naechstes Atom in Kette   */

/* ============================================================================
 *  ENERGIEBUCHHALTUNG
 * ---------------------------------------------------------------------------
 *  E_broken: Wenn eine Feder reisst, geht ihre potenzielle Energie nicht
 *  einfach "verloren".  Sie wird in E_broken akkumuliert.  Damit gilt:
 *     E_total = E_kin + E_spring + E_rep + E_broken = const  (ohne Daempfung)
 *  Das ist der physikalisch korrekte Erhaltungssatz!
 * ========================================================================== */
static double E_broken = 0.0;        /* Kumulierte verlorene Bindungsenergie      */
static double E_damp   = 0.0;        /* Kumulierte Daempfungsarbeit (derzeit unused)*/
static double DAMPING  = 0.0;        /* Daempfungskoeffizient (0 = keine Daempfung)*/

/* ============================================================================
 *  HEALING  (Frenkel-Paar-Rekombination)
 * ---------------------------------------------------------------------------
 *  In echten Metallen koennen Gitterdefekte (Frenkel-Paare: Vakanz + Interstitial)
 *  wieder rekombinieren.  Dieses Modell erlaubt gerissene Federn, sich wieder
 *  zu schliessen, wenn die Atome nah genug zusammenkommen und langsam genug sind.
 * ========================================================================== */
static int    HEALING       = 0;     /* 0 = Healing aus, 1 = Healing an           */
static double HEALING_DIST  = 1.10;  /* Rekombination wenn r < HEALING_DIST * L0  */
static double HEALING_VREL  = 0.5;   /* ...und Relativgeschwindigkeit < Schwelle   */

/* ============================================================================
 *  LAUFZEIT-KONFIGURATION  (wird aus params.ini ueberschrieben)
 * ========================================================================== */
static int      NSTEPS     = 12000;  /* Anzahl Zeitschritte                       */
static double   DT         = 1e-3;   /* Zeitschrittweite                          */
static unsigned SEED       = 1;      /* Seed fuer reproduzierbare RNG-Sequenz     */
static int      LOG_EVERY  = 0;      /* 0 -> auto (NSTEPS/40), >0 -> explizit     */
static int      DUMP_EVERY = 0;      /* 0 -> kein XYZ-Dump, >0 -> alle N Schritte */
static int      REORDER    = 0;      /* 1 -> Morton-Sortierung (Cache-Optimierung) */
static double   PKA_X      = 0.5;    /* Rel. x-Startposition des PKA (0..1)      */
static double   PKA_Y      = 0.5;    /* Rel. y-Startposition des PKA (0..1)      */
static double   PKA_ENERGY = 1000.0; /* Kinetische Startenergie des PKA           */
static double   PKA_ANGLE  = 15.0;   /* Abschusswinkel in Grad (0=oben, 90=rechts)*/
static double   PKA_MASS   = 1.0;    /* Masse des PKA (kann von MASS abweichen)   */
static double   ABSORB_BORDER = 0.0; /* Absorbierende Randzone (Platzhalter)      */
static char     OUT_PREFIX[128] = "run"; /* Dateinamen-Praefix fuer CSV-Ausgaben  */
static char     MODELSEL[8]     = "A";   /* Modell-Auswahl (nur "A" implementiert) */

/* ============================================================================
 *  MULTI-PKA UNTERSTUETZUNG
 * ---------------------------------------------------------------------------
 *  Standardmaessig gibt es nur einen PKA.  Man kann mehrere starten indem man
 *  n_pka > 1 setzt (automatische Verteilung) oder [pka1],[pka2],... Sektionen
 *  in der ini-Datei definiert (manuell, mit individuellen Parametern).
 *
 *  PkaDef.hx/.hy/.he/.ha/.hm sind "has"-Flags:
 *    1 = dieses Feld wurde explizit in [pkaN] angegeben
 *    0 = aus dem globalen [pka]-Block erben
 * ========================================================================== */
static int      N_PKA = 1;            /* Anzahl gleichzeitiger PKAs               */
#define MAX_PKA 64
typedef struct {
    double x, y, energy, angle, mass;  /* PKA-Parameter                          */
    int hx, hy, he, ha, hm;            /* has-Flags: wurde angegeben? (1=ja)      */
} PkaDef;
static PkaDef   MANUAL_PKA[MAX_PKA];  /* Array manuell definierter PKAs           */
static int      N_MANUAL = 0;         /* Anzahl definierter [pkaN]-Sektionen      */

/* ============================================================================
 *  DETERMINISTISCHER ZUFALLSZAHLENGENERATOR  (LCG, 64-Bit)
 * ---------------------------------------------------------------------------
 *  Liefert bei gleichem SEED immer dieselbe Sequenz -> reproduzierbare Laeufe.
 *  Wichtig: IDENTISCH zur MPI-Version, damit auto-verteilte PKAs auf beiden
 *  Implementierungen exakt dieselben Positionen haben (Vergleichbarkeit).
 *
 *  Algorithmus: Linear Congruential Generator (LCG)
 *    RNG_S = a * RNG_S + c  (mod 2^64)
 *  mit den Knuth-Konstanten a=6364136223846793005, c=1442695040888963407.
 *  Gibt eine Gleichverteilung in [0,1) als double zurueck.
 * ========================================================================== */
static unsigned long RNG_S;
static double rng_u(void){
    RNG_S = RNG_S * 6364136223846793005UL + 1442695040888963407UL;
    /* Die oberen 52 Bit des Zustands bilden die Mantisse eines double in [0,1) */
    return ((RNG_S >> 11) & ((1UL << 52) - 1)) / (double)(1UL << 52);
}

/* ============================================================================
 *  CONFIG-PARSER  (liest params.ini)
 * ---------------------------------------------------------------------------
 *  Liest Zeile fuer Zeile.  Format:  schluessel = wert
 *  Kommentare (# und ;) werden abgeschnitten.
 *  Sektions-Header [name] werden erkannt um [pka1],[pka2],...-Blocks zu routen.
 *  Unbekannte Schluessel werden still ignoriert (Vorwaertskompatibilitaet).
 * ========================================================================== */
static void read_config(const char*fn){
    FILE*f=fopen(fn,"r");
    if(!f){ fprintf(stderr,"WARN: Config '%s' nicht gefunden -> Defaults\n",fn); return; }
    char line[256];
    char section[64]="";
    while(fgets(line,sizeof line,f)){
        char*h=strpbrk(line,";#"); if(h)*h=0;        /* Kommentar abschneiden */
        char*lb=strchr(line,'[');                    /* Sektions-Header?      */
        if(lb){ section[0]=0; sscanf(lb,"[%63[^]]]",section); continue; }
        char*eq=strchr(line,'='); if(!eq) continue;  /* keine Zuweisung       */
        *eq=0;
        char key[64], val[128];
        if(sscanf(line,"%63s",key)!=1) continue;
        if(sscanf(eq+1,"%127s",val)!=1) continue;
        #define KV(k) (!strcmp(key,k))
        /* [pka1],[pka2],... in das Manual-Array routen:
         * sscanf prueft ob section="pka1","pka2",... und liest den Index. */
        { int pidx;
          if(sscanf(section,"pka%d",&pidx)==1 && pidx>=1 && pidx<=MAX_PKA){
            PkaDef*p=&MANUAL_PKA[pidx-1]; if(pidx>N_MANUAL)N_MANUAL=pidx;
            if      (KV("pka_x"))     {p->x=atof(val);     p->hx=1;}
            else if (KV("pka_y"))     {p->y=atof(val);     p->hy=1;}
            else if (KV("pka_energy")){p->energy=atof(val);p->he=1;}
            else if (KV("pka_angle")) {p->angle=atof(val); p->ha=1;}
            else if (KV("pka_mass"))  {p->mass=atof(val);  p->hm=1;}
            continue;
          } }
        /* Alle globalen Parameter */
        if      (KV("model"))        strncpy(MODELSEL,val,7);
        else if (KV("healing"))      HEALING=(!strcmp(val,"true")||!strcmp(val,"1"));
        else if (KV("healing_dist")) HEALING_DIST=atof(val);
        else if (KV("healing_vrel")) HEALING_VREL=atof(val);
        else if (KV("NX"))           NX=atoi(val);
        else if (KV("NY"))           NY=atoi(val);
        else if (KV("L0"))           L0=atof(val);
        else if (KV("K_SPRING"))     K_SPRING=atof(val);
        else if (KV("MAX_STRETCH"))  MAX_STRETCH=atof(val);
        else if (KV("RCUT"))         RCUT=atof(val);
        else if (KV("K_REP"))        K_REP=atof(val);
        else if (KV("REP_N"))        REP_N=atof(val);
        else if (KV("pka_x"))        PKA_X=atof(val);
        else if (KV("pka_y"))        PKA_Y=atof(val);
        else if (KV("pka_energy"))   PKA_ENERGY=atof(val);
        else if (KV("pka_angle"))    PKA_ANGLE=atof(val);
        else if (KV("pka_mass"))     PKA_MASS=atof(val);
        else if (KV("n_pka"))        N_PKA=atoi(val);
        else if (KV("dt"))           DT=atof(val);
        else if (KV("n_steps"))      NSTEPS=atoi(val);
        else if (KV("damping"))      DAMPING=atof(val);
        else if (KV("v_thresh"))     V_THRESH=atof(val);
        else if (KV("absorb_border"))ABSORB_BORDER=atof(val);
        else if (KV("log_every"))    LOG_EVERY=atoi(val);
        else if (KV("dump_every"))   DUMP_EVERY=atoi(val);
        else if (KV("reorder"))      REORDER=(!strcmp(val,"true")||!strcmp(val,"1"));
        else if (KV("out_prefix"))   strncpy(OUT_PREFIX,val,127);
        else if (KV("seed"))         SEED=(unsigned)atoi(val);
        /* unbekannte Schluessel (z.B. "lattice", "enable_vtk") still ignoriert */
        #undef KV
    }
    fclose(f);
}

/* ============================================================================
 *  GITTERAUFBAU
 * ---------------------------------------------------------------------------
 *  Erzeugt das Dreiecksgitter und alle Bindungen.
 *
 *  Dreiecksgitter-Geometrie:
 *    Gerade Zeilen   (j=0,2,4,...): Atome bei x = i*L0
 *    Ungerade Zeilen (j=1,3,5,...): Atome bei x = i*L0 + L0/2  (nach rechts versetzt)
 *    Alle Zeilen:                   y = j * L0*sqrt(3)/2
 *
 *  Dadurch hat jedes Innenatom genau 6 gleich weit entfernte Nachbarn.
 *
 *  Bindungen: Jede Bindung wird NUR EINMAL gespeichert.
 *  Wir erzeugen fuer Atom (j,i) nur Bindungen NACH RECHTS und NACH OBEN,
 *  damit keine Doppelzaehlung passiert.
 * ========================================================================== */
static int idx_of(int j, int i){ return j*NX + i; }  /* linearisierter Gitterindex */

/* Traegt a und b als gegenseitige Nachbarn ein (fuer Abstoss-Ausschluss) */
static void add_neighbor(int a, int b){
    nbr[a*MAXNBR + n_nbr[a]++] = b;
    nbr[b*MAXNBR + n_nbr[b]++] = a;
}

/* ============================================================================
 *  PKA-PLATZIERUNG
 * ---------------------------------------------------------------------------
 *  place_one_pka(): Setzt einem Atom die PKA-Masse und Anfangsgeschwindigkeit.
 *    (xr, yr) = relative Koordinaten (0..1) -> werden in Gitterindizes umgerechnet
 *    energy   = kinetische Startenergie  -> v = sqrt(2*E/m)
 *    ang      = Abschusswinkel in Grad (0=oben/neg.y, 90=rechts/pos.x)
 *    Formel:  vx = v*sin(ang),  vy = -v*cos(ang)
 *    (Negatives vy, weil y im Gitter nach oben zeigt aber 0 Grad = nach oben)
 * ========================================================================== */
static void place_one_pka(double xr, double yr, double energy, double ang, double m){
    /* Relative Koordinaten (0..1) in Gitter-Zeilenindex/Spaltenindex umrechnen */
    int pj=(int)(yr*(NY-1)+0.5), pi=(int)(xr*(NX-1)+0.5);
    /* Sicherstellen dass der Index im Gitter bleibt */
    if(pj<0)pj=0; if(pj>=NY)pj=NY-1; if(pi<0)pi=0; if(pi>=NX)pi=NX-1;
    int p=idx_of(pj,pi);
    mass[p]=m;
    if(energy>0.0){
        double v=sqrt(2.0*energy/mass[p]);  /* E_kin = 0.5*m*v^2 -> v = sqrt(2E/m) */
        double th=ang*M_PI/180.0;           /* Grad -> Bogenmass fuer sin/cos        */
        vx[p]=v*sin(th); vy[p]=-v*cos(th);
    }
}

/* Setzt ALLE PKAs -- manueller Modus ([pkaN] vorhanden) hat Vorrang.
 * Auto-Modus: N_PKA PKAs deterministisch verteilt (identisch zur MPI-Version).
 * Rueckgabe: tatsaechlich gesetzte Anzahl. */
static int place_pkas(void){
    if(N_MANUAL>0){                   /* --- MANUELLER MODUS: aus [pka1],[pka2]... ---*/
        int placed=0;
        for(int k=0;k<N_MANUAL && placed<N_PKA;k++){
            PkaDef*p=&MANUAL_PKA[k];
            /* Leere Luecke (kein Feld angegeben) ueberspringen */
            if(!(p->hx||p->hy||p->he||p->ha||p->hm)) continue;
            /* Fehlende Felder aus globalen Standardwerten erben */
            place_one_pka(p->hx?p->x:PKA_X, p->hy?p->y:PKA_Y,
                          p->he?p->energy:PKA_ENERGY,
                          p->ha?p->angle:PKA_ANGLE,
                          p->hm?p->mass:PKA_MASS);
            placed++;
        }
        return placed;
    }
    /* --- AUTO-MODUS: N_PKA PKAs deterministisch verteilen --- */
    RNG_S=(unsigned long)SEED*2654435761UL+12345UL;  /* RNG initialisieren      */
    for(int kp=0;kp<N_PKA;kp++){
        double xr, yr, ang;
        if(kp==0 && N_PKA==1){
            /* Erster (einziger) PKA: genau wie in params.ini angegeben */
            xr=PKA_X; yr=PKA_Y; ang=PKA_ANGLE;
        } else {
            /* Weitere PKAs: zufaellig im mittleren Bereich (Rand vermeiden) */
            xr=0.15+0.70*rng_u(); yr=0.10+0.80*rng_u(); ang=360.0*rng_u();
        }
        place_one_pka(xr,yr,PKA_ENERGY,ang,PKA_MASS);
    }
    return N_PKA;
}

/* Speicher allokieren + Dreiecksgitter aufbauen */
static void build_lattice(void){
    DY = L0 * sqrt(3.0)/2.0;   /* Zeilenabstand: Hoehe eines gleichseitigen Dreiecks */
    N  = NX*NY;                 /* Gesamtzahl der Atome */

    /* Speicher fuer alle per-Atom-Arrays allokieren */
    px=malloc(N*sizeof*px); py=malloc(N*sizeof*py);
    vx=calloc(N,sizeof*vx); vy=calloc(N,sizeof*vy);  /* calloc: mit 0 initialisiert */
    fx=calloc(N,sizeof*fx); fy=calloc(N,sizeof*fy);
    px0=malloc(N*sizeof*px0); py0=malloc(N*sizeof*py0);
    mass=malloc(N*sizeof*mass);
    nbr=malloc((size_t)N*MAXNBR*sizeof*nbr);
    n_nbr=calloc(N,sizeof*n_nbr);

    /* Atome auf dem Dreiecksgitter positionieren */
    for(int j=0;j<NY;j++)
        for(int i=0;i<NX;i++){
            int id=idx_of(j,i);
            /* (j&1) ist 1 fuer ungerade Zeilen -> Versatz um L0/2 nach rechts */
            px[id]= i*L0 + (j&1)*(L0*0.5);
            py[id]= j*DY;
            px0[id]=px[id]; py0[id]=py[id];  /* Anfangspositionen merken */
            mass[id]=MASS;
        }

    /* Bindungen erzeugen: pro Atom nur RECHT und OBEN -> jede Feder genau 1x.
     * Maximal 3 Bindungen pro Atom in diese Richtungen -> cap=3*N genuegt. */
    int cap = 3*N;
    bond_a=malloc(cap*sizeof*bond_a);
    bond_b=malloc(cap*sizeof*bond_b);
    bond_intact=malloc(cap*sizeof*bond_intact);
    n_bonds=0;

    for(int j=0;j<NY;j++)
        for(int i=0;i<NX;i++){
            int a=idx_of(j,i);
            /* Bindung nach rechts: (j,i) -- (j,i+1) */
            if(i+1<NX){ int b=idx_of(j,i+1);
                bond_a[n_bonds]=a; bond_b[n_bonds]=b; bond_intact[n_bonds]=1;
                add_neighbor(a,b); n_bonds++; }
            /* Bindungen nach oben (zwei, je nach gerader/ungerader Zeile) */
            if(j+1<NY){
                int li, ri;
                /* gerade Zeile j: obere Nachbarn links=(j+1,i-1), rechts=(j+1,i)  */
                /* ungerade Zeile: obere Nachbarn links=(j+1,i),   rechts=(j+1,i+1)*/
                if((j&1)==0){ li=i-1; ri=i; } else { li=i; ri=i+1; }
                if(li>=0 && li<NX){ int b=idx_of(j+1,li);
                    bond_a[n_bonds]=a; bond_b[n_bonds]=b; bond_intact[n_bonds]=1;
                    add_neighbor(a,b); n_bonds++; }
                if(ri>=0 && ri<NX){ int b=idx_of(j+1,ri);
                    bond_a[n_bonds]=a; bond_b[n_bonds]=b; bond_intact[n_bonds]=1;
                    add_neighbor(a,b); n_bonds++; }
            }
        }
}

/* Prueft ob Atom p und q durch eine intakte Bindung verbunden sind.
 * Lineare Suche in der Nachbarliste (max. 6 Eintraege -> sehr schnell). */
static int is_bonded(int p, int q){
    int base=p*MAXNBR, n=n_nbr[p];
    for(int k=0;k<n;k++) if(nbr[base+k]==q) return 1;
    return 0;
}

/* Entfernt q aus der Nachbarliste von p.
 * "Swap-with-last": letztes Element ueberschreibt das geloeschte -> O(1). */
static void drop_neighbor(int p, int q){
    int base=p*MAXNBR, n=n_nbr[p];
    for(int k=0;k<n;k++) if(nbr[base+k]==q){
        nbr[base+k]=nbr[base+n-1];  /* letzten Eintrag nach vorne kopieren */
        n_nbr[p]--;
        return;
    }
}

/* ============================================================================
 *  CELL LIST AUFBAU
 * ---------------------------------------------------------------------------
 *  cell_setup(): Einmalig die Gitterdimensionen und Speicher berechnen.
 *  cell_build(): Jeden Zeitschritt: alle Atome in ihre aktuelle Zelle eintragen.
 *
 *  Grosszuegiger Rand (5*L0): Das PKA schlaegt Atome aus dem Gitterbereich
 *  heraus.  Ohne Rand wuerden diese Atome ihre Zelle nicht finden.
 * ========================================================================== */
static void cell_setup(void){
    /* Bounding Box der Atome bestimmen */
    double x0=1e30,y0=1e30,x1=-1e30,y1=-1e30;
    for(int p=0;p<N;p++){
        if(px[p]<x0)x0=px[p]; if(px[p]>x1)x1=px[p];
        if(py[p]<y0)y0=py[p]; if(py[p]>y1)y1=py[p];
    }
    double margin = 5.0*L0;  /* Rand: Atome koennen herausgeschlagen werden */
    dom_x0=x0-margin; dom_y0=y0-margin;
    double w=(x1-x0)+2*margin, h=(y1-y0)+2*margin;
    cell_size = RCUT;  /* Zellenbreite = Cutoff: Nachbarn nur in den 9 Nachbarzellen */
    ncx=(int)(w/cell_size)+1; ncy=(int)(h/cell_size)+1;
    cell_head=malloc((size_t)ncx*ncy*sizeof*cell_head);
    cell_next=malloc((size_t)N*sizeof*cell_next);
}

/* Berechnet den Zell-Index fuer Position (x,y). Klampert auf gueltige Grenzen. */
static inline int cell_of(double x, double y){
    int cx=(int)((x-dom_x0)/cell_size);
    int cy=(int)((y-dom_y0)/cell_size);
    if(cx<0)cx=0; if(cx>=ncx)cx=ncx-1;  /* Ausreisser klampen */
    if(cy<0)cy=0; if(cy>=ncy)cy=ncy-1;
    return cy*ncx+cx;
}

/* Baut die verkettete Cell-List neu auf (O(N)).
 * Alle Zellen auf leer (-1) setzen, dann jedes Atom als Listenkopf eintragen. */
static void cell_build(void){
    for(int c=0;c<ncx*ncy;c++) cell_head[c]=-1;   /* alle Zellen leeren */
    for(int p=0;p<N;p++){
        int c=cell_of(px[p],py[p]);
        /* Prepend in verkettete Liste: p zeigt auf alten Kopf, wird neuer Kopf */
        cell_next[p]=cell_head[c]; cell_head[c]=p;
    }
}

/* ============================================================================
 *  ABSTOSS-KRAFT UND -POTENTIAL  (isoliert, leicht austauschbar)
 * ---------------------------------------------------------------------------
 *  Diese beiden Funktionen sind bewusst herausgezogen: man kann spaeter
 *  Morse- oder Lennard-Jones-Potenzial als Drop-in-Ersatz einsetzen.
 *
 *  rep_force(r): Kraftbetrag (positiv = abstossend) fuer Abstand r < RCUT.
 *    F(r) = K_REP * ((RCUT/r)^REP_N - 1)
 *    F(RCUT) = K_REP * (1-1) = 0  -> kraftstetig am Cutoff!
 *    F->inf fuer r->0              -> harter Kern.
 *
 *  rep_pot(r): Potenzial mit V(RCUT)=0 (Stammfunktion von -F, normiert).
 *    Wird nur fuer die Energiebuchhaltung benoetigt, NICHT fuer die Kraft.
 * ========================================================================== */
static inline double rep_force(double r){
    return K_REP * ( pow(RCUT/r, REP_N) - 1.0 );
}
static inline double rep_pot(double r){
    double n=REP_N;
    return K_REP * ( (RCUT - pow(RCUT,n)*pow(r,1.0-n))/(1.0-n) - (RCUT - r) );
}

/* ============================================================================
 *  KRAFTBERECHNUNG  (Herzstuck der Simulation, ~90% der Laufzeit)
 * ---------------------------------------------------------------------------
 *  Berechnet fuer jedes Atom die Kraft aus drei Beitraegen:
 *  (1) Federn:    Harmonische Kraft entlang intakter Bindungen + Bruchdetektion
 *  (2) Abstossung: Kurzreichweits-Kraft zwischen nicht-gebundenen Paaren (Cell List)
 *  (3) Daempfung: Geschwindigkeitsproportionale Bremskraft (nur wenn DAMPING>0)
 *
 *  "Gather-Formulation" fuer die Abstossung:
 *  Jedes Paar (p,q) wird genau einmal betrachtet.  Beide Atome bekommen
 *  ihre Kraft direkt zugewiesen (Newton 3: F_pq = -F_qp).
 *  Das funktioniert serial, ist aber nicht direkt GPU-parallelisierbar
 *  (wegen Schreibkonflikten auf fx[q]).  Die CUDA-Version nutzt stattdessen
 *  "Gather": jedes Atom p summiert nur seine eigene Kraft, ohne in fx[q] zu
 *  schreiben -- sicher parallelisierbar ohne atomicAdd.
 * ========================================================================== */
static void compute_forces(void){
    memset(fx,0,N*sizeof*fx);  /* alle Kraefte auf 0 zuruecksetzen */
    memset(fy,0,N*sizeof*fy);

    /* ---------------------------------------------------------------------- */
    /* (1) FEDERKRAEFTE + BRUCHDETEKTION                                       */
    /* ---------------------------------------------------------------------- */
    for(int b=0;b<n_bonds;b++){
        if(!bond_intact[b]) continue;  /* gerissene Feder ueberspringen */
        int a=bond_a[b], c=bond_b[b];
        double dx=px[c]-px[a], dy=py[c]-py[a];
        double r=sqrt(dx*dx+dy*dy);

        /* Bruchbedingung: Feder zu weit gedehnt? */
        if(r > L0*MAX_STRETCH){
            bond_intact[b]=0;                           /* dauerhaft reissen */
            E_broken += 0.5*K_SPRING*(r-L0)*(r-L0);    /* Energie buchen */
            drop_neighbor(a,c); drop_neighbor(c,a);     /* Nachbarlisten aktualisieren */
            continue;
        }

        /* Harmonische Federkraft: F = K*(r-L0), gerichtet entlang (dx,dy)
         * Effiziente Form: f = K*(r-L0)/r, dann Fx = f*dx, Fy = f*dy */
        double f = K_SPRING*(r-L0)/r;
        double fxc=f*dx, fyc=f*dy;
        fx[a]+=fxc; fy[a]+=fyc;  /* Atom a: wird von c angezogen/abstossen */
        fx[c]-=fxc; fy[c]-=fyc;  /* Atom c: entgegengesetzte Kraft (Newton 3) */
    }

    /* ---------------------------------------------------------------------- */
    /* (2) ABSTOSS-KRAEFTE via Cell List (nur nicht-gebundene Paare)           */
    /* ---------------------------------------------------------------------- */
    double rc2=RCUT*RCUT;  /* Quadrat des Cutoffs (Vergleich ohne sqrt effizienter) */

    /* Halb-Stern der Nachbarzellen:
     * Wir durchsuchen nur die 4 "rechten" Nachbarzellen plus die gleiche Zelle.
     * Damit wird jedes Paar genau EINMAL betrachtet (kein Doppelzaehlen).
     * Offsets: (1,0)=rechts, (1,1)=rechts-oben, (0,1)=oben, (-1,1)=links-oben */
    static const int off[4][2]={{1,0},{1,1},{0,1},{-1,1}};
    for(int cy=0;cy<ncy;cy++)
    for(int cx=0;cx<ncx;cx++){
        int c=cy*ncx+cx;
        for(int p=cell_head[c]; p!=-1; p=cell_next[p]){
            /* Paare innerhalb der gleichen Zelle:
             * q iteriert nur ueber Atome die in der Liste NACH p kommen -> kein Doppelzaehlen */
            for(int q=cell_next[p]; q!=-1; q=cell_next[q]){
                double dx=px[p]-px[q], dy=py[p]-py[q];
                double r2=dx*dx+dy*dy;
                if(r2>=rc2||r2==0.0) continue;  /* zu weit oder gleich */
                if(is_bonded(p,q)) continue;     /* gebunden -> Feder wirkt, keine Abstossung */
                double r=sqrt(r2), f=rep_force(r)/r;
                fx[p]+=f*dx; fy[p]+=f*dy;   /* p: von q weggedruckt */
                fx[q]-=f*dx; fy[q]-=f*dy;   /* q: von p weggedruckt (Newton 3) */
            }
            /* Paare mit den 4 Nachbarzellen */
            for(int o=0;o<4;o++){
                int nxc=cx+off[o][0], nyc=cy+off[o][1];
                if(nxc<0||nxc>=ncx||nyc<0||nyc>=ncy) continue;
                int nc=nyc*ncx+nxc;
                for(int q=cell_head[nc]; q!=-1; q=cell_next[q]){
                    double dx=px[p]-px[q], dy=py[p]-py[q];
                    double r2=dx*dx+dy*dy;
                    if(r2>=rc2||r2==0.0) continue;
                    if(is_bonded(p,q)) continue;
                    double r=sqrt(r2), f=rep_force(r)/r;
                    fx[p]+=f*dx; fy[p]+=f*dy;
                    fx[q]-=f*dx; fy[q]-=f*dy;
                }
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* (3) DAEMPFUNG / ELECTRONIC STOPPING  (nur wenn DAMPING > 0)            */
    /* ---------------------------------------------------------------------- */
    /* Schnelle Atome verlieren Energie an das Elektronensystem (electronic stopping).
     * Modell: geschwindigkeitsproportionale Bremskraft F_damp = -gamma * v.
     * Schwelle V_THRESH: nur Atome oberhalb dieser Geschwindigkeit werden gedaempft,
     * damit nicht die thermischen Gitterschwingungen abgebremst werden. */
    if(DAMPING>0.0){
        for(int p=0;p<N;p++){
            double v=sqrt(vx[p]*vx[p]+vy[p]*vy[p]);
            if(v>V_THRESH){
                fx[p]-=DAMPING*vx[p];
                fy[p]-=DAMPING*vy[p];
            }
        }
    }
}

/* ============================================================================
 *  ENERGIEBERECHNUNG  (nur fuer Diagnosezwecke, alle LOG_EVERY Schritte)
 * ---------------------------------------------------------------------------
 *  Berechnet E_kin, E_spring, E_rep separat fuer die Ausgabe-Tabelle.
 *  Mit DAMPING=0 muss E_kin + E_spring + E_rep + E_broken = E0 (konstant) gelten.
 *  Erlaubte Abweichung: < 0.01% (Velocity-Verlet: symplektisch, kein Drift).
 * ========================================================================== */
static void energies(double*ek, double*es, double*er){
    double sk=0, ss=0, sr=0;

    /* Kinetische Energie: Summe ueber alle Atome von 0.5 * m * v^2 */
    for(int p=0;p<N;p++) sk+=0.5*mass[p]*(vx[p]*vx[p]+vy[p]*vy[p]);

    /* Federpotenzial: Summe ueber alle intakten Bindungen von 0.5*K*(r-L0)^2 */
    for(int b=0;b<n_bonds;b++){
        if(!bond_intact[b])continue;
        int a=bond_a[b],c=bond_b[b];
        double dx=px[c]-px[a],dy=py[c]-py[a],r=sqrt(dx*dx+dy*dy);
        ss+=0.5*K_SPRING*(r-L0)*(r-L0);
    }

    /* Abstoss-Potenzial: gleiche Cell-List-Topologie wie in compute_forces */
    double rc2=RCUT*RCUT;
    static const int off[4][2]={{1,0},{1,1},{0,1},{-1,1}};
    for(int cy=0;cy<ncy;cy++)for(int cx=0;cx<ncx;cx++){
        int c=cy*ncx+cx;
        for(int p=cell_head[c];p!=-1;p=cell_next[p]){
            for(int q=cell_next[p];q!=-1;q=cell_next[q]){
                double dx=px[p]-px[q],dy=py[p]-py[q],r2=dx*dx+dy*dy;
                if(r2>=rc2||r2==0.0)continue; if(is_bonded(p,q))continue;
                sr+=rep_pot(sqrt(r2));
            }
            for(int o=0;o<4;o++){int nxc=cx+off[o][0],nyc=cy+off[o][1];
                if(nxc<0||nxc>=ncx||nyc<0||nyc>=ncy)continue;
                int nc=nyc*ncx+nxc;
                for(int q=cell_head[nc];q!=-1;q=cell_next[q]){
                    double dx=px[p]-px[q],dy=py[p]-py[q],r2=dx*dx+dy*dy;
                    if(r2>=rc2||r2==0.0)continue; if(is_bonded(p,q))continue;
                    sr+=rep_pot(sqrt(r2));
                }
            }
        }
    }
    *ek=sk; *es=ss; *er=sr;
}

/* ============================================================================
 *  HEALING / FRENKEL-PAAR-REKOMBINATION  (nur wenn HEALING=1)
 * ---------------------------------------------------------------------------
 *  Eine gerissene Bindung kann sich WIEDER schliessen, wenn:
 *    (a) Die Atome nah genug sind: RCUT <= r < HEALING_DIST * L0
 *        Unterhalb RCUT wirkt die Abstossung -> Energiebilanz unvollstaendig,
 *        daher dort kein Healing.  Physikalisch: kollidierende Atome heilen nicht.
 *    (b) Relativgeschwindigkeit klein: |v_a - v_b| < HEALING_VREL
 *        Schnelle Atome wuerden sofort wieder auseinanderfliegen -> kein Healing.
 *
 *  Energiebuchhaltung: Beim urspruenglichen Bruch wurde 0.5*K*(r-L0)^2 in
 *  E_broken eingebucht.  Beim Heilen wird der neue Federwert wieder abgezogen.
 *  So bleibt E_total = E_kin + E_spring + E_rep + E_broken stetig.
 * ========================================================================== */
static void heal_bonds(void){
    if(!HEALING) return;
    double hd = HEALING_DIST*L0;
    for(int b=0;b<n_bonds;b++){
        if(bond_intact[b]) continue;    /* intakt -> kein Healing noetig */
        int a=bond_a[b], c=bond_b[b];
        double dx=px[c]-px[a], dy=py[c]-py[a];
        double r=sqrt(dx*dx+dy*dy);
        /* Heilfenster: [RCUT, hd) */
        if(r<RCUT || r>=hd) continue;
        /* Relativgeschwindigkeit pruefen */
        double dvx=vx[c]-vx[a], dvy=vy[c]-vy[a];
        if(sqrt(dvx*dvx+dvy*dvy)>=HEALING_VREL) continue;
        /* Alle Bedingungen erfuellt: Bindung reaktivieren */
        bond_intact[b]=1;
        E_broken -= 0.5*K_SPRING*(r-L0)*(r-L0);  /* Energie zurueckgeben */
        add_neighbor(a,c);  /* Nachbarliste aktualisieren (Abstoss ausschliessen) */
    }
}

/* ============================================================================
 *  VELOCITY-VERLET ZEITSCHRITT  (kick-drift-kick Formulation)
 * ---------------------------------------------------------------------------
 *  Ablauf:
 *    1. Halber Kick:  v += (F/m) * dt/2     (Geschwindigkeit halb aktualisieren)
 *    2. Drift:        x += v * dt           (Position mit neuer Geschw. bewegen)
 *    3. Cell List:    neu aufbauen           (Positionen haben sich geaendert)
 *    4. Healing:      pruefen               (vor Kraftberechnung)
 *    5. Kraefte:      F = F(x_neu)          (mit neuen Positionen)
 *    6. Halber Kick:  v += (F/m) * dt/2     (Geschwindigkeit mit neuen Kraeften)
 *
 *  Warum "symplektisch"?  Der Integrator erhalt eine diskrete Version des
 *  Phasenraumvolumens (Liouville-Theorem).  Das bedeutet keine akkumulierende
 *  Energiedrift, sondern nur eine begrenzte Schwankung um E0.
 * ========================================================================== */
static void step(double dt){
    /* (1+2) Halber Kick + Drift */
    for(int p=0;p<N;p++){
        vx[p]+=0.5*dt*fx[p]/mass[p];   /* Kick */
        vy[p]+=0.5*dt*fy[p]/mass[p];
        px[p]+=dt*vx[p];               /* Drift */
        py[p]+=dt*vy[p];
    }
    /* (3) Cell List neu aufbauen (Positionen haben sich veraendert) */
    cell_build();
    /* (4) Healing pruefen (muss vor Kraftberechnung passieren) */
    heal_bonds();
    /* (5) Neue Kraefte aus aktualisierten Positionen */
    compute_forces();
    /* (6) Zweiter halber Kick mit den neuen Kraeften */
    for(int p=0;p<N;p++){
        vx[p]+=0.5*dt*fx[p]/mass[p];
        vy[p]+=0.5*dt*fy[p]/mass[p];
    }
}

/* ============================================================================
 *  AUSGABE-FUNKTIONEN
 * ========================================================================== */

/* Endzustand als CSV: Position, Verschiebung, gerissene Bindungen.
 * "broken" = Anzahl der verlorenen Nachbarn vs. Anfangsgitter.
 * Randatome haben von vornherein weniger Nachbarn!  Im Python-Auswerteskript
 * werden nur Innenatome als echte Defekte gezaehlt. */
static void dump_state(const char*fn){
    FILE*f=fopen(fn,"w");
    fprintf(f,"x,y,disp,broken\n");
    for(int p=0;p<N;p++){
        double dx=px[p]-px0[p], dy=py[p]-py0[p];
        double disp=sqrt(dx*dx+dy*dy);     /* Verschiebung vom Gitterplatz */
        int broken = (MAXNBR - n_nbr[p]);  /* fehlende Bindungen */
        fprintf(f,"%.4f,%.4f,%.4f,%d\n",px[p],py[p],disp,broken);
    }
    fclose(f);
}

/* Animations-Frame im Extended-XYZ-Format (OVITO/ParaView).
 * Spalten: Element, x, y, z=0, Verschiebung, Geschwindigkeit, gerissene Bindungen.
 * "vmag" zeigt die Schockwellenfront: hohe Geschwindigkeit = Welle gerade hier.
 * "disp" zeigt bleibende Verformung: gross = Atom weit vom Gitterplatz. */
static void dump_frame(FILE*f, double t){
    fprintf(f,"%d\n",N);
    fprintf(f,"Properties=species:S:1:pos:R:3:disp:R:1:vmag:R:1:broken:I:1 Time=%.4f\n",t);
    for(int p=0;p<N;p++){
        double dx=px[p]-px0[p], dy=py[p]-py0[p];
        double disp=sqrt(dx*dx+dy*dy);
        double vmag=sqrt(vx[p]*vx[p]+vy[p]*vy[p]);
        int broken = nbr0[p]-n_nbr[p];  /* nbr0[p] = Nachbarzahl zu Beginn */
        fprintf(f,"A %.4f %.4f 0.0 %.4f %.4f %d\n",px[p],py[p],disp,vmag,broken);
    }
}

/* Positionen der gerissenen Bindungen (ANFANGSPOSITIONEN der Atom-Paare).
 * Zeigt die Schaden-Topologie: wo im urspruenglichen Gitter sind die Risse? */
static void dump_broken_bonds(const char*fn){
    FILE*f=fopen(fn,"w");
    fprintf(f,"xa,ya,xb,yb\n");
    for(int b=0;b<n_bonds;b++){
        if(bond_intact[b]) continue;
        int a=bond_a[b], c=bond_b[b];
        fprintf(f,"%.4f,%.4f,%.4f,%.4f\n",px0[a],py0[a],px0[c],py0[c]);
    }
    fclose(f);
}

/* ============================================================================
 *  CACHE-OPTIMIERUNG: MORTON-SORTIERUNG  (optionales Reorder)
 * ---------------------------------------------------------------------------
 *  PROBLEM: Nach dem Gitteraufbau sind Atom-Indizes row-major (Zeile 0 zuerst).
 *  Die Cell-List-Suche greift auf raeumlich benachbarte Atome zu.  Mit row-major-
 *  Ordnung sind "nah benachbarte" Atome oft WEIT voneinander im Speicher
 *  (Atom oben-rechts: hoher Index; Atom unten-links: niedriger Index) -> hohe
 *  Cache-Miss-Rate, da bei jedem Zugriff eine neue Cache-Zeile geladen wird.
 *
 *  LOESUNG: Z-Order (Morton-Kurve).
 *  Die Morton-Kurve verschraenkt die Bits der x- und y-Zellenkoordinaten:
 *    x = (x3 x2 x1 x0)_2,  y = (y3 y2 y1 y0)_2
 *    -> morton = (y3 x3 y2 x2 y1 x1 y0 x0)_2
 *  Raeumlich benachbarte Zellen haben aehnliche Morton-Zahlen und liegen daher
 *  auch im Speicher nah beieinander.
 *
 *  WICHTIG: Die Physik aendert sich NICHT.
 *  Alle Index-Verweise (Bindungslisten, Nachbarlisten) werden mit dem
 *  Permutationsvektor o2n umgeschrieben.  Der Energieverlauf muss
 *  exakt gleich bleiben -- das ist der Korrektheitstest.
 * ========================================================================== */
static unsigned morton2(unsigned x, unsigned y){
    unsigned z=0;
    for(int i=0;i<16;i++){
        z|=((x>>i)&1u)<<(2*i);    /* x-Bit i an Position 2*i */
        z|=((y>>i)&1u)<<(2*i+1);  /* y-Bit i an Position 2*i+1 */
    }
    return z;
}
static unsigned *g_key;
static int cmp_key(const void*a, const void*b){
    unsigned ka=g_key[*(const int*)a], kb=g_key[*(const int*)b];
    return (ka>kb)-(ka<kb);
}
static void reorder_morton(void){
    unsigned *key=malloc(N*sizeof*key);
    int *perm=malloc(N*sizeof*perm);  /* perm[neuer_idx] = alter_idx */
    int *o2n =malloc(N*sizeof*o2n);   /* o2n[alter_idx]  = neuer_idx */
    /* Morton-Schluessel fuer jeden Atom berechnen */
    for(int p=0;p<N;p++){
        int cx=(int)((px[p]-dom_x0)/cell_size); if(cx<0)cx=0;
        int cy=(int)((py[p]-dom_y0)/cell_size); if(cy<0)cy=0;
        key[p]=morton2((unsigned)cx,(unsigned)cy); perm[p]=p;
    }
    g_key=key; qsort(perm,N,sizeof(int),cmp_key);  /* nach Morton-Key sortieren */
    for(int np=0;np<N;np++) o2n[perm[np]]=np;      /* Umkehrpermutation */
    /* Alle per-Atom-Double-Arrays gemaess perm umordnen */
    #define PD(arr){ double*t=malloc(N*sizeof*t); for(int i=0;i<N;i++)t[i]=arr[perm[i]]; memcpy(arr,t,N*sizeof*t); free(t);}
    PD(px);PD(py);PD(vx);PD(vy);PD(fx);PD(fy);PD(mass);PD(px0);PD(py0);
    #undef PD
    #define PI(arr){ int*t=malloc(N*sizeof*t); for(int i=0;i<N;i++)t[i]=arr[perm[i]]; memcpy(arr,t,N*sizeof*t); free(t);}
    PI(n_nbr); PI(nbr0);
    #undef PI
    /* Nachbarliste (2D-Array [N*MAXNBR]) umordnen */
    { int*t=malloc((size_t)N*MAXNBR*sizeof*t);
      for(int i=0;i<N;i++) memcpy(&t[i*MAXNBR],&nbr[perm[i]*MAXNBR],MAXNBR*sizeof(int));
      memcpy(nbr,t,(size_t)N*MAXNBR*sizeof*t); free(t); }
    /* Alle gespeicherten Atom-Indizes auf neue Nummerierung umschreiben */
    for(int p=0;p<N;p++){ int base=p*MAXNBR; for(int k=0;k<n_nbr[p];k++) nbr[base+k]=o2n[nbr[base+k]]; }
    for(int b=0;b<n_bonds;b++){ bond_a[b]=o2n[bond_a[b]]; bond_b[b]=o2n[bond_b[b]]; }
    free(key); free(perm); free(o2n);
}

/* ============================================================================
 *  MAIN  --  Hauptprogramm
 * ---------------------------------------------------------------------------
 *  Ablauf:
 *    1. Konfiguration lesen
 *    2. Gitter + Bindungen aufbauen
 *    3. Cell List initialisieren
 *    4. Anfangskraefte berechnen
 *    5. PKA(s) platzieren (Anfangsgeschwindigkeit setzen)
 *    6. Optional: Morton-Sortierung
 *    7. Startenergie messen
 *    8. Hauptschleife: NSTEPS Zeitschritte + Energie-Log
 *    9. Ausgabe: Endzustand, Bindungsbrueche
 * ========================================================================== */
int main(int argc, char**argv){
    const char*cfg = (argc>1)? argv[1] : "params.ini";
    read_config(cfg);
    if(strcmp(MODELSEL,"A")!=0)
        fprintf(stderr,"HINWEIS: model=%s, dieses Binary rechnet Modell A "
                       "(Morse = morse_md).\n",MODELSEL);
    srand(SEED);

    build_lattice();
    /* Anfangs-Nachbarzahl speichern (Randatome haben weniger als 6) */
    nbr0=malloc(N*sizeof*nbr0); memcpy(nbr0,n_nbr,N*sizeof*nbr0);
    cell_setup();
    cell_build();
    compute_forces();

    /* PKA(s) platzieren: manuell aus [pkaN] oder auto ueber N_PKA */
    int npka_set = place_pkas();

    /* Optionale Cache-Optimierung (aendert Speicherlayout, nicht die Physik) */
    if(REORDER) reorder_morton();

    /* Startenergie messen (Referenz fuer Drift-Check) */
    double ek,es,er,e0=0;
    energies(&ek,&es,&er);
    e0 = ek+es+er+E_broken;

    /* Simulationsparameter auf stdout ausgeben */
    printf("# Modell A | Config: %s\n",cfg);
    printf("# N=%d bonds=%d  NX=%d NY=%d  K_SPRING=%.1f MAX_STRETCH=%.3f\n",
           N,n_bonds,NX,NY,K_SPRING,MAX_STRETCH);
    printf("# RCUT=%.3f K_REP=%.1f REP_N=%.1f  healing=%d (dist=%.2f vrel=%.2f)\n",
           RCUT,K_REP,REP_N,HEALING,HEALING_DIST,HEALING_VREL);
    printf("# PKA: n=%d (%s) | Default ort=(%.2f,%.2f) E=%.1f angle=%.1f mass=%.2f | dt=%g steps=%d seed=%u\n",
           npka_set, N_MANUAL>0?"manuell":"auto",
           PKA_X,PKA_Y,PKA_ENERGY,PKA_ANGLE,PKA_MASS,DT,NSTEPS,SEED);
    printf("# step      t       E_kin     E_spring    E_rep     E_broken    E_total     drift%%\n");

    /* Ausgabedateien oeffnen */
    char fn_e[256],fn_s[256],fn_b[256];
    snprintf(fn_e,sizeof fn_e,"%s_energy.csv",OUT_PREFIX);
    snprintf(fn_s,sizeof fn_s,"%s_state_final.csv",OUT_PREFIX);
    snprintf(fn_b,sizeof fn_b,"%s_broken_bonds.csv",OUT_PREFIX);
    FILE*elog=fopen(fn_e,"w");
    fprintf(elog,"step,t,E_kin,E_spring,E_rep,E_broken,E_total\n");

    /* XYZ-Animationsdatei (nur wenn dump_every > 0) */
    char fn_xyz[256]; FILE*fxyz=NULL;
    if(DUMP_EVERY>0){ snprintf(fn_xyz,sizeof fn_xyz,"%s.xyz",OUT_PREFIX);
                      fxyz=fopen(fn_xyz,"w"); }

    /* Haupt-Schleife: 0 bis NSTEPS (inklusive) -- erst ausgeben, dann Schritt */
    int logevery = (LOG_EVERY>0)? LOG_EVERY : NSTEPS/40; if(logevery<1)logevery=1;
    for(int s=0;s<=NSTEPS;s++){
        /* Energie ausgeben und Drift pruefen */
        if(s%logevery==0){
            energies(&ek,&es,&er);
            double et=ek+es+er+E_broken;
            double drift = (e0!=0)? 100.0*(et-e0)/fabs(e0) : 0.0;
            printf("%6d %8.3f %10.3f %10.3f %10.3f %10.3f %11.3f %8.4f\n",
                   s, s*DT, ek,es,er,E_broken,et,drift);
            fprintf(elog,"%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                   s,s*DT,ek,es,er,E_broken,et);
        }
        /* Animations-Frame speichern */
        if(fxyz && s%DUMP_EVERY==0) dump_frame(fxyz, s*DT);
        /* Zeitschritt (nicht nach dem letzten Schritt) */
        if(s<NSTEPS) step(DT);
    }
    fclose(elog);
    if(fxyz) fclose(fxyz);

    /* Ergebnis: Anzahl gerissener Bindungen */
    int nbroken=0; for(int b=0;b<n_bonds;b++) if(!bond_intact[b]) nbroken++;
    printf("# RESULT broken_bonds=%d total_bonds=%d healing=%d seed=%u\n",
           nbroken,n_bonds,HEALING,SEED);

    /* Endzustand als CSV fuer Python-Auswertung speichern */
    dump_state(fn_s);
    dump_broken_bonds(fn_b);
    return 0;
}
