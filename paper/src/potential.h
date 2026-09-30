/* ============================================================================
 *  potential.h  --  Kurzreichweitige Abstossung (austauschbar)
 * ----------------------------------------------------------------------------
 *  Gemeinsame Definition der Paar-Abstossung, damit Serial/MPI/CUDA nicht je
 *  eine eigene Kopie pflegen.  Aktuell nur von cascade_serial.c benutzt.
 *
 *  Modelle (Schalter rep_model in params.ini, Sektion [potential]):
 *    r12  -- bisheriges Modell (Standard):
 *              F(r) = K_REP * ((RCUT/r)^n - 1)            fuer r < RCUT
 *            kraftstetig am Cutoff, V(RCUT) = 0.  Wirkt NUR auf ungebundene
 *            Paare (gebundene Paare spueren nur die Feder).
 *    zbl  -- abgeschirmte Coulomb-Abstossung (Ziegler-Biersack-Littmark)
 *            fuer ALLE Paare (auch gebundene), glatt abgeschaltet zwischen
 *            rs1 und rs2 < RCUT.  Siehe docs/plan_zbl.md.
 *
 *  POT_FN erlaubt spaeter __host__ __device__ fuer die CUDA-Version.
 * ========================================================================== */
#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <math.h>

#ifndef POT_FN
#define POT_FN static inline
#endif

enum { REP_R12 = 0, REP_ZBL = 1 };

/* r12: Kraftbetrag (positiv = abstossend) fuer r < rcut. */
POT_FN double r12_force(double r, double rcut, double krep, double n){
    return krep * ( pow(rcut/r, n) - 1.0 );
}

/* r12: Potenzial mit V(rcut) = 0 (Stammfunktion von -F, normiert).
 * Nur fuer die Energiebuchhaltung, nicht fuer die Kraft. */
POT_FN double r12_pot(double r, double rcut, double krep, double n){
    return krep * ( (rcut - pow(rcut,n)*pow(r,1.0-n))/(1.0-n) - (rcut - r) );
}

/* ============================================================================
 *  ZBL  (Ziegler, Biersack, Littmark 1985: universelles abgeschirmtes Coulomb)
 * ---------------------------------------------------------------------------
 *    E(r) = A/r * phi(r/a)
 *    phi(x) = 0.18175 e^-3.1998x + 0.50986 e^-0.94229x
 *           + 0.28022 e^-0.4029x + 0.028171 e^-0.20162x
 *    A = Z1*Z2*e^2,   a = 0.46850 A / (Z1^0.23 + Z2^0.23)
 *
 *  Einheiten: Das Modell ist dimensionslos (L0 = 1, Energie in Einheiten eps).
 *  zbl_setup() rechnet aus Z, L0 [Angstrom] und eps [eV] um:
 *    A_mod = Z^2 * 14.399645 eV*A / (eps * L0_A)     [Energie * L0]
 *    a_mod = 0.46850 A / (2 Z^0.23) / L0_A           [L0]
 *
 *  Abschaltung wie LAMMPS "pair zbl": Fuer rs1 < r < rs2 wird zur Kraft
 *  ein Polynom addiert, sodass E, E' und E'' bei rs2 gleich 0 sind:
 *    E_sw(r) = E(r) + SA/3 t^3 + SB/4 t^4 + SC,   t = r - rs1
 *  Fuer r < rs1 nur die Konstante SC (Kraft unveraendert).
 * ========================================================================== */
typedef struct {
    double A, a;            /* Vorfaktor [Energie*L0], Abschirmlaenge [L0]      */
    double rs1, rs2;        /* Beginn und Ende der Abschaltung [L0]             */
    double SA, SB, SC;      /* Koeffizienten des Abschaltpolynoms               */
} ZblParams;

static const double ZBL_C[4] = {0.18175, 0.50986, 0.28022, 0.028171};
static const double ZBL_D[4] = {3.19980, 0.94229, 0.40290, 0.201620};

/* Rohes ZBL-Potenzial und seine ersten beiden Ableitungen nach r. */
POT_FN void zbl_raw(const ZblParams*z, double r, double*e, double*de, double*d2e){
    double ph=0, dph=0, d2ph=0;
    for(int i=0;i<4;i++){
        double k=ZBL_D[i]/z->a, t=ZBL_C[i]*exp(-k*r);
        ph+=t; dph-=k*t; d2ph+=k*k*t;
    }
    double ir=1.0/r;
    *e   = z->A*ph*ir;
    *de  = z->A*(dph*ir - ph*ir*ir);
    *d2e = z->A*(d2ph*ir - 2.0*dph*ir*ir + 2.0*ph*ir*ir*ir);
}

/* Parameter aus physikalischen Groessen, Abschaltpolynom bestimmen. */
POT_FN void zbl_setup(ZblParams*z, double Zat, double L0_A, double eps_eV,
                      double rs1, double rs2){
    z->A   = Zat*Zat*14.399645/(eps_eV*L0_A);
    z->a   = 0.46850/(2.0*pow(Zat,0.23))/L0_A;
    z->rs1 = rs1; z->rs2 = rs2;
    double e,de,d2e, tc=rs2-rs1;
    zbl_raw(z,rs2,&e,&de,&d2e);
    z->SA = (-3.0*de + tc*d2e)/(tc*tc);
    z->SB = ( 2.0*de - tc*d2e)/(tc*tc*tc);
    z->SC = -e + 0.5*tc*de - tc*tc*d2e/12.0;
}

/* Abgeschaltetes Potenzial, 0 fuer r >= rs2. */
POT_FN double zbl_pot(const ZblParams*z, double r){
    if(r>=z->rs2) return 0.0;
    double e,de,d2e; zbl_raw(z,r,&e,&de,&d2e);
    e += z->SC;
    if(r>z->rs1){ double t=r-z->rs1; e += z->SA/3.0*t*t*t + z->SB/4.0*t*t*t*t; }
    return e;
}

/* Abgeschalteter Kraftbetrag -dE/dr (positiv = abstossend), 0 fuer r >= rs2. */
POT_FN double zbl_force(const ZblParams*z, double r){
    if(r>=z->rs2) return 0.0;
    double e,de,d2e; zbl_raw(z,r,&e,&de,&d2e);
    if(r>z->rs1){ double t=r-z->rs1; de += z->SA*t*t + z->SB*t*t*t; }
    return -de;
}

#endif /* POTENTIAL_H */
