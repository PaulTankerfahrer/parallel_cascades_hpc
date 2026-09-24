/* zbl_scatter_test.c -- Zwei-Koerper-Streuung am ZBL-Potenzial (potential.h).
 *
 * Prueft Potenzial + adaptive Zeitschrittregel gegen die klassische Theorie:
 *   - Mindestabstand r_min:  1 - V(r)/E_c - b^2/r^2 = 0
 *   - Streuwinkel im Schwerpunktsystem
 *       theta_cm = pi - 2 b INT_{r_min}^inf dr / (r^2 sqrt(1 - V/E_c - b^2/r^2))
 *     (Substitution u = r_min/r, u = 1-w^2 entfernt die Wurzelsingularitaet)
 *   Gleiche Massen m=1: E_c = E_lab/2, Laborwinkel des Projektils = theta_cm/2.
 *
 * Simulation: Velocity-Verlet, dt wie choose_dt() in cascade_serial.c.
 * Bauen:  gcc -O2 -Isrc -o zbl_scatter scripts/ci/zbl_scatter_test.c -lm
 * Exit-Code 1, wenn ein Fall die Toleranz verfehlt.
 */
#include <stdio.h>
#include <math.h>
#include "potential.h"

static ZblParams Z;

static double g_of(double r, double b, double Ec){
    return 1.0 - zbl_pot(&Z,r)/Ec - b*b/(r*r);
}

/* r_min per Bisektion: g < 0 innen, g > 0 aussen */
static double theory_rmin(double b, double Ec){
    double lo=1e-6, hi=Z.rs2+b+1.0;
    for(int i=0;i<200;i++){ double m=0.5*(lo+hi); if(g_of(m,b,Ec)<0) lo=m; else hi=m; }
    return hi;
}

static double theory_theta_cm(double b, double Ec, double rmin){
    int n=200000; double s=0;
    for(int i=0;i<n;i++){
        double w=(i+0.5)/n, r=rmin/(1.0-w*w);
        double g=g_of(r,b,Ec); if(g<=0) g=1e-300;
        s += 2.0*w/sqrt(g);
    }
    s/=n;
    return M_PI - 2.0*b/rmin*s;
}

/* Projektil (0) mit E_lab auf ruhendes Target (1) im Ursprung, Stossparameter b */
static void simulate(double Elab, double b, double dxmax, double*rmin, double*th_lab, double*derr){
    double x[2]={-(Z.rs2+0.2), 0}, y[2]={b,0}, vx[2]={sqrt(2*Elab),0}, vy[2]={0,0};
    double fx[2], fy[2];
    #define FORCES() { double dx=x[0]-x[1], dy=y[0]-y[1], r=sqrt(dx*dx+dy*dy); \
        double f=zbl_force(&Z,r)/r; fx[0]=f*dx; fy[0]=f*dy; fx[1]=-fx[0]; fy[1]=-fy[0]; }
    FORCES();
    double e0=Elab; *rmin=1e30;
    for(long it=0; it<50000000; it++){
        double dx=x[0]-x[1], dy=y[0]-y[1], r=sqrt(dx*dx+dy*dy);
        if(r<*rmin) *rmin=r;
        if(it>10 && r>Z.rs2+0.2 && (dx*(vx[0]-vx[1])+dy*(vy[0]-vy[1]))>0) break;
        double v=0,a=0;
        for(int k=0;k<2;k++){ double vv=hypot(vx[k],vy[k]), aa=hypot(fx[k],fy[k]);
            if(vv>v)v=vv; if(aa>a)a=aa; }
        double dt = 2*dxmax/(v+sqrt(v*v+2*a*dxmax));   /* wie choose_dt() */
        for(int k=0;k<2;k++){ vx[k]+=0.5*dt*fx[k]; vy[k]+=0.5*dt*fy[k]; x[k]+=dt*vx[k]; y[k]+=dt*vy[k]; }
        FORCES();
        for(int k=0;k<2;k++){ vx[k]+=0.5*dt*fx[k]; vy[k]+=0.5*dt*fy[k]; }
    }
    *th_lab = atan2(vy[0],vx[0]);
    double e1=0.5*(vx[0]*vx[0]+vy[0]*vy[0]+vx[1]*vx[1]+vy[1]*vy[1]);
    *derr=fabs(e1-e0)/e0;
}

int main(void){
    zbl_setup(&Z, 74.0, 2.74, 10.0, 0.50, 0.85);
    double dxmax=0.0025;   /* Standard in cascade_serial.c */
    double Es[]={10,100,1000,1e4,1e5}, bs[]={0.0,0.02,0.05,0.1,0.2,0.4};
    int fail=0;
    printf("# ZBL W-W, eps=10 eV, dx_max=%g   (theta in Grad, Labor)\n",dxmax);
    printf("#   E_lab      b   rmin_th  rmin_sim   th_th   th_sim   |dE|/E\n");
    for(int i=0;i<5;i++) for(int j=0;j<6;j++){
        double E=Es[i], b=bs[j], Ec=E/2;
        double rt=theory_rmin(b,Ec), tt=0.5*theory_theta_cm(b,Ec,rt)*180/M_PI;
        double rs,ts,de; simulate(E,b,dxmax,&rs,&ts,&de); ts*=180/M_PI;
        if(b==0.0) tt=0.0;  /* zentraler Stoss: Projektil bleibt liegen, Winkel undefiniert */
        double er=fabs(rs-rt)/rt, et=fabs(ts-tt);
        int bad = er>1e-2 || (b>0 && et>0.05*fabs(tt)+0.05) || de>2e-4;
        fail|=bad;
        printf("%9.0f %6.2f %9.4f %9.4f %7.3f %8.3f %8.1e %s\n",E,b,rt,rs,tt,b>0?ts:0.0,de,bad?"FEHLER":"");
    }
    printf(fail? "== ZBL-Streutest FEHLGESCHLAGEN ==\n" : "== ZBL-Streutest bestanden ==\n");
    return fail;
}
