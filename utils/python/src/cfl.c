#include "evolve.h"
double cfl(void) {
    int i,j,k;
    i=j=k=0;
    double soundspeed2,soundspeed,visc;
    double cfl1_a, cfl1_b, cfl1;
    double cfl7_a, cfl7_b, cfl7;
    double cfl5_a,cfl5_b,cfl5;
    double cfl2, cfl3;
    double res,fac,vxx,vxxp;
    cfl5_a = 0;
    cfl5_b = 0;
    res = 1e99;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,vxx,vxxp,soundspeed2,visc,soundspeed,cfl1_a,cfl1_b,cfl1,cfl2,cfl3,cfl7_a,cfl7_b,cfl7,cfl5_a,cfl5_b,cfl5) reduction(min:res)
#endif
    for(k=0;k<size_z;k++) {
        for(j=NGHY;j<size_y-NGHY;j++) {
            for(i=0;i<size_x;i++) {
#ifdef FARGO
                vxx = vx[l] - vmed[l2D];
                vxxp = vx[lxp] - vmed[l2D];
#else
                vxx = vx[l];
                vxxp = vx[lxp];
#endif
                soundspeed2 = energy[l]*energy[l];
                visc = params.alpha*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j));
                soundspeed = sqrt(soundspeed2);
                cfl1_a = soundspeed/zone_size_x(j,k);
                cfl1_b = soundspeed/zone_size_y(j,k);
                cfl1 = fmax(cfl1_a,cfl1_b);

                cfl2 =fmax(fabs(vxx),fabs(vxxp))/zone_size_x(j,k);
                cfl3 = fmax(fabs(vy[l]),fabs(vy[lyp]))/zone_size_y(j,k);

                cfl7_a = 1.0/zone_size_x(j,k);	
                cfl7_b = 1.0/zone_size_y(j,k);
                cfl7 = 4.0*visc*pow(fmax(cfl7_a,cfl7_b),2);
#ifdef ARTIFICIALVISCOSITY
                cfl5_a = fabs(vx[lxp]-vx[l])/zone_size_x(j,k);
                cfl5_b = fabs(vy[lyp]-vy[l])/zone_size_y(j,k);
                cfl5 = fmax(cfl5_a, cfl5_b)*4.0*CVNR;
#else
                cfl5 = 0;
#endif
	            fac = CFL/sqrt(cfl1*cfl1 + cfl2*cfl2 + 
                        cfl3*cfl3 + cfl5*cfl5 +  cfl7*cfl7);

                if (fac < res) {
                    res = fac;
                }
            }
        }
    }
    
    return res;

}
