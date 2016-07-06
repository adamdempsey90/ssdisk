#include "evolve.h"


void artificial_visc(void) {
    int i,j,k;
    i=j=k=0;
    double dvx,dvy;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,dvx,dvy)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {

	            dvx = vx[lxp]-vx[l];
	            if (dvx < 0.0) {
	                Piym[l] = CVNR*CVNR*dens[l]*dvx*dvx;
	            }
	            else {
	                Piym[l] = 0.0;
	            }
                dvy = vy[lyp]-vy[l];
                if (dvy < 0.0) {
                  Piyp[l] = CVNR*CVNR*dens[l]*dvy*dvy;
                }
                else {
                  Piyp[l] = 0.0;
                }
            }
        }
    }


    double res,fac,resd;
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) private(i,j,k,res,fac,resd)
#endif
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            res = 0;
            resd = 0;
            fac = 0;
            for(i=0;i<size_x;i++) {
                fac =  - 2.0*(Piym[l]-Piym[lxm])/(dens[l]+dens[lxm])/zone_size_x(j,k);
                res += fac;
                resd += dens[l];
	            vx_temp[l] += fac*dt; 
                vy_temp[l] += - 2.0*(Piyp[l]-Piyp[lym])/(dens[l]+dens[lym])*dt/zone_size_y(j,k);           
            }
            res /= (double)nx;
            resd /= (double)nx;
            LamdepS[j + size_y*6] += ymed(j)*res*resd*dt;
        }
    }

    return;

}
