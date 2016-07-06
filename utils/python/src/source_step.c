#include "evolve.h"
void source_step(void) {
    int i,j,k;
    i=j=k=0;
    double res,resd;
    double vxc;
    compute_Pres();


#ifdef _OPENMP
    #pragma omp parallel for collapse(2) private(i,j,k,res,resd)
#endif
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            res = 0;
            resd = 0;
            for(i=0;i<size_x;i++) {
                resd += dens[l];
                res +=  -2.0/(dens[l]+dens[lxm]) *(Pres[l]-Pres[lxm])/zone_size_x(j,k);
                //res += -dt/(dens[lym] + dens[lxm - pitch]) *(Pres[lym] - Pres[lxm-pitch])/zone_size_x(j-1,k);
            }
            res /=(double)nx;
            resd /= (double)nx;
            dtLd_rhs[j] += dt*res*ymed(j);
            LamdepS[j + size_y*4] += dt*res*ymed(j)*resd;
        }
    }


#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,vxc)
#endif
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {
                // X 
                vx_temp[l] = vx[l] +  -2*dt/(dens[l]+dens[lxm]) *(Pres[l]-Pres[lxm])/zone_size_x(j,k);
                vx_temp[l] -= dt*(Pot[l]-Pot[lxm])/zone_size_x(j,k);
                if (IndirectTerm) {
                    vx_temp[l] -= dt*(indPot[l]-indPot[lxm])/zone_size_x(j,k);
                }
                // Y
                vy_temp[l] = vy[l] -2*dt/(dens[l]+dens[lym])*(Pres[l]-Pres[lym])/(ymed(j)-ymed(j-1));
                vxc = .25*(vx[l]+vx[lxp]+vx[lym]+vx[lxp-pitch])+omf*ymin(j);
                vy_temp[l] += dt*vxc*vxc/ymin(j);
                vy_temp[l] -= dt*(Pot[l]-Pot[lym])/(ymed(j)-ymed(j-1));
                if (IndirectTerm) {
                    vy_temp[l] -= dt*(indPot[l]-indPot[lym])/(ymed(j)-ymed(j-1));
                }
            }
        }
    }


    return;
}
