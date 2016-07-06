#include "evolve.h"



void advect_shift(double *q, int *nshift) {
    /* Assume NGHX = 0 */
    int i,j,k;
    i=j=k=0;
    
    int itarget,ltarget;

#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,itarget,ltarget)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
            

	            itarget = i-nshift[l2D_int]; // l2D_int
	            while (itarget <  0)  {
                    itarget += nx;
                }
	            while (itarget >= nx+0) {
                    itarget -= nx;
                }
	            ltarget = l-i+itarget;
	            Pres[l] = q[ltarget];
            }
        }
    }

    memcpy(q,Pres,sizeof(double)*size_x*size_y*size_z);
    return;

}
void compute_residuals(double dt) {
    int i,j,k;
    i=j=k=0;
    double ntilde, nround;

#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,ntilde,nround)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                ntilde = vmed[l2D]*dt/zone_size_x(j,k);
                nround = floor(ntilde+0.5);
                if(i == 0) {
                    nshift[l2D_int] = (int)nround;
                }
                vx[l] = vx_temp[l] - vmed[l2D];
                vx_temp[l] = (ntilde-nround)*zone_size_x(j,k)/dt;
            }
        }
    }




    return;
}
/*
void fargo_transport(void) {

  compute_residual(dt);
  transportX (vx, dt); // Vx => variable residual
  vanleer_ppa (vx_temp, dt); // Vx_temp => fixed residual @ given r. This one only is done with PPA
  advect_shift(Pixp, nshift);
  advect_shift(Pixm, nshift);
  advect_shift(Piyp, nshift);
  advect_shift(Piym, nshift);
  advect_shift(dens, nshift);


    return;
}
*/
