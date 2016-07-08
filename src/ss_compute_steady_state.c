#include "fargo3d.h"

void compute_steady_state(void) {
    int i,j,k;
    i = j = k = 0;
    
    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;

    real y_min = YMIN;
    real y_max = YMAX;


    real *energy = Energy->field_cpu;
    real *rho = Density->field_cpu;

    real *densityss = DensitySS->field_cpu;
    real *vyss = VySS->field_cpu;
    real *vxss = VxSS->field_cpu;
    real *drfd = drFluxDisk->field_cpu;
    real *drfnu = drFluxVisc->field_cpu;
    real *dtld = dtLDisk->field_cpu;

    real mdot = MDOT;
    real fac,norm,lamdep,ld,visc;
    real intlam = 0;
    real wkzin = 0;
    real wkzout = 1e99;
#ifdef STOCKHOLM
    wkzin = y_min + .0476*(y_max-y_min);
    wkzout = y_max - .19*(y_max-y_min);
#endif

    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
#ifdef ALPHAVISCOSITY
#ifdef ISOTHERMAL
            visc= ALPHA*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j)/(G*MSTAR));
#else
                // Not good for non-isothermal
            visc= ALPHA*GAMMA*(GAMMA-1.0)*(energy[l])/(rho[l])*sqrt(ymed(j)*ymed(j)*ymed(j)/(G*MSTAR));
#endif
#else
            visc =  NU;
#endif
            norm = mdot/(3*M_PI*visc);
            fac = 2*M_PI*ymed(j);
            // Exclude wavekilling zones
            if ((ymed(j) >= wkzin) && (ymed(j) <= wkzout)) {
                lamdep += fac*zone_size_y(j,k)*(drfnu[l2D] + drfd[l2D] + dtld[l2D]);
            }
                ld = (vxss[l2D] + OMEGAFRAME*ymed(j))*ymed(j);                
                densityss[l2D] = norm*(1 + lamdep/(ld * mdot)); // lamdep = 0 in inner wkz and constant in outer wkz
                vyss[l2D] = -mdot/(fac*densityss[l2D]);
            
        }
    }
            

    return;
}   
