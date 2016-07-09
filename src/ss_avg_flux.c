#include "fargo3d.h"

void avg_flux(real dt, Field *Qs) {
    
    int i,j,k;
    i = j = k =0;
    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;

    real res1;


#ifdef GPU
    Dev2Host3D(Qs);
    Dev2Host3D(DensStar);
    Dev2Host3D(Vy_temp);
#endif



    real* denstar= DensStar->field_cpu;
    real* vy = Vy_temp->field_cpu;
    real* qs = Qs->field_cpu;

    real* pibarstar = Pibarstar->field_cpu;
    real* drfd = drFluxDisk->field_cpu;
    

// Calculate avg momentum
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res1 = 0;
            for(i=0;i<size_x;i++) {
                res1 += qs[l];
            }
            res1 /= (real)Nx;
            pibarstar[l2D] = res1;
        }
    }

// Now we calculate the flux for deposited torque

    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y-1;j++) {
            
            res1 = 0;

            for(i=0;i<size_x;i++) {
           
                res1 += (vy[l]*denstar[l]*pibarstar[l2D]*SurfY(j,k) - vy[lyp]*denstar[lyp]*pibarstar[l2D+1]*SurfY(j+1,k))*InvVol(j,k);
        
            }
            res1 /= (real)Nx;

            drfd[l2D] -= res1*dt*.5;

        }
    }

    return;

}
