#include "fargo3d.h"

void avg_flux(real dt) {
    
    int i,j,k;
    i = j = k =0;
    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;

    real res1,res2;

    real dqm, dqp;

    real *rho = Density->field_cpu;
    real *denstar= DensStar->field_cpu;
    real *vx  = Vx_temp->field_cpu;
    real *vy = Vy_temp->field_cpu;

    real *pibar = Pibar->field_cpu;
    real *pibarstar = Pibarstar->field_cpu;
    real *slopebar = Slopebar->field_cpu;
    real *drfd = drFluxDisk->field_cpu;

// Calculate avg momentum
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res1 = 0;
            res2 = 0;
            for(i=0;i<size_x;i++) {
                res1 += rho[l]*ymed(j)*(vx[l] + OMEGAFRAME*ymed(j));
                res2 += rho[l];
            }
            res1 /= (real)Nx;
            res2 /= (real)Nx;
            pibar[l2D] = res1/res2;
        }
    }

// Calculte slope

    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            dqm = (pibar[l2D] - pibar[l2D-1])/zone_size_y(j,k);
            dqp = (pibar[l2D+1] - pibar[l2D])/zone_size_y(j+1,k);
            
            if (dqp*dqm <= 0) {
                slopebar[l2D] = 0;
            }
            else {
                slopebar[l2D] = 2*dqp*dqm/(dqp + dqm);
            }
        }
    }

// Transport momentum to edge 
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            
            res1 = 0;

            for(i=0;i<size_x;i++) {
                res1 += vy[l];
            }
            res1 /= (real)Nx;

            if (res1 > 0.0) {
                pibarstar[l2D] = pibar[l2D-1] + .5*(zone_size_y(j-1,k) - res1*dt)*slopebar[l2D-1];
            }
            else {
                pibarstar[l2D] = pibar[l2D] - .5*(zone_size_y(j,k)+res1*dt)*slopebar[l2D];
            }
            

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

            drfd[l2D] -= res1*dt;

        }
    }

    return;

}
