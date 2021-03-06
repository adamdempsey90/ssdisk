#include "fargo3d.h"

void start_Ld_avg(void) {
    int i,j,k;
    i = j = k = 0;

    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;
    
    real res,res1;
#ifdef GPU
    Dev2Host3D(Vx);
    Dev2Host3D(Density);
#endif
    real* vx = Vx->field_cpu;
    real* rho = Density->field_cpu;
    real* vxss = VxSS->field_cpu;
    real* dtld = dtLDisk->field_cpu;

    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res = 0;
            res1 = 0;
            for(i=0;i<size_x;i++) {
                res += vx[l];
                res1 += rho[l];
            }
            res /= (real)Nx;
            res1 /= (real)Nx;
            vxss[l2D] = res;
            dtld[l2D] -= res1*ymed(j)*(res + OMEGAFRAME*ymed(j));
        }
    }

    return;

}
void end_Ld_avg(void) {
    int i,j,k;
    i = j = k = 0;
    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;

    real res,res1;
    
#ifdef GPU
    Dev2Host3D(Vx);
    Dev2Host3D(Density);
#endif

    real *vx = Vx->field_cpu;
    real *rho = Density->field_cpu;
    real *dtld = dtLDisk->field_cpu;

    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res = 0;
            res1 = 0;
            for(i=0;i<size_x;i++) {
                res1 += rho[l];
                res += vx[l];
            }
            res /= (real)Nx;
            res1 /= (real)Nx;
            dtld[l2D] += res1*ymed(j)*(res + OMEGAFRAME*ymed(j));
        }
    }

    return;

}

