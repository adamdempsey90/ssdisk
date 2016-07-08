#include "fargo3d.h"

void start_Ld_avg(void) {
    int i,j,k;
    i = j = k = 0;

    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    real res;
    real *vxss = VxSS->field_cpu;

    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res = 0;
            for(i=0;i<size_x;i++) {
                res += vx[l];
            }
            res /= (real)nx;
            vxss[l2D] = res;
            dtLDisk[l2D] -= ymed(j)*(res + omf*ymed(j));
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

    real res;

    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res = 0;
            for(i=0;i<size_x;i++) {
                res += vx[l];
            }
            res /= (real)nx;
            dtLDisk[l2D] += ymed(j)*(res + omf*ymed(j));
        }
    }

    return;

}
