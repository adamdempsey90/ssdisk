#include "fargo3d.h"


void add_avgs(void) {
    int i,j,k;
    i = j = k = 0;
    
    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;

    real *dens = Density->field_cpu;
    real *vx = Vx->field_cpu;
    real *vy = Vy->field_cpu;

    real *dens_ss = DensitySS->field_cpu;
    real *vx_ss = VxSS->field_cpu;
    real *vy_ss = VySS->field_cpu;


    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                dens[l] = dens_ss[l2D];
                vx[l] = vx_ss[l2D];
                vy[l] = vy_ss[l2D];

            }
        }
    }

#ifdef GPU
    Host2Dev3D(Density);
    Host2Dev3D(Vy);
    Host2Dev3D(Vx);
#endif

    return;
}
