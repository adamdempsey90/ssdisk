#include "fargo3d.h"
void clear_averages(void) {
    int i,j,k;
    i = j = k =0;
    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;

    real *drfnu = drFluxVisc->field_cpu;
    real *drfd = drFluxDisk->field_cpu;
    real *dtld = dtLDisk->field_cpu;


    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            drfnu[l2D] = 0;
            drfd[l2D] = 0;
            dtld[l2D] = 0;
        }
    }
    return;
}
