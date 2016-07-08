#include "fargo3d.h"



void take_average(real dt) {

    int i,j,k;

    i = j = k =0;

    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;



    real *drfd = drFluxDisk->field_cpu;
    real *drfnu = drFluxVisc->field_cpu;
    real *dtld = dtLDisk->field_cpu;


    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            drfd[l2D] /= dt;
            drfnu[l2D] /= dt;
            dtld[l2D] /= dt;
        }
    }

    return;
}
