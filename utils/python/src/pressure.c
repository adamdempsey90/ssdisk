#include "evolve.h"
void compute_Pres(void) {
    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                Pres[l] = energy[l]*energy[l]*dens[l];
            }
        }
    }
    return;
}
void compute_energy(void) {
    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                energy[l] = params.h * pow(ymed(j),params.flaringindex) * sqrt(1./ymed(j));
            }
        }
    }
    return;
}
