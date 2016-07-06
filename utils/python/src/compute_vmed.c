#include "evolve.h"


void compute_vmed(double *vt) {
    int i,j,k;
    i=j=k=0;
    double res;

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) private(i,j,k,res)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res = 0;
            for(i=0;i<size_x;i++) {
                res += vt[l];
            }
            vmed[l2D] = res/(double)nx;
        }
    }

    return;


}
