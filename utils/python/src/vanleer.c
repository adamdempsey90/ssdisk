#include "evolve.h"
void vanleer_y_a(double *q) {
    int i,j,k;
    i=j=k=0;
    double dqm, dqp;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,dqm,dqp)
#endif
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {

                dqm = (q[l] - q[lym])/zone_size_y(j,k);
                dqp = (q[lyp]-q[l])/zone_size_y(j+1,k);

                if (dqp*dqm <= 0) {
                    slope[l] = 0;
                }
                else {
                    slope[l] = 2*dqp*dqm/(dqp+dqm);
                }
            }
        }
    }

    return;
}
void vanleer_y_a_avg(double *q) {
    int j;
    double dqm, dqp;
#ifdef _OPENMP
    #pragma omp parallel for private(j,dqm,dqp)
#endif
    for(j=1;j<size_y-1;j++) {

        dqm = (q[j] - q[j-1])/zone_size_y(j,k);
        dqp = (q[j+1]-q[j])/zone_size_y(j+1,k);

        if (dqp*dqm <= 0) {
            slope[j] = 0;
        }
        else {
            slope[j] = 2*dqp*dqm/(dqp+dqm);
        }
    }

    return;
}
void vanleer_x_a(double *q) {
    int i,j,k;
    i=j=k=0;
    double dqm, dqp;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,dqm,dqp)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {

                dqm = (q[l] - q[lxm]);
                dqp = (q[lxp]-q[l]);

                if (dqp*dqm <= 0) {
                    slope[l] = 0;
                }
                else {
                    slope[l] = 2*dqp*dqm/((dqp+dqm)*zone_size_x(j,k));
                }
            }
        }
    }

    return;
}
void vanleer_y_b(double *q, double *qs, double dt) {
    int i,j,k;
    i=j=k=0;

#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {

                if (vy_temp[l] > 0.) {
                    qs[l] = q[lym] + .5*(zone_size_y(j-1,k)-vy_temp[l]*dt)*slope[lym];
                }
                else {
                    qs[l] = q[l] - .5*(zone_size_y(j,k)+vy_temp[l]*dt)*slope[l];
                }

            }
        }
    }
    return;
}
void vanleer_y_b_avg(double *q, double *qs, double dt) {
    int i,j,k;
    i=j=k=0;
    double res;
#ifdef _OPENMP
    #pragma omp parallel for private(i,j,res)
#endif
    for(j=1;j<size_y-1;j++) {
        res =0 ;
        for(i=0;i<size_x;i++) {
            res += vy_temp[l];
        }
        res /=(double)nx;
        if (res > 0.) {
            qs[j] = q[j-1] + .5*(zone_size_y(j-1,k)-res*dt)*slope[j-1];
        }
        else {
            qs[j] = q[j] - .5*(zone_size_y(j,k)+res*dt)*slope[j];
        }

    }
    return;
}
void vanleer_x_b(double *q, double *qs, double dt,double *vxt) {
    int i,j,k;
    i=j=k=0;

#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {

                if (vxt[l] > 0.) {
                    qs[l] = q[lxm] + .5*(zone_size_x(j,k)-vxt[l]*dt)*slope[lxm];
                }
                else {
                    qs[l] = q[l] - .5*(zone_size_x(j,k)+vxt[l]*dt)*slope[l];
                }

            }
        }
    }
    return;
}
