#include "evolve.h"
void updateX(double *q, double *qs,double dt,double *vxt) {
    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {

                q[l] += ((vxt[l]*qs[l]*denstar[l]-vxt[lxp]*qs[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

            }
        }
    }

    return;
}
void update_flux_avg(double *qs, double *q) {
    int i,j,k;
    i=j=k=0;
    double res, res1,res2,fac;
    double vravg, vravgp;
    double davg, davgp;
#ifdef _FFTW
    for(j=0;j<size_y-1;j++) {
        conv_prefac[j] = -dt*qs[j]*SurfY(j,k)*InvVol(j,k);
    }
    convolution_2d(vy_temp,denstar,drFd,conv_prefac,0,size_y,size_y-1);

    for(j=0;j<size_y-1;j++) {
        conv_prefac[j] = dt*qs[j+1]*SurfY(j+1,k)*InvVol(j,k);
    }
    convolution_2d(&vy_temp[size_x],&denstar[size_x],drFd,conv_prefac,0,size_y,size_y-1);
#endif
#ifdef _OPENMP
    //#pragma omp parallel for private(i,j,k,res,fac,res2,res1,vravg,vravgp,davg,davgp)
#endif
    for(j=0;j<size_y-1;j++) {
        i = 0;
     //   convolution(&vy_temp[l],&denstar[l],drFd,-dt*qs[j]*SurfY(j,k)*InvVol(j,k),j,size_y);
    //    convolution(&vy_temp[lyp],&denstar[lyp],drFd,dt*qs[j+1]*SurfY(j+1,k)*InvVol(j,k),j,size_y);
    
        
        res = 0;
        res1 = 0;
        res2 = 0;
        fac = 0;
        vravg = 0;
        vravgp = 0;
        davg = 0;
        davgp = 0;
        for(i=0;i<size_x;i++) {    
//            fac += (denstar[lyp]*vy_temp[lyp]*SurfY(j+1,k)*(qs[j+1] - q[j]) - denstar[l]*vy_temp[l]*SurfY(j,k)*(qs[j] - q[j]))*InvVol(j,k);
            vravg += vy_temp[l];
            vravgp += vy_temp[lyp];
            davg += denstar[l];
            davgp += denstar[lyp];

            fac = ((vy_temp[l]*qs[j]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[j+1]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
            res += fac;
            //res2 += fac - .5*(qs[j] + qs[j+1])*((vy_temp[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
            
            res2 += dens[l]*.5*(vy_temp[l]+vy_temp[lyp])*(qs[j+1]-qs[j])/(ymin(j+1)-ymin(j));

            //res2 += (vy_temp[lyp]*denstar[lyp]*(qs[j+1]-q[j]) - vy_temp[l]*denstar[l]*(qs[j]-q[j]))/(ymin(j+1)-ymin(j));
        //    res2 += fac - q[j]*((vy_temp[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
            //res += ((vy_temp[l]*qs[j]*denstar[l]*ymin(j)-vy_temp[lyp]*qs[j+1]*denstar[lyp]*ymin(j+1)))/(ymin(j+1)-ymin(j));
        }
        res /=(double)nx;
        res2 /= (double)nx;
        fac /=(double)nx;
        vravg /= (double)nx;
        vravgp /= (double)nx;
        davg /= (double)nx;
        davgp /= (double)nx;
        
        res1 = ((vravg*qs[j]*davg*SurfY(j,k)-vravgp*qs[j+1]*davgp*SurfY(j+1,k))*InvVol(j,k));

        drFd[j] -= dt*res;
        mdotl[j] -= dt*res;
        drFdB[j] -= dt*res1;
        LamdepS[j] += dt*res2;
    }
    return;
}
void update_flux(double *qs,double *q) {
    int i,j,k;
    i=j=k=0;
    double res,fac,facd,res2;
#ifdef _OPENMP
    #pragma omp parallel for private(i,j,k,res,fac,facd,res2)
#endif
    for(j=0;j<size_y-1;j++) {
        res = 0;
        res2 = 0;
        fac  = 0;
        facd = 0;
        for(i=0;i<size_x;i++) {    
            facd += dens[l];
            //fac += (SurfY(j+1,k)*vy_temp[lyp]*(qs[lyp]-q[l]) - SurfY(j,k)*vy_temp[l]*(qs[l]-q[l]))*InvVol(j,k);
            fac = ((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
            res += fac;

            //res2 += (( vy_temp[lyp]*(qs[lyp]-q[l]) - vy_temp[l]*(qs[l]-q[l]))/(ymin(j+1)-ymin(j)));
            
            res2 += .5*(vy_temp[l]+vy_temp[lyp])*(qs[lyp]-qs[l])/(ymin(j+1)-ymin(j));


            //res2 += ((vy_temp[l]*qs[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*SurfY(j+1,k))*InvVol(j,k));
            //res2 -= q[l]*((vy_temp[l]*SurfY(j,k) - vy_temp[lyp]*SurfY(j+1,k)))*InvVol(j,k);
        }
        res /=(double)nx;
        res2 /= (double)nx;
        fac /=(double)nx;
        facd /=(double)nx;
        drFt[j] -= dt*res*.5;
        LamdepS[j + size_y] -= .5*res2*dt*facd; 
        dtLd_rhs[j] -= fac*dt*.5;
    }
    return;
}
void updateY(double *q, double *qs,double dt) {
    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {
                
                q[l] += dt*((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
            }
        }
    }
    return;
}
void update_density_Y(double dt) {
    int i,j,k;
    i=j=k=0;
    double res;
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) private(i,j,k,res)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y-1;j++) {
            res = 0;
            for(i=0;i<size_x;i++) {
                res += vy_temp[l]*denstar[l];
                dens[l] += dt*((vy_temp[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));

            }
            res /= (double)nx;
            mdotavg[j] += res * -2*M_PI*ymin(j)*dt;
        }
    }
    return;
}
void update_density_X(double dt,double *vxt) {
    int i,j,k;
    i = j= k =0;

#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {

                dens[l] += ((vxt[l]*denstar[l]-vxt[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

            }
        }
    }
    return;
}
