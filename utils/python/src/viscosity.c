#include "evolve.h"
void viscosity(void) {
    int i,j,k;
    i=j=k=0;
    double visc, viscm,div_v;
    double res,fac,facp;
    double res2,fac2;
    double res3,fac3;
    double resd,facd;
    double res4;
    visc = 0;
    viscm = 0;
    compute_energy();
#ifdef _OPENMP
    #pragma omp parallel 
    {
    #pragma omp for collapse(3) private(i,j,k,visc,viscm,div_v)
#endif
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {

                visc = params.alpha*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j));
                viscm = params.alpha*.5*(energy[l]*energy[l]+energy[lym]*energy[lym])*sqrt(ymin(j)*ymin(j)*ymin(j));
                div_v = 0.0;
                div_v += (vx[lxp]-vx[l])*SurfX(j,k);
                div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
                div_v *= 2.0/3.0*InvVol(j,k);

                tauxx[l] = visc*dens[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
                tauxx[l] += visc*dens[l]*(vy[lyp]+vy[l])/ymed(j);
                tauyy[l] = visc*dens[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
                tauxy[l] = viscm*.25*(dens[l]+dens[lxm]+dens[lym]+dens[lxm-pitch])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z


            }
        }
    }



    i = j = k = 0;
#ifdef _OPENMP
#pragma omp for private(j,viscm)
#endif
    for(j=1;j<size_y-1;j++) {
        viscm = params.alpha*.5*(energy[l]*energy[l]+energy[lym]*energy[lym])*sqrt(ymin(j)*ymin(j)*ymin(j));
        tauxyavg[j] = viscm*.5*(dbar[j]+dbar[j-1])*( (vxbar[j]-vxbar[j-1])/(ymed(j)-ymed(j-1))-.5*(vxbar[j]+vxbar[j-1])/ymin(j)); 
    }
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for collapse(2) private(i,j,k,res,fac,facp,resd,res2,res3,fac2,fac3,facd,res4)
#endif
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-2;j++) {
            res = 0;
            resd = 0;
            res3 = 0;
            res2 = 0;
            res4 = 0;
            for(i=0;i<size_x;i++) {
                // X
                resd += dens[l];
                
                fac3 = 2.0*(tauxx[l]-tauxx[lxm])/(zone_size_x(j,k)*(dens[l]+dens[lxm]));
                res3 += fac3;
                
                vx_temp[l] += fac3*dt;
                
                fac =  (ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j))*ymed(j));

                
                //facp =  (ymin(j+1)*ymin(j+1)*tauxy[lxp+pitch]-ymin(j)*ymin(j)*tauxy[lxp])/((ymin(j+1)-ymin(j))*ymed(j));

                fac2 = fac*2.0/(ymed(j)*(dens[l]+dens[lxm]));

                //res2 += .5*fac2 + .5*facp*2.0/(ymed(j)*(dens[lxp]+dens[l]));

                res2 += fac2;

                vx_temp[l] += fac2*dt; 

                res += fac;

                //res4 += fac;
                // Y
                vy_temp[l] += 2.0*(ymed(j)*tauyy[l]-ymed(j-1)*tauyy[lym])/((ymed(j)-ymed(j-1))*(dens[l]+dens[lym])*ymin(j))*dt;
                vy_temp[l] += 2.0*(tauxy[lxp]-tauxy[l])/(dx*ymin(j)*(dens[l]+dens[lym]))*dt;
                vy_temp[l] -= (tauxx[l]+tauxx[lym])/(ymin(j)*(dens[l]+dens[lym]))*dt;
                //res += .5*(fac+facp);
            }
            res2 /= (double)nx;
            res3 /= (double)nx;
            resd /= (double)nx;
            res /= (double)nx;
            res4 /= (double)nx;

            LamdepS[j + size_y*3] += dt*res3*ymed(j)*resd;
            LamdepS[j + size_y*2] += dt*res2*ymed(j)*resd;
             drFt[j] += -dt*res;
             dtLd_rhs[j] += res*ymed(j);
            //drFd[j] = -(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]*SurfY(j+1,k) - ymin(j)*ymin(j)*tauxyavg[j]*SurfY(j,k))*InvVol(j,k);
            
            facd  =  -dt*(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]-ymin(j)*ymin(j)*tauxyavg[j])/((ymin(j+1)-ymin(j))*ymed(j));

            LamdepS[j + size_y*5] += facd; //-dt*res;
            drFd[j]  +=  facd; //-dt*res;             
            drFnu[j] += facd;
            drFdB[j]  +=  facd; //-dt*res;             

        }
    }
#ifdef _OPENMP
    }
#endif
    return;
}
