#include "evolve.h"
void set_Lamdep(void) {
    int j;
#ifdef _OPENMP
    #pragma omp parallel for private(j)
#endif
    for(j=NGHY;j<size_y - NGHY;j++) {
        dtLt[j] += (Lt[j+size_y]-Lt[j]);
        dtLd[j] += (Ld[j+size_y]-Ld[j]);
        dtLdS[j] += (LdS[j+size_y]-LdS[j]);
        dtdbar[j] += (dbarS[j+size_y]-dbarS[j]);
        dtLw[j] = dtLt[j] - dtLd[j];
        Lamdep[j] = dtLd[j] + drFd[j];
        LamdepB[j] = dtLd[j] + drFdB[j];
        drFw[j] = drFt[j] - drFd[j];
        drFwB[j] = drFt[j] - drFdB[j];

    }
    int mi;
    for(mi=1;mi<MMAX+2;mi++) {
        for(j=NGHY;j<size_y-NGHY;j++) {
            dtLt[j + mi*size_y] = Lt[j+size_y + mi*size_y*2] - Lt[j+mi*size_y*2]; 
        }
    }

    return;
}
void set_Lamex(void) {
    int i,j,k;
    i=j=k=0;
    
    double res; //,resi;

#ifdef _FFTW
    for(j=NGHY;j<size_y-NGHY;j++) {
        conv_prefac[j] = dt;
    }
    convolution_deriv_2d(&dens[NGHY*nx],&Pot[NGHY*nx],&Lamex[NGHY],&conv_prefac[NGHY],0,size_y,size_y-2*NGHY);
#endif

#ifdef _OPENMP
    #pragma omp parallel for private(i,j,k,res)
#endif
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        //resi = 0;
    //    i = k = 0;
   //     convolution_deriv(&dens[l],&Pot[l],Lamex,-dt,j,size_y);

        for(i=0;i<size_x;i++) {
            res -= dens[l]*(Pot[lxp]-Pot[lxm])/(2*dx);
          //  resi -= dens[l]*(indPot[lxp]-indPot[lxm])/(2*dx);
        }
        res /=(double)nx;
        //resi /=(double)nx;
        Lamex[j] += dt*res;
       // Lamex[j + size_y] += dt*resi;

    }

    return;

}
