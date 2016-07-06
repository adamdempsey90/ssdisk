#include "evolve.h"
void set_avg(int p) {
    int i,j,k;
    i=j=k=0;
    double resx,resy,resd;
    double resL;
/*
    for(j=0;i<size_y;j++) {
        conv_prefac[j] = ymed(j);
    }

    convolution_2d(dens,vx,Lt,conv_prefac,p*size_y,2*size_y,size_y);
    */

    for(j=0;j<size_y;j++) {
        resx = 0;
        resy = 0;
        resd = 0;
        resL = 0;

        i = k = 0;
        //convolution(&dens[l],&vx[l],Lt,ymed(j),j+p*size_y,size_y*2);
        for(i=0;i<size_x;i++) {
            resx += vx[l];
            resy += vy[l];
            resd += dens[l];
            resL += ymed(j)*dens[l]*(.5*(vx[l] + vx[lxp]) + omf*ymed(j));
        }
        vxbar[j] = resx/(double)nx;
        vybar[j] = resy/(double)nx;
        dbar[j] = resd/(double)nx;

        dbarS[j + p*size_y] = dbar[j];
        Lt[j + p*size_y] = resL/(double)nx;
        if (p==1) dbart[j] += dbar[j]*dt;

        Ld[j + p*size_y] = ymed(j)*( vxbar[j] + omf*ymed(j))*dbar[j];
        LdS[j + p*size_y] = ymed(j)*( vxbar[j] + omf*ymed(j));
        Lw[j + p*size_y] = Lt[j + p*size_y] - Ld[j + p*size_y];
    }
    return;
}

void set_dtLt(double fac) {
    int j;
#ifdef _FFTW
    for(j=0;j<size_y;j++) {
        conv_prefac[j] = ymed(j)*fac;
    }

    convolution_2d(dens,vx,dtLt,conv_prefac,0,size_y,size_y);
#endif
    return;
}
