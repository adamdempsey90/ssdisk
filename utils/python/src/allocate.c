#include "evolve.h"
void allocate_all(void) {
    
    MALLOC_SAFE((Ymin = (double *)malloc(sizeof(double)*(size_y+1))));
    MALLOC_SAFE((Xmin = (double *)malloc(sizeof(double)*(size_x+1))));

    MALLOC_SAFE((Ymed = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Xmed = (double *)malloc(sizeof(double)*(size_x))));
    
    MALLOC_SAFE((dens0 = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((vmed = (double *)malloc(sizeof(double)*(size_y*size_z))));
    MALLOC_SAFE((vx0 = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((vy0 = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((qR = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((qL = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((dens = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vx = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vx_temp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vy_temp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pres = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((energy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pot = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((indPot = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pixm = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pixp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Piym = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Piyp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((denstar = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Qs = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((slope = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((divrho = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauxx = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauxy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauyy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));

    MALLOC_SAFE((nshift = (int *)malloc(sizeof(int)*(size_y*size_z))));

    MALLOC_SAFE((Lt = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((Ld = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((LdS = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((dbarS = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((Lw = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((drFt = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFd = (double *)malloc(sizeof(double)*(size_y*(MMAX+2)))));
    MALLOC_SAFE((drFdB = (double *)malloc(sizeof(double)*(size_y*(MMAX+2)))));
    MALLOC_SAFE((mdotl = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFnu = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFw = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFwB = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Lamdep = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((LamdepB = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((LamdepS = (double *)malloc(sizeof(double)*(size_y*7))));
    MALLOC_SAFE((Lamex = (double *)malloc(sizeof(double)*(size_y*(MMAX+2)))));
    MALLOC_SAFE((tauxyavg = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dbar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((mdotavg = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dbart = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((vxbar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((vybar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dbarstar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dtLt = (double *)malloc(sizeof(double)*(size_y*(MMAX+2)))));
    MALLOC_SAFE((dtLd = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dtLdS = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dtdbar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dtLd_rhs = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dtLw = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((conv_prefac = (double *)malloc(sizeof(double)*(size_y))));
    int i,j,k;
    i = j = k = 0;

    for(j=0;j<size_y;j++) {
        Ymed[j] = 0;
        Ymin[j] = 0;
        Lt[j] = 0;
        Ld[j] = 0;
        Lw[j] = 0;
        Lt[j+size_y] = 0;
        Ld[j+size_y] = 0;
        Lw[j+size_y] = 0;
        drFt[j] = 0;
        drFd[j] = 0;
        drFdB[j] = 0;
        mdotl[j] = 0;
        drFnu[j] = 0;
        drFw[j] = 0;
        drFwB[j] = 0;
        Lamex[j] = 0;
        Lamex[j+size_y] = 0;
        LdS[j] = 0;
        LdS[j+size_y] = 0;
        dbarS[j] = 0;
        dbarS[j+size_y] = 0;

        Lamdep[j] = 0;
        LamdepB[j] = 0;
        LamdepS[j + size_y*0] = 0;
        LamdepS[j + size_y*1] = 0;
        LamdepS[j + size_y*2] = 0;
        LamdepS[j + size_y*3] = 0;
        LamdepS[j + size_y*4] = 0;
        LamdepS[j + size_y*5] = 0;
        LamdepS[j + size_y*6] = 0;
        tauxyavg[j] = 0;
        vxbar[j] = 0;
        vybar[j] = 0;
        dbar[j] = 0;
        mdotavg[j] = 0;
        dbart[j] = 0;
        dbarstar[j] = 0;
        dtLt[j] = 0;
        dtLd_rhs[j] = 0;
        dtLd[j] = 0;
        dtLdS[j] = 0;
        dtdbar[j] = 0;
        dtLw[j] = 0;
        conv_prefac[j] = 0;
    }
    Ymin[size_y] = 0;
    for(i=0;i<size_x;i++) {
        Xmed[i] = 0;
        Xmin[i] = 0;
    }
    Xmin[size_x] = 0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            vmed[l2D] = 0;
            for(i=0;i<size_x;i++) {
        
                
                dens[l]=0;
                vx[l]=0;
                vy[l]=0;
                vx_temp[l]=0;
                vy_temp[l]=0;
                Pres[l]=0;
                energy[l]=0;
                Pot[l]=0;
                indPot[l]=0;
                Pixm[l]=0;
                Pixp[l]=0;
                Piym[l]=0;
                Piyp[l]=0;
                denstar[l]=0;
                Qs[l]=0;
                slope[l]=0;
                divrho[l]=0;
                tauxx[l]=0;
                tauxy[l]=0;
                tauyy[l]=0;
            }
        }
    }
    int mi;
    for(mi=0;mi<(MMAX+2);mi++) {
        for(j=0;j<size_y;j++) {
            drFd[j + mi*size_y] = 0;
            Lamex[j + mi*size_y] = 0;
            dtLt[j + mi*size_y] = 0;
        }
    }


    return;
}
void free_all(void) {
    
    free(Ymin);
    free(Xmin);

    free(Ymed);
    free(Xmed);
    
    free(dens);
    free(vx);
    free(vy);
    free(vx_temp);
    free(vy_temp);
    free(Pres);
    free(energy);
    free(Pot);
    free(indPot);
    free(Pixm);
    free(Pixp);
    free(Piym);
    free(Piyp);
    free(denstar);
    free(Qs);
    free(slope);
    free(divrho);
    free(tauxx);
    free(tauxy);
    free(tauyy);

    free(Lt);
    free(Ld);
    free(Lw);
    free(drFt);
    free(drFd);
    free(drFw);
    free(Lamdep);
    free(Lamex);
    free(tauxyavg);
    free(dbar);
    free(mdotavg);
    free(dbart);
    free(vxbar);
    free(vybar);
    free(dbarstar);
    free(dtLt);
    free(dtLd);
    free(dtLw);
    return;
}
