#include "evolve.h"
void potential(void) {
    int i,j,k,n;
    i=j=k=0;
    double xpl, ypl,mp;
    double smoothing,rad,distp;

    /*
    double resx = 0;
    double resy = 0;
    double cellmass;
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++)  {
            cellmass = Vol(j,k)*dens[l];
            rad = XC*XC + YC*YC;
            rad = pow(rad+smoothing,-1.5);
            resx += G * cellmass * XC * rad;
            resy += G * cellmass * YC * rad;
        }
    }
    */
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,n,mp,xpl,ypl,distp,smoothing,rad)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                Pot[l] = 0;
                for(n=0;n<nb;n++) {
                    mp = psys[n].mp;
                    xpl = psys[n].x;
                    ypl = psys[n].y;
                    distp = sqrt(xpl*xpl + ypl*ypl);
                    smoothing = params.h*pow(distp,params.flaringindex)*distp*params.soft;
                    smoothing *= smoothing;
                    rad = (XC-xpl)*(XC-xpl) + (YC-ypl)*(YC-ypl);
                    Pot[l] -= G*mp/sqrt(rad + smoothing);
                }
                indPot[l] = 0;
                /*
                indPot[l] = G*planet.mp*(XC*xpl+YC*ypl)/(distp*distp*distp); 
                indPot[l]  -= resx*XC + resy*YC ;
                */
            }
        }
    }
    return;
}
