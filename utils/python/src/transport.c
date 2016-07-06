#include "evolve.h"
void vel_to_temp(void) {
    memcpy(vx_temp,vx,sizeof(double)*size_x*size_y*size_z);
    memcpy(vy_temp,vy,sizeof(double)*size_x*size_y*size_z);
/*
    int i,j,k;
    i = j = k =0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                vx_temp[l] = vx[l];
                vy_temp[l] = vy[l];
            }
        }
    }
*/
    return;
}
void temp_to_vel(void) {
    memcpy(vx,vx_temp,sizeof(double)*size_x*size_y*size_z);
    memcpy(vy,vy_temp,sizeof(double)*size_x*size_y*size_z);
/*
    int i,j,k;
    i = j = k = 0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                vx[l] = vx_temp[l];
                vy[l] = vy_temp[l];
            }
        }
    }
*/
    return;
}
void transport_step(void) {
    set_momenta();    
    transportY();
#ifdef FARGO
    compute_residuals(dt);
    transportX(vx,FALSE);
    transportX(vx_temp,TRUE);
    advect_shift(Pixp, nshift);
    advect_shift(Pixm, nshift);
    advect_shift(Piyp, nshift);
    advect_shift(Piym, nshift);
    
    advect_shift(dens, nshift);
   
#else
    transportX(vx_temp,FALSE);
#endif
    set_vel();
    return;
}

void set_momenta(void) {
    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z; k++) {
        for(j=0;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {
                Pixm[l] = ymed(j)*(vx_temp[l]+omf*ymed(j))*dens[l];
                Pixp[l] = ymed(j)*(vx_temp[lxp]+omf*ymed(j))*dens[l];
                Piym[l] = vy_temp[l]*dens[l];
                Piyp[l] = vy_temp[lyp]*dens[l];
            }
        }
    }

    return;
}
void transportY(void) {
    // Y direction

    vanleer_y_a(dens);
    vanleer_y_b(dens,denstar,dt);

    DividebyRho(Pixm);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    update_flux(Qs,divrho);
    updateY(Pixm,Qs,dt);

    DividebyRho(Pixp);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    update_flux(Qs,divrho);
    updateY(Pixp,Qs,dt);

    DividebyRho(Piym);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    updateY(Piym,Qs,dt);

    DividebyRho(Piyp);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    updateY(Piyp,Qs,dt);


    vanleer_y_a_avg(dbar);
    vanleer_y_b_avg(dbar,dbarstar,dt);

    DividebyRhoavg(Ld);
    vanleer_y_a_avg(divrho);
    vanleer_y_b_avg(divrho,Qs,dt);
    update_flux_avg(Qs,divrho);

    
    update_density_Y(dt);



    return;

}
void DividebyRho(double *q) {
    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                divrho[l] = q[l]/dens[l];
            }
        }
    }
    return;
}
void DividebyRhoavg(double *q) {
    int i,j,k;
    double resv;
    i=j=k=0;
    /*
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) private(i,j,k,resv)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            resv  = 0;
            for(i=0;i<size_x;i++) {
                resv += ( .5*(vx_temp[l]+vx_temp[lxp]) + ymed(j)*omf)*ymed(j);
            }
            resv /= (double)nx;
            divrho[j] = resv;
        }
    }
    */
    
    for(j=0;j<size_y;j++) {
        divrho[j] = q[j]/dbar[j];
    }

    return;
}
void transportX(double *vxt, int ppa) {
    // X direction

    if (ppa) {
        vanleer_ppa(dt,dens,denstar,vxt);
    }
    else {
        vanleer_x_a(dens);
        vanleer_x_b(dens,denstar,dt,vxt);
    }
/* Pixm */
    DividebyRho(Pixm);
    if (ppa) {
        vanleer_ppa(dt,divrho,Qs,vxt);
    }
    else {
        vanleer_x_a(divrho);
        vanleer_x_b(divrho,Qs,dt,vxt);
    }
    updateX(Pixm,Qs,dt,vxt);

/* Pixp */
    DividebyRho(Pixp);
    if (ppa) {
        vanleer_ppa(dt,divrho,Qs,vxt);
    }
    else {
        vanleer_x_a(divrho);
        vanleer_x_b(divrho,Qs,dt,vxt);
    }
    updateX(Pixp,Qs,dt,vxt);


/* Piym */
    DividebyRho(Piym);
    if (ppa) {
        vanleer_ppa(dt,divrho,Qs,vxt);
    }
    else {
        vanleer_x_a(divrho);
        vanleer_x_b(divrho,Qs,dt,vxt);
    }
    updateX(Piym,Qs,dt,vxt);
/* Piyp */
    DividebyRho(Piyp);
    if (ppa) {
        vanleer_ppa(dt,divrho,Qs,vxt);
    }
    else {
        vanleer_x_a(divrho);
        vanleer_x_b(divrho,Qs,dt,vxt);
    }
    updateX(Piyp,Qs,dt,vxt);

/* Density */
    update_density_X(dt,vxt);

    return;

}
void set_vel(void) {

    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                vy[l] = (Piym[l] + Piyp[lym])/(dens[l]+dens[lym]);
                vx[l] = (Pixm[l] + Pixp[lxm])/(ymed(j)*(dens[l]+dens[lxm])) - omf*ymed(j);
            }
        }
    }
/*
    for(j=1;j<size_y;j++) {
        resv = 0;
        for(i=0;i<size_x;i++) {
            resv += .5*(Pixm[l] + Pixp[l]);
            //resv += dens[l]*(.5*(vx[lxp]+vx[l])+omf*ymed(j))*ymed(j);
        }
        resv /= (double)nx;
        dtLt[j] = -resv;
        dtLd[j] = -dbar[j]*(vx[j] + omf*ymed(j))*ymed(j);
    }
*/
    return;
}
