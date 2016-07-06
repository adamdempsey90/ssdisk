#include "evolve.h"
double *k1,*k2,*k3,*k4,*k5,*k6,*q0,*q1;

void get_accel(double *q, double *k, double dtn) {
    int i,j;
    double xpl,ypl,zpl,coeff,distp;
    for(i=0;i<nb;i++) {
        k[0+i*6] = q[3+i*6];
        k[1+i*6] = q[4+i*6];
        k[2+i*6] = q[5+i*6];
        k[3+i*6] = 0;
        k[4+i*6] = 0;
        k[5+i*6] = 0;
        for(j=0;j<nb;j++) {
            if (i != j) {
                xpl = q[0 + i*6] - q[0 + j*6];
                ypl = q[1 + i*6] - q[1 + j*6];
                zpl = q[2 + i*6] - q[2 + j*6];
                distp = xpl*xpl + ypl*ypl + zpl*zpl; 
                coeff = -G*psys[j].mp *pow(distp,-1.5);
                k[3 + i*6] += coeff*xpl;
                k[4 + i*6] += coeff*ypl;
                k[5 + i*6] += coeff*zpl;

            }
        }
    }
    for(i=0;i<nb*6;i++) {
        k[i] *= dtn;
    }
    return;
}
void move_planet_step(double dtn) {
    int i,n;
    double cvals1[5] =  {0.2, 0.0, 0.0, 0.0, 0.0};
    double cvals2[5] = {0.075, 0.225, 0.0, 0.0, 0.0};
    double cvals3[5] = { 0.3, -0.9, 1.2, 0.0, 0.0};
    double cvals4[5] = {-11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0, 0.0 };
    double cvals5[5] = {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0 };
    


    for(n=0;n<nb;n++) { 
        q0[0+ n*6] = psys[n].x;
        q0[1+ n*6] = psys[n].y;
        q0[2+ n*6] = psys[n].z;
        q0[3+ n*6] = psys[n].vx;
        q0[4+ n*6] = psys[n].vy;
        q0[5+ n*6] = psys[n].vz;

        for(i=0;i<6;i++) {
            k1[i + n*6] = 0;
            k2[i + n*6] = 0;
            k3[i + n*6] = 0;
            k4[i + n*6] = 0;
            k5[i + n*6] = 0;
            k6[i + n*6] = 0;
            q1[i + n*6] = 0;
        }
    }

    get_accel(q0, k1, dtn);
    for(i=0;i<6*nb;i++) {
        q1[i] = q0[i] + cvals1[0]*k1[i] + cvals1[1]*k2[i] + cvals1[2]*k3[i] + cvals1[3]*k4[i] + cvals1[4]*k5[i];
    }

    get_accel(q1, k2, dtn);
    for(i=0;i<6*nb;i++) {
        q1[i] = q0[i] + cvals2[0]*k1[i] + cvals2[1]*k2[i] + cvals2[2]*k3[i] + cvals2[3]*k4[i] + cvals2[4]*k5[i];
    }


    get_accel(q1, k3, dtn);
    for(i=0;i<6*nb;i++) {
        q1[i] = q0[i] + cvals3[0]*k1[i] + cvals3[1]*k2[i] + cvals3[2]*k3[i] + cvals3[3]*k4[i] + cvals3[4]*k5[i];
    }

    get_accel(q1, k4, dtn);
    for(i=0;i<6*nb;i++) {
        q1[i] = q0[i] + cvals4[0]*k1[i] + cvals4[1]*k2[i] + cvals4[2]*k3[i] + cvals4[3]*k4[i] + cvals4[4]*k5[i];
    }

    get_accel(q1, k5, dtn);
    for(i=0;i<6*nb;i++) {
        q1[i] = q0[i] + cvals5[0]*k1[i] + cvals5[1]*k2[i] + cvals5[2]*k3[i] + cvals5[3]*k4[i] + cvals5[4]*k5[i];
    }

    get_accel(q1, k6, dtn);

    for (i = 0; i < 6*nb; i++) {
        q1[i] = (q0[i]+
                    37.0/378.0*k1[i]  + 
                    250.0/621.0*k3[i] + 
                    125.0/594.0*k4[i] + 
                    512.0/1771.0*k6[i]);
    }
    for(n=0;n<nb;n++) {
        psys[n].x =  q1[0+n*6];
        psys[n].y =  q1[1+n*6];
        psys[n].z =  q1[2+n*6];
        psys[n].vx = q1[3+n*6];
        psys[n].vy = q1[4+n*6];
        psys[n].vz = q1[5+n*6];
    }
    return;
}
void rotate_sys(double angle) {
    int n;
    double xt, yt,vxt,vyt;

    double ca = cos(angle);
    double sa = sin(angle);
    for(n=0;n<nb;n++) {
        xt = psys[n].x;
        yt = psys[n].y;
        vxt = psys[n].vx;
        vyt = psys[n].vy;

        psys[n].x = xt * ca  + yt * sa;
        psys[n].y = -xt * sa  + yt * ca;

        psys[n].vx = vxt * ca  + vyt * sa;
        psys[n].vy = -vxt * sa  + vyt * ca;
    }

    return;
}
void move_planet(void) {
    int i,j,k;
    int subcycles = 5;
    double dt_frac = 1./(double)subcycles;
    double omf_new;
    double xpl0,ypl0;
    double xpl1,ypl1;
    double d1, d2,cross,domega;

    xpl0 = psys[0].x;
    ypl0 = psys[0].y;

    for(i=0;i<subcycles;i++) {
        move_planet_step(dt*dt_frac);
    }
    xpl1 = psys[0].x;
    ypl1 = psys[0].y;
    
    d2 = sqrt(xpl1*xpl1 + ypl1*ypl1);
    d1 = sqrt(xpl0*xpl0+ypl0*ypl0);
    cross = xpl0*ypl1-xpl1*ypl0;
    omf_new = asin(cross/(d1*d2))/dt;
    domega = omf_new - omf;
    omf = omf_new;
    rotate_sys(omf*dt);

    i = j = k = 0;
    double res=0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res = 0;
            for(i=0;i<size_x;i++) {
                vx[l] -= domega * ymed(j);
                res += vx[l];
            }
            res /=(double)nx;
            vxbar[j] = res;
        }
    }


    return;
}
void init_rk5(void) {

    k1 = (double *)malloc(sizeof(double)*6*nb);
    k2 = (double *)malloc(sizeof(double)*6*nb);
    k3 = (double *)malloc(sizeof(double)*6*nb);
    k4 = (double *)malloc(sizeof(double)*6*nb);
    k5 = (double *)malloc(sizeof(double)*6*nb);
    k6 = (double *)malloc(sizeof(double)*6*nb);
    q0 = (double *)malloc(sizeof(double)*6*nb);
    q1 = (double *)malloc(sizeof(double)*6*nb);

    int i;
    for(i=0;i<6*nb;i++) {
        k1[i] = 0;
        k2[i] = 0;
        k3[i] = 0;
        k4[i] = 0;
        k5[i] = 0;
        k6[i] = 0;
        q0[i] = 0;
        q1[i] = 0;
    }
    return;
}
void free_rk5(void) {

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
    free(q0);
    free(q1);
    return;
}
void move_to_com(void) {
    psys[1].x = -psys[0].x*psys[0].mp/psys[1].mp;
    psys[1].y = -psys[0].y*psys[0].mp/psys[1].mp;
    psys[1].vx = -psys[0].vx*psys[0].mp/psys[1].mp;
    psys[1].vy = -psys[0].vy*psys[0].mp/psys[1].mp;

    return;
}
void read_planet_file(int n, char *directory) {
    char filename[512],filename2[512];
    FILE *f,*f1;
    int i,j;
    double xpl,ypl,zpl,vxpl,vypl,vzpl,mpl,tpl,omfpl;
    double epl,apl,vpl,psipl,phipl,ipl,wpl,alphapl;
    int scanned = 0;

    for(j=0;j<nb;j++) {
        sprintf(filename,"%splanet%d.dat",directory,j);
        f = fopen(filename,"r");
        while ( (scanned = fscanf(f,"%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    &i,&xpl,&ypl,&zpl,&vxpl,&vypl,&vzpl,&mpl,&tpl,&omfpl)) != EOF) {

                if (i == n) {
                    psys[j].x = xpl; 
                    psys[j].y = ypl; 
                    psys[j].z = zpl; 
                    psys[j].vx = vxpl; 
                    psys[j].vy = vypl; 
                    psys[j].vz = vzpl; 
                    psys[j].mp = mpl; 
                    psys[j].omf = omfpl; 
                    psys[j].t = tpl;
                    break;
                }

        }


        fclose(f);
        sprintf(filename2,"%sorbit%d.dat",directory,j);
        f1 = fopen(filename2,"r");

        while ( (scanned = fscanf(f1,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    &tpl,&epl,&apl,&mpl,&vpl,&psipl,&phipl,&ipl,&wpl,&alphapl)) != EOF) {

                if (fabs(tpl-psys[j].t)<1e-8) {
                    psys[j].orbit.e = epl; 
                    psys[j].orbit.a = apl; 
                    psys[j].orbit.M = mpl; 
                    psys[j].orbit.V = vpl; 
                    psys[j].orbit.psi = psipl; 
                    psys[j].orbit.phi = phipl; 
                    psys[j].orbit.i = ipl; 
                    psys[j].orbit.w = wpl; 
                    psys[j].orbit.alpha = alphapl; 
                    break;
                }

        }
        fclose(f1);
    }
    omf = psys[0].omf;


    return;

}
