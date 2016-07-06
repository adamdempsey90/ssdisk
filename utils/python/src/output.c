#include "evolve.h"
void output_init(char *directory) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%soutput_init.dat",directory);
    f = fopen(fname,"w");
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fclose(f);
    return;
}
void output_stock(char *directory) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%soutput_stock.dat",directory);
    f = fopen(fname,"w");
    fwrite(&dens0[0],sizeof(double),size_y*size_z,f);
    fwrite(&vy0[0],sizeof(double),size_y*size_z,f);
    fwrite(&vx0[0],sizeof(double),size_y*size_z,f);
    fclose(f);
    return;
}
void output_psys(char *directory,int n) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%spsys.dat",directory);
    if (n==0) {
        f = fopen(fname,"w");
    }
    else {
        f = fopen(fname,"a");
    }
    fprintf(f,"%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",n,
            psys[0].x, psys[0].y, psys[0].vx, psys[0].vy,
            psys[1].x, psys[1].y, psys[1].vx, psys[1].vy);
    fclose(f);
    return;
}
void output(char *directory) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%soutput.dat",directory);
    f = fopen(fname,"w");
    fwrite(&Ymin[0],sizeof(double),size_y+1,f);
    fwrite(&Xmin[0],sizeof(double),size_x+1,f);
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&Pot[0],sizeof(double),size_x*size_y*size_z,f);

    fclose(f);
    return;
}
void output_torque(char *directory,int n) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%storque%d.dat",directory,n);
    f = fopen(fname,"w");
    fwrite(&Ymed[NGHY],sizeof(double),ny,f);
    fwrite(&dbart[NGHY],sizeof(double),ny,f);
    fwrite(&Lt[NGHY],sizeof(double),ny,f);
    fwrite(&Lt[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&Ld[NGHY],sizeof(double),ny,f);
    fwrite(&Ld[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&Lw[NGHY],sizeof(double),ny,f);
    fwrite(&Lw[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&LdS[NGHY],sizeof(double),ny,f);
    fwrite(&LdS[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&drFt[NGHY],sizeof(double),ny,f);
    fwrite(&drFd[NGHY],sizeof(double),ny,f);
    fwrite(&mdotl[NGHY],sizeof(double),ny,f);
    fwrite(&drFnu[NGHY],sizeof(double),ny,f);
    fwrite(&drFdB[NGHY],sizeof(double),ny,f);
    fwrite(&drFw[NGHY],sizeof(double),ny,f);
    fwrite(&drFwB[NGHY],sizeof(double),ny,f);
    fwrite(&Lamex[NGHY],sizeof(double),ny,f);
//    fwrite(&Lamex[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&Lamdep[NGHY],sizeof(double),ny,f);
    fwrite(&LamdepB[NGHY],sizeof(double),ny,f);
    fwrite(&dtLt[NGHY],sizeof(double),ny,f);
    fwrite(&dtLd[NGHY],sizeof(double),ny,f);
    fwrite(&dtLw[NGHY],sizeof(double),ny,f);
    fwrite(&dtLdS[NGHY],sizeof(double),ny,f);
    fwrite(&dtdbar[NGHY],sizeof(double),ny,f);
    fwrite(&mdotavg[NGHY],sizeof(double),ny,f);
    fwrite(&LamdepS[NGHY + size_y*0],sizeof(double),ny,f);
    fwrite(&LamdepS[NGHY + size_y*1],sizeof(double),ny,f);
    fwrite(&LamdepS[NGHY + size_y*2],sizeof(double),ny,f);
    fwrite(&LamdepS[NGHY + size_y*3],sizeof(double),ny,f);
    fwrite(&LamdepS[NGHY + size_y*4],sizeof(double),ny,f);
    fwrite(&LamdepS[NGHY + size_y*5],sizeof(double),ny,f);
    fwrite(&LamdepS[NGHY + size_y*6],sizeof(double),ny,f);
    fwrite(&dtLd_rhs[NGHY],sizeof(double),ny,f);

    fclose(f);

    sprintf(fname,"%storque_m%d.dat",directory,n);
    f = fopen(fname,"w");
    int mi;
    for(mi=0;mi<MMAX+2;mi++) {
        fwrite(&drFd[NGHY+size_y*mi],sizeof(double),ny,f);
    }
    for(mi=0;mi<MMAX+2;mi++) {
        fwrite(&Lamex[NGHY+size_y*mi],sizeof(double),ny,f);
    }
    for(mi=0;mi<MMAX+2;mi++) {
        fwrite(&dtLt[NGHY+size_y*mi],sizeof(double),ny,f);
    }


    fclose(f);


    return;
}
