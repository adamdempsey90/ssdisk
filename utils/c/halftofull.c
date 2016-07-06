#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif


#define DIMSLENGTH 14
#define DIMSNXINDEX 6
#define NGHZ 3

int nx,ny,nz,stride;
void convert_grid(char *directory);
void convert_field(char *name, char *directory, int n, int flip);
void read_dims(char *directory);

int main(int argc, char *argv[]) {
    
    int n = atoi(argv[1]);
    printf("%d\t%s\n",n,argv[2]);

    read_dims(argv[2]);
    printf("Detected grid with nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
    printf("Converting grid\n");
    convert_grid(argv[2]);
    printf("Converting density\n");
    convert_field("dens",argv[2],n,1);
    printf("Converting energy\n");
    convert_field("energy",argv[2],n,1);
    printf("Converting vx\n");
    convert_field("vx",argv[2],n,1);
    printf("Converting vy\n");
    convert_field("vy",argv[2],n,1);
    printf("Converting vz\n");
    convert_field("vz",argv[2],n,-1);
    printf("Done\n");


    return 1;
}


void convert_grid(char *directory) {
    FILE *f;
    char fname[512];
    int k;
    double *fld;
    double dz,zp;
    sprintf(fname,"%sdomain_z.dat",directory);

    printf("opening %s\n", fname);
    fld = (double *)malloc(sizeof(double)*(nz+1+NGHZ));


    f = fopen(fname,"r");
    if (f == NULL) printf("Error finding grid file %s\n",fname);

    for(k=0;k<nz+1 + NGHZ;k++) {
        fscanf(f,"%lg\n",&fld[k]);
    }
    fclose(f);

    f = fopen(fname,"w");

    for(k=0;k< nz+1 + NGHZ;k++) {
        fprintf(f,"%.16f\n",fld[k]);
    }
    for(k=nz+NGHZ-1;k>=0;k--) {
        dz = fld[k] - fld[nz+NGHZ];
        zp = fld[nz+NGHZ] - dz;
        fprintf(f,"%.16f\n",zp);
    }
    fclose(f);
    free(fld);
    return;
}

void convert_field(char *name, char *directory, int n, int flip) {
    FILE *f;
    char fname[512];
    int i,j,k;
    double *fld;
    sprintf(fname,"%sgas%s%d.dat",directory,name,n);

    fld = (double *)malloc(sizeof(double)*nx*ny*nz*2);


    f = fopen(fname,"r");

    fread(fld,sizeof(double),nx*ny*nz,f);
    fclose(f);
#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=nz;k<nz*2;k++) {
        for(j=0;j<ny;j++) {
            for(i=0;i<nx;i++) {
                fld[i + nx*j + stride*k] = flip*fld[i + nx*j + stride*( nz - 1 - (k-nz))];
            }
        }
    }

    f = fopen(fname,"w");

    fwrite(fld,sizeof(double),nx*ny*nz*2,f);
    fclose(f);

    free(fld);
    return;
}   

void read_dims(char *directory) {
    char fname[512];
    FILE *f;
    char tempstr[512];
    double temp;
    int i;

    sprintf(fname,"%sdimensions.dat",directory);

    f = fopen(fname,"r");

    for(i=0;i<DIMSLENGTH;i++) {
        fscanf(f,"%s\n",tempstr);

    }

    for(i=0;i<DIMSLENGTH;i++) {
        fscanf(f,"%lg\t",&temp);
        if (i == DIMSNXINDEX) {
            nx = (int)temp;
        }
        if (i == DIMSNXINDEX+1) {
            ny = (int)temp;
        }
        if (i == DIMSNXINDEX+2) {
            nz = (int)temp;
        }
    }
    fclose(f);
    stride = nx*ny;
    return;
}
