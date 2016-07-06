#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define DIMSLENGTH 14
#define DIMSNXINDEX 6
#define NGHY 3
#define NGHZ 3
#define TRUE 1
#define FALSE 0
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define l ( i + pitch*j + stride*k)
#define ls ( i + pitch_s*j + stride_s*k)
#define lxp (((i)<(nx-1)) ? ((l)+1) : ((l)-(nx-1)))
#define lxm (((i)>0) ? ((l)-1) : ((l)+(nx-1)))
#define lyp ((l)+nx)
#define lym ((l)-nx)

#define lzp (l+stride)
#define lzm (l-stride)


int nx,ny,nz;
int nx_s, ny_s, nz_s;
int  pitch, stride;
int  pitch_s, stride_s;

double *xmed, *ymed, *zmed;
double *xm, *ym, *zm;

double *xmed2d, *ymed2d;
double *xm2d, *ym2d;
double *fld, *fld2d;

double YMIN,YMAX;
double ZMIN;
double h = .05;


void read_dims(char *directory);
void dump_domain(char *directory);
void read_domain(char *directory);
void interpolate(char *name, char *directory, int n);
double bilinear( double x, double y);
double scaleH(double x);

int main(int argc, char *argv[]) {
    int n = atoi(argv[1]);
    printf("%d\t%s\n",n,argv[2]);
    
    nx_s = atoi(argv[2]);
    ny_s = atoi(argv[3]);
    nz_s = atoi(argv[4]);

    read_dims(argv[5]);
    nz = 1;
    pitch = nx; stride = nx*ny;
    pitch_s = nx_s; stride_s = nx_s*ny_s;

    xm = (double *)malloc(sizeof(double)*(nx_s+1));
    ym = (double *)malloc(sizeof(double)*(ny_s+1));
    zm = (double *)malloc(sizeof(double)*(nz_s+1));
    xmed = (double *)malloc(sizeof(double)*(nx_s));
    ymed = (double *)malloc(sizeof(double)*(ny_s));
    zmed = (double *)malloc(sizeof(double)*(nz_s));

    xm2d = (double *)malloc(sizeof(double)*(nx+1));
    ym2d = (double *)malloc(sizeof(double)*(ny+1));
    xmed2d = (double *)malloc(sizeof(double)*(nx));
    ymed2d = (double *)malloc(sizeof(double)*(ny));

    fld = (double *)malloc(sizeof(double)*nx_s*ny_s*nz_s);
    fld2d = (double *)malloc(sizeof(double)*(nx)*(ny+2));

    printf("Detected grid with nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
    printf("Converting to  grid with nx = %d, ny = %d, nz = %d\n",nx_s,ny_s,nz_s);

    read_domain(argv[5]);

    YMIN = ym2d[0];
    YMAX = ym2d[ny];
    ZMIN = atan(.25);

    printf("Domain size, r=(%lg,%lg), Delta z = %lg\n",YMIN,YMAX,ZMIN);

    dump_domain(argv[5]);



    printf("Converting density\n");
    interpolate("dens",argv[5],n);
    printf("expanding energy\n");
    interpolate("energy",argv[5],n);



    free(xm); free(xmed);
    free(ym); free(ymed);
    free(zm); free(zmed);
    free(xm2d); free(xmed2d);
    free(ym2d); free(ymed2d);
    free(fld); free(fld2d); 
    return 1;
}

double scaleH(double x) {
    return h*x;
}


void interpolate(char *name, char *directory, int n) {
    int i,j,k;

    int ii,jj;
    int iR, iL, jR, jL;

    double phi,R;
    double zval,norm;

    FILE *f;
    char fname[512];

    sprintf(fname,"%sgas%s%d.dat",directory,name,n);

    f = fopen(fname,"r");

    for(j=1;j<ny+1;j++) {
        fread(&fld2d[nx*j],sizeof(double),nx,f);
    }
    fclose(f);

    for(i=0;i<nx;i++) {
        fld2d[i] = fld2d[i + nx];
    }
    for(i=0;i<nx;i++) {
        fld2d[i + nx*(ny+1)] = fld2d[i + nx*(ny)];
    }




    for(k=0;k<nz_s;k++) {
        for(j=0;j<ny_s;j++) {
            for(i=0;i<nx_s;i++) {
                R = sin(zmed[k]) * ymed[j];
                phi = xmed[i];
                
                for(ii=0;ii<nx;ii++) {
                    if ( xmed2d[ii]>= phi ) {
                        iR = ii;
                        if (xmed2d[ii] == phi) {
                            iL = ii;
                        }
                        else {
                            if (iR == 0) {
                                iL = nx-1;
                            }
                            else {
                                iL = iR - 1;
                            }
                        }
                        break;
                    }
                }
                
                if (R < ymed2d[0]) {
                    jR = 1;
                    jL = 0;
                }
                else if (R > ymed2d[ny-1]) {
                    jL = ny-1;
                    jR = ny;
                }
                else {
                    for(jj=0;jj<ny;jj++) {
                        if ( ymed2d[jj]>= R ) {
                            jR = jj+1;
                            if (ymed2d[jj] == R) {
                                jL = jR;
                            }
                            else {
                                jL = jR-1;
                            }
                            break;
                        }
                    }
                }

                zval = 1/tan(zmed[k]);
                zval *= zval/2;
                norm = sqrt(2*M_PI)*scaleH(R);

                fld[ls] = (fld2d[iL + jL*nx] + fld2d[iL + jR*nx] + fld2d[iR + jL*nx] + fld2d[iR + jR*nx])*.25 * exp(-zval)/norm;

            }
        }
    }

    f = fopen(fname,"w");
    fwrite(fld,sizeof(double),nx_s*ny_s*nz_s,f);
    fclose(f);

    return;
}


/*
void expand_field( char *name, char *directory, int n) {
    FILE *f;
    char fname[512];
    int i,j,k;
    double R,phi,temp,zval,norm;

    sprintf(fname,"%sgas%s%d.dat",directory,name,n);

    f = fopen(fname,"r");

    fread(&fld2d[nx],sizeof(double),nx*ny,f);
    fclose(f);


    for(i=0;i<nx;i++) {
        fld2d[i] = fld2d[i + nx];
    }
    for(i=0;i<nx;i++) {
        fld2d[i+nx*ny] = fld2d[i + nx*(ny-1)];
    }

#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(i,j,k,R,phi,temp,zval,norm)
#endif
    for(k=0;k<nz_s;k++) {
        for(j=0;j<ny_s;j++) {
            for(i=0;i<nx_s;i++) {
                R = ymed[j]*sin(zmed[k]);
                phi = xmed[i];
                printf("Interpolating to %lg\t%lg\n",R,phi);

                temp = bilinear(phi,R);
                printf("Result %lg\n",temp);

                zval = 1/tan(zmed[k]);
                zval *= zval/2;
                norm = sqrt(2*M_PI)*scaleH(R);

                fld[i + nx_s*j + nx_s*ny_s*k] = temp * exp(-zval)/norm;

            }
        }
    }


    f = fopen(fname,"w");

    fwrite(fld,sizeof(double),nx_s*ny_s*nz_s,f);
    fclose(f);

    return;
}
*/
double bilinear(double x, double y) {
    int i,j,k,ii,jj;
    double xL, xR, yL, yR;
    double tempx, tempy,res;

    i = j = k = 0;

    i = nx-2;
    for(ii=0;ii<nx;ii++) {
        if (xmed2d[ii] >= x) {
            i = ii-1;
            break;
        }
    }
    if (i == nx-1) {
        xL = xmed2d[i];
        xR = xmed2d[i] + xmed2d[i]-xmed2d[i-1];
    }
    else {
        if (i == -1) {
            i = nx-1;
            xL = xmed2d[0];
            xR = xmed2d[i+1];
        }

    }

    j = ny-2;
    for(jj=0;jj<ny;jj++) {
        if (ymed2d[jj] >= y) {
            j = jj-1;
            break;
        }
    }
    if (j==ny-1) {
        yL = ymed2d[j];
        yR = ymed2d[j] + ymed2d[j]-ymed2d[j-1];
    }
    else {
        if (j == -1) {
            j = 0;
            yL = ymed2d[0]-(ymed2d[1]-ymed2d[0]);
            yR = ymed2d[0];
        }
        else {
            yL = ymed2d[j];
            yR = ymed2d[j+1];
        }
    }

    tempx  = fld2d[l]*(xR-x)/(xR-xL) + fld2d[lxp]*(x - xL)/(xR-xL); 
    tempy = fld2d[lyp]*(xR-x)/(xR-xL) + fld2d[lxp + pitch]*(x-xL)/(xR-xL);


    res = tempx*(yR-y)/(yR-yL) + tempy*(y-yL)/(yR-yL);
    printf("\n%lg %lg %lg %lg\n",xL,xR,yL,yR);
    printf("%lg %lg %lg %lg\n",fld2d[l],fld2d[lxp],fld2d[lyp],fld2d[lxp+pitch]);
    printf("%lg\n",res);

    return res;
}


void dump_domain(char *directory) {
    FILE *fx,*fy,*fz;
    char fnamex[512];
    char fnamey[512];
    char fnamez[512];
    int i,j,k;
    double temp;

    sprintf(fnamex,"%sdomain_x.dat",directory);
    fx = fopen(fnamex,"w");
    for(i=0;i<nx_s+1;i++) {
        xm[i] = i*2*M_PI/nx_s - M_PI;
        fprintf(fx,"%.16f\n",xm[i]);
    }
    for(i=0;i<nx_s;i++) {
        xmed[i] = .5*(xm[i] + xm[i+1]);
    }
    fclose(fx);

    sprintf(fnamey,"%sdomain_y.dat",directory);
    fy = fopen(fnamey,"w");
    for(j=-NGHY;j<ny_s+1+NGHY;j++) {
        temp  = exp(log(YMIN) + j*log(YMAX/YMIN)/ny_s);
        fprintf(fy,"%.16f\n",temp);
        if ( (j>=0) && (j < ny_s + 1)) {
            ym[j] = temp;
        }
    }
    for(j=0;j<ny_s;j++) {
        ymed[j] = .5*(ym[j] + ym[j+1]);
    }
    fclose(fy);

    sprintf(fnamez,"%sdomain_z.dat",directory);
    fz = fopen(fnamez,"w");
    for(k=-NGHZ;k<nz_s+1+NGHZ;k++) {
        temp = M_PI/2 + -ZMIN + k*ZMIN/(nz_s/2);
        fprintf(fz,"%.16f\n",temp);
        if ( (k>=0) && (k < nz_s + 1)) {
            zm[k] = temp;
        }
    }
    for(k=0;k<nz_s;k++) {
        zmed[k] = .5*(zm[k] + zm[k+1]);
    }
    fclose(fz);
    return;
}

void read_domain(char *directory) {
    FILE *f;
    char fname[512];
    int i,j;
    double temp;

    sprintf(fname,"%sdomain_y.dat",directory);

    f = fopen(fname,"r");
    if (f == NULL) printf("Error finding grid file %s\n",fname);

    for(j=0;j<ny+1 +NGHY;j++) {
        fscanf(f,"%lg\n",&temp);
        if (j>=NGHY) {
            ym2d[j-NGHY] = temp;
        }
    }
    fclose(f);


    for(i=0;i<nx+1;i++) {
        xm2d[i] = i*2*M_PI/(float)nx - M_PI;
    }
    for(j=0;j<ny;j++) {
        ymed2d[j] = .5*(ym2d[j] + ym2d[j+1]);
    }
    for(i=0;i<nx;i++) {
        xmed2d[i] = .5*(xm2d[i] + xm2d[i+1]);
    }

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
    }
    fclose(f);
    return;
}



