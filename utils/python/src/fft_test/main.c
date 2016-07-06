#include "evolve.h"

#define l (i + j*nx)

int main(void) {

    int i,j;
    ny = 4;
    nx = 10;
    size_y = ny;
    size_x = nx;
    if (size_x % 2 == 0) {
        NEVEN = TRUE;
    }
    else {
        NEVEN = FALSE;
    }

    allocate_conv();

    double *fld1  = (double *)malloc(sizeof(double)*size_x*size_y);
    double *fld2  = (double *)malloc(sizeof(double)*size_x*size_y);
    double *res  = (double *)malloc(sizeof(double)*size_y*(MMAX+2));
    double *ans  = (double *)malloc(sizeof(double)*size_y*(MMAX+2));
    double *fhat  = (double *)malloc(sizeof(double)*size_x);
    double *fac  = (double *)malloc(sizeof(double)*size_y);

    for(i=0;i<size_y*(MMAX+2);i++) {
        res[i] = 0;
    }
    for(i=0;i<size_y;i++) {
        fac[i] = 1.0;
    }


    
    FILE *f = fopen("input1.dat","r");

    fread(fld1,sizeof(double),size_x*size_y,f);
    fclose(f);
    f = fopen("input2.dat","r");

    fread(fld2,sizeof(double),size_x*size_y,f);
    fclose(f);





    i = 0;
    double temp;
    for(j=0;j<size_y;j++) {
        i=0;
        convolution(&fld1[l],&fld2[l],res,1.0,j,size_y);
        temp = 0;
        for(i=0;i<size_x;i++) {
            temp += fld1[i+j*nx]*fld2[i+j*nx];
        }
        temp /= (double)nx;
        res[j] = temp;
    }
    f = fopen("answer.dat","r");
    fread(ans,sizeof(double),size_y*(MMAX+2),f);
    fclose(f);
    

    double err = 0;
    double tol = 1e-5;

    int mi;
    for(mi=0;mi<MMAX+2;mi++) {
        for(j=0;j<size_y;j++) {
            err += fabs((res[j + size_y*mi] - ans[j+size_y*mi]));
        }
    }
    if (err < tol) {
        printf("convolution PASSED! \n");
        printf("Error %.5e\n",err);
    }
    else {
        printf("Failed, error: %.15e\n",err);

        for(mi=0;mi<MMAX+2;mi++) {
            printf("m = %d\n",mi);
            for(j=0;j<ny;j++) {
                printf("%.16f\t%.16f\n",res[j + mi*size_y],ans[j+mi*size_y]);
            }
            printf("\n\n");
        }
    }

    for(i=0;i<size_y*(MMAX+2);i++) res[i] = 0;

    for(j=0;j<size_y;j++) {
        i=0;
        temp = 0;
        for(i=0;i<size_x;i++) {
            temp += fld1[i+j*nx]*fld2[i+j*nx];
        }
        temp /= (double)nx;
        res[j] = temp;
    }
    convolution_2d(fld1,fld2,res,fac,size_y); 
    err = 0;

    for(mi=0;mi<MMAX+2;mi++) {
        for(j=0;j<size_y;j++) {
            err += fabs((res[j + size_y*mi] - ans[j+size_y*mi]));
        }
    }
    if (err < tol) {
        printf("convolution_2d PASSED! \n");
        printf("Error %.5e\n",err);
    }
    else {
        printf("Failed, error: %.15e\n",err);

        for(mi=0;mi<MMAX+2;mi++) {
            printf("m = %d\n",mi);
            for(j=0;j<ny;j++) {
                printf("%.16f\t%.16f\n",res[j + mi*size_y],ans[j+mi*size_y]);
            }
            printf("\n\n");
        }
    }


    for(i=0;i<size_y*(MMAX+2);i++) res[i] = 0;

    for(j=0;j<size_y;j++) {
        i=0;
        convolution_deriv(&fld1[l],&fld2[l],res,1.0,j,size_y);
        for(i=1;i<MMAX+2;i++) {
            res[j] += res[j + i*size_y];
        }
    }
    f = fopen("answerd.dat","r");
    fread(ans,sizeof(double),size_y*(MMAX+2),f);
    fclose(f);
    

    err = 0;
    for(mi=1;mi<MMAX+2;mi++) {
        for(j=0;j<size_y;j++) {
                err += fabs((res[j + size_y*mi] - ans[j+size_y*mi]));

        }
    }
    if (err < tol) {
        printf("convolution_deriv PASSED!\n");
        printf("Error %.5e\n",err);
    }
    else {
        printf("Failed, error: %.15e\n",err);

        for(mi=0;mi<MMAX+2;mi++) {
            printf("m = %d\n",mi);
            for(j=0;j<ny;j++) {
                printf("%.16f\t%.16f\n",res[j + mi*size_y],ans[j+mi*size_y]);
            }
            printf("\n\n");
        }
    }


    for(i=0;i<size_y*(MMAX+2);i++) res[i] = 0;

    convolution_deriv_2d(fld1,fld2,res,fac,size_y);
    for(j=0;j<size_y;j++) {
        for(i=1;i<MMAX+2;i++) {
            res[j] += res[j + i*size_y];
        }
    }

    err = 0;
    for(mi=1;mi<MMAX+2;mi++) {
        for(j=0;j<size_y;j++) {
                err += fabs((res[j + size_y*mi] - ans[j+size_y*mi]));

        }
    }
    if (err < tol) {
        printf("convolution_deriv_2d PASSED!\n");
        printf("Error %.5e\n",err);
    }
    else {
        printf("Failed, error: %.15e\n",err);

        for(mi=0;mi<MMAX+2;mi++) {
            printf("m = %d\n",mi);
            for(j=0;j<ny;j++) {
                printf("%.16f\t%.16f\n",res[j + mi*size_y],ans[j+mi*size_y]);
            }
            printf("\n\n");
        }
    }

    free_conv();
    free(ans);
    free(res);
    free(fld1);
    free(fld2);
    free(fhat);

    return 1;
}




