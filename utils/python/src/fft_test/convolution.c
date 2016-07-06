#include "evolve.h"
#include "rfftw.h"


int size_f;
double fft_norm_fac;
rfftw_plan fplan;

/*
void fft1d(const double *fld1, double *res) {
    memcpy(&in1[0],&fld1[0],sizeof(double)*size_x);
    rfftw_one(fplan,in1,out1);
    memcpy(&res[0],&out1[0],sizeof(double)*size_x);
    return;
}
*/

void convolution_2d( double *fld1,  double *fld2, double *res, double *fac, int ncols) {
    fftw_real *out1, *out2;
    int mi;
    double temp;


    out1 = (fftw_real *)malloc(sizeof(fftw_real)*size_x*ncols);
    out2 = (fftw_real *)malloc(sizeof(fftw_real)*size_x*ncols);

    rfftw(fplan,ncols,fld1,1,size_x,out1,1,size_x);
    rfftw(fplan,ncols,fld2,1,size_x,out2,1,size_x);
    

    int j;

    for(j=0;j<ncols;j++) {

        res[j + ncols] += fac[j]*out1[0 + j*size_x]*out2[0 + j*size_x]/fft_norm_fac; 
    
        for(mi=1;mi<MMAX;mi++) {
            res[j + (mi+1)*ncols] += fac[j]*2*(out1[mi + j*size_x]*out2[mi + j*size_x] + out1[size_x-mi + j*size_x]*out2[size_x-mi + j*size_x])/fft_norm_fac;
        }
        temp = 0;
        for(mi=MMAX; mi<size_f;mi++) {
            temp += 2*(out1[mi + j*size_x]*out2[mi + j*size_x] + out1[size_x-mi + j*size_x]*out2[size_x-mi + j*size_x]);
        }
        if ( NEVEN ) {
            temp += out1[size_x/2 + j*size_x]*out2[size_x/2 + j*size_x];
        }
        res[j + (MMAX+1)*ncols] += fac[j]*temp/fft_norm_fac;
    }
    

    fftw_free(out1);
    fftw_free(out2);
    return;
}

void convolution( double *fld1,  double *fld2, double *res, double fac, int jres, int ncols) {
    fftw_real *out1, *out2;
    int mi;
    double temp;


    out1 = (fftw_real *)malloc(sizeof(fftw_real)*size_x);
    out2 = (fftw_real *)malloc(sizeof(fftw_real)*size_x);
    


    rfftw_one(fplan,fld1,out1);
    rfftw_one(fplan,fld2,out2);

    res[jres + ncols] += fac*out1[0]*out2[0]/fft_norm_fac; 
    
    for(mi=1;mi<MMAX;mi++) {
        res[jres + (mi+1)*ncols] += fac*2*(out1[mi]*out2[mi] + out1[size_x-mi]*out2[size_x-mi])/fft_norm_fac;
    }
    temp = 0;
    for(mi=MMAX; mi<size_f;mi++) {
        temp += 2*(out1[mi]*out2[mi] + out1[size_x-mi]*out2[size_x-mi]);
    }
    if ( NEVEN ) {
        temp += out1[size_x/2]*out2[size_x/2];
    }
    res[jres + (MMAX+1)*ncols] += fac*temp/fft_norm_fac;
    

    fftw_free(out1);
    fftw_free(out2);
    return;
}
void convolution_deriv_2d( double *fld1,  double *fld2, double *res, double *fac, int ncols) {
    int mi;
    fftw_real *out1, *out2;

    out1 = (fftw_real *)malloc(sizeof(fftw_real)*size_x*ncols);
    out2 = (fftw_real *)malloc(sizeof(fftw_real)*size_x*ncols);

    rfftw(fplan,ncols,fld1,1,size_x,out1,1,size_x);
    rfftw(fplan,ncols,fld2,1,size_x,out2,1,size_x);
    

    int j;

    for(j=0;j<ncols;j++) {
        for(mi=1;mi<MMAX;mi++) {
            res[j + (mi+1)*ncols] += mi*fac[j]*2*(out1[mi+j*size_x]*out2[size_x-mi+j*size_x] - out1[size_x-mi+j*size_x]*out2[mi+j*size_x])/fft_norm_fac;
        }
        double temp = 0;
        for(mi=MMAX;mi< size_f;mi++) {
            temp += mi*2*(out1[mi + j*size_x]*out2[size_x-mi + j*size_x] - out1[size_x-mi + j*size_x]*out2[mi + j*size_x]);
        }
        res[j + (MMAX+1)*ncols] += fac[j]*temp/fft_norm_fac;
    }


    fftw_free(out1);
    fftw_free(out2);
    return;
}
void convolution_deriv( double *fld1,  double *fld2, double *res, double fac, int jres,int ncols) {
    int mi;
    fftw_real *out1, *out2;

    out1 = (fftw_real *)malloc(sizeof(fftw_real)*size_x);
    out2 = (fftw_real *)malloc(sizeof(fftw_real)*size_x);

    //memcpy(&in1[0],&fld1[0],sizeof(double)*size_x);
    //memcpy(&in2[0],&fld2[0],sizeof(double)*size_x);
    rfftw_one(fplan,fld1,out1);
    rfftw_one(fplan,fld2,out2);

    for(mi=1;mi<MMAX;mi++) {
        res[jres + (mi+1)*ncols] += mi*fac*2*(out1[mi]*out2[size_x-mi] - out1[size_x-mi]*out2[mi])/fft_norm_fac;
    }
    double temp = 0;
    for(mi=MMAX;mi< size_f;mi++) {
        temp += mi*2*(out1[mi]*out2[size_x-mi] - out1[size_x-mi]*out2[mi]);
    }
    res[jres + (MMAX+1)*ncols] += fac*temp/fft_norm_fac;


    fftw_free(out1);
    fftw_free(out2);
    return;
}

void allocate_conv(void) {
    size_f = (size_x + 1)/2;
    
    fplan = rfftw_create_plan(size_x,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);

    fft_norm_fac = (double)size_x*size_x;

    return;
}
void free_conv(void) {
    rfftw_destroy_plan(fplan);
    return;
}
