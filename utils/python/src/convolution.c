#include "evolve.h"

#ifdef _FFTW
#ifdef _OPENMP
#include "rfftw_threads.h"
#else
#include "rfftw.h"
#endif

int size_f;
double fft_norm_fac;
rfftw_plan fplan;
#endif

void convolution_2d( double *fld1,  double *fld2, double *res, double *fac, int jstart,int ncols,int ntrans) {
#ifdef _FFTW
    fftw_real *out1, *out2;
    int mi;
    double temp;


    out1 = (fftw_real *)malloc(sizeof(fftw_real)*size_x*ntrans);
    out2 = (fftw_real *)malloc(sizeof(fftw_real)*size_x*ntrans);

#ifdef _OPENMP
    rfftw_threads(num_threads,fplan,ntrans,fld1,1,size_x,out1,1,size_x);
    rfftw_threads(num_threads,fplan,ntrans,fld2,1,size_x,out2,1,size_x);
#else
    rfftw(fplan,ntrans,fld1,1,size_x,out1,1,size_x);
    rfftw(fplan,ntrans,fld2,1,size_x,out2,1,size_x);
#endif
    

    int j;

    for(j=jstart;j<ntrans;j++) {

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
#endif // _FFTW
    return;
}

void convolution( double *fld1,  double *fld2, double *res, double fac, int jres, int ncols) {
#ifdef _FFTW
    fftw_real *out1, *out2;
    int mi;
    double temp;


    out1 = (fftw_real *)malloc(sizeof(fftw_real)*size_x);
    out2 = (fftw_real *)malloc(sizeof(fftw_real)*size_x);
    

#ifdef _OPENMP
    rfftw_threads_one(num_threads,fplan,fld1,out1);
    rfftw_threads_one(num_threads,fplan,fld2,out2);
#else
    rfftw_one(fplan,fld1,out1);
    rfftw_one(fplan,fld2,out2);
#endif

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

#endif // _FFTW
    return;
}
void convolution_deriv_2d( double *fld1,  double *fld2, double *res, double *fac, int jstart,int ncols,int ntrans) {
#ifdef _FFTW
    int mi;
    fftw_real *out1, *out2;

    out1 = (fftw_real *)malloc(sizeof(fftw_real)*size_x*ntrans);
    out2 = (fftw_real *)malloc(sizeof(fftw_real)*size_x*ntrans);

#ifdef _OPENMP
    rfftw_threads(num_threads,fplan,ntrans,fld1,1,size_x,out1,1,size_x);
    rfftw_threads(num_threads,fplan,ntrans,fld2,1,size_x,out2,1,size_x);
#else
    rfftw(fplan,ntrans,fld1,1,size_x,out1,1,size_x);
    rfftw(fplan,ntrans,fld2,1,size_x,out2,1,size_x);
#endif
    

    int j;

    for(j=jstart;j<ntrans;j++) {
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
#endif
    return;
}
void convolution_deriv( double *fld1,  double *fld2, double *res, double fac, int jres,int ncols) {
#ifdef _FFTW
    int mi;
    fftw_real *out1, *out2;

    out1 = (fftw_real *)malloc(sizeof(fftw_real)*size_x);
    out2 = (fftw_real *)malloc(sizeof(fftw_real)*size_x);

    //memcpy(&in1[0],&fld1[0],sizeof(double)*size_x);
    //memcpy(&in2[0],&fld2[0],sizeof(double)*size_x);
#ifdef _OPENMP
    rfftw_threads_one(num_threads,fplan,fld1,out1);
    rfftw_threads_one(num_threads,fplan,fld2,out2);
#else
    rfftw_one(fplan,fld1,out1);
    rfftw_one(fplan,fld2,out2);
#endif

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
#endif // _FFTW
    return;
}

void allocate_conv(void) {
#ifdef _FFTW
    size_f = (size_x + 1)/2;

#ifdef _OPENMP
    int fftw_status = fftw_threads_init();
    if (fftw_status != 0) printf("Problem initializing fftw threads\n"); 
#endif
    fplan = rfftw_create_plan(size_x,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);

    fft_norm_fac = (double)size_x*size_x;

#endif // _FFTW
    return;
}
void free_conv(void) {
#ifdef _FFTW
    rfftw_destroy_plan(fplan);
#endif // _FFTW
    return;
}
