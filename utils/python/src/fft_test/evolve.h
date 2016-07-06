#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define TRUE 1
#define FALSE 0
#define MMAX 4
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif



int ny,nx,size_x,size_y;
int NEVEN;


void convolution_deriv_2d( double *fld1,  double *fld2, double *res, double *fac, int ncols);
void convolution_2d( double *fld1,  double *fld2, double *res, double *fac,int ncols);
void convolution( double *fld1,  double *fld2, double *res, double fac, int jres, int ncols);
void convolution_deriv( double *fld1,  double *fld2, double *res, double fac, int jres,int ncols);
void free_conv(void);
void allocate_conv(void);
