#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif




#define NGHY 3
#define TRUE 1
#define FALSE 0
#define G 1.0
#define MSTAR 1.0
#define MAXSTEPS 1000000
#define CVNR 1.41
#define CVNL 0.05
#define MINDT 1e-10
#define MMAX 30
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define ARTIFICIALVISCOSITY
#define FARGO
#define STOCKHOLMACC
#define NOWAVEKILLRHO


#define MALLOC_SAFE(ptr) if (ptr == NULL) printf("Malloc error at line %d!\n",__LINE__);
//#define FREE_SAFE(ptr) free(ptr); ptr=NULL; 

#define ymed(j) Ymed[(j)]
#define xmed(i) Xmed[(i)]
#define ymin(j) Ymin[(j)]
#define xmin(i) Xmed[(i)]
#define zone_size_x(j,k) (dx*ymed(j))
#define zone_size_y(j,k) (ymin(j+1)-ymin(j))
#define SurfY(j,k) (ymin(j)*dx)
#define SurfX(j,k) (ymin(j+1)-ymin(j))
#define InvVol(j,k) (2/(dx*(ymin(j+1)*ymin(j+1) - ymin(j)*ymin(j))))
#define Vol(j,k) (0.5*(dx*(ymin(j+1)*ymin(j+1) - ymin(j)*ymin(j))))
#define XC ymed(j)*cos(xmed(i))
#define YC ymed(j)*sin(xmed(i))

#define l   ((i)+(j)*(nx)+((k)*stride))
#define l_f(ii,jj,kk)   ((ii)+(jj)*(nx)+((kk)*stride))
#define lact   ((i)+(jact)*(nx)+((k)*stride))
#define lxp (((i)<(nx-1)) ? ((l)+1) : ((l)-(nx-1)))
#define lxm (((i)>0) ? ((l)-1) : ((l)+(nx-1)))
#define l2D ((j)+((k)*pitch2d))
#define l2D_int ((j)+((k)*pitch2d))

#define ixm ((i)>0 ? ((i)-1) : nx-1)
#define ixp ((i)<nx-1 ? ((i)+1) : 0)

#define lyp ((l)+nx)
#define lym ((l)-nx)

#define lzp (l+stride)
#define lzm (l-stride)


typedef struct Parameters {
    double wkzout;
    double xmin;
    double dt;
    double ymin;
    double mdot;
    double flaringindex;
    double xmax;
    double wkzin;
    double ymax;
    double soft;
    double alpha;
    double h;
    int ninterm;
    int ntot;
    int nx;
    int ny;
    int log;
    int corotate;
    int indirect;
    double cfl;
    int nz;
    double mp;
    double a;
    double omf;
    double nuindex;
    double vrindex;

} Parameters;

typedef struct Orbit {

    double a; // semi-major axis
    double e; // eccentricity
    double M; // mean anomaly
    double V; // true anomaly
    double psi; // arg of periastron from ascending ndoe
    double phi; // angle b/w actual and initial position of x-axis
    double i; // inclination
    double w; // longitude of ascending node w.r.t actual x-axis
    double alpha; // projection of perihelion w.r.t actual x-axis

} Orbit;

typedef struct Planet {

    double x,y,z;
    double vx,vy,vz;
    double mp, omf;
    double dist;
    double t;
    Orbit orbit;
} Planet;


double *dens, *vx, *vy, *Pres, *indPot,*Pot, *energy;
double *vmed;
double *qR, *qL;
double *dens0, *vx0, *vy0;
double *dbar,*dbart, *vxbar, *vybar, *dbarstar;
double *mdotavg;
double *conv_prefac;
double *vx_temp, *vy_temp;
double *Pixp, *Pixm, *Piym, *Piyp;
double *slope, *divrho, *denstar, *Qs;
double *Ymed, *Xmed, *Ymin, *Xmin;
double *tauxx, *tauxy, *tauyy, *tauxyavg;
double *Lt, *Ld, *Lw, *drFw, *drFwB, *drFd, *drFdB, *drFt, *Lamdep,*LamdepB, *Lamex;
double *mdotl,*drFnu;
double *LdS, *dtLdS, *dbarS, *dtdbar;
double *LamdepS, *dtLd_rhs;
double *dtLt, *dtLd, *dtLw;
int *nshift;

double dt,omf,dx,time_step,omf0;
double CFL;
int nx, ny, nz, size_x, size_y, size_z,stride,pitch,pitch2d,pitch2d_int,nsteps;
int NEVEN;
int nb;
int IndirectTerm;
int num_threads;

Parameters params;
Planet psys[2];

double Nu(double x);
double Cs(double X);


void set_ang(void);
void allocate_all(void);
void free_all(void);
void compute_Pres(void);
void compute_energy(void);
void viscosity(void);
void potential(void);
void temp_to_vel(void);
void vel_to_temp(void);
void source_step(void);
void transport_step(void);
void set_momenta(void);
void transportY(void);
void DividebyRho(double *q);
void DividebyRhoavg(double *q);
void transportX(double *vxt, int ppa);
void set_vel(void);
void vanleer_y_a(double *q);
void vanleer_y_a_avg(double *q);
void vanleer_x_a(double *q);
void vanleer_y_b(double *q, double *qs, double dt);
void vanleer_y_b_avg(double *q, double *qs, double dt);
void vanleer_x_b(double *q, double *qs, double dt,double *vxt);
void updateX(double *q, double *qs,double dt,double *vxt);
void updateY(double *q, double *qs,double dt);
void update_flux(double *qs,double *q);
void update_flux_avg(double *qs, double *q);
void update_density_X(double dt,double *vxt);
void update_density_Y(double dt);
void set_bc(void);
void ymax_bound(void);
void ymin_bound(void);
void ymin_bound_acc(void);
void ymax_bound_acc(void);
void read_domain(char *directory);
void read_single_file(int n,int i, char *directory);
void read_files(int n, char *directory);
void output(char *directory);
void output_init(char *directory);
void set_Lamdep(void);
void set_waves(void);
void set_Lamex(void);
void output_torque(char *directory,int n);
void set_avg(int p);
double cfl(void);
void read_planet_file(int n, char *directory);
void get_accel(double *q, double *k, double dtn);
void move_planet_step(double dtn);
void rotate_sys(double angle);
void move_planet(void);
void time_avg(void);
void init_rk5(void);
void free_rk5(void);
void artificial_visc(void);
void move_to_com(void);
void read_param_file(char *directory);
void stockholm(void);
void fargo_transport(void);
void vanleer_ppa_b(double dt, double *q, double *qs, double *vxt);
void vanleer_ppa_a(double *q);
void vanleer_ppa(double dt, double *q, double *qs, double *vxt);
void advect_shift(double *q, int *nshift);
void compute_residuals(double dt);
void output_stock(char *directory);
void read_stockholm(char *directory);
void output_psys(char *directory, int n);
void compute_vmed(double *vt);
void convolution( double *fld1,  double *fld2, double *result, double fac, int jres,int ncols);
void convolution_deriv( double *fld1,  double *fld2, double *res, double fac, int jres,int ncols);
void free_conv(void);
void allocate_conv(void);
void convolution_2d( double *fld1,  double *fld2, double *res, double *fac, int jstart, int ncols,int ntrans);
void convolution_deriv_2d( double *fld1,  double *fld2, double *res, double *fac, int jstart, int ncols,int ntrans);
void set_dtLt(double fac);
