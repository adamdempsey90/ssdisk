//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void StockholmBoundary_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT2D(Density0);
  INPUT2D(Density0_avg);
  OUTPUT(Density);
#ifdef ADIABATIC
  INPUT(Energy);
  INPUT2D(Energy0);
  OUTPUT(Energy);
#endif
#ifdef X
  INPUT(Vx);
  INPUT2D(Vx0);
  INPUT2D(Vx0_avg);
  OUTPUT(Vx);
#endif
#ifdef Y
  INPUT(Vy);
  INPUT2D(Vy0);
  INPUT2D(Vy0_avg);
  OUTPUT(Vy);
#endif
#ifdef Z
  INPUT(Vz);
  INPUT2D(Vz0);
  OUTPUT(Vz);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho  = Density->field_cpu;
  real* rho0 = Density0->field_cpu;
  real* rho0_avg = Density0_avg->field_cpu;
#ifdef X
  real* vx  = Vx->field_cpu;
  real* vx0 = Vx0->field_cpu;
  real* vx0_avg = Vx0_avg->field_cpu;
#endif
#ifdef Y
  real* vy  = Vy->field_cpu;
  real* vy0 = Vy0->field_cpu;
  real* vy0_avg = Vy0_avg->field_cpu;
#endif
#ifdef Z
  real* vz  = Vz->field_cpu;
  real* vz0 = Vz0->field_cpu;
#endif
#ifdef ADIABATIC
  real* e    = Energy->field_cpu;
  real* e0   = Energy0->field_cpu;
#endif
  int pitch   = Pitch_cpu;
  int stride  = Stride_cpu;
  int size_x  = Nx+2*NGHX;
  int size_y  = Ny+2*NGHY;
  int size_z  = Nz+2*NGHZ;
  int pitch2d = Pitch2D;
  int nghy = NGHY;
  real y_min = YMIN;
  real y_max = YMAX;
  real z_min = ZMIN;
  real z_max = ZMAX;
  real of    = OMEGAFRAME;
  real of0   = OMEGAFRAME0;
  real r0 = R0;
#ifdef USERWKZ
  real wkzin = WKZIN;
  real wkzout = WKZOUT;
#else
  real wkzin = .0476;
  real wkzout = .19;
#endif
  int periodic_z = PERIODICZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  // Below: we match the standard Stockholm's prescription
  // see De Val Borro et al. (2006), section 3.2
  //  real Y_inf = y_min + (y_max-y_min)*0.0476;
  //  real Y_sup = y_max - (y_max-y_min)*0.19;
//  real Y_inf = exp( log(y_min) + log(y_max/y_min)*0.05);
//  real Y_sup = exp( log(y_max) - log(y_max/y_min)*0.10);
  real Y_inf = y_min + (y_max-y_min)*wkzin;
  real Y_sup = y_max - (y_max-y_min)*wkzout;
  real Z_inf = z_min - (z_max-z_min); // Here we push Z_inf & Z_sup
  real Z_sup = z_max + (z_max-z_min); // out of the mesh
#ifdef CYLINDRICAL
  Z_inf = z_min + (z_max-z_min)*0.1; 
  Z_sup = z_max - (z_max-z_min)*0.1;
  if (periodic_z) { // Push Z_inf & Z_sup out of mesh if periodic in Z
    Z_inf = z_min-r0;
    Z_sup = z_max+r0;
  }
#endif
#ifdef SPHERICAL
  Z_inf = M_PI/2.0-(M_PI/2.0-z_min)*1.2;
  Z_sup = M_PI/2.0+(M_PI/2.0-z_min)*0.8; // Avoid damping in ghost zones
  // if only half upper disk is covered by the mesh
#endif
  real radius;
  real vx0_target;
  real rampy;
  real rampz;
  real rampzz;
  real rampi;
  real ramp;
  real rho_target;
  real vr_target;
  real vp_target;
  real tau;
  real taud;
  real ds = 0.03333;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++) {
#endif
//<#>
	rampy = 0.0;
	rampz = 0.0;
	rampzz = 0.0;
#ifdef Y
	if(ymed(j) > Y_sup) {
	  rampy   = (ymed(j)-Y_sup)/(y_max-Y_sup);
#ifdef NOSTOCKHOLMRIGHT
	  rampy   = 0 ;
#endif
#ifdef STOCKHOLMACC
      vr_target = vy0_avg[l2D];
      vp_target = vx0_avg[l2D];
      rho_target = rho0_avg[l2D];
#else
      vr_target = vy0[l2D];
      vp_target = vx0[l2D];
      rho_target = rho0[l2D];
#endif
	}
	if(ymed(j) < Y_inf) {
	  rampy   = (Y_inf-ymed(j))/(Y_inf-y_min);
      vr_target = vy0[l2D];
      vp_target = vx0[l2D];
      rho_target = rho0[l2D];
/*
#ifdef STOCKHOLMACC
      rho_target = rho0_avg[l2D];
#else
      rho_target = rho0[l2D];
#endif
*/
	}
	rampy *= rampy;		/* Parabolic ramp as in De Val Borro et al (2006) */

#endif
#ifdef Z
	if(zmed(k) > Z_sup) {
	  rampz   = (zmed(k)-Z_sup)/(z_max-Z_sup);
	}
	if(zmed(k) < Z_inf) {
	  rampz   = (Z_inf-zmed(k))/(Z_inf-z_min);
	}
	rampz = rampz * rampz;		/* vertical ramp in X^2 */
	if(zmin(k) > Z_sup) {
	  rampzz  = (zmin(k)-Z_sup)/(z_max-Z_sup);
	}
	if(zmin(k) < Z_inf) {
	  rampzz  = (Z_inf-zmin(k))/(Z_inf-z_min);
	}
	rampzz= rampzz * rampzz;		/* vertical ramp in X^2 */
#endif
	if (periodic_z) {
	  rampz = 0.0;
	  rampzz = 0.0;
	}
	ramp = rampy+rampz;
	rampi= rampy+rampzz;
	tau = ds*sqrt(ymed(j)*ymed(j)*ymed(j)/G/MSTAR);
	if(ramp>0.0) {
	  taud = tau/ramp;
#ifndef NOWAVEKILLRHO
      rho[l] = (rho[l]*taud+rho_target*dt)/(dt+taud);
#endif
#ifndef NOWAVEKILLVPHI
#ifdef X
	  vx0_target = vx0[l2D];
	  radius = ymed(j);
#ifdef SPHERICAL
	  radius *= sin(zmed(k));
#endif
	  vx0_target -= (of-of0)*radius;
	  vx[l] = (vx[l]*taud+vx0_target*dt)/(dt+taud);
#endif
#endif
#ifdef Y
	  vy[l] = (vy[l]*taud+vr_target*dt)/(dt+taud);
#endif
	}
#ifdef Z
	if(rampi>0.0) {
	  taud = tau/rampi;
	  vz[l] = (vz[l]*taud+vz0[l2D]*dt)/(dt+taud);
	}
#endif
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
