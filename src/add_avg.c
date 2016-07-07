//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void AddAvgs_cpu () {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Vx);
  INPUT(Vy);
  INPUT(VxMedSS);
  INPUT(VyMedSS);
  INPUT(RhoMedSS);
  INPUT(VxMed);
  INPUT(VyMed);
  INPUT(RhoMed);
//<\USER_DEFINED>


//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* vx_avg = VxMed->field_cpu;
  real* vy_avg = VyMed->field_cpu;
  real* dens_avg = RhoMed->field_cpu;
  real* vx_avg_ss = VxMed->field_cpu;
  real* vy_avg_ss = VyMed->field_cpu;
  real* dens_avg_ss = RhoMed->field_cpu;

  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  int ll2d;
//<\INTERNAL>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
    ll2d = l2d;
    dens[ll] += dens_avg_ss[ll2d] - dens_avg[ll2d];
    vx[ll] += vx_avg_ss[ll2d] - vx_avg[ll2d];
    vy[ll] += vy_avg_ss[ll2d] - vy_avg[ll2d];
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
   
