//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void visctensor_cyl_ss_cpu(){

//<USER_DEFINED>
  INPUT(DensityMed);
  INPUT(LamDep);
#ifdef ALPHAVISCOSITY
  INPUT(Energy);
#endif
#ifdef X
  INPUT(VxMed);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho = Density->field_cpu;
#ifdef ALPHAVISCOSITY
  real* energy = Energy->field_cpu;
#endif
#ifdef X
  real* vx = Vx->field_cpu;
#endif
  real* lamdep = LamDep->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP; 
  int size_y = Ny+NGHY;
  int size_z = Nz+2*NGHZ-1;
  int begin_j = NGHY;
  real dx = Dx;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  real viscosityp;
  real viscositym;
  real pixyp;
  real pixym;
//<\INTERNAL>

//<CONSTANT>
// real NU(1);
// real GAMMA(1);
// real ALPHA(1);
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvVj(Ny+2*NGHY);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=begin_j; j<size_y; j++) {
#endif
//<#>
#ifdef ALPHAVISCOSITY
#ifdef ISOTHERMAL
	viscosityp= ALPHA*.5*(energy[lyp]*energy[lyp]+energy[l]*energy[l])*sqrt(ymin(j+1)*ymin(j+1)*ymin(j+1)/(G*MSTAR));
	viscositym= ALPHA*.5*(energy[l]*energy[l]+energy[lym]*energy[lym])*sqrt(ymin(j)*ymin(j)*ymin(j)/(G*MSTAR));
#else
	viscosityp= ALPHA*GAMMA*(GAMMA-1.0)*(energy[lyp]+energy[l])/(rho[lyp]+rho[l])*sqrt(ymin(j+1)*ymin(j+1)*ymin(j+1)/(G*MSTAR));
	viscositym= ALPHA*GAMMA*(GAMMA-1.0)*(energy[l]+energy[lym])/(rho[l]+rho[lym])*sqrt(ymin(j)*ymin(j)*ymin(j)/(G*MSTAR));
#endif
#else
	viscositym = viscosity = NU;
#endif

#if defined(Y) && defined(X)

	pixy = viscositym*.25*(rho[l2D]+rho[l2D-1])*( (vx[l2D]-vx[l2D-1])/(ymed(j)-ymed(j-1))-.5*(vx[l2D]+vx[l2D-1])/ymin(j)); //centered on left, inner vertical edge in z
	pixyp = viscosityp*.25*(rho[l2D+1]+rho[l2D])*( (vx[l2D+1]-vx[l2D])/(ymed(j+1)-ymed(j))-.5*(vx[l2D+1]+vx[l2D])/ymin(j+1)); //centered on left, inner vertical edge in z
#endif

	lamdep[l2D] += (ymin(j+1)*ymin(j+1)*pixyp - ymin(j)*ymin(j)*pixy)/((ymin(j+1)-ymin(j))*ymed(j));
//<\#>
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
