#include "fargo3d.h"

void init_stockholm() {
  
  static boolean init = TRUE;

  boolean error_density = TRUE;
  boolean error_vx = TRUE;
  boolean error_vy = TRUE;
  boolean error_vz = TRUE;
  boolean error_energy = TRUE;
  
  if(init) {

  INPUT(Density);
  OUTPUT2D(Density0);
#ifdef ADIABATIC
  INPUT(Energy);
  OUTPUT2D(Energy0);
#endif
#ifdef X
  INPUT(Vx);
  OUTPUT2D(Vx0);
#endif
#ifdef Y
  INPUT(Vy);
  OUTPUT2D(Vy0);
#endif
#ifdef Z
  INPUT(Vz);
  OUTPUT2D(Vz0);
#endif

  int i,j,k;

  i = j = k = 0;
  
#ifdef X
  real* vx  = Vx->field_cpu;
  real* vx0 = Vx0->field_cpu;
#endif
#ifdef Y
  real* vy  = Vy->field_cpu;
  real* vy0 = Vy0->field_cpu;
#endif
#ifdef Z
  real* vz  = Vz->field_cpu;
  real* vz0 = Vz0->field_cpu;
#endif
#ifdef ADIABATIC
  real* e    = Energy->field_cpu;
  real* e0   = Energy0->field_cpu;
#endif
  real* rho  = Density->field_cpu;
  real* rho0 = Density0->field_cpu;

  if ((Restart == YES) || (Restart_Full == YES)) {
    error_density = Read2D(Density0, "density0_2d.dat", OUTPUTDIR, GHOSTINC);
    error_vx = Read2D(Vx0, "vx0_2d.dat", OUTPUTDIR, GHOSTINC);
    error_vy = Read2D(Vy0, "vy0_2d.dat", OUTPUTDIR, GHOSTINC);
    error_vz = Read2D(Vz0, "vz0_2d.dat", OUTPUTDIR, GHOSTINC);
#ifdef ADIABATIC
    error_energy = Read2D(Energy0, "e0_2d.dat", OUTPUTDIR, GHOSTINC);
#endif
  }
  
#ifdef Z
    for (k=0; k<Nz+2*NGHZ; k++) {
#endif
#ifdef Y
      for (j=0; j<Ny+2*NGHY; j++) {
#endif
#ifdef ADIABATIC
	if (error_energy)
	  e0[l2D]   = e[l];
#endif
#ifdef X
	if (error_vx)
	  vx0[l2D]  = vx[l];
#endif
#ifdef Y
	if (error_vy)
	  vy0[l2D]  = vy[l];
#endif
#ifdef Z
	if (error_vz)
	  vz0[l2D]  = vz[l];
#endif
	if (error_density)
	  rho0[l2D] = rho[l];

#ifdef Y
      }
#endif
#ifdef Z
    }
#endif
    Write2D(Density0, "density0_2d.dat", OUTPUTDIR, GHOSTINC);
#ifdef X
    Write2D(Vx0,  "vx0_2d.dat", OUTPUTDIR, GHOSTINC);
#endif
#ifdef Y
    Write2D(Vy0,  "vy0_2d.dat", OUTPUTDIR, GHOSTINC);
#endif
#ifdef Z
    Write2D(Vz0,  "vz0_2d.dat", OUTPUTDIR, GHOSTINC);
#endif
#ifdef ADIABATIC
    Write2D(Energy0,   "e0_2d.dat", OUTPUTDIR, GHOSTINC);
#endif
    init = FALSE;
  }
}
void init_stockholm_accretion() {

  INPUT(Density);
  INPUT(Vx);
  INPUT(Vy);
  OUTPUT2D(Density0_avg);
  OUTPUT2D(Vx0_avg);
  OUTPUT2D(Vy0_avg);
  FARGO_SAFE(ComputeVxmed(Vx));
  FARGO_SAFE(ComputeVymed(Vy));
  FARGO_SAFE(ComputeRhomed(Density));

}
/*
  INPUT(Density);
  OUTPUT2D(Density0);
  OUTPUT2D(Density0_avg);
#ifdef ADIABATIC
  INPUT(Energy);
  OUTPUT2D(Energy0);
  OUTPUT2D(Energy0_avg);
#endif
#ifdef X
  INPUT(Vx);
  OUTPUT2D(Vx0);
  OUTPUT2D(Vx0_avg);
#endif
#ifdef Y
  INPUT(Vy);
  OUTPUT2D(Vy0);
  OUTPUT2D(Vy0_avg);
#endif
#ifdef Z
  INPUT(Vz);
  OUTPUT2D(Vz0);
  OUTPUT2D(Vz0_avg);
#endif

  int i,j,k;

  int size_x = Nx + 2*NGHX;
  i = j = k = 0;
  
#ifdef X
  real* vx  = Vx->field_cpu;
  real* vx0 = Vx0->field_cpu;
  real* vx0_avg = Vx0_avg->field_cpu;
  real vx_temp;
#endif
#ifdef Y
  real* vy  = Vy->field_cpu;
  real* vy0 = Vy0->field_cpu;
  real* vy0_avg = Vy0_avg->field_cpu;
  real vy_temp;
#endif
#ifdef Z
  real* vz  = Vz->field_cpu;
  real* vz0 = Vz0->field_cpu;
  real* vz0_avg = Vz0_avg->field_cpu;
  real vz_temp;
#endif
#ifdef ADIABATIC
  real* e    = Energy->field_cpu;
  real* e0   = Energy0->field_cpu;
  real* e0_avg   = Energy0_avg->field_cpu;
  real e_temp;
#endif
  real* rho  = Density->field_cpu;
  real* rho0 = Density0->field_cpu;
  real* rho0_avg = Density0_avg->field_cpu;
  real rho_temp;

   
#ifdef Z
    for (k=0; k<Nz+2*NGHZ; k++) {
#endif
#ifdef Y
      for (j=0; j<Ny+2*NGHY; j++) {
#endif

#ifdef ADIABATIC
          e_temp   = 0;
#endif
#ifdef X
          vx_temp = 0;
#endif
#ifdef Y
          vy_temp  = 0;
#endif
#ifdef Z
          vz_temp  = 0;
#endif
          rho_temp = 0;
#ifdef X
      for (i=0; i<size_x; i++) {
#endif

#ifdef ADIABATIC
          e_temp += e[l];
#endif
#ifdef X
          vx_temp += vx[l];
#endif
#ifdef Y
          vy_temp += vy[l];
#endif
#ifdef Z
          vz_temp += vz[l];
#endif
          rho_temp += rho[l];
#ifdef X
      }
#endif
#ifdef ADIABATIC
	  e0_avg[l2D]   = e_temp/size_x;
#endif
#ifdef X
	  vx0_avg[l2D]  =  vx_temp/size_x;
#endif
#ifdef Y
	  vy0_avg[l2D]  =  vy_temp/size_x;
#endif
#ifdef Z
	  vz0_avg[l2D]  =  vz_temp/size_x;
#endif
	  rho0_avg[l2D] =  rho_temp/size_x;

#ifdef Y
      }
#endif
#ifdef Z
    }
#endif
}
*/
