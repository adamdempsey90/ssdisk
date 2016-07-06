//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Fill_GhostsX_erad_cpu () {

//<USER_DEFINED>
  INPUT(Energyrad);
  OUTPUT(Energyrad);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  int lghost1;
  int lcopy1;
  int lghost2;
  int lcopy2;
//<\INTERNAL>

//<EXTERNAL>
  real* erad = Energyrad->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int nx = Nx;
  int size_x = NGHX; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>

//<MAIN_LOOP>
    i = j = k = 0;
#ifdef Z
    for (k = 0; k < size_z; k++) {
#endif
#ifdef Y
      for (j = 0; j < size_y; j++) {
#endif	 
	for (i = 0; i < size_x; i++) {
//<#>
	  lghost1 = l;
          lcopy1 = l + nx; 
          lcopy2 = lghost1 + NGHX;
          lghost2 = lcopy2 + nx; 
	  
	  erad[lghost1]    = erad[lcopy1];
	  erad[lghost2]    = erad[lcopy2];
//<\#>
	}
#ifdef Y
      }
#endif	 
#ifdef Z
    }
#endif
//<\MAIN_LOOP>
}
