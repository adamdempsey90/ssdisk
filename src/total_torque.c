//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

real total_torque_cpu () {

//<USER_DEFINED>
  INPUT(Density);
  real rplanet = sqrt(Xplanet*Xplanet+Yplanet*Yplanet+Zplanet*Zplanet);
  real rsmoothing = THICKNESSSMOOTHING*ASPECTRATIO*pow(rplanet/R0,FLARINGINDEX)*rplanet;
//<\USER_DEFINED>


//<EXTERNAL>
  real* dens = Density->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real rsm2 = rsmoothing*rsmoothing;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real dx;
  real dy;
  real dz=0.0;
  real InvDist3;
  real cellmass;
  real dist2;
  real distance;
  real fxi;
  real fyi;
  real tottorq;
  real Gtottorq;
//<\INTERNAL>

//<CONSTANT>
// real Xplanet(1);
// real Yplanet(1);
// real Zplanet(1);
// real VXplanet(1);
// real VYplanet(1);
// real VZplanet(1);
// real MplanetVirtual(1);
// real Syk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;
tottorq = 0;
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
	cellmass = Vol(j,k)*dens[ll];
#ifdef CARTESIAN
	dx = xmed(i)-Xplanet;
	dy = ymed(j)-Yplanet;
#ifdef Z
	dz = zmed(k)-Zplanet;
#endif
#endif
#ifdef CYLINDRICAL
	dx = ymed(j)*cos(xmed(i))-Xplanet;
	dy = ymed(j)*sin(xmed(i))-Yplanet;
#ifdef Z
	dz = zmed(k)-Zplanet;
#endif
#endif
#ifdef SPHERICAL
	dx = ymed(j)*cos(xmed(i))*sin(zmed(k))-Xplanet;
	dy = ymed(j)*sin(xmed(i))*sin(zmed(k))-Yplanet;
#ifdef Z
	dz = ymed(j)*cos(zmed(k))-Zplanet;
#endif
#endif
	dist2 = dx*dx+dy*dy+dz*dz;
	dist2 += rsm2;
	distance = sqrt(dist2);
	InvDist3 = 1.0/(dist2*distance);
	InvDist3 *= G*cellmass;
	
	fxi  = dx*InvDist3;
	fyi  = dy*InvDist3;
	tottorq += Xplanet*fyi-Yplanet*fxi;
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
    

  MPI_Reduce (tottorq, Gtottorq, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  return Gtottorq;
}
