#include "fargo3d.h"

void InitDensPlanet() {
   int i,j,k;
  real *field;
  real viscosity,fac1,fac2,fac3,nuind,nu0,r;
  field = Density->field_cpu;
  nuind = 0.5+2*FLARINGINDEX;
  nu0 = ASPECTRATIO*ASPECTRATIO*ALPHA;
  boolean GhostInclude = TRUE;

  

    int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;
   for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
        r = Ymed(j);
      //  fac1 = -1.5 * nu0 * pow(r/YMAX,-0.5+2*FLARINGINDEX);
        fac2 = MDOT/(3*M_PI*nu0*pow(r,nuind));
//        viscosity = nu0*pow(r,nuind);
//        fac1=SIGMAOUT*YMAX/r-SIGMAIN*YMIN/r;
//        fac3=sqrt(YMAX) - sqrt(YMIN);
//        fac2 = sqrt(r)-sqrt(YMIN);
        //fac1=pow(YMAX,.5-nuind)*YMIN*sig_in-pow(YMIN,.5-nuind)*YMAX*SIGMAOUT;
        //fac2=-pow(YMAX,1-nuind)*YMIN*SIGMAIN+pow(YMIN,1-nuind)*YMAX*SIGMAOUT;
        //fac3 = pow(YMIN,1-nuind)*pow(YMAX,0.5-nuind)-pow(YMAX,1-nuind)*pow(YMIN,0.5-nuind);
        //fac1 *= pow(r,-nuind);
        //fac2 *= pow(r,-.5*-nuind);
      for (i = begin_i; i<end_i; i++) {
          //
          //
//         field[l] = SIGFLOOR + SIGMA0 * exp( - (Ymed(j)-1.8)*(Ymed(j)-1.8)/.01);
//        field[l]= SIGMAIN*YMIN/r + fac1*fac2/fac3;//(fac1-fac2)/fac3;
          field[l] = fac2;

      }
    }
  }
}

void InitSoundSpeedPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t, omega, vk;
  real rho;
  FILE *fo;
  real *d;
  real *e;

  field = Energy->field_cpu;
  d = Density->field_cpu;

  boolean GhostInclude = TRUE;

  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
	r = Ymed(j);
	vk = sqrt(G*MSTAR/r);
	field[l] = ASPECTRATIO * pow(Ymed(j)/R0, FLARINGINDEX) * vk; //sqrt(G*MSTAR/Ymed(j))
#ifdef ADIABATIC
	field[l] = field[l]*field[l]*d[l]/(GAMMA-1.0);
#endif
      }
    }
  }
}

void InitVazimPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t;
  real rho;
  FILE *fo;
  real vt, omega;
  real *vr;
  real *cs;
  real md0,nuind,nu0;
  real vr0, fac1, fac2;
  field = Vx->field_cpu;
  vr = Vy->field_cpu;
  cs = Energy->field_cpu;

  boolean GhostInclude = TRUE;
  nuind = 0.5+2*FLARINGINDEX;
  nu0 = ALPHA*ASPECTRATIO*ASPECTRATIO;

  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
       r = Ymed(j);

   //     vr0=-1.5*nu0*pow(r,nuind);
   //     fac1 = pow(YMAX,.5-nuind)*SIGMAIN*YMIN-pow(YMIN,.5-nuind)*SIGMAOUT*YMAX;
   //
   //     fac2 = fac1 + (-pow(YMAX,1-nuind)*SIGMAIN*YMIN+pow(YMIN,1-nuind)*SIGMAOUT*YMAX)/sqrt(r);


//        fac1=SIGMAOUT*YMAX/r-SIGMAIN*YMIN/r;
//        fac1 /= (sqrt(YMAX) - sqrt(YMIN));
//        fac2 = SIGMAIN*YMIN/r + fac1*(sqrt(r)-sqrt(YMIN));
//        md0 = 1.5*nu0*fac1;

        fac1 = -1.5 * nu0 * pow(r,nuind-1);

      for (i = begin_i; i<end_i; i++) {

	omega = sqrt(G*MSTAR/r/r/r);
	vt = omega*r*sqrt(1.0+pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
			  (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));

	vt -= OMEGAFRAME*r;

	field[l] = vt*(1.+ASPECTRATIO*NOISE*(drand48()-.5));
	//vr[l] = r*omega*ASPECTRATIO*NOISE*(drand48()-.5);
//    vr[l] = -1.5*ALPHA*ASPECTRATIO*ASPECTRATIO*pow(r,-0.5+2*FLARINGINDEX);
    vr[l] = -1.5 * ALPHA*ASPECTRATIO*ASPECTRATIO*pow(Ymin(j),2*FLARINGINDEX-0.5);
      }
    }
  }
}

void CondInit() {

  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);

  int i,j,k;
  int index;
  real vt;
  real *field;
  real *rho;
  real *v1;
  real *v2;
  real *e;
#ifdef PLANETS
  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  if(COROTATING)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);
  else
#endif
    OMEGAFRAME = OMEGAFRAME;

  InitDensPlanet ();
  InitSoundSpeedPlanet ();
  InitVazimPlanet ();
}
