//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void compute_potential(real dt) {

  real OmegaNew, domega;
  int i;
  int subcycling = 5;
  
  if (Corotating) GetPsysInfo (MARK);
  
#ifdef GPU
  //Copy all the planetary data to device
  DevMemcpyH2D(Sys->x_gpu, Sys->x_cpu, sizeof(real)*(Sys->nb+1));
  DevMemcpyH2D(Sys->y_gpu, Sys->y_cpu, sizeof(real)*(Sys->nb+1));
  DevMemcpyH2D(Sys->z_gpu, Sys->z_cpu, sizeof(real)*(Sys->nb+1));
  DevMemcpyH2D(Sys->mass_gpu, Sys->mass_cpu, sizeof(real)*(Sys->nb+1));
#endif
  
  DiskOnPrimaryAcceleration = ComputeAccel(0.0, 0.0, 0.0, 0.0, 0.0);
  FARGO_SAFE(ComputeIndirectTerm());
  FARGO_SAFE(Potential()); // Gravitational potential from star and planet(s)
#ifdef FTPOTENTIAL
  if (InitPotential) {
        FARGO_SAFE(init_potential_m());
        InitPotential = NO;
  }
#endif

#ifndef FIXEDPSYS
  FARGO_SAFE(AdvanceSystemFromDisk(dt));

  if (ThereIsACentralBinary)
    subcycling = 30;		/* Arbitrary number of subcycles which
				   should fit most needs */
  for (i = 0; i < subcycling; i++)
    FARGO_SAFE(AdvanceSystemRK5(1.0/((double)(subcycling))*dt));
  
  if (Corotating) {
    OmegaNew = GetPsysInfo(GET)/dt;
    domega = OmegaNew-OMEGAFRAME;
    FARGO_SAFE(CorrectVtheta(domega));
    OMEGAFRAME = OmegaNew;
  }
  RotatePsys(OMEGAFRAME*dt);
#endif
}

void Potential_cpu() {
  
//<USER_DEFINED>
#ifdef FTPOTENTIAL
  INPUT(FTPot)
#endif
  OUTPUT(Pot);
  real planetmass_taper;
  if (MASSTAPER == 0.0)
    planetmass_taper = 1.0;
  else
    planetmass_taper = (PhysicalTime >= MASSTAPER ? 1.0 : .5*(1.0-cos(M_PI*PhysicalTime/MASSTAPER)));
//<\USER_DEFINED>

//<EXTERNAL>
  real* pot  = Pot->field_cpu;
#ifdef FTPOTENTIAL
  real* ftpot = FTPot->field_cpu;
#endif
  real* xplanet = Sys->x_cpu;
  real* yplanet = Sys->y_cpu;
  real* zplanet = Sys->z_cpu;
  real* mplanet = Sys->mass_cpu;
  int nb        = Sys->nb;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int indirect_term = INDIRECTTERM;
  int indirectx = IndirectTerm.x;
  int indirecty = IndirectTerm.y;
  int indirectz = IndirectTerm.z;
  real taper = planetmass_taper;
  int istar1 = BinaryStar1;
  int istar2 = BinaryStar2;
  int binary_true = ThereIsACentralBinary;
//<\EXTERNAL>
  
//<INTERNAL>
  int i;
  int j;
  int k;
  int n;
  real smoothing;
  real dist;
  real rroche;
  real planetdistance;
  real mp;
  real invd3;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real ASPECTRATIO(1);
// real ROCHESMOOTHING(1);
// real FLARINGINDEX(1);
// real THICKNESSSMOOTHING(1);
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
#ifndef NODEFAULTSTAR
#ifdef SPHERICAL
	pot[l] =  -G*MSTAR/ymed(j); //Potential from star
#endif
#ifdef CYLINDRICAL
	pot[l] =  -G*MSTAR/sqrt(ymed(j)*ymed(j)+ZC*ZC); //Potential from star
#endif
#ifdef CARTESIAN
	pot[l] = -G*MSTAR/sqrt(XC*XC+YC*YC+ZC*ZC);
#endif
#else
	pot[l] = 0.0; // No default star
#endif

	for(n=0; n<nb; n++) {
	  mp = mplanet[n]*taper;

	  planetdistance = sqrt(xplanet[n]*xplanet[n]+
				yplanet[n]*yplanet[n]+
				zplanet[n]*zplanet[n]);
	  rroche = planetdistance*pow((1.0/3.0*mp/MSTAR),1.0/3.0);

	  if (ROCHESMOOTHING != 0)
	    smoothing = rroche*ROCHESMOOTHING;
	  else
	    smoothing = ASPECTRATIO*
	      pow(planetdistance/R0,FLARINGINDEX)*
	      planetdistance*THICKNESSSMOOTHING;

	  smoothing*=smoothing;

	  dist = ((XC-xplanet[n])*(XC-xplanet[n])+
		  (YC-yplanet[n])*(YC-yplanet[n])+
		  (ZC-zplanet[n])*(ZC-zplanet[n]));
#ifndef NODEFAULTSTAR
	  if (indirect_term == YES) {
	    /* Indirect term due to planets */
	    pot[l] += G*mp*(XC*xplanet[n]+YC*yplanet[n]+ZC*zplanet[n])/(planetdistance*
										planetdistance*
										planetdistance); 
#ifndef NOINDIRECTDISK
	    pot[l] -= indirectx*XC + indirecty*YC + indirectz*ZC; /* Indirect term due to gas */
#endif
	  }
#endif
#ifdef NODEFAULTSTAR
	  if (binary_true && (indirect_term == YES)) {
	    if ((n != istar1) && (n != istar2)) { /* For all non-stellar objects */
	      planetdistance = sqrt((xplanet[n]-xplanet[istar1])*(xplanet[n]-xplanet[istar1])+
				    (yplanet[n]-yplanet[istar1])*(yplanet[n]-yplanet[istar1])+
				    (zplanet[n]-zplanet[istar1])*(zplanet[n]-zplanet[istar1]));
	      invd3 = 1.0/(planetdistance*planetdistance*planetdistance);
	      pot[l] += G*mp*invd3*mplanet[1]*((xplanet[n]-xplanet[istar1])*XC+	\
					       (yplanet[n]-yplanet[istar1])*YC+	\
					       (zplanet[n]-zplanet[istar1])*ZC)/\
		(mplanet[1]+mplanet[2]);
	      
	      planetdistance = sqrt((xplanet[n]-xplanet[istar2])*(xplanet[n]-xplanet[istar2])+
				    (yplanet[n]-yplanet[istar2])*(yplanet[n]-yplanet[istar2])+
				    (zplanet[n]-zplanet[istar2])*(zplanet[n]-zplanet[istar2]));
	      invd3 = 1.0/(planetdistance*planetdistance*planetdistance);
	      pot[l] += G*mp*invd3*mplanet[2]*((xplanet[n]-xplanet[istar2])*XC+	\
					       (yplanet[n]-yplanet[istar2])*YC+	\
					       (zplanet[n]-zplanet[istar2])*ZC)/\
		(mplanet[1]+mplanet[2]);
	      
	    }
	  }
#endif
#ifndef FTPOTENTIAL
	  pot[l] += -G*mp/sqrt(dist+smoothing); //Potential from planets
#else
      pot[l] += ftpot[l]*taper;
#endif
//      }

	}
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

void init_potential_m(void) {
    OUTPUT(Pot);
    OUTPUT(FTPot);
    OUTPUT(PotP);
    int i;
    int j;
    int k;
    int n;
    int m;
    int size_x;
    int size_y;
    int size_z;
    real facr;
    real faci;
    int mstart;
    int mend;
    
    int nb        = Sys->nb;
    real* xplanet = Sys->x_cpu;
    real* yplanet = Sys->y_cpu;
    real* zplanet = Sys->z_cpu;
    real* mplanet = Sys->mass_cpu;
    real smoothing;
    real dist;
    real rroche;
    real planetdistance;
    real mp;
    real* pot = Pot->field_cpu;
    real* fpot = FTPot->field_cpu;
    real* potp = PotP->field_cpu;
#if defined(MSTART) && defined(MEND)
    mstart = MSTART;
    mend = MEND;
#else
    mstart = 0;
    mend = 0;
#endif
    size_x = Nx + 2*NGHX;
    size_y = Ny + 2*NGHY;
    size_z = Nz + 2*NGHZ;
    i=j=k=m=0;
#ifdef Z
        for (k=0; k<size_z; k++) {
#endif
#ifdef Y
            for (j=0; j<size_y; j++) {
#endif
                for(i=0; i<size_x; i++) {

                for(n=0; n<nb; n++) {
                  mp = mplanet[n];

                  planetdistance = sqrt(xplanet[n]*xplanet[n]+
                            yplanet[n]*yplanet[n]+
                            zplanet[n]*zplanet[n]);
                  rroche = planetdistance*pow((1.0/3.0*mp/MSTAR),1.0/3.0);

                  if (ROCHESMOOTHING != 0)
                    smoothing = rroche*ROCHESMOOTHING;
                  else
                    smoothing = ASPECTRATIO*
                      pow(planetdistance/R0,FLARINGINDEX)*
                      planetdistance*THICKNESSSMOOTHING;

                  smoothing*=smoothing;

                  dist = ((XC-xplanet[n])*(XC-xplanet[n])+
                      (YC-yplanet[n])*(YC-yplanet[n])+
                      (ZC-zplanet[n])*(ZC-zplanet[n]));
                  potp[l] += -G*mp/sqrt(smoothing + dist);
                }
                }
#ifdef Y 
            }
#endif
#ifdef Z
        }
#endif
    WriteField(PotP,0);
    printf("%d\t%d\n",mstart,mend);
    for (m=mstart; m<mend+1; m++) {                
        printf("Adding the m=%d mode\n",m);
#ifdef Z
        for (k=0; k<size_z; k++) {
#endif
#ifdef Y
            for (j=0; j<size_y; j++) {
#endif
                i=0;
                facr = .5*cos(m*xmed(i))*potp[l];
                faci = -.5*sin(m*xmed(i))*potp[l];
                i=Nx + 2*NGHX-1;
                facr += .5*cos(m*xmed(i))*potp[l];
                faci += -.5*sin(m*xmed(i))*potp[l];
                for( i=1; i<size_x-1; i++) {
                    facr += cos(m*xmed(i))*potp[l];
                    faci += -sin(m*xmed(i))*potp[l];
                }
                facr /= (real)size_x;
                faci /= (real)size_x;

                for(i=0;i<size_x;i++) {
                    fpot[l] += 2*(facr*cos(m*xmed(i)) - faci*sin(m*xmed(i)));
        
                }
#ifdef Y 
            }
#endif
#ifdef Z
        }
#endif
    }


#ifdef Z
        for (k=0; k<size_z; k++) {
#endif
#ifdef Y
            for (j=0; j<size_y; j++) {
#endif
                for(i=0; i<size_x; i++) {
                    pot[l] += fpot[l];
                }
#ifdef Y 
            }
#endif
#ifdef Z
        }
#endif

        FARGO_SAFE(WriteField(Pot,0));
        FARGO_SAFE(WriteField(FTPot,0));
}




