#include "fargo3d.h"


void viscous_flux(real dt) {
    
    int i,j,k;
    i = j = k =0;

    int size_x = Nx + 2*NGHX; 
    int size_y = Ny + 2*NGHY; 
    int size_z = Nz + 2*NGHZ; 
    int pitch = Pitch_cpu;
    real viscositym, viscosityp,resc, resv,rescp,resvp;

    real *rho = Density->field_cpu;
    real *vx  = Vx->field_cpu;
    real *energy = Energy->field_cpu;
    real *drfnu = drFluxVisc->field_cpu;

// Calculate avg momentum
    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            resc = 0;
            rescp = 0;
            resv = 0;
            resvp = 0;
            for(i=0;i<size_x;i++) {
#ifdef ALPHAVISCOSITY
#ifdef ISOTHERMAL
                viscositym= ALPHA*.5*(energy[l]*energy[l]+energy[lym]*energy[lym])*sqrt(ymin(j)*ymin(j)*ymin(j)/(G*MSTAR));
                viscosityp= ALPHA*.5*(energy[lyp]*energy[lyp]+energy[l]*energy[l])*sqrt(ymin(j+1)*ymin(j+1)*ymin(j+1)/(G*MSTAR));
#else
                viscositym= ALPHA*GAMMA*(GAMMA-1.0)*(energy[l]+energy[lym])/(rho[l]+rho[lym])*sqrt(ymin(j)*ymin(j)*ymin(j)/(G*MSTAR));
                viscosityp= ALPHA*GAMMA*(GAMMA-1.0)*(energy[lyp]+energy[l])/(rho[lyp]+rho[l])*sqrt(ymin(j+1)*ymin(j+1)*ymin(j+1)/(G*MSTAR));
#endif
#else
            	viscositym =  NU;
            	viscosityp =  NU;
#endif
                resc += viscositym*.25*(rho[l]+rho[lxm]+rho[lym]+rho[lxm-pitch]) ;
	            resv +=(vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j); //centered on left, inner vertical edge in z
                rescp += viscosityp*.25*(rho[lyp]+rho[lxm+pitch]+rho[l]+rho[lxm]) ;
	            resvp +=(vx[lyp]-vx[l])/(ymed(j+1)-ymed(j))-.5*(vx[lyp]+vx[l])/ymin(j+1); //centered on left, inner vertical edge in z

            }
            resc /= (real)Nx;
            resv /= (real)Nx;
            rescp /= (real)Nx;
            resvp /= (real)Nx;
	        drfnu[l2D] -= (ymin(j+1)*ymin(j+1)*rescp*resvp -ymin(j)*ymin(j)*resc*resv)/((ymin(j+1)-ymin(j))*ymed(j))*dt;
        }
    }

    return;

}
