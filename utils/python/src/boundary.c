#include "evolve.h"
void set_bc(void) {
    
    ymin_bound_acc();
    ymax_bound_acc();
    return;

}

void ymax_bound(void) {
    int i,j,k;

  int jact;
  int jgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs;
  int lactbs_null;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,lgh,lghs,lactb,lactbs,lactbs_null,jgh,jact)
#endif
  for(k=0; k<size_z; k++) {
    for(j=0; j<NGHY; j++) {
      for(i=0; i<size_x; i++) {

	lgh = i + (ny+NGHY+j)*pitch + k*stride;
	lghs = i + (ny+NGHY+1+j)*pitch + k*stride;
	lactb = i + (ny+NGHY-1-j)*pitch + k*stride;
	lactbs = i + (ny+NGHY-1-j)*pitch + k*stride;
	lactbs_null = i + (ny+NGHY)*pitch + k*stride;
	jgh = (ny+NGHY+j);
	jact = (ny+NGHY-1-j);

	dens[lgh] = (dens[i+(ny+NGHY-1)*pitch]-dens[i+(ny+NGHY-2)*pitch])/(ymed(ny+NGHY-1)-ymed(ny+NGHY-2))*(ymed(jgh)-ymed(ny+NGHY-1))+dens[i+(ny+NGHY-1)*pitch];
	vx[lgh] = (vx[lactb]+ymed(jact)*omf)*sqrt(ymed(jact)/ymed(jgh)) - ymed(jgh)*omf;
	if (j<size_y-1)
		vy[lghs] = -vy[lactbs];
	vy[lactbs_null] = 0;
      }
    }
  }

    return;
}
void ymin_bound_acc(void) {

  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs_null;
  int nghy = NGHY;
  double sig1;
 double vr1;
 double ri1;
 double rm1;
  double omegaframe = omf;
  double nu_index = 0.5 + 2*params.flaringindex;
  double vr_index = -0.5 + 2*params.flaringindex;


  i = j = k = 0;

#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,sig1,ri1,rm1,vr1,lgh,lghs,lactb,lactbs_null,jgh,jact)
#endif
  for(k=0; k<size_z; k++) {
    for(j=0; j<NGHY; j++) {
      for(i=0; i<size_x; i++) {

	lgh = l;
	lghs = l;
	lactb = i + (2*nghy-j-1)*pitch + k*stride;
	lactbs_null = i + nghy*pitch + k*stride;
	jgh = j;
	jact = (2*nghy-j-1);

	sig1 = dens[ i + (nghy)*pitch + k*stride];
	ri1 = ymed(nghy);
	rm1 = ymin(nghy);
	vr1 = vy[ i + (nghy)*pitch + k*stride];
	dens[lgh] = sig1*pow(ri1/ymed(jgh),nu_index);
	vx[lgh] = (vx[lactb]+ymed(jact)*omegaframe)*sqrt(ymed(jact)/ymed(jgh))-ymed(jgh)*omegaframe;
	vy[lghs] = vr1*pow(ymed(jgh)/ri1,vr_index);
	vy[lactbs_null] = vr1*pow(rm1/ri1,vr_index);
      }
    }
  }
  return;
}

void ymax_bound_acc(void) {
  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs_null;
  double sig1;
 double ri1;
 double fh1;
  int nghy = NGHY;
  double mdot = params.mdot;
  double omegaframe = omf;
  double nu_0 = params.alpha*params.h*params.h;
  double nu_index = 0.5 + 2*params.flaringindex;
  double vnorm = -1.5*params.alpha*params.h*params.h;
  double vr_index = -0.5 + 2*params.flaringindex;
  double pi = M_PI;

  i = j = k = 0;

#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,sig1,ri1,fh1,lgh,lghs,lactb,lactbs_null,jgh,jact)
#endif
  for(k=0; k<size_z; k++) {
    for(j=0; j<NGHY; j++) {
      for(i=0; i<size_x; i++) {

	lgh = i + (ny+nghy+j)*pitch + k*stride;
	lghs = i + (ny+nghy+1+j)*pitch + k*stride;
	lactb = i + (ny+nghy-1-j)*pitch + k*stride;
	lactbs_null = i + (ny+nghy)*pitch + k*stride;
	jgh = (ny+nghy+j);
	jact = (ny+nghy-1-j);

	sig1 = dens[ i + (ny+nghy-1)*pitch + k*stride];
	ri1 = ymed(ny+nghy-1);
	fh1 = 3*pi*nu_0*pow(ri1,nu_index+0.5)*sig1;
	dens[lgh] = (fh1+mdot*(sqrt(ymed(jgh))-sqrt(ri1)))/(3*pi*nu_0*pow(ymed(jgh),nu_index+0.5));
	vx[lgh] = (vx[lactb]+ymed(jact)*omegaframe)*sqrt(ymed(jact)/ymed(jgh))-ymed(jgh)*omegaframe;
	if (j<size_y-1)
		vy[lghs] = vnorm*pow(ymin(jgh),vr_index)*mdot*sqrt(ymin(jgh))/(fh1+mdot*(sqrt(ymin(jgh))-sqrt(ri1)));
	vy[lactbs_null] = vnorm*pow(ymin(jgh),vr_index)*mdot*sqrt(ymin(jgh))/(fh1+mdot*(sqrt(ymin(jgh))-sqrt(ri1)));
      }
    }
  }
  return;
}
void ymin_bound(void) {

  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs;
  int lactbs_null;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k,lgh,lghs,lactb,lactbs,lactbs_null,jgh,jact)
#endif
for(k=0; k<size_z; k++) {
    for(j=0; j<NGHY; j++) {
      for(i=0; i<size_x; i++) {

	lgh = l;
	lghs = l;
	lactb = i + (2*NGHY-j-1)*pitch + k*stride;
	lactbs = i + (2*NGHY-j)*pitch + k*stride;
	lactbs_null = i + NGHY*pitch + k*stride;
	jgh = j;
	jact = (2*NGHY-j-1);

	dens[lgh] = (dens[i+(NGHY+1)*pitch]-dens[i+NGHY*pitch])/(ymed(NGHY+1)-ymed(NGHY))*(ymed(jgh)-ymed(NGHY))+dens[i+NGHY*pitch];
	vx[lgh] = (vx[lactb]+ymed(jact)*omf)*sqrt(ymed(jact)/ymed(jgh)) - ymed(jgh)*omf;
	vy[lghs] = -vy[lactbs];
	vy[lactbs_null] = 0;
//<\#>
      }
    }
  }
    return;
}
