#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NGHY 3

#define ymed(j) Ymed[(j)]
#define xmed(i) Xmed[(i)]
#define ymin(j) Ymin[(j)]
#define xmin(i) Xmed[(i)]
#define zone_size_x(j,k) (dx*ymed(j))
#define zone_size_y(j,k) (ymin(j+1)-ymin(j))
#define SurfY(j,k) ymin(j)*dx
#define SurfX(j,k) (ymin(j+1)-ymin(j))
#define InvVol(j,k) 2/(dx*(ymin(j+1)*ymin(j+1) - ymin(j)*ymin(j)))
#define Vol(j,k) 0.5*(dx*(ymin(j+1)*ymin(j+1) - ymin(j)*ymin(j)))

#define l   ((i)+(j)*(nx)+((k)*stride))
#define l_f(ii,jj,kk)   ((ii)+(jj)*(nx)+((kk)*stride))
#define lact   ((i)+(jact)*(nx)+((k)*stride))
#define lxp (((i)<(nx-1)) ? ((l)+1) : ((l)-(nx-1)))
#define lxm (((i)>0) ? ((l)-1) : ((l)+(nx-1)))

#define ixm ((i)>0 ? ((i)-1) : nx-1)
#define ixp ((i)<nx-1 ? ((i)+1) : 0)

#define lyp ((l)+nx)
#define lym ((l)-nx)

#define lzp (l+stride)
#define lzm (l-stride)

typedef struct Parameters {
    double alpha;
    double mp;
    double a;
    double omf;
    double h;
    double flaringindex;
    double nuindex;
    double vrindex;
    double mdot;
    double soft;

} Parameters;


int nx, ny, nz;
int size_x, size_y, size_z;
int stride, pitch;
double *Ymin, *Xmin, *Ymed, *Xmed;
double *rho, *vx, *vy;
double *momp, *momm;
double *rhos, *mdot,*mdotmed;
double *lstar,*lmed;
double *tauxx, *tauyy, *tauxy;
double *work, *work1;
double *mdotavg, *lavg,*lbar, *mdotbar,*rhoavg, *tauxyavg;
double *rhobar,*vxbar,*vybar;
double *divp, *divpavg;
double dx;

Parameters params;

double Nu(double x);
double Cs(double x);

void allocate_all(void);
void free_all(void);
void set_lstar(char *directory);
void set_indpot(char *directory);
void set_lmed(char *directory);
void set_rhostar(char *directory);
void set_mdot(char *directory);
void set_Ld(char *directory);
void set_Fd(char *directory);
void set_lam_ex(char *directory);
void set_lam_dep(char *directory);
void set_div(void);
void set_tensor(char *directory);
void set_slopes_x(double *q);
void transport_x(double *qs,double *res);
void set_slopes(double *q);
void transport(double *qs,double *res);
void clear_work(void);
void add_boundary(void);
void read_param_file(char *filename);
void read_domain(char *directory);
void read_files(int n, char *directory);


int main(int arc, char *argv[]) { 
  int n;
  int i,j,k;
  char directory[256];
  char param_fname[100];

  n = atoi(argv[1]);
  strcpy(directory,argv[2]);
  sprintf(param_fname,"%sparam_file.txt",directory);
//strcpy(param_fname,argv[2]);

  read_param_file(param_fname);

  

  size_x = nx;
  size_y = ny+2*NGHY;
  size_z = nz;
  stride = size_x*size_y;
  pitch = size_x;
  dx = 2*M_PI/nx;


  printf("Allocating.\n");
  allocate_all();


// Read in fluid variables
  printf("Reading files.\n");
  read_domain(directory); 
  read_files(n,directory);
  add_boundary();
  
// Set basic derived properties
  printf("Setting rhostar.\n");
  set_rhostar(directory);
  printf("Setting lstar.\n");
  set_lmed(directory);
  set_lstar(directory);
  printf("Setting mdot.\n");
  set_mdot(directory);
  set_indpot(directory);
  printf("Setting tensor.\n");
  set_tensor(directory);
  set_div();


// Calculate fluxes and torques and output

  printf("Setting lam ex\n");
  set_lam_ex(directory);
  printf("Setting L\n");
  set_Ld(directory);
  printf("Setting F\n");
  set_Fd(directory);
  printf("Setting lam dep\n");

  set_lam_dep(directory);


// Done
  printf("Freeing.\n");
  free_all();

  return 0;

}
void allocate_all(void) {
  int i,j,k;
  Xmin  =     (double*)malloc(sizeof(double)*(size_x+1));
  Ymin  =     (double*)malloc(sizeof(double)*(size_y+1));
  Xmed  =     (double*)malloc(sizeof(double)*(size_x));
  Ymed  =     (double*)malloc(sizeof(double)*(size_y));

  lavg  =     (double*)malloc(sizeof(double)*(size_y));
  lbar  =     (double*)malloc(sizeof(double)*(size_y));
  rhoavg  =     (double*)malloc(sizeof(double)*(size_y));
  rhobar  =     (double*)malloc(sizeof(double)*(size_y));
  vxbar  =     (double*)malloc(sizeof(double)*(size_y));
  vybar  =     (double*)malloc(sizeof(double)*(size_y));
  mdotavg  =     (double*)malloc(sizeof(double)*(size_y));
  mdotbar  =     (double*)malloc(sizeof(double)*(size_y));
  tauxyavg  =     (double*)malloc(sizeof(double)*(size_y));
  divpavg  =     (double*)malloc(sizeof(double)*(size_y));

  vx =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  vy =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  rho = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  rhos = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  lstar = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  lmed = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  mdot = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  mdotmed = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  momm = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  momp = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  
  tauxx =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  tauyy =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  tauxy = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  divp = (double*)malloc(sizeof(double)*size_x*size_y*size_z);


  work = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  work1 = (double*)malloc(sizeof(double)*size_x*size_y*size_z);

  for(k=0; k<size_z; k++) {
    for(j=0; j<size_y; j++) {
        divpavg[j] = 0;
        tauxyavg[j] = 0;
        mdotavg[j] = 0;
        mdotbar[j] = 0;
      for(i=0; i<size_x; i++) {
            vx[l] = 0;
            vy[l] = 0;
            rho[l] = 0;
            tauxx[l] = 0;
            tauyy[l] = 0;
            tauxy[l] = 0;
            work[l] = 0;
            work1[l] = 0;
            lstar[l] = 0;
            lmed[l] = 0;
            mdot[l] = 0;
            mdotmed[l] = 0;
            rhos[l] = 0;
            divp[l] = 0;
            momm[l] = 0;
            momp[l] = 0;
      }
    }
  }


    return;
}
void free_all(void) {


  free(Xmin);
  free(Ymin);
  free(Xmed);
  free(Ymed);
  free(lavg);
  free(lbar);
  free(rhoavg);
  free(rhobar);
  free(vxbar);
  free(vybar);
  free(mdot);
  free(mdotavg);
  free(mdotbar);
  free(tauxyavg);
  free(divp);
  free(divpavg);
  free(vx);
  free(momm);
  free(momp);
  free(vy);
  free(rho);
  free(rhos);
  free(tauxx);
  free(tauyy);
  free(tauxy);
  free(work);
  free(work1);
  free(lstar);
  free(lmed);


    return;
}
double Nu(double x) {
    return params.alpha*params.h*params.h*pow(x,params.nuindex);
}
double Cs(double x) {
    return params.h*pow(x,params.flaringindex-0.5);
}


void set_lstar(char *directory) {
    int i,j,k;
    k=0;
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/lstar.dat",directory);
    f = fopen(fname,"w");

/*
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            momm[l] = rho[l]*vx[l]*ymed(j);
            momp[l] = rho[l]*vx[lxp]*ymed(j);
        }
    }
    clear_work();
    transport(momm,NULL);
    clear_work();
    transport(momp,NULL);

    for(j=1;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            lstar[l] = (momm[l] + momp[lym])/(rhos[lym]+rhos[l]);
        }
    }
*/
    clear_work();
    transport(lmed,lstar);
    fwrite(&lstar[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fclose(f);
    double res;
    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += lstar[l];
        }
        lavg[j] = res/(double)nx;
    }
    return;
}


void set_lmed(char *directory) {
    int i,j,k;
    k=0;
/*
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            momm[l] = rho[l]*vx[l]*ymed(j);
            momp[l] = rho[l]*vx[lxp]*ymed(j);
        }
    }
    clear_work();
    transport_x(momm,NULL);
    clear_work();
    transport_x(momp,NULL);

    for(j=1;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            lmed[l] = (momm[l] + momp[lxm])/(rho[lxm]+rho[l]);
        }
    }

    transport_x(vx,lmed);
*/
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            lmed[l] = .5*(vx[l] + vx[lxp])*ymed(j); 
        }
    }


    double res;
    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += lmed[l];
        }
        lbar[j] = res/(double)nx;
    }

    return;

}
void set_rhostar(char *directory) {
    int i,j,k;
    k=0;
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/rhostar.dat",directory);
    f = fopen(fname,"w");
    clear_work();
    transport(rho,rhos);

    fwrite(&rhos[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fclose(f);
    double res;
    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += rhos[l];
        }
        rhoavg[j] = res/(double)nx;
    }
    return;

}
void set_mdot(char *directory) {
    int i,j,k;
    k = 0;

    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/mdot.dat",directory);
    f = fopen(fname,"w");
    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {        
            mdot[l] = rhos[l] * - 2 * M_PI*ymin(j)*vy[l];
            mdotmed[l] = rho[l]*-2*M_PI*ymed(j)*.5*(vy[l]+vy[lyp]);
        }
    }
    fwrite(&mdot[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fclose(f);
    double res,res1;
    for(j=0;j<size_y;j++) {
        res = 0;
        res1=0;
        for(i=0;i<size_x;i++) {
            res += mdot[l];
            res1 += mdotmed[l];
        }
        mdotavg[j] = res/(double)nx;
        mdotbar[j] = res1/(double)nx;
    }


    return;
}

void set_Ld(char *directory) {
    int i,j,k;
    k=0;
    FILE *ft,*fw,*fd;
    double Ld,Lw,Ltot;
    double rest, resw;
    double fac;
    char fnamet[256],fnamew[256],fnamed[256];
    sprintf(fnamet,"%stemp_files/Lt.dat",directory);
    sprintf(fnamew,"%stemp_files/Lw.dat",directory);
    sprintf(fnamed,"%stemp_files/Ld.dat",directory);
    
    ft = fopen(fnamet,"w");
    fw = fopen(fnamew,"w");
    fd = fopen(fnamed,"w");

    
    for(j=NGHY;j<size_y-NGHY;j++) {
        Ld = 2*M_PI*ymed(j)*rhobar[j]*lbar[j];
        rest = 0;
        for(i=0;i<size_x;i++) {

            rest += 2*M_PI*ymed(j)*rho[l]*lmed[l];

        }
        rest /= (double)nx;
        resw = rest - Ld;
        fwrite(&rest,sizeof(double),1,ft);
        fwrite(&resw,sizeof(double),1,fw);
        fwrite(&Ld,sizeof(double),1,fd);
    }


    fclose(ft);
    fclose(fw);
    fclose(fd);
    return;
}

void set_Fd(char *directory) {
    int i,j,k;
    k=0;
    FILE *ft,*fw,*fd;
    double Fd,Fw,Ftot;
    double resw, rest;
    char fnamet[256],fnamew[256],fnamed[256];
    sprintf(fnamet,"%stemp_files/Ft.dat",directory);
    sprintf(fnamew,"%stemp_files/Fw.dat",directory);
    sprintf(fnamed,"%stemp_files/Fd.dat",directory);
    
    ft = fopen(fnamet,"w");
    fw = fopen(fnamew,"w");
    fd = fopen(fnamed,"w");

    
    for(j=NGHY;j<size_y-NGHY+1;j++) {
        Fd = -mdotavg[j] * lavg[j] - 2*M_PI*ymin(j)*ymin(j)*tauxyavg[j];
        rest = 0;
        for(i=0;i<size_x;i++) {
            rest += -mdot[l]*lstar[l] - M_PI*ymin(j)*ymin(j)*(tauxy[l]+tauxy[lxp]);

        }
        rest /= (double)nx;
        resw = rest - Fd;
        fwrite(&rest,sizeof(double),1,ft);
        fwrite(&resw,sizeof(double),1,fw);
        fwrite(&Fd,sizeof(double),1,fd);
    }


    fclose(ft);
    fclose(fw);
    fclose(fd);


    return;
}

void set_indpot(char *directory) {
    int i,j,k;
    k=0;

    double res,resj;

    double resx, resy;
    double dm;
    resx=0;
    resy=0;
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
           
            dm= Vol(j,k)*rho[l];
            resx += dm*cos(xmed(i))/(ymed(j)*ymed(j));
            resy += dm*sin(xmed(i))/(ymed(j)*ymed(j));
            
        }

    }

    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/lamind.dat",directory);
    f = fopen(fname,"w");
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            
            resj = params.mp*ymed(j)*sin(xmed(i))/(params.a*params.a); 
            resj += ymed(j)*(cos(xmed(i))*resy - sin(xmed(i))*resx); 
            resj *= -2*M_PI*ymed(j)*rho[l];
            res += resj;
            
        }
        res /= (double)nx;
        fwrite(&res,sizeof(double),1,f);

    }
    fclose(f);
    return;
}

void set_lam_ex(char *directory) {
    int i,j,k;
    k = 0;
    FILE *fl;
    char fname[256];
    sprintf(fname,"%stemp_files/lamex.dat",directory);

    double smoothing = params.soft*params.h*pow(params.a,1 + params.flaringindex);
    smoothing *= smoothing;
    double rad,res;
    double pot;

    fl = fopen(fname,"w");
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {    
            rad = ymed(j)*ymed(j) + params.a*params.a -2*params.a*ymed(j)*cos(xmed(i))+smoothing;
            rad *= sqrt(rad);
            pot = -params.mp * ymed(j)*params.a*sin(xmed(i))/rad;
            res += -2*M_PI*ymed(j)*rho[l]*pot;
        }
        res /= (double)nx;
        fwrite(&res,sizeof(double),1,fl);
    }

    fclose(fl);
    return;
}

void set_lam_dep(char *directory) {
    int i,j,k;
    k=0;
    FILE *f;
    double vyc;
    char fname[256];
    sprintf(fname,"%stemp_files/lamdep.dat",directory);
    f = fopen(fname,"w");
    if (f==NULL) printf("File not opened.\n");
    double fac,res,sfac;
    for(j=NGHY;j<size_y-NGHY;j++) {
        fac = 2*M_PI*ymed(j);
        res=0;
        for(i=0;i<size_x;i++) {
            //sfac = 1 + (rho[l] - rhobar[j])/rhobar[j];
            //res += .5*(divp[l] + divp[lxp])/sfac - divpavg[j];
            res -= (lavg[j+1]*mdot[lyp]-lavg[j]*mdot[l] - lbar[j]*(mdot[lyp]-mdot[l]))/(ymin(j+1)-ymin(j));
            //res -= mdotmed[l]*(lavg[j+1]-lavg[j])/(ymin(j+1)-ymin(j));
            res -= fac*rhobar[j]*(vy[lyp]*lstar[lyp]-vy[l]*lstar[l]  - lmed[l]*(vy[lyp]-vy[l]) )/(ymin(j+1)-ymin(j));
        }
        res /= (double)nx;
        res += divpavg[j];
        fwrite(&res,sizeof(double),1,f);
    }

    fclose(f);
    return;
}
void set_div(void) {
    int i,j,k;
    k=0;
    double res,res1,res2;
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
       for(i=0;i<size_x;i++) {

	        res += 2.0*(tauxx[l]-tauxx[lxm])/(zone_size_x(j,k)*(rho[l]+rho[lxm]));
	        res += 2.0*(ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j))*ymed(j)*ymed(j)*(rho[lxm]+rho[l]));
        }
       res /= (double)nx;
       res *= rhobar[l]*2*M_PI*ymed(j)*ymed(j);
       res -= 2*M_PI*(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]-ymin(j)*ymin(j)*tauxyavg[j])/(ymin(j+1)-ymin(j));
       divpavg[j] = res;
    }


    return;
}
/*
void set_div(void) {
    int i,j,k;
    k=0;

//    clear_work();
    double res,sigfac;
    for(j=NGHY;j<size_y-NGHY;j++) {
        divpavg[j] = 2*M_PI*(tauxyavg[j+1]*ymin(j+1)*ymin(j+1) - tauxyavg[j]*ymin(j)*ymin(j))/(ymin(j+1)-ymin(j));
       for(i=0;i<size_x;i++) {
            //sigfac = .5*(rho[l]+rho[lxm]);
            res += sigfac*2*M_PI*ymed(j)*(tauxx[l]-tauxx[lxm])/(dx);
            res += sigfac*2*M_PI*(tauxy[lyp]*ymin(j+1)*ymin(j+1) - tauxy[l]*ymin(j)*ymin(j))/(ymin(j+1)-ymin(j));
       }
        res /= (double)nx;
        divpavg[j] = res - divpavg[j];
    }

    return;

}
*/
void set_tensor(char *directory) {
    int i,j,k;
    double visc, viscm, div_v;
    double cs, csm;
    double rhoc, rhocb;
    k = 0;


    for(j=1;j<size_y-1;j++) {
        cs = Cs(ymed(j));
        csm = Cs(ymin(j));
        cs *= cs; csm *= csm;
        visc = Nu(ymed(j));
        viscm = Nu(ymin(j));
        tauxyavg[j] = viscm*.5*(rhobar[j]+rhobar[j-1])*( 
                (vxbar[j]-vxbar[j-1])/(ymed(j)-ymed(j-1))
                -.5*(vxbar[j]+vxbar[j-1])/ymin(j));

        for(i=0;i<size_x;i++) {
            rhoc=.25*(rho[l]+rho[lxm]+rho[lym]+rho[lxm-pitch]);
            div_v = (vx[lxp]-vx[l])*SurfX(j,k);
            div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
            div_v *= 2.0/3.0*InvVol(j,k);

            tauxx[l] = visc*rho[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
            tauxx[l] += visc*rho[l]*(vy[lyp]+vy[l])/ymed(j);
            tauxx[l] += -cs*rho[l]; // Isothermal Pressure
            tauyy[l] = visc*rho[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
            tauyy[l] += -cs*rho[l]; // Isothermal Pressure
            tauxy[l] = viscm*rhoc*((vy[l]-vy[lxm])/(dx*ymin(j)) + ymin(j)*(vx[l]/ymed(j)-vx[lym]/ymed(j-1))/(ymed(j)-ymed(j-1)));//-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z
            

        }
    }
    return;
}

void set_slopes_x(double *q) {
    int i,j,k;
    k=0;
    double dqm, dqp;

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {

        dqm = (q[l]-q[lxm]);
        dqp = (q[lxp]-q[l]);
        if(dqp*dqm<=0) work[l] = 0;
        else  work[l] = 2.*dqp*dqm/(zone_size_x(j,k)*dqm+dqp);
      }
    }
    j=0;
    for(i=0;i<size_x;i++) {
        work[l] = 0;
    }
    j=size_y-1;
    for(i=0;i<size_x;i++) {
        work[l] = 0;
    }

    
    return;
}
void transport_x(double *qs,double *res) {
    int i,j,k;
    k = 0;
   
    if (res == NULL) {
        res = qs;
    }

    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            work1[l] = qs[l];
        }
    }
    set_slopes_x(qs);   // Slope is in work array

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {

        res[l] = work1[lxm] + 0.5 * zone_size_x(j,k)*work[lxm];

      }

    }
    return;
}

void set_slopes(double *q) {
    int i,j,k;
    k=0;
    double dqm, dqp;

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {

        dqm = (q[l]-q[lym])/zone_size_y(j,k);
        dqp = (q[lyp]-q[l])/zone_size_y(j+1,k);
        if(dqp*dqm<=0) work[l] = 0;
        else  work[l] = 2.*dqp*dqm/(dqm+dqp);
      }
    }
    j=0;
    for(i=0;i<size_x;i++) {
        work[l] = 0;
    }
    j=size_y-1;
    for(i=0;i<size_x;i++) {
        work[l] = 0;
    }

    
    return;
}
void transport(double *qs,double *res) {
    int i,j,k;
    k = 0;
   
    if (res == NULL) {
        res = qs;
    }

    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            work1[l] = qs[l];
        }
    }
    set_slopes(qs);   // Slope is in work array

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {

    	if (vy[l]>0.) {
            res[l] = work1[lym] + 0.5 * (zone_size_y(j-1,k))*work[lym];
        }
	    else {
	        res[l] = work1[l] - 0.5 * (zone_size_y(j,k))*work[l];
        }

      }

    }
    return;
}
void clear_work(void){
    int i,j,k;
    k=0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                work[l] = 0;
                work1[l] =0;
            }
        }
    }
    return;
}


void add_boundary(void) {
    double fh1;
    double fac,rhom;
    double facm;
    int i,j,k,jact;
    k = 0;
    
    jact = NGHY;
    for(j=0;j<NGHY;j++) {
            for(i=0;i<size_x;i++) {
                rho[l] = rho[lact]*pow(ymed(j)/ymed(jact),-params.nuindex);
                vx[l] = vx[lact]*pow(ymed(j)/ymed(jact),-.5);
                vy[l] = vy[lact]*pow(ymin(j)/ymin(jact),params.vrindex);

            }
    }
    jact=size_y-NGHY-1;
    for(j=size_y-NGHY;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                fh1 = rho[lact]*3*M_PI*Nu(ymed(jact))*sqrt(ymed(jact));
                fac = 3*M_PI*Nu(ymed(j))*sqrt(ymed(j));
                facm = 3*M_PI*Nu(ymin(j))*sqrt(ymin(j));
                
                rho[l] = (fh1 + params.mdot*(sqrt(ymed(j))-sqrt(ymed(jact)))/(3*M_PI))/fac;
                rhom = (fh1 + params.mdot*(sqrt(ymin(j))-sqrt(ymed(jact)))/(3*M_PI))/facm;
                vy[l] = params.mdot/(-2*M_PI*ymin(j)*rhom);
                vx[l] = vx[lact]*pow(ymed(j)/ymed(jact),-.5);
                //vy[l] = vy[lact]*pow(ymin(j)/ymin(jact),params.vrindex);

            }
    }
    double res1,res2,res3;

    for(j=0;j<size_y;j++) {
        res1 = 0;
        res2 = 0;
        res3 = 0;
        for(i=0;i<size_x;i++) {
            res1 += rho[l];
            res2 += vx[l];
            res3 += vy[l];
        }
        rhobar[j] = res1/(double)nx;
        vxbar[j] = res2/(double)nx;
        vybar[j] = res3/(double)nx;
    }
    return;

}

void read_param_file(char *filename) {
    FILE *f;
    
    printf("Reading %s\n",filename);
    f = fopen(filename,"r");
    if (f == NULL) {
        printf("Can't find parameter file, %s\n",filename);
        exit(0);
    }
    fscanf(f,"%d\n",&nx);
    fscanf(f,"%d\n",&ny);
    fscanf(f,"%d\n",&nz);
    fscanf(f,"%lg\n",&params.alpha);
    fscanf(f,"%lg\n",&params.mp);
    fscanf(f,"%lg\n",&params.a);
    fscanf(f,"%lg\n",&params.omf);
    fscanf(f,"%lg\n",&params.h);
    fscanf(f,"%lg\n",&params.flaringindex);
    fscanf(f,"%lg\n",&params.mdot);
    fscanf(f,"%lg\n",&params.soft);

    params.nuindex = 2*params.flaringindex + 0.5;
    params.vrindex = params.nuindex - 1.0;

    fclose(f);
    printf("nx=%d\tny=%d\tnz=%d\n",nx,ny,nz);
    printf("alpha=%.1e\tmp=%.1e\ta=%lg\n",params.alpha,params.mp,params.a);
    printf("omf=%lg\th=%lg\tflaring=%lg\n",params.omf,params.h,params.flaringindex);
    printf("mdot=%.2e\tsoft=%lg\n",params.mdot,params.soft);
    return;
}

void read_domain(char *directory) {
    FILE *fx, *fy;
    char filename[512];
    char filename2[512];
    double temp;
    int i,j;

    sprintf(filename,"%sdomain_x.dat",directory);
    printf("Reading %s\n",filename);
    fx = fopen(filename,"r");
    if (fx == NULL) printf("Error reading %s\n",filename);

    for(i=0;i<size_x+1;i++) {
        fscanf(fx,"%lg\n",&Xmin[i]);
    }
    fclose(fx);
    sprintf(filename2,"%sdomain_y.dat",directory);
    printf("Reading %s\n",filename2);
    fy = fopen(filename2,"r");
    if (fy == NULL) printf("Error reading %s\n",filename);
    for(j=0;j<size_y+1;j++) {
            fscanf(fy,"%lg\n",&Ymin[j]);
    }
    fclose(fy);
   

    for(i=0;i<size_x;i++) {
        Xmed[i] = .5*(Xmin[i] + Xmin[i+1]);
    }
    for(j=0;j<size_y;j++) {
        Ymed[j] = .5*(Ymin[j] + Ymin[j+1]);
    }
    return;
}

void read_files(int n, char *directory) {
    char filename[512];
    FILE *fd,*fx,*fy;
    int i,j,k;
    
    sprintf(filename,"%sgasdens%d.dat",directory,n);
    printf("Reading %s\n",filename);
    fd = fopen(filename,"r");
    if (fd == NULL) printf("Error loading %s\n",filename); 
    sprintf(filename,"%sgasvx%d.dat",directory,n);
    printf("Reading %s\n",filename);
    fx = fopen(filename,"r");
    if (fx == NULL) printf("Error loading %s\n",filename); 
    sprintf(filename,"%sgasvy%d.dat",directory,n);
    printf("Reading %s\n",filename);
    fy = fopen(filename,"r");
    if (fy == NULL) printf("Error loading %s\n",filename); 
    for(k=0;k<size_z;k++) {
        for(j =NGHY; j<size_y-NGHY;j++) {
            for(i=0;i<size_x;i++) {
                fread(&rho[l],sizeof(double),1,fd);
                fread(&vx[l],sizeof(double),1,fx);
                fread(&vy[l],sizeof(double),1,fy);
                vx[l] += params.omf*ymed(j);
            }
        }
    }

    fclose(fd);
    fclose(fx);
    fclose(fy);
    return;

}




