#include "evolve.h"

int main(int argc, char *argv[]) {
    int j;
    int n;
    int status;
    char directory[256],outputdir[256];



    if (argc > 4) {
        printf("Too many arguments.Only using first 3.\n");
    }
    n = atoi(argv[1]);
    time_step = atof(argv[2]);
    strcpy(directory,argv[3]);

    nb = 2;

    sprintf(outputdir,"%stemp_files/",directory);
    status = mkdir(outputdir,S_IRWXU | S_IRWXG | S_IRWXO);
    if (status == -33)  {
        printf("Status is -33 for some reason.\n");
    }
    //sprintf(param_fname,"%sparam_file.txt",directory);


    
    read_param_file(directory);
    nx = params.nx;
    ny = params.ny;
    nz = params.nz;
    CFL = params.cfl;
    IndirectTerm = params.indirect;

    size_x = nx;
    size_y = ny+2*NGHY;
    size_z = nz;
    stride = size_x*size_y;
    pitch = size_x;
    pitch2d = 0;
    //size_t double_to_int = sizeof(double)/sizeof(int);
    pitch2d_int = 0 ;//pitch*(int)double_to_int;
    dx = 2*M_PI/nx;

    if (size_x % 2 == 0) NEVEN = TRUE;
    else NEVEN = FALSE;

    num_threads = 1;

#ifdef _OPENMP
    //omp_set_num_threads(8);
    num_threads = omp_get_max_threads();
    printf("Using %d threads\n",num_threads);
#endif
    allocate_all();
#ifdef _FFTW
    allocate_conv();
#endif
    read_domain(directory);
    
    
    read_files(n,directory);
    //read_single_file(n,1,directory);
    /*
    if (cfl_dt < dt) {
        printf("Using cfl limited timestep of %.4e instead of %.4e\n",cfl_dt,dt);
        dt = cfl_dt;
    }
    */

    
    //move_to_com();


    read_stockholm(directory);
/*
    for(j=0;j<size_y;j++) {
        dens0[j] = params.mdot/(3*M_PI*Nu(ymed(j)));
        vy0[j] = -1.5*Nu(ymin(j))/ymin(j);
        vx0[j] = pow(ymed(j),-.5);
	    vx0[j] *= sqrt(1.0+pow(params.h,2)*pow(ymed(j),2*params.flaringindex)*
			  (2.0*params.flaringindex - 1.5));
        vx0[j] -= omf*ymed(j);

    }
*/
    omf0 = omf;

    output_stock(outputdir);

    output_init(outputdir);
   
    init_rk5();
    nsteps = 0;
    double tstart = psys[0].t;
    double tend = psys[0].t + time_step;
    double time = tstart;


#ifdef _FFTW
    printf("Using FFTW to compute fourier transforms\n");
#endif
#ifdef ARTIFICIALVISCOSITY
    printf("Using ARTIFICIAL VISCOSITY\n");
#endif
#ifdef FARGO
    printf("Using FARGO TRANSPORT\n");
#endif
#ifdef STOCKHOLMACC
    printf("Damping to current average in stockholm boundaries.\n");
#endif
#ifdef NOWAVEKILLRHO
    printf("Not damping the density.\n");
#endif
#ifdef FIXEDPSYS
    printf("Planetary system is FIXED.\n");
#endif
    printf("Starting time = %.3f\nEnding time = %.3f\n",tstart,tend);
    set_bc();

    set_dtLt(-1.0);

    while ((time < tend) && nsteps<MAXSTEPS) {
#ifdef FARGO
        compute_vmed(vx);
#endif
        dt = cfl();
        if (dt <= MINDT) {
            printf("Timestep has fallen below minimum!\n");
            break;
        }
        if (time + dt > tend) {
            dt = tend-time;
        }
        printf("\r t=%.3f, dt = %.8f,  %02d%% complete...",time,dt,(int)(100*(time-tstart)/(tend-tstart)));
        fflush(stdout);
        set_avg(0);
    
        potential();
        
#ifndef FIXEDPSYS
        move_planet();
#endif
    
       
        source_step();
       
        set_Lamex();
        viscosity();
#ifdef ARTIFICIALVISCOSITY
        
        artificial_visc();
#endif
    
        temp_to_vel();
        set_bc();
        
        vel_to_temp();
       
#ifdef FARGO
        compute_vmed(vx_temp);
#endif        
        transport_step();
        
        time += dt;
      
        set_avg(1);

        set_Lamdep();

        set_bc();
        stockholm();
        //output_psys(outputdir,nsteps);
        nsteps++;
    }
//    temp_to_vel();   

    set_dtLt(1.0);

    
    dt = tend - tstart;
    for(j=0;j<size_y;j++) {
        dbart[j]/=dt;
        //Lt[j]/=dt;
        //Lt[j+size_y]/=dt;
        //Ld[j]/=dt;
        //Ld[j+size_y]/=dt;
        //Lw[j]/=dt;
        //Lw[j+size_y]/=dt;
        Lw[j] = Lt[j] - Ld[j];
        Lw[j+size_y] = Lt[j+size_y] - Ld[j+size_y];
        drFt[j]/=dt;
        drFd[j]/=dt;
        drFdB[j]/=dt;
        mdotl[j]/=dt;
        drFnu[j]/=dt;
        drFw[j]/=dt;
        drFwB[j]/=dt;
        Lamex[j]/=dt;
        //Lamex[j+size_y]/=dt;
        Lamdep[j]/=dt;
        LamdepB[j]/=dt;
        dtLt[j]/=dt;
        dtLd[j]/=dt;
        dtLw[j]/=dt;
        dtLd_rhs[j]/=dt;
        mdotavg[j] /= dt;
        LamdepS[j + size_y*0] /= dt;
        LamdepS[j + size_y*1] /= dt;
        LamdepS[j + size_y*2] /= dt;
        LamdepS[j + size_y*3] /= dt;
        LamdepS[j + size_y*4] /= dt;
        LamdepS[j + size_y*5] /= dt;
        LamdepS[j + size_y*6] /= dt;
    }
    int mi;
    for(mi=1;mi<MMAX+2;mi++) {
        for(j=0;j<size_y;j++) {
            Lamex[j + size_y*mi] /= dt;
            drFd[j + size_y*mi] /= dt;
            dtLt[j + size_y*mi] /= dt;

        }
    }

    output(outputdir);
    output_torque(directory,n);
    output_torque(outputdir,n);
    free_rk5();
#ifdef _FFTW
    free_conv();
#endif
    //free_all();
    return 0;
}
