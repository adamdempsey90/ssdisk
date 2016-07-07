/** \file main.c

Main file of the distribution. 
Manages the call to initialization
functions, then the main loop.

*/
#include "fargo3d.h"

int   begin_i = 0, NbRestart = 0;
int   InnerOutputCounter=0, StillWriteOneOutput;
real dt;
real dtemp = 0.0;

int main(int argc, char *argv[]) {
  
  int   i, OutputNumber = 0, d;
  int   ni;
  char  ParameterFile[MAXLINELENGTH];
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i+=d) {
    d=1;
    if (*(argv[i]) == '+') {
      if (strspn (argv[i], \
		  "+S#D") \
	  != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], '#')) {
	d=2;
	ArrayNb = atoi(argv[i+1]);
	EarlyOutputRename = YES;
	if (ArrayNb <= 0) {
	  masterprint ("Incorrect Array number after +# flag\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'D')) {
	d=2;
	strcpy (DeviceFile, argv[i+1]);
	DeviceFileSpecified = YES;
      }
      if (strchr (argv[i], 'S')) {
	d=2;
	StretchNumber = atoi(argv[i+1]);
	StretchOldOutput = YES;
      }
    }
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], \
		  "-tofCmkspSVBD0#") \
	  != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strchr (argv[i], 'f'))
	ForwardOneStep = YES;
      if (strchr (argv[i], '0'))
	OnlyInit = YES;
      if (strchr (argv[i], 'C')) {
	EverythingOnCPU = YES;
#ifdef GPU
	mastererr ("WARNING: Forcing execution of all functions on CPU\n");
#else
	mastererr ("WARNING: Flag -C meaningless for a CPU built\n");
#endif
      }
      if (strchr (argv[i], 'm')) {
	Merge = YES;
      }
      if (strchr (argv[i], 'k')) {
	Merge = NO;
      }
      if (strchr (argv[i], 'o')) {
	RedefineOptions = YES;
	ParseRedefinedOptions (argv[i+1]) ;
	d=2;
      }
      if (strchr (argv[i], 's')) {
	Restart = YES;
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], '#')) {
	d=2;
	ArrayNb = atoi(argv[i+1]);
	if (ArrayNb <= 0) {
	  masterprint ("Incorrect Array number after -# flag\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'p')) {
	PostRestart = YES;
      }
      if (strchr (argv[i], 'S')) {
	Restart_Full = YES;
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'V')) {
	Dat2vtk = YES;
	Restart_Full = YES;
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect output number\n");
	  PrintUsage (argv[0]);	  
	}
      }
      if (strchr (argv[i], 'B')) {
	Vtk2dat = YES;
	Restart_Full = YES;
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect output number\n");
	  PrintUsage (argv[0]);	  
	}
      }
      if (strchr (argv[i], 'D')) {
	d=2;
	DeviceManualSelection = atoi(argv[i+1]);
      }
    }
    else strcpy (ParameterFile, argv[i]);
  }

#ifdef WRITEGHOSTS
  if (Merge == YES) {
	mastererr ("Cannot merge outputs when dumping ghost values.\n");
	mastererr ("'make nofulldebug' could fix this problem.\n");
	mastererr ("Using the -k flag could be another solution.\n");
	prs_exit (1);
  }
#endif
  
  if (ParameterFile[0] == 0) PrintUsage (argv[0]);

#ifdef MPICUDA
  EarlyDeviceSelection();
#endif
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);
  

#ifndef MPICUDA
  SelectDevice(CPU_Rank);
#endif
  InitVariables ();
  MPI_Barrier(MPI_COMM_WORLD);
  ReadDefaultOut ();
  ReadVarFile (ParameterFile);
  if (strcmp (PLANETCONFIG, "NONE") != 0)
    ThereArePlanets = YES;


  if (ORBITALRADIUS > 1.0e-30){
    YMIN *= ORBITALRADIUS;
    YMAX *= ORBITALRADIUS;
    DT   *= sqrt(ORBITALRADIUS*ORBITALRADIUS*ORBITALRADIUS);
  }


  SubsDef (OUTPUTDIR, DefaultOut);

  /* This must be placed ***BEFORE*** reading the input files in case of a restart */
  if ((ArrayNb) && (EarlyOutputRename == YES)) {
    i = strlen(OUTPUTDIR);
    if (OUTPUTDIR[i-1] == '/') OUTPUTDIR[i-1] = 0;//Remove trailing slash if any
    sprintf (OUTPUTDIR, "%s%06d/", OUTPUTDIR, ArrayNb); //Append numerical suffix
    /* There is no need to perform the wildcard (@) substitution. This has already been done */
    printf ("\n\n***\n\nNew Output Directory is %s\n\n***\n\n", OUTPUTDIR);
    MakeDir(OUTPUTDIR); /*Create the output directory*/
  }


  MakeDir(OUTPUTDIR); /*Create the output directory*/

#if !defined(X)
  NX = 1;
#endif
#if !defined(Y)
  NY = 1;
#endif
#if !defined(Z)
  NZ = 1;
#endif

  SelectWriteMethod();

#if !defined(Y) && !defined(Z)
  if (CPU_Rank==1){
    prs_error ("You cannot split a 1D mesh in x. Sequential runs only!");
  }
  if (CPU_Number > 1) {
    MPI_Finalize();
    prs_exit(EXIT_FAILURE);
  }
#endif
    
  ListVariables ("variables.par"); //Writes all variables defined in set up
  ListVariablesIDL ("IDL.var");
  ChangeArch(); /*Changes the name of the main functions
		  ChangeArch adds _cpu or _gpu if GPU is activated.*/
  split(&Gridd); /*Split mesh over PEs*/
  InitSpace();
  WriteDim();
  InitSurfaces();
  LightGlobalDev(); /* Copy light arrays to the device global memory */
  CreateFields(); // Allocate all fields.

  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  if(Corotating)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);
  OMEGAFRAME0 = OMEGAFRAME;
  /* We need to keep track of initial azimuthal velocity to correct
the target velocity in Stockholm's damping prescription. We copy the
value above *after* rescaling, and after any initial correction to
OMEGAFRAME (which is used afterwards to build the initial Vx field. */

  
  if(Restart == YES || Restart_Full == YES) {
    CondInit (); //Needed even for restarts: some setups have custom
		 //definitions (eg potential for setup MRI) or custom
		 //scaling laws (eg. setup planetesimalsRT).
    begin_i = RestartSimulation(NbRestart);
    if (ThereArePlanets) {
      PhysicalTime  = GetfromPlanetFile (NbRestart, 9, 0);
      OMEGAFRAME  = GetfromPlanetFile (NbRestart, 10, 0);
      RestartPlanetarySystem (NbRestart, Sys);
    }
  }
  else {
    if (ThereArePlanets)
      EmptyPlanetSystemFiles ();
    CondInit(); // Initialize set up
    // Note: CondInit () must be called only ONCE (otherwise some
    // custom scaling laws may be applied several times).
  }

  if (StretchOldOutput == YES) {
    StretchOutput (StretchNumber);
  }

  FARGO_SAFE(comm(ENERGY)); //Very important for isothermal cases!

  /* This must be placed ***after*** reading the input files in case of a restart */
  if ((ArrayNb) && (EarlyOutputRename == NO)) {
    i = strlen(OUTPUTDIR);
    if (OUTPUTDIR[i-1] == '/') OUTPUTDIR[i-1] = 0;//Remove trailing slash if any
    sprintf (OUTPUTDIR, "%s%06d/", OUTPUTDIR, ArrayNb); //Append numerical suffix
    /* There is no need to perform the wildcard (@) substitution. This has already been done */
    printf ("\n\n***\n\nNew Output Directory is %s\n\n***\n\n", OUTPUTDIR);
    MakeDir(OUTPUTDIR); /*Create the output directory*/
    ListVariables ("variables.par"); //Writes all variables defined in set up
    ListVariablesIDL ("IDL.var");
    InitSpace();
    WriteDim ();
  }
  
  DumpToFargo3drc(argc, argv);

  FillGhosts(PrimitiveVariables()); 
#ifdef STOCKHOLM 
  FARGO_SAFE(init_stockholm());
#ifdef STOCKHOLMACC
  FARGO_SAFE(ComputeVymed(Vy));
  FARGO_SAFE(ComputeRhomed(Density));
  
  Write2D(Density0_avg, "density0_2d_avg.dat", OUTPUTDIR, GHOSTINC);
  Write2D(Vy0_avg, "vy0_2d_avg.dat", OUTPUTDIR, GHOSTINC);
#endif
#endif

#ifdef GHOSTSX
  masterprint ("\n\nNew version with ghost zones in X activated\n");
#else
  masterprint ("Standard version with no ghost zones in X\n");
#endif
#ifdef TIMER
    clock_t begin_timer_time, end_timer_time;
    real timer_time_elapsed;
#endif
  for (ni = 0; ni<NITER; ni++) { // Iteration loop
      masterprintf("Start of %d iteration\n", ni);
      for (i = begin_i; i<=NTOT; i++) { // MAIN LOOP
    #ifdef TIMER
        if (i==begin_i) {
            begin_timer_time = clock();
        }
    #endif
        if (NINTERM * (TimeStep = (i / NINTERM)) == i) {

    #if defined(MHD) && defined(DEBUG)
          FARGO_SAFE(ComputeDivergence(Bx, By, Bz));
    #endif
          if (ThereArePlanets)
        WritePlanetSystemFile(TimeStep, NO);
          
    #ifndef NOOUTPUTS
          WriteOutputsAndDisplay(ALL);


          if(CPU_Master) printf("OUTPUTS %d at date t = %f OK\n", TimeStep, PhysicalTime);
    #endif

          if (TimeInfo == YES)
        GiveTimeInfo (TimeStep);
        }

        if (NSNAP != 0) {
          if (NSNAP * (TimeStep = (i / NSNAP)) == i) {
        WriteOutputsAndDisplay(SPECIFIC);
          }
        }
        
        AlgoGas();
        MonitorGlobal (MONITOR2D      | MONITORY | MONITORY_RAW|	\
               MONITORSCALAR  | MONITORZ | MONITORZ_RAW);
        
        if (ThereArePlanets) {
          WriteTorqueAndWork(TimeStep, 0);
          WritePlanetSystemFile(TimeStep, YES);
          SolveOrbits (Sys);
        }
    #ifdef TIMER
        if (i==begin_i) {
            end_timer_time = clock();
            timer_time_elapsed =( (double)(end_timer_time-begin_timer_time))/CLOCKS_PER_SEC;
            masterprint("time for time_step was %g s\n",timer_time_elapsed);
        }
    #endif

      }
      masterprintf("End of %d iteration\n", ni);
      compute_avgs();
      compute_steady_state();
      add_avg();
      output_steady_state();
  }

  MPI_Finalize();  
  printf("End of simulation!\n");
  return 0;  
}