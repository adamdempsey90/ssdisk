#include "fargo3d.h"

void output_steady_state(int n) {
    FILE *f;
    char fname[256];

    real *drfnu = drFluxVisc->field_cpu;
    real *drfd = drFluxDisk->field_cpu;
    real *dtld = dtLDisk->field_cpu;
    real *density_ss = DensitySS->field_cpu;
    real *vy_ss = VySS->field_cpu;
    real *vx_ss= VxSS->field_cpu;

    


    sprintf(fname,"%sss_fluxes%d.dat",OUTPUTDIR,n);

    f = fopen(fname,"w");

    fwrite(&drfnu[NGHY],sizeof(real),Ny,f);
    fwrite(&drfd[NGHY],sizeof(real),Ny,f);
    fwrite(&dtld[NGHY],sizeof(real),Ny,f);
    fwrite(&density_ss[NGHY],sizeof(real),Ny,f);
    fwrite(&vy_ss[NGHY],sizeof(real),Ny,f);
    fwrite(&vx_ss[NGHY],sizeof(real),Ny,f);

    fclose(f);

    sprintf(fname,"%sfinal%d.dat",OUTPUTDIR,n);

    f = fopen(fname,"w");

    real *dens = Density->field_cpu;
    fwrite(&dens[NGHY*Nx],sizeof(real),Ny*Nx,f);

    fclose(f);

    DumpAllFields(0,n+1);


    

    return;
}   
