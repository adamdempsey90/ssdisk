#include "evolve.h"

void set_var(char *name,int int_val, double double_val, int bool_val, char *strval) {
    if (strcmp(name,"WKZOUT") == 0) {
        params.wkzout = double_val;
    }
    else if (strcmp(name,"XMIN") == 0) {	
        params.xmin = double_val;

    }
    else if (strcmp(name,"DT") == 0) {	
        params.dt = double_val;
    }
    else if (strcmp(name,"YMIN") == 0) {
        params.ymin = double_val;

    }
    else if (strcmp(name,"MDOT") == 0) {
        params.mdot = double_val;
    }
    else if (strcmp(name,"FLARINGINDEX") == 0) {
        params.flaringindex = double_val;
    }
    else if (strcmp(name,"XMAX") == 0) {
        params.xmax = double_val;
    }
    else if (strcmp(name,"WKZIN") == 0) {
        params.wkzin = double_val;

    }
    else if (strcmp(name,"YMAX") == 0) {
        params.ymax = double_val;

    }
    else if (strcmp(name,"THICKNESSSMOOTHING") == 0) {	
        params.soft = double_val;

    }
    else if (strcmp(name,"ALPHA") == 0) {
        params.alpha = double_val;

    }
    else if (strcmp(name,"ASPECTRATIO") == 0) {
        params.h = double_val;

    }
    else if (strcmp(name,"NINTERM") == 0) {
        params.ninterm = int_val;

    }
    else if (strcmp(name,"NTOT") == 0) {
        params.ntot = int_val;

    }
    else if (strcmp(name,"NX") == 0) {	
        params.nx = int_val;

    }
    else if (strcmp(name,"NY") == 0) {	
        params.ny= int_val;

    }
    else if (strcmp(name,"SPACING") == 0) {
        if (strval[0] == 'L') {
            if (strval[1] == 'O') {
                params.log = TRUE;
            }
            else {
                params.log = FALSE;
            }
        }
        else {
            params.log = FALSE;
        }

    }
    else if (strcmp(name,"FRAME") == 0) {
        if (strval[0] == 'F') {
            params.corotate = FALSE;
        }
        else {
            params.corotate = TRUE;
        }

    }
    else if (strcmp(name,"INDIRECTTERM") == 0) {
        params.indirect = bool_val;

    }
    else if (strcmp(name,"CFL") == 0) {
        params.cfl = double_val;
    }
    else if (strcmp(name,"NZ") == 0) {	
        params.nz = double_val;

    }

    return;
}
void read_param_file(char *directory) {
    FILE *f;

    char tok[20] = "\t :=>";

    char line[100],name[100],strval[100];
    char *data;
    double temp;
    int status;
    int int_val;
    int bool_val;
    char testbool;
    unsigned int i;
    char fname[256];
    sprintf(fname,"%svariables.par",directory);

    f= fopen(fname,"r");

    while (fgets(line,100,f)) {
       // printf("%s\n",line);
        status = sscanf(line,"%s",name);
      //  printf("%s\n",name);
        if (name[0] != '#' && status == 1) {
        
             data = line + (int)strlen(name);
             sscanf(data + strspn(data,tok),"%lf",&temp);
             sscanf(data + strspn(data,tok),"%s",strval);
             //printf("%lf\t%s\n",temp,strval);
            int_val = (int)temp;
            testbool = toupper(strval[0]);
            if (testbool == 'Y') bool_val = TRUE;
            else bool_val = FALSE;
            
            for (i = 0; i<strlen(name); i++) name[i] = (char)toupper(name[i]);
            
            set_var(name,int_val,temp,bool_val,strval);

        }
    }

    printf("%lg\t%s\n",params.wkzout,"WKZOUT");
    printf("%lg\t%s\n",params.xmin,"XMIN");
    printf("%lg\t%s\n",params.dt,"DT");
    printf("%lg\t%s\n",params.ymin,"YMIN");
    printf("%lg\t%s\n",params.mdot,"MDOT");
    printf("%lg\t%s\n",params.flaringindex,"FLARE");
    printf("%lg\t%s\n",params.xmax,"XMAX");
    printf("%lg\t%s\n",params.wkzin,"WKZIN");
    printf("%lg\t%s\n",params.ymax,"YMAX");
    printf("%lg\t%s\n",params.soft,"EPS");
    printf("%lg\t%s\n",params.alpha,"ALPHA");
    printf("%lg\t%s\n",params.h,"H");
    printf("%d\t%s\n",params.ninterm,"NINTERM");
    printf("%d\t%s\n",params.ntot,"NTOT");
    printf("%d\t%s\n",params.nx,"NX");
    printf("%d\t%s\n",params.ny,"NY");
    printf("%d\t%s\n",params.log,"LOG");
    printf("%d\t%s\n",params.corotate,"COROTATE");
    printf("%d\t%s\n",params.indirect,"INDIRECT");
    printf("%lg\t%s\n",params.cfl,"CFL");
    printf("%d\t%s\n",params.nz,"NZ");
    return;
}
