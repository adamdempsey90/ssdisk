#include <stdio.h>
#include <stdlib.h>

void WriteBinFile(int n1, int n2, int n3,	\
		  double *var1, char *name) {
  int i,j,k;
  int ntemp; 
  FILE *F;
  char filename[512];
  float *var;

  sprintf(filename, "%s%06d.bin", name, 0);
  F = fopen(filename,"w"); 
    
  var = (float*)malloc(sizeof(float)*n1*n2*n3);

  for(i=0; i<n1*n2*n3; i++) {
    var[i] = (float)var1[i];
  }
  
  ntemp = 12;
  fwrite(&ntemp,4,1,F);
  fwrite(&n1,4,1,F);
  fwrite(&n2,4,1,F);
  fwrite(&n3,4,1,F);
  fwrite(&ntemp,4,1,F);
  
  ntemp = n1*n2*n3*sizeof(float);
  fwrite(&ntemp,4,1,F); fwrite(var,sizeof(float)*n1*n2*n3,1,F); fwrite(&ntemp,4,1,F);
  fclose(F);
}

void read_subsample(double *q, int nx, int ny, int nz, FILE *file) {
  int i,j,k;
  for(k = 0; k<nz; k++) {
    for(j = 0; j<ny; j++) {
      fread(q+i+j*nx+k*nx*ny, sizeof(double), 1, file);
    }
  }
}

FILE *open_file(char *name, int n, char *m) {
  char filename[512];
  FILE *f;
  sprintf(filename, "%s%06d.dat", name, n);
  printf("Abriendo %s\n",filename);
  f = fopen(filename, m);
  if(f == NULL) {
    printf("Error leyendo %s\n", filename);
    return 0;
  }
  else printf("%s OK\n",filename);
  return f;
}

void div_mesh(double *q, int nx, int ny, int nz) {
  int i,j,k;
  double temp;
  for(k = 0; k<nz; k++) {
    for(j = 0; j<ny; j++) {
      temp = q[j*nx+k*nx*ny];
      for (i=1; i<nx; i++) {
	q[i+j*nx+k*nx*ny] = temp;
      }
    }
  }
}

void output(double *q, int nx, int ny, int nz, char *name) {
  FILE *f;
  char filename[512];
  int i,j,k;  
  sprintf(filename, "%s%06d.dat", name, 0);
  f = open_file(name, 0, "w");
  fwrite(q, sizeof(double)*nx*ny*nz, 1, f);
  printf("Subdivision de %s OK\n", filename);
  fclose(f);
}

int main(int argc, char *argv[]) {

  FILE *density;
  FILE *vx;
  FILE *vy;
  FILE *vz;
  FILE *energy;
  char name[512];
  
  int n;
  int nx;
  int ny;
  int nz;
  int stride;
  int dim;

  int i,j,k;

  double *Vx, *Vy, *Vz;
  double *Density;
  double *Energy;

  if(argc == 1) {
    printf("usage: ./div <n> <ny> <nz> <new nx dim> \n");
    exit(1);
  }
  printf("%d\n", argc);

  n = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  nx = atoi(argv[4]);
  
  printf("Your mesh will has a final size:\n nx = %d, ny = %d, nz = %d\n", nx, ny, nz);

  stride = nx*ny;

  Vx =     (double*)malloc(sizeof(double)*nx*ny*nz);
  Vy =     (double*)malloc(sizeof(double)*nx*ny*nz);
  Vz =     (double*)malloc(sizeof(double)*nx*ny*nz);
  Density= (double*)malloc(sizeof(double)*nx*ny*nz);
  Energy = (double*)malloc(sizeof(double)*nx*ny*nz);
  
  for(k=0; k<nz; k++) {
    for(j=0; j<ny; j++) {
      for(i=0; i<nx; i++) {
	Vx[i+j*nx+k*stride] =			\
	  Vy[i+j*nx+k*stride] =			\
	  Vz[i+j*nx+k*stride] =			\
	  Density[i+j*nx+k*stride] =		\
	  Energy[i+j*nx+k*stride] = 0;
      }
    }
  }
  
  density = open_file("Density", n, "r");
  vx = open_file("Vx", n, "r");
  vy = open_file("Vy", n, "r");
  vz = open_file("Vz", n, "r");
  energy = open_file("Energy", n, "r");
  
  printf("Comenzando la lectura...\n");
  
  read_subsample(Density, nx, ny, nz, density);
  printf("Lectura del archivo %s%.06d.dat OK\n", "Density", n);
  read_subsample(Vx, nx, ny, nz, vx);
  printf("Lectura del archivo %s%.06d.dat OK\n", "Vx", n);
  read_subsample(Vy, nx, ny, nz, vy);
  printf("Lectura del archivo %s%.06d.dat OK\n", "Vy", n);
  read_subsample(Vz, nx, ny, nz, vz);
  printf("Lectura del archivo %s%.06d.dat OK\n", "Vz", n);
  read_subsample(Energy, nx, ny, nz, energy);
  printf("Lectura del archivo %s%.06d.dat OK\n", "Energy", n);

  fclose(density);
  fclose(vx);
  fclose(vy);
  fclose(vz);
  fclose(energy);

  div_mesh(Density, nx, ny, nz);
  div_mesh(Vx, nx, ny, nz);
  div_mesh(Vy, nx, ny, nz);
  div_mesh(Vz, nx, ny, nz);
  div_mesh(Energy, nx, ny, nz);

  output(Density, nx, ny, nz, "Density");
  output(Vx, nx, ny, nz, "Vx");
  output(Vy, nx, ny, nz, "Vy");
  output(Vz, nx, ny, nz, "Vz");
  output(Energy, nx, ny, nz, "Energy");
  
  WriteBinFile(nx, ny, nz, Density, "Density");
  WriteBinFile(nx, ny, nz, Vx, "Vx");
  WriteBinFile(nx, ny, nz, Vy, "Vy");
  WriteBinFile(nx, ny, nz, Vz, "Vz");
  WriteBinFile(nx, ny, nz, Energy, "Energy");

  free(Vx);
  free(Vy);
  free(Vz);
  free(Energy);
  free(Density);

  printf("Final mesh size: nx=%d, ny=%d, nz=%d\n", nx, ny, nz);
  
  return 1;
}

