#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef __APPLE__
#include <malloc.h>
#endif

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

#ifndef __FLOAT__
typedef double real;
#else
typedef float real;
#endif

struct quadruplet {
  real min;
  real max;
  real avg;
  real rms;
};

typedef struct quadruplet Quadruplet;

Quadruplet stat_info (real *buf, long size)
{
  real min=1e30;
  real max=-min;
  long i;
  real sum=0.0, avg;
  real rms=0.0;
  Quadruplet res;
  for (i = 0; i < size; i++) {
    if (buf[i] < min)
      min = buf[i];
    if (buf[i] > max)
      max = buf[i];
    sum += buf[i];
  }
  res.max = max;
  res.min = min;
  avg = res.avg = sum/(real)size;
  for (i = 0; i < size; i++) {
    rms += (buf[i]-avg)*(buf[i]-avg);
  }
  res.rms = sqrt(rms/(real)size);
  return res;
}

int main (int argc, char *argv[])
{
  char file1[1024];
  char file2[1024];
  FILE *f1, *f2;
  long sz1, sz2, i;
  real *buf1, *buf2, b1, b2;
  char *c1, *c2;
  int differ=0;
  size_t foo;
  Quadruplet res;
  int identical=0, flatratio=0, zerorequired=0;
  if (argc != 4) {
    fprintf (stderr, "%s admits exactly two arguments: -flag (z or f) file1 and file2\n", argv[0]);
    exit (1);
  }
  if (*(argv[1]+1) == 'z') zerorequired=1;
  snprintf (file1, 1020, "%s", argv[2]);
  snprintf (file2, 1020, "%s", argv[3]);

  f1 = fopen (file1, "r");
  if (f1 == NULL) {
    fprintf (stderr, "I cannot open %s.\n", file1);
    exit (1);
  }
  fseek (f1, 0, SEEK_END);
  sz1 = ftell(f1);
  fseek (f1, 0, SEEK_SET);

  f2 = fopen (file2, "r");
  if (f2 == NULL) {
    fprintf (stderr, "I cannot open %s.\n", file2);
    exit (1);
  }
  fseek (f2, 0, SEEK_END);
  sz2 = ftell(f2);
  fseek (f2, 0, SEEK_SET);
  if (sz1 != sz2) {
    fprintf (stderr, "Error : %s and %s have different sizes.\n", file1, file2);
    exit (1);
  }
  buf1 = (real *)malloc (sz1);
  buf2 = (real *)malloc (sz2);
  if ((buf1 == NULL) || (buf2 == NULL)) {
    fprintf (stderr, "Error : cannot allocate memory.\n");
    exit (1);
  }
  foo = fread (buf1, sz1, 1, f1);
  foo = fread (buf2, sz2, 1, f2);
  fclose (f1);
  fclose (f2);
  c1 = (char *)buf1;
  c2 = (char *)buf2;
  for (i = 0; i < sz1; i++) {
    if (c1[i] != c2[i])
      differ = 1;
  }
  if (differ) 
    printf ("%s and %s differ.\n", file1, file2);
  else {
    printf ("%s and %s are identical.\n", file1, file2);
    identical = 1;
  }
  sz1 /= sizeof(real);
  sz2 /= sizeof(real);
  res = stat_info (buf1, sz1);
  printf ("Size of file : %ld bytes, ie %ld ", sz1*sizeof(real), sz1);
  if (sizeof(real) == 4)
    printf ("single precision floating point values\n");
  else
    printf ("double precision floating point values\n");
  printf ("File1: min = %g\tmax = %g\tavg = %g\trms = %g\n", res.min, res.max, res.avg, res.rms);
  res = stat_info (buf2, sz1);
  printf ("File2: min = %g\tmax = %g\tavg = %g\trms = %g\n", res.min, res.max, res.avg, res.rms);
  for (i = 0; i < sz1; i++) {
    b1 = buf1[i];
    b2 = buf2[i];
    buf1[i] = b1-b2;
    buf2[i] = b2/b1;
  }
  res = stat_info (buf1, sz1);
  printf ("(1)-(2): min = %g\tmax = %g\tavg = %g\trms = %g\n", res.min, res.max, res.avg, res.rms);
  res = stat_info (buf2, sz1);
  printf ("Ratio 2/1: min = %g\tmax = %g\tavg = %g\trms = %g\n", res.min, res.max, res.avg, res.rms);
  if (fabs(res.rms/res.avg) < 1e-12) {
    printf ("flat ratio\n");
    flatratio = 1;
  }
  else
    printf ("variable ratio\n");
  for (i = 0; i < sz1; i++)
    buf2[i] = 1.0-buf2[i];
  res = stat_info (buf2, sz1);
  printf ("(1-2)/(1): min = %g\tmax = %g\tavg = %g\trms = %g\n", res.min, res.max, res.avg, res.rms);
  if ((flatratio && !zerorequired) || identical) {
    printf ("%sTEST PASSED (SUCCESS) : files ",BOLDGREEN);
    if (identical)
      printf ("are identical%s\n",RESET);
    else
      printf ("have a flat ratio%s\n",RESET);
  } else
    printf ("%sTEST FAILED : files differ%s\n",BOLDRED,RESET);
  free (buf1);
  free (buf2);
}
