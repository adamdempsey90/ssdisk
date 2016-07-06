#include <Python.h>
#include <numpy/ndarraytypes.h>

//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>

typedef double real;

struct field {
  real *data;
  real *x;
  real *y;
  int nx;
  int ny;
};

typedef struct field Field;

double interp2d(real x1, real x2,
		real y1, real y2,
		real f11, real f12, real f21, real f22,
		real x, real y) {

  double t = (x-x1)/(x2-x1);
  double u = (y-y1)/(y2-y1);

  return (1.0-t)*(1.0-u)*f11 + t*(1.0-u)*f12 + t*u*f22 + u*(1.0-t)*f21;
}

double get_value(Field *f, real x0, real y0) {
  
  int i, j;
  real value;
  int nx, ny;

  nx = f->nx;  ny = f->ny;

  i = (int)((x0-f->x[0])/(f->x[nx-1]-f->x[0])*nx);
  j = (int)((y0-f->y[0])/(f->y[ny-1]-f->y[0])*ny);

  if (i<0 || j<0 || i>(nx - 2) || j>(ny - 2)) { //This line is for the D-plots. 
    //    fprintf(stderr,"i=%d\tj=%d\tnx=%d\tny=%d\tYou are trying to get a value outside the mesh...\n",i,j,nx,ny);
    return -1.0;
  }

  value = interp2d(f->x[i], f->x[i+1],
		   f->y[j], f->y[j+1],
		   f->data[i+j*nx], f->data[(i+1)+j*nx],
		   f->data[i+(j+1)*nx], f->data[(i+1)+(j+1)*nx],
		   x0,y0);
  return value;
}

int euler(Field *Vx, Field *Vy, real *x, real *y, int sign) {

  real vphi;
  real vrad;
  real l, h;

  vphi = get_value(Vx,*x,*y);
  vrad = get_value(Vy,*x,*y);

  if (vphi == -1.0 || vrad == -1.0) {
    //    fprintf(stderr,"There was an error in euler, vphi or vrad is -1\n");
    return -1;
  }
  
//  if ( (Vx->x[Vx->nx-1]-Vx->x[0])/Vx->nx < (Vx->y[Vx->ny-1]-Vx->y[0])/Vx->ny )
//    l = (Vx->x[Vx->nx-1]-Vx->x[0])/Vx->nx;
//  else
  l = (Vx->y[Vx->ny-1]-Vx->y[0])/Vx->ny;

  h = 0.1*l/sqrt(vphi*vphi+vrad*vrad);

  *x += sign*h*vphi/(*y);
  *y += sign*h*vrad;

  return 0;
}

int get_stream(Field *Vx, Field *Vy, real x0, real y0, real *x, real *y, int nmax, real lmax, int sign) {

  int i,k;
  real xold, yold, l;
  real dx, dy;
  int control;

  k = 0;
  l = 0.0;
  xold = y0*cos(x0);
  yold = y0*sin(x0);
  for (i=0;i<nmax;i++) {
    x[i] = x0;
    y[i] = y0;
    dx = y0*cos(x0)-xold;
    dy = y0*sin(x0)-yold;
    l += sqrt(dx*dx+dy*dy);
    xold = y0*cos(x0);
    yold = y0*sin(x0);
    control = euler(Vx,Vy,&x0,&y0,sign);
    if (control == -1)
      break;
    if (l>lmax) {
      break;
    }
    k+=1;
  }
  return k;
}

Field *create_fields(real *field, real *domx, real *domy, int nx, int ny) {

  Field *f;
  
  f = (Field*)malloc(sizeof(Field));

  f->data = field;
  f->x = domx;
  f->y = domy;  
  f->nx = nx;
  f->ny = ny;

  return f;
}

static PyObject *stream(PyObject* self, PyObject* args) {
  
  PyArrayObject *Vx, *Vx_x, *Vx_y;
  PyArrayObject *Vy, *Vy_x, *Vy_y;

  real *vx, *vx_x, *vx_y;
  real *vy, *vy_x, *vy_y;

  int nx_x, ny_x;
  int nx_y, ny_y;

  real *x, *y;

  int nmax = 10000;
  int sign = 1;
  real lmax = 2.0*M_PI;

  int i,k;

  Field *Field_Vx;
  Field *Field_Vy;

  real x0,y0;

  if (!PyArg_ParseTuple(args, "OOOOOOdd|idi", &Vx, &Vx_x, &Vx_y, &Vy, &Vy_x, &Vy_y, &x0, &y0, &nmax, &lmax, &sign)) //Getting the numpy arrays
    return NULL;

  x = (real*)malloc(sizeof(real)*nmax);
  y = (real*)malloc(sizeof(real)*nmax);

  nx_x = (int)Vx_x->dimensions[0];
  ny_x = (int)Vx_y->dimensions[0];

  nx_y = (int)Vy_x->dimensions[0];
  ny_y = (int)Vy_y->dimensions[0];

  vx = (real*)Vx->data;
  vy = (real*)Vy->data;

  vx_x = (real*)Vx_x->data;
  vy_x = (real*)Vy_x->data;

  vx_y = (real*)Vx_y->data;
  vy_y = (real*)Vy_y->data;

  Field_Vx = create_fields(vx, vx_x, vx_y, nx_x, ny_x);
  Field_Vy = create_fields(vy, vy_x, vy_y, nx_y, ny_y);

  k = get_stream(Field_Vx, Field_Vy, x0, y0, x, y, nmax, lmax, sign);

  // Output region

  PyObject *x_out = PyList_New(k);
  PyObject *y_out = PyList_New(k);

  for (i=0;i<k;i++) {
    PyList_SET_ITEM(x_out, i, Py_BuildValue("d",(real)x[i]));
    PyList_SET_ITEM(y_out, i, Py_BuildValue("d",(real)y[i]));
  }

  return Py_BuildValue("OO", x_out, y_out);
}

static PyMethodDef StreamMethods[] = {
  {"stream", stream, METH_VARARGS, "Method for computing a streamline."},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initstream(void) {
  (void) Py_InitModule("stream", StreamMethods);
}
