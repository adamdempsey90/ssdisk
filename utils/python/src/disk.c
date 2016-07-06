#include "evolve.h"
double Nu(double x) {
    return params.alpha*params.h*params.h*pow(x,2*params.flaringindex+.5);
}
double Cs(double x) {
    return params.h*pow(x,params.flaringindex-0.5);
}
