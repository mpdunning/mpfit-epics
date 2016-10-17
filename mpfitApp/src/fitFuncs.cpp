/* fitFuncs.cpp
 * Collection of fitting routines to be used with the mpfit function.
 *
 * 07/07/2014 zms, md
 *---------------------------------------------------------------------------*/

#include <math.h>
#include "fitFuncs.h"

int polyFit(int m, int n, double* p, double* dy, double** dv, void* pv) {
/*-----------------------------------------------------------------------------
 * n'th order polynomial fit function, where:
 * m - number of data points
 * n - number of parameters (2..10)
 * p - array of fit parameters
 * dy - array of residuals to be returned
 * dv - user computed derivatives, not used
 * pv - private data, in this case it is xydata_t type
 *---------------------------------------------------------------------------*/
    int i;
    fitdata_t* v = (fitdata_t*)pv;
    double* x; double* y; double* e; double* f; double fn = 0.0;

    x=v->x; y=v->y; e=v->e; f=v->f;
    
    for(i=0; i<m; i++){
      fn = polyFunc(x[i], p, n);
      dy[i] = (y[i] - fn)/e[i];
      f[i] = fn;
    }

    return 0;
}

double polyFunc(double x, double* p, int n) {
/*-----------------------------------------------------------------------------
 * Evaluates polynomial function value at x using n parameters in p.
 *---------------------------------------------------------------------------*/
    double fn, xx; 
    int i;

    for(i=0, fn=0.0, xx=1.0; i<n; i++){
      fn += p[i]*xx;
      xx *= x;
    }
    
    return(fn);
}

int gausFit(int m, int n, double *p, double* dy, double** dv, void* pv) {
/*-----------------------------------------------------------------------------
 * Gaussian fit function, where:
 * m - number of data points
 * n - number of parameters (4)
 * p - array of fit parameters:
 *     p[0] - constant offset
 *     p[1] - peak y value
 *     p[2] - x centroid position
 *     p[3] - gaussian sigma width
 * dy - array of residuals to be returned
 * dv - user computed derivatives, not used
 * pv - private data, in this case it is xydata_t type
 *---------------------------------------------------------------------------*/
    int i; 
    fitdata_t* v = (fitdata_t*)pv;
    double* x; double* y; double* e; double* f;
    double fx;

    x=v->x; y=v->y; e=v->e; f=v->f;
    
    for(i=0; i<m; i++){
      f[i] = fx = gausFunc(x[i], p);
      dy[i] = (y[i]-fx)/e[i];
    }
  
    return 0;
}

double gausFunc(double x, double* p) {
/*-----------------------------------------------------------------------------
 * Returns a y=p0+p1*exp(-0.5*(x-p2)*(x-p2)/(p3*p3))
 *---------------------------------------------------------------------------*/
  double f, xc, sig2;
  
    sig2 = p[3]*p[3];
    xc = x - p[2];
    f = p[0] + p[1]*exp(-0.5*xc*xc/sig2);
    
    return(f);
}


