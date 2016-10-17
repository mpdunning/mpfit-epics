/* file: fitFuncs.h
 * Standard fitting functions to be used with the mpfit routine.
 * 07/07/2014 zms, md
 *---------------------------------------------------------------------------*/

#ifndef _fitFuncs_h
#define _fitFuncs_h

typedef struct {
    double*	x;
    double*	y;
    double*	e;
    double*	f;
} fitdata_t;

int polyFit(int m, int n, double* p, double* dy, double** dv, void* pv);
int gausFit(int m, int n, double* p, double* dy, double** dv, void* pv);
double polyFunc(double x, double* p, int n);
double gausFunc(double x, double* p);

#endif	/* _fitFuncs_h */
