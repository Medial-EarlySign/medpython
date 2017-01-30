#ifndef __BLAS__H__
#define __BLAS__H__

#include <stdlib.h>
#include <math.h>

int dqrdc2_(double *x, int *ldx, int *n, int *p, double *tol, int *rank, double *qraux, int *pivot, double *work);

int dqrsl_(double *x, int *ldx, int *n, int *k, double *qraux, double *y, double *qy, double *qty, double *b, double *rsd, double *xb, int *job, int *info);

int daxpy_(const int *n, const double *da, const double *dx, const int *incx, double *dy, const int *incy);

double ddot_(const int *n, double *dx, const int *incx, double *dy, const int *incy);

double dnrm2_(int *n, double *x, int *incx);

double d_sign(double *a, double *b);

int dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);

int dqrls_(double *x, int *n, int *p, double *y, int *ny, double *tol, double *b, double *rsd, double *qty, int *k, int *jpvt, double *qraux, double *work);

int dscal_(int *n, double *da, double *dx, int *incx);

int dtrsl_(double *t, int *ldt, int *n, double *b, int *job, int *info);


#endif
