/* daxpy.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#include "f2c.hpp"
#include "../blas.hpp"
/* Subroutine */ int daxpy_(const int *n, const double *da, const double *dx, const int *incx, double *dy, const int *incy)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i, m, ix, iy, mp1;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     DAXPY constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*  =====================================================================
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0)
		return 0;

    if (*da == 0.)
		return 0;

    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

	    m = *n % 4;
		if (m != 0) {
			i__1 = m;
			for (i = 1; i <= i__1; ++i) {
			dy[i] += *da * dx[i];
			}
		}

		if (*n < 4) {
			return 0;
		}
		mp1 = m + 1;
		i__1 = *n;
		for (i = mp1; i <= i__1; i += 4) {
			dy[i] += *da * dx[i];
			dy[i + 1] += *da * dx[i + 1];
			dy[i + 2] += *da * dx[i + 2];
			dy[i + 3] += *da * dx[i + 3];
		}
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
        dy[iy] += *da * dx[ix];
        ix += *incx;
        iy += *incy;
    }
    }
    return 0;
} /* daxpy_ */
