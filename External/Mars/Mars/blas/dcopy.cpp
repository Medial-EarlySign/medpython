/* dcopy.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#include "f2c.hpp"

/* Subroutine */ int dcopy_(int *n, double *dx, int *incx, double *dy, int *incy)
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

/*     DCOPY copies a vector, x, to a vector, y. */
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
    if (*n <= 0) {
		return 0;
    }

    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

    m = *n % 7;
    if (m != 0) {
        i__1 = m;
        for (i = 1; i <= i__1; ++i) {
        dy[i] = dx[i];
        }
        if (*n < 7) {
        return 0;
        }
    }
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 7) {
        dy[i] = dx[i];
        dy[i + 1] = dx[i + 1];
        dy[i + 2] = dx[i + 2];
        dy[i + 3] = dx[i + 3];
        dy[i + 4] = dx[i + 4];
        dy[i + 5] = dx[i + 5];
        dy[i + 6] = dx[i + 6];
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
        dy[iy] = dx[ix];
        ix += *incx;
        iy += *incy;
    }
    }
    return 0;
} /* dcopy_ */
