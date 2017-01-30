/* ddot.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#include "f2c.hpp"

double ddot_(const int *n, double *dx, const int *incx, double *dy, const int *incy)
{
    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    static int i, m;
    static double dtemp;
    static int ix, iy, mp1;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     DDOT forms the dot product of two vectors. */
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
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
    return ret_val;
    }
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

    m = *n % 5;
    if (m != 0) {
        i__1 = m;
        for (i = 1; i <= i__1; ++i) {
        dtemp += dx[i] * dy[i];
        }
        if (*n < 5) {
        ret_val = dtemp;
        return ret_val;
        }
    }
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
        dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2]
            * dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i +
            4];
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
        dtemp += dx[ix] * dy[iy];
        ix += *incx;
        iy += *incy;
    }
    }
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */
