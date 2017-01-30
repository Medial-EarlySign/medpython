/* dscal.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#include "f2c.hpp"

/* Subroutine */ int dscal_(int *n, double *da, double *dx, int *incx)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i, m, nincx, mp1;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*     DSCAL scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*  =====================================================================
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
    return 0;
    }
    if (*incx == 1) {

/*        code for increment equal to 1 */


/*        clean-up loop */

    m = *n % 5;
    if (m != 0) {
        i__1 = m;
        for (i = 1; i <= i__1; ++i) {
        dx[i] = *da * dx[i];
        }
        if (*n < 5) {
        return 0;
        }
    }
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
        dx[i] = *da * dx[i];
        dx[i + 1] = *da * dx[i + 1];
        dx[i + 2] = *da * dx[i + 2];
        dx[i + 3] = *da * dx[i + 3];
        dx[i + 4] = *da * dx[i + 4];
    }
    } else {

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
        dx[i] = *da * dx[i];
    }
    }
    return 0;
} /* dscal_ */
