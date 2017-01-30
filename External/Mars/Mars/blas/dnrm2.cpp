/* dnrm2.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#include "f2c.hpp"

//#include <math.h>

double dnrm2_(int *n, double *x, int *incx)
{
    /* System generated locals */
    int i__1, i__2;
    double ret_val, d__1;

    /* Builtin functions */
    //double sqrt();

    /* Local variables */
    static double norm, scale, absxi;
    static int ix;
    static double ssq;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DNRM2 returns the euclidean norm of a vector via the function */
/*  name, so that */

/*     DNRM2 := sqrt( x'*x ) */

/*  Further Details */
/*  =============== */

/*  -- This version written on 25-October-1982. */
/*     Modified on 14-October-1993 to inline the call to DLASSQ. */
/*     Sven Hammarling, Nag Ltd. */

/*  =====================================================================
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n < 1 || *incx < 1) {
	   norm = 0.;
    } else if (*n == 1) {
	  norm = fabs(x[1]);
    } else {
	    scale = 0.;
	    ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK
*/
/*        auxiliary routine: */
/*        CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

	    i__1 = (*n - 1) * *incx + 1;
	    i__2 = *incx;
	   for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	      if (x[ix] != 0.) {
	       absxi = (d__1 = x[ix], fabs(d__1));
	       if (scale < absxi) {
/* Computing 2nd power */
	            d__1 = scale / absxi;
	            ssq = ssq * (d__1 * d__1) + 1.;
	           scale = absxi;
	        } else {
/* Computing 2nd power */
	            d__1 = absxi / scale;
	            ssq += d__1 * d__1;
	        }
        }
/* L10: */
    }
	    norm = scale * sqrt((double)ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of DNRM2. */

} /* dnrm2_ */
