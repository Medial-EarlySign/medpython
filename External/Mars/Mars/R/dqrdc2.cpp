/* dqrdc2.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#include "../blas/f2c.hpp"

/* Table of constant values */

//#include <math.h>

static int c__1 = 1;


/*     dqrdc2 uses householder transformations to compute the qr */
/*     factorization of an n by p matrix x.  a limited column */
/*     pivoting strategy based on the 2-norms of the reduced columns */
/*     moves columns with near-zero norm to the right-hand edge of */
/*     the x matrix.  this strategy means that sequential one */
/*     degree-of-freedom effects can be computed in a natural way. */

/*     i am very nervous about modifying linpack code in this way. */
/*     if you are a computational linear algebra guru and you really */
/*     understand how to solve this problem please feel free to */
/*     suggest improvements to this code. */

/*     Another change was to compute the rank. */

/*     on entry */

/*        x       double precision(ldx,p), where ldx .ge. n. */
/*                x contains the matrix whose decomposition is to be */
/*                computed. */

/*        ldx     integer. */
/*                ldx is the leading dimension of the array x. */

/*        n       integer. */
/*                n is the number of rows of the matrix x. */

/*        p       integer. */
/*                p is the number of columns of the matrix x. */

/*        tol     double precision */
/*                tol is the nonnegative tolerance used to */
/*                determine the subset of the columns of x */
/*                included in the solution. */

/*        jpvt    integer(p). */
/*                integers which are swapped in the same way as the */
/*                the columns of x during pivoting.  on entry these */
/*                should be set equal to the column indices of the */
/*                columns of the x matrix (typically 1 to p). */

/*        work    double precision(p,2). */
/*                work is a work array. */

/*     on return */

/*        x       x contains in its upper triangle the upper */
/*                triangular matrix r of the qr factorization. */
/*                below its diagonal x contains information from */
/*                which the orthogonal part of the decomposition */
/*                can be recovered.  note that if pivoting has */
/*                been requested, the decomposition is not that */
/*                of the original matrix x but that of x */
/*                with its columns permuted as described by jpvt. */

/*        k       integer. */
/*                k contains the number of columns of x judged */
/*                to be linearly independent. */

/*        qraux   double precision(p). */
/*                qraux contains further information required to recover */
/*                the orthogonal part of the decomposition. */

/*        jpvt    jpvt(k) contains the index of the column of the */
/*                original matrix that has been interchanged into */
/*                the k-th column. */

/*     original (dqrdc.f) linpack version dated 08/14/78 . */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     this version dated 22 august 1995 */
/*     ross ihaka */

/*     bug fixes 29 September 1999 BDR (p > n case, inaccurate ranks) */


/*     dqrdc2 uses the following functions and subprograms. */

/*     blas daxpy,ddot,dscal,dnrm2 */
/*     fortran dabs,dmax1,min0,dsqrt */

/* Subroutine */ int dqrdc2_(double *x, int *ldx, int *n, int *p, double *tol, int *k, double *qraux, int *jpvt, double *work)
{
    /* System generated locals */
    int x_dim1, x_offset, work_dim1, work_offset, i__1, i__2, i__3;
    double d__1, d__2;

    /* Builtin functions */
    //double d_sign(), sqrt();

    /* Local variables */
    //extern double ddot_(), dnrm2_();
    static int i, j, l;
    static double t;
    //extern /* Subroutine */ int dscal_(), daxpy_();
    static double nrmxl, tt;
    static int lp1, lup;
    static double ttt;

/*     internal variables */

/*     compute the norms of the columns of x. */

    /* Parameter adjustments */
    work_dim1 = *p;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    --jpvt;
    --qraux;
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
    qraux[j] = dnrm2_(n, &x[j * x_dim1 + 1], &c__1);
    work[j + work_dim1] = qraux[j];
    work[j + (work_dim1 << 1)] = qraux[j];
    if (work[j + (work_dim1 << 1)] == 0.) {
        work[j + (work_dim1 << 1)] = 1.;
    }
/* L70: */
    }

/*     perform the householder reduction of x. */

    lup = min(*n,*p);
    *k = *p + 1;
    i__1 = lup;
    for (l = 1; l <= i__1; ++l) {

/*     previous version only cycled l to lup */

/*     cycle the columns from l to p left-to-right until one */
/*     with non-negligible norm is located.  a column is considered */
/*     to have become negligible if its norm has fallen below */
/*     tol times its original norm.  the check for l .le. k */
/*     avoids infinite cycling. */

L80:
    if (l >= *k || qraux[l] >= work[l + (work_dim1 << 1)] * *tol) {
        goto L120;
    }
    lp1 = l + 1;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i) {
        t = x[i + l * x_dim1];
        i__3 = *p;
        for (j = lp1; j <= i__3; ++j) {
        x[i + (j - 1) * x_dim1] = x[i + j * x_dim1];
/* L90: */
        }
        x[i + *p * x_dim1] = t;
/* L100: */
    }
    i = jpvt[l];
    t = qraux[l];
    tt = work[l + work_dim1];
    ttt = work[l + (work_dim1 << 1)];
    i__2 = *p;
    for (j = lp1; j <= i__2; ++j) {
        jpvt[j - 1] = jpvt[j];
        qraux[j - 1] = qraux[j];
        work[j - 1 + work_dim1] = work[j + work_dim1];
        work[j - 1 + (work_dim1 << 1)] = work[j + (work_dim1 << 1)];
/* L110: */
    }
    jpvt[*p] = i;
    qraux[*p] = t;
    work[*p + work_dim1] = tt;
    work[*p + (work_dim1 << 1)] = ttt;
    --(*k);
    goto L80;
L120:
    if (l == *n) {
        goto L190;
    }

/*           compute the householder transformation for column l. */

    i__2 = *n - l + 1;
    nrmxl = dnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
    if (nrmxl == 0.) {
        goto L180;
    }
    if (x[l + l * x_dim1] != 0.) {
        nrmxl = d_sign(&nrmxl, &x[l + l * x_dim1]);
    }
    i__2 = *n - l + 1;
    d__1 = 1. / nrmxl;
    dscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
    x[l + l * x_dim1] += 1.;

/*              apply the transformation to the remaining columns, */
/*              updating the norms. */

    lp1 = l + 1;
    if (*p < lp1) {
        goto L170;
    }
    i__2 = *p;
    for (j = lp1; j <= i__2; ++j) {
        i__3 = *n - l + 1;
        t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &c__1) / x[l + l * x_dim1];
        i__3 = *n - l + 1;
        daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &c__1);
        if (qraux[j] == 0.) {
        goto L150;
        }
/* Computing 2nd power */
        d__2 = (d__1 = x[l + j * x_dim1], fabs(d__1)) / qraux[j];
        tt = 1. - d__2 * d__2;
        tt = max(tt,0.);
        t = tt;

/* modified 9/99 by BDR. Re-compute norms if there is large reduct
ion */
/* The tolerance here is on the squared norm */
/* In this version we need accurate norms, so re-compute often. */
/* work(j,1) is only updated in one case: looks like a bug -- no l
onger used*/

/*                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j,1))
**2 */
/*                     if (tt .eq. 1.0d0) go to 130 */
        if (fabs(t) < 1e-6) {
        goto L130;
        }
        qraux[j] *= sqrt(t);
        goto L140;
L130:
        i__3 = *n - l;
        qraux[j] = dnrm2_(&i__3, &x[l + 1 + j * x_dim1], &c__1);
        work[j + work_dim1] = qraux[j];
L140:
L150:
/* L160: */
        ;
    }
L170:

/*              save the transformation. */

    qraux[l] = x[l + l * x_dim1];
    x[l + l * x_dim1] = -nrmxl;
L180:
L190:
/* L200: */
    ;
    }
/* Computing MIN */
    i__1 = *k - 1;
    *k = min(i__1,*n);
    return 0;
} /* dqrdc2_ */
