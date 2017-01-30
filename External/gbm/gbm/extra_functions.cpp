// Extra functions from R //

#include "stdafx.h"
#include <R_ext/Random.h>

#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static Int32 dummy[625];
static Int32 *mt = dummy+1; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

static void
MT_sgenrand(Int32 seed)
{
    int i;

    for (i = 0; i < N; i++) {
	mt[i] = seed & 0xffff0000;
	seed = 69069 * seed + 1;
	mt[i] |= (seed & 0xffff0000) >> 16;
	seed = 69069 * seed + 1;
    }
    mti = N;
}

static double MT_genrand(void)
{
    Int32 y;
    static Int32 mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

 //   mti = dummy[0];

    if (mti >= N) { /* generate N words at one time */
	int kk;

	if (mti == N+1)   /* if sgenrand() has not been called, */
	    MT_sgenrand(4357); /* a default initial seed is used   */

	for (kk = 0; kk < N - M; kk++) {
	    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
	    mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	for (; kk < N - 1; kk++) {
	    y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
	    mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

	mti = 0;
    }

    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
 //   dummy[0] = mti;

    return ( (double)y * 2.3283064365386963e-10 ); /* reals: [0,1)-interval */
}

#define i2_32m1 2.328306437080797e-10/* = 1/(2^32 - 1) */

static double fixup(double x)
{
    /* ensure 0 and 1 are never returned */
    if(x <= 0.0) return 0.5*i2_32m1;
    if((1.0 - x) <= 0.0) return 1.0 - 0.5*i2_32m1;
    return x;
}

double unif_rand(void)
{
	return fixup(MT_genrand());
}

// Sorting ..
static int rcmp(double x, double y, Rboolean nalast)
{
//    int nax = ISNAN(x), nay = ISNAN(y);
//    if (nax && nay)     return 0;
//    if (nax)            return nalast?1:-1;
//    if (nay)            return nalast?-1:1;

	if (x < y)          return -1;
    if (x > y)          return 1;
    return 0;
}

void my_rsort_with_index(double *x, int *indx, int n)
{
    double v;
    int i, j, h, iv;

    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
        for (i = h; i < n; i++) {
            v = x[i]; iv = indx[i];
            j = i;
            while (j >= h && rcmp(x[j - h], v, TRUE) > 0)
                 { x[j] = x[j - h]; indx[j] = indx[j-h]; j -= h; }
            x[j] = v; indx[j] = iv;
        }
}

double my_fmax2(double x, double y)
{
        return (x < y) ? y : x;
}

double my_fmin2(double x, double y)
{
        return (x < y) ? x : y;
}