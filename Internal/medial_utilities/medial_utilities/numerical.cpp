// Utilities from numerical recipies

#include "medial_utilities.h"

// FFT.
#define SWAP(aa,bb) {tempr = (aa); (aa)=(bb); (bb) = tempr;} 

// Fast Fourier Transform on data
void four1 (double data[], unsigned long nn, int isign) {

	unsigned long n,mmax,m,j,istep,i ;
	double wtemp,wr,wpr,wpi,wi,theta ;
	double tempr,tempi ;

	n = nn << 1;
	j=1 ; 

	// Bit Reversal
	for (i=1; i<n; i+=2) {
		if (j>i) {
			SWAP(data[j],data[i]) ;
			SWAP(data[j+1],data[i+1]) ;
		}

		m=nn ;
		while (m>=2 && j>m) {
			j -= m ;
			m >>= 1 ;
		}
		j += m ;
	}

	// Danielson-Lanczos
	mmax = 2 ;
	while (n>mmax) {
		istep = mmax << 1 ;
		theta = isign * (6.28318530717959/mmax) ;
		wtemp = sin(0.5*theta) ;
		wpr = -2.0*wtemp*wtemp ;
		wpi = sin(theta) ;
		wr=1.0 ;
		wi=0.0 ;

		for (m=1; m<mmax; m+=2) {
			for (i=m; i<=n; i+=istep) {
				j=i+mmax; 
				tempr = wr*data[j] - wi*data[j+1] ;
				tempi = wr*data[j+1] + wi*data[j] ;
				data[j] = data[i] - tempr ;
				data[j+1] = data[i+1] - tempi ;
				data[i] += tempr ;
				data[i+1] += tempi ;
			}

			wr = (wtemp=wr)*wpr - wi*wpi+wr ;
			wi = wi*wpr + wtemp*wpi + wi ;
		}
		mmax=istep ;
	}
}

void get_r_phi(double *fft_data,int size, int i, double *r, double *phi) {

	double ap = fft_data[2*i+2] ;
	double bp = fft_data[2*i+3] ;
	double am = fft_data[2*size - 2*i - 2] ;
	double bm = fft_data[2*size - 2*i - 1] ;
	double x = ap+am ;
	double y = bm-bp ;
	
	*r = sqrt(x*x+y*y); 
	
	if(x==0)
		*phi = 3.1415/2 ;
	else
		*phi = atan(y/x) * 180/3.1415 ;

	return ;
}

//Gauss Jordan elimination to solve AX = b ;
void gaussj (double **a, int n, double **b, int m) {

	int *indxc = (int *) malloc (n*sizeof(int)) ;
	int *indxr = (int *) malloc (n*sizeof(int)) ;
	int *ipiv = (int *) malloc (n*sizeof(int)) ;
	assert (indxc != NULL && indxr != NULL && ipiv != NULL) ;

	for (int i=0; i<n; i++)
		ipiv[i] = 0 ;

	int irow,icol ;
	double big,pivinv,dum,tempr ;

	for (int i=0; i<n; i++) {
		big = 0.0 ;
		for (int j=0; j<n; j++) {
			if (ipiv[j] != 1) {
				for (int k=0; k<n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = fabs(a[j][k]) ;
							irow = j; 
							icol = k;
						}
					}
				}
			}
		}
		++(ipiv[icol]) ;

		if (irow != icol) {
			for (int l=0; l<n; l++)
				SWAP(a[irow][l],a[icol][l]) ;
			for (int l=0; l<m; l++)
				SWAP(b[irow][l],b[icol][l]) ;
		}

		indxr[i] = irow ;
		indxc[i] = icol ;
		assert(a[icol][icol] != 0.0) ;

		pivinv = 1.0/a[icol][icol] ;
		a[icol][icol] = 1.0 ;
		for (int l=0; l<n; l++)
			a[icol][l] *= pivinv ;
		for (int l=0; l<m; l++)
			b[icol][l] *= pivinv ;

		for (int ll=0; ll<n; ll++) {
			if (ll != icol) {
				dum = a[ll][icol] ;
				a[ll][icol] = 0.0 ;
				for (int l=0; l<n; l++)
					a[ll][l] -= a[icol][l]*dum ;
				for (int l=0; l<m; l++)
					b[ll][l] -= b[icol][l]*dum ;
			}
		}
	}

	for (int l=n-1; l>=0; l--) {
		if (indxr[l] != indxc[l])  {
			for (int k=0; k<n; k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]) ;
		}
	}
	
	free(ipiv) ;
	free(indxr) ;
	free(indxc) ;
	
	return ;
}