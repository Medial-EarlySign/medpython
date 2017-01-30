//Linear-regression functions

#include "classifiers.h"
#include "ROCR/ROCR/ROCR.h"
#include "medial_utilities/medial_utilities/medial_utilities.h"
#include "classifiers.h"

#define _CRT_SECURE_NO_WARNINGS

#define MY_EPSILON	1e-5
//-----------------------------------------------------------------------------
int lm_wrapper(double *x, double *y, int nsamples, int nftrs, int niter, double eiter , double rfactor, double *b, double *err)
{
	double *xn=NULL, *yn=NULL, *bn=NULL;
	double *avg = NULL, *std = NULL;
	double avgy, stdy;

	xn = (double *)new double[nftrs*nsamples];
	yn = (double *)new double[nsamples];
	bn = (double *)new double[nftrs];
	avg = (double *)new double[nftrs];
	std = (double *)new double[nftrs];

	if (xn == NULL || yn == NULL || bn == NULL) {
		if (xn != NULL) delete xn;
		if (yn != NULL) delete yn;
		if (bn != NULL) delete bn;
		if (avg != NULL) delete avg;
		if (std != NULL) delete std;
		fprintf(stderr,"lm_wrapper: Can't allocate working area\n");
		return -1;
	}

	int i,j;
	double sx,sxx,sn;
	sn = 1.0/(double)nsamples;
	for (i=0; i<nftrs; i++) {
		sx = 0;
		sxx = 0;
		for (j=0; j<nsamples; j++) {
			sx += x[i*nsamples+j]*sn;
			sxx += x[i*nsamples+j]*x[i*nsamples+j]*sn;
		}
		avg[i] = sx; //sx/(double)nftrs;
		std[i] = sxx - avg[i]*avg[i]; //sxx/(double)nftrs - avg[i]*avg[i];
		std[i] = sqrt(abs(std[i]));
		//fprintf(stderr,"i=%d nsamples=%d nftrs %d sn %f sx %f sxx %f avg %f std %f\n",i,nsamples,nftrs,sn,sx,sxx,avg[i],std[i]);
		if (std[i] < MY_EPSILON)
			std[i] = MY_EPSILON;
		//fprintf(stderr,"i=%d nsamples=%d nftrs %d sn %f sx %f sxx %f avg %f std %f\n",i,nsamples,nftrs,sn,sx,sxx,avg[i],std[i]);
		for (j=0; j<nsamples; j++) {
			xn[i*nsamples+j] = (x[i*nsamples + j] - avg[i])/std[i];
		}
	}

	sx = 0;
	sxx = 0;
	for (j=0; j<nsamples; j++) {
		sx += y[j];
		sxx += y[j]*y[j];
	}
	avgy = sx/(double)nsamples;
	stdy = sxx/(double)nsamples - avgy*avgy;
	stdy = sqrt(stdy);
	if (stdy < MY_EPSILON)
		stdy = MY_EPSILON;
	for (j=0; j<nsamples; j++)
		yn[j] = (y[j] - avgy)/stdy;

//	fprintf(stderr,"before internal lm\n"); fflush(stderr);
	int rclm = lm (xn, yn, nsamples, nftrs, niter, eiter , rfactor, bn, err) ;
//	fprintf(stderr,"after internal lm rclm=%d err=%f nsamples=%d nftrs=%d\n",rclm,*err,nsamples,nftrs); fflush(stderr);

	if (rclm >= 0) {
		b[nftrs] = avgy;
		for (i=0; i<nftrs; i++) {
//			fprintf(stderr,"i=%d bn %f avgy %f stdy %f avgi %f stdi %f\n",i,bn[i],avgy,stdy,avg[i],std[i]); fflush(stderr);
			b[i] = bn[i]*stdy/std[i];
			b[nftrs] += -bn[i]*(avg[i]/std[i])*stdy;
		}

		//for (i=0; i<nftrs; i++) {
		//	fprintf(stderr,"i=%d bn %f b %f avgy %f stdy %f avgi %f stdi %f\n",i,bn[i],b[i],avgy,stdy,avg[i],std[i]);
		//}


	}

	if (xn != NULL) delete xn;
	if (yn != NULL) delete yn;
	if (bn != NULL) delete bn;
	if (avg != NULL) delete avg;
	if (std != NULL) delete std;

	return rclm;
}


//-----------------------------------------------------------------------------
int lm_wrapper(float *x, float *y, int nsamples, int nftrs, int niter, double eiter , double rfactor, float *b, double *err)
{
	double *xn=NULL, *yn=NULL, *bn=NULL;
	double *avg = NULL, *std = NULL;
	double avgy, stdy;

	xn = (double *)new double[nftrs*nsamples];
	yn = (double *)new double[nsamples];
	bn = (double *)new double[nftrs];
	avg = (double *)new double[nftrs];
	std = (double *)new double[nftrs];

	if (xn == NULL || yn == NULL || bn == NULL) {
		if (xn != NULL) delete xn;
		if (yn != NULL) delete yn;
		if (bn != NULL) delete bn;
		if (avg != NULL) delete avg;
		if (std != NULL) delete std;
		fprintf(stderr,"lm_wrapper: Can't allocate working area\n");
		return -1;
	}



	int i,j;

	//for (i=0; i<10; i++) {
	//	fprintf(stderr,"## ");
	//	for (j=0; j<nftrs; j++)
	//		fprintf(stderr,"%7.2f ",x[i*nsamples+j]);
	//	fprintf(stderr," :: %7.2f\n",y[i]);

	//}

	double sx,sxx,sn;
	sn = 1.0/(double)nsamples;
	for (i=0; i<nftrs; i++) {
		sx = 0;
		sxx = 0;
		for (j=0; j<nsamples; j++) {
			sx += (double)x[i+nftrs*j]*sn;
			sxx += (double)x[i+nftrs*j]*(double)x[i+nftrs*j]*sn;
		}
		avg[i] = sx; //sx/(double)nftrs;
		std[i] = sxx - avg[i]*avg[i]; //sxx/(double)nftrs - avg[i]*avg[i];
		std[i] = sqrt(abs(std[i]));
		//fprintf(stderr,"i=%d nsamples=%d nftrs %d sn %f sx %f sxx %f avg %f std %f\n",i,nsamples,nftrs,sn,sx,sxx,avg[i],std[i]);
		if (std[i] < MY_EPSILON)
			std[i] = MY_EPSILON;
		//fprintf(stderr,"i=%d nsamples=%d nftrs %d sn %f sx %f sxx %f avg %f std %f\n",i,nsamples,nftrs,sn,sx,sxx,avg[i],std[i]);
		for (j=0; j<nsamples; j++) {
			xn[i*nsamples+j] = ((double)x[i+nftrs*j] - avg[i])/std[i];
		}
	}

	sx = 0;
	sxx = 0;
	for (j=0; j<nsamples; j++) {
		sx += (double)y[j];
		sxx += (double)y[j]*(double)y[j];
	}
	avgy = sx/(double)nsamples;
	stdy = sxx/(double)nsamples - avgy*avgy;
	stdy = sqrt(stdy);
	if (stdy < MY_EPSILON)
		stdy = MY_EPSILON;
	for (j=0; j<nsamples; j++)
		yn[j] = ((double)y[j] - avgy)/stdy;

//	fprintf(stderr,"before internal lm\n"); fflush(stderr);
	int rclm = lm (xn, yn, nsamples, nftrs, niter, eiter , rfactor, bn, err) ;
//	fprintf(stderr,"after internal lm rclm=%d err=%f nsamples=%d nftrs=%d\n",rclm,*err,nsamples,nftrs); fflush(stderr);

	if (rclm >= 0) {
		b[nftrs] = (float)avgy;
		for (i=0; i<nftrs; i++) {
//			fprintf(stderr,"i=%d bn %f avgy %f stdy %f avgi %f stdi %f\n",i,bn[i],avgy,stdy,avg[i],std[i]); fflush(stderr);
			b[i] = (float)(bn[i]*stdy/std[i]);
			b[nftrs] += (float)(-bn[i]*(avg[i]/std[i])*stdy);
		}

		//for (i=0; i<nftrs; i++) {
		//	fprintf(stderr,"i=%d bn %f b %f avgy %f stdy %f avgi %f stdi %f\n",i,bn[i],b[i],avgy,stdy,avg[i],std[i]);
		//}


	}

	if (xn != NULL) delete xn;
	if (yn != NULL) delete yn;
	if (bn != NULL) delete bn;
	if (avg != NULL) delete avg;
	if (std != NULL) delete std;

	return rclm;
}
// Learn		
int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double rfactor, double *b, double *err) {

	// Prepare
	double *y ;
	if ((y = (double *) malloc (nsamples*sizeof (double))) == NULL)
		return -1 ;

	memcpy(y,_y,nsamples*sizeof(double)) ;	
	memset(b,0,nftrs*sizeof(double)) ;

	double *sumxx ;
	if ((sumxx = (double *) malloc (nftrs*sizeof(double)))==NULL)
		return -1 ;
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += x[XIDX(j,i,nsamples)]*x[XIDX(j,i,nsamples)];
	}

	// Do the iterations
	int it = 0 ;
	double  prev_err,alpha,oldb ;
	double sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {

		it++ ;
		prev_err = *err ;
		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+=y[i]*y[i];

//		fprintf(stderr,"it = %d : Err = %f ,  PrevErr = %f\n",it,*err,prev_err) ;

		for(int j=0; j<nftrs; j++) {
			sumxy = 0 ;

			for(int i=0; i<nsamples; i++)
				sumxy+= x[XIDX(j,i,nsamples)]*y[i];

			alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;
			oldb=b[j];

			b[j] = rfactor * (b[j] + alpha) ;
			alpha = b[j] - oldb;

			for(int i=0; i<nsamples; i++)
				y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
		}
	}	

	free(y) ;
	free(sumxx) ;
	return 0 ;
}

int lm_positive(double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double rfactor, double *b, double *err) {

	// Prepare
	double *y ;
	if ((y = (double *) malloc (nsamples*sizeof (double))) == NULL)
		return -1 ;

	memcpy(y,_y,nsamples*sizeof(double)) ;	
	memset(b,0,nftrs*sizeof(double)) ;

	double *sumxx ;
	if ((sumxx = (double *) malloc (nftrs*sizeof(double)))==NULL)
		return -1 ;
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += x[XIDX(j,i,nsamples)]*x[XIDX(j,i,nsamples)];
	}

	// Do the iterations
	int it = 0 ;
	double  prev_err,alpha,oldb ;
	double sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {

		it++ ;
		prev_err = *err ;
		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+=y[i]*y[i];

//		fprintf(stderr,"it = %d : Err = %f ,  PrevErr = %f\n",it,*err,prev_err) ;

		for(int j=0; j<nftrs; j++) {
			sumxy = 0 ;

			for(int i=0; i<nsamples; i++)
				sumxy+= x[XIDX(j,i,nsamples)]*y[i];

			alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;
			oldb=b[j];

			b[j] =  (b[j] + alpha) ;
			if (b[j]<0) b[j]=0;
			b[j]*=0.99;
			alpha = b[j] - oldb;

			for(int i=0; i<nsamples; i++)
				y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
		}
	}	

	free(y) ;
	free(sumxx) ;
	return 0 ;
}

int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double rfactor, double *ws, double *b, double *err) {

	// Prepare
	double *y ;
	if ((y = (double *) malloc (nsamples*sizeof (double))) == NULL)
		return -1 ;

	memcpy(y,_y,nsamples*sizeof(double)) ;	
	memset(b,0,nftrs*sizeof(double)) ;

	double *sumxx ;
	if ((sumxx = (double *) malloc (nftrs*sizeof(double)))==NULL)
		return -1 ;
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += ws[i] * x[XIDX(j,i,nsamples)]  *x[XIDX(j,i,nsamples)];
	}

	// Do the iterations
	int it = 0 ;
	double  prev_err,alpha,oldb ;
	double sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {

		it++ ;
		prev_err = *err ;
		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+= ws[i] * y[i] * y[i];

//		fprintf(stderr,"it = %d : Err = %f ,  PrevErr = %f\n",it,*err,prev_err) ;

		for(int j=0; j<nftrs; j++) {
			sumxy = 0 ;

			for(int i=0; i<nsamples; i++)
				sumxy+= ws[i] * x[XIDX(j,i,nsamples)]*y[i];

			alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;
			oldb=b[j];

			b[j] = rfactor * (b[j] + alpha) ;
			alpha = b[j] - oldb;

			for(int i=0; i<nsamples; i++)
				y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
		}
	}	

	free(y) ;
	free(sumxx) ;
	return 0 ;
}

int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *ws, double *b, double *err) {

	// Prepare
	double *y ;
	if ((y = (double *) malloc (nsamples*sizeof (double))) == NULL)
		return -1 ;

	memcpy(y,_y,nsamples*sizeof(double)) ;	
	memset(b,0,nftrs*sizeof(double)) ;

	double *sumxx ;
	if ((sumxx = (double *) malloc (nftrs*sizeof(double)))==NULL)
		return -1 ;
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += ws[i] * x[XIDX(j,i,nsamples)]  *x[XIDX(j,i,nsamples)];
	}

	// Do the iterations
	int it = 0 ;
	double  prev_err,alpha,oldb ;
	double sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {

		it++ ;
		prev_err = *err ;
		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+= ws[i] * y[i] * y[i];

//		fprintf(stderr,"it = %d : Err = %f ,  PrevErr = %f\n",it,*err,prev_err) ;

		for(int j=0; j<nftrs; j++) {
			sumxy = 0 ;

			for(int i=0; i<nsamples; i++)
				sumxy+= ws[i] * x[XIDX(j,i,nsamples)]*y[i];

			alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;
			oldb=b[j];

			b[j] = rfactors[j] * (b[j] + alpha) ;
			alpha = b[j] - oldb;

			for(int i=0; i<nsamples; i++)
				y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
		}
	}	

	free(y) ;
	free(sumxx) ;
	return 0 ;
}

int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore) {

	// Prepare
	double *y ;
	if ((y = (double *) malloc (nsamples*sizeof (double))) == NULL)
		return -1 ;

	memcpy(y,_y,nsamples*sizeof(double)) ;	
	memset(b,0,nftrs*sizeof(double)) ;

	double *sumxx ;
	if ((sumxx = (double *) malloc (nftrs*sizeof(double)))==NULL)
		return -1 ;
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += x[XIDX(j,i,nsamples)]*x[XIDX(j,i,nsamples)];
	}

	// Do the iterations
	int it = 0 ;
	double  prev_err,alpha,oldb ;
	double sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {

		it++ ;
		prev_err = *err ;
		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+=y[i]*y[i];

//		fprintf(stderr,"it = %d : Err = %f ,  PrevErr = %f\n",it,*err,prev_err) ;

		for(int j=0; j<nftrs; j++) {
			if (! ignore[j]) {
				sumxy = 0 ;

				for(int i=0; i<nsamples; i++)
					sumxy+= x[XIDX(j,i,nsamples)]*y[i];

				alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;
				oldb=b[j];

				b[j] = rfactors[j] * (b[j] + alpha) ;
				alpha = b[j] - oldb;

				for(int i=0; i<nsamples; i++)
					y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
			}
		}
	}	

	free(y) ;
	free(sumxx) ;
	return 0 ;
}

int lm (double *x, double *_y, double *w, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore, double *corrs) {

	// Prepare
	double *y ;
	if ((y = (double *) malloc (nsamples*sizeof (double))) == NULL)
		return -1 ;

	memcpy(y,_y,nsamples*sizeof(double)) ;	
	memset(b,0,nftrs*sizeof(double)) ;

	double *sumxx ;
	if ((sumxx = (double *) malloc (nftrs*sizeof(double)))==NULL)
		return -1 ;
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += w[i]*x[XIDX(j,i,nsamples)]*x[XIDX(j,i,nsamples)];
	}


	// Do the iterations
	int it = 0 ;
	double  prev_err,alpha,oldb ;
	double sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {

		it++ ;
		prev_err = *err ;
		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+=w[i]*y[i]*y[i];

//		fprintf(stderr,"it = %d : Err = %f ,  PrevErr = %f\n",it,*err,prev_err) ;

		for(int j=0; j<nftrs; j++) {
			if (! ignore[j]) {
				sumxy = 0 ;

				for(int i=0; i<nsamples; i++)
					sumxy+= w[i]*x[XIDX(j,i,nsamples)]*y[i];

				alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;

				if (alpha*corrs[j] < 0)
					alpha = 0 ;

				oldb=b[j];

				b[j] = rfactors[j] * (b[j] + alpha) ;
				alpha = b[j] - oldb;

				for(int i=0; i<nsamples; i++)
					y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
			}
		}
	}	

	free(y) ;
	free(sumxx) ;
	return 0 ;
}

int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore, double *corrs) {

	double *sumxx ;
	if ((sumxx = (double *) malloc (nftrs*sizeof(double)))==NULL)
		return -1 ;
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += x[XIDX(j,i,nsamples)]*x[XIDX(j,i,nsamples)];
	}

	int rc = lm(x,_y,nsamples,nftrs,niter,eiter,rfactors,b,err,ignore,corrs,sumxx) ;

	free(sumxx) ;
	return rc ;

}

int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore, double *corrs, double *sumxx) {

	// Prepare
	double *y ;
	if ((y = (double *) malloc (nsamples*sizeof (double))) == NULL)
		return -1 ;

	memcpy(y,_y,nsamples*sizeof(double)) ;	
	memset(b,0,nftrs*sizeof(double)) ;

	// Do the iterations
	int it = 0 ;
	double  prev_err,alpha,oldb ;
	double sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {

		it++ ;
		prev_err = *err ;
		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+=y[i]*y[i];

//		fprintf(stderr,"it = %d : Err = %f ,  PrevErr = %f\n",it,*err,prev_err) ;

		for(int j=0; j<nftrs; j++) {
			if (! ignore[j]) {
				sumxy = 0 ;

				for(int i=0; i<nsamples; i++)
					sumxy+= x[XIDX(j,i,nsamples)]*y[i];

				alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;

				if (alpha*corrs[j] < 0)
					alpha = 0 ;

				oldb=b[j];

				b[j] = rfactors[j] * (b[j] + alpha) ;
				alpha = b[j] - oldb;

				for(int i=0; i<nsamples; i++)
					y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
			}
		}
	}	

	free(y) ;
	return 0 ;
}

int label_changing_lm (double *x, double *y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore, double *corrs) {

	// Prepare
	memset(b,0,nftrs*sizeof(double)) ;

	double *sumxx ;
	if ((sumxx = (double *) malloc (nftrs*sizeof(double)))==NULL) {
		fprintf(stderr,"Sumxx allocation failed\n") ;
		return -1 ;
	}
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += x[XIDX(j,i,nsamples)]*x[XIDX(j,i,nsamples)];
	}

	// Do the iterations
	int it = 0 ;
	double  prev_err,alpha,oldb ;
	double sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {

		it++ ;
		prev_err = *err ;
		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+=y[i]*y[i];


		for(int j=0; j<nftrs; j++) {
			if (! ignore[j]) {
				sumxy = 0 ;

				for(int i=0; i<nsamples; i++)
					sumxy+= x[XIDX(j,i,nsamples)]*y[i];

				alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;
				if (alpha*corrs[j] < 0)
					alpha = 0 ;

				oldb=b[j];

				b[j] = rfactors[j] * (b[j] + alpha) ;
				alpha = b[j] - oldb;

				for(int i=0; i<nsamples; i++)
					y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
			}
		}
	}	

	free(sumxx) ;
	return 0 ;
}

// Predict
void lm_predict(double *x, double *avg, int nsamples, int nftrs, double *b, double *preds) {

	memset(preds,0,nsamples*sizeof(double)) ;

	for (int j=0; j<nftrs; j++) {
		for (int i=0; i<nsamples; i++) 
			preds[i] += b[j] * (x[XIDX(j,i,nsamples)] - avg[j]) ;
	}

	return ;
}

void lm_predict(double *x, int nsamples, int nftrs, double *b, double *preds) {

	memset(preds,0,nsamples*sizeof(double)) ;

	for (int j=0; j<nftrs; j++) {
		for (int i=0; i<nsamples; i++)
			preds[i] += b[j] * x[XIDX(j,i,nsamples)] ;
	}

	return ;
}

void lm_predict(double *x, int nsamples, int nftrs, double *b, double norm, double *preds) {

	memset(preds,0,nsamples*sizeof(double)) ;

	for (int j=0; j<nftrs; j++) {
		for (int i=0; i<nsamples; i++)
			preds[i] += b[j] * x[XIDX(j,i,nsamples)] ;
	}

	for (int i=0; i<nsamples; i++)
		preds[i] += norm ;

	return ;
}


// Complete features in (non-transposed) matrix according to linear regression params
void complete_ftrs(double *x, int *missing, int nsamples, int nftrs, double *complete_x, double **ftr_bs) {

	double value ;
	for (int i=0; i<nsamples; i++) {
		for (int j=0; j<nftrs; j++) {

			if (missing[XIDX(i,j,nftrs)]) {
				value = 0 ;
				for (int k=0; k<nftrs; k++)
					value += ftr_bs[j][k] * x[XIDX(i,k,nftrs)] ;
			} else
				value = x[XIDX(i,j,nftrs)] ;

			complete_x[XIDX(i,j,nftrs)] = value ;
		}
	}

	return ;
}

// Linear regression to predict missing features
int features_lm(double *x, int nsamples , int nftrs,int niter ,double eiter, double *rfactors, double **ftr_bs) {

	int *ignore ;
	double *ftr_labels ;
	if ((ignore = (int *) malloc (nftrs*sizeof(int)))==NULL || (ftr_labels = (double *) malloc (nsamples*sizeof(double)))==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}
	
	double err ;

	for (int i=0; i<nftrs; i++) {
		if ((ftr_bs[i] = (double *) malloc (nftrs*sizeof(double)))==NULL) {
			fprintf(stderr,"Allocation %d failed\n",i) ;
			return -1 ;
		}

		memset(ignore,0,nftrs*sizeof(int)) ;
		ignore[i] = 1 ;

		for (int j=0; j<nsamples; j++)
			ftr_labels[j] = x[XIDX(i,j,nsamples)] ;

		if (lm(x,ftr_labels,nsamples,nftrs,niter,eiter,rfactors,ftr_bs[i],&err,ignore)==-1) {
			fprintf(stderr,"lm failed\n") ;
			return -1 ;
		}

		fprintf(stderr,"Prediction of label %d : Err = %.3f\n",i,err) ;
	}

	return 0 ;
}
