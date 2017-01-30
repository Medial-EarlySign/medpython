// KNN functions

#include "data_handler/data_handler/data_handler.h"
#include "classifiers/classifiers/classifiers.h"
#include "medial_utilities/medial_utilities/globalRNG.h"

#define INF_DIST  99999999.0
#define MIN_DIST 1e-5
#define EPSILON  1e-5
#define NCHECK -1
#define NRAND 500

// Find k nearest neighbors and their squared-distances. test2learn is the index of the sample in question which is ignored.
void find_nbrs(double *test_x, int ind, int k, double *learn_x, int nlearn, int nftrs, int *order, double *ws, int *nbrs, double *dists, int norm, int test2learn) {

	// Initialize with INF
	for (int i=0; i<k; i++) {
		dists[i] = INF_DIST ;
		nbrs[i] = -1 ;
	  }
	
	double target_dist = INF_DIST ;
	int target_nbr = 0 ;

	// Go over learn-table and find nearest neighbours
	for (int i=0; i<nlearn; i++) {

		if (i == test2learn)
			continue ;

		double dist = 0 ;
		for (int j=0; j<nftrs; j++) {
			if (ws[order[j]] == 0)
				break ;

			double d = ws[order[j]] * (learn_x[XIDX(i,order[j],nftrs)] - test_x[XIDX(ind,order[j],nftrs)]) ;
			if (norm == 2)
				dist += d*d ;
			else
				dist += fabs(d) ;

			if (dist >= target_dist)
				break ;
		}

		// Replace
		if (dist < target_dist) {
			nbrs[target_nbr] = i ;
			target_dist = dists[target_nbr] = dist ;

			// Find new target (most distance current neighbor)
			for (int j=0; j<k; j++) {
				if (dists[j] > target_dist) {
					target_dist = dists[j] ;
					target_nbr = j ;
				}
			}
		}
	}
}


// Find mean distance by random sampling
double get_mean_dist(double *test_x, int ind, double *learn_x, int nlearn, int nftrs, int nrand, double *ws, int norm) {
				   
	double sum = 0 ;

	for (int i=0; i<nrand; i++) {
		int irand = (int) (nlearn * globalRNG::rand()/(globalRNG::max() + 1.0)) ;

		for (int j=0; j<nftrs; j++) {
			double d = ws[j] * (learn_x[XIDX(irand,j,nftrs)] - test_x[XIDX(ind,j,nftrs)]) ;
			if (norm == 2)
				sum += d*d ;
			else 
				sum += fabs(d) ;
		}
	}

	return sum/nrand ;
}

// Get score from neihgbors
double nbrs_score(int *nbrs, double *dists, double *y, int k, double mean_dist, int type) {

	double pred = 0 ;

	if (type == 1) {
		double sumw = 0 ;
		for (int i=0; i<k; i++) {
			if (dists[i] != -1) {
				double w = (dists[i] > mean_dist) ? 0 : ( 1 - sqrt(dists[i]/mean_dist) ) ;
				pred += y[nbrs[i]] * w ;
				sumw += w ;
			}
		}
		pred /= sumw ;
	} else if (type == 2) {
		double sumw = 0 ;
		for (int i=0; i<k; i++) {
			if (nbrs[i] != -1) {
				if (dists[i] < 1e-5)
					dists[i] = 1e-5 ;
				double w = 1/(dists[i]) ;
				pred += y[nbrs[i]] * w ;
				sumw += w ;
			}
		}
		pred /= sumw ;
	} 

	return pred ;
}

// Get Score from weighted linear model on neighbours
int nbrs_ls(double *testx, int ind, int *nbrs, double *dists, double *x, double *y, int k, int nftrs, double mean_dist, double *localx, double *localy, double *localw,
			double *localb, double *localr, double *pred) {

	// Create
	int nrows = 0 ;
	for (int i=0; i<k; i++) {
		if (nbrs[i] != -1)
			nrows++ ;
	}

	int irow = 0 ;
	for (int i=0; i<k; i++) {
		if (nbrs[i] != -1) {
			localw[irow] = (dists[i] > mean_dist) ? 0 : ( 1 - sqrt(dists[i]/mean_dist) ) ;
			localy[irow] = y[nbrs[i]] ;

			for (int j=0; j<nftrs; j++)
				localx[XIDX(j,irow,nrows)] = x[XIDX(nbrs[i],j,nftrs)] ;
			irow ++ ;
		}
	}

	if (nrows <= nftrs) {
		fprintf(stderr,"Couldn't find enough neihbors for nbrs-ls (found %d with %d features)\n",nrows,nftrs) ;
		return -1 ;
	}

	// Normalize
	double yavg,*xavg,*xstd ;
	if (tcalc_stats(localx,localy,localw,nrows,nftrs,&xavg,&xstd,&yavg) == -1)
		return -1 ;

	tnormalize_data(localx,localy,nrows,nftrs,xavg,yavg) ;

	// Learn
	for (int i=0; i<nftrs; i++) {
		double corr = tget_corr(localx,localy,nrows,i) ;
		localr[i] = sqrt(fabs(corr)) ;
	}

	double err ;
	if (lm(localx,localy,nrows,nftrs,NITER,EITER,localr,localw,localb,&err)==-1)
		return -1 ;

	// Predict
	*pred = yavg ;
	for (int j=0; j<nftrs; j++)
		*pred += localb[j] * (testx[XIDX(ind,j,nftrs)] - xavg[j]) ;

	free(xavg) ;
	free(xstd) ;

	return 0 ;
}

// Get orders of features by stdv x weight
int order_ftrs(double *x, int nsamples, int nftrs, double *weights, int *order) {

	for (int i=0; i<nftrs; i++)
		order[i] = i ;


	double *avg,*std ;
	if ((avg = (double *) malloc (nftrs*sizeof(double)))==NULL || (std = (double *) malloc (nftrs*sizeof(double)))==NULL) {
		fprintf(stderr,"Stats allocation failed\n") ;
		return -1 ;
	}

	calc_xstats(x,nsamples,nftrs,avg,std) ;

	// sort by weights x std
	for (int i=0; i<nftrs-1; i++) {
		for (int j=i+1; j<nftrs; j++) {
			if (weights[order[i]]*std[order[i]] < weights[order[j]]*std[order[j]]) {
				int temp = order[i] ;
				order[i] = order[j] ;
				order[j] = temp ;
			}
		}
	}

	free(avg) ;
	free(std) ;

	return 0 ;
}

// Predict a single instance using KNN
int knn_predict (double *test_x, int ind, double *learn_x, double *learn_y, int nlearn, int nftrs, double *weights, int *order, int *nbrs, double *dists, int k,
					int type, double *nbrs_x, double *nbrs_y, double *nbrs_w, double *nbrs_b, double *nbrs_r, double *pred) {

	int norm = 2 ;
	if (type == 4) {
		type = 3 ;
		norm = 1 ;
	}

	find_nbrs(test_x,ind,k,learn_x,nlearn,nftrs,order,weights,nbrs,dists,norm,-1) ;
	double mean = get_mean_dist(test_x,ind,learn_x,nlearn,nftrs,NRAND,weights,norm) ;

	if (type == 1 || type == 2) {
		*pred = nbrs_score(nbrs,dists,learn_y,k,mean,type) ;
		return 0 ;
	} else if (type == 3) {
		return nbrs_ls(test_x,ind,nbrs,dists,learn_x,learn_y,k,nftrs,mean,nbrs_x,nbrs_y,nbrs_w,nbrs_b,nbrs_r,pred) ;
	} else
		return -1 ;
}

// Predict a single instance using KNN. test2learn is the index of the sample in question which is ignored.
int knn_predict (double *test_x, int ind, double *learn_x, double *learn_y, int nlearn, int nftrs, double *weights, int *order, int *nbrs, double *dists, int k,
					int type, double *nbrs_x, double *nbrs_y, double *nbrs_w, double *nbrs_b, double *nbrs_r, int test2learn, double *pred) {

	int norm = 2 ;
	if (type == 4) {
		type = 3 ;
		norm = 1 ;
	}

	find_nbrs(test_x,ind,k,learn_x,nlearn,nftrs,order,weights,nbrs,dists,norm,test2learn) ;
	double mean = get_mean_dist(test_x,ind,learn_x,nlearn,nftrs,NRAND,weights,norm) ;

	if (type == 1 || type == 2) {
		*pred = nbrs_score(nbrs,dists,learn_y,k,mean,type) ;
		return 0 ;
	} else if (type == 3) {
		return nbrs_ls(test_x,ind,nbrs,dists,learn_x,learn_y,k,nftrs,mean,nbrs_x,nbrs_y,nbrs_w,nbrs_b,nbrs_r,pred) ;
	} else
		return -1 ;
}


// Predict using KNN :
// Type 1 = KNN using Eucleadean distance, weighing each neighbor by dist/mean
// Type 2 = KNN using Eucleadean distance, weighing each neighbor by 1/dist
// Type 3 == KNN using Eucleadean distance + Weighted LS
// Type 4 == KNN Using L1 distance + Weighted LS

int knn_predict(double *test_x, int ntest, double *learn_x, double *learn_y, int nlearn, int nftrs, double *weights, double *preds, int k, int type) {

	// Checks ...
	if (type != 1 && type != 2 && type != 3 && type != 4) {
		fprintf(stderr,"Unknown knn type %d\n",type) ;
		return -1 ;
	}

	if ((type == 3 || type == 4) && k < nftrs) {
		fprintf(stderr,"k (%d) must be larger than nftrs (%d) in KNN+LS\n",k,nftrs) ;
		return -1 ;
	}

	// OK, lets go ...
	fprintf(stderr,"Running knn : K = %d , Data = (%d + %d) x %d\n",k,ntest,nlearn,nftrs) ;

	// Allocation
	int *order,*nbrs ;
	double *dists ;

	if ((order = (int *) malloc (nftrs*sizeof(int)))==NULL || (nbrs = (int *) malloc (k*sizeof(int)))==NULL || 
		(dists = (double *) malloc (k*sizeof(double)))==NULL) {
			fprintf(stderr,"Allocation failed\n") ;
			return -1 ;
	}

	double *nbrs_x = NULL, *nbrs_y = NULL, *nbrs_w = NULL, *nbrs_b = NULL, *nbrs_r = NULL;
	if (type == 3 || type == 4) {
		if ((nbrs_x = (double *) malloc (k*nftrs*sizeof(double)))==NULL || (nbrs_y = (double *) malloc(k*sizeof(double)))==NULL ||
			(nbrs_w = (double *) malloc (k*sizeof(double)))==NULL || (nbrs_b = (double *) malloc(nftrs*sizeof(double)))==NULL || 
			(nbrs_r = (double *) malloc (nftrs*sizeof(double)))==NULL) {
				fprintf(stderr,"nbrs data allocation failed\n") ;
				return -1 ;
		}
	}

	// Order features
	if (order_ftrs(learn_x,nlearn,nftrs,weights,order)==-1)
		return -1 ;

	for (int i=0; i<ntest; i++) {
		if (i%1000 == 1)
			fprintf(stderr,"Predicting %d/%d\n",i,ntest) ;
		if (knn_predict(test_x,i,learn_x,learn_y,nlearn,nftrs,weights,order,nbrs,dists,k,type,nbrs_x,nbrs_y,nbrs_w,nbrs_b,nbrs_r,&(preds[i]))==-1) {
			fprintf(stderr,"knn prediction failed\n") ;
			return -1 ;
		}
	}

	free(order) ;
	free(dists) ;
	free(nbrs) ;

	if (type == 3 || type ==4) {
		free(nbrs_x) ;
		free(nbrs_y) ;
		free(nbrs_w) ;
		free(nbrs_b) ;
		free(nbrs_r) ;
	}

	return 0 ;
}

// Predict using KNN :
// Type 1 = KNN using Eucleadean distance, weighing each neighbor by dist/mean
// Type 2 = KNN using Eucleadean distance, weighing each neighbor by 1/dist
// Type 3 == KNN using Eucleadean distance + Weighted LS
// Type 4 == KNN Using L1 distance + Weighted LS

int knn_predict(double *test_x, int ntest, double *learn_x, double *learn_y, int nlearn, int nftrs, double *weights, double *preds, int k, int *finds, int type) {

	// Checks ...
	if (type != 1 && type != 2 && type != 3 && type != 4) {
		fprintf(stderr,"Unknown knn type %d\n",type) ;
		return -1 ;
	}

	if ((type == 3 || type == 4) && k < nftrs) {
		fprintf(stderr,"k (%d) must be larger than nftrs (%d) in KNN+LS\n",k,nftrs) ;
		return -1 ;
	}

	// OK, lets go ...
	fprintf(stderr,"Running knn : K = %d , Data = (%d + %d) x %d\n",k,ntest,nlearn,nftrs) ;

	// Allocation
	int *order,*nbrs ;
	double *dists ;

	if ((order = (int *) malloc (nftrs*sizeof(int)))==NULL || (nbrs = (int *) malloc (k*sizeof(int)))==NULL || 
		(dists = (double *) malloc (k*sizeof(double)))==NULL) {
			fprintf(stderr,"Allocation failed\n") ;
			return -1 ;
	}

	double *nbrs_x = NULL, *nbrs_y = NULL, *nbrs_w = NULL, *nbrs_b = NULL, *nbrs_r = NULL;
	if (type == 3 || type == 4) {
		if ((nbrs_x = (double *) malloc (k*nftrs*sizeof(double)))==NULL || (nbrs_y = (double *) malloc(k*sizeof(double)))==NULL ||
			(nbrs_w = (double *) malloc (k*sizeof(double)))==NULL || (nbrs_b = (double *) malloc(nftrs*sizeof(double)))==NULL || 
			(nbrs_r = (double *) malloc (nftrs*sizeof(double)))==NULL) {
				fprintf(stderr,"nbrs data allocation failed\n") ;
				return -1 ;
		}
	}

	// Order features
	if (order_ftrs(learn_x,nlearn,nftrs,weights,order)==-1)
		return -1 ;

	for (int i=0; i<ntest; i++) {
		if (i%1000 == 1)
			fprintf(stderr,"Predicting %d/%d\n",i,ntest) ;
		if (knn_predict(test_x,i,learn_x,learn_y,nlearn,nftrs,weights,order,nbrs,dists,k,type,nbrs_x,nbrs_y,nbrs_w,nbrs_b,nbrs_r,finds[i],&(preds[i]))==-1) {
			fprintf(stderr,"knn prediction failed\n") ;
			return -1 ;
		}
	}

	free(order) ;
	free(dists) ;
	free(nbrs) ;

	return 0 ;
}