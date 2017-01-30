#include "data_handler.h"


#define DEF_NRUNS 1
#define NCLN_ITER 15
#define MAX_STD 15
#define RESOLUTION 100
#define RATIO 3

// Split (Test+Learn)
int tsplit(double *xtable, double *labels, int nsamples, int nftrs, double *xtable1, double *labels1, double *xtable2, double *labels2, int nlearn) {

	// Random ordering
	int *order ;
	if ((order = randomize(nsamples))==NULL)
		return -1 ;

	// Learn
	for (int i=0; i<nlearn; i++) {
		labels1[i] = labels[order[i]] ;
		for (int j=0; j<nftrs; j++)
			xtable1[XIDX(j,i,nlearn)] = xtable[XIDX(j,order[i],nsamples)] ;
	}

	// Test
	int ntest = nsamples-nlearn ;
	for (int i=0; i<ntest; i++) {
		labels2[i] = labels[order[ntest+i]] ;
		for (int j=0; j<nftrs; j++)
			xtable2[XIDX(j,i,ntest)] = xtable[XIDX(j,order[ntest+i],nsamples)] ;
	}


	free(order) ;
	return 0 ;
}

int tsplit(double *xtable, double *labels, int nsamples, int nftrs, double *xtable1, double *labels1, double *xtable2, double *labels2, int *order, int test_start, 
		   int ntest) {

	int nlearn = nsamples - ntest ;
	
	// Learn
	for (int i=0; i<test_start; i++) {
		labels1[i] = labels[order[i]] ;
		for (int j=0; j<nftrs; j++) 
			xtable1[XIDX(j,i,nlearn)] = xtable[XIDX(j,order[i],nsamples)] ;
	}

	// Test
	for (int i=0; i<ntest; i++) {
		labels2[i] = labels[order[test_start+i]] ;
		for (int j=0; j<nftrs; j++)
			xtable2[XIDX(j,i,ntest)] = xtable[XIDX(j,order[test_start+i],nsamples)] ;
	}

	// Learn
	for (int i=0; i<(nsamples-test_start-ntest);  i++) {
		labels1[test_start+i] = labels[order[test_start+ntest+i]] ;
		for (int j=0; j<nftrs; j++)
			xtable1[XIDX(j,test_start+i,nlearn)] = xtable[XIDX(j,order[test_start+ntest+i],nsamples)] ;
	}

	return 0 ;
}


// Split, keeping all samples from a each patient in one of the groups.
int pid_tsplit(double *x, double *y, char *ids, int nrows, int ncols, int *order, double learn_ratio, double *learn_x, double *learn_y, double *test_x, double *test_y,
			  int *nlearn, char *test_ids) {

	// count
	int pid,prev_pid = -1 ;

	int *pid_cnts = (int *) malloc(nrows*sizeof(int)) ;

	if (pid_cnts==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	int npids=0 ;
	for (int i=0; i<nrows; i++) {
		if (sscanf(ids+HIDX(i),"%d",&pid)!= 1) {
			fprintf(stderr,"Cannot parse id %s\n",ids+HIDX(i)) ;
			return -1 ;
		}

		if (pid != prev_pid) {
			npids++ ;
			pid_cnts[npids-1] = 1 ;
		} else
			pid_cnts[npids-1] ++ ;

		prev_pid = pid ;
	}

	// Mark
	int *mark = (int *) malloc (npids*sizeof(int)) ;
	if (mark==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	*nlearn = 0 ;
	memset(mark,0,npids*sizeof(int)) ;
	for (int i=0; i< ((int) (learn_ratio*npids)); i++) {
		mark[order[i]] = 1 ;
		(*nlearn) += pid_cnts[order[i]] ;
	}

	int ntest = nrows - (*nlearn) ;

	// Another pass
	int ilearn = 0 ;
	int itest = 0 ;
	int ipid = -1 ;
	prev_pid = -1 ;

	for (int i=0; i<nrows; i++) {
		sscanf(ids+HIDX(i),"%d",&pid) ;
		if (pid != prev_pid)
			ipid++ ;

		if (mark[ipid] == 1) {
			for (int j=0; j<ncols; j++)
				learn_x[XIDX(j,ilearn,*nlearn)] = x[XIDX(j,i,nrows)] ;
			learn_y[ilearn++] = y[i] ;
		} else {
			for (int j=0; j<ncols; j++)
				test_x[XIDX(j,itest,ntest)] = x[XIDX(j,i,nrows)] ;
			strcpy(test_ids+HIDX(itest),ids+HIDX(i)) ;
			test_y[itest++] = y[i] ;
		}

		prev_pid = pid ;
	}

	// Clean
	free(mark) ;
	free(pid_cnts) ;
	free(order) ;

	return 0 ;
}

// Get Flags for splitting
int get_flags(char *ids, int nrows, double lratio, int *flags, int *nlearn) {

	// count
	int pid,prev_pid = -1 ;

	int *pid_cnts = (int *) malloc(nrows*sizeof(int)) ;

	if (pid_cnts==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	int npids=0 ;
	for (int i=0; i<nrows; i++) {
		if (sscanf(ids+HIDX(i),"%d",&pid)!= 1) {
			fprintf(stderr,"Cannot parse id %s\n",ids+HIDX(i)) ;
			return -1 ;
		}

		if (pid != prev_pid) {
			npids++ ;
			pid_cnts[npids-1] = 1 ;
		} else
			pid_cnts[npids-1] ++ ;

		prev_pid = pid ;
	}

	// Randomize ...
	int *order = randomize(npids) ;
	if (order == NULL) {
		fprintf(stderr,"randomization oreder failed\n") ;
		return -1 ;
	}

	// Mark
	int *mark = (int *) malloc (npids*sizeof(int)) ;
	if (mark==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	*nlearn = 0 ;
	memset(mark,0,npids*sizeof(int)) ;
	for (int i=0; i< ((int) (lratio*npids)); i++) {
		mark[order[i]] = 1 ;
		(*nlearn) += pid_cnts[order[i]] ;
	}

	// Another pass for flags
	int ipid = -1 ;
	prev_pid = -1 ;

	for (int i=0; i<nrows; i++) {
		sscanf(ids+HIDX(i),"%d",&pid) ;
		if (pid != prev_pid)
			ipid++ ;

		flags[i] = 1 - mark[ipid] ;
		prev_pid = pid ;
	}

	// Clean
	free(mark) ;
	free(pid_cnts) ;
	free(order) ;

	return 0 ;
}

int get_flags(char *ids, int nrows, int *order, double lratio, int *flags, int *nlearn) {

	// count
	int pid,prev_pid = -1 ;

	int *pid_cnts = (int *) malloc(nrows*sizeof(int)) ;

	if (pid_cnts==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	int npids=0 ;
	for (int i=0; i<nrows; i++) {
		if (sscanf(ids+HIDX(i),"%d",&pid)!= 1) {
			fprintf(stderr,"Cannot parse id %s\n",ids+HIDX(i)) ;
			return -1 ;
		}

		if (pid != prev_pid) {
			npids++ ;
			pid_cnts[npids-1] = 1 ;
		} else
			pid_cnts[npids-1] ++ ;

		prev_pid = pid ;
	}

	// Mark
	int *mark = (int *) malloc (npids*sizeof(int)) ;
	if (mark==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1 ;
	}

	*nlearn = 0 ;
	memset(mark,0,npids*sizeof(int)) ;
	for (int i=0; i< ((int) (lratio*npids)); i++) {
		mark[order[i]] = 1 ;
		(*nlearn) += pid_cnts[order[i]] ;
	}

	// Another pass for flags
	int ipid = -1 ;
	prev_pid = -1 ;

	for (int i=0; i<nrows; i++) {
		sscanf(ids+HIDX(i),"%d",&pid) ;
		if (pid != prev_pid)
			ipid++ ;

		flags[i] = 1 - mark[ipid] ;
		prev_pid = pid ;
	}

	// Clean
	free(mark) ;
	free(pid_cnts) ;

	return 0 ;
}

// Split (Test+Learn) - untransposed data, with missing flags
int split(double *xtable, int *missing, double *labels, int nsamples, int nftrs, double *xtable1, int *missing1, double *labels1, 
		  double *xtable2, int *missing2, double *labels2, int nlearn) {

	// Random ordering
	int *order ;
	if ((order = randomize(nsamples))==NULL)
		return -1 ;

	// Learn
	for (int i=0; i<nlearn; i++) {
		labels1[i] = labels[order[i]] ;
		for (int j=0; j<nftrs; j++)  {
			xtable1[XIDX(i,j,nftrs)] = xtable[XIDX(order[i],j,nftrs)] ;
			missing1[XIDX(i,j,nftrs)] = missing[XIDX(order[i],j,nftrs)] ;
		}
	}

	// Test
	for (int i=nlearn; i<nsamples; i++) {
		labels2[i-nlearn] = labels[order[i]] ;
		for (int j=0; j<nftrs; j++) {
			xtable2[XIDX(i-nlearn,j,nftrs)] = xtable[XIDX(order[i],j,nftrs)] ;	
			missing1[XIDX(i-nlearn,j,nftrs)] = missing[XIDX(order[i],j,nftrs)] ;
		}
	}

	return 0 ;
}


// Handle missing values - make all true values positive and set missing values to -1.
void adjust_matrix (double *xtable, int *missing, int nsamples, int nftrs) {

	for (int i=0; i<nftrs; i++) {
		double min = 0 ;
		for (int j=0; j<nsamples; j++) {
			if (!missing[XIDX(j,i,nftrs)] && xtable[XIDX(j,i,nftrs)]<min)
				min = xtable[XIDX(j,i,nftrs)] ;
		}

		for (int j=0; j<nsamples; j++) {
			if (missing[XIDX(j,i,nftrs)])
				xtable[XIDX(j,i,nftrs)] = -1.0 ;
			else
				xtable[XIDX(j,i,nftrs)] -= min ;
		}
	}

	return ;
}

// Calculate statistics of xtable 
void tcalc_xstats(double *x, int nsamples, int nftrs, double *avg, double *std, double missing)
{

	memset(avg,0,nftrs*sizeof(double)) ;
	memset(std,0,nftrs*sizeof(double)) ;

	for (int j=0; j<nftrs; j++) {
		double sum = 0 , sum2 = 0 ;
		int cnt = 0 ;
		for (int i=0; i<nsamples; i++) {
			double val = x[XIDX(j,i,nsamples)] ;

			if (val != missing) {
				sum += val ;
				sum2 += val*val ;
				cnt++ ;
			}
		}

		if (cnt != 0) {
			avg[j] = sum/cnt ;
			std[j] = sqrt ((sum2 - avg[j]*sum)/(cnt-1)) ;
		} else {
			avg[j] = 0.0 ;
			std[j] = 1.0 ;
		}
	}

	return ;
}

void tcalc_xstats(double *x, double *w, int nsamples, int nftrs, double *avg, double *std, double missing)
{

	memset(avg,0,nftrs*sizeof(double)) ;
	memset(std,0,nftrs*sizeof(double)) ;

	for (int j=0; j<nftrs; j++) {
		double sum = 0 ;
		double norm = 0 ;
		for (int i=0; i<nsamples; i++) {
			double val = x[XIDX(j,i,nsamples)] ;
			double weight = w[i] ;

			if (val != missing) {
				sum += val*weight ;
				norm += weight ;
			}
		}

		if (norm != 0) {
			avg[j] = sum/norm ;

			sum = 0 ;
			for (int i=0; i<nsamples; i++) {
				double val = x[XIDX(j,i,nsamples)] ;
				double weight = w[i] ;

				if (val != missing)
					sum += weight*(val - avg[j])*(val - avg[j]) ;
			}

			std[j] = sqrt (sum/norm) ;
			if (std[j] == 0)
				std[j] = 1.0 ;

		} else {
			avg[j] = 0.0 ;
			std[j] = 1.0 ;
		}
	}

	return ;
}

// Calculate statistics of xtable and ytable
int tcalc_stats(double *x, double *y, int nsamples, int nftrs, double **avg, double **std, double *yavg, double missing)
{
	if (((*avg) = (double *) malloc(nftrs*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate averages for %d\n", nftrs) ;
		return -1 ;
	}

	if (((*std) = (double *) malloc(nftrs*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate stds for %d\n", nftrs) ;
		free(*avg) ;
		return -1 ;
	}

	tcalc_xstats(x,nsamples,nftrs,*avg,*std, missing);

	double sum=0; 
	for (int i=0; i<nsamples; i++)
		sum += y[i] ;
	(*yavg) = sum/nsamples ;

	return 0;
}

int tcalc_stats(double *x, double *y, double *w, int nsamples, int nftrs, double **avg, double **std, double *yavg, double missing)
{
	if (((*avg) = (double *) malloc(nftrs*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate averages for %d\n", nftrs) ;
		return -1 ;
	}

	if (((*std) = (double *) malloc(nftrs*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate stds for %d\n", nftrs) ;
		free(*avg) ;
		return -1 ;
	}

	tcalc_xstats(x,w,nsamples,nftrs,*avg,*std,missing);

	double sum=0; 
	double norm=0 ;
	for (int i=0; i<nsamples; i++) {
		sum += y[i]*w[i] ;
		norm += w[i] ;
	}
	(*yavg) = sum/norm ;

	return 0;
}

// Calculate statistics of central part of each vector
int tcalc_partial_stats(double *xtable, int npatient, int nvar, double central_p, double *avg, double *std, double missing) {

	double *vec = (double *) malloc(npatient*sizeof(double)) ;
	if (vec == NULL) {
		fprintf(stderr,"Allocation Failed\n") ;
		return -1 ;
	}

	for (int i=0; i<nvar; i++) {
		int n=0 ;
		for (int j=0; j<npatient; j++) {
			if (xtable[XIDX(i,j,npatient)] != missing)
				vec[n++] = xtable[XIDX(i,j,npatient)] ;
		}

		qsort(vec,n,sizeof(double),double_compare) ;

		int start = (int) (n*(1-central_p)/2) ;
		int neff = (int) (n*central_p) ;

		if (get_moments(vec + start,neff,avg+i,std+i,missing) == -1) {
			avg[i] = -1.0 ;
			std[i] = -1.0 ;
		}
	}

	return 0 ;
}

// Data Normalization - mean + standard deviation
void tnormalize_data(double *x, double *y, int nsamples, int nftrs, double *avg, double *std, double yavg, double missing)
{
	// Reduce average from each column 

	for(int i=0; i<nsamples; i++) {
		for(int j=0; j<nftrs; j++) {
			if (x[XIDX(j,i,nsamples)]==missing)
				x[XIDX(j,i,nsamples)] = 0;		
			else
				x[XIDX(j,i,nsamples)]= (x[XIDX(j,i,nsamples)] - avg[j])/std[j];
		}

		y[i]-=yavg;
	}
}

// Normalize data - set means to zero.
void tnormalize_data(double *x, double *y, int nsamples, int nftrs, double *avg, double yavg, double missing)
{
	// Reduce average from each column 


	for(int i=0; i<nsamples; i++) {
		for(int j=0; j<nftrs; j++) {
			if (x[XIDX(j,i,nsamples)]==missing)
				x[XIDX(j,i,nsamples)] = 0;		
			else
				x[XIDX(j,i,nsamples)]-=avg[j];
		}

		y[i]-=yavg;
	}
}

void tnormalize_data(double *x, int nsamples, int nftrs, double *avg, double missing)
{
	// Reduce average from each column 


	for(int i=0; i<nsamples; i++) {
		for(int j=0; j<nftrs; j++) {
			if (x[XIDX(j,i,nsamples)]==missing)
				x[XIDX(j,i,nsamples)] = 0;		
			else
				x[XIDX(j,i,nsamples)]-=avg[j];
		}
	}
}

void tnormalize_data(double *x, int nsamples, int nftrs, double *avg, double *std, double missing)
{
	// Reduce average from each column 


	for(int i=0; i<nsamples; i++) {
		for(int j=0; j<nftrs; j++) {
			if (x[XIDX(j,i,nsamples)]==missing)
				x[XIDX(j,i,nsamples)] = 0;		
			else
				x[XIDX(j,i,nsamples)] = (x[XIDX(j,i,nsamples)] - avg[j])/std[j] ;
		}
	}
}

void tnormalize_data_keeping_missing(double *x, double *y, int nsamples, int nftrs, double *avg, double yavg, double missing)
{
	// Reduce average from each column 


	for(int i=0; i<nsamples; i++) {
		for(int j=0; j<nftrs; j++) {
			if (x[XIDX(j,i,nsamples)]!=missing)
				x[XIDX(j,i,nsamples)]-=avg[j];
		}

		y[i]-=yavg;
	}
}

void tnormalize_data_keeping_missing(double *x, int nsamples, int nftrs, double *avg, double missing)
{
	// Reduce average from each column 


	for(int i=0; i<nsamples; i++) {
		for(int j=0; j<nftrs; j++) {
			if (x[XIDX(j,i,nsamples)]!=missing)
				x[XIDX(j,i,nsamples)]-=avg[j];
		}
	}
}



// Transpose matrix (for linear regression efficieny)
int transpose(double *tx, double **x, int nsamples, int nftrs, int reverse) {

	if (((*x) = (double *) malloc (nsamples * nftrs * sizeof(double)))==NULL)
		return -1 ;

	for (int i=0; i<nsamples; i++) {
		for (int j=0; j<nftrs; j++) {
			if (reverse)
				(*x)[XIDX(i,j,nftrs)] = tx[XIDX(j,i,nsamples)] ;
			else
				(*x)[XIDX(j,i,nsamples)] = tx[XIDX(i,j,nftrs)] ;
		}
	}

	return 0 ;
}

int transpose(int *tx, int **x, int nsamples, int nftrs, int reverse) {

	if (((*x) = (int *) malloc (nsamples * nftrs * sizeof(int)))==NULL)
		return -1 ;

	for (int i=0; i<nsamples; i++) {
		for (int j=0; j<nftrs; j++) {
			if (reverse)
				(*x)[XIDX(i,j,nftrs)] = tx[XIDX(j,i,nsamples)] ;
			else
				(*x)[XIDX(j,i,nsamples)] = tx[XIDX(i,j,nftrs)] ;
		}
	}

	return 0 ;
}

// Normalize x1,y1 and x2 according to means and standard-deviation of x1,y1
int normalize_inout(double *x1, double *y1, double *w1, int nrows1, double *x2, int nrows2, int ncols, double missing_val) {

	double *xavg,*xstd,yavg ;
	if (tcalc_stats(x1,y1,w1,nrows1,ncols,&xavg,&xstd,&yavg,missing_val)==-1) {
		fprintf(stderr,"TCalc_stats failed\n") ;
		return -1 ;
	}

	tnormalize_data(x1,y1,nrows1,ncols,xavg,xstd,yavg,missing_val) ;
	tnormalize_data(x2,nrows2,ncols,xavg,xstd,missing_val) ;
	free(xavg) ; free(xstd) ;

	return 0 ;
}

int normalize_inout(double *x1, double *y1, int nrows1, double *x2, int nrows2, int ncols,  double missing_val) {

	double *xavg,*xstd,yavg ;
	if (tcalc_stats(x1,y1,nrows1,ncols,&xavg,&xstd,&yavg,missing_val)==-1) {
		fprintf(stderr,"TCalc_stats failed\n") ;
		return -1 ;
	}

	tnormalize_data(x1,y1,nrows1,ncols,xavg,xstd,yavg,missing_val) ;
	tnormalize_data(x2,nrows2,ncols,xavg,xstd,missing_val) ;
	free(xavg) ; free(xstd) ;

	return 0 ;
}

// Correlation to label
double tget_corr(double *x, double *y, int nsamples, int ind, double missing) {

	double *vec1,*vec2 ;
	if ((vec1 = (double *) malloc (nsamples*sizeof(double)))==NULL || (vec2 = (double *) malloc (nsamples*sizeof(double)))==NULL) {
		fprintf(stderr,"Allocation of vector failed\n") ;
		return -1 ;
	}

	int n=0 ;
	for (int i=0; i<nsamples; i++) {
		if (x[XIDX(ind,i,nsamples)] != missing) {
			vec1[n] = x[XIDX(ind,i,nsamples)] ;
			vec2[n++] = y[i] ;
		}
	}

	double r2 = pearson(vec1,vec2,n) ;

	free(vec1) ;
	free(vec2) ;

	return r2 ;
}

double get_corr(double *x, double *y, int nsamples, int ind, int nftrs, double missing) {

	double *vec1,*vec2 ;
	if ((vec1 = (double *) malloc (nsamples*sizeof(double)))==NULL || (vec2 = (double *) malloc (nsamples*sizeof(double)))==NULL) {
		fprintf(stderr,"Allocation of vector failed\n") ;
		return -1 ;
	}

	int n=0 ;
	for (int i=0; i<nsamples; i++) {
		if (x[XIDX(i,ind,nftrs)] != missing) {
			vec1[n] = x[XIDX(i,ind,nftrs)] ;
			vec2[n++] = y[i] ;
		}
	}

	double r2 = pearson(vec1,vec2,n) ;

	free(vec1) ;
	free(vec2) ;

	return r2 ;
}

// Filter to a given neg/pos ratio
int tfilter(double *x, double *y, int nsamples, int nftrs, double ratio, double *fx, double *fy, int *fnsamples) {

	// Do nothing if ratio <= 0
	if (ratio <= 0) {
		*fnsamples = nsamples ;
		memcpy(fx,x,nsamples*nftrs*sizeof(double)) ;
		memcpy(fy,y,nsamples*sizeof(double)); 
		return 0 ;
	}

	// Cnt
	int npos = 0 ;
	for (int i=0; i<nsamples; i++) {
		if (y[i] > 0)
			npos++ ;
	}
	int nneg = nsamples - npos ;

	int onpos,onneg ;
	if ((nneg+0.0)/npos > ratio) {
		onpos = npos ;
		onneg = (int) (ratio*npos) ;
	} else {
		onneg = nneg ;
		onpos = (int) (nneg/ratio) ;
	}
	
	*fnsamples = onpos+onneg ;
	fprintf(stderr,"Will filter to %d pos + %d neg samples (Ratio = %f -> %f)\n",onpos,onneg,(nneg+0.0)/npos,ratio) ;

	// Order
	int *order ;
	if ((order = randomize(nsamples))==NULL)
		return -1 ;

	npos = nneg = 0 ;
	int nout = 0 ;
	for (int i=0; i<nsamples; i++) {
		if ((y[order[i]] > 0 && npos++ < onpos) || (y[order[i]] <= 0 && nneg++ < onneg)) {
			for (int j=0; j<nftrs; j++)
				fx[XIDX(j,nout,*fnsamples)] = x[XIDX(j,order[i],nsamples)] ;
			fy[nout++] = y[order[i]] ;
		}
	}

	free(order) ;
	return 0 ;
}

int filter(double *x, double *y, int nsamples, int nftrs, double ratio, double *fx, double *fy, int *fnsamples) {

	// Do nothing if ratio <= 0
	if (ratio <= 0) {
		*fnsamples = nsamples ;
		memcpy(fx,x,nsamples*nftrs*sizeof(double)) ;
		memcpy(fy,y,nsamples*sizeof(double)); 
		return 0 ;
	}

	// Cnt
	int npos = 0 ;
	for (int i=0; i<nsamples; i++) {
		if (y[i] > 0)
			npos++ ;
	}
	int nneg = nsamples - npos ;

	int onpos,onneg ;
	if ((nneg+0.0)/npos > ratio) {
		onpos = npos ;
		onneg = (int) (ratio*npos) ;
	} else {
		onneg = nneg ;
		onpos = (int) (nneg/ratio) ;
	}
	
	*fnsamples = onpos+onneg ;
	fprintf(stderr,"Will filter to %d pos + %d neg samples (Ratio = %f -> %f)\n",onpos,onneg,(nneg+0.0)/npos,ratio) ;

	// Order
	int *order ;
	if ((order = randomize(nsamples))==NULL)
		return -1 ;

	npos = nneg = 0 ;
	int nout = 0 ;
	for (int i=0; i<nsamples; i++) {
		if ((y[order[i]] > 0 && npos++ < onpos) || (y[order[i]] <= 0 && nneg++ < onneg)) {
			for (int j=0; j<nftrs; j++)
				fx[XIDX(nout,j,nftrs)] = x[XIDX(order[i],j,nftrs)] ;
			fy[nout++] = y[order[i]] ;
		}
	}

	free(order) ;
	return 0 ;
}

int filter(double *x, double *y, char *ids, int nsamples, int nftrs, double ratio, double **fx, double **fy, char **fids, int *fnsamples) {

	// Do nothing if ratio <= 0
	if (ratio <= 0) {
		*fnsamples = nsamples ;

		*fx = (double *) malloc(nsamples * nftrs * sizeof(double)) ;
		*fy = (double *) malloc(nsamples * sizeof(double)) ;
		*fids = (char *) malloc(nsamples * MAX_STRING_LEN) ;
		if (*fx == NULL || *fy == NULL || *fids == NULL)
			return -1 ;

		memcpy(fx,x,nsamples*nftrs*sizeof(double)) ;
		memcpy(fy,y,nsamples*sizeof(double)); 
		memcpy(fids,ids,nsamples*MAX_STRING_LEN) ;
		return 0 ;
	}

	// Cnt
	int npos = 0 ;
	for (int i=0; i<nsamples; i++) {
		if (y[i] > 0)
			npos++ ;
	}
	int nneg = nsamples - npos ;

	int onpos,onneg ;
	if ((nneg+0.0)/npos > ratio) {
		onpos = npos ;
		onneg = (int) (ratio*npos) ;
	} else {
		onneg = nneg ;
		onpos = (int) (nneg/ratio) ;
	}
	
	*fnsamples = onpos+onneg ;
	fprintf(stderr,"Will filter to %d pos + %d neg samples (Ratio = %f -> %f)\n",onpos,onneg,(nneg+0.0)/npos,ratio) ;

	*fx = (double *) malloc((*fnsamples) * nftrs * sizeof(double)) ;
	*fy = (double *) malloc((*fnsamples) * sizeof(double)) ;
	*fids = (char *) malloc((*fnsamples) * MAX_STRING_LEN) ;
	if (*fx == NULL || *fy == NULL || *fids == NULL)
		return -1 ;

	// Order
	int *order ;
	if ((order = randomize(nsamples))==NULL)
		return -1 ;

	npos = nneg = 0 ;
	int nout = 0 ;
	for (int i=0; i<nsamples; i++) {
		if ((y[order[i]] > 0 && npos++ < onpos) || (y[order[i]] <= 0 && nneg++ < onneg)) {
			for (int j=0; j<nftrs; j++)
				(*fx)[XIDX(nout,j,nftrs)] = x[XIDX(order[i],j,nftrs)] ;
			(*fy)[nout] = y[order[i]] ;
			strcpy((*fids)+HIDX(nout),ids+HIDX(order[i])) ;
			nout ++ ;
		}
	}

	free(order) ;
	return 0 ;
}

int filter(double *x, double *y, int nsamples, int nftrs, double ratio, double *fx, double *fy, int *fnsamples, int *finds) {


	// Do nothing if ratio <= 0
	if (ratio <= 0) {
		*fnsamples = nsamples ;
		memcpy(fx,x,nsamples*nftrs*sizeof(double)) ;
		memcpy(fy,y,nsamples*sizeof(double)); 

		for (int i=0; i<nsamples; i++)
			finds[i] = i ;

		return 0 ;
	}

	for (int i=0; i<nsamples; i++)
		finds[i] = -1 ;

	// Cnt
	int npos = 0 ;
	for (int i=0; i<nsamples; i++) {
		if (y[i] > 0)
			npos++ ;
	}
	int nneg = nsamples - npos ;

	int onpos,onneg ;
	if ((nneg+0.0)/npos > ratio) {
		onpos = npos ;
		onneg = (int) (ratio*npos) ;
	} else {
		onneg = nneg ;
		onpos = (int) (nneg/ratio) ;
	}
	
	*fnsamples = onpos+onneg ;
	fprintf(stderr,"Will filter to %d pos + %d neg samples (Ratio = %f -> %f)\n",onpos,onneg,(nneg+0.0)/npos,ratio) ;

	// Order
	int *order ;
	if ((order = randomize(nsamples))==NULL)
		return -1 ;

	npos = nneg = 0 ;
	int nout = 0 ;
	for (int i=0; i<nsamples; i++) {
		if ((y[order[i]] > 0 && npos++ < onpos) || (y[order[i]] <= 0 && nneg++ < onneg)) {
			for (int j=0; j<nftrs; j++)
				fx[XIDX(nout,j,nftrs)] = x[XIDX(order[i],j,nftrs)] ;
			finds[order[i]] = nout ;
			fy[nout++] = y[order[i]] ;
		}
	}

	free(order) ;
	return 0 ;
}

int filter(double **x, double **y, int *nsamples, int nftrs, double ratio) {

	// Do nothing if ratio <= 0
	if (ratio <= 0)
		return 0 ;

	// Cnt
	int npos = 0 ;
	for (int i=0; i<*nsamples; i++) {
		if ((*y)[i] > 0)
			npos++ ;
	}
	int nneg = *nsamples - npos ;

	int onpos,onneg ;
	if (nneg/npos > ratio) {
		onpos = npos ;
		onneg = (int) (ratio*npos) ;
	} else {
		onneg = nneg ;
		onpos = (int) (nneg/ratio) ;
	}
	
	int onsamples = onpos+onneg ;
	fprintf(stderr,"Will filter to %d pos + %d neg samples (Ratio = %f)\n",onpos,onneg,ratio) ;

	// Allocate
	double *fx,*fy ;
	if ((fx = (double *) malloc (nftrs*onsamples*sizeof(double)))==NULL || (fy = (double *) malloc (onsamples*sizeof(double)))==NULL)
		return -1 ;

	// Order
	int *order ;
	if ((order = randomize(*nsamples))==NULL)
		return -1 ;

	npos = nneg = 0 ;
	int nout = 0 ;
	for (int i=0; i<(*nsamples); i++) {
		if (((*y)[order[i]] > 0 && npos++ < onpos) || ((*y)[order[i]] <= 0 && nneg++ < onneg)) {
			for (int j=0; j<nftrs; j++)
				fx[XIDX(nout,j,nftrs)] = (*x)[XIDX(order[i],j,nftrs)] ;
			fy[nout++] = (*y)[order[i]] ;
		}
	}

	// Finish
	free(*x) ;
	free(*y) ;

	*x = fx ;
	*y = fy ;
	*nsamples = onsamples ;

	return 0 ;
}

int tfilter(double **x, double **y, int *nsamples, int nftrs, double ratio) {

	// Do nothing if ratio <= 0
	if (ratio <= 0)
		return 0 ;

	// Cnt
	int npos = 0 ;
	for (int i=0; i<*nsamples; i++) {
		if ((*y)[i] > 0)
			npos++ ;
	}
	int nneg = *nsamples - npos ;

	int onpos,onneg ;
	if (nneg/npos > ratio) {
		onpos = npos ;
		onneg = (int) (ratio*npos) ;
	} else {
		onneg = nneg ;
		onpos = (int) (nneg/ratio) ;
	}
	
	int onsamples = onpos+onneg ;
	fprintf(stderr,"Will filter to %d pos + %d neg samples (Ratio = %f)\n",onpos,onneg,ratio) ;

	// Allocate
	double *fx,*fy ;
	if ((fx = (double *) malloc (nftrs*onsamples*sizeof(double)))==NULL || (fy = (double *) malloc (onsamples*sizeof(double)))==NULL)
		return -1 ;

	// Order
	int *order ;
	if ((order = randomize(*nsamples))==NULL)
		return -1 ;

	npos = nneg = 0 ;
	int nout = 0 ;
	for (int i=0; i<(*nsamples); i++) {
		if (((*y)[order[i]] > 0 && npos++ < onpos) || ((*y)[order[i]] <= 0 && nneg++ < onneg)) {
			for (int j=0; j<nftrs; j++)
				fx[XIDX(j,nout,onsamples)] = (*x)[XIDX(j,order[i],(*nsamples))] ;
			fy[nout++] = (*y)[order[i]] ;
		}
	}

	// Finish
	free(*x) ;
	free(*y) ;
	free(order) ;

	*x = fx ;
	*y = fy ;
	*nsamples = onsamples ;

	return 0 ;
}

// Take a subset of columns and rows
int get_submatrix(double *xtable, double *ytable, char *header, char *ids, int *nrows, int *ncols, int *take_rows, int ntake_rows, int *take_cols, int ntake_cols) {

	for (int i=0; i<ntake_cols; i++) {
		if (take_cols[i] < i)
			return -1 ;		
	}

	for (int i=0; i<ntake_rows; i++) {
		if (take_rows[i] < i)
			return -1 ;
	}

	for (int i=0; i<ntake_cols; i++)
		strncpy(header+HIDX(i),header+HIDX(take_cols[i]),MAX_STRING_LEN) ;

	for (int i=0; i<ntake_rows; i++)
		strncpy(ids+HIDX(i),ids+HIDX(take_rows[i]),MAX_STRING_LEN) ;

	for (int i=0; i<ntake_rows; i++) {
		ytable[i] = ytable[take_rows[i]] ;
		for (int j=0; j<ntake_cols; j++)
			xtable[XIDX(i,j,ntake_cols)] = xtable[XIDX(take_rows[i],take_cols[j],*ncols)] ;
	}

	*nrows = ntake_rows ;
	*ncols = ntake_cols ;

	return 0 ;
}

// Take a subset of columns
int get_cols_submatrix(double *xtable, char *header, int nrows, int *ncols, int *take_cols, int ntake_cols) {

	for (int i=0; i<ntake_cols; i++) {
		if (take_cols[i] < i)
			return -1 ;		
	}

	for (int i=0; i<ntake_cols; i++)
		strncpy(header+HIDX(i),header+HIDX(take_cols[i]),MAX_STRING_LEN) ;


	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ntake_cols; j++)
			xtable[XIDX(i,j,ntake_cols)] = xtable[XIDX(i,take_cols[j],*ncols)] ;
	}

	*ncols = ntake_cols ;

	return 0 ;
}

void tsplit_by_flags(double *xtable, double *labels, int *flags, int nrows, int ncols, double *xtable1, double *labels1, double *xtable2, double *labels2,
					int nlearn) {

	int ntest = nrows - nlearn ;
	int ilearn = 0 ;
	int itest = 0 ;

	for (int i=0; i<nrows; i++) {
		if (flags[i] == 0) {
			for (int j=0; j<ncols; j++)
				xtable1[XIDX(j,ilearn,nlearn)] = xtable[XIDX(j,i,nrows)] ;
			labels1[ilearn++] = labels[i] ;
		} else {
			for (int j=0; j<ncols; j++)
				xtable2[XIDX(j,itest,ntest)] = xtable[XIDX(j,i,nrows)] ;
			labels2[itest++] = labels[i] ;
		}
	}

	return ;
}

void tsplit_by_flags(double *xtable, double *labels, char *ids, int *flags, int nrows, int ncols, double *xtable1, double *labels1, char *ids1, double *xtable2, 
					 double *labels2, char *ids2, int nlearn) {

	int ntest = nrows - nlearn ;
	int ilearn = 0 ;
	int itest = 0 ;

	for (int i=0; i<nrows; i++) {
		if (flags[i] == 0) {
			for (int j=0; j<ncols; j++)
				xtable1[XIDX(j,ilearn,nlearn)] = xtable[XIDX(j,i,nrows)] ;
			strcpy(ids1+HIDX(ilearn),ids+HIDX(i)) ;
			labels1[ilearn++] = labels[i] ;
		} else {
			for (int j=0; j<ncols; j++)
				xtable2[XIDX(j,itest,ntest)] = xtable[XIDX(j,i,nrows)] ;
			strcpy(ids2+HIDX(itest),ids+HIDX(i)) ;
			labels2[itest++] = labels[i] ;
		}
	}

	return ;
}