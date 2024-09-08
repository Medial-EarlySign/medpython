#include "medial_utilities.h"
#include "globalRNG.h"

// Utilitities for vector calculations in a matrix

#define VIDX(i,j,ncol) ((i)*(ncol)+(j))

// QuickSort comparison functions
int float_compare (const void *el1, const void* el2) {
	if (*(float *)el1 > *(float *)el2) return 1;
	if (*(float *)el1 < *(float *)el2) return -1 ;
	return 0 ;
}

// Find Average of a matrix column
float calc_col_avg(float *table, int col, int nrow, int ncol, float missing)
{
	float sum = 0;
	float cntr= 0;
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=1.0;
		sum += table[VIDX(i,col,ncol)];
	}
	float result = sum / cntr;
	return result;
}

// Find Average of a subset of the rows in a matrix column
float calc_col_sub_avg(float *table, int col, int *inds, int ninds, int ncol, float missing)
{
	float sum = 0;
	float cntr= 0;
	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=1.0;
		sum += table[VIDX(inds[i],col,ncol)];
	}
	float result = sum / cntr;
	return result;
}

// Find standard deviation of a matrix column
float calc_col_std(float *table, int col, int nrow, int ncol, float missing)
{
	float result = 0;
	float cntr = 0;
	float avg = calc_col_avg(table,col,nrow,ncol,missing);
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=1.0;
		result += (table[VIDX(i,col,ncol)] - avg) * (table[VIDX(i,col,ncol)] - avg);
	}
	result = result / cntr;
	result = sqrt( result );
	return result;
}

float calc_col_std(float *table, int col, int nrow, int ncol, float avg, float missing)
{
	float result = 0;
	float cntr = 0;

	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=1.0;
		result += (table[VIDX(i,col,ncol)] - avg) * (table[VIDX(i,col,ncol)] - avg);
	}
	result = result / cntr;
	result = sqrt( result );
	return result;
}

// Find standard deviation of a subset of the rows in a matrix column
float calc_col_sub_std(float *table, int col, int *inds, int ninds, int ncol, float missing)
{
	float result = 0;
	float cntr = 0;
	float avg = calc_col_sub_avg(table,col,inds,ninds,ncol,missing);
	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=1.0;
		result += (table[VIDX(inds[i],col,ncol)] - avg) * (table[VIDX(inds[i],col,ncol)] - avg);
	}
	result = result / cntr;
	result = sqrt( result );
	return result;
}

float calc_col_sub_std(float *table, int col, int *inds, int ninds, int ncol,float avg, float missing)
{
	float result = 0;
	float cntr = 0;

	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=1.0;
		result += (table[VIDX(inds[i],col,ncol)] - avg) * (table[VIDX(inds[i],col,ncol)] - avg);
	}
	result = result / cntr;
	result = sqrt( result );
	return result;
}

// Find weighted-average of a matrix column
float weighted_calc_col_avg(float *table, int col, float *weights, int nrow, int ncol, float missing)
{
	float sum = 0;
	float cntr= 0;
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=weights[i];
		sum += weights[i]*table[VIDX(i,col,ncol)];
	}
	float result = sum / cntr;
	return result;
}

// Find weighted-average of a subset of the lines in a matrix column
float weighted_calc_col_sub_avg(float *table, int col, float *weights, int *inds, int ninds, int ncol, float missing)
{
	float sum = 0;
	float cntr= 0;
	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=weights[inds[i]];
		sum += weights[inds[i]]*table[VIDX(inds[i],col,ncol)];
	}
	float result = sum / cntr;
	return result;
}

// Find weighted-standard-deviation of matrix column
float weighted_calc_col_std(float *table, int col, float *weights, int nrow, int ncol, float missing)
{
	float result = 0;
	float cntr = 0;
	float avg = weighted_calc_col_avg(table,col,weights,nrow,ncol,missing);
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=weights[i];
		result += weights[i]*(table[VIDX(i,col,ncol)] - avg) * (table[VIDX(i,col,ncol)] - avg);
	}
	result = result / cntr;
	result = sqrt( result );
	return result;
}

float weighted_calc_col_std(float *table, int col, float *weights, int nrow, int ncol, float avg, float missing)
{
	float result = 0;
	float cntr = 0;

	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=weights[i];
		result += weights[i]*(table[VIDX(i,col,ncol)] - avg) * (table[VIDX(i,col,ncol)] - avg);
	}
	result = result / cntr;
	result = sqrt( result );
	return result;
}

// Find weighted-standard-deviation of a subset of the lines in a matrix column
float weighted_calc_col_sub_std(float *table, int col, float *weights, int *inds, int ninds, int ncol, float missing)
{
	float result = 0;
	float cntr = 0;
	float avg = weighted_calc_col_sub_avg(table,col,weights,inds,ninds,ncol,missing);
	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=weights[inds[i]];
		result += weights[inds[i]]*(table[VIDX(inds[i],col,ncol)] - avg) * (table[VIDX(inds[i],col,ncol)] - avg);
	}
	result = result / cntr;
	result = sqrt( result );
	return result;
}

float weighted_calc_col_sub_std(float *table, int col, float *weights, int *inds, int ninds, int ncol, float avg, float missing)
{
	float result = 0;
	float cntr = 0;

	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=weights[inds[i]];
		result += weights[inds[i]]*(table[VIDX(inds[i],col,ncol)] - avg) * (table[VIDX(inds[i],col,ncol)] - avg);
	}
	result = result / cntr;
	result = sqrt( result );
	return result;
}

// Find sum of a matrix column
float calc_col_sum(float *table, int col, int nrow, int ncol , float missing)
{
	float result = 0;
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		result += table[VIDX(i,col,ncol)];
	}
	return result;
}

// Shuffle a vector
int shuffle(float *vec, int n, float **new_vec) {

	if (((*new_vec) = (float *) malloc(n*sizeof(float))) == NULL) {
		printf ("error : cannot allocate vetor for %d floats", n) ;
		return -1 ;
	}

	for (int i=0; i<n; i++) 
		(*new_vec)[i] = vec[i] ;

	for (int i=0; i<n; i++) {
		float r = (float) (((float) globalRNG::rand())/(globalRNG::max() + 1.0)) ;
		int id = (int) (r*(n-i)) ;

		float temp = (*new_vec)[i+id] ;
		(*new_vec)[i+id] = (*new_vec)[i] ;
		(*new_vec)[i] = temp ;
	}

	return 0;
}

// Find Pearson correlation of two vectors
float pearson(float *vec1, float *vec2, int n) {

	float sx = 0 ;
	float sy = 0 ;
	for (int j=0; j<n; j++) {
		sx += vec1[j] ;
		sy += vec2[j];	
	}
	float meanx = sx/n ;
	float meany = sy/n ;

	float sxx = 0 ;
	float sxy = 0 ;
	float syy = 0 ;

	for (int j=0; j<n; j++) {
		sxx += (vec1[j]-meanx)*(vec1[j]-meanx) ;
		syy += (vec2[j]-meany)*(vec2[j]-meany) ;
		sxy += (vec1[j]-meanx)*(vec2[j]-meany) ;
	}

//	fprintf(stderr,"Pearson Debug : %f %f %f %f %f\n",meanx,meany,sxx,sxy,syy) ;

	if (sxx==0) sxx=1 ;
	if (syy==0) syy=1 ;

	float r2 = sxy/sqrt(sxx*syy) ;
	return r2 ;
}

// Get indices of vec2 elements in vec1
int get_indices(int *vec1, int n1, int *vec2, int n2, int **indices, int *n3) {

	if (((*indices) = (int *) malloc(n2*sizeof(int))) == NULL) {
		fprintf(stderr,"error : cannot allocate %d indices\n",n2) ;
		return -1 ;
	}

	*n3 = 0 ;
	for (int i2=0; i2<n2; i2++) {
		int index = -1 ;
		for (int i1=0; i1<n1; i1++) {
			if (vec2[i2] == vec1[i1]) {
				index = i1 ;
				break ;
			}
		}

		if (index >= 0)
			(*indices)[(*n3)++] = index ;
	}

	return 0 ;
}

// Get the (tied) order of a vector.
int get_order(float *vec, int n, float **order) {

	// Prepare
	if (((*order) = (float *) malloc (n*sizeof(float))) == NULL) {
		fprintf(stderr,"cannot allocate %d indices\n",n) ;
		return -1 ;
	}

	if (n==1) {
		(*order)[0] = 0.0 ;
		return 0 ;
	}

	struct val_idx *temp ;
	if ((temp = (struct val_idx *) malloc (n*sizeof(struct val_idx))) == NULL) {
		fprintf(stderr,"cannot allocate %d val_idx pairs\n",n) ;
		return -1 ;
	}
	
	for (int i=0; i<n; i++) {
		(temp[i]).val = vec[i] ;
		(temp[i]).idx = i ;
	}

	// Sort
	qsort(temp,n,sizeof(struct val_idx),val_idx_compare) ;

	// Parse
	int min_idx = 0 ;
	float curr_val = (temp[0]).val ;

	for (int i=1; i<n; i++) {
		if ((temp[i]).val != curr_val) {
			for (int j=min_idx; j<i; j++)
				(*order)[(temp[j]).idx] =  ((float)(min_idx + i-1))/2 ;
			min_idx = i ;
			curr_val = (temp[i]).val ;
		}
	}

	for (int j=min_idx; j<n; j++)
		(*order)[(temp[j]).idx] = ((float)(min_idx + n-1))/2 ;

	free(temp) ;
	return 0 ;
}



// Find Spearman correlation of two vectors
float spearman(float *vec1, float *vec2, int n) {

	float *order1 ;
	if (get_order(vec1,n,&order1) == -1)
		return -1 ;

	float *order2 ;
	if (get_order(vec2,n,&order2) == -1)
		return -1 ;

	float tau =  pearson(order1,order2,n) ;
	free(order1) ;
	free(order2) ;

	return tau ;
}

// Create a random permutation
int *randomize (int nrows) {

	int *order ;
	if ((order = (int *) malloc(nrows*sizeof(int))) == NULL) 
		return order ;

	for (int i=0; i<nrows ; i++)
		order[i] = i ;

	for (int i=0; i<nrows; i++) {
		float r = (float) (((float) globalRNG::rand())/(globalRNG::max() + 1.0)) ;
		int id = (int) (r*(nrows-i)) ;
	
		int temp = order[i+id] ;
		order[i+id] = order[i] ;
		order[i] = temp ;
	}

	return order ;
}

// Moments
int get_moments (double *v, int n, double *mean, double *sdv, double missing) {

	int effn = 0;
	double sum = 0 ;
	double sum2 = 0 ;

	for (int i=0; i<n; i++) {
		if (v[i] != missing) {
			effn++ ;
			sum += v[i] ;
			sum2 += v[i]*v[i] ;
		}
	}

	if (effn==0)
		return -1 ;

	(*mean) = sum/effn ;
	if (effn==1)
		(*sdv) = 0.0 ;
	else
		(*sdv) = sqrt((sum2 - (*mean)*sum)/(effn-1)) ;

	return 0 ;
}

int get_mean (double *v, int n, double *mean, double missing) {

	int effn = 0;
	double sum = 0 ;

	for (int i=0; i<n; i++) {
		if (v[i] != missing) {
			effn++ ;
			sum += v[i] ;
		}
	}

	if (effn==0)
		return -1 ;

	(*mean) = sum/effn ;
	return 0 ;
}

int get_sdv (double *v, int n, double mean, double *sdv, double missing) {

	int effn = 0;
	double sum2 = 0 ;

	for (int i=0; i<n; i++) {
		if (v[i] != missing) {
			effn++ ;
			sum2 += (v[i]-mean)*(v[i]-mean) ;
		}
	}

	if (effn==0)
		return -1 ;

	if (effn==1)
		(*sdv) = 0 ;
	else
		(*sdv) = sqrt(sum2/(effn-1)) ;

	return 0 ;
}

int get_median (double *v, int n, double *median, double missing) {

	double *tempv = (double *) malloc(n*sizeof(double)) ;
	if (tempv==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1;
	}

	int effn=0 ;
	for (int i=0; i<n; i++) {
		if (v[i] != missing)
			tempv[effn++] = v[i] ;
	}

	qsort(tempv,effn,sizeof(double),double_compare) ;
	(*median) = (tempv[effn/2] + tempv[(effn-1)/2])/2 ;

	return 0 ;
}

int get_quantiles (double *v, int n, double *qs, int nqs, double *vals, double missing) {

	double *tempv = (double *) malloc(n*sizeof(double)) ;
	if (tempv==NULL) {
		fprintf(stderr,"Allocation failed\n") ;
		return -1;
	}

	int effn=0 ;
	for (int i=0; i<n; i++) {
		if (v[i] != missing)
			tempv[effn++] = v[i] ;
	}

	qsort(tempv,effn,sizeof(double),double_compare) ;

	for (int iq=0; iq<nqs; iq++) {
		double ptr = qs[iq]*effn ;
		double p = ptr - ((int) ptr) ;
		vals[iq] = (1-p)*(tempv[(int)ptr]) + p*(tempv[1+(int)ptr]) ;
	}

	return 0 ;
}
