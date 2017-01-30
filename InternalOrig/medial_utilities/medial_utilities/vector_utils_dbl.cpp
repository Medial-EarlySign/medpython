// Utilitities for vector calculations in a matrix

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "medial_utilities.h"
#include "globalRNG.h"

#define VIDX(i,j,ncol) ((i)*(ncol)+(j))

// QuickSort comparison functions
int dbl_idx_compare (const void *el1, const void* el2) {
	if ((*(struct dbl_idx *)el1).val > (*(struct dbl_idx *)el2).val) return 1;
	if ((*(struct dbl_idx *)el1).val < (*(struct dbl_idx *)el2).val) return -11;
	return 0 ;
}

int val_idx_compare (const void *el1, const void* el2) {
	if ((*(struct val_idx *)el1).val > (*(struct val_idx *)el2).val) return 1;
	if ((*(struct val_idx *)el1).val < (*(struct val_idx *)el2).val) return -11;
	return 0 ;
}

int double_compare (const void *el1, const void* el2) {
	if (*(double *)el1 > *(double *)el2) return 1;
	if (*(double *)el1 < *(double *)el2) return -1 ;
	return 0 ;
}

int int_compare (const void *el1, const void* el2) {
	return (*(int *) el1 - *(int *) el2) ;
}

// Find Average of a matrix column
double calc_col_avg(double *table, int col, int nrow, int ncol, double missing)
{
	double sum = 0;
	double cntr= 0;
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=1.0;
		sum += table[VIDX(i,col,ncol)];
	}
	
	double result = (cntr == 0) ? 0 : sum / cntr;
	return result;
}

// Find Average of a subset of the rows in a matrix column
double calc_col_sub_avg(double *table, int col, int *inds, int ninds, int ncol, double missing)
{
	double sum = 0;
	double cntr= 0;
	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=1.0;
		sum += table[VIDX(inds[i],col,ncol)];
	}
	double result = (cntr == 0) ? 0 : sum / cntr;
	return result;
}

// Find standard deviation of a matrix column
double calc_col_std(double *table, int col, int nrow, int ncol, double missing)
{
	double result = 0;
	double cntr = 0;
	double avg = calc_col_avg(table,col,nrow,ncol,missing);
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=1.0;
		result += (table[VIDX(i,col,ncol)] - avg) * (table[VIDX(i,col,ncol)] - avg);
	}
	result = (cntr == 0) ? 0 : sqrt(result / cntr);
	return result;
}

double calc_col_std(double *table, int col, int nrow, int ncol, double avg, double missing)
{
	double result = 0;
	double cntr = 0;

	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=1.0;
		result += (table[VIDX(i,col,ncol)] - avg) * (table[VIDX(i,col,ncol)] - avg);
	}
	result = (cntr == 0) ? 0 : sqrt(result / cntr);
	return result;
}

// Find standard deviation of a subset of the rows in a matrix column
double calc_col_sub_std(double *table, int col, int *inds, int ninds, int ncol, double missing)
{
	double result = 0;
	double cntr = 0;
	double avg = calc_col_sub_avg(table,col,inds,ninds,ncol,missing);
	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=1.0;
		result += (table[VIDX(inds[i],col,ncol)] - avg) * (table[VIDX(inds[i],col,ncol)] - avg);
	}
	result = (cntr == 0) ? 0 : sqrt(result / cntr);
	return result;
}

double calc_col_sub_std(double *table, int col, int *inds, int ninds, int ncol,double avg, double missing)
{
	double result = 0;
	double cntr = 0;

	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=1.0;
		result += (table[VIDX(inds[i],col,ncol)] - avg) * (table[VIDX(inds[i],col,ncol)] - avg);
	}
	result = (cntr == 0) ? 0 : sqrt(result / cntr);
	return result;
}

// Find weighted-average of a matrix column
double weighted_calc_col_avg(double *table, int col, double *weights, int nrow, int ncol, double missing)
{
	double sum = 0;
	double cntr= 0;
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=weights[i];
		sum += weights[i]*table[VIDX(i,col,ncol)];
	}
	double result = (cntr == 0) ? 0 : (sum / cntr) ;
	return result;
}

// Find weighted-average of a subset of the lines in a matrix column
double weighted_calc_col_sub_avg(double *table, int col, double *weights, int *inds, int ninds, int ncol, double missing)
{
	double sum = 0;
	double cntr= 0;
	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=weights[inds[i]];
		sum += weights[inds[i]]*table[VIDX(inds[i],col,ncol)];
	}
	double result = (cntr == 0) ? 0 : (sum / cntr) ;
	return result;
}

// Find weighted-standard-deviation of matrix column
double weighted_calc_col_std(double *table, int col, double *weights, int nrow, int ncol, double missing)
{
	double result = 0;
	double cntr = 0;
	double avg = weighted_calc_col_avg(table,col,weights,nrow,ncol,missing);
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=weights[i];
		result += weights[i]*(table[VIDX(i,col,ncol)] - avg) * (table[VIDX(i,col,ncol)] - avg);
	}
	result = (cntr == 0) ? 0 : sqrt(result / cntr);
	return result;
}

double weighted_calc_col_std(double *table, int col, double *weights, int nrow, int ncol, double avg, double missing)
{
	double result = 0;
	double cntr = 0;

	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		cntr+=weights[i];
		result += weights[i]*(table[VIDX(i,col,ncol)] - avg) * (table[VIDX(i,col,ncol)] - avg);
	}
	result = (cntr == 0) ? 0 : sqrt(result / cntr);
	return result;
}

// Find weighted-standard-deviation of a subset of the lines in a matrix column
double weighted_calc_col_sub_std(double *table, int col, double *weights, int *inds, int ninds, int ncol, double missing)
{
	double result = 0;
	double cntr = 0;
	double avg = weighted_calc_col_sub_avg(table,col,weights,inds,ninds,ncol,missing);
	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=weights[inds[i]];
		result += weights[inds[i]]*(table[VIDX(inds[i],col,ncol)] - avg) * (table[VIDX(inds[i],col,ncol)] - avg);
	}
	result = (cntr == 0) ? 0 : sqrt(result / cntr);
	return result;
}

double weighted_calc_col_sub_std(double *table, int col, double *weights, int *inds, int ninds, int ncol, double avg, double missing)
{
	double result = 0;
	double cntr = 0;

	for( int i = 0; i < ninds; i++ ) {
		if (table[VIDX(inds[i],col,ncol)]==missing) continue;
		cntr+=weights[inds[i]];
		result += weights[inds[i]]*(table[VIDX(inds[i],col,ncol)] - avg) * (table[VIDX(inds[i],col,ncol)] - avg);
	}
	result = (cntr == 0) ? 0 : sqrt(result / cntr);
	return result;
}

// Find sum of a matrix column
double calc_col_sum(double *table, int col, int nrow, int ncol, double missing)
{
	double result = 0;
	for( int i = 0; i < nrow; i++ ) {
		if (table[VIDX(i,col,ncol)]==missing) continue;
		result += table[VIDX(i,col,ncol)];
	}
	return result;
}

// Shuffle a vector
int shuffle(double *vec, int n, double **new_vec) {

	if (((*new_vec) = (double *) malloc(n*sizeof(double))) == NULL) {
		printf ("error : cannot allocate vetor for %d doubles", n) ;
		return -1 ;
	}

	for (int i=0; i<n; i++) 
		(*new_vec)[i] = vec[i] ;

	for (int i=0; i<n; i++) {
		double r = globalRNG::rand()/(globalRNG::max() + 1.0) ;
		int id = (int) (r*(n-i)) ;

		double temp = (*new_vec)[i+id] ;
		(*new_vec)[i+id] = (*new_vec)[i] ;
		(*new_vec)[i] = temp ;
	}

	return 0;
}

// Find Pearson correlation of two vectors
double pearson(double *vec1, double *vec2, int n) {

	double sx = 0 ;
	double sy = 0 ;
	for (int j=0; j<n; j++) {
		sx += vec1[j] ;
		sy += vec2[j];	
	}
	double meanx = sx/n ;
	double meany = sy/n ;

	double sxx = 0 ;
	double sxy = 0 ;
	double syy = 0 ;

	for (int j=0; j<n; j++) {
		sxx += (vec1[j]-meanx)*(vec1[j]-meanx) ;
		syy += (vec2[j]-meany)*(vec2[j]-meany) ;
		sxy += (vec1[j]-meanx)*(vec2[j]-meany) ;
	}


	if (sxx==0) sxx=1 ;
	if (syy==0) syy=1 ;

	double r2 = sxy/sqrt(sxx*syy) ;
	return r2 ;
}

// Find Pearson correlation of two vectors (-2 if all are missing)
double pearson(double *vec1, double *vec2, int n, double missing) {

	double sx = 0 ;
	double sy = 0 ;
	int neff = 0 ;
	for (int j=0; j<n; j++) {
		if (vec1[j] != missing && vec2[j] != missing) {
			sx += vec1[j] ;
			sy += vec2[j];	
			neff ++ ;
		}
	}

	if (neff==0)
		return -2 ;

	double meanx = sx/neff ;
	double meany = sy/neff ;

	double sxx = 0 ;
	double sxy = 0 ;
	double syy = 0 ;

	for (int j=0; j<n; j++) {
		if (vec1[j] != missing && vec2[j] != missing) {
			sxx += (vec1[j]-meanx)*(vec1[j]-meanx) ;
			syy += (vec2[j]-meany)*(vec2[j]-meany) ;
			sxy += (vec1[j]-meanx)*(vec2[j]-meany) ;
		}
	}


	if (sxx==0) sxx=1 ;
	if (syy==0) syy=1 ;

	double r2 = sxy/sqrt(sxx*syy) ;
	return r2 ;
}

// Get the (tied) order of a vector.
int get_order(double *vec, int n, double **order) {

	// Prepare
	if (((*order) = (double *) malloc (n*sizeof(double))) == NULL) {
		fprintf(stderr,"cannot allocate %d indices\n",n) ;
		return -1 ;
	}

	if (n==1) {
		(*order)[0] = 0.0 ;
		return 0 ;
	}

	struct dbl_idx *temp ;
	if ((temp = (struct dbl_idx *) malloc (n*sizeof(struct dbl_idx))) == NULL) {
		fprintf(stderr,"cannot allocate %d dbl_idx pairs\n",n) ;
		return -1 ;
	}
	
	for (int i=0; i<n; i++) {
		(temp[i]).val = vec[i] ;
		(temp[i]).idx = i ;
	}

	// Sort
	qsort(temp,n,sizeof(struct dbl_idx),dbl_idx_compare) ;

	// Parse
	int min_idx = 0 ;
	double curr_val = (temp[0]).val ;

	for (int i=1; i<n; i++) {
		if ((temp[i]).val != curr_val) {
			for (int j=min_idx; j<i; j++)
				(*order)[(temp[j]).idx] =  ((double)(min_idx + i-1))/2 ;
			min_idx = i ;
			curr_val = (temp[i]).val ;
		}
	}

	for (int j=min_idx; j<n; j++)
		(*order)[(temp[j]).idx] = ((double)(min_idx + n-1))/2 ;

	free(temp) ;
	return 0 ;
}



// Find Spearman correlation of two vectors
double spearman(double *vec1, double *vec2, int n) {

	double *order1 ;
	if (get_order(vec1,n,&order1) == -1)
		return -1 ;

	double *order2 ;
	if (get_order(vec2,n,&order2) == -1)
		return -1 ;

	double tau =  pearson(order1,order2,n) ;
	free(order1) ;
	free(order2) ;

	return tau ;
}
