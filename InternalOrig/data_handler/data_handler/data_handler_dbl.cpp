// Data handler : Reading, Processing, and Printing data matrices (using doubles)

#include "data_handler.h"

// Work space
char dbl_buf[BUF_SIZE] ;
char dbl_field[MAX_STRING_LEN] ;
char dbl_line[DH_MAX_COLS][MAX_STRING_LEN] ;

#define _CRT_SECURE_NO_WARNINGS

using namespace std ;

// Convert text table to clean data.
int convert_data(char *text_table, int nrow, int ncol, double **xtable, double **ytable, char **headers, int *cols_to_read, int npatient, int nvar, int ycol)
{
	if (((*xtable) = (double *) malloc(npatient*nvar*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", npatient,nvar) ;
		return -1 ;
	}
	memset(*xtable,0,npatient*nvar*sizeof(double)) ;

	if (((*ytable) = (double *) malloc(npatient*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate ytable for %d", nvar) ;
		free(*xtable); 
		return -1 ;
	}
	memset(*ytable,0,npatient*sizeof(double)) ;

	if (((*headers) = (char *) malloc(nvar*MAX_STRING_LEN)) == NULL) {
		fprintf (stderr,"error : cannot allocate headers for %d vars",nvar) ;
		free(*xtable) ;
		free(*ytable) ;
		return -1 ;
	}
	memset(*headers,0,nvar*MAX_STRING_LEN) ;

	
	// variable (text, xtable)

	int   xtable_cols = 0;  
	fprintf(stderr, "Working with Label : %s\n",text_table + IDX(0,ycol,ncol)) ;

	for(int j = 0; j < nvar; j++ )
	{
		int text_col = cols_to_read[j];

		if (text_col>1000)
		{
			sprintf((*headers)+HIDX(xtable_cols), "%s%d", text_table + IDX(0,(text_col-1000)/10,ncol), (text_col-1000)%10);

		} else
		{
			strcpy_s((*headers)+HIDX(xtable_cols), MAX_STRING_LEN, text_table + IDX(0,text_col,ncol));
		}

		for( int i = 0; i < npatient; i++ )
		{
			if (text_col>1000)
			{
				int real_col = (text_col-1000)/10;
				int val = (text_col-1000)%10;
				 
				if (atoi(text_table + IDX(i+1,real_col,ncol))==val) 
				{
					(*xtable)[XIDX(i,xtable_cols,nvar)] = 1; 
				} else 
				{
					(*xtable)[XIDX(i,xtable_cols,nvar)] = 0;
				}
			} else if (text_table[IDX(i+1,text_col,ncol)] =='æ') {
				(*xtable)[XIDX(i,xtable_cols,nvar)]=1; 
			} else if (text_table[IDX(i+1,text_col,ncol)] =='ð') {
				(*xtable)[XIDX(i,xtable_cols,nvar)]=0; 
			} else 	(*xtable)[XIDX(i,xtable_cols,nvar)] = atof( text_table + IDX(i+1,text_col,ncol) ) ;
		}

		xtable_cols++;
	}

	for( int i = 0; i < npatient; i++ ) {
		(*ytable)[i] = atof( text_table + IDX(i+1,ycol,ncol) );	
	}	

	return 0;
}

// Calculate statistics of xtable with weights
void weighted_calc_xstats(double *xtable, double *weights, int npatient, int nvar, double *avg, double *std)
{
	for( int i = 0; i < nvar; i++ ) {
		avg[i] = weighted_calc_col_avg(xtable,i, weights, npatient, nvar,-1.0);
		std[i] = weighted_calc_col_std(xtable,i, weights, npatient, nvar, avg[i],-1.0);
	}
}

// Calculate statistics of xtable with weights
void calc_xstats(double *xtable, int npatient, int nvar, double *avg, double *std, double missing)
{
	for( int i = 0; i < nvar; i++ ) {
		avg[i] = calc_col_avg(xtable,i, npatient, nvar, missing);
		std[i] = calc_col_std(xtable,i, npatient, nvar, avg[i], missing);
		if (std[i] == 0)
			std[i] = 1.0 ;
	}
}

// Calculate statistics of xtable and ytable
int calc_stats(double *xtable, double *ytable, int npatient, int nvar,double **avg, double **std, double *yavg, double missing)
{
	if (((*avg) = (double *) malloc(nvar*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate averages for %d\n", nvar) ;
		return -1 ;
	}
	memset(*avg,0,nvar*sizeof(double)) ;

	if (((*std) = (double *) malloc(nvar*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate stds for %d\n", nvar) ;
		free(*avg) ;
		return -1 ;
	}
	memset(*std,0,nvar*sizeof(double)) ;

	calc_xstats(xtable,npatient,nvar,*avg,*std,missing) ;

	double ysum=0; 
	for (int i=0; i<npatient; i++)
		ysum += ytable[i] ;
	(*yavg) = ysum/npatient ;

	return 0;
}

// Calculate statistics of central part of each vector
int calc_partial_stats(double *xtable, int npatient, int nvar, double central_p, double *avg, double *std, double missing) {

	double *vec = (double *) malloc(npatient*sizeof(double)) ;
	if (vec == NULL) {
		fprintf(stderr,"Allocation Failed\n") ;
		return -1 ;
	}

	for (int i=0; i<nvar; i++) {
		int n=0 ;
		for (int j=0; j<npatient; j++) {
			if (xtable[XIDX(j,i,nvar)] != missing)
				vec[n++] = xtable[XIDX(j,i,nvar)] ;
		}

		qsort(vec,n,sizeof(double),double_compare) ;

		int start = (int) (n*(1-central_p)/2) ;
		int neff = (int) (n*central_p) ;

		if ((vec[start]==vec[start+neff-1]) || get_moments(vec + start,neff,avg+i,std+i,missing) == -1) {
			avg[i] = -1.0 ;
			std[i] = -1.0 ;
		}
	}

	return 0 ;
}

// Calculate statistics of xtable and ytable with weights
int weighted_calc_stats(double *xtable, double *ytable, double *weights, int npatient, int nvar,double **avg, double **std, double *yavg)
{

	if (((*avg) = (double *) malloc(nvar*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate averages for %d\n", nvar) ;
		return -1 ;
	}
	memset(*avg,0,nvar*sizeof(double)) ;

	if (((*std) = (double *) malloc(nvar*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate stds for %d\n", nvar) ;
		free(*avg) ;
		return -1 ;
	}
	memset(*std,0,nvar*sizeof(double)) ;

	weighted_calc_xstats(xtable,weights,npatient,nvar,*avg,*std) ;

	double ysum=0; 
	double ynum = 0 ;	
	for (int i=0; i<npatient; i++) {
		ynum += weights[i] ;
		ysum += weights[i]*ytable[i] ;
	}

	(*yavg) = ysum/ynum ;

	return 0;
}

// Check outliers and replace with max/min allowed values
int outliers(double *xtable, double *avg, double *std, int npatients, int nvar, int ivar, double max_std, double missing) {
	 
	int nclear = 0 ;
	for( int i = 0; i < npatients; i++ ) {
		if (xtable[XIDX(i,ivar,nvar)]==missing) continue;
		if( xtable[XIDX(i,ivar,nvar)] > avg[ivar] +  max_std * std[ivar] ) {
			xtable[XIDX(i,ivar,nvar)] = avg[ivar] +  (max_std-1) * std[ivar] ;
			nclear++ ;
		}
		if( xtable[XIDX(i,ivar,nvar)] < avg[ivar] - max_std * std[ivar] ) {
			xtable[XIDX(i,ivar,nvar)] = avg[ivar] - (max_std-1) * std[ivar] ;
			nclear++ ;
		}
	}

	fprintf(stderr,"Outliers ; Adjusted %d entries\n",nclear) ;

	avg[ivar] = calc_col_avg(xtable,ivar, npatients, nvar,missing);
	std[ivar] = calc_col_std(xtable,ivar, npatients, nvar, avg[ivar],missing);
	if (std[ivar] == 0)
		std[ivar] = 1 ;

	return nclear ;
}

int outliers(double *xtable, double *avg, double *std, int npatient, int nvar, int *check, int check_num, double max_std, double missing)
{
	 
	int nclear = 0 ;
	for( int i = 0; i < npatient; i++ ) {
		for( int k = 0; k < check_num; k++ ) {
			
			int j = check[k] ;
			 
			if (xtable[XIDX(i,j,nvar)]==missing) continue;
			if( xtable[XIDX(i,j,nvar)] > avg[j] +  max_std * std[j] ) {
				xtable[XIDX(i,j,nvar)] = avg[j] +  (max_std-1) * std[j] ;
				nclear++ ;
			}
			if( xtable[XIDX(i,j,nvar)] < avg[j] - max_std * std[j] ) {
				xtable[XIDX(i,j,nvar)] = avg[j] - (max_std-1) * std[j] ;
				nclear++ ;
			}
		}
	}

	fprintf(stderr,"Outliers ; Adjusted %d entries\n",nclear) ;

	calc_xstats(xtable,npatient,nvar,avg,std,missing);
	return nclear ;
}

int outliers(double *xtable, double *avg, double *std, int npatient, int nvar, int *fcheck, double max_std, double missing)
{
	 
	int nclear = 0 ;
	for( int i = 0; i < npatient; i++ ) {
		for( int j=0; j<nvar; j++) {
			if (fcheck[j]) {
			 
				if (xtable[XIDX(i,j,nvar)]==missing) continue;
				if( xtable[XIDX(i,j,nvar)] > avg[j] +  max_std * std[j] ) {
					xtable[XIDX(i,j,nvar)] = avg[j] +  (max_std-1) * std[j] ;
					nclear++ ;
				}
				if( xtable[XIDX(i,j,nvar)] < avg[j] - max_std * std[j] ) {
					xtable[XIDX(i,j,nvar)] = avg[j] - (max_std-1) * std[j] ;
					nclear++ ;
				}
			}
		}
	}

	fprintf(stderr,"Outliers ; Adjusted %d entries\n",nclear) ;

	calc_xstats(xtable,npatient,nvar,avg,std,missing);
	return nclear ;
}

// Check outliers and remove (adjusting weights)
int remove_outliers_with_weights(double **xtable, double **ytable, double **weights, double *avg, double *std, int *npatient, int nvar, int *check, int check_num, double max_std)
{
	 
	int out_i = 0 ;
	for( int i = 0; i < (*npatient); i++ ) {

		int remove = 0 ;
		for( int k = 0; k < check_num; k++ ) {
			int j = check[k] ;

			if( (*xtable)[XIDX(i,j,nvar)] > avg[j] +  max_std * std[j] ) {
				remove = 1 ;
			}
			if( (*xtable)[XIDX(i,j,nvar)] < avg[j] - max_std * std[j] ) {
				remove = 1 ;
			}
		}

		if (remove == 0) {
			if (out_i < i) {
				for (int j = 0; j<nvar; j++)
					(*xtable)[XIDX(out_i,j,nvar)] = (*xtable)[XIDX(i,j,nvar)] ;
				(*ytable)[out_i] = (*ytable)[i] ;
				(*weights)[out_i] = (*weights)[i] ;
			}
			out_i++ ;
		}
	}

	if (out_i == 0)
		return 0 ;

	if (out_i <= *npatient) {
		*npatient = out_i ;
		if (((*xtable) = (double *) realloc(*xtable,(*npatient)*nvar*sizeof(double))) == NULL) {
			fprintf(stderr,"Cannot reallocate data for %d samples\n",*npatient) ;
			return -1 ;
		}

		if (((*ytable) = (double *) realloc(*ytable,(*npatient)*sizeof(double))) == NULL) {
			fprintf(stderr,"Cannot reallocate labels for %d samples\n",*npatient) ;
			return -1 ;
		}

		if (((*weights) = (double *) realloc(*weights,(*npatient)*sizeof(double))) == NULL) {
			fprintf(stderr,"Cannot reallocate weights for %d samples\n",*npatient) ;
			return -1 ;
		}
	}

//	fprintf(stderr,"After removing outliers - left with %d samples\n",*npatient) ;
	calc_xstats(*xtable,*npatient,nvar,avg,std);
	return 0 ;
}

// Check outliers and remove
int remove_outliers(double **xtable, double **ytable, double *avg, double *std, int *npatient, int nvar, int *check, int check_num, double max_std)
{

	double *weights ;
	if ((weights = (double *) malloc((*npatient) * sizeof(double))) == NULL) {
		printf ("error : cannot allocate dummy weights\n") ;
		return -1 ;
	}

	for (int i=0; i<*npatient; i++)
		weights[i] = 1.0 ;

	int rc = remove_outliers_with_weights(xtable,ytable,&weights,avg,std,npatient,nvar,check,check_num,max_std) ;

	free(weights) ;
	return rc ;
}		

// Clear data - run several iterations of checking and marking outliers with -1
void clear_data(double *xtable, double *avg, double *std, int npatient, int nvar, int niter, int *check, int check_num, double max_std, double missing)
{

	// Prepare mean+std
	calc_xstats(xtable,npatient,nvar,avg,std,missing) ;

	// Iteratively clear all required variables
	for (int i=0; i<check_num; i++) {
		int ivar = check[i] ;
		fprintf(stderr,"Clearing var %d (%d/%d) with %f/%f\n",ivar,i,check_num,avg[ivar],std[ivar]) ;

		for(int it=0; it<niter-1; it++) {
			if (outliers(xtable,avg,std,npatient,nvar,ivar,max_std,missing)==0)
				break ;
		}

		outliers(xtable,avg,std,npatient,nvar,ivar,max_std,missing) ;
	}
}

void clear_data(double *xtable, double *avg, double *std, int npatient, int nvar, int niter, int *fcheck, double max_std)
{
	for(int it=0; it<niter-1; it++) {
		if (outliers(xtable,avg,std,npatient,nvar,fcheck,max_std)==0)
			return ;
	}

// 	fprintf(stderr,"--------------- LAST ONE ---------\n");
	outliers(xtable,avg,std,npatient,nvar,fcheck,max_std) ;
}

// Clear data - Check marked features
int clear_data(double *xtable, double *avg, double *std, int npatient, int nvar, int niter, double max_std, bool *mask, double missing_val)
{

	int ncheck = 0 ;
	for (int i=0; i<nvar; i++) {
		if (mask[i])
			ncheck ++ ;
	}

	if (ncheck == 0)
		return 0 ;

	int *check ;
	if ((check = (int *) malloc (ncheck*sizeof(int)))==NULL) {
		fprintf(stderr,"Cannot allocate dummy check_cols\n") ;
		return -1 ;
	}

	int icheck = 0;
	for (int i=0; i<nvar; i++) {
		if (mask[i])
			check[icheck++] = i ;
	}

	clear_data(xtable,avg,std,npatient,nvar,niter,check,ncheck,max_std,missing_val) ;
	return 0 ;
}
// Clear data - Check all features
int clear_data(double *xtable, double *avg, double *std, int npatient, int nvar, int niter, double max_std)
{

	int ncheck = nvar ;
	int *check ;

	if ((check = (int *) malloc (ncheck*sizeof(int)))==NULL) {
		fprintf(stderr,"Cannot allocate dummy check_cols\n") ;
		return -1 ;
	}

	for (int i=0; i<nvar; i++)
		check[i] = i ;

	clear_data(xtable,avg,std,npatient,nvar,niter,check,ncheck,max_std) ;
	return 0 ;
}

// Clear data - run several iterations of checking and removing outliers
int aggresive_clear_data(double **xtable, double **ytable, double *avg, double *std, int *npatient, int nvar, int niter, int *check, int check_num,double max_std)
{
	for(int it=0; it<niter-1; it++) {
		if (remove_outliers(xtable,ytable,avg,std,npatient,nvar,check,check_num,max_std) == -1) 
			return -1 ;
	}

	if ((*npatient == 0)) {
		fprintf(stderr,"Aggressive Data Clearing : Nothing left. Quitting\n") ;
		return 0 ;
	}

// 	fprintf(stderr,"--------------- LAST ONE ---------\n");

	if (remove_outliers(xtable,ytable,avg,std,npatient,nvar,check,check_num,max_std) == -1)
		return -1 ;

	return 0 ;
}

// Clear data - run several iterations of checking and removing outliers (adjusting weights)
int aggresive_clear_data_with_weights(double **xtable, double **ytable, double **weights,double *avg, double *std, int *npatient, int nvar, int niter, int *check, int check_num,double max_std)
{
	for(int it=0; it<niter-1; it++) {
		if (remove_outliers_with_weights(xtable,ytable,weights,avg,std,npatient,nvar,check,check_num,max_std) == -1) 
			return -1 ;
	}

	if ((*npatient == 0)) {
		fprintf(stderr,"Aggressive Data Clearing : Nothing left. Quitting\n") ;
		return 0 ;
	}

 //	fprintf(stderr,"--------------- LAST ONE ---------\n");

	if (remove_outliers_with_weights(xtable,ytable,weights,avg,std,npatient,nvar,check,check_num,max_std) == -1)
		return -1 ;

	return 0 ;
}

// Print summary of data
int print_data(char *headers, double *avg, double *std, int npatient, int nvar)
{
	printf("\nPatients number = %d\n\n", npatient );

	printf( "    <%15s> (%8s, %8s)\n", "variable", "avg", "std" );
	
	for( int i = 0; i < nvar; i++ ) {
		printf( "%02d. <%15s> (%.3f, %.3f)\n", i, headers+HIDX(i), avg[i], std[i]);
	}

	return 0;
}

// Print data in a matrix format into file
int print_matrix(double *xtable, double *ytable, char *headers, int npatient, int nvar, const char *file_name, double *weights) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++)
		fprintf(fp,"%s\t",headers+HIDX(i)) ;
	fprintf(fp,"Label\n") ;

	for (int i=0; i<npatient; i++) {
		if (weights[i] > 0) {
			for (int j=0; j<nvar; j++)
				fprintf(fp,"%f\t",xtable[XIDX(i,j,nvar)]) ;
			fprintf(fp,"%f\n",ytable[i]) ;
		}
	}

	fclose(fp) ;
	return 0 ;
}

int print_matrix(double *xtable, char *headers, int npatient, int nvar, const char *file_name, double *weights) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar-1; i++)
		fprintf(fp,"%s\t",headers+HIDX(i)) ;
	fprintf(fp,"%s\n",headers+HIDX(nvar-1)) ;

	for (int i=0; i<npatient; i++) {
		if (weights[i] > 0) {
			for (int j=0; j<nvar - 1; j++)
				fprintf(fp,"%f\t",xtable[XIDX(i,j,nvar)]) ;
			fprintf(fp,"%f\n",xtable[XIDX(i,nvar-1,nvar)]) ;
		}
	}

	fclose(fp) ;
	return 0 ;
}

int print_matrix(double *xtable, double *ytable, char *headers, int npatient, int nvar, const char *file_name) {

	double *weights ;
	if ((weights = (double *) malloc (npatient*sizeof(double))) == NULL) {
		fprintf(stderr,"weight allocation failed\n") ;
		return -1 ;
	}

	for (int i=0; i<npatient; i++)
		weights[i]=1.0 ;

	return print_matrix(xtable,ytable,headers,npatient,nvar,file_name,weights) ;
}	

int print_matrix(double *xtable, char *headers, int npatient, int nvar, const char *file_name) {

	double *weights ;
	if ((weights = (double *) malloc (npatient*sizeof(double))) == NULL) {
		fprintf(stderr,"weight allocation failed\n") ;
		return -1 ;
	}

	for (int i=0; i<npatient; i++)
		weights[i]=1.0 ;

	return print_matrix(xtable,headers,npatient,nvar,file_name,weights) ;
}	

// Printing Matrix without header.
int print_matrix(double *xtable, double *ytable, int npatient, int nvar, const char *file_name) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++)
		fprintf(fp,"Var%d\t",i) ;
	fprintf(fp,"Label\n") ;

	for (int i=0; i<npatient; i++) {
		for (int j=0; j<nvar; j++)
			fprintf(fp,"%f\t",xtable[XIDX(i,j,nvar)]) ;
		fprintf(fp,"%f\n",ytable[i]) ;
	}

	fclose(fp) ;
	return 0 ;
}

int tprint_matrix(double *xtable, double *ytable, int npatient, int nvar, const char *file_name) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++)
		fprintf(fp,"Var%d\t",i) ;
	fprintf(fp,"Label\n") ;

	for (int i=0; i<npatient; i++) {
		for (int j=0; j<nvar; j++)
			fprintf(fp,"%f\t",xtable[XIDX(j,i,npatient)]) ;
		fprintf(fp,"%f\n",ytable[i]) ;
	}

	fclose(fp) ;
	return 0 ;
}

int tprint_matrix(double *xtable, double *ytable, int npatient, int nvar, int *ignore, const char *file_name) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++) {
		if (! ignore[i])
			fprintf(fp,"Var%d\t",i) ;
	}
	fprintf(fp,"Label\n") ;

	for (int i=0; i<npatient; i++) {
		for (int j=0; j<nvar; j++) {
			if (!ignore[j])
				fprintf(fp,"%f\t",xtable[XIDX(j,i,npatient)]) ;
		}
		fprintf(fp,"%f\n",ytable[i]) ;
	}

	fclose(fp) ;
	return 0 ;
}

int tprint_matrix(double *xtable, double *ytable, int npatient, int nvar, int *ignore, const char *file_name, double missing) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++) {
		if (! ignore[i])
			fprintf(fp,"Var%d\t",i) ;
	}
	fprintf(fp,"Label\n") ;

	for (int i=0; i<npatient; i++) {
		for (int j=0; j<nvar; j++) {
			if (!ignore[j]) {
				if (xtable[XIDX(j,i,npatient)]==missing)
					fprintf(fp,"NA\t") ;
				else
					fprintf(fp,"%f\t",xtable[XIDX(j,i,npatient)]) ;
			}
		}
		fprintf(fp,"%f\n",ytable[i]) ;
	}

	fclose(fp) ;
	return 0 ;
}


int print_matrix(double *xtable, int npatient, int nvar, const char *file_name) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++) {
		fprintf(fp,"Var%d",i) ;
		if (i != nvar -1)
			fprintf(fp,"\t") ;
	}
	fprintf(fp,"\n") ;

	for (int i=0; i<npatient; i++) {
		for (int j=0; j<nvar; j++) {
			fprintf(fp,"%f",xtable[XIDX(i,j,nvar)]) ;
			if (j != nvar -1)
				fprintf(fp,"\t") ;
		}
		fprintf(fp,"\n") ;
	}

	fclose(fp) ;
	return 0 ;
}

int tprint_matrix(double *xtable, int npatient, int nvar, const char *file_name) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++) {
		fprintf(fp,"Var%d",i) ;
		if (i != nvar -1)
			fprintf(fp,"\t") ;
	}
	fprintf(fp,"\n") ;

	for (int i=0; i<npatient; i++) {
		for (int j=0; j<nvar; j++) {
			fprintf(fp,"%f",xtable[XIDX(j,i,npatient)]) ;
			if (j != nvar -1)
				fprintf(fp,"\t") ;
		}
		fprintf(fp,"\n") ;
	}

	fclose(fp) ;
	return 0 ;
}

int tprint_matrix(double *xtable, int npatient, int from, int to, int nvar, const char *file_name) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++) {
		fprintf(fp,"Var%d",i) ;
		if (i != nvar -1)
			fprintf(fp,"\t") ;
	}
	fprintf(fp,"\n") ;

	for (int i=from; i<=to; i++) {
		for (int j=0; j<nvar; j++) {
			fprintf(fp,"%f",xtable[XIDX(j,i,npatient)]) ;
			if (j != nvar -1)
				fprintf(fp,"\t") ;
		}
		fprintf(fp,"\n") ;
	}

	fclose(fp) ;
	return 0 ;
}

int tprint_matrix(double *xtable, int npatient, int nvar, const char *file_name, double missing) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nvar; i++) {
		fprintf(fp,"Var%d",i) ;
		if (i != nvar -1)
			fprintf(fp,"\t") ;
	}
	fprintf(fp,"\n") ;

	for (int i=0; i<npatient; i++) {
		for (int j=0; j<nvar; j++) {
			if (xtable[XIDX(j,i,npatient)]==missing)
				fprintf(fp,"NA") ;
			else
				fprintf(fp,"%f",xtable[XIDX(j,i,npatient)]) ;

			if (j != nvar -1)
				fprintf(fp,"\t") ;
		}
		fprintf(fp,"\n") ;
	}

	fclose(fp) ;
	return 0 ;
}

int tprint_matrix(double *xtable, int npatient, int nvar, int *ignore, const char *file_name) {

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) writing to file %s", errno, strerror(errno), file_name );
		return -1;
	}

	int last  ;
	for (int i=nvar-1; i>=0; i--) {
		if (! ignore[i]) {
			last = i ;
			break ;
		}
	}

	for (int i=0; i<nvar; i++) {
		if (! ignore[i]) {
			fprintf(fp,"Var%d",i) ;
			if (i != last)
				fprintf(fp,"\t") ;
		}
	}
	fprintf(fp,"\n") ;

	for (int i=0; i<npatient; i++) {
		for (int j=0; j<nvar; j++) {
			if (! ignore[j]) {
				fprintf(fp,"%f",xtable[XIDX(j,i,npatient)]) ;
				if (j != last)
					fprintf(fp,"\t") ;
			}
		}
		fprintf(fp,"\n") ;
	}

	fclose(fp) ;
	return 0 ;
}

// Normalize data - set means to zero.
void normalize_data(double *xtable, double *ytable, int npatient, int nvar, double *avg, double yavg, double missing)
{
	// Reduce average from each column 


	for(int i=0; i<npatient; i++)
	{
		for(int j=0; j<nvar; j++)
		{
			if (xtable[XIDX(i,j,nvar)]==missing)
			{
				xtable[XIDX(i,j,nvar)] = 0;
			} else
			{
				xtable[XIDX(i,j,nvar)]-=avg[j];
			}
		}

		ytable[i]-=yavg;
	}
}

void normalize_data_keeping_missing(double *xtable, double *ytable, int npatient, int nvar, double *avg, double yavg, double missing)
{
	// Reduce average from each column 


	for(int i=0; i<npatient; i++)
	{
		for(int j=0; j<nvar; j++)
		{
			if (xtable[XIDX(i,j,nvar)]!=missing)
			{
				xtable[XIDX(i,j,nvar)]-=avg[j];
			}
		}

		ytable[i]-=yavg;
	}
}

// Normalize data - set means to zero and sdv to ~1
void normalize_data(double *xtable, double *ytable, int npatient, int nvar, double *avg, double *sdv, double yavg)
{
	// Reduce average from each column 


	for(int i=0; i<npatient; i++)
	{
		for(int j=0; j<nvar; j++)
		{
			if (xtable[XIDX(i,j,nvar)]==-1)
			{
				xtable[XIDX(i,j,nvar)] = 0;
			} else
			{
				if (sdv[j] == 0) {
					xtable[XIDX(i,j,nvar)] -= avg[j] ;
				} else {
					xtable[XIDX(i,j,nvar)] = (xtable[XIDX(i,j,nvar)] - avg[j])/sdv[j];
				}
			}
		}

		ytable[i]-=yavg;
	}
}

// Normalize data - set means to zero and sdv to ~1
void normalize_data(double *xtable, int npatient, int nvar, double *avg, double *sdv)
{
	// Reduce average from each column 


	for(int i=0; i<npatient; i++)
	{
		for(int j=0; j<nvar; j++)
		{
			if (xtable[XIDX(i,j,nvar)]==-1)
			{
				xtable[XIDX(i,j,nvar)] = 0;
			} else
			{
				if (sdv[j] == 0) {
					xtable[XIDX(i,j,nvar)] -= avg[j] ;
				} else {
					xtable[XIDX(i,j,nvar)] = (xtable[XIDX(i,j,nvar)] - avg[j])/sdv[j];
				}
			}
		}
	}
}

// Normalize data - set means to zero 
void normalize_data(double *xtable, int npatient, int nvar, double *avg)
{
	// Reduce average from each column 


	for(int i=0; i<npatient; i++)
	{
		for(int j=0; j<nvar; j++)
		{
			if (xtable[XIDX(i,j,nvar)]==-1)
			{
				xtable[XIDX(i,j,nvar)] = 0;
			} else
			{
				xtable[XIDX(i,j,nvar)] -= avg[j] ;
			}
		}
	}
}

void normalize_data_keeping_missing(double *xtable, int npatient, int nvar, double *avg)
{
	// Reduce average from each column 


	for(int i=0; i<npatient; i++)
	{
		for(int j=0; j<nvar; j++)
		{
			if (xtable[XIDX(i,j,nvar)]!=-1)
			{
				xtable[XIDX(i,j,nvar)] -= avg[j] ;
			}
		}
	}
}

// Normalize x data - set means to zero and sdv to ~1
void normalize_xdata(double *xtable, double *ytable, int npatient, int nvar, double *avg, double *sdv)
{
	// Reduce average from each column 


	for(int i=0; i<npatient; i++)
	{
		for(int j=0; j<nvar; j++)
		{
			if (xtable[XIDX(i,j,nvar)]==-1)
			{
				xtable[XIDX(i,j,nvar)] = 0;
			} else
			{
				if (sdv[j] == 0) {
					xtable[XIDX(i,j,nvar)] -= avg[j] ;
				} else {
					xtable[XIDX(i,j,nvar)] = (xtable[XIDX(i,j,nvar)] - avg[j])/sdv[j];
				}
			}
		}
	}
}


// Split data into two files - given a list of rows.
int split_data(double *xtable, double *ytable, int nrows, int ncols, int *indices, int nindices, double **in_xtable, double **in_ytable, double **ex_xtable, double **ex_ytable) {

	// Allocation
	if (((*in_xtable) = (double *) malloc(nindices*ncols*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", nindices,ncols) ;
		return -1 ;
	}

	if (((*ex_xtable) = (double *) malloc((nrows-nindices)*ncols*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", nrows-nindices,ncols) ;
		free(*in_xtable) ;
		return -1 ;
	}

	if (((*in_ytable) = (double *) malloc(nindices*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate ytable for %d ", nindices) ;
		free(*in_xtable) ; free(*ex_xtable) ;
		return -1 ;
	}

	if (((*ex_ytable) = (double *) malloc((nrows-nindices)*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate ytable for %d ", nrows-nindices) ;
		free(*in_xtable) ; free(*ex_xtable) ; free(*in_ytable) ;
		return -1 ;
	}

	//Temporary Vector
	int *in_flags ;
	if ((in_flags =  (int *) malloc(nrows*sizeof(int))) == NULL) {
		fprintf(stderr,"error : cannot allocate flags for %d",nindices) ;
		return -1 ;
	}
	
	memset(in_flags,0,nrows*sizeof(int)) ;
	for (int i=0; i<nindices; i++)
		in_flags[indices[i]] = 1 ;

	// Split
	int indx = 0 ;
	int exdx = 0 ;

	for (int i=0; i<nrows; i++) {
		if (in_flags[i] == 1) {
			(*in_ytable)[indx] = ytable[i] ;
			for (int j=0; j<ncols; j++)
				(*in_xtable)[XIDX(indx,j,ncols)] = xtable[XIDX(i,j,ncols)] ;
			indx++ ;
		} else {
			(*ex_ytable)[exdx] = ytable[i] ;
			for (int j=0; j<ncols; j++)
				(*ex_xtable)[XIDX(exdx,j,ncols)] = xtable[XIDX(i,j,ncols)] ;
			exdx++ ;
		}
	}

	return 0 ;
}

// Create a new xtable with part of the columns
int get_cols_subset(double *all_xtable, double **xtable, char *all_header, char **header, int nrows, int ncols, int *cols_subset, int subset_size) 
{
	// Select indicated columns

	if (((*xtable) = (double *) malloc(nrows*subset_size*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", nrows,subset_size) ;
		return -1 ;
	}

	if (((*header) = (char *) malloc(subset_size*MAX_STRING_LEN)) == NULL) {
		fprintf (stderr,"error : cannot allocate headers for %d vars",subset_size) ;
		free(*xtable) ;
		return -1 ;
	}

	for (int i=0; i<subset_size; i++)
		strcpy_s((*header)+HIDX(i),MAX_STRING_LEN,all_header+HIDX(cols_subset[i])) ;
	
	for (int i=0; i<nrows; i++) {
		for (int j=0; j<subset_size; j++) {
			(*xtable)[XIDX(i,j,subset_size)] = all_xtable[XIDX(i,cols_subset[j],ncols)] ;
		}
	}

	return 0 ;
}

// Create a new xtable with product featuures
int add_products(double *xtable, char *header, int nrows, int ncols, double **out_xtable, char **out_header, int *prods, int nprods) {
	
	int ntot_cols = ncols + (nprods*(nprods-1))/2 ;

	if (((*out_xtable) = (double *) malloc(nrows*ntot_cols*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", nrows,ntot_cols) ;
		return -1 ;
	}

	if (((*out_header) = (char *) malloc(ntot_cols*MAX_STRING_LEN)) == NULL) {
		fprintf (stderr,"error : cannot allocate headers for %d vars",ntot_cols) ;
		free(*out_xtable) ;
		return -1 ;
	}

	for (int j=0; j<ncols; j++)
		strcpy_s((*out_header)+HIDX(j),MAX_STRING_LEN,header+HIDX(j)) ;

	int col = ncols ;
	for (int j=0; j<(nprods-1); j++) {
		for (int k=j+1; k<nprods; k++) {
			sprintf((*out_header)+HIDX(col++),"%sX%s",header+HIDX(prods[j]),header+HIDX(prods[k])) ;
		}
	}

	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ncols; j++)
			(*out_xtable)[XIDX(i,j,ntot_cols)] = xtable[XIDX(i,j,ncols)] ;

		int col = ncols ;
		for (int j=0; j<(nprods-1); j++) {
			for (int k=j+1; k<nprods; k++) {
				if (xtable[XIDX(i,prods[j],ncols)] == -1 || xtable[XIDX(i,prods[k],ncols)] == -1) {
					(*out_xtable)[XIDX(i,col,ntot_cols)] = -1 ;
				} else {
					(*out_xtable)[XIDX(i,col,ntot_cols)] = xtable[XIDX(i,prods[j],ncols)] * xtable[XIDX(i,prods[k],ncols)] ;
				}
				col++ ;
			}
		}
	}

	return 0 ;
}

// Create a new xtable with product featuures
int add_products(double *xtable, char *header, int nrows, int ncols, double **out_xtable, char **out_header, int *prods1, int nprods1, int *prods2, int nprods2) {
	
	int ntot_cols = ncols + (nprods1*nprods2);

	if (((*out_xtable) = (double *) malloc(nrows*ntot_cols*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", nrows,ntot_cols) ;
		return -1 ;
	}

	if (((*out_header) = (char *) malloc(ntot_cols*MAX_STRING_LEN)) == NULL) {
		fprintf (stderr,"error : cannot allocate headers for %d vars",ntot_cols) ;
		free(*out_xtable) ;
		return -1 ;
	}

	for (int j=0; j<ncols; j++)
		strcpy_s((*out_header)+HIDX(j),MAX_STRING_LEN,header+HIDX(j)) ;

	int col = ncols ;
	for (int j=0; j<nprods1; j++) {
		for (int k=0; k<nprods2; k++) {
			sprintf((*out_header)+HIDX(col++),"%sX%s",header+HIDX(prods1[j]),header+HIDX(prods2[k])) ;
		}
	}

	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ncols; j++)
			(*out_xtable)[XIDX(i,j,ntot_cols)] = xtable[XIDX(i,j,ncols)] ;

		int col = ncols ;
		for (int j=0; j<nprods1; j++) {
			for (int k=0; k<nprods2; k++) {
				if (xtable[XIDX(i,prods1[j],ncols)] == -1 || xtable[XIDX(i,prods2[k],ncols)] == -1) {
					(*out_xtable)[XIDX(i,col,ntot_cols)] = -1 ;
				} else {
					(*out_xtable)[XIDX(i,col,ntot_cols)] = xtable[XIDX(i,prods1[j],ncols)] * xtable[XIDX(i,prods2[k],ncols)] ;
				}
				col++ ;
			}
		}
	}

	return 0 ;
}

// Create a new xtable with square featuures
int add_squares(double *xtable, char *header, int nrows, int ncols, double **out_xtable, char **out_header, int *squares, int nsquares) {
	
	int ntot_cols = ncols + nsquares ;

	if (((*out_xtable) = (double *) malloc(nrows*ntot_cols*sizeof(double))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", nrows,ntot_cols) ;
		return -1 ;
	}

	if (((*out_header) = (char *) malloc(ntot_cols*MAX_STRING_LEN)) == NULL) {
		fprintf (stderr,"error : cannot allocate headers for %d vars",ntot_cols) ;
		free(*out_xtable) ;
		return -1 ;
	}

	for (int j=0; j<ncols; j++)
		strcpy_s((*out_header)+HIDX(j),MAX_STRING_LEN,header+HIDX(j)) ;

	int col = ncols ;
	for (int j=0; j<nsquares; j++)
		sprintf((*out_header)+HIDX(col++),"%s^2",header+HIDX(squares[j])) ;

	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ncols; j++)
			(*out_xtable)[XIDX(i,j,ntot_cols)] = xtable[XIDX(i,j,ncols)] ;

		int col = ncols ;
		for (int j=0; j<nsquares; j++) {
			if (xtable[XIDX(i,squares[j],ncols)] == -1) {
				(*out_xtable)[XIDX(i,col,ntot_cols)] = -1 ;
			} else {
				(*out_xtable)[XIDX(i,col,ntot_cols)] = xtable[XIDX(i,squares[j],ncols)] * xtable[XIDX(i,squares[j],ncols)] ;
			}
			col++ ;
		}
	}

	return 0 ;
}

// Read clean matrix of unknown size (with header, tab delimeted, C++ code)
int read_matrix(char *input_file_name, double **xtable, double **ytable, char **header, int *nrows, int *nvars) {

	ifstream ifs(input_file_name) ;
	if (!ifs.is_open()) {
		fprintf(stderr,"Cannot open %s for reading",input_file_name) ;
		return -1 ;
	}

	string line ;
	bool header_flag = true ;
	vector<double> data ;
	vector<string> head ;
	stringstream sline ;
	string element ;
	(*nrows) = 0 ;

	data.clear() ;
	head.clear() ;
	int blockn ;
	int blocki ;

	while (getline(ifs,line)) {
		sline << line ;

		if (header_flag) {
			header_flag = false ;
			while (getline(sline,element,'\t'))
				head.push_back(element) ;
			(*nvars) = (int) head.size() - 1;
			if (head.back() != "Label") {
				fprintf(stderr,"Last column must be a label ... Not Aborting\n") ;
//				return -1 ;
			}

			if (((*header) = (char *) malloc((*nvars)*MAX_STRING_LEN*sizeof(char))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d string\n",*nvars) ;
				return -1 ;
			}
			
			for (int i=0; i<(*nvars); i++)
				sprintf((*header)+HIDX(i),"%s",head[i].c_str()) ;

			// Allocate initial block_size ;
			if (((*ytable) = (double *) malloc(BLOCK_SIZE*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d labales\n",*nrows) ;
				return -1 ;
			}

			if (((*xtable) = (double *) malloc(BLOCK_SIZE*(*nvars)*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",*nrows,*nvars) ;
				return -1 ;
			}
			
			blockn=0 ;
			blocki=1 ;

		} else {
			(*nrows)++ ;
			if ((*nrows)%1000 == 0)
				fprintf(stderr,"Read %d lines\n",(*nrows)) ;

			blockn++ ;
			if (blockn>BLOCK_SIZE) {
				blocki++ ;
				if (((*ytable) = (double *) realloc(*ytable,BLOCK_SIZE*blocki*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d values\n",BLOCK_SIZE*blocki) ;
					return -1 ;
				}

				if (((*xtable) = (double *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}

			data.clear() ;

			while (getline(sline,element,'\t')) {
				stringstream selement(element) ;
				double value ;
				if ((selement >> value).fail()) {
					fprintf(stderr,"Cannot convert %s to double at \'%s\' (%d,%zd) after %f\n",element.c_str(),line.c_str(),*nrows,data.size(),data.back()) ;
					return -1 ;
				}
				data.push_back(value) ;
			}
			if (data.size() != (*nvars) + 1) {
				fprintf(stderr,"Inconssistent number of columsn at row %d (%zd vs. %d)\n",(*nrows),data.size(),(*nvars)+1) ;
				return -1 ;
			}

			(*ytable)[(*nrows)-1] = data.back() ;
			for (int j=0; j<(*nvars); j++) {
				(*xtable)[XIDX((*nrows)-1,j,(*nvars))] = data[j] ;
			}
		}

		sline.clear() ;
	}

	fprintf(stderr,"Finished reading %d rows of %d vars + label\n",*nrows,*nvars) ;
	return 0 ;
}

// Read clean matrix of unknown size (with header, tab delimeted, C++ code) n_label_cols last columns are labels, and we take the "label_col"-th one.
int read_matrix(char *input_file_name, double **xtable, double **ytable, char **header, int *nrows, int *nvars, int n_label_cols, int label_col, int max_samples) {

	ifstream ifs(input_file_name) ;
	if (!ifs.is_open()) {
		fprintf(stderr,"Cannot open %s for reading\n",input_file_name) ;
		return -1 ;
	}

	string line ;
	bool header_flag = true ;
	vector<double> data ;
	vector<string> head ;
	stringstream sline ;
	string element ;
	(*nrows) = 0 ;

	data.clear() ;
	head.clear() ;
	int blockn ;
	int blocki ;

	while (getline(ifs,line)) {
		sline << line ;

		if (header_flag) {
			header_flag = false ;
			while (getline(sline,element,'\t'))
				head.push_back(element) ;
			(*nvars) = (int) head.size() - n_label_cols;

			if (((*header) = (char *) malloc((*nvars)*MAX_STRING_LEN*sizeof(char))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d string\n",*nvars) ;
				return -1 ;
			}
			
			for (int i=0; i<(*nvars); i++)
				sprintf((*header)+HIDX(i),"%s",head[i].c_str()) ;

			// Allocate initial block_size ;
			if (((*ytable) = (double *) malloc(BLOCK_SIZE*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d labales\n",*nrows) ;
				return -1 ;
			}

			if (((*xtable) = (double *) malloc(BLOCK_SIZE*(*nvars)*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",*nrows,*nvars) ;
				return -1 ;
			}
			
			blockn=0 ;
			blocki=1 ;

		} else {
			(*nrows)++ ;
			if ((*nrows)%1000 == 0)
				fprintf(stderr,"Read %d lines\n",(*nrows)) ;

			blockn++ ;
			if (blockn>BLOCK_SIZE) {
				blocki++ ;
				if (((*ytable) = (double *) realloc(*ytable,BLOCK_SIZE*blocki*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d values\n",BLOCK_SIZE*blocki) ;
					return -1 ;
				}

				if (((*xtable) = (double *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}

			data.clear() ;

			while (getline(sline,element,'\t')) {
				stringstream selement(element) ;
				double value ;
				if ((selement >> value).fail()) {
					fprintf(stderr,"Cannot convert %s to double at \'%s\' (%d,%zd) after %f\n",element.c_str(),line.c_str(),*nrows,data.size(),data.back()) ;
					return -1 ;
				}
				data.push_back(value) ;
			}
			if (data.size() != (*nvars) + n_label_cols) {
				fprintf(stderr,"Inconssistent number of columsn at row %d (%zd vs. %d)\n",(*nrows),data.size(),(*nvars)+1) ;
				return -1 ;
			}

			(*ytable)[(*nrows)-1] = data[(*nvars) + label_col] ;
			for (int j=0; j<(*nvars); j++) {
				(*xtable)[XIDX((*nrows)-1,j,(*nvars))] = data[j] ;
			}
				
			if ((*nrows) == max_samples)
				break ;
		}

		sline.clear() ;
	}

	fprintf(stderr,"Finished reading %d rows of %d vars + label\n",*nrows,*nvars) ;
	return 0 ;
}


// Read clean matrix of unknown size (with header, tab delimeted) n_label_cols last columns are labels, and we take the "label_col"-th one.
int fast_read_matrix(char *input_file_name, double **xtable, double **ytable, char **header, int *nrows, int *nvars, int n_label_cols, int label_col, int max_samples) {

	FILE *fin = safe_fopen(input_file_name, "rb", false);
	if (fin == NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",input_file_name) ;
		return -1 ;
	}

	char *startbuf,*endbuf ;

	int blockn ;
	int blocki ;

	int n = DH_MAX_COLS ;
	int header_flag = 1 ;
	*nrows = 0 ;

	while(!(feof(fin))) {      
		fgets(dbl_buf, sizeof(dbl_buf), fin);
		if (feof(fin))
			break ;
			                
		startbuf = dbl_buf;             
		endbuf = dbl_buf;

		int ivars = 0 ;
		for( ;  ;  ) {
			if ((*endbuf == '\n') || (*endbuf == '\r') || (*endbuf == '\t')) {
				strncpy(dbl_field, startbuf, endbuf-startbuf);                     
				dbl_field[endbuf-startbuf]='\0';

				ivars++ ;
				if (ivars > n) {
					fprintf(stderr,"Problem reading file\n") ;
					return -1 ;
				}
					
				strncpy(dbl_line[ivars-1],dbl_field,MAX_STRING_LEN) ;
				dbl_line[ivars-1][MAX_STRING_LEN-1] = '\0' ;
				startbuf=endbuf+1 ;
			}

			if ((*endbuf == '\n') || (*endbuf == '\r'))
				break ;                      
			endbuf++;
		}

		if (header_flag) {
			header_flag = 0 ;
			n = ivars ;
			(*nvars) = n - n_label_cols ;

			if (((*header) = (char *) malloc((*nvars)*MAX_STRING_LEN*sizeof(char))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d string\n",*nvars) ;
				return -1 ;
			}
			
			for (int i=0; i<(*nvars); i++)
				strcpy((*header)+HIDX(i),dbl_line[i]) ;

			// Allocate initial block_size ;
			if (((*ytable) = (double *) malloc(BLOCK_SIZE*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d labales\n",*nrows) ;
				return -1 ;
			}

			if (((*xtable) = (double *) malloc(BLOCK_SIZE*(*nvars)*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",*nrows,*nvars) ;
				return -1 ;
			}
			
			blockn=0 ;
			blocki=1 ;

		} else {
			(*nrows)++ ;
			if ((*nrows)%1000 == 0)
				fprintf(stderr,"Read %d lines\n",(*nrows)) ;

			blockn++ ;
			if (blockn>BLOCK_SIZE) {
				blocki++ ;
				if (((*ytable) = (double *) realloc(*ytable,BLOCK_SIZE*blocki*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d values\n",BLOCK_SIZE*blocki) ;
					return -1 ;
				}

				if (((*xtable) = (double *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}


			for (int j=0; j<(*nvars); j++)
				(*xtable)[XIDX((*nrows)-1,j,(*nvars))] =  atof(dbl_line[j]) ;
			(*ytable)[(*nrows)-1] =  atof(dbl_line[(*nvars) + label_col]) ;
				
			if ((*nrows) == max_samples)
				break ;
		}
	}


	fprintf(stderr,"Finished reading %d rows of %d vars + label\n",*nrows,*nvars) ;
	return 0 ;
}

// Read clean matrix of unknown size (with header, tab delimeted) 
int fast_read_matrix(char *input_file_name, double **xtable, char **header, int *nrows, int *nvars, int max_samples) {

	FILE *fin = safe_fopen(input_file_name, "rb", false);
	if (fin == NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",input_file_name) ;
		return -1 ;
	}

	char *startbuf,*endbuf ;

	int blockn ;
	int blocki ;

	int n = DH_MAX_COLS ;
	int header_flag = 1 ;
	*nrows = 0 ;

	while(!(feof(fin))) {      
		fgets(dbl_buf, sizeof(dbl_buf), fin);
		if (feof(fin))
			break ;
			                
		startbuf = dbl_buf;             
		endbuf = dbl_buf;

		int ivars = 0 ;
		for( ;  ;  ) {
			if ((*endbuf == '\n') || (*endbuf == '\r') || (*endbuf == '\t')) {
				strncpy(dbl_field, startbuf, endbuf-startbuf);                     
				dbl_field[endbuf-startbuf]='\0';

				ivars++ ;
				if (ivars > n) {
					fprintf(stderr,"Problem reading file\n") ;
					return -1 ;
				}
					
				strncpy(dbl_line[ivars-1],dbl_field,MAX_STRING_LEN) ;
				dbl_line[ivars-1][MAX_STRING_LEN-1] = '\0' ;
				startbuf=endbuf+1 ;
			}

			if ((*endbuf == '\n') || (*endbuf == '\r'))
				break ;                      
			endbuf++;
		}

		if (header_flag) {
			header_flag = 0 ;
			(*nvars) = n = ivars ;

			if (((*header) = (char *) malloc((*nvars)*MAX_STRING_LEN*sizeof(char))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d string\n",*nvars) ;
				return -1 ;
			}
			
			for (int i=0; i<(*nvars); i++)
				strcpy((*header)+HIDX(i),dbl_line[i]) ;

			// Allocate initial block_size ;
			if (((*xtable) = (double *) malloc(BLOCK_SIZE*(*nvars)*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",*nrows,*nvars) ;
				return -1 ;
			}
			
			blockn=0 ;
			blocki=1 ;

		} else {
			(*nrows)++ ;
			if ((*nrows)%1000 == 0)
				fprintf(stderr,"Read %d lines\n",(*nrows)) ;

			blockn++ ;
			if (blockn>BLOCK_SIZE) {
				blocki++ ;

				if (((*xtable) = (double *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}


			for (int j=0; j<(*nvars); j++)
				(*xtable)[XIDX((*nrows)-1,j,(*nvars))] =  atof(dbl_line[j]) ;
				
			if ((*nrows) == max_samples)
				break ;
		}
	}


	fprintf(stderr,"Finished reading %d rows of %d vars\n",*nrows,*nvars) ;
	return 0 ;
}

int fast_read_matrix(char *input_file_name, double **xtable, int *nrows, int *nvars, int max_samples) {

	FILE *fin = safe_fopen(input_file_name, "rb", false);
	if (fin == NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",input_file_name) ;
		return -1 ;
	}

	char *startbuf,*endbuf ;

	int blockn = BLOCK_SIZE ;
	int blocki = 0;

	int n = DH_MAX_COLS ;
	*nrows = 0 ;
	*nvars  ;

	*xtable = NULL ;

	while(!(feof(fin))) {      
		fgets(dbl_buf, sizeof(dbl_buf), fin);
		if (feof(fin))
			break ;
			                
		startbuf = dbl_buf;             
		endbuf = dbl_buf;

		int ivars = 0 ;
		for( ;  ;  ) {
			if ((*endbuf == '\n') || (*endbuf == '\r') || (*endbuf == '\t')) {
				strncpy(dbl_field, startbuf, endbuf-startbuf);                     
				dbl_field[endbuf-startbuf]='\0';

				ivars++ ;
				if (ivars > n) {
					fprintf(stderr,"Problem reading file\n") ;
					return -1 ;
				}
					
				strncpy(dbl_line[ivars-1],dbl_field,MAX_STRING_LEN) ;
				dbl_line[ivars-1][MAX_STRING_LEN-1] = '\0' ;
				startbuf=endbuf+1 ;
			}

			if ((*endbuf == '\n') || (*endbuf == '\r'))
				break ;                      
			endbuf++;
		}

		if ((*nrows) == 0)
			*nvars = ivars ;
		else if (ivars != *nvars) {
			fprintf(stderr,"nvars inconsistency\n") ;
			return -1 ;
		}

		(*nrows)++ ;
		if ((*nrows)%50000 == 0)
			fprintf(stderr,"Read %d lines\n",(*nrows)) ;

		blockn++ ;
		if (blockn>BLOCK_SIZE) {
			blocki++ ;

			if (((*xtable) = (double *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
				return -1 ;
			}
			blockn=1 ;
		}


		for (int j=0; j<(*nvars); j++)
			(*xtable)[XIDX((*nrows)-1,j,(*nvars))] = (double) atof(dbl_line[j]) ;
				
		if ((*nrows) == max_samples)
			break ;
	}


	fprintf(stderr,"Finished reading %d rows of %d vars\n",*nrows,*nvars) ;
	return 0 ;
}


// Read clean matrix of unknown size (with header, tab delimeted) first col is Id, last col is label
int read_matrix(char *input_file_name, double **xtable, char **ids, double **ytable, char **header, int *nrows, int *nvars, int max_samples) {

	FILE *fin = safe_fopen(input_file_name, "rb", false);
	if (fin == NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",input_file_name) ;
		return -1 ;
	}

	char *startbuf,*endbuf ;

	int blockn ;
	int blocki ;

	int n = DH_MAX_COLS ;
	int header_flag = 1 ;
	*nrows = 0 ;

	while(!(feof(fin))) {      
		fgets(dbl_buf, sizeof(dbl_buf), fin);
		if (feof(fin))
			break ;
			                
		startbuf = dbl_buf;             
		endbuf = dbl_buf;

		int ivars = 0 ;
		for( ;  ;  ) {
			if ((*endbuf == '\n') || (*endbuf == '\r') || (*endbuf == '\t')) {
				strncpy(dbl_field, startbuf, endbuf-startbuf);                     
				dbl_field[endbuf-startbuf]='\0';

				ivars++ ;
				if (ivars > n) {
					fprintf(stderr,"Problem reading file (%d vs %d)\n",ivars,n) ;
					return -1 ;
				}
					
				strncpy(dbl_line[ivars-1],dbl_field,MAX_STRING_LEN) ;
				dbl_line[ivars-1][MAX_STRING_LEN-1] = '\0' ;
				startbuf=endbuf+1 ;
			}

			if ((*endbuf == '\n') || (*endbuf == '\r'))
				break ;                      
			endbuf++;
		}

		if (header_flag) {
			header_flag = 0 ;
			n = ivars ;
			(*nvars) = n - 2;

			if (((*header) = (char *) malloc((*nvars)*MAX_STRING_LEN*sizeof(char))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d string\n",*nvars) ;
				return -1 ;
			}
			
			// Ignore first column
			for (int i=0; i<(*nvars); i++)
				strcpy((*header)+HIDX(i),dbl_line[i+1]) ;

			// Allocate initial block_size ;
			if (((*ytable) = (double *) malloc(BLOCK_SIZE*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d labales\n",BLOCK_SIZE) ;
				return -1 ;
			}

			if (((*xtable) = (double *) malloc(BLOCK_SIZE*(*nvars)*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE,*nvars) ;
				return -1 ;
			}

			if (((*ids) = (char *) malloc(BLOCK_SIZE*MAX_STRING_LEN)) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d ids\n",BLOCK_SIZE) ;
				return -1 ;
			}
			
			blockn=0 ;
			blocki=1 ;

		} else {
			(*nrows)++ ;
			if ((*nrows)%1000 == 0)
				fprintf(stderr,"Read %d lines\n",(*nrows)) ;

			blockn++ ;
			if (blockn>BLOCK_SIZE) {
				blocki++ ;
				if (((*ytable) = (double *) realloc(*ytable,BLOCK_SIZE*blocki*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d values\n",BLOCK_SIZE*blocki) ;
					return -1 ;
				}

				if (((*xtable) = (double *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}

				if (((*ids) = (char *) realloc(*ids,BLOCK_SIZE*blocki*MAX_STRING_LEN)) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d flags\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}

			strcpy((*ids) + HIDX((*nrows)-1),dbl_line[0]) ;
			for (int j=0; j<(*nvars); j++)
				(*xtable)[XIDX((*nrows)-1,j,(*nvars))] =  atof(dbl_line[j+1]) ;
			(*ytable)[(*nrows)-1] =  atof(dbl_line[(*nvars) + 1]) ;
				
			if ((*nrows) == max_samples)
				break ;
		}
	}


	fprintf(stderr,"Finished reading %d rows of %d vars + label\n",*nrows,*nvars) ;
	return 0 ;
}

// Read clean matrix of unknown size (no header, no y, tab delimeted) first col is Id
int read_matrix(char *input_file_name, double **xtable, char **ids, int *nrows, int *nvars, int max_samples) {

	FILE *fin = safe_fopen(input_file_name, "rb", false);
	if (fin == NULL) {
		fprintf(stderr,"Cannot open %s for reading\n",input_file_name) ;
		return -1 ;
	}

	char *startbuf,*endbuf ;

	int blockn ;
	int blocki ;

	int n = DH_MAX_COLS ;
	int header_flag = 1 ;
	*nrows = 0 ;
	int first = 1 ;

	while(!(feof(fin))) {      
		fgets(dbl_buf, sizeof(dbl_buf), fin);
		if (feof(fin))
			break ;
			                
		startbuf = dbl_buf;             
		endbuf = dbl_buf;

		int ivars = 0 ;
		for( ;  ;  ) {
			if ((*endbuf == '\n') || (*endbuf == '\r') || (*endbuf == '\t')) {
				strncpy(dbl_field, startbuf, endbuf-startbuf);                     
				dbl_field[endbuf-startbuf]='\0';

				ivars++ ;
				if (ivars > n) {
					fprintf(stderr,"Problem reading file (%d vs %d)\n",ivars,n) ;
					return -1 ;
				}
					
				strncpy(dbl_line[ivars-1],dbl_field,MAX_STRING_LEN) ;
				dbl_line[ivars-1][MAX_STRING_LEN-1] = '\0' ;
				startbuf=endbuf+1 ;
			}

			if ((*endbuf == '\n') || (*endbuf == '\r'))
				break ;                      
			endbuf++;
		}

		if (first) {
			first = 0 ;
			n = ivars ;
			(*nvars) = n - 1;
			(*nrows) = 1 ;

			// Allocate initial block_size ;
			if (((*xtable) = (double *) malloc(BLOCK_SIZE*(*nvars)*sizeof(double))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE,*nvars) ;
				return -1 ;
			}

			if (((*ids) = (char *) malloc(BLOCK_SIZE*MAX_STRING_LEN)) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d ids\n",BLOCK_SIZE) ;
				return -1 ;
			}
			
			blockn=0 ;
			blocki=1 ;

		} else {
			(*nrows)++ ;
			if ((*nrows)%1000 == 0)
				fprintf(stderr,"Read %d lines\n",(*nrows)) ;

			blockn++ ;
			if (blockn>BLOCK_SIZE) {
				blocki++ ;

				if (((*xtable) = (double *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(double))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}

				if (((*ids) = (char *) realloc(*ids,BLOCK_SIZE*blocki*MAX_STRING_LEN)) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d flags\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}
		}

		strcpy((*ids) + HIDX((*nrows)-1),dbl_line[0]) ;
		for (int j=0; j<(*nvars); j++)
			(*xtable)[XIDX((*nrows)-1,j,(*nvars))] =  atof(dbl_line[j+1]) ;
				
		if ((*nrows) == max_samples)
			break ;
	}


	fprintf(stderr,"Finished reading %d rows of %d vars\n",*nrows,*nvars) ;
	return 0 ;
}

// Split (Test+Learn)
int split(double *xtable, double *labels, int nsamples, int nftrs, double *xtable1, double *labels1, double *xtable2, double *labels2, int *order, int test_start, 
		   int ntest) {

	int nlearn = nsamples - ntest ;
 

	// Learn
	for (int i=0; i<test_start; i++) {
		labels1[i] = labels[order[i]] ;
		for (int j=0; j<nftrs; j++) 
			xtable1[XIDX(i,j,nftrs)] = xtable[XIDX(order[i],j,nftrs)] ;
	}

	// Test
	for (int i=0; i<ntest; i++) {
		labels2[i] = labels[order[test_start+i]] ;
		for (int j=0; j<nftrs; j++)
			xtable2[XIDX(i,j,nftrs)] = xtable[XIDX(order[test_start+i],j,nftrs)] ;
	}

	// Learn
	for (int i=0; i<(nsamples-test_start-ntest);  i++) {
		labels1[test_start+i] = labels[order[test_start+ntest+i]] ;
		for (int j=0; j<nftrs; j++)
			xtable1[XIDX(test_start+i,j,nftrs)] = xtable[XIDX(order[test_start+ntest+i],j,nftrs)] ;
	}

	return 0 ;
}

int split(double *xtable, double *labels, int nsamples, int nftrs, double *xtable1, double *labels1, double *xtable2, double *labels2, int nlearn) {

	// Random ordering
	int *order ;
	if ((order = randomize(nsamples))==NULL)
		return -1 ;

	// Learn
	for (int i=0; i<nlearn; i++) {
		labels1[i] = labels[order[i]] ;
		for (int j=0; j<nftrs; j++)
			xtable1[XIDX(i,j,nftrs)] = xtable[XIDX(order[i],j,nftrs)] ;
	}

	// Test
	for (int i=nlearn; i<nsamples; i++) {
		labels2[i-nlearn] = labels[order[i]] ;
		for (int j=0; j<nftrs; j++)
			xtable2[XIDX(i-nlearn,j,nftrs)] = xtable[XIDX(order[i],j,nftrs)] ;
	}

	free(order) ;
	return 0 ;
}

// Split, keeping all samples from a each patient in one of the groups.
int pid_split(double *x, double *y, char *ids, int nrows, int ncols, double learn_ratio, double **learn_x, double **learn_y, double **test_x, double **test_y,
			  int *nlearn, char **test_ids) {

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
	for (int i=0; i< ((int) (learn_ratio*npids)); i++) {
		mark[order[i]] = 1 ;
		(*nlearn) += pid_cnts[order[i]] ;
	}

	// Allocate
	int ntest = nrows - (*nlearn) ;
	(*learn_x) = (double *) malloc((*nlearn) * ncols * sizeof(double)) ;
	(*learn_y) = (double *) malloc((*nlearn) * sizeof(double)) ;
	(*test_x) = (double *) malloc(ntest * ncols * sizeof(double)) ;
	(*test_y) = (double *) malloc(ntest * sizeof(double)) ;
	(*test_ids) = (char *) malloc(ntest * MAX_STRING_LEN) ;
	if ((*learn_x)==NULL || (*learn_y)==NULL || (*test_x)==NULL || (*test_y)==NULL || (test_ids)==NULL) {
		fprintf(stderr,"Cannot allocate datasets\n") ;
		return -1 ;
	}

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
				(*learn_x)[XIDX(ilearn,j,ncols)] = x[XIDX(i,j,ncols)] ;
			(*learn_y)[ilearn++] = y[i] ;
		} else {
			for (int j=0; j<ncols; j++)
				(*test_x)[XIDX(itest,j,ncols)] = x[XIDX(i,j,ncols)] ;
			strcpy((*test_ids)+HIDX(itest),ids+HIDX(i)) ;
			(*test_y)[itest++] = y[i] ;
		}

		prev_pid = pid ;
	}

	// Clean
	free(mark) ;
	free(pid_cnts) ;
	free(order) ;

	return 0 ;
}