// Data handler : Reading, Processing, and Printing data matrices (using floats)

#include "data_handler.h"

#define _CRT_SECURE_NO_WARNINGS

// Work space
char buf[BUF_SIZE] ;
char field[MAX_STRING_LEN] ;
char line[DH_MAX_COLS][MAX_STRING_LEN] ;

using namespace std ;

// Read data from test file
int read_xl_file(char *file_name, int max_nrow, int ncol, char **text_table, int *nrow)
{
	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "rb", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) reading file %s\n", errno, strerror(errno), file_name );
		return -1;
	}

	if (((*text_table) = (char *) malloc(ncol*max_nrow*MAX_STRING_LEN)) == NULL) {
		fprintf (stderr,"error : cannot allocate text table for %s\n", file_name) ;
		return -1 ;
	}
	
	memset(*text_table,0,ncol*max_nrow*MAX_STRING_LEN) ;

	//read variable values
	(*nrow) = 0 ;

	for( int i = 0; i < max_nrow + 1; i++ ) {
		for( int j = 0; j < ncol; j++ ) {
			char value[MAX_STRING_LEN] = "";
			ret = fscanf( fp, "%s\t", value );
			if (ret == EOF && j == 0)
				break ;

			if (i == max_nrow) {
				fprintf(stderr,"Error : Too many input lines at %s\n",file_name) ;
				return -1 ;
			}

			if( ret < 1 ) {
				fprintf( stderr,"error : ret = %d reading file (%s) at (%d/%d, %d/%d)\n", ret, file_name, i, max_nrow, j , ncol);
				return -1;
			}
			strcpy_s( (*text_table) + IDX(i,j,ncol), MAX_STRING_LEN, value);
		}
		if (ret == EOF)
			break ;

		(*nrow)++ ;
	}
		
	ret = fclose( fp );
	if( ret != 0 ) {
		fprintf( stderr, "error (%d) closing file %s", errno, file_name );
		return -2;
	}

	fprintf( stderr, "Done reading file!\n\n" );
	return ret;
}

// Read part of data file (reading line shift,shift+part,shift+2*part,....)
int read_part_of_xl_file(char *file_name, int nrow, int ncol, char **text_table, int part, int shift)
{
	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "rb", false)) == NULL) {
		fprintf( stderr,"error (%d) (%s) reading file %s\n", errno, strerror(errno), file_name );
		return -1;
	}

	int part_nrow = (nrow-1)/part + 1 ;
	if ((nrow-1)%part >= shift)
		part_nrow++ ;

	if (((*text_table) = (char *) malloc(ncol*part_nrow*MAX_STRING_LEN)) == NULL) {
		fprintf (stderr,"error : cannot allocate text table for %s\n", file_name) ;
		return -1 ;
	}
	
	memset(*text_table,0,ncol*part_nrow*MAX_STRING_LEN) ;

	//read variable values
	int irow = 0 ;

	for( int i = 0; i < nrow ; i++ ) {
		if (i==0 || (i-1)%part == shift) {
			for( int j = 0; j < ncol; j++ ) {
				char value[MAX_STRING_LEN] = "";
				ret = fscanf( fp, "%s\t", value );
				if( ret < 1 ) {
					fprintf( stderr, "error (%d) (%s) reading file (%s) at (%d, %d)\n",
						errno, strerror(errno), file_name, i, j );
					return -1;
				}
				strcpy_s( (*text_table) + IDX(irow,j,ncol), MAX_STRING_LEN, value);
			}
			irow++ ;
		}
	}
		
	ret = fclose( fp );
	if( ret != 0 ) {
		fprintf( stderr, "error (%d) closing file %s", errno, file_name );
		return -2;
	}

	fprintf( stderr, "Done reading file!\n\n" );
	return ret;
}

// Fill blanks in data file
int fill_blanks(char *text_table, int nrow, int ncol )
{
	int count = 0;
	for( int i = 0; i < nrow ; i++ ) {
		for( int j = 0; j < ncol ; j++ ) {
			if( text_table[IDX(i,j,ncol)] == '.' ) {
				count++;
				strcpy( text_table + IDX(i,j,ncol), "-1" );
			}
		}
	}

	fprintf( stderr, "Filled %d missing values.\n\n", count );
	return count;
}

// Filter table - require the column 'col' to be 'value'
int filter_table(char *text_table, int nrow, int ncol, char **filtered_table, int col, char *value, int *onrows) {

	if (((*filtered_table) = (char *) malloc(nrow*ncol*MAX_STRING_LEN)) == NULL) {
		fprintf (stderr, "error : cannot allocate filtered table for %d rows\n", nrow) ;
		return -1 ;
	}

	int irow = 0 ;
	for (int i=0; i<nrow; i++) {
		if (i==0 || strcmp(text_table + IDX(i,col,ncol),value) != 0) {
			for (int j=0; j<ncol; j++)
				strcpy((*filtered_table)+IDX(irow,j,ncol),text_table+IDX(i,j,ncol)) ;
			irow++ ;
		}
	}

	(*onrows) = irow ;
	return 0 ;
}

// Print table into file
int print_table(char *file_name, char *table, int nrows, int ncols) { 

	FILE *fp = NULL;
	int ret = 0;

	if ((fp = safe_fopen(file_name, "w", false)) == NULL) {
		fprintf( stderr, "error (%d) (%s) writing file %s\n", errno, strerror(errno), file_name );
		return -1;
	}

	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ncols-1; j++)
			fprintf(fp,"%s\t",table + IDX(i,j,ncols)) ;
		fprintf(fp,"%s\n",table + IDX(i,ncols-1,ncols)) ;
	}

	return 0 ;
}

// Convert text table to clean data.
int convert_data(char *text_table, int nrow, int ncol, float **xtable, float **ytable, char **headers, int *cols_to_read, int npatient, int nvar, int ycol)
{
	if (((*xtable) = (float *) malloc(npatient*nvar*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", npatient,nvar) ;
		return -1 ;
	}
	memset(*xtable,0,npatient*nvar*sizeof(float)) ;

	if (((*ytable) = (float *) malloc(npatient*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate ytable for %d", nvar) ;
		free(*xtable); 
		return -1 ;
	}
	memset(*ytable,0,npatient*sizeof(float)) ;

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
			} else 	(*xtable)[XIDX(i,xtable_cols,nvar)] = (float)atof( text_table + IDX(i+1,text_col,ncol) ) ;
		}

		xtable_cols++;
	}

	for( int i = 0; i < npatient; i++ ) {
		(*ytable)[i] = (float)atof( text_table + IDX(i+1,ycol,ncol) );	
	}	

	return 0;
}

// Calculate statistics of central part of each vector
int calc_partial_stats(float *xtable, int npatient, int nvar, double central_p, float *avg, float *std, float missing) {

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

		qsort(vec,n,sizeof(float),float_compare) ;

		int start = (int) (n*(1-central_p)/2) ;
		int neff = (int) (n*central_p) ;

		double davg,dstd ;
		if (get_moments(vec + start,neff,&davg,&dstd,missing) == -1) {
			avg[i] = -1.0 ;
			std[i] = -1.0 ;
		} else {
			avg[i] = (float) davg ;
			std[i] = (float) dstd ;
		}
	}

	return 0 ;
}

// Calculate statistics of xtable with weights
void weighted_calc_xstats(float *xtable, float *weights, int npatient, int nvar, float *avg, float *std)
{
	for( int i = 0; i < nvar; i++ ) {
		avg[i] = weighted_calc_col_avg(xtable,i, weights, npatient, nvar, -1.0);
		std[i] = weighted_calc_col_std(xtable,i, weights, npatient, nvar, avg[i], -1.0);
		if (std[i] == 0)
			std[i] = 1 ;
	}
}

// Calculate statistics of xtable with weights
void calc_xstats(float *xtable, int npatient, int nvar, float *avg, float *std)
{
	for( int i = 0; i < nvar; i++ ) {
		avg[i] = calc_col_avg(xtable,i, npatient, nvar,-1.0);
		std[i] = calc_col_std(xtable,i, npatient, nvar, avg[i],-1.0);
		if (std[i] == 0)
			std[i] = 1 ;
	}
}

// Calculate statistics of xtable and ytable
int calc_stats(float *xtable, float *ytable, int npatient, int nvar,float **avg, float **std, float *yavg)
{
	if (((*avg) = (float *) malloc(nvar*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate averages for %d\n", nvar) ;
		return -1 ;
	}
	memset(*avg,0,nvar*sizeof(float)) ;

	if (((*std) = (float *) malloc(nvar*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate stds for %d\n", nvar) ;
		free(*avg) ;
		return -1 ;
	}
	memset(*std,0,nvar*sizeof(float)) ;

	calc_xstats(xtable,npatient,nvar,*avg,*std) ;

	float ysum=0; 
	for (int i=0; i<npatient; i++)
		ysum += ytable[i] ;
	(*yavg) = ysum/npatient ;

	return 0;
}

// Calculate statistics of xtable and ytable with weights
int weighted_calc_stats(float *xtable, float *ytable, float *weights, int npatient, int nvar,float **avg, float **std, float *yavg)
{

	if (((*avg) = (float *) malloc(nvar*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate averages for %d\n", nvar) ;
		return -1 ;
	}
	memset(*avg,0,nvar*sizeof(float)) ;

	if (((*std) = (float *) malloc(nvar*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate stds for %d\n", nvar) ;
		free(*avg) ;
		return -1 ;
	}
	memset(*std,0,nvar*sizeof(float)) ;

	weighted_calc_xstats(xtable,weights,npatient,nvar,*avg,*std) ;

	float ysum=0; 
	float ynum = 0 ;	
	for (int i=0; i<npatient; i++) {
		ynum += weights[i] ;
		ysum += weights[i]*ytable[i] ;
	}

	(*yavg) = ysum/ynum ;

	return 0;
}

// Check outliers and replace with max/min allowed values
int outliers(float *xtable, float *avg, float *std, int npatient, int nvar, int *check, int check_num, float max_std)
{
	 
	int nclear = 0 ;
	for( int i = 0; i < npatient; i++ ) {
		for( int k = 0; k < check_num; k++ ) {
			
			int j = check[k] ;
			 
			if (xtable[XIDX(i,j,nvar)]==-1) continue;
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

	calc_xstats(xtable,npatient,nvar,avg,std);
	return nclear ;
}

// Check outliers and remove (adjusting weights)
int remove_outliers_with_weights(float **xtable, float **ytable, float **weights, float *avg, float *std, int *npatient, int nvar, int *check, int check_num, float max_std)
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
		if (((*xtable) = (float *) realloc(*xtable,(*npatient)*nvar*sizeof(float))) == NULL) {
			fprintf(stderr,"Cannot reallocate data for %d samples\n",*npatient) ;
			return -1 ;
		}

		if (((*ytable) = (float *) realloc(*ytable,(*npatient)*sizeof(float))) == NULL) {
			fprintf(stderr,"Cannot reallocate labels for %d samples\n",*npatient) ;
			return -1 ;
		}

		if (((*weights) = (float *) realloc(*weights,(*npatient)*sizeof(float))) == NULL) {
			fprintf(stderr,"Cannot reallocate weights for %d samples\n",*npatient) ;
			return -1 ;
		}
	}

//	fprintf(stderr,"After removing outliers - left with %d samples\n",*npatient) ;
	calc_xstats(*xtable,*npatient,nvar,avg,std);
	return 0 ;
}

// Check outliers and remove
int remove_outliers(float **xtable, float **ytable, float *avg, float *std, int *npatient, int nvar, int *check, int check_num, float max_std)
{

	float *weights ;
	if ((weights = (float *) malloc((*npatient) * sizeof(float))) == NULL) {
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
void clear_data(float *xtable, float *avg, float *std, int npatient, int nvar, int niter, int *check, int check_num, float max_std)
{
	for(int it=0; it<niter-1; it++) {
		if (outliers(xtable,avg,std,npatient,nvar,check,check_num,max_std)==0)
			return ;
	}

// 	fprintf(stderr,"--------------- LAST ONE ---------\n");
	outliers(xtable,avg,std,npatient,nvar,check,check_num,max_std) ;
}

// Clear data - run several iterations of checking and removing outliers
int aggresive_clear_data(float **xtable, float **ytable, float *avg, float *std, int *npatient, int nvar, int niter, int *check, int check_num,float max_std)
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
int aggresive_clear_data_with_weights(float **xtable, float **ytable, float **weights,float *avg, float *std, int *npatient, int nvar, int niter, int *check, int check_num,float max_std)
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
int print_data(char *headers, float *avg, float *std, int npatient, int nvar)
{
	printf("\nPatients number = %d\n\n", npatient );

	printf( "    <%15s> (%8s, %8s)\n", "variable", "avg", "std" );
	
	for( int i = 0; i < nvar; i++ ) {
		printf( "%02d. <%15s> (%.3f, %.3f)\n", i, headers+HIDX(i), avg[i], std[i]);
	}

	return 0;
}

// Print data in a matrix format into file
int print_matrix(float *xtable, float *ytable, char *headers, int npatient, int nvar, const char *file_name, float *weights) {

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

int print_matrix(float *xtable, float *ytable, char *headers, int npatient, int nvar, const char *file_name) {

	float *weights ;
	if ((weights = (float *) malloc (npatient*sizeof(float))) == NULL) {
		fprintf(stderr,"weight allocation failed\n") ;
		return -1 ;
	}

	for (int i=0; i<npatient; i++)
		weights[i]=1.0 ;

	return print_matrix(xtable,ytable,headers,npatient,nvar,file_name,weights) ;
}	

// Normalize data - set means to zero.
void normalize_data(float *xtable, float *ytable, int npatient, int nvar, float *avg, float yavg)
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
				xtable[XIDX(i,j,nvar)]-=avg[j];
			}
		}

		ytable[i]-=yavg;
	}
}

// Normalize data - set means to zero and sdv to ~1
void normalize_data(float *xtable, float *ytable, int npatient, int nvar, float *avg, float *sdv, float yavg)
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


// Split data into two files - given a list of rows.
int split_data(float *xtable, float *ytable, int nrows, int ncols, int *indices, int nindices, float **in_xtable, float **in_ytable, float **ex_xtable, float **ex_ytable) {

	// Allocation
	if (((*in_xtable) = (float *) malloc(nindices*ncols*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", nindices,ncols) ;
		return -1 ;
	}

	if (((*ex_xtable) = (float *) malloc((nrows-nindices)*ncols*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate xtable for %d x %d", nrows-nindices,ncols) ;
		free(*in_xtable) ;
		return -1 ;
	}

	if (((*in_ytable) = (float *) malloc(nindices*sizeof(float))) == NULL) {
		fprintf (stderr,"error : cannot allocate ytable for %d ", nindices) ;
		free(*in_xtable) ; free(*ex_xtable) ;
		return -1 ;
	}

	if (((*ex_ytable) = (float *) malloc((nrows-nindices)*sizeof(float))) == NULL) {
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
int get_cols_subset(float *all_xtable, float **xtable, char *all_header, char **header, int nrows, int ncols, int *cols_subset, int subset_size) 
{
	// Select indicated columns

	if (((*xtable) = (float *) malloc(nrows*subset_size*sizeof(float))) == NULL) {
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
int add_products(float *xtable, char *header, int nrows, int ncols, float **out_xtable, char **out_header, int *prods, int nprods) {
	
	int ntot_cols = ncols + (nprods*(nprods-1))/2 ;

	if (((*out_xtable) = (float *) malloc(nrows*ntot_cols*sizeof(float))) == NULL) {
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
int add_products(float *xtable, char *header, int nrows, int ncols, float **out_xtable, char **out_header, int *prods1, int nprods1, int *prods2, int nprods2) {
	
	int ntot_cols = ncols + (nprods1*nprods2);

	if (((*out_xtable) = (float *) malloc(nrows*ntot_cols*sizeof(float))) == NULL) {
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
int add_squares(float *xtable, char *header, int nrows, int ncols, float **out_xtable, char **out_header, int *squares, int nsquares) {
	
	int ntot_cols = ncols + nsquares ;

	if (((*out_xtable) = (float *) malloc(nrows*ntot_cols*sizeof(float))) == NULL) {
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
int read_matrix(char *input_file_name, float **xtable, float **ytable, char **header, int *nrows, int *nvars) {

	ifstream ifs(input_file_name) ;
	if (!ifs.is_open()) {
		fprintf(stderr,"Cannot open %s for reading",input_file_name) ;
		return -1 ;
	}

	string line ;
	bool header_flag = true ;
	vector<float> data ;
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
			if (((*ytable) = (float *) malloc(BLOCK_SIZE*sizeof(float))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d labales\n",*nrows) ;
				return -1 ;
			}

			if (((*xtable) = (float *) malloc(BLOCK_SIZE*(*nvars)*sizeof(float))) == NULL) {
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
				if (((*ytable) = (float *) realloc(*ytable,BLOCK_SIZE*blocki*sizeof(float))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d values\n",BLOCK_SIZE*blocki) ;
					return -1 ;
				}

				if (((*xtable) = (float *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(float))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}

			data.clear() ;

			while (getline(sline,element,'\t')) {
				stringstream selement(element) ;
				float value ;
				if ((selement >> value).fail()) {
					fprintf(stderr,"Cannot convert %s to float at \'%s\' (%d,%zd) after %f\n",element.c_str(),line.c_str(),*nrows,data.size(),data.back()) ;
					return -1 ;
				}
				data.push_back((float) value) ;
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
int read_matrix(char *input_file_name, float **xtable, float **ytable, char **header, int *nrows, int *nvars, int n_label_cols, int label_col, int max_samples) {

	ifstream ifs(input_file_name) ;
	if (!ifs.is_open()) {
		fprintf(stderr,"Cannot open %s for reading\n",input_file_name) ;
		return -1 ;
	}

	string line ;
	bool header_flag = true ;
	vector<float> data ;
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
			if (((*ytable) = (float *) malloc(BLOCK_SIZE*sizeof(float))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d labales\n",*nrows) ;
				return -1 ;
			}

			if (((*xtable) = (float *) malloc(BLOCK_SIZE*(*nvars)*sizeof(float))) == NULL) {
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
				if (((*ytable) = (float *) realloc(*ytable,BLOCK_SIZE*blocki*sizeof(float))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d values\n",BLOCK_SIZE*blocki) ;
					return -1 ;
				}

				if (((*xtable) = (float *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(float))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}

			data.clear() ;

			while (getline(sline,element,'\t')) {
				stringstream selement(element) ;
				float value ;
				if ((selement >> value).fail()) {
					fprintf(stderr,"Cannot convert %s to float at \'%s\' (%d,%zd) after %f\n",element.c_str(),line.c_str(),*nrows,data.size(),data.back()) ;
					return -1 ;
				}
				data.push_back((float) value) ;
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

// Set learning and testing sets for cross validation
void set_cv_indices(int nsamples, int nfold, int ifold, int *order, int *learn, int *nlearn, int *test, int *ntest) {
	
	// Nfold == 1 means that we don't really do CV ...
	if (nfold == 1) {
		(*nlearn) = (*ntest) = nsamples ;
		for (int j=0; j<nsamples; j++)
			learn[j] = test[j] = order[j] ;
	} else {
		float size = (float) (nsamples + 0.0)/nfold ;

		int test_from = (int) (size*ifold) ;
		int test_to = (int) (size*(ifold+1) - 1) ;
		if (test_to + size >= nsamples)
			test_to = nsamples-1 ;

		*ntest = (test_to-test_from)+1 ;
		*nlearn = nsamples-(*ntest) ;

		for (int j=0; j<test_from; j++)
			learn[j]=order[j] ;

		for (int j=test_from; j<=test_to; j++) 
			test[j-test_from] = order[j] ;
	
		for (int j=test_to+1; j<nsamples; j++) 
			learn[j-(test_to+1)+test_from] = order[j] ;

	}
	
	return ;
}

// Read clean matrix of unknown size (with header, tab delimeted) n_label_cols last columns are labels, and we take the "label_col"-th one.
int fast_read_matrix(char *input_file_name, float **xtable, float **ytable, char **header, int *nrows, int *nvars, int n_label_cols, int label_col, int max_samples) {

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
		fgets(buf, sizeof(buf), fin);
		if (feof(fin))
			break ;
			                
		startbuf = buf;             
		endbuf = buf;

		int ivars = 0 ;
		for( ;  ;  ) {
			if ((*endbuf == '\n') || (*endbuf == '\r') || (*endbuf == '\t')) {
				strncpy(field, startbuf, endbuf-startbuf);                     
				field[endbuf-startbuf]='\0';

				ivars++ ;
				if (ivars > n) {
					fprintf(stderr,"Problem reading file\n") ;
					return -1 ;
				}
					
				strncpy(line[ivars-1],field,MAX_STRING_LEN) ;
				line[ivars-1][MAX_STRING_LEN-1] = '\0' ;
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
				strcpy((*header)+HIDX(i),line[i]) ;

			// Allocate initial block_size ;
			if (((*ytable) = (float *) malloc(BLOCK_SIZE*sizeof(float))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d labales\n",*nrows) ;
				return -1 ;
			}

			if (((*xtable) = (float *) malloc(BLOCK_SIZE*(*nvars)*sizeof(float))) == NULL) {
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
				if (((*ytable) = (float *) realloc(*ytable,BLOCK_SIZE*blocki*sizeof(float))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d values\n",BLOCK_SIZE*blocki) ;
					return -1 ;
				}

				if (((*xtable) = (float *) realloc(*xtable,BLOCK_SIZE*blocki*(*nvars)*sizeof(float))) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*nvars)) ;
					return -1 ;
				}
				blockn=1 ;
			}


			for (int j=0; j<(*nvars); j++)
				(*xtable)[XIDX((*nrows)-1,j,(*nvars))] = (float) atof(line[j]) ;
			(*ytable)[(*nrows)-1] = (float) atof(line[(*nvars) + label_col]) ;
				
			if ((*nrows) == max_samples)
				break ;
		}
	}


	fprintf(stderr,"Finished reading %d rows of %d vars + label\n",*nrows,*nvars) ;
	return 0 ;
}


// Read a text table of unknown size with header (max_nsamples = optional maximal size) ;
int read_text_table(char *input_file_name, char **table, char **header, int *nrows, int *ncols, int max_samples) {
	

	fprintf(stderr,"Read %s\n",input_file_name) ;


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
		memset(buf,0,sizeof(buf)) ;
		fgets(buf, sizeof(buf), fin);
		
		if (feof(fin))
			break ;

		if (!strlen(buf))
			continue ;
			                
		startbuf = buf;             
		endbuf = buf;
		
		int ivars = 0 ;
		for( ;  ;  ) {
			if ((*endbuf == '\n') || (*endbuf == '\r') || (*endbuf == '\t')) {
				strncpy(field, startbuf, endbuf-startbuf);                     
				field[endbuf-startbuf]='\0';
				
				ivars++ ;
				if (ivars > n) {
					fprintf(stderr,"Problem reading file\n") ;
					return -1 ;
				}
					
				strncpy(line[ivars-1],field,MAX_STRING_LEN) ;
				line[ivars-1][MAX_FIELD_SIZE-1] = '\0' ; 
				startbuf=endbuf+1 ;
			}

			if ((*endbuf == '\n') || (*endbuf == '\r'))
				break ;                      
			endbuf++;
		}
		
		if (header_flag) {
			header_flag = 0 ;
			(*ncols) = ivars ;

			if (((*header) = (char *) malloc((*ncols)*MAX_STRING_LEN*sizeof(char))) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d string\n",*ncols) ;
				return -1 ;
			}
			
			for (int i=0; i<(*ncols); i++)
				strcpy((*header)+HIDX(i),line[i]) ;

			// Allocate initial block_size ;
			if (((*table) = (char *) malloc(BLOCK_SIZE*(*ncols)*MAX_FIELD_SIZE)) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",*nrows,*ncols) ;
				return -1 ;
			}
			
			blockn=0 ;
			blocki=1 ;

		} else {
			(*nrows)++ ;
			if ((*nrows)%50000 == 0)
				fprintf(stderr,"Read %d lines\n",(*nrows)) ;

			blockn++ ;
			if (blockn>BLOCK_SIZE) {
				blocki++ ;

				if (((*table) = (char *) realloc(*table,BLOCK_SIZE*blocki*(*ncols)*MAX_FIELD_SIZE)) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*ncols)) ;
					return -1 ;
				}
				blockn=1 ;
			}


			for (int j=0; j<(*ncols); j++) {
				if (strlen(line[j]) > MAX_FIELD_SIZE) {
					fprintf(stderr,"Cannot store field %s\n",line[j]) ;
					return -1 ;
				}
				strcpy_s( (*table) + FIDX((*nrows)-1,j,(*ncols)), MAX_FIELD_SIZE, line[j]);
			}
				
			if ((*nrows) == max_samples)
				break ;
		}
	}


	fprintf(stderr,"Finished reading %d rows of %d vars\n",*nrows,*ncols) ;
//	fprintf(stderr,"Vars are : \n") ;
//	for (int i=0; i<(*ncols); i++)
//		fprintf(stderr,"%d. %s\n",i,*header+HIDX(i)) ;

	fclose(fin) ;
	return 0 ;
}

// Read a text table of unknown size without header (max_nsamples = optional maximal size) ;
int read_text_table_wo_header(char *input_file_name, char **table, int *nrows, int *ncols, int max_samples) {

	fprintf(stderr,"Read %s\n",input_file_name) ;

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
		memset(buf,0,sizeof(buf)) ;
		fgets(buf, sizeof(buf), fin);
		if (feof(fin))
			break ;

		if (!strlen(buf))
			continue ;
			                
		startbuf = buf;             
		endbuf = buf;

		int ivars = 0 ;
		for( ;  ;  ) {
			if ((*endbuf == '\n') || (*endbuf == '\r') || (*endbuf == '\t')) {
				strncpy(field, startbuf, endbuf-startbuf);                     
				field[endbuf-startbuf]='\0';

				ivars++ ;
				if (ivars > n) {
					fprintf(stderr,"Problem reading file\n") ;
					return -1 ;
				}
					
				strncpy(line[ivars-1],field,MAX_STRING_LEN) ;
				line[ivars-1][MAX_FIELD_SIZE-1] = '\0' ;
				startbuf=endbuf+1 ;
			}

			if ((*endbuf == '\n') || (*endbuf == '\r'))
				break ;                      
			endbuf++;
		}

		if (header_flag) {
			header_flag = 0 ;
			(*ncols) = ivars ;
			(*nrows)++ ;

			// Allocate initial block_size ;
			if (((*table) = (char *) malloc(BLOCK_SIZE*(*ncols)*MAX_FIELD_SIZE)) == NULL) {
				fprintf(stderr,"Cannot allocate memory for %d x %d values\n",*nrows,*ncols) ;
				return -1 ;
			}
			
			blockn=0 ;
			blocki=1 ;

		} else {
			(*nrows)++ ;
			if ((*nrows)%50000 == 0)
				fprintf(stderr,"Read %d lines\n",(*nrows)) ;

			blockn++ ;
			if (blockn>BLOCK_SIZE) {
				blocki++ ;

				if (((*table) = (char *) realloc(*table,BLOCK_SIZE*blocki*(*ncols)*MAX_FIELD_SIZE)) == NULL) {
					fprintf(stderr,"Cannot allocate memory for %d x %d values\n",BLOCK_SIZE*blocki,(*ncols)) ;
					return -1 ;
				}
				blockn=1 ;
			}
		}


		for (int j=0; j<(*ncols); j++) {
			if (strlen(line[j]) > MAX_FIELD_SIZE) {
				fprintf(stderr,"Cannot store field %s\n",line[j]) ;
				return -1 ;
			}
			strcpy_s( (*table) + FIDX((*nrows)-1,j,(*ncols)), MAX_FIELD_SIZE, line[j]);
		}
			
		if ((*nrows) == max_samples)
			break ;
	}


	fprintf(stderr,"Finished reading %d rows of %d vars\n",*nrows,*ncols) ;
	return 0 ;
}

