// Data handler : Reading, Processing, and Printing data matrices
#ifndef __DATA_HANDLER__
#define __DATA_HANDLER__

#define _CRT_SECURE_NO_WARNINGS

#pragma once

#ifndef _WIN32
#define strcpy_s(DST, SZ, SRC) strncpy((DST), (SRC), (SZ))
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "medial_utilities/medial_utilities/medial_utilities.h"

#define BLOCK_SIZE 500000
#define BUF_SIZE 200000
#define MAX_FIELDS 20

#define DH_MAX_COLS 4400 

#define DEF_NRUNS 1
#define NCLN_ITER 15
#define MAX_STD 15
#define RESOLUTION 100
#define RATIO 3

// Functions //

// Reading Data
// Read data from test file
int read_xl_file(char *file_name, int max_nrow, int ncol, char **text_table, int *nrow) ;

// Read part of data file (reading line shift,shift+part,shift+2*part,....)
int read_part_of_xl_file(char *file_name, int nrow, int ncol, char **text_table, int part, int shift) ;

// Read a text table of unknown size with header (max_nsamples = optional maximal size) ;
int read_text_table(char *input_file_name, char **table, char **header, int *nrows, int *ncols, int max_samples = 0);

// Read clean matrix of unknown size (with header, tab delimeted, C++ code)
int read_matrix(char *input_file_name, float **xtable, float **ytable, char **header, int *nrows, int *nvars) ;
int read_matrix(char *input_file_name, double **xtable, double **ytable, char **header, int *nrows, int *nvars) ;

// Read clean matrix of unknown size (with header, tab delimeted, C++ code) n_label_cols last columns are labels, and we take the "label_col"-th one.
int read_matrix(char *input_file_name, float **xtable, float **ytable, char **header, int *nrows, int *nvars, int n_label_cols, int label_col, int max_samples = 0) ;
int read_matrix(char *input_file_name, double **xtable, double **ytable, char **header, int *nrows, int *nvars, int n_label_cols, int label_col, int max_samples = 0) ;

// Read clean matrix of unknown size (with header, tab delimeted) n_label_cols last columns are labels, and we take the "label_col"-th one.
int fast_read_matrix(char *input_file_name, float **xtable, float **ytable, char **header, int *nrows, int *nvars, int n_label_cols, int label_col, int max_samples = 0) ;
int fast_read_matrix(char *input_file_name, double **xtable, double **ytable, char **header, int *nrows, int *nvars, int n_label_cols, int label_col, int max_samples = 0) ;
int fast_read_matrix(char *input_file_name, double **xtable, char **header, int *nrows, int *nvars, int max_samples = 0) ;
int fast_read_matrix(char *input_file_name, double **xtable, int *nrows, int *nvars, int max_samples = 0) ;
// Read clean matrix of unknown size (no header, no y, tab delimeted) first col is Id
int read_matrix(char *input_file_name, double **xtable, char **ids, int *nrows, int *nvars, int max_samples=0) ;

// Read clean matrix of unknown size (with header, tab delimeted) first col is Id, last col is label
int read_matrix(char *input_file_name, double **xtable, char **ids, double **ytable, char **header, int *nrows, int *nvars, int max_samples=0) ;

// Read a text table of unknown size without header (max_nsamples = optional maximal size) ;
int read_text_table_wo_header(char *input_file_name, char **table, int *nrows, int *ncols, int max_samples = 0) ;

// Convert text table to clean data.
int convert_data(char *text_table, int nrow, int ncol, float **xtable, float **ytable, char **headers, int *cols_to_read, int npatient, int nvar, int ycol) ;
int convert_data(char *text_table, int nrow, int ncol, double **xtable, double **ytable, char **headers, int *cols_to_read, int npatient, int nvar, int ycol) ;

// Fill blanks in data file
int fill_blanks(char *text_table, int nrow, int ncol) ;

// Filter table - require the column 'col' to be 'value'
int filter_table(char *text_table, int nrow, int ncol, char **filtered_table, int col, char *value, int *onrows) ;

// Printing
// Print table into file
int print_table(char *file_name, char *table, int nrows, int ncols) ;

// Print summary of data
int print_data(char *headers, float *avg, float *std, int npatient, int nvar) ;
int print_data(char *headers, double *avg, double *std, int npatient, int nvar) ;

// Print data in a matrix format into file
int print_matrix(float *xtable, float *ytable, char *headers, int npatient, int nvar, const char *file_name, float *weights)  ;
int print_matrix(double *xtable, double *ytable, char *headers, int npatient, int nvar, const char *file_name, double *weights) ;
int print_matrix(float *xtable, float *ytable, char *headers, int npatient, int nvar, const char *file_name) ;
int print_matrix(double *xtable, char *headers, int npatient, int nvar, const char *file_name, double *weights) ;
int print_matrix(double *xtable, double *ytable, char *headers, int npatient, int nvar, const char *file_name)  ;
int print_matrix(double *xtable, char *headers, int npatient, int nvar, const char *file_name)  ;

// Printing Matrix without header
int print_matrix(double *xtable, double *ytable, int npatient, int nvar, const char *file_name) ; // Basic
int tprint_matrix(double *xtable, double *ytable, int npatient, int nvar, const char *file_name) ; // Transposed input
int tprint_matrix(double *xtable, double *ytable, int npatient, int nvar, int *ignore, const char *file_name)  ; // as above, Ignoring marked variable
int tprint_matrix(double *xtable, double *ytable, int npatient, int nvar, int *ignore, const char *file_name, double missing) ; // as above, Marking missing values as NA
int print_matrix(double *xtable, int npatient, int nvar, const char *file_name) ; // Basic - no label
int tprint_matrix(double *xtable, int npatient, int nvar, const char *file_name) ; // Transposed input
int tprint_matrix(double *xtable, int npatient, int from, int to, int nvar, const char *file_name) ; // as above, printing only part of the lines
int tprint_matrix(double *xtable, int npatient, int nvar, const char *file_name, double missing) ; // Marking missing values as NA
int tprint_matrix(double *xtable, int npatient, int nvar, int *ignore, const char *file_name) ; // Ignoring marked variable

// Memory handling
// Transpose Matrix. If reverse=0, start with all features of a given sample in a line.
int transpose(double *tx, double **x, int nsamples, int nftrs, int reverse = 0) ;
int transpose(int *tx, int **x, int nsamples, int nftrs, int reverse = 0) ;

// Stats and Cleaning
// Calculate statistics of xtable with weights
void weighted_calc_xstats(float *xtable, float *weights, int npatient, int nvar, float *avg, float *std) ;
void weighted_calc_xstats(double *xtable, double *weights, int npatient, int nvar, double *avg, double *std) ;

// Calculate statistics of xtable with weights
void calc_xstats(float *xtable, int npatient, int nvar, float *avg, float *std) ;
void calc_xstats(double *xtable, int npatient, int nvar, double *avg, double *std, double missing = -1.0) ;
void tcalc_xstats(double *x, int nsamples, int nftrs, double *avg, double *std, double missing =-1.0) ;
void tcalc_xstats(double *x, double *w, int nsamples, int nftrs, double *avg, double *std, double missing=-1.0) ;

// Calculate statistics of xtable and ytable
int calc_stats(float *xtable, float *ytable, int npatient, int nvar,float **avg, float **std, float *yavg) ;
int calc_stats(double *xtable, double *ytable, int npatient, int nvar,double **avg, double **std, double *yavg, double missing = -1.0) ;
int tcalc_stats(double *x, double *y, int nsamples, int nftrs, double **avg, double **std, double *yavg, double missing=-1.0) ;
int tcalc_stats(double *x, double *y, double *w, int nsamples, int nftrs, double **avg, double **std, double *yavg, double missing=-1.0) ;

// Calculate moments on central part of data
int calc_partial_stats(float *xtable, int npatient, int nvar, double central_p, float *avg, float *std, float missing = -1.0) ;
int calc_partial_stats(double *xtable, int npatient, int nvar, double central_p, double *avg, double *std, double missing = -1.0) ;
int tcalc_partial_stats(float *xtable, int npatient, int nvar, double central_p, float *avg, float *std, float missing = -1.0) ;

// Calculate statistics of xtable and ytable with weights
int weighted_calc_stats(float *xtable, float *ytable, float *weights, int npatient, int nvar,float **avg, float **std, float *yavg) ;
int weighted_calc_stats(double *xtable, double *ytable, double *weights, int npatient, int nvar,double **avg, double **std, double *yavg) ;

// Check outliers and replace with max/min allowed values
int outliers(float *xtable, float *avg, float *std, int npatient, int nvar, int *check, int check_num, float max_std) ;
int outliers(double *xtable, double *avg, double *std, int npatient, int nvar, int *check, int check_num, double max_std, double missing = -1.0) ;
int outliers(double *xtable, double *avg, double *std, int npatient, int nvar, int *fcheck, double max_std, double missing=-1.0) ;

// Check outliers and remove (adjusting weights)
int remove_outliers_with_weights(float **xtable, float **ytable, float **weights, float *avg, float *std, int *npatient, int nvar, int *check, int check_num, float max_std) ;
int remove_outliers_with_weights(double **xtable, double **ytable, double **weights, double *avg, double *std, int *npatient, int nvar, int *check, int check_num, double max_std) ;

// Check outliers and remove
int remove_outliers(float **xtable, float **ytable, float *avg, float *std, int *npatient, int nvar, int *check, int check_num, float max_std) ;
int remove_outliers(double **xtable, double **ytable, double *avg, double *std, int *npatient, int nvar, int *check, int check_num, double max_std) ;

// Clear data - run several iterations of checking and marking outliers with -1
void clear_data(float *xtable, float *avg, float *std, int npatient, int nvar, int niter, int *check, int check_num, float max_std) ;
void clear_data(double *xtable, double *avg, double *std, int npatient, int nvar, int niter, int *check, int check_num, double max_std,double missing = -1.0);
int clear_data(double *xtable, double *avg, double *std, int npatient, int nvar, int niter, double max_std, bool *mask, double missing = -1.0) ;
void clear_data(double *xtable, double *avg, double *std, int npatient, int nvar, int niter, int *fcheck, double max_std) ;
int clear_data(double *xtable, double *avg, double *std, int npatient, int nvar, int niter, double max_std) ;

// Clear data - run several iterations of checking and removing outliers
int aggresive_clear_data(float **xtable, float **ytable, float *avg, float *std, int *npatient, int nvar, int niter, int *check, int check_num,float max_std);
int aggresive_clear_data(double **xtable, double **ytable, double *avg, double *std, int *npatient, int nvar, int niter, int *check, 
						 int check_num,double max_std) ;

// Clear data - run several iterations of checking and removing outliers (adjusting weights)
int aggresive_clear_data_with_weights(float **xtable, float **ytable, float **weights,float *avg, float *std, int *npatient, int nvar, int niter, int *check, int check_num,float max_std) ;
int aggresive_clear_data_with_weights(double **xtable, double **ytable, double **weights ,double *avg, double *std, int *npatient, int nvar, int niter, int *check, int check_num,double max_std) ;

// Normalize data - set means to zero.
void normalize_data(float *xtable, float *ytable, int npatient, int nvar, float *avg, float yavg) ;
void normalize_data(double *xtable, int npatient, int nvar, double *avg) ;
void normalize_data(double *xtable, double *ytable, int npatient, int nvar, double *avg, double yavg, double missing = -1.0) ;
void normalize_data_keeping_missing(double *xtable, double *ytable, int npatient, int nvar, double *avg, double yavg, double missing = -1.0) ;
void normalize_data_keeping_missing(double *xtable, int npatient, int nvar, double *avg) ;
void tnormalize_data(double *x, double *y, int nsamples, int nftrs, double *avg, double yavg, double missing = -1.0) ;
void tnormalize_data(double *x, int nsamples, int nftrs, double *avg, double missing = -1.0) ;
void tnormalize_data_keeping_missing(double *x, double *y, int nsamples, int nftrs, double *avg, double yavg, double missing = -1.0) ;
void tnormalize_data_keeping_missing(double *x, int nsamples, int nftrs, double *avg, double missing = -1.0) ;

// Normalize data - set means to zero and sdv to ~1
void normalize_data(float *xtable, float *ytable, int npatient, int nvar, float *avg, float *sdv, float yavg) ;
void normalize_data(double *xtable, double *ytable, int npatient, int nvar, double *avg, double *sdv, double yavg) ;
void normalize_data(double *xtable, int npatient, int nvar, double *avg, double *sdv) ;
void normalize_xdata(double *xtable, double *ytable, int npatient, int nvar, double *avg, double *sdv) ;
void tnormalize_data(double *x, double *y, int nsamples, int nftrs, double *avg, double *std, double yavg, double missing = -1.0) ;
void tnormalize_data(double *x, int nsamples, int nftrs, double *avg, double *std, double missing = -1.0) ;

// Normalize x1,y1 and x2 according to means and standard-deviation of x1,y1
int normalize_inout(double *x1, double *y1, double *w1, int nrows1, double *x2, int nrows2, int ncols, double missing_val=-1) ;
int normalize_inout(double *x1, double *y1, int nrows1, double *x2, int nrows2, int ncols,  double missing_val=-1) ;

// Cross Validation ...
// Split data into two files - given a list of rows.
int split_data(float *xtable, float *ytable, int nrows, int ncols, int *indices, int nindices, float **in_xtable, float **in_ytable, float **ex_xtable, float **ex_ytable) ;
int split_data(double *xtable, double *ytable, int nrows, int ncols, int *indices, int nindices, double **in_xtable, double **in_ytable, 
			   double **ex_xtable, double **ex_ytable) ;

// Split (Test+Learn)
int split(double *xtable, double *labels, int nsamples, int nftrs, double *xtable1, double *labels1, double *xtable2, double *labels2, int nlearn) ;
int tsplit(double *xtable, double *labels, int nsamples, int nftrs, double *xtable1, double *labels1, double *xtable2, double *labels2, int nlearn) ;
int split(double *xtable, double *labels, int nsamples, int nftrs, double *xtable1, double *labels1, double *xtable2, double *labels2, int *order, int test_start, 
		   int ntest) ;
int tsplit(double *xtable, double *labels, int nsamples, int nftrs, double *xtable1, double *labels1, double *xtable2, double *labels2, int *order, int test_start, 
		   int ntest) ;

// Split, keeping all samples from a each patient in one of the groups.
int pid_split(double *x, double *y, char *ids, int nrows, int ncols, double learn_ratio, double **learn_x, double **learn_y, double **test_x, double **test_y,
			  int *nlearn, char **test_ids) ;
int pid_tsplit(double *x, double *y, char *ids, int nrows, int ncols, int *order, double learn_ratio, double *learn_x, double *learn_y, double *test_x, double *test_y,
			  int *nlearn, char *test_ids) ;

// Split (Test+Learn) - untransposed data, with missing flags
int split(double *xtable, int *missing, double *labels, int nsamples, int nftrs, double *xtable1, int *missing1, double *labels1, 
		  double *xtable2, int *missing2, double *labels2, int nlearn) ;

// Set learning and testing sets for cross validation
void set_cv_indices(int nsamples, int nfold, int ifold, int *order, int *learn, int *nlearn, int *test, int *ntest) ;

// Splitting according to flags
int get_flags(char *ids, int nrows, double lratio, int *flags, int *nlearn) ;
int get_flags(char *ids, int nrows, int *order, double lratio, int *flags, int *nlearn) ;

// Splitting according to flags
void tsplit_by_flags(double *xtable, double *labels, int *flags, int nrows, int ncols, double *xtable1, double *labels1, double *xtable2, double *labels2,
					int nlearn) ;
void tsplit_by_flags(double *xtable, double *labels, char *ids, int *flags, int nrows, int ncols, double *xtable1, double *labels1, char *ids1, double *xtable2, 
					 double *labels2, char *ids2, int nlearn) ;

// Extending Data
// Create a new xtable with part of the columns
int get_cols_subset(float *all_xtable, float **xtable, char *all_header, char **header, int nrows, int ncols, int *cols_subset, int subset_size) ;

// Create a new xtable with product featuures
int add_products(float *xtable, char *header, int nrows, int ncols, float **out_xtable, char **out_header, int *prods, int nprods) ;
int add_products(float *xtable, char *header, int nrows, int ncols, float **out_xtable, char **out_header, int *prods1, int nprods1, int *prods2, int nprods2) ;
int add_products(double *xtable, char *header, int nrows, int ncols, double **out_xtable, char **out_header, int *prods, int nprods) ;
int add_products(double *xtable, char *header, int nrows, int ncols, double **out_xtable, char **out_header, int *prods1, int nprods1, int *prods2, int nprods2) ;

// Create a new xtable with square featuures
int add_squares(float *xtable, char *header, int nrows, int ncols, float **out_xtable, char **out_header, int *squares, int nsquares) ;
int add_squares(double *xtable, char *header, int nrows, int ncols, double **out_xtable, char **out_header, int *squares, int nsquares) ;

// Filtering Data
// Create a new xtable with part of the columns
int get_cols_subset(double *all_xtable, double **xtable, char *all_header, char **header, int nrows, int ncols, int *cols_subset, int subset_size) ;
int get_submatrix(double *xtable, double *ytable, char *header, char *ids, int *nrows, int *ncols, int *take_rows, int ntake_rows, int *take_cols, int ntake_cols) ;
int get_cols_submatrix(double *xtable, char *header, int nrows, int *ncols, int *take_cols, int ntake_cols) ;

// Set pos/neg ratio
int tfilter(double *x, double *y, int nsamples, int nftrs, double ratio, double *fx, double *fy, int *fnsamples) ;
int filter (double *x, double *y, int nsamples, int nftrs, double ratio, double *fx, double *fy, int *fnsamples) ;
int filter (double *x, double *y, char *ids, int nsamples, int nftrs, double ratio, double **fx, double **fy, char **fids, int *fnsamples) ;
int filter (double *x, double *y, int nsamples, int nftrs, double ratio, double *fx, double *fy, int *fnsamples, int *finds) ;
int filter (double **x, double **y, int *nsamples, int nftrs, double ratio) ;
int tfilter(double **x, double **y, int *nsamples, int nftrs, double ratio) ;

// Handle missing values - make all true values positive and set missing values to -1.
void adjust_matrix (double *xtable, int *missing, int nsamples, int nftrs) ;

// Query Data
// Correlations
double tget_corr(double *x, double *y, int nsamples, int ind, double missing=-1.0) ;
double get_corr(double *x, double *y, int nsamples, int ind, int nftrs, double missing=-1.0)  ;

#endif