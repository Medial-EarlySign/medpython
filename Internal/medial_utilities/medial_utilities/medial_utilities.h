// medial_utilities: cluster-analysis and numerical utilities
#ifndef __MED_UTIL_H__
#define __MED_UTIL_H__ 
#pragma once

#include "zlib/zlib/zlib.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <string>
#include <sstream> 
#include <iostream>
#include <map>
#include <random>
#include <algorithm>

#define MAX_STRING_LEN     1024
#define MAX_FIELD_SIZE 50

#define IDX(i,j,ncol) ((i)*(ncol)*(MAX_STRING_LEN) + (j)*(MAX_STRING_LEN))
#define XIDX(i,j,ncol) ((i)*(ncol) + (j))
#define HIDX(i) ((i))*(MAX_STRING_LEN)
#define SIDX(i,j,ncol,size) ((i)*(ncol)*(size+1) + (j)*(size+1))
#define FIDX(i,j,ncol) ((i)*(ncol)*(MAX_FIELD_SIZE) + (j)*(MAX_FIELD_SIZE))

// Utilities for sorting

struct val_idx {
	int idx ;
	float val ;
} ;

struct dbl_idx {
	int idx ;
	double val ;
} ;

// IO
int read_blob(const char *file_name, unsigned char **data) ;

// QuickSort comparison functions
int dbl_idx_compare (const void *el1, const void* el2) ;
int val_idx_compare (const void *el1, const void* el2) ;
int double_compare (const void *el1, const void* el2) ;
int int_compare (const void *el1, const void* el2) ;
int float_compare (const void *el1, const void* el2) ;

// Utilitities for vector calculations in a matrix
// Find Average of a matrix column
float calc_col_avg(float *table, int col, int nrow, int ncol, float missing) ;
double calc_col_avg(double *table, int col, int nrow, int ncol, double missing) ;

// Find Average of a subset of the rows in a matrix column
float calc_col_sub_avg(float *table, int col, int *inds, int ninds, int ncol, float missing) ;
double calc_col_sub_avg(double *table, int col, int *inds, int ninds, int ncol, double missing) ;

// Find standard deviation of a matrix column
float calc_col_std(float *table, int col, int nrow, int ncol, float missing) ;
double calc_col_std(double *table, int col, int nrow, int ncol, double missing) ;
float calc_col_std(float *table, int col, int nrow, int ncol, float avg, float missing) ;
double calc_col_std(double *table, int col, int nrow, int ncol, double avg, double missing) ;

// Find standard deviation of a subset of the rows in a matrix column
float calc_col_sub_std(float *table, int col, int *inds, int ninds, int ncol, float missing) ;
double calc_col_sub_std(double *table, int col, int *inds, int ninds, int ncol, double missing) ;
float calc_col_sub_std(float *table, int col, int *inds, int ninds, int ncol, float avg, float missing) ;
double calc_col_sub_std(double *table, int col, int *inds, int ninds, int ncol,double avg, double missing) ;

// Find weighted-average of a matrix column
float weighted_calc_col_avg(float *table, int col, float *weights, int nrow, int ncol, float missing) ;
double weighted_calc_col_avg(double *table, int col, double *weights, int nrow, int ncol, double missing) ;

// Find weighted-average of a subset of the lines in a matrix column
float weighted_calc_col_sub_avg(float *table, int col, float *weights, int *inds, int ninds, int ncol, float missing) ;
double weighted_calc_col_sub_avg(double *table, int col, double *weights, int *inds, int ninds, int ncol, double missing) ;

// Find weighted-standard-deviation of matrix column
float weighted_calc_col_std(float *table, int col, float *weights, int nrow, int ncol, float missing) ;
double weighted_calc_col_std(double *table, int col, double *weights, int nrow, int ncol, double missing) ;
float weighted_calc_col_std(float *table, int col, float *weights, int nrow, int ncol, float avg, float missing) ;
double weighted_calc_col_std(double *table, int col, double *weights, int nrow, int ncol, double avg, double missing) ;

// Find weighted-standard-deviation of a subset of the lines in a matrix column
float weighted_calc_col_sub_std(float *table, int col, float *weights, int *inds, int ninds, int ncol, float missing) ;
double weighted_calc_col_sub_std(double *table, int col, double *weights, int *inds, int ninds, int ncol, double missing) ;
float weighted_calc_col_sub_std(float *table, int col, float *weights, int *inds, int ninds, int ncol, float avg, float missing) ;
double weighted_calc_col_sub_std(double *table, int col, double *weights, int *inds, int ninds, int ncol, double avg, double missing) ;

// Find sum of a matrix column
float calc_col_sum(float *table, int col, int nrow, int ncol, float missing) ;
double calc_col_sum(double *table, int col, int nrow, int ncol, double missing) ;

// Shuffle a vector
int shuffle(float *vec, int n, float **new_vec) ;
int shuffle(double *vec, int n, double **new_vec) ;

// Utilities for independent vector calculations
// Pearson correlation of two vectors
float pearson(float *vec1, float *vec2, int n)  ;
double pearson(double *vec1, double *vec2, int n) ;
double pearson(double *vec1, double *vec2, int n, double missing) ;

// Spearman correlation of two vectors
float spearman(float *vec1, float *vec2, int n) ;
double spearman(double *vec1, double *vec2, int n) ;

// Get indices of vec2 elements in vec1
int get_indices(int *vec1, int n1, int *vec2, int n2, int **indices, int *n3) ;

// Get the (tied) order of a vector.
int get_order(float *vec, int n, float **order) ;
int get_order(double *vec, int n, double **order) ;

// Create a random permutation
int *randomize (int nrows) ;

// Mathematical Utilities
// Fast Fourier Transform on data
void four1 (double data[], unsigned long nn, int isign) ;

// Get R and Phi from FFT data
void get_r_phi(double *fft_data,int size, int i, double *r, double *phi) ;

//Gauss Jordan elimination to solve AX = b ;
void gaussj (double **a, int n, double **b, int m) ;

// Clustering
// K-Means
int kmeans (double *x, int nrows, int ncols, int k, int *clusters, double *means) ;

// Identifier closest cluster
int get_closest(double *x,int nrows, int ncols, double *mean) ;

// General utilities
// check if a file exists
bool file_exists (const char *filename) ;
// open files with error message in case of failure
FILE* safe_fopen(const char* filename, const char* mode, bool exit_on_fail = true) ;
// gzopen files with error message in case of failure
gzFile safe_gzopen(const char* filename, const char* mode, bool print_msg = true, bool exit_on_fail = true) ;
// read line from a gzFile (wrapper of gzgets)
char * gzGetLine(gzFile file, std::string& str);
// bring various path formats ("W:/path/to/file" and "/cygdrive/w/path/to/file") to the Windows format as accepted by condor ("\\nas1\Work\path\to\file")
int fix_path(const std::string& in, std::string& out) ;

// time
double get_hour(const char *time, int format = 0) ;
void hours2time (double hours, char *time) ;
double get_day(const char *time, int format = 0) ;
double get_min(char *time, int format = 0) ;
double min2day (double mins) ;

int get_day(int date) ;
int get_date(int days) ;

// Vector Utils
int get_moments (double *v, int n, double *mean, double *sdv, double missing=-1) ;
int get_mean (double *v, int n, double *mean, double missing=-1) ;
int get_sdv (double *v, int n, double mean, double *sdv, double missing=-1) ;
int get_median (double *v, int n, double *median, double missing=-1) ;
int get_quantiles (double *v, int n, double *qs, int nqs, double *vals, double missing=-1) ;

#endif
