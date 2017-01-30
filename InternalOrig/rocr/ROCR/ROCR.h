// ROCR utilities : Analyze regression/classification results

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Prediction Values and Comparison
struct pred_val {
	double pred ;
	int ind ;
} ;
int compare_pred_vals( const void *arg1, const void *arg2 ) ;


// Generate graph with x = xmode = rec/sens/fpr/rpp/scr ; and y = ymode = tpr/ppv/prec/lift
int generate_graph(float *labels, float *predictions, int nsamples, char *xmode, char *ymode, int resolution, float **x, float **y) ;
int generate_graph(double *labels, double *predictions, int nsamples, char *xmode, char *ymode, int resolution, double **x, double **y) ;

// Find optimal confusion table
int get_optimal_confusion(float *labels, float *predictions, int nsamples, int confusion[2][2]) ;
int get_optimal_confusion(double *labels, double *predictions, int nsamples, int confusion[2][2]) ;

// Generate a performance histogram
int generate_hist(float *labels, float *predictions, int nsamples, char *mode, int resolution, float **x, float **y, int **count) ;
int generate_hist(double *labels, double *predictions, int nsamples, char *mode, int resolution, double **x, double **y, int **count) ;

// Calculate AUC. remove-part indicates there are samples with y>1.0 that should be removed
int get_auc(float *labels, float *predictions, int nsamples, int resolution, float *auc) ;
int get_auc(double *ilabels, double *predictions, int nsamples, int resolution, double *auc, int remove_part = 0) ;

// Calculate AUC using a subset of MARKED samples
int get_auc(float *labels, float *predictions, int tot_nsamples, int resolution, int *marked, float *auc) ;
int get_auc(double *labels, double *predictions, int tot_nsamples, int resolution, int *marked, double *auc)  ;