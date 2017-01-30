#pragma once

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#define EXTERNAL_HEADER

#include<vector>
#include "dataset.h"
#include "distribution.h"
#include "bernoulli.h"
#include "adaboost.h"
#include "poisson.h"
#include "gaussian.h"
//#include "coxph.h"
#include "laplace.h"
#include "quantile.h"
#include "gbm_engine.h"
#include "gbm_predict.h"

// Management of trees and predictions.
struct gbm_tree ;
struct cat_splits ;

enum GBM_LossFunctions {
	GBM_Loss_AdaBoost,
	GBM_Loss_Bernulli,
	GBM_Loss_Gaussian,
	GBM_Loss_Laplace,
	GBM_Loss_Poisson,
	GBM_Loss_Quantile,
	GBM_Loss_Last
} ;


typedef struct {
	double shrinkage ;
	double bag_p ;
	int ntrees ;
	int depth  ;
	int min_obs_in_node ;
	bool take_all_pos  ;
} gbm_parameters ;

typedef struct {	 
	int ntrees; 
	gbm_tree *trees;
	double initF;
	int **rcSplits;
	int *rcSplitSizes;     
	int ncSplits;
} full_gbm_learn_info_t ;

#define MIN_OBS_IN_NODE 500

// Learn 
int *get_order(double *x, int nrows, int ncols) ;
int get_gbm_predictor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, double bag_p, 
					  bool take_all_pos, int ntrees, int depth,
						int min_obs_in_node, gbm_tree *trees, double *initF, int ***rcSplits, 
						int **rcSplitSizes, int *ncSplits, GBM_LossFunctions lossFunction = GBM_Loss_AdaBoost, double alpha = 0.25) ;
int get_gbm_regressor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, double bag_p, int ntrees, int depth, int min_obs_in_node, 
					  gbm_tree *trees, double *initF, int ***rcSplits, int **rcSplitSizes, int *ncSplits) ;
int get_gbm_predictor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, double bag_p, bool take_all_pos, int ntrees, int depth, int min_obs_in_node, full_gbm_learn_info_t *gbm_info,
					  GBM_LossFunctions lossFunction = GBM_Loss_AdaBoost, double alpha = 0.25) ;
int get_gbm_regressor(double *x, double *y, double *w, int nrows, int ncols, double shrinkage, double bag_p, int ntrees, int depth, int min_obs_in_node, full_gbm_learn_info_t *gbm_info) ;

// Predict
int gbm_predict(double *x, int nrows, int ncols, gbm_tree *trees, int ntrees, double init_f, int **rcSplits, double *preds) ;
int gbm_predict(double *x, int nrows, int ncols, int ntrees, full_gbm_learn_info_t *gbm_info, double *preds) ;

// Clear
void clear_gbm_info(full_gbm_learn_info_t *gbm_info) ;

// Read from File
int read_full_gbm_info(full_gbm_learn_info_t *gbm_info, char *fname) ;
int read_full_gbm_info(full_gbm_learn_info_t *gbm_info, FILE *fp) ;

// Write to File
int write_full_gbm_info(full_gbm_learn_info_t *gbm_info, char *fname) ;
void write_full_gbm_info(full_gbm_learn_info_t *gbm_info, FILE *fp) ;
void write_full_gbm_info(const string& prefix, full_gbm_learn_info_t *gbm_info, FILE *fp) ;

// Serialization
size_t get_gbm_info_size(full_gbm_learn_info_t *gbm_info) ;
int gbm_serialize(full_gbm_learn_info_t *gbm_info, unsigned char *blob) ;
int gbm_deserialize(unsigned char *blob, full_gbm_learn_info_t *gbm_info) ;