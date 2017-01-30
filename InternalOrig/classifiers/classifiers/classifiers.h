// classifiers : All Sorts of Classifiers 

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <tchar.h>
#endif
#include <string.h>

#include "time.h"
#include "gbm/gbm/gbm_utils.h"
#include "medial_utilities/medial_utilities/medial_utilities.h"
#include "data_handler/data_handler/data_handler.h"

#define NITER 200
#define EITER 0.00001
#define RFACTOR 0.99

#define LRN_FILE "R_LrnData"
#define PRD_FILE "R_PrdData"
#define PRED_FILE "R_Preds"
#define OOB_FILE "R_OobData"
#define ERR_FILE "R_Err"
#define R_MDL_FILE "R_Model"
#define LEARN_BATCH "R_Learner.r"
#define IMPORTANCE_BATCH "R_getImportance.r"
#define PRED_BATCH "R_Predictor.r"
#define IMP_BATCH "R_ImpLearner.r"
#define OOB_BATCH "R_Oob.r"
#define STD_FILE "R_StdErr"

#ifdef _WIN32
#define STD_DIR "\\\\nas1\\temp\\R_stderr"
#define R_EXEC "\\\\nas1\\Work\\Applications\\R\\R-latest\\bin\\x64\\R CMD BATCH --silent --no-timing"
#else
#define STD_DIR "/nas1/Temp/R_stderr"
#define R_EXEC "R CMD BATCH --silent --no-timing"
#endif

#define MDL_FILE "R_Model"
#define TREE_FILE "R_Trees"
#define IMP_FILE "R_Importance"
#define MAX_R_ROWS 10000

#define RF_SEED 19740301

// Linear Plus
#define AUCPARTS 3
#define ITERATION_FOR_ONE_CHOICE 30000
#define NUMBER_OF_FEATURES_PER_RUN  20
#define ENS_SIZE 5
#define INTERNAL_TRAIN_PERCENT 1.0
#define ACTIONS_NUM 1
#define LP_MAXLINES 12000000

/*********************************************************/
/************************** KNN **************************/
/*********************************************************/

// Find k nearest neighbors and their squared-distances . test2learn is the index of the sample in question which is ignored.
void find_nbrs(double *test_x, int ind, int k, double *learn_x, int nlearn, int nftrs, int *order, double *ws, int *nbrs, double *dists, int norm, int test2learn) ;

// Find mean distance by random sampling
double get_mean_dist(double *test_x, int ind, double *learn_x, int nlearn, int nftrs, int nrand, double *ws, int norm = 2) ;

// Get score from neihgbors
double nbrs_score(int *nbrs, double *dists, double *y, int k, double mean_dist, int type) ;

// Get Score from weighted linear model on neighbours
int nbrs_ls(double *testx, int ind, int *nbrs, double *dists, double *x, double *y, int k, int nftrs, double mean_dist, double *localx, double *localy, double *localw,
			double *localb, double *localr, double *pred) ;

// Get orders of features by stdv x weight
int order_ftrs(double *x, int nsamples, int nftrs, double *weights, int *order) ;

// Predict a single instance using KNN
int knn_predict (double *test_x, int ind, double *learn_x, double *learn_y, int nlearn, int nftrs, double *weights, int *order, int *nbrs, double *dists, int k,
					int type, double *nbrs_x, double *nbrs_y, double *nbrs_w, double *nbrs_b, double *nbrs_r, double *pred) ;

// Predict a single instance using KNN. test2learn is the index of the sample in question which is ignored.
int knn_predict (double *test_x, int ind, double *learn_x, double *learn_y, int nlearn, int nftrs, double *weights, int *order, int *nbrs, double *dists, int k,
					int type, double *nbrs_x, double *nbrs_y, double *nbrs_w, double *nbrs_b, double *nbrs_r, int test2learn, double *pred) ;

// Predict using KNN : Types - 1=Eucleadean dist, nbr-w=dist/mean ; 2=Eucleadean dist, nbr-w=1/dist; 3=Eucleadean dist + Weighted LS; 4=L1 dist + Weighted LS;
int knn_predict(double *test_x, int ntest, double *learn_x, double *learn_y, int nlearn, int nftrs, double *weights, double *preds, int k, int type=1) ;

// Predict using KNN. finds = indices of test samples in learn samples.
int knn_predict(double *test_x, int ntest, double *learn_x, double *learn_y, int nlearn, int nftrs, double *weights, double *preds, int k, int *finds, int type=1) ;

/*********************************************************/
/******************* Linear Regression *******************/
/*********************************************************/

// Ax = y ; X(nsamples*nftrs) ; niter Iterations up to relative error eiter, Ridge-factor(s) = rfactor(s) ; w(s) = weights ; ignore = columns to ignore ; corrs = used
//    to force correct sign to coeffieicnts. Output = b,err 
// sumxx = pre-calculated sum(X^2) for efficiency
int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double rfactor, double *b, double *err) ;
int lm_positive (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double rfactor, double *b, double *err) ;
int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double rfactor, double *ws, double *b, double *err) ;
int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *ws, double *b, double *err) ;
int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore) ;
int lm (double *x, double *_y, double *w, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore, double *corrs) ;
int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore, double *corrs) ;
int lm (double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore, double *corrs, double *sumxx) ;


// wrapper for first lm func, doing internally all x,y normalizations and returning a b vec with affine coefficient in b[nftrs] (b[0...nftrs-1] are the features coefficients)
int lm_wrapper(double *x, double *_y, int nsamples, int nftrs, int niter, double eiter , double rfactor, double *b, double *err);
int lm_wrapper(float *x, float *y, int nsamples, int nftrs, int niter, double eiter , double rfactor, float *b, double *err);

// lm that CHANGES y
int label_changing_lm (double *x, double *y, int nsamples, int nftrs, int niter, double eiter , double *rfactors, double *b, double *err, int *ignore, double *corrs) ;

// Predict
void lm_predict(double *x, double *avg, int nsamples, int nftrs, double *b, double *preds) ;
void lm_predict(double *x, int nsamples, int nftrs, double *b, double *preds) ;
void lm_predict(double *x, int nsamples, int nftrs, double *b, double norm, double *preds) ;

// Linear regression to predict missing features
int features_lm(double *x, int nsamples , int nftrs,int niter ,double eiter, double *rfactors, double **ftr_bs) ;

// Complete features in (non-transposed) matrix according to linear regression params
void complete_ftrs(double *x, int *missing, int nsamples, int nftrs, double *complete_x, double **ftr_bs)  ;

/*********************************************************/
/********************* Random Forest *********************/
/*********************************************************/

typedef struct {
	int left,right,status,svar,pred ;
	double sval ;
} rf_node ;

typedef struct {
	rf_node *nodes ;
	int nnodes ;
} rf_tree ;

typedef struct {
	rf_tree *trees ;
	int ntrees ;
} random_forest ;

// Single tree prediction
int tree_pred(rf_tree tr, double *values) ;

// Forest prediction
void forest_pred(random_forest *forest, double *values, int *counts, int nclass) ;

// Read from file
int read_forest(char *file_name, random_forest *forest) ;
int write_forest(random_forest *forest, char *file_name) ;
void print_forest(random_forest *forest, FILE *fp) ;
int read_text_forest(char *file_name, random_forest *forest) ;

// predict (default ftree = TREE_FILE)
int rf_predict(double *x, double *preds, int nrows, int ncols, int trans_flag, random_forest *rf) ;

int rf_predict(double *tx, double *preds, int nrows, int ncols, random_forest *rf)  ;
int rf_predict(double *tx, double *preds, int nrows, int ncols, unsigned char *model) ;
int rf_predict(double *tx, double *preds, int nrows, int ncols) ;
int rf_predict(double *tx, double *preds, int nrows, int ncols, char *ftree) ;
int rf_predict_from_text_forest(double *tx, double *preds, int nrows, int ncols, char *ftree) ;
int rf_predict_from_text_forest(double *tx, double *preds, int nrows, int ncols) ;
int get_prediction_matrix(double *tx, int nrows, int ncols, double **preds, int *ntrees) ;

// R-Envelope: learn
int R_learn(double *tx, double *y, int nrows, int ncols) ;

// R-Envelope: check importance
int R_get_classifiction_importance(double *tx, double *y, int nrows, int ncols, int *ignore, double *importance) ;

// R-Envelope: learn classification
int R_learn_classification(double *tx, double *y, int nrows, int ncols, char *type, int eq_size = 1) ;
int R_learn_classification(double *tx, double *y, int nrows, int ncols, int *ignore, char *type, int eq_size = 1) ;
int R_learn_classification(double *tx, double *y, int nrows, int ncols, char *type, char *rf_file, int eq_size = 1) ;

// R-Envelope: predict
int R_predict(double *tx, double *preds, int nrows, int ncols) ;

// R-Envelope : Get out-of-bag error
int R_get_oob_err (double *tx, double *y, int nrows, int ncols, int *ignore, int ntrees, double *err) ;

// Serialization
size_t get_rf_size(random_forest *forest) ;
int rf_serialize(random_forest *forest, unsigned char *rf_data) ;
int rf_deserialize(unsigned char *rf_data, random_forest *forest) ;

// Cleaning
void clear_random_forest(random_forest *forest) ;

/*********************************************************/
/********************* QRF Forest ************************/
/*********************************************************/
int qrf_predict(double *tx, double *preds, int nrows, int ncols, unsigned char *model) ;
int qrf_predict(double *tx, double *preds, int nrows, int ncols, unsigned char *model, int nthreads) ;

/*********************************************************/
/********************** R Predictor **********************/
/*********************************************************/

// R-Envelope: predict class (default mdl_file = MDL_FILE). type = rf/lm/....
int R_predict_class(double *tx, double *preds, int nrows, int ncols, char *type) ;
int R_predict_class(double *tx, double *preds, int nrows, int ncols, int *ignore, char *type) ;
int R_predict_class(double *tx, double *preds, int nrows, int ncols, char *type, char *mdl_file) ;
int R_predict_class_with_NA(double *tx, double *preds, int nrows, int ncols, char *type, double missing) ;
int R_predict_class_with_NA(double *tx, double *preds, int nrows, int ncols, int *ignore, char *type, double missing) ;


/*********************************************************/
/********************** Barak-Trees **********************/
/*********************************************************/

#define NLOOP 3
#define BTREE_IT  1500 // 2000

struct barak_trees {
	char  **stree_leaf[NLOOP] ;
	float ***stree_val_final[NLOOP] ;
	int   **stree_x[NLOOP] ;
	char  **stree_sign[NLOOP] ;
	float **stree_saf[NLOOP] ;
	int   **stree_next[NLOOP] ;
	float *avrx ;
	float *lm_avgx[NLOOP] ;
	float lm_avgy[NLOOP] ;
	float *streeb[NLOOP] ;
	int n_trees ;
} ;

// Initialize
int initialize_brk_trees(barak_trees *bt, int ntrees, int ncols) ;

// Free
void free_brk_trees(barak_trees *bt) ;

// Learn
int learn_brk(int n_samples ,int n_feature ,float *p_trainx ,float *_trainy, int data_use ,int n_leaf ,float step_gap ,barak_trees  *my_tree_struct, 
			  float missing_val = -65536.0) ;

// Predict
int brk_predict (int n_sample ,int n_feature ,  float *p_testx , barak_trees  *test_tree_struc  , double *test_pred , float missing_val = -65536.0) ;

// Write to file
int make_trees_file(barak_trees  *my_tree_struc , char *fname) ;

// Read from file
int read_trees_file(barak_trees  *my_tree_struc , char *fname) ;

/*********************************************************/
/*********************** Nir-Trees ***********************/
/*********************************************************/

// x = predictors matrix (transposed)
// y = response vector
// (nrows,ncols) = dimensions of x
// training = flags for internal-training set

// Learn
int learn_nir_trees(double *x, double *y, int nrows, int ncols, char *file_name) ;

// Learn : Internal step
int make_trees(double *x, double *y, int nrows, int ncols, int *training, int gloop) ;

// Predict
int nir_trees_predict(double *x, double *preds, int nrows, int ncols, char *file_name) ;

// Predict : Internal step
int runontest(double *x, int nrows, int ncols, int gloop) ;
void finalresult(double *x, int nrows, int ncols, double *predictions) ;

// Print to file
int print_nir_trees(char *file_name) ;

// Read from file
int read_nir_trees(char *file_name) ;

/*********************************************************/
/******************** Nir-GBM Enemble ********************/
/*********************************************************/

#define GBM_ENS_SIZE 100 

typedef struct  {
	int ens_size ;
    double *alpha;
	full_gbm_learn_info_t *full_info_tree;

} full_gbm_ens_learn_info_t ;

int gbm_predict_ens(double *x, int nrows, int ncols, full_gbm_ens_learn_info_t *full_info_ens, double *preds);
int get_gbm_predictor_ens(double *x1, double *y1, int train_size, int nftrs, int ens_size, full_gbm_ens_learn_info_t *full_info_ens, gbm_parameters *gbm_params);
int get_gbm_predictor_ens(double *x1, double *y1, double *w1, int train_size, int nftrs, int ens_size, full_gbm_ens_learn_info_t *full_info_ens, gbm_parameters *gbm_params);
void clear_gbm_predictor_ens(full_gbm_ens_learn_info_t *full_info_ens) ;
int read_full_gbm_ens_info(full_gbm_ens_learn_info_t *ens_info, char *fname) ;
int write_full_gbm_ens_info(full_gbm_ens_learn_info_t *ens_info, char *fname) ;
void write_full_gbm_ens_info(full_gbm_ens_learn_info_t *ens_info, FILE *fp) ;

// Serialization
size_t get_gbm_ens_info_size(full_gbm_ens_learn_info_t *ens_info) ;
int gbm_ens_serialize(full_gbm_ens_learn_info_t *ens_info, unsigned char *ens_data) ;
int gbm_ens_deserialize(unsigned char *ens_data, full_gbm_ens_learn_info_t *ens_info) ;

/*********************************************************/
/********************** Linear Plus **********************/
/*********************************************************/

typedef struct 
{
	int l_rand_ftr;
	int l_rand_ftrb;
	int l_rand_ftrc;
	int l_rand_ftrd;
	int l_rand_act;
	int l_rand_dir;
	int l_rand_ab;
	int l_rand_abc;
	int l_rand_abcd;
	double l_first_saf;
	double l_second_saf;
	double l_first_safb;
	double l_second_safb;
	double l_first_safc;
	double l_second_safc;
	double l_first_safd;
	double l_second_safd;
	double l_rand_size;
}full_linear_plus_info_line_t;

typedef struct
{	 	
	int nftrs;
	double *avrgs;
	double *stdevs;
	int linesno;
	full_linear_plus_info_line_t *lines; 
} full_linear_plus_info_t ;

// Learn a LinearPlus predictor
void get_linear_plus(double *x1, double *y1, int train_size, int nftrs, full_linear_plus_info_t *full_info_lin_plus );
// Learn a LinearPlus predictor and write to a file (textually)
int get_linear_plus_f(double *x, double *y, double *w, int nrows, int ncols, char *fname);
// Learn a LinearPlus predictor and serialize
int learn_linear_plus_predictor(double *x, double *y, int nrows, int ncols, unsigned char **model) ;
// Write LinearPlus object to a textual file
int write_full_linear_plus(full_linear_plus_info_t *linear_info, char *fname) ;
void write_full_linear_plus(full_linear_plus_info_t *linear_info, FILE *fname)  ;
// Read LinearPlus object from a textual file
int read_full_linear_plus(full_linear_plus_info_t *linear_info, char *fname);
// Predict using a LinearPlus object
int predict_linear_plus(double *x, int nrows, int ncols, full_linear_plus_info_t *linear_info, double *preds);
// Predict using a serialized LinearPlus objet
int linear_plus_predict(double *x, double *preds, int nrows, int ncols, unsigned char *model) ;
// Predict using a LinearPlus object from a textual file
int linear_plus_predict(double *x, double *preds, int nrows, int ncols, char *fname);
// Linear Plus cross validation - learn on trainig-set and predict on test-set.
int get_linear_plus_predictions(double *x1, double *y1, int nrows1, double *x2, int nrows2, int ncols, double *preds) ;
// Serialization
size_t get_linear_plus_size(full_linear_plus_info_t *linear_info) ;
int linear_plus_serialize(full_linear_plus_info_t *linear_info, unsigned char *linear_data) ;
int linear_plus_deserialize(unsigned char *linear_data, full_linear_plus_info_t *linear_info) ;
// Cleaning
void clear_linear_plus(full_linear_plus_info_t *linear_info) ;
