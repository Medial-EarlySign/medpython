// Management of trees and predictions.
struct gbm_tree ;
struct cat_splits ;
#include <string>

using namespace std;

// Prediction
GBMRESULT gbm_pred (double *radX, int rcRows, int rcCols, int cTrees, double rdInitF, const gbm_tree *rTrees, int **rCSplits, int *raiVarType, double *radPredF, double missing_val) ;
GBMRESULT gbm_pred (double *radX, int irow, int rcRows, int rcCols, const gbm_tree *rTrees, int itree, int **rCSplits, int *raiVarType, double *radPredF, double missing_val) ;

GBMRESULT gbm_pred (double *radX, int rcRows, int rcCols, int cTrees, double rdInitF, const gbm_tree *rTrees, int **rCSplits, int *raiVarType, double *radPredF) ;
GBMRESULT gbm_pred (double *radX, int irow, int rcRows, int rcCols, const gbm_tree *rTrees, int itree, int **rCSplits, int *raiVarType, double *radPredF) ;


// structures management memory
void print_gbm_tree(const string& prefix, FILE *file, gbm_tree *trees, int i) ;
void print_gbm_tree(const string& prefix, FILE *file, gbm_tree &tree) ;
void print_gbm_tree(FILE *file, gbm_tree *trees, int i) ;
void print_gbm_tree(FILE *file, gbm_tree &tree) ;
int read_gbm_tree(FILE *file, gbm_tree *trees, int i) ;
int read_gbm_tree(FILE *file, gbm_tree &tree) ;
gbm_tree *allocate_gbm_trees(int ntrees) ;
GBMRESULT set_tree(CGBM *pGBM, int nnodes, VEC_VEC_CATEGORIES&vecSplitCodes, gbm_tree *trees, int itree) ; 
void free_gbm_trees(gbm_tree *trees) ;