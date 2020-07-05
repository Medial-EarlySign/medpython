//
// Quantized Random Forest - an attempt to ultra fast RF version
//
// Currently - only binary categories (y values - 0/1)
//

#ifndef __QRF__H__
#define __QRF__H__

#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <vector>
#include <map>
#include <random>

#define QRF_MODE_QRF		0

enum QRF_TreeType {
	QRF_BINARY_TREE = 0,
	QRF_REGRESSION_TREE = 1,
	QRF_CATEGORICAL_CHI2_TREE = 2,
	QRF_CATEGORICAL_ENTROPY_TREE = 3,
	QRF_MULTILABEL_ENTROPY_TREE = 4,
	QRF_LAST = 4
};

#define REGRESSION_SPLIT_IMPROVE 0.999999

//#define REGRESSION_SPLIT_IMPROVE 0.8

#define MIN_SPLIT_NODE_SIZE 50
#define MIN_SPLIT_SPREAD	0.1

#define PREDS_REGRESSION_AVG			0
#define PREDS_REGRESSION_WEIGHTED_AVG	1
#define PREDS_REGRESSION_QUANTILE		2
#define PREDS_REGRESSION_WEIGHTED_QUANTILE		3
#define PREDS_REGRESSION_SAMPLE			4
#define PROBS_CATEG_MAJORITY_AVG		0
#define PROBS_CATEG_AVG_PROBS			1
#define PROBS_CATEG_AVG_COUNTS			2
#define PREDS_CATEG_MAJORITY_AVG		4
#define PREDS_CATEG_AVG_PROBS			5
#define PREDS_CATEG_AVG_COUNTS			6

#define COLLECTED_VALUES_SIZE	200000

using namespace std;

struct ValInd {
	float val;
	int idx;
};

class QRF_Node {

public:
	// next used for both learn & predict
	short split_feat;
	float split_val;
	int counts[2];
	char is_leaf;
	int r_ind;
	int l_ind;
	float pred; // #1 / (#1 + #0)
	float s2;   // avg sq err for regression
	int depth;

	// next used only for learn
	short split_q_idx;
	int from_sample;
	int to_sample;
	int size() { return (to_sample - from_sample + 1); }


	size_t estimated_size() { return 12; }
};

class QRF_Tree {
public:
	// next used for both learn & predict
	vector<QRF_Node> nodes;
	int n_nodes;
	int max_nodes;

	// next used only for learn
	vector<int> sample_ids;
	map<short, short> feat_chosen;	// used in random feature choices within a node while building the tree.
	vector<int> hist[2];			// histogram for binary case (here to  avoid reallocation)
	vector<ValInd> qy;				// holding and sorting y values in a node when searching a split
	vector<int> inds;

	vector<double> histr_sum;
	vector<int>	histr_num;


	size_t estimated_size() { size_t s = histr_num.size() + histr_num.size() + inds.size() + qy.size() * 2 + hist[0].size() + hist[1].size() + feat_chosen.size() * 4 + sample_ids.size(); for (auto &n : nodes) s += n.estimated_size(); return s; }

	mt19937 rand_gen;




	void print(FILE *f);
	void init_rand_state();
};

struct OOB_Result {
	float sum;
	int n_tests;
	float mean;
	float std;
	vector<int> cnts;
	vector<float> probs;
};

class QRF_ResNode : public SerializableObject {
	// this class is designed to hold the minimal needed information in a node in order to use a tree as a predictor
public:
	int mode;
	short ifeat;
	float split_val;
	int is_leaf;
	int left;
	int right;
	float pred;
	int n_size;
	vector<int> counts;		// none for regression, 2 for binary, ncateg (defined in QRF_Forest) for categorized case
	vector<int> values; // Counting of values in learning set in regression learning
	vector<pair<int, unsigned int>> value_counts; // Counting of values in learning set in regression learning when working in sparse mode
	int tot_n_values;
	int majority;			// for categories cases


	size_t estimated_size() { return 10 + counts.size() + values.size() + 2*value_counts.size(); }

	void get_scores(int mode, int get_counts_flag, int n_categ, vector<float> &scores) const;

	ADD_CLASS_NAME(QRF_ResNode)
		ADD_SERIALIZATION_FUNCS(n_size, mode, ifeat, split_val, is_leaf, left, right, pred, counts, values, value_counts, tot_n_values, majority)

};


class QRF_ResTree : public SerializableObject {
public:
	vector<QRF_ResNode> qnodes;
	
	size_t estimated_size() { size_t s = 0; for (auto &q : qnodes) s += q.estimated_size(); return s; }

	ADD_CLASS_NAME(QRF_ResTree)
		ADD_SERIALIZATION_FUNCS(qnodes)
};

class QuantizedRF {

public:

	int tree_mode; // QRF_BINARY_TREE or QRF_REGRESSION_TREE or QRF_CATEGORICAL_CHI2_TREE
	int NSamples, NFeat, MaxQ;
	vector<short> max_q;
	vector<char> y;
	vector<int> ids[2];
	vector<vector<float>> quant_values;
	vector<vector<short>> q_data;
	vector<float> yr; // y for regression trees
	vector<vector<int>> yr_multilabel; // y for multilabel splitting
	vector<float> w; // weights

	vector<double> log_table;

	vector<OOB_Result> cv;

	vector<char> test_s;

	vector<int> groups;

	int n_groups;


	int init_all(float *X, int *Y, float *Yr, const float *W, int nfeat, int nsamples, int maxq);

	void init_groups(vector<int> &groups_in);
	// binary problem related
	int init(float *X, int *Y, int nfeat, int nsamples, int maxq);
	int get_Tree(int *sampsize, int ntry, QRF_Tree &tree);
	int collect_Tree_oob_scores(float *x, int nfeat, QRF_ResTree &resTree, vector<int>& sample_ids);
	void complete_oob_cv();
	void score_tree_by_index(float *x, int nfeat, QRF_ResTree &tree, int id, float& score, int& majority, vector<int> &counts);


	double get_cross_validation_auc();

	// regresssion problem
	int init_regression(float *X, float *Y, const float *W, int nfeat, int nsamples, int maxq);
	int get_regression_Tree(int *sampsize, int ntry, QRF_Tree &tree);

	// categorized, split by chi square problem
	int n_categ;

	// stopping criteria related
	int min_split_node_size;
	float min_split_spread;


	int n_called;

	void clear();

	bool take_all_samples;

	int max_depth;

private:
	//int quantize_feature(vector<ValInd> &x, int nsamples, int maxq, vector<float> &quant_val, vector<short> &qd);
	int quantize_no_loss(vector<ValInd> &vals, int nsamples, int maxq, vector<float> &quant_val, vector<short> &qd);

	// binary problem
	int find_best_split(QRF_Tree &tree, int node, int ntry);
	int split_node(QRF_Tree &tree, int node);

	// regression problem
	int find_best_regression_split(QRF_Tree &tree, int node, int ntry);
	int split_regression_node(QRF_Tree &tree, int node);

	// categorized, split by chi square problem
	int find_best_categories_chi2_split(QRF_Tree &tree, int node, int ntry);
	int find_best_categories_entropy_split(QRF_Tree &tree, int node, int ntry);
	int find_best_categories_entropy_split_multilabel(QRF_Tree &tree, int node, int ntry);

};

struct qrf_scoring_thread_params {
	const float *x;
	const vector<QRF_ResTree> *trees;
	int from;
	int to;
	int nfeat;
	int nsamples;
	float *res;

	int serial;
	int state;

	int mode;
	int n_categ;
	vector<float> *quantiles;
	const vector<float> *sorted_values;
	bool sparse_values;

	int get_counts; // for CATEGORICAL runs there's such an option, 0 - don't get, 1 - sum counts 2 - sum probs
					//	thread th_handle;
};

//============================================================================================================================
class QRF_Forest : public SerializableObject {
public:
	int mode;				// 0 - regular mode (only QRF) 
	int min_node_size;		// for stop criteria - min size of a node above which we keep on splitting (for categorical or regression)
	float min_spread;		// for stop criteria - min spread of a node above which we keep on splitting (for regression)
	int n_categ;			// categories in y must be 0,1,...,n_categ-1 in this case
	int get_counts_flag;	// how to get results
	int get_only_this_categ; // in case of categorical model, return only probs for this category in predictions (def is -1: return all)

	bool keep_all_values; // Keep all values in each node in regression mode
	bool sparse_values; // Keep all values in sparse mode (as histogram)
	vector<float> quantiles; // Quantiles to predict in quantile-regression mode
	vector<float> sorted_values; // sorted prediction values actually appearing in learning set in regression mode

	vector<int> groups;		// if given should be at length nsamples. It maps each sample to a group (for example samples from the same patient at different times)
							// this is important when we randomize elements to each tree and when we test oob.

	vector<QRF_ResTree> qtrees;		// actual collection of trees built or deserialized

	// out of bag related arrays
	int collect_oob;
	vector<vector<float> > oob_scores; // used to keep oob scores which we use in some cases - if regression : (mean,std) ; if classification : vector of probs

	// constructor: 
	QRF_Forest() {
		qtrees.clear(); mode = 0; collect_oob = 0; oob_scores.clear(); keep_all_values = false; sparse_values = true; nthreads = 1; min_node_size = MIN_SPLIT_NODE_SIZE; get_only_this_categ = -1;
		min_spread = 0; n_categ = -1; get_counts_flag = 0; take_all_samples = false; max_depth = 0;
	}

	int nthreads; // number of threads to use in building a forest, and in scoring it.

	//full samples to take - similar algorithm to knn where k=min_node. use with n_trees=1
	bool take_all_samples;

	//max depth  features to take in branch of tree
	int max_depth;

	// Learn Methods :
	//-----------------

	// binary problem :: 
	// y must be 0/1 values
	// sampsize is either NULL (all taken for bagging) or of size 2 (sampsize[0], sampsize[1] - bagging with these sizes)
	int get_forest(double *x, double *y, int nfeat, int nsamples, int *sampsize, int ntry, int ntrees, int maxq);
	int get_forest(float *x, int *y, int nfeat, int nsamples, int *sampsize, int ntry, int ntrees, int maxq);

	// regression problem
	// sampsize - 0: take all, other: bag using given size.
	int get_forest_regression_trees(float *x, float *y, int nfeat, int nsamples, int sampsize, int ntry, int ntrees, int maxq, int min_node, float spread);

	// categorical chi2/entropy problem - Splitting method: QRF_CATEGORICAL_ENTROPY_TREE/QRF_CATEGORICAL_CHI2_TREE
	// y must be in 0 ... (ncateg-1) values.
	// sampsize : if NULL all bagging size is the number of samples, otherwise sampsize for each category is expected (sampsize[0]....sampsize[ncateg-1])
	int get_forest_categorical(float *x, float *y, const float *w, int nfeat, int nsamples, int *sampsize, int ntry, int ntrees, int maxq, int min_node, int ncateg, int splitting_method);

	// Predict Methods :
	//------------------

	int score_samples(float *x, int nfeat, int nsamples, float *&res) const; // same as next one but get_counts = 0;
	// res should be allocated outside.
	// for classification nsamples x n_categ. get_counts : 0 - average on majority vote of each tree ; 1 - average on probabilities in each tree ; 2 (or otherwise) - average on counts in each tree
	//													   4 - prediction (based on majority)
	// for regression nsamples x 1. get_counts : 0 - average prediction ; 1 (or otherwise) - weighted average.
	int score_samples(float *x_in, int nfeat, int nsamples, float *&res, int get_counts) const;

	int score_samples_t(double *x, int nfeat, int nsamples, double *&res); // tailored for predictor formats

	void get_single_score_fast(qrf_scoring_thread_params &params, const vector<float> &x, vector<float> &preds) const;

	// OOB Methods:
	//--------------
	//int collect_Tree_oob_scores_threaded(float *x, int nfeat, QRF_ResTree &resTree, vector<int>& sample_ids);


	// serialization
	ADD_CLASS_NAME(QRF_Forest)
		ADD_SERIALIZATION_FUNCS(qtrees, mode, min_node_size, min_spread, n_categ, get_counts_flag, get_only_this_categ, keep_all_values, sparse_values, quantiles, sorted_values, nthreads, take_all_samples, max_depth)

		// IO Methods :
		//-------------
		void write(FILE *fp);

	// Variable Importance:
	//---------------------

	void variableImportance(vector<pair<short, double> >& rankedFeatures, unsigned int nFeatures);

private:
	// learn for all modes
	int get_forest_trees_all_modes(float *x, void *y, const float *w, int nfeat, int nsamples, int *sampsize, int ntry, int ntrees, int maxq, int mode);

	// scoring with threads for all modes
	void score_with_threads(float *x, int nfeat, int nsamples, float *res) const;

	// transferring from inner algo data structure to exposed on
	int transfer_to_forest(vector<QRF_Tree> &trees, QuantizedRF &qrf, int mode);
	int transfer_tree_to_res_tree(QuantizedRF &qrf, QRF_Tree &tree, QRF_ResTree &qt, int mode, map<float, int> &all_values);
	int init_keep_all_values(QuantizedRF &qrf, int mode, map<float, int> &all_values);
};

MEDSERIALIZE_SUPPORT(QRF_ResNode)
MEDSERIALIZE_SUPPORT(QRF_ResTree)
MEDSERIALIZE_SUPPORT(QRF_Forest)

#endif