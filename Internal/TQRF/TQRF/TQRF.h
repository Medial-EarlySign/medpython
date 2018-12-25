#ifndef __TQRF_H__
#define __TQRF_H__
//
// TQRF
//

#include <Logger/Logger/Logger.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <MedMat/MedMat/MedMat.h>
#include <MedTime/MedTime/MedTime.h>
#include <queue>

using namespace std;

enum TQRF_TreeTypes {
	TQRF_TREE_ENTROPY = 0,
	TQRF_TREE_REGRESSION = 1,
	TQRF_TREE_LIKELIHOOD = 2,
	TQRF_TREE_WEIGHTED_LIKELIHOOD = 3,
	TQRF_TREE_DEV = 4, // free place to use when developing new score ideas
	TQRF_TREE_UNDEFINED = 5
};

enum TQRF_Missing_Value_Method {
	TQRF_MISSING_VALUE_MEAN = 0,
	TQRF_MISSING_VALUE_MEDIAN = 1,
	TQRF_MISSING_VALUE_LARGER_NODE = 2,
	TQRF_MISSING_VALUE_LEFT = 3,
	TQRF_MISSING_VALUE_RIGHT = 4,
	TQRF_MISSING_VALUE_RAND_ALL = 5,
	TQRF_MISSING_VALUE_RAND_EACH_SAMPLE = 6,
	TQRF_MISSING_VALUE_NOTHING = 7
};


enum TQRF_Node_Working_State {
	TQRF_Node_State_Initiated = 0,
	TQRF_Node_State_In_Progress = 1,
	TQRF_Node_State_Done = 2
};

enum TQRF_Missing_Direction {
	TQRF_MISSING_DIRECTION_LEFT = 0,
	TQRF_MISSING_DIRECTION_RIGHT = 1,
	TQRF_MISSING_DIRECTION_RAND_EACH_SAMPLE = 2
};

#define TQRF_MAX_TIME_SLICE			10000000
#define MIN_ELEMENTS_IN_TIME_SLICE	100
#define UNSET_BETA	((float)-1e-10)


//==========================================================================================================================
class TQRF_Params : public SerializableObject {
public:

	//========================================================================================================================
	// params list
	//========================================================================================================================
	string init_string = "";			/// sometimes it helps to keep it for debugging
	string samples_time_unit = "Date";
	int samples_time_unit_i; /// calculated upon init

	int ncateg = 2;			/// number of categories (1 for regression)

	// time slices
	string time_slice_unit = "Days";
	int time_slice_unit_i; /// calculated upon init
	int time_slice_size = -1;			/// the size of the basic time slice, -1: is like infinity: a single time slice like a regular QRF
	int n_time_slices = 1;				/// if time_slices vector is not given, one will be created using time_slice_size and this parameter. 
	vector<int> time_slices = {};		/// if not empty: defines the borders of all the time lines. Enables a very flexible time slicing strategy
	
	vector<float> time_slices_wgts ={}; /// default is all 1.0 , but can be assigned by the user, will be used to weight the scores from different time windows
	int censor_cases = 0;				/// when calclating the time slices distributions we have an option to NOT count the preciding 0's of non 0 cases.

	// quantization
	int max_q = 200;					/// maximal quantization
	//int max_q_sample = 100000;			/// the max number of values to use when deciding q limits

	// trees and stopping criteria
	string tree_type = "";				/// options: regression, entropy, logrank
	int tree_type_i = -1;				/// tree type code : calulated from tree_type the string
	int	ntrees = 50;					/// number of trees to learn
	int max_depth = 100;				/// maximal depth of tree
	int min_node_last_slice = 10;		/// stopping criteria : minimal number of samples in a node in the last time slice
	int min_node = 10;					/// stopping criteria : minimal number of samples in a node in the first time slice
	float random_split_prob = 0;		/// at this probability we will split a node in a random manner, in order to add noise to the tree.

	// feature sampling
	int ntry = -1;						/// -1: use the ntry_prob rule, > 0 : choose this number of features.
	float ntry_prob = (float)0.1;		/// choose ntry_prob * nfeat features each time
	int nsplits = -1;					/// -1: check all splits for each feature , then split the max, > 0: choose this number of split points at random and choose best

	// speedup by subsample control
	int max_node_test_samples = 50000;	/// when a node is bigger than this number : choose this number of random samples to make decisions

	// bagging control
	int single_sample_per_pid = 1;		/// when bagging select a single sample per pid (which in itself can be repeated)
	int bag_with_repeats = 1;			/// weather to bag with repeats or not
	float bag_prob = (float)0.5;		/// random choice of samples for each tree prob
	float bag_ratio = -1;				/// control ratio of #0 : #NonZero of labels, if < -1 , leave as is.
	float bag_feat = (float)1.0;		/// proportion of random features chosen for each tree
	int qpoints_per_split = 0;			/// if > 0 : will only choose this random number of points to test split points at, otherwise will test all of them

	// categorial featues
	int nvals_for_categorial = 0;			/// features with number of different values below nvals_for_categ will be assumed categorial
	vector<string> categorial_str;		/// all features containing one of the strings defined here in their name will be assumed categorial
	vector<string> categorial_tags;		/// all features containing these tags will be assumed categorial
	//vector<int> categorial;				/// calculated from the above (in learning, once train data is given. In testing - already ready)

	// missing value
	float missing_val = MED_MAT_MISSING_VALUE;	/// missing value
	string missing_method_str = "median";		/// how to handle missing values: median , left, right, mean, rand
	int missing_method = -1;					/// to be initialized from missing_method_str

	// sanities
	int test_for_inf = 1;				/// will fail on non finite values in input data	
	int test_for_missing = 0;			/// will fail if missing value found in data

	// prediction configuration
	int only_this_categ = -1;			/// relavant only to categorial predictions: -1: give all categs, 0 and above: give only those categs
										/// remember that currently 0 is a special category in TQRF : the control category (nothing happens, healthy, etc...)
	int predict_from_slice = -1;		/// will give predictions for slices [predict_from_slice,...,predict_to_slice]. if negative: all slices.
	int predict_to_slice = -1;
	int predict_sum_times = 0;			/// will sum predictions over different times


	// weights
	float case_wgt = 1;					/// the weight to use for cases with y!=0 in a weighted case

	// ada boost mode
	int nrounds = 1;					/// a single round means simply running TQRF as defined with no boosting applied
	float min_p = (float)0.01;			/// minimal probability to trim to when recalculating weights
	float max_p = (float)0.99;			/// maximal probability to trip to when recalculating weights
	float alpha = 1;					/// shrinkage factor
	float wgts_pow = 2;					/// power for the pow(-log(p), wgts_pow) used for adaboost weights

	// lists
	float tuning_size = 0;				/// size of group to tune tree weights by.
	int tune_max_depth = 0;				/// max depth of a node to get a weight for. 0 means 1 weight per tree.
	int tune_min_node_size = 0;			/// min node size for a node to have a weight

	// tuning gradient descent parameters
	float gd_rate = (float)0.01;		/// gradient descent step size
	int gd_batch = 1000;				/// gradient descent batch size
	float gd_momentum = (float)0.95;	/// gradient descent momentum
	float gd_lambda = 0;				/// regularization
	int gd_epochs = 0;					/// 0 : stop automatically , Otherwise: do this number of epochs

	// verbosity
	int verbosity = 0;					/// for debug prints
	int ids_to_print = 30;				/// control debug prints in certain places
	int debug = 0;						/// extra param for use when debugging

	//========================================================================================================================


	/// initialization from string
	int init(map<string, string>& map);

	// next are non serialized helpers we keep here as they are common to ALL the forest
	vector<double> log_table;


	// Serialization
	ADD_CLASS_NAME(TQRF_Params)
	ADD_SERIALIZATION_FUNCS(init_string, samples_time_unit, samples_time_unit_i, ncateg, time_slice_unit, time_slice_unit_i, time_slice_size, n_time_slices,
		time_slices, time_slices_wgts, censor_cases, max_q, tree_type, tree_type_i, ntrees, max_depth, min_node_last_slice, min_node, random_split_prob, ntry, ntry_prob,
		nsplits, max_node_test_samples, single_sample_per_pid, bag_with_repeats, bag_prob, bag_ratio, bag_feat, qpoints_per_split, nvals_for_categorial, categorial_str,
		categorial_tags, missing_val, missing_method_str, missing_method, test_for_inf, test_for_missing, only_this_categ, predict_from_slice, predict_to_slice,
		predict_sum_times, case_wgt, nrounds, min_p, max_p, alpha, wgts_pow, tuning_size, tune_max_depth, tune_min_node_size, gd_rate, gd_batch, gd_momentum,
		gd_lambda, gd_epochs, verbosity, ids_to_print, debug)

	//ADD_SERIALIZATION_FUNCS(init_string, samples_time_unit, time_slice_unit, time_slice_size, time_slices, max_q, max_q_sample, tree_type, ntrees, max_depth, min_node_last_slice, min_node, )
};


//==========================================================================================================================
// contains all the needed data for training including all quantizations (features, time slices) that are needed
//==========================================================================================================================
class Quantized_Feat : public SerializableObject {

public:
	vector<vector<short>> qx;		/// a vector of features that mimics the input x_in features matrix, but coded into quantized values
	vector<vector<float>> q_to_val; /// from a q value to float value : q=0 is reserved for missing value
									/// the range for q>0 is : [q_to_val[q], q_to_val[q+1])
	const MedFeatures *orig_medf;			   /// pointer to the original MedFeatures

	int nfeat = 0;					/// just an easy helper that = qx.size()
	vector<string> feature_names;	/// useful for debugging
	int ncateg = 0;					/// ncateg 0 is regression, otherwise categories are assumed to be 0 ... ncateg-1
	vector<const vector<float> *> orig_data; /// pointers to the original data given
	vector<string> feat_names;		   /// as given in train
	vector<float> y;
	vector<int> y_i;
	vector<int> last_time_slice;	/// when there's more than 1 time slice there may be censoring involved and the last_time_slice is the last uncensored one.
	int n_time_slices;				/// 1 time slice is simply the regular case of a label for the whole future
	vector<int> slice_counts[2];	/// counts of elements in slices (in case of non regression trees). slices with no variability are not interesting.

	vector<vector<int>> lists;		/// lists[0] is always the lines used for training the trees in round 1
									/// the others can be used for later stages, for example :
									/// lists[1] could be used for early stopping measurements or for estimating weights for trees/nodes

	vector<int> is_categorial_feat;

	// next are pre computed for bagging purposes
	vector<vector<vector<int>>> time_categ_pids;
	vector<vector<vector<int>>> time_categ_idx;
	vector<vector<int>> categ_pids;
	vector<vector<int>> categ_idx;
	unordered_map<int, vector<vector<vector<int>>>> pid2time_categ_idx;

	// next are helper arrays used when doind adaboost
	vector<float> wgts;
	vector<float> orig_wgts;
	vector<float> probs;
	vector<float> w_to_sum;
	vector<vector<float>> sum_over_trees;
	float alpha0;

	int init(const MedFeatures &medf, TQRF_Params &params);

	~Quantized_Feat() {  pid2time_categ_idx.clear(); }

private:
	int quantize_feat(int i_feat, TQRF_Params &params);
	int init_time_slices(const MedFeatures &medf, TQRF_Params &params);
	int init_pre_bagging(TQRF_Params &params);
	int init_lists(const MedFeatures &medf, TQRF_Params &params);

};


//==========================================================================================================================
// a basic node class : currently a single node type serves all trees .... could be changed to 
//==========================================================================================================================
class TQRF_Node : public SerializableObject {
public:
	// Next are must for every node and are ALWAYS serialized
	int node_idx = -1;			/// for debugging and prints
	int i_feat = -1;				/// index of feature used in this node
	float bound = (float)-1e10;		/// samples with <= bound go to Left , the other to Right
	int is_terminal = 0;
	int left_node = -1;
	int right_node = -1;
	int depth = -1;
	int missing_direction = TQRF_MISSING_DIRECTION_RAND_EACH_SAMPLE;	/// 0: left , 1: right , 2: randomize each sample

	// next are needed while learning , and if asked to keep samples in nodes - we keep them always for now
	int from_idx = -1;		/// the node elements are those given in its tree indexes from place from_idx, to to_idx.
	int to_idx = -1;
	int size() { return to_idx-from_idx+1; }

	int node_serialization_mask = 0x1; /// choose which of the following to serialize
	int beta_idx = -1;

	// categorical : mask |= 0x1 , time_categ_count[t][c] : how many counts in this node are in timeslice t and category c
	vector<vector<float>> time_categ_count;

	// regression : mask |= 0x2
	//float pred_mean = (float)-1e10;
	//float pred_std = (float)1;

	// quantiles: mask |= 0x4
	//vector<pair<float, float>> quantiles;


	// following are never serialized - only for learn time
	int state = TQRF_Node_State_Initiated; // 0 - created 1 - in process 2 - done with

	ADD_CLASS_NAME(TQRF_Node)
	ADD_SERIALIZATION_FUNCS(node_idx, i_feat, bound, is_terminal, left_node, right_node, depth, missing_direction, from_idx, to_idx,
		node_serialization_mask, beta_idx, time_categ_count, state)
};

//----------------------------------------------------------------------------------------------------------------------
class TQRF_Node_Categorial : TQRF_Node {

};

//==========================================================================================================================
// Split_Stat contains the quantized data structures used in order to make a split decision
// Basically : 
// for categorial outcomes:
// for each time slot, for each quanta -> counts for each category
// for regression outcomes:
// for each time slot, for each quanta -> nvals and sum (??)
//==========================================================================================================================
class TQRF_Split_Stat {

public:
	// categorial case
	// vector<vector<vector<int>>> counts; /// counts[t][q][c] = how many counts were in time slot t, quanta q, and category c.

	// suggestion (categorial case):
	// TQRF_Split will get a node with indexes, and then:
	// (1) Go over all the samples in the node and for each one add +1 to the relevant counts[t][q][c]
	// (2) Will enable going over counts and choose the best q for splitting given some params (score type, minimal size, minimal change, etc)
	//
	// This may allow a very elegant code for a tree in which all the hard stuff is implemented inside.
	// 
	// issues to think about:
	// parallelism over nodes
	// memory allocation - we want it single time
	// tricks for efficient calculation of counts and the scores

	virtual ~TQRF_Split_Stat() {};

	virtual int init(Quantized_Feat &qf, TQRF_Params &params) {	return 0;};
	virtual int prep_histograms(int i_feat, TQRF_Node &node, vector<int> &indexes, Quantized_Feat &qf, TQRF_Params &params) { return 0; };
	virtual int get_best_split(TQRF_Params &params, int &best_q, double &best_score) { return 0; };

	int get_q_test_points(int feat_max_q, TQRF_Params &params, vector<int> &qpoints);


	// helper vector for qpoints
	vector<int> qpoints;
	
	// debug
	virtual void print_histograms() { return; };

	// the actual number of q values used (full or after qpoints squeeze if it was done)
	int counts_q = 0;


	static TQRF_Split_Stat *make_tqrf_split_stat(int tree_type);
};

//==========================================================================================================================
class TQRF_Split_Categorial : public TQRF_Split_Stat {
public:
	// categorial case : counts[t][q][c] : time_slice , quanta, category : number of elements
	vector<vector<vector<int>>> counts;

	// sums[t][c] = number of samples in time slice t and category c summed on all q vals
	// this is needed for a more efficient computation of scores later
	vector<vector<int>> sums;

	// sums_t[t] = number of samples in time slice t (needed later for more efficient calculations)
	vector<int> sums_t;

	int total_sum = 0; // sum of the sum_t vector

	~TQRF_Split_Categorial() {};

	// next are for easy access
	int ncateg = 0;
	int nslices = 0;
	int maxq = 0; // overall

	// API's
	int init(Quantized_Feat &qf, TQRF_Params &params);
	int prep_histograms(int i_feat, TQRF_Node &node, vector<int> &indexes, Quantized_Feat &qf, TQRF_Params &params);
	//virtual int get_best_split(TQRF_Params &params, int &best_q, float &best_score);

	void print_histograms();
};


//==========================================================================================================================
class TQRF_Split_Likelihood : public TQRF_Split_Categorial {
public:
	int get_best_split(TQRF_Params &params, int &best_q, double &best_score);
};


//==========================================================================================================================
class TQRF_Split_Entropy : public TQRF_Split_Categorial {
public:
	int get_best_split(TQRF_Params &params, int &best_q, double &best_score);
};

//==========================================================================================================================
class TQRF_Split_Weighted_Categorial : public TQRF_Split_Stat {
public:
	// categorial case : counts[t][q][c] : time_slice , quanta, category : number of elements
	vector<vector<vector<float>>> counts;

	// sums[t][c] = number of samples in time slice t and category c summed on all q vals
	// this is needed for a more efficient computation of scores later
	vector<vector<float>> sums;

	// sums_t[t] = number of samples in time slice t (needed later for more efficient calculations)
	vector<float> sums_t;

	float total_sum = 0; // sum of the sum_t vector

	~TQRF_Split_Weighted_Categorial() {};

	// next are for easy access
	int ncateg = 0;
	int nslices = 0;
	int maxq = 0; // overall

				  // API's
	int init(Quantized_Feat &qf, TQRF_Params &params);
	int prep_histograms(int i_feat, TQRF_Node &node, vector<int> &indexes, Quantized_Feat &qf, TQRF_Params &params);
	//virtual int get_best_split(TQRF_Params &params, int &best_q, float &best_score);

	void print_histograms();
};

//==========================================================================================================================
class TQRF_Split_Weighted_Likelihood : public TQRF_Split_Weighted_Categorial {
public:
	int get_best_split(TQRF_Params &params, int &best_q, double &best_score);
};

//==========================================================================================================================
class TQRF_Split_Dev : public TQRF_Split_Categorial {
public:
	~TQRF_Split_Dev() {};
	int get_best_split(TQRF_Params &params, int &best_q, double &best_score) { return 0; };
};

//==========================================================================================================================
class TQRF_Split_Regression : public TQRF_Split_Stat {
public:
	// categorial case
	vector<vector<vector<pair<float,int>>>> sum_num;

	int init(Quantized_Feat &qf, TQRF_Params &params) { return 0; };
	int prep_histograms(int i_feat, TQRF_Node &node, vector<int> &indexes, Quantized_Feat &qf, TQRF_Params &params) { return 0; };
	int get_best_split(TQRF_Params &params, int &best_q, double &best_score) { return 0; };
};

//==========================================================================================================================
// A tree base class
//==========================================================================================================================
class TQRF_Tree : public SerializableObject {

public:
	// next are needed also for predictions, and hence should be serialized
	int tree_type;
	int id;						// for debug prints - a specific tree identifier
	int keep_indexes = 0;
	vector<int> indexes;		// indexes[i] = an index of a sample in the given Quantized_Feat
	vector<TQRF_Node> nodes;	// this node supports currently all possible nodes for all trees... to save ugly templated code

	// next variables are no-need-to-serialize helpers
	vector<int> i_feats;		// feature indexes to be used in this tree (they can be bagged as well)

	TQRF_Tree() {};

	void init(Quantized_Feat &qfeat, TQRF_Params &params) { _qfeat = &qfeat; _params = &params; }
	int Train(Quantized_Feat &qfeat, TQRF_Params &params) {	init(qfeat, params); return Train(); }

	const TQRF_Node *Get_Node_for_predict(MedMat<float> &x, int i_row, float missing_val, int &beta_idx) const;
	TQRF_Node *Get_Node(MedMat<float> &x, int i_row, float missing_val);

	int Train();

	// helpers inside Train:

	// get indexes vector ready
	int get_bagged_indexes();

	// initialize root node
	int init_root_node();

	// get the next node to work on
	int get_next_node(int curr_node);

	// get the list of features to work on
	int	get_feats_to_test(vector<int> &feats_to_test);

	// init the vector for splits: to nfeat, right sizes, and right type + a free mem api
	int init_split_stats(vector<TQRF_Split_Stat *> &tqs);
	void free_split_stats(vector<TQRF_Split_Stat *> &tqs);

	// close work on current node and make the split if needed
	int node_splitter(int i_curr_node, int i_best, int q_best);

	// once a node is finalized : prepares its internal counts (with or without taking weights into account)
	float prep_node_counts(int i_curr_node, int use_wgts_flag);


	void pre_serialization() { if (keep_indexes == 0) indexes.clear(); }
	ADD_CLASS_NAME(TQRF_Tree)
	ADD_SERIALIZATION_FUNCS(tree_type, id, keep_indexes, indexes, nodes)

private:
	Quantized_Feat *_qfeat;
	TQRF_Params *_params;

	// next used to manage nodes while building
	int n_nodes_in_process = 0;
	int i_last_node_in_process = 0;

	int bag_chooser(float p, int _t, int _c, /* OUT APPEND */ vector<int> &_indexes);
	int bag_chooser(int choose_with_repeats, int single_sample_per_id, float p, vector<int> &pids, vector<int> &idx, unordered_map<int, vector<int>> &pid2idx, /* OUT APPEND */ vector<int> &_indexes);


};



//==========================================================================================================================
class TQRF_Forest : public SerializableObject {

public:

	TQRF_Params params;
	vector<TQRF_Tree> trees;
	vector<float> alphas;
	vector<float> betas;

	int init(map<string, string>& map) { return params.init(map); }
	int init_from_string(string init_string) { params.init_string = init_string; return SerializableObject::init_from_string(init_string); }

	void init_tables(Quantized_Feat &qfeat);

	/// The basic train matrix for TQRF is MedFeatures (!!) the reason is that it contains everything in one place:
	/// that is: the X features, the Y outcome, the weights and the samples for each row.
	/// All of these are needed when calculating a logrank score for example
	/// The y matrix is added since we may want to use regression with y values given for every time slice ...
	int Train(const MedFeatures &medf, const MedMat<float> &Y);
	int Train(const MedFeatures &medf) { MedMat<float> dummy; return Train(medf, dummy); }

	int Train_AdaBoost(const MedFeatures &medf, const MedMat<float> &Y);
	int update_counts(vector<vector<float>> &sample_counts, MedMat<float> &x, Quantized_Feat &qf, int zero_counts, int round);

	/// tuning : solving a gd problem of finding the optimal betas for nodes at some certain chosen depth in the trees on a kept-a-side 
	/// set of samples.
	int tune_betas(Quantized_Feat &qfeat);
	int solve_betas_gd(MedMat<float>& C, MedMat<float>& S, vector<float> &b);

	/// However - the basic predict for this model is MedMat !! , as here it is much simpler :
	/// we only need to find the terminal nodes in the trees and calculate our scores
	int Predict(MedMat<float> &x, vector<float> &preds) const;
	int n_preds_per_sample() const;

	int Predict_Categorial(MedMat<float> &x, vector<float> &preds) const; // currently like this... with time should consider inheritance to do it right.

	// print average bagging reports
	void print_average_bagging(int _n_time_slices, int _n_categ);

	// simple helpers
	static int get_tree_type(const string &str);
	static int get_missing_value_method(const string &str);

	ADD_CLASS_NAME(TQRF_Forest)
	ADD_SERIALIZATION_FUNCS(params, trees, alphas, betas)

private:

};

// helper struct
struct TreeNodeIdx {
	int i_tree = -1;
	int i_node = -1;
	TreeNodeIdx(int i_t, int i_n) { i_tree = i_t; i_node = i_n; }
};


/// next is for debugging
class AllIndexes : public SerializableObject {

public:

	vector<vector<int>> all_indexes;

	void init_all_indexes(vector<TQRF_Tree> &trees) {
		for (auto &tree : trees)
			all_indexes.push_back(tree.indexes);
	}

	ADD_CLASS_NAME(AllIndexes)
	ADD_SERIALIZATION_FUNCS(all_indexes)
};

//========================================
// Join the serialization Waggon
//========================================
MEDSERIALIZE_SUPPORT(TQRF_Params);
MEDSERIALIZE_SUPPORT(TQRF_Node);
MEDSERIALIZE_SUPPORT(TQRF_Tree);
MEDSERIALIZE_SUPPORT(TQRF_Forest);
MEDSERIALIZE_SUPPORT(AllIndexes)

#endif