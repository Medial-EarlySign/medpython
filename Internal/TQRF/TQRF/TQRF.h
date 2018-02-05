#pragma once
//
// TQRF
//

#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <MedUtils/MedUtils/MedMat.h>
#include <MedTime/MedTime/MedTime.h>
#include <queue>

using namespace std;

enum TQRF_TreeTypes {
	TQRF_TREE_ENTROPY = 0,
	TQRF_TREE_REGRESSION = 1,
	TQRF_TREE_LOGRANK = 3,
	TQRF_TREE_UNDEFINED = 4
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

#define TQRF_MAX_TIME_SLICE			10000000


//==========================================================================================================================
class TQRF_Params : SerializableObject {
public:

	//========================================================================================================================
	// params list
	//========================================================================================================================
	string init_string = "";			/// sometimes it helps to keep it for debugging
	string samples_time_unit = "Date";
	
	// time slices
	string time_slice_unit = "Days";
	int time_slice_size = -1;			/// the size of the basic time slice, -1: is like infinity: a single time slice like a regular QRF
	int n_time_slices = 1;				/// if time_slices vector is not given, one will be created using time_slice_size and this parameter. 
	vector<int> time_slices = {};		/// if not empty: defines the borders of all the time lines. Enables a very flexible time slicing strategy
	
	// quantization
	int max_q = 200;					/// maximal quantization
	//int max_q_sample = 100000;			/// the max number of values to use when deciding q limits

	// trees and stopping criteria
	string tree_type = "";				/// options: regression, entropy, logrank
	int tree_type_i = -1;				/// tree type code : calulated from tree_type the string
	int	ntrees = 50;					/// number of trees to learn
	int max_depth = 100;				/// maximal depth of tree
	int min_node_last_slice = 10;		/// stopping criteria : minimal number of samples in a node in the last time slice
	int min_node = 5000;				/// stopping criteria : minimal number of samples in a node in the first time slice

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
	int test_for_missing = 0;			/// will fail on if missing value found in data

	// verbosity
	int verbosity = 0;					/// for debug prints

	//========================================================================================================================


	/// initialization from string
	int init(map<string, string>& map);

	// Serialization
	//ADD_SERIALIZATION_FUNCS(init_string, samples_time_unit, time_slice_unit, time_slice_size, time_slices, max_q, max_q_sample, tree_type, ntrees, max_depth, min_node_last_slice, min_node, )
};


//==========================================================================================================================
// contains all the needed data for training including all quantizations (features, time slices) that are needed
//==========================================================================================================================
class Quantized_Feat : SerializableObject {

public:
	vector<vector<short>> qx;		/// a vector of features that mimics the input x_in features matrix, but coded into quantized values
	vector<vector<float>> q_to_val; /// from a q value to float value : q=0 is reserved for missing value
									/// the range for q>0 is : [q_to_val[q], q_to_val[q+1])
	MedFeatures *orig_medf;			   /// pointer to the original MedFeatures

	int nfeat = 0;					/// just an easy helper that = qx.size()
	int ncateg = 0;					/// ncateg 0 is regression, otherwise categories are assumed to be 0 ... ncateg-1
	vector<vector<float> *> orig_data; /// pointers to the original data given
	vector<string> feat_names;		   /// as given in train
	vector<float> y;
	vector<int> last_time_slice;	/// when there's more than 1 time slice there may be censoring involved and the last_time_slice is the last uncensored one.
	int n_time_slices;				/// 1 time slice is simply the regular case of a label for the whole future
	vector<int> is_categorial_feat;

	int init(MedFeatures &medf, TQRF_Params &params);

private:
	int quantize_feat(int i_feat, TQRF_Params &params);

};


//==========================================================================================================================
// a basic node class : currently a single node type serves all trees .... could be changed to 
//==========================================================================================================================
class TQRF_Node : SerializableObject {
public:
	// Next are must for every node and are ALWAYS serialized
	int node_idx;			/// for debugging and prints
	int i_feat = -1;				/// index of feature used in this node
	float bound = (float)-1e10;		/// samples with <= bound go to Left , the other to Right
	int is_terminal = 0;
	int left_idx = -1;
	int right_idx = -1;
	int depth = -1;

	// next are needed while learning , and if asked to keep samples in nodes - we keep them always for now
	int from_idx = -1;		/// the node elements are those given in its tree indexes from place from_idx, to to_idx.
	int to_idx = -1;
	int size() { return to_idx-from_idx+1; }

	int node_serialization_mask = 0x1; /// choose which of the following to serialize

	// categorical : mask |= 0x1
	vector<vector<int>> time_categ_count;

	// regression : mask |= 0x2
	float pred_mean = (float)-1e10;
	float pred_std = (float)1;

	// quantiles: mask |= 0x4
	vector<pair<float, float>> quantiles;


	// following are never serialized - only for learn time
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

};

//==========================================================================================================================
// A tree base class
//==========================================================================================================================
class TQRF_Tree : SerializableObject {

public:
	int tree_type;
	int id;						// for debug prints - a specific tree identifier
	int keep_indexes = 0;
	vector<int> indexes;		// indexes[i] = an index of a sample in the given Quantized_Feat
	vector<TQRF_Node> nodes;	// this node supports currently all possible nodes for all trees... to save ugly templated code
	vector<int> i_feats;		// feature indexes to be used in this tree (they can be bagged as well)

	void init(Quantized_Feat &qfeat, TQRF_Params &params) { _qfeat = qfeat; _params = params; }
	int Train(Quantized_Feat &qfeat, TQRF_Params &params) {	init(qfeat, params); return Train(); }

	int Train();

	// Step(1) in Train: get indexes vector ready
	int get_bagged_indexes();

	// Step(2) in Train: initialize root
	int init_root_node();

	// Step(3) : train and add new nodes as long as needed
	int train_node(int i_node);

	// within train_node()

	// Step(3.1) : get the relevant histograms 
	int get_feature_stat();

	// static TQRF_Tree* make_tree(const string &type);

private:
	Quantized_Feat &_qfeat;
	TQRF_Params &_params;

	// next used to manage nodes while building
	int n_nodes_in_process = 0;
	int i_last_node_in_process = 0;

	int bag_chooser(int choose_with_repeats, int single_sample_per_id, float p, vector<int> &pids, vector<int> &idx, unordered_map<int, vector<int>> &pid2idx, /* OUT APPEND */ vector<int> &_indexes);


};



//==========================================================================================================================
class TQRF_Forest : SerializableObject {

public:

	TQRF_Params params;
	vector<TQRF_Tree> trees;

	/// The basic train matrix for TQRF is MedFeatures (!!) the reason is that it contains everything in one place:
	/// that is: the X features, the Y outcome, the weights and the samples for each row.
	/// All of these are needed when calculating a logrank score for example
	/// The y matrix is added since we may want to use regression with y values given for every time slice ...
	int Train(const MedFeatures &medf, const MedMat<float> &Y);
	int Train(const MedFeatures &medf) { MedMat<float> dummy; return Train(medf, dummy); }



	// simple helpers
	static int get_tree_type(const string &str);
	static int get_missing_value_method(const string &str);
private:

};

//========================================
// Join the serialization Waggon
//========================================