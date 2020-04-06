//#pragma once
#ifndef _ITERATIVEIMPUTER_H_
#define _ITERATIVEIMPUTER_H_
#include <MedProcessTools/MedProcessTools/FeatureProcess.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <unordered_set>

//=========================================================================
// IterativeImputer :: Inner Class : the actual implementation
//=========================================================================
class IterativeImputerParams : public SerializableObject {
public:
	// params
	vector<string> features_to_impute = {}; // if empty : impute all
	string regressor = "qrf";
	string regressor_params = "type=regression;ntrees=30;min_node=100;spread=0.1;learn_nthreads=40;predict_nthreads=40";
	//	string regressor = "linear_model";
	//	string regressor_params = "rfactor=0.9";
	string multi_categ_classifier = "qrf";
	string multi_categ_classifier_params = "type=categorical_entropy;ntrees=30;min_node=10;learn_nthreads=40;predict_nthreads=40";
	string add_ncateg_var_name = "n_categ"; // for algorithms that need it
	string round1_strata = "Age,40,80,5:Gender,1,2,1";
	int do_round1 = 1;
	int round1_moment = 0; // 0 mean , 1 median
	int categorial_bound = 0;  // if n diff values for a feature < this bound --> we declare it categorial
	int max_iterations = 1;
	float p_validation = (float)0.1; // helping to print intermediate results
	int min_vals_for_training = 10000;
	float missing_value = MED_MAT_MISSING_VALUE;
	int round_to_resolution = 1;
	int verbose = 1;
	float missing_bound = 0.5;

	// a few constant params
	int min_vals_for_impute = 1000; // min number of samples to run impute rounds.


	int init(map<string, string>& mapper);

	ADD_CLASS_NAME(IterativeImputerParams)
	ADD_SERIALIZATION_FUNCS(features_to_impute, regressor, regressor_params, multi_categ_classifier, multi_categ_classifier_params, add_ncateg_var_name,
		round1_strata, do_round1, round1_moment, categorial_bound, max_iterations, p_validation, min_vals_for_training, missing_value,
		round_to_resolution, verbose, missing_bound);

};


class feature_info : public SerializableObject {
public:
	string name = "";
	string full_name = "";
	int n_diff_vals = 0;
	int is_categorial = 0;
	float min = (float)1e10;
	float max = (float)-1e10;
	float resolution = 0;
	int predictor_type = 0;
	vector<int> inds_for_pred;

	// no need to serialize helpers
	int n_missing = 0;
	int n_with_values = 0;
	int n_with_non_zero_values = 0;
	float *data = NULL;
	int data_len = 0;
	vector<char> is_missing;
	vector<string> feats_for_pred; // full names of features to be used in the prediction matrix
	vector<int> train_idx;
	vector<int> test_idx;
	vector<int> pred_idx;

	int prep_feats_for_pred(MedFeatures &mfd);
	int prep_indexes(const vector<int> &external_train_idx, const vector<int> &external_test_idx, float missing_value);


	void print() {
		fprintf(stderr, "Feature Info :: %s :: %s :: data_len %d : n_missing %d ( %5.2f ): n_with %d ( non zero %d ): n_diff_vals %d : categorial %d : min %f : max %f : resolution %f\n",
			name.c_str(), full_name.c_str(), data_len, n_missing, (float)100 * n_missing / (float)data_len, n_with_values, n_with_non_zero_values, n_diff_vals, is_categorial, min, max, resolution);
	}

	ADD_CLASS_NAME(feature_info)
	ADD_SERIALIZATION_FUNCS(name, full_name, n_diff_vals, is_categorial, min, max, resolution, predictor_type, inds_for_pred);
};

/// \n
/// IterativeImputer\n
///\n
/// A general strong imputer that does the following:\n
/// (1) Runs a simple stratified imputer\n
/// (2) Runs iterations completing values (from the least missing to the max missing) where:\n
///     (a) continuous values are calculated using a regressor\n
///     (b) categorial values (less than some bound) are calculated using a multi category classifier\n
/// (3) Repeats the process several times until it converges or until max_iters is reached.\n
///
class IterativeImputer : public SerializableObject {
public:

	IterativeImputerParams params; // needs serialization

	// helpers not for serialization
	vector<int> train_idx, test_idx;
	vector<char> is_train;
	MedFeatures learn_features_map; // in the learn stage we must keep a copy of our MedFeatures since we change it with each step and iteration
	unordered_set<int> train_ids, test_ids;


	// helpers for (partial) serialization
	vector<feature_info> feats;

	// First round imputers (needs serialization)
	vector<FeatureImputer> first_round_imputers;

	// Iterations imputers for each split (needs serialization)
	vector<int> predictors_order;
	vector<vector<MedPredictor *>> predictors; // predictors[i][j] = the predictor at iteration i , for feature j

	/// The parsed fields from init command.
	/// @snippet IterativeImputer.cpp IterativeImputerParams::init
	int init(map<string, string>& mapper) { return params.init(mapper); }

	int init_internals(MedFeatures &mfd);
	int init_feature_info(MedFeatures &mfd, string feat_name);
	int init_feature_info_update(MedFeatures &mfd, feature_info &fi);
	int round_to_resolution(MedFeatures &mfd);
	int round_arr(float *arr, int len, float resolution, float _min, float _max);
	int learn_first_round(MedFeatures &mfd);
	int learn_iteration(MedFeatures &mfd, int iter);
	int apply_first_round(MedFeatures &mfd, bool learning);
	int apply_iteration(MedFeatures &mfd, int iter);

	int find_feats_to_learn_from(int f_idx);
	int feats_for_pred_inds_to_names(feature_info &fi);


	//int Learn(MedFeatures &mfd) { init_internals(mfd); }
	//int Apply(MedFeatures &mfd) { fprintf(stderr, "IterativeImputer::Apply() NOT IMPLEMENTED YET\n"); }
	int Learn(MedFeatures &mfd);
	int Apply(MedFeatures &mfd, bool learning);

	// Serialization
	ADD_CLASS_NAME(IterativeImputer)
	ADD_SERIALIZATION_FUNCS(params, feats, first_round_imputers, predictors_order, predictors)

};


//=============================================================================
// FeatureIterativeImputer :: Wrapper of IterativeImputer as a FeatureProcessor
//=============================================================================
class FeatureIterativeImputer : public FeatureProcessor {
public:

	// holding an IterativeImputer instance
	IterativeImputer imputer;

	// Constructor
	FeatureIterativeImputer() { processor_type = FTR_PROCESS_ITERATIVE_IMPUTER; }

	/// The parsed fields from init command.
	/// @snippet IterativeImputer.cpp IterativeImputerParams::init
	int init(map<string, string>& mapper) { return imputer.init(mapper); }

	// Learn imputing model
	int Learn(MedFeatures& features, unordered_set<int>& ids);
	int Learn(MedFeatures& features) { return imputer.Learn(features); }

	// Apply cleaning model
	int _apply(MedFeatures& features, unordered_set<int>& ids, bool learning);
	int _apply(MedFeatures& features, bool learning) { return imputer.Apply(features, learning); }

	/// check if a set of features is affected by the current processor
	bool are_features_affected(unordered_set<string>& out_req_features);

	/// update sets of required as input according to set required as output to processor
	void update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features);

	// Serialization
	ADD_CLASS_NAME(FeatureIterativeImputer)
	ADD_SERIALIZATION_FUNCS(processor_type, imputer)

};


//=======================================================
// Join the MedSerialize Wagon
//=======================================================
MEDSERIALIZE_SUPPORT(IterativeImputerParams)
MEDSERIALIZE_SUPPORT(feature_info)
MEDSERIALIZE_SUPPORT(IterativeImputer)
MEDSERIALIZE_SUPPORT(FeatureIterativeImputer)
#endif