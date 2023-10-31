/// @file
#ifndef __EXPLAIN_WRAPPER_H__
#define __EXPLAIN_WRAPPER_H__

#include <vector>
#include <string>
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedProcessTools/MedProcessTools/PostProcessor.h>
#include <MedStat/MedStat/GibbsSampler.h>
#include <MedAlgo/MedAlgo/tree_shap.h>
#include <MedAlgo/MedAlgo/SamplesGenerator.h>

using namespace std;

/**
* Parameters for filtering explanations
*/
class ExplainFilters : public SerializableObject {
public:
	int sort_mode; ///< 0 - both pos and negative (sorted by abs), -1 - only negatives, +1 - only positives
	int max_count; ///< maximal limit to take as explain features. 0 - no limit
	float sum_ratio; ///< percentage of sum of explain values to take from sort_mode. [0 - 1]

	ExplainFilters();

	int init(map<string, string> &map);

	/// commit filterings
	void filter(map<string, float> &explain_list) const;

	ADD_CLASS_NAME(ExplainFilters)
		ADD_SERIALIZATION_FUNCS(sort_mode, max_count, sum_ratio)

};

/**
* Processings of explanations - grouping, Using covariance matrix for taking feature correlations into account.
*/
class ExplainProcessings : public SerializableObject {
private:
	bool postprocessing_cov = false; ///< should covariance correction be done in the postprocessing stage (this is not the case when working iteratively)
public:
	bool group_by_sum = false; ///< If true will do grouping by sum of each feature, otherwise will use internal special implementation
	bool learn_cov_matrix = false; ///< If true will learn cov_matrix
	int normalize_vals = 0; ///< If != 0 will normalize contributions. 1: normalize by sum of (non b0) abs of all contributions 2: same, but also corrects for groups
	int zero_missing = 0; ///<  if != 0 will throw bias terms and zero all contributions of missing values and groups of missing values
	bool keep_b0 = false; ///< if true will keep b0 prior
	bool iterative = false; ///< if true will add explainers iteratively, conditioned on those already selected
	int iteration_cnt = 0; ///< if >0 the maximal number of iterations
	bool use_max_cov = false; ///< If true will use max cov logic

	bool use_mutual_information; ///< if true will use mutual information instead of covariance
	BinSettings mutual_inf_bin_setting; ///< the bin setting for mutual information

	MedMat<float> abs_cov_features; ///< absolute values of covariance features for matrix.either read from file (and then apply absolute value), or learn if learn_cov_matrix is on , 

	string grouping; ///< grouping file or "BY_SIGNAL" keyword to group by signal or "BY_SIGNAL_CATEG" - for category signal to split by values (aggreagates time windows) or "BY_SIGNAL_CATEG_TREND" - also splitby TRENDS
	vector<vector<int>> group2Inds;
	vector<string> groupNames;
	map<string, vector<int>> groupName2Inds;

	ExplainProcessings();

	int init(map<string, string> &map);

	/// Learns process - for example cov matrix
	void learn(const MedFeatures &train_mat);

	/// commit processings
	void process(map<string, float> &explain_list) const;
	// same as process but zero-ing all contributions of missing values features and groups with all the participants inside missing
	void process(map<string, float> &explain_list, unsigned char *missing_value_mask) const;

	/// helper func: returns the normalized contribution for a specific group given original contributions
	float get_group_normalized_contrib(const vector<int> &group_inds, vector<float> &contribs, float total_normalization_factor) const;

	void post_deserialization();

	///Creates the feature groups from the argument file_name and by existing features
	static void read_feature_grouping(const string &file_name, const vector<string>& features, vector<vector<int>>& group2index,
		vector<string>& group_names, bool verbose = true);

	ADD_CLASS_NAME(ExplainProcessings)
		ADD_SERIALIZATION_FUNCS(group_by_sum, abs_cov_features, normalize_vals, zero_missing, groupNames, group2Inds, keep_b0,
			iterative, iteration_cnt, postprocessing_cov, use_mutual_information, mutual_inf_bin_setting, use_max_cov)
};

/**
* A wrapper class to hold all global arguments needed for ModelExplainer
*/
class GlobalExplainerParams : public SerializableObject {
public:
	string attr_name = ""; ///< attribute name for explainer
	bool store_as_json = false; ///< If true will store ButWhy output as json in string attributes
	bool denorm_features = true; ///< If true will save feature values denorm

	//No init - will be initialized directly in ModelExplainer::init

	ADD_CLASS_NAME(GlobalExplainerParams)
		ADD_SERIALIZATION_FUNCS(attr_name, store_as_json, denorm_features)
};

/**
* An abstract class API for explainer
*/
class ModelExplainer : public PostProcessor {
private:
	/// init function for specific explainer
	virtual void _init(map<string, string> &mapper) = 0;
	/// map from feature names to thier normalizer if exsits. Get filled after init_post_processor. no need to serialize
	unordered_map<string, const FeatureNormalizer *> feats_to_norm;
public:
	MedPredictor * original_predictor = NULL; ///< predictor we're trying to explain
	ExplainFilters filters; ///< general filters of results
	ExplainProcessings processing; ///< processing of results, like groupings, COV
	GlobalExplainerParams global_explain_params;
	
	/// Global init for general args in all explainers. initialize directly all args in GlobalExplainerParams
	virtual int init(map<string, string> &mapper);

	virtual int update(map<string, string>& mapper);
	/// overload function for ModelExplainer - easier API
	virtual void _learn(const MedFeatures &train_mat) = 0;

	///Learns from predictor and train_matrix (PostProcessor API)
	virtual void Learn(const MedFeatures &train_mat);
	void Apply(MedFeatures &matrix) { explain(matrix); } ///< alias for explain

	void get_input_fields(vector<Effected_Field> &fields) const;
	void get_output_fields(vector<Effected_Field> &fields) const;

	///Init ModelExplainer from MedModel - copies predictor pointer, might save normalizers pointers
	void init_post_processor(MedModel& model);

	///Virtual - return explain results in sample_feature_contrib
	virtual void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const = 0;

	/// Stores explain results in matrix
	virtual void explain(MedFeatures &matrix) const; //stores _explain results in MedFeatures

	static void print_explain(MedSample &smp, int sort_mode = 0);

	void dprint(const string &pref) const;

	virtual ~ModelExplainer() {};
};

/** @enum
* Tree Explainer modes - will be choosen by explainer
*/
enum TreeExplainerMode {
	ORIGINAL_IMPL = 0,
	CONVERTED_TREES_IMPL = 1,
	PROXY_IMPL = 2
};

/**
* A generic tree explainer:
* 1. Reads tree model into structure to calc SHAP values - QRF, XGB, lightGBM (future - BART ?)
* 2. for xgboost/lightGBM where conversion fails, use internal implementation of TreeShap
* 3. train LightGBM/Xgboost proxy to explain non-tree Predictor
*/
class TreeExplainer : public ModelExplainer {
private:
	MedPredictor * proxy_predictor = NULL; //uses this if model has no tree implementation
	//Tree structure of generic ensamble trees
private:
	bool convert_qrf_trees();
	bool convert_lightgbm_trees();
	bool convert_xgb_trees();
	void _init(map<string, string> &mapper);
public:
	bool try_convert_trees();
	TreeEnsemble generic_tree_model;
	string proxy_model_type = ""; ///< proxy predictor type to relearn original predictor output with tree models
	string proxy_model_init = ""; ///< proxy predictor arguments
	bool interaction_shap = false; ///< If true will calc interaction_shap values (slower)
	int approximate = false; ///< if true will run SAABAS alg - which is faster
	float missing_value = MED_MAT_MISSING_VALUE; ///< missing value
	bool verbose = false;

	TreeExplainer() { processor_type = FTR_POSTPROCESS_TREE_SHAP; }

	void init_post_processor(MedModel& model);

	TreeExplainerMode get_mode() const;

	void _learn(const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	void post_deserialization();

	~TreeExplainer();

	ADD_CLASS_NAME(TreeExplainer)
		ADD_SERIALIZATION_FUNCS(proxy_predictor, interaction_shap, filters, processing, global_explain_params, verbose)
};

/**
* Shapely Explainer - Based on learning training data to handle missing_values as "correct" input.
* so predicting on missing_values mask would give E[F(X)], whwn X has missing values, will allow much faster compution.
* All we need to do is reweight or manipulate missing_values(erase/add missing values) that each sample would look like:
* First we sample uniformally how many missing values should be - than randomally remove those value and set them as missing values.
* This will cause the weight of each count of missing values to be equal in train - same as weights in SHAP values calculation
*/
class MissingShapExplainer : public ModelExplainer {
private:
	MedPredictor * retrain_predictor = NULL; //the retrain model

	void _init(map<string, string> &mapper);

	float avg_bias_score;
public:
	int add_new_data; ///< how many new data data points to add for train according to sample masks
	bool no_relearn; ///< If true will use original model without relearn. assume original model is good enough for missing vals (for example LM model)

	int max_test; ///< max number of samples in SHAP
	float missing_value; ///< missing value 
	bool sample_masks_with_repeats; ///< Whether or not to sample masks with repeats
	float select_from_all; ///< If max_test is beyond this percentage of all options than sample from all options (to speed up runtime)
	bool uniform_rand; ///< it True will sample masks uniformlly
	bool use_shuffle; ///< if not sampling uniformlly, If true will use shuffle (to speed up runtime)
	string predictor_args; ///< arguments to change in predictor - for example to change it into regression
	string predictor_type;
	bool verbose_learn; ///< If true will print more in learn
	string verbose_apply; ///< If has value - output file
	float max_weight; ///< the maximal weight number. if < 0 no limit
	int subsample_train; ///< if not zero will use this to subsample original train sampels to this number
	int limit_mask_size; ///< if set will limit mask size in the train - usefull for minimal_set

	// parameters for minimal_set usage if use_minimal_set is true will do different thing
	bool use_minimal_set; ///< If true will use different method to find minimal set
	float sort_params_a; ///< weight for minimal distance from original score importance
	float sort_params_b; ///< weight for variance in prediction using imputation. the rest is change from prev
	float sort_params_k1; ///< weight for minimal distance from original score importance
	float sort_params_k2; ///< weight for variance in prediction using imputation. the rest is change from prev
	int max_set_size; ///< the size to look for to explain
	float override_score_bias; ///< when given will use it as score bias it train is very different from test
	float split_to_test; ///< to report RMSE on this ratio > 0 and < 1


	MissingShapExplainer();

	void _learn(const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	~MissingShapExplainer();

	ADD_CLASS_NAME(MissingShapExplainer)
		ADD_SERIALIZATION_FUNCS(retrain_predictor, max_test, missing_value, sample_masks_with_repeats,
			select_from_all, uniform_rand, use_shuffle, no_relearn, avg_bias_score, filters, processing, global_explain_params,
			predictor_type, predictor_args, max_weight, use_minimal_set, sort_params_a, sort_params_b,
			sort_params_k1, sort_params_k2, max_set_size, override_score_bias, verbose_apply, subsample_train,
			limit_mask_size, split_to_test)
};

/**
* shapley explainer with gibbs, GAN or other samples generator
*/
class ShapleyExplainer : public ModelExplainer {
private:
	unique_ptr<SamplesGenerator<float>> _sampler = NULL;
	void *sampler_sampling_args = NULL;

	GibbsSampler<float> _gibbs;
	GibbsSamplingParams _gibbs_sample_params;

	float avg_bias_score;

	void init_sampler(bool with_sampler = true);

	void _init(map<string, string> &mapper);
public:
	GeneratorType gen_type = GeneratorType::GIBBS; ///< generator type
	string generator_args = ""; ///< for learn
	string sampling_args = ""; ///< args for sampling
	int n_masks = 100; ///< how many test to conduct from shapley
	bool use_random_sampling = true; ///< If True will use random sampling - otherwise will sample mask size and than create it
	float missing_value = MED_MAT_MISSING_VALUE; ///< missing value

	ShapleyExplainer() { processor_type = FTR_POSTPROCESS_SHAPLEY; avg_bias_score = 0; }

	void _learn(const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	void post_deserialization();

	void load_GIBBS(MedPredictor *original_pred, const GibbsSampler<float> &gibbs, const GibbsSamplingParams &sampling_args);
	void load_GAN(MedPredictor *original_pred, const string &gan_path);
	void load_MISSING(MedPredictor *original_pred);
	void load_sampler(MedPredictor *original_pred, unique_ptr<SamplesGenerator<float>> &&generator);

	void dprint(const string &pref) const;

	ADD_CLASS_NAME(ShapleyExplainer)
		ADD_SERIALIZATION_FUNCS(_sampler, gen_type, generator_args, n_masks, missing_value, sampling_args,
			use_random_sampling, avg_bias_score, filters, processing, global_explain_params)
};

/**
* shapley-Lime explainer with gibbs, GAN or other sampler generator
*/
class LimeExplainer : public ModelExplainer {
private:
	unique_ptr<SamplesGenerator<float>> _sampler = NULL;
	void *sampler_sampling_args = NULL;

	//just for gibbs memory hold when init & learn
	GibbsSampler<float> _gibbs;
	GibbsSamplingParams _gibbs_sample_params;

	void init_sampler(bool with_sampler = true);
	void _init(map<string, string> &mapper);
	medial::shapley::LimeWeightMethod get_weight_method(string method_s);
public:
	GeneratorType gen_type = GeneratorType::GIBBS; ///< generator type
	string generator_args = ""; ///< for learn
	string sampling_args = ""; ///< args for sampling
	float missing_value = MED_MAT_MISSING_VALUE; ///< missing value
	float p_mask = 0; ///< prob for 1 in mask, if 0 - mask generation done by first selecting # of 1's in mask (uniformly) and then selecting the 1's
	medial::shapley::LimeWeightMethod weighting = medial::shapley::LimeWeightSum;
	int n_masks = 1250; ///< number of masks

	LimeExplainer() { processor_type = FTR_POSTPROCESS_LIME_SHAP; }

	void _learn(const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	void load_GIBBS(MedPredictor *original_pred, const GibbsSampler<float> &gibbs, const GibbsSamplingParams &sampling_args);
	void load_GAN(MedPredictor *original_pred, const string &gan_path);
	void load_MISSING(MedPredictor *original_pred);
	void load_sampler(MedPredictor *original_pred, unique_ptr<SamplesGenerator<float>> &&generator);

	void post_deserialization();

	void dprint(const string &pref) const;

	ADD_CLASS_NAME(LimeExplainer)
		ADD_SERIALIZATION_FUNCS(_sampler, gen_type, generator_args, missing_value, sampling_args, p_mask, n_masks, weighting,
			filters, processing, global_explain_params)
};

/**
* KNN explainer
*/
class KNN_Explainer : public ModelExplainer {
private:
	MedFeatures trainingMap;
	vector<float> average, std;

	// do the calculation for a single sample after normalization
	void computeExplanation(vector<float> thisRow, map<string, float> &sample_explain_reasons, vector <vector<int>> knnGroups, vector<string> knnGroupNames)const;

	void _init(map<string, string> &mapper);
public:

	int numClusters = -1; ///< how many samples (randomly chosen) represent the training space  -1:all. If larger than size of matrix, size of matrix will be used and warning generated.
	float fraction = (float)0.02; ///<fraction of points that is considered neighborhood to a point
	float chosenThreshold = MED_MAT_MISSING_VALUE; ///< Threshold to use on scores. If missing use thresholdQ to define threshold
	float thresholdQ = MED_MAT_MISSING_VALUE;///< defines threshold by positive ratio  on training set  ( when chosenThreshold missing). If this one is missing too, no thresholding. Explain by raw scoes.

	KNN_Explainer() { processor_type = FTR_POSTPROCESS_KNN_EXPLAIN; }

	void _learn(const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	ADD_CLASS_NAME(KNN_Explainer)
		ADD_SERIALIZATION_FUNCS(numClusters, trainingMap, average, std, fraction, chosenThreshold, filters, processing, global_explain_params)
};

/**
* Simple Linear Explainer - puts zeros for each feature and measures change in score
*/
class LinearExplainer : public ModelExplainer {
private:
	void _init(map<string, string> &mapper);

	float avg_bias_score;
public:
	LinearExplainer() { processor_type = FTR_POSTPROCESS_LINEAR; avg_bias_score = 0; }

	void _learn(const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	ADD_CLASS_NAME(LinearExplainer)
		ADD_SERIALIZATION_FUNCS(avg_bias_score, filters, processing, global_explain_params)
};

/**
* iterative set explainer with (gibbs, GAN or other samples generator) or proxy predictor algorithm
* to get as close as we can to final prediction score with lowest variance with the smallest set as possible of varaibles
*/
class IterativeSetExplainer : public ModelExplainer {
private:
	unique_ptr<SamplesGenerator<float>> _sampler = NULL;
	void *sampler_sampling_args = NULL;

	GibbsSampler<float> _gibbs;
	GibbsSamplingParams _gibbs_sample_params;

	float avg_bias_score;

	void init_sampler(bool with_sampler = true);

	void _init(map<string, string> &mapper);
public:
	GeneratorType gen_type = GeneratorType::GIBBS; ///< generator type
	string generator_args = ""; ///< for learn
	string sampling_args = ""; ///< args for sampling
	int n_masks = 100; ///< how many test to conduct from shapley
	bool use_random_sampling = true; ///< If True will use random sampling - otherwise will sample mask size and than create it
	float missing_value = MED_MAT_MISSING_VALUE; ///< missing value

	float sort_params_a; ///< weight for minimal distance from original score importance
	float sort_params_b; ///< weight for variance in prediction using imputation. the rest is change from prev
	float sort_params_k1; ///< weight for minimal distance from original score importance
	float sort_params_k2; ///< weight for variance in prediction using imputation. the rest is change from prev
	int max_set_size; ///< the size to look for to explain

	IterativeSetExplainer() { processor_type = FTR_POSTPROCESS_ITERATIVE_SET; avg_bias_score = 0; }

	void _learn(const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	void post_deserialization();

	void load_GIBBS(MedPredictor *original_pred, const GibbsSampler<float> &gibbs, const GibbsSamplingParams &sampling_args);
	void load_GAN(MedPredictor *original_pred, const string &gan_path);
	void load_MISSING(MedPredictor *original_pred);
	void load_sampler(MedPredictor *original_pred, unique_ptr<SamplesGenerator<float>> &&generator);

	void dprint(const string &pref) const;

	ADD_CLASS_NAME(IterativeSetExplainer)
		ADD_SERIALIZATION_FUNCS(_sampler, gen_type, generator_args, n_masks, missing_value, sampling_args,
			use_random_sampling, avg_bias_score, filters, processing, global_explain_params, max_set_size,
			sort_params_a, sort_params_b, sort_params_k1, sort_params_k2)
};

MEDSERIALIZE_SUPPORT(ExplainFilters)
MEDSERIALIZE_SUPPORT(ExplainProcessings)
MEDSERIALIZE_SUPPORT(GlobalExplainerParams)
MEDSERIALIZE_SUPPORT(TreeExplainer)
MEDSERIALIZE_SUPPORT(MissingShapExplainer)
MEDSERIALIZE_SUPPORT(ShapleyExplainer)
MEDSERIALIZE_SUPPORT(LimeExplainer)
MEDSERIALIZE_SUPPORT(LinearExplainer)
MEDSERIALIZE_SUPPORT(KNN_Explainer)
MEDSERIALIZE_SUPPORT(IterativeSetExplainer)



#endif
