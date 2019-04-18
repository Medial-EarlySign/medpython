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
* An abstract class API for explainer
*/
class ModelExplainer : public PostProcessor {
public:
	MedPredictor * original_predictor = NULL; //uses this if model has implementation of SHAP (like xgboost, lightGBM)

	/// overload function for ModelExplainer - easier API
	virtual void Learn(MedPredictor *original_pred, const MedFeatures &train_mat) = 0; 

	///Learns from predictor and train_matrix (PostProcessor API)
	virtual void Learn(MedModel &model, MedPidRepository& rep, const MedFeatures &train_mat);
	void Apply(MedFeatures &matrix) { explain(matrix); } //alias for explain

														 ///Virtual - return explain results in sample_feature_contrib
	virtual void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const = 0;

	/// Stores explain results in matrix
	virtual void explain(MedFeatures &matrix) const; //stores _explain results in MedFeatures

	static void print_explain(MedSample &smp);

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
* 1. Uses xgboost/lightGBM feature contirbution
* 2. Reads tree model into structure to calc SHAP values - QRF, BART..?
* 3. train LightGBM/Xgboost proxy to explain non-tree Predictor
*/
class TreeExplainer : public ModelExplainer {
private:
	MedPredictor * proxy_predictor = NULL; //uses this if model has no tree implementation
	TreeEnsemble generic_tree_model;
	//Tree structure of generic ensamble trees
private:
	bool try_convert_trees();
public:
	string proxy_model_type;
	string proxy_model_init;
	bool interaction_shap = false;
	int approximate = false;
	float missing_value = MED_MAT_MISSING_VALUE;

	int init(map<string, string> &mapper);

	TreeExplainerMode get_mode() const;

	void Learn(MedPredictor *original_pred, const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	void post_deserialization();

	~TreeExplainer();

	ADD_CLASS_NAME(TreeExplainer)
		ADD_SERIALIZATION_FUNCS(original_predictor, proxy_predictor, generic_tree_model, interaction_shap)
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
public:
	int add_new_data; ///< how many new data data points to add for train according to sample masks

	int max_test; ///< max number of samples in SHAP
	float missing_value; ///< missing value 
	bool sample_masks_with_repeats; ///< Whether or not to sample masks with repeats
	float select_from_all; ///< If max_test is beyond this percentage of all options than sample from all options (to speed up runtime)
	bool uniform_rand; ///< it True will sample masks uniformlly
	bool use_shuffle; ///< if not sampling uniformlly, If true will use shuffle (to speed up runtime)
	string change_learn_args; ///< arguments to change in predictor - for example to change it into regression
	bool verbose_learn; ///< If true will print more in learn

	MissingShapExplainer();

	int init(map<string, string> &mapper);

	void Learn(MedPredictor *original_pred, const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	~MissingShapExplainer();

	ADD_CLASS_NAME(MissingShapExplainer)
		ADD_SERIALIZATION_FUNCS(original_predictor, retrain_predictor, max_test, missing_value,
			sample_masks_with_repeats, select_from_all, uniform_rand, use_shuffle)
};

/// @enum
/// Generator Types options
enum GeneratorType
{
	GIBBS = 0, ///< "GIBBS" - to use GibbsSampler
	GAN = 1, ///< "GAN" to use GAN generator, accepts GAN path
	MISSING = 2 ///< "MISSING" to use no generator, just puts missing values where mask[i]==0
};

/// convert function for generator type to string
string GeneratorType_toStr(GeneratorType type);
/// convert function for generator
GeneratorType GeneratorType_fromStr(const string &type);

/**
* shapley explainer with gibbs, GAN or other sampler generator
*/
class ShapleyExplainer : public ModelExplainer {
private:
	unique_ptr<SamplesGenerator<float>> _sampler = NULL;
	void *sampler_sampling_args = NULL;

	GibbsSampler<float> _gibbs;
	GibbsSamplingParams _gibbs_sample_params;

	void init_sampler();
public:
	GeneratorType gen_type = GeneratorType::GIBBS; ///< generator type
	string generator_args = ""; ///< for learn
	string sampling_args = ""; ///< args for sampling
	int max_test = 100; ///< how many test to conduct from shapley
	float missing_value = MED_MAT_MISSING_VALUE; ///< missing value

	int init(map<string, string> &mapper);

	void Learn(MedPredictor *original_pred, const MedFeatures &train_mat);

	void explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const;

	void post_deserialization();

	void load_GIBBS(MedPredictor *original_pred, const GibbsSampler<float> &gibbs, const GibbsSamplingParams &sampling_args);
	void load_GAN(MedPredictor *original_pred, const string &gan_path);
	void load_MISSING(MedPredictor *original_pred);

	ADD_CLASS_NAME(ShapleyExplainer)
		ADD_SERIALIZATION_FUNCS(original_predictor, _sampler, gen_type, generator_args, max_test, missing_value, sampling_args)
};

MEDSERIALIZE_SUPPORT(TreeExplainer)
MEDSERIALIZE_SUPPORT(MissingShapExplainer)
MEDSERIALIZE_SUPPORT(ShapleyExplainer)

#endif
