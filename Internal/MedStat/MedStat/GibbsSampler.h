#ifndef __GIBBS_SAMPLER_H__
#define __GIBBS_SAMPLER_H__
#include <vector>
#include <string>
#include <map>
#include <random>
#include <unordered_map>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedAlgo/MedAlgo/BinSplitOptimizer.h>
#include "MedProcessTools/MedProcessTools/Calibration.h"

using namespace std;

/**
* A wrapper class to store same predictor trained on random selected samples to return prediction dist
*/
template<typename T> class PredictorOrEmpty : public SerializableObject {
public:
	int input_size;
	vector<T> sample_cohort; ///< all data points of feature
	MedPredictor *predictor; ///< predictors for each feature and probability to see Y (logloss function)
	vector<Calibrator> calibrators; ///< calibrator for probability for each pred
	vector<float> bin_vals; ///< the value of feature for each pred 

	vector<float> cluster_centers; ///< for kMeans centers
	vector<vector<float>> clusters_y; ///< for kMeans centers

	PredictorOrEmpty();
	~PredictorOrEmpty();

	/// retrieves random sample for feature based on all other features
	T get_sample(vector<T> &x, mt19937 &gen);

	ADD_CLASS_NAME(PredictorOrEmpty)
		ADD_SERIALIZATION_FUNCS(input_size, sample_cohort, cluster_centers, clusters_y, predictor, calibrators, bin_vals)
};

/**
* Parameters fo Gibbs Sampling
*/
class Gibbs_Params : public SerializableObject {
public:
	//learn args - kmeans or predictor
	int kmeans; ///< If > 0 will use kmeans to find clusters and look on each cluster y distribution - select 1 randomly and learn that
	float selection_ratio; ///< selection_ratio for kMeans - down sample
	int selection_count; ///< selection down sample count
	bool select_with_repeats; ///< If true will selct with repeats
	int max_iters; ///< max_iters for kmeans

	string predictor_type; ///< predictor args for multi-class
	string predictor_args; ///< predictor args for multi-class
	string num_class_setup; ///< param to control number of classes if needed in predictor
	BinSettings bin_settings; ///< binning method for each signal

	float calibration_save_ratio; ///< if given will use calibrate each prediction score on the saved_ratio. [0, 1]
	string calibration_string; ///< if calibration_save_ratio > 0 will use this init for calibration string
	
	Gibbs_Params();

	int init(map<string, string>& map);

	ADD_CLASS_NAME(Gibbs_Params)
		ADD_SERIALIZATION_FUNCS(kmeans, selection_ratio, max_iters, select_with_repeats, predictor_type, predictor_args,
			num_class_setup, bin_settings, calibration_save_ratio, calibration_string)
};

/**
* A class that contains all sampling arguments
*/
class GibbsSamplingParams : public SerializableObject {
public:
	//sample args
	int burn_in_count; ///< how many rounds in the start to ignore
	int jump_between_samples; ///< how many rounds to ignore between taking samples
	int samples_count; ///< how many samples to output
	bool find_real_value_bin; ///< If true will find closet real value to result - to be in same resolution, real value from train

	GibbsSamplingParams();
	int init(map<string, string>& map);

	ADD_CLASS_NAME(GibbsSamplingParams)
	ADD_SERIALIZATION_FUNCS(burn_in_count, jump_between_samples, samples_count, find_real_value_bin)
};

/**
* A gibbs sampler - has learn and create sample based on mask
*/
template<typename T> class GibbsSampler : public SerializableObject {
private:
	mt19937 _gen;
	bool done_prepare;
public:
	Gibbs_Params params; ///< gibbs params
	vector<PredictorOrEmpty<T>> feats_predictors; ///< gibbs_feature generators based on predictors
	vector<string> all_feat_names; ///< all features names (saved in learn)
	vector<string> impute_feat_names; ///< all features names (saved in learn)
	vector<vector<T>> uniqu_value_bins; ///< to round samples to those resoultions! - important for no leak!

	GibbsSampler();

	/// <summary>
	/// learn gibbs sample - for each feature creates predictors
	/// </summary>
	void learn_gibbs(const map<string, vector<T>> &cohort_data);

	/// <summary>
	/// learn gibbs sample - for each feature creates predictors
	/// </summary>
	void learn_gibbs(const map<string, vector<T>> &cohort_data, const vector<string> &learn_features, bool skip_missing);

	/// <summary>
	/// Should be called before first get_samples when used in parallel manner
	/// </summary>
	void prepare_predictors();

	/// <summary>
	/// generates samples based on gibbs sampling process
	/// </summary>
	void get_samples(map<string, vector<T>> &results, const GibbsSamplingParams &sampling_params,
		const vector<bool> *mask = NULL, const vector<T> *mask_values = NULL, bool print_progress = false);

	/// <summary>
	/// generates samples based on gibbs sampling process. const and can be called parallel
	/// </summary>
	void get_samples(map<string, vector<T>> &results, const GibbsSamplingParams &sampling_params, mt19937 &rnd_gen,
		const vector<bool> *mask = NULL, const vector<T> *mask_values = NULL, bool print_progress = false);

	/// <summary>
	/// generates samples based on gibbs sampling process - uses only burn rate and creates one sample and exits
	/// </summary>
	void get_parallel_samples(map<string, vector<T>> &results, const GibbsSamplingParams &sampling_params,
		const vector<bool> *mask = NULL, const vector<T> *mask_values = NULL);

	/// <summary>
	/// takes original cohort and results samples - filters and keep only samples that are similar to original population
	/// </summary>
	void filter_samples(const map<string, vector<float>> &cohort_data,
		map<string, vector<T>> &results, const string &predictor_type, const string &predictor_args,
		float filter_sens);

	int init(map<string, string>& map); ///< initialized params init function. reffer to that

	virtual ~GibbsSampler();

	ADD_CLASS_NAME(GibbsSampler<T>)
		ADD_SERIALIZATION_FUNCS(params, feats_predictors, uniqu_value_bins, all_feat_names, impute_feat_names)
};

MEDSERIALIZE_SUPPORT(PredictorOrEmpty<float>)
MEDSERIALIZE_SUPPORT(PredictorOrEmpty<double>)
MEDSERIALIZE_SUPPORT(Gibbs_Params)
MEDSERIALIZE_SUPPORT(GibbsSamplingParams)
MEDSERIALIZE_SUPPORT(GibbsSampler<float>)
MEDSERIALIZE_SUPPORT(GibbsSampler<double>)

#endif
