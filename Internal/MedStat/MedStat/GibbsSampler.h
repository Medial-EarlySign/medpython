#ifndef __GIBBS_SAMPLER_H__
#define __GIBBS_SAMPLER_H__
#include <vector>
#include <string>
#include <map>
#include <random>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>

using namespace std;

/**
* A wrapper class to store same predictor trained on random selected samples to return prediction dist
*/
class PredictorOrEmpty : public SerializableObject {
private:
	mt19937 gen;
public:
	int input_size;
	vector<float> sample_cohort; ///< all data points of feature

	vector<float> cluster_centers; ///< for kMeans centers
	vector<vector<float>> clusters_y; ///< for kMeans centers

	PredictorOrEmpty();

	/// retrieves random sample for feature based on all other features
	float get_sample(vector<float> &x);

	ADD_CLASS_NAME(PredictorOrEmpty)
		ADD_SERIALIZATION_FUNCS(input_size, sample_cohort, cluster_centers, clusters_y)
};

/**
* Parameters fo Gibbs Sampling
*/
class Gibbs_Params : public SerializableObject {
public:
	//learn args
	int kmeans; ///< If > 0 will use kmeans to find clusters and look on each cluster y distribution - select 1 randomly and learn that
	float selection_ratio; ///< selection_ratio for kMeans - down sample
	bool select_with_repeats; ///< If true will selct with repeats
	int max_iters; ///< max_iters for kmeans

	//sample args
	int burn_in_count; ///< how many rounds in the start to ignore
	int jump_between_samples; ///< how many rounds to ignore between taking samples
	int samples_count; ///< how many samples to output
	bool find_real_value_bin; ///< If true will find closet real value to result - to be in same resolution, real value from train

	Gibbs_Params();

	int init(map<string, string>& map);

	ADD_CLASS_NAME(Gibbs_Params)
		ADD_SERIALIZATION_FUNCS(kmeans, selection_ratio, max_iters, select_with_repeats, burn_in_count,
			jump_between_samples, samples_count, find_real_value_bin)
};

/**
* A gibbs sampler - has learn and create sample based on mask
*/
class GibbsSampler : public SerializableObject {
public:
	Gibbs_Params params; ///< gibbs params
	vector<PredictorOrEmpty> feats_predictors; ///< gibbs_feature generators based on predictors
	vector<string> all_feat_names; ///< all features names (saved in learn)
	vector<vector<float>> uniqu_value_bins; ///< to round samples to those resoultions! - important for no leak!

	/// <summary>
	/// learn gibbs sample - for each feature creates predictors
	/// </summary>
	void learn_gibbs(const map<string, vector<float>> &cohort_data);

	/// <summary>
	/// generates samples based on gibbs sampling process
	/// </summary>
	void get_samples(map<string, vector<float>> &results,
		const vector<bool> *mask = NULL, const vector<float> *mask_values = NULL);

	/// <summary>
	/// generates samples based on gibbs sampling process - uses only burn rate and creates one sample and exits
	/// </summary>
	void get_parallel_samples(map<string, vector<float>> &results, uniform_real_distribution<> &real_dist,
		const vector<bool> *mask = NULL);

	/// <summary>
	/// takes original cohort and results samples - filters and keep only samples that are similar to original population
	/// </summary>
	void filter_samples(const map<string, vector<float>> &cohort_data,
		map<string, vector<float>> &results, const string &predictor_type, const string &predictor_args,
		float filter_sens);

	int init(map<string, string>& map); ///< initialized params init function. reffer to that

	ADD_CLASS_NAME(GibbsSampler)
		ADD_SERIALIZATION_FUNCS(params, feats_predictors, uniqu_value_bins, all_feat_names)
};

MEDSERIALIZE_SUPPORT(PredictorOrEmpty)
MEDSERIALIZE_SUPPORT(Gibbs_Params)
MEDSERIALIZE_SUPPORT(GibbsSampler)

#endif
