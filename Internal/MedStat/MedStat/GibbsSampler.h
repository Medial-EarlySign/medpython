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
class PredictorOrEmpty {
private:
	mt19937 gen;
public:
	vector<MedPredictor *> predictors; ///< list of all predictors for feature
	vector<float> sample_cohort; ///< all data points of feature

	PredictorOrEmpty();

	/// retrieves random sample for feature based on all other features
	float get_sample(vector<float> &x);
};

/**
* Parameters fo Gibbs Sampling
*/
class Gibbs_Params : public SerializableObject {
public:
	//learn args
	string predictor_type; ///< predictor type to learn
	string predictor_args; ///< predictor arg to learn
	int predictors_counts; ///< how many random predictors for feature
	float selection_ratio; ///< which ratio of data to take for train in each predictor for randomness

	//sample args
	int burn_in_count; ///< how many rounds in the start to ignore
	int jump_between_samples; ///< how many rounds to ignore between taking samples
	int samples_count; ///< how many samples to output

	Gibbs_Params();

	int init(map<string, string>& map);
};

/**
* A gibbs sampler - has learn and create sample based on mask
*/
class GibbsSampler : public SerializableObject {
public:
	Gibbs_Params params; ///< gibbs params
	vector<PredictorOrEmpty> feats_predictors; ///< gibbs_feature generators based on predictors
	vector<string> all_feat_names; ///< all features names (saved in learn)

	/// <summary>
	/// learn gibbs sample - for each feature creates predictors
	/// </summary>
	void learn_gibbs(const map<string, vector<float>> &cohort_data);

	/// <summary>
	/// generates samples based on gibbs sampling process
	/// </summary>
	void get_samples(map<string, vector<float>> &results,
		const vector<bool> *mask = NULL, const vector<float> *mask_values = NULL);

	int init(map<string, string>& map); ///< initialized params init function. reffer to that

	ADD_CLASS_NAME(GibbsSampler)
		ADD_SERIALIZATION_FUNCS(params, feats_predictors, all_feat_names)
};

#endif
