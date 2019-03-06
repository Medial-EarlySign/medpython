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
	vector<MedPredictor *> predictors;
	vector<float> sample_cohort;

	PredictorOrEmpty();

	float get_sample(vector<float> &x);
};

/**
* Parameters fo Gibbs Sampling
*/
class Gibbs_Params
	: public SerializableObject {
public:
	int burn_in_count; //how many rounds in the start to ignore
	int jump_between_samples; // how many rounds to ignore between taking samples
	int samples_count; //how many samples to output
	string predictor_type;
	string predictor_args;
	int predictors_counts;
	float selection_ratio;

	Gibbs_Params();

	int init(map<string, string>& map);
};

namespace medial {

	/*!
	*  \brief statsitical namespace
	*/
	namespace stats {
		/// \brief gibb_sampling process on cohort
		void gibbs_sampling(const map<string, vector<float>> &cohort_data, const Gibbs_Params &params,
			map<string, vector<float>> &results, const vector<bool> *mask = NULL, const vector<float> *mask_values = NULL);


	}
}



#endif
