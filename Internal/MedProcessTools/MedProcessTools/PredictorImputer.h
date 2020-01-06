#ifndef _FTR_PREDICTOR_IMPUTER_H_
#define _FTR_PREDICTOR_IMPUTER_H_

#include <string>
#include <random>
#include <MedAlgo/MedAlgo/BinSplitOptimizer.h>
#include <MedProcessTools/MedProcessTools/Calibration.h>
#include <MedAlgo/MedAlgo/SamplesGenerator.h>

using namespace std;

/**
* Predictor Imputer - use all features in the matrix to predict value to impute
* selects randomly a value based on probability to get that value (similar to our gibbs)
*/
class PredictorImputer : public FeatureProcessor {
private:
	unique_ptr<SamplesGenerator<float>> _sampler = NULL;
	void *sampler_sampling_args = NULL;

	GibbsSampler<float> _gibbs;
	GibbsSamplingParams _gibbs_sample_params;
	int n_masks = 1;
	vector<string> impute_features;

	void init_sampler(bool with_sampler = true);

	mt19937 gen;
public:
	float missing_value; ///< missing value to look for to impute
	bool verbose_learn; ///< if true will output more info when learning
	bool verbose_apply; ///< if true will output verbose output in apply
	string tag_search; ///< feature tag search

	GeneratorType gen_type; ///< generator type
	string generator_args; ///< for learn
	string sampling_args; ///< args for sampling

	PredictorImputer() : FeatureProcessor() { init_defaults(); }

	// Copy
	//void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<PredictorImputer *>(processor)); }

	void init_defaults();

	void post_deserialization();

	void load_GIBBS(const GibbsSampler<float> &gibbs, const GibbsSamplingParams &sampling_args);
	void load_GAN(const string &gan_path);
	void load_MISSING();
	void load_sampler(unique_ptr<SamplesGenerator<float>> &&generator);

	/// The parsed fields from init command.
	/// @snippet PredictorImputer.cpp PredictorImputer::init
	int init(map<string, string>& mapper);

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int _apply(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	ADD_CLASS_NAME(PredictorImputer)
		ADD_SERIALIZATION_FUNCS(processor_type, tag_search, missing_value, gen_type, _sampler,
			generator_args, sampling_args, verbose_apply, impute_features)
};

MEDSERIALIZE_SUPPORT(PredictorImputer)


#endif