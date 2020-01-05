#ifndef _FTR_PREDICTOR_IMPUTER_H_
#define _FTR_PREDICTOR_IMPUTER_H_

#include <string>
#include <random>
#include <MedAlgo/MedAlgo/BinSplitOptimizer.h>
#include <MedProcessTools/MedProcessTools/Calibration.h>

using namespace std;

/**
* Parameters fo Predictor Imputer
*/
class Predictor_Imputer_Params : public SerializableObject {
public:
	string predictor_type; ///< predictor args for multi-class
	string predictor_args; ///< predictor args for multi-class
	string num_class_setup; ///< param to control number of classes if needed in predictor
	BinSettings bin_settings; ///< binning method for each signal
	int sub_sample; ///< if given > 0 will sun_sample training to increase speed of train

	float calibration_save_ratio; ///< if given will use calibrate each prediction score on the saved_ratio. [0, 1]
	string calibration_string; ///< if calibration_save_ratio > 0 will use this init for calibration string

	Predictor_Imputer_Params();

	int init(map<string, string>& map);

	ADD_CLASS_NAME(Predictor_Imputer_Params)
		ADD_SERIALIZATION_FUNCS(predictor_type, predictor_args,
			num_class_setup, bin_settings, calibration_save_ratio, calibration_string, sub_sample)
};

/**
* Predictor Imputer - use all features in the matrix to predict value to impute
* selects randomly a value based on probability to get that value (similar to our gibbs)
*/
class PredictorImputer : public FeatureProcessor {
private:
	MedPredictor * predictor;
	vector<float> sorted_uniq_vals;
	vector<float> sorted_bin_vals;
	vector<Calibrator> calibrators; ///< calibrator for probability for each pred
	vector<string> predictor_features;
	int num_classes;
	mt19937 gen;
public:
	float missing_value; ///< missing value to look for to impute
	Predictor_Imputer_Params params; ///< parameters for the predictor
	bool verbose_learn; ///< if true will output more info when learning
	bool find_real_value; ///< if true will round to most similar origianl value
	bool debug; ///< if true will output verbose output in apply

	PredictorImputer() : FeatureProcessor() { init_defaults(); }

	~PredictorImputer();
	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<PredictorImputer *>(processor)); }

	void init_defaults();

	/// The parsed fields from init command.
	/// @snippet PredictorImputer.cpp PredictorImputer::init
	int init(map<string, string>& mapper);

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int _apply(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	ADD_CLASS_NAME(PredictorImputer)
		ADD_SERIALIZATION_FUNCS(processor_type, feature_name, resolved_feature_name, missing_value, params,
			find_real_value, predictor, predictor_features, sorted_bin_vals, sorted_uniq_vals, calibrators, num_classes,
			debug)
};

MEDSERIALIZE_SUPPORT(Predictor_Imputer_Params)
MEDSERIALIZE_SUPPORT(PredictorImputer)


#endif