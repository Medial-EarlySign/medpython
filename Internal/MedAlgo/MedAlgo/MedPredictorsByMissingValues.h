#ifndef __MEDPREDICTOR_BY_MISSING_VALUE_H__
#define __MEDPREDICTOR_BY_MISSING_VALUE_H__

#include "MedAlgo.h"
#include <MedProcessTools/MedProcessTools/Calibration.h>

class MedPredictorsByMissingValues : public MedPredictor {
public:
	vector<MedPredictor *> predictors;
	vector<Calibrator *> calibrators;
	string predictor_type;
	string predictor_params;
	string masks_params;
	string masks_tw;
	string calibrate_predictions;

	/// <summary>
	/// an initialization for model
	/// @snippet MedPredictorsByMissingValues.cpp MedPredictorsByMissingValues::init
	/// </summary>
	int init(map<string, string>& mapper);
	//int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	vector<int> get_cols_of_mask(string &full_mask, vector<string> &signals) const;
	vector<int> get_cols_of_predictor(int i, vector<vector<int>> &cols_per_mask, vector<int> &all_cols) const;
	int learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts);
	int predict(MedFeatures& features) const;

	~MedPredictorsByMissingValues() {
		for (size_t i = 0; i < predictors.size(); ++i)
		{
			delete predictors[i];
			predictors[i] = NULL;
		}
		predictors.clear();
	}

	MedPredictorsByMissingValues() {
		classifier_type = MODEL_BY_MISSING_VALUES_SUBSET;
	}

	ADD_CLASS_NAME(MedPredictorsByMissingValues)
		ADD_SERIALIZATION_FUNCS(classifier_type, predictors, predictor_type, predictor_params, masks_params, masks_tw, calibrate_predictions, calibrators)
};

MEDSERIALIZE_SUPPORT(MedPredictorsByMissingValues)

#endif