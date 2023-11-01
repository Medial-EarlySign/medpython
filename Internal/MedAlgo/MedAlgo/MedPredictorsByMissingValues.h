#ifndef __MEDPREDICTOR_BY_MISSING_VALUE_H__
#define __MEDPREDICTOR_BY_MISSING_VALUE_H__

#include "MedAlgo.h"

class MedPredictorsByMissingValues : public MedPredictor {
public:
	vector<MedPredictor *> predictors;
	string predictor_type;
	string predictor_params;
	string masks_params;

	/// <summary>
	/// an initialization for model
	/// @snippet MedPredictorsByMissingValues.cpp MedPredictorsByMissingValues::init
	/// </summary>
	int init(map<string, string>& mapper);
	//int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts);
	int predict(MedMat<float> &x, vector<float> &preds) const;

	~MedPredictorsByMissingValues() {
		for (size_t i = 0; i < predictors.size(); ++i)
		{
			delete predictors[i];
			predictors[i] = NULL;
		}
		predictors.clear();
	}

	ADD_CLASS_NAME(MedPredictorsByMissingValues)
		ADD_SERIALIZATION_FUNCS(classifier_type, predictors, predictor_type, predictor_params)
};

MEDSERIALIZE_SUPPORT(MedPredictorsByMissingValues)

#endif