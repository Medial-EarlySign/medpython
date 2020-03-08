#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>

//================================================================
// MultiClass
//================================================================
enum MedMultiClassType {
	MULTI_CLASS_ONE_VS_ALL = 1,
	MULTI_CLASS_LAST
};

struct MedMultiClassParams {

	MedPredictorTypes method;
	MedMultiClassType multi_class_type;

	vector<float> class_values;
	void *internal_params;

};


struct MedMultiClass : public MedPredictor {

	MedMultiClassParams params;
	vector<MedPredictor *> internal_predictors;

	// Function
	MedMultiClass();
	MedMultiClass(void *params);
	MedMultiClass(MedMultiClassParams& params);

	int init(void *params);
	void set_internal_method(MedPredictorTypes type);
	void init_defaults();
	~MedMultiClass();

	int init_classifiers();
	int init_classifier(int index);

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	// Print
	void print(FILE *fp, const string& prefix, int level = 0) const;

	// Predictions per sample
	int n_preds_per_sample() const;

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	ADD_CLASS_NAME(MedMultiClass)
		ADD_SERIALIZATION_FUNCS(classifier_type, params.method, params.multi_class_type, internal_predictors)
};

MEDSERIALIZE_SUPPORT(MedMultiClass)
