#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <Logger/Logger/Logger.h>


class MedSimpleEnsemble : public MedPredictor {
public:
	// Model
	vector<MedPredictor *> predictors;
	MedPredictor *combiner = NULL;
	vector<float>	weights;
	int n_preds = 1;
	float p_combine = 0;

	/// Parameters
	vector<string> predictor_names;
	vector<string> predictor_params;
	string combiner_name = "";
	string combiner_params = "";


	// Function
	MedSimpleEnsemble() { classifier_type = MODEL_SIMPLE_ENSEMBLE; };

	~MedSimpleEnsemble() {};

	int set_params(map<string, string>& mapper);

	void init_defaults() { predictors.clear(); weights.clear(); predictor_names.clear(); predictor_params.clear(); combiner_name = ""; combiner_params = ""; p_combine = 0; };

	// learn simply calls init from file
	int learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts);
	int Learn(float *x, float *y, int nsamples, int nftrs) { HMTHROW_AND_ERR("SimpleEnsemble: Learn(float *,...) not implemented, used the MedMat API instead\n"); };
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs) { HMTHROW_AND_ERR("SimpleEnsemble: Learn(float *,...) not implemented, used the MedMat API instead\n"); };

	// predict - we only have the medmat option
	int predict_pre_combine(MedMat<float> &x, MedMat<float> &preds) const;
	int predict(MedMat<float> &x, vector<float> &preds) const;	
	int Predict(float *x, float *&preds, int nsamples, int nftrs) { HMTHROW_AND_ERR("SimpleEnsemble: Predict(float *,...) not implemented, used the MedMat API instead\n"); };
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) { HMTHROW_AND_ERR("SimpleEnsemble: Predict(float *,...) not implemented, used the MedMat API instead\n"); };

	int n_preds_per_sample() { return n_preds; }

	ADD_CLASS_NAME(MedSimpleEnsemble)
	ADD_SERIALIZATION_FUNCS(classifier_type, predictors, weights, predictor_names, predictor_params, combiner, combiner_name, combiner_params, p_combine, n_preds)
};

MEDSERIALIZE_SUPPORT(MedSimpleEnsemble)

