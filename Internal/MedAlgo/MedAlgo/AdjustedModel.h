#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <Logger/Logger/Logger.h>

class priors_info : public SerializableObject{
public:
	vector<string> names;
	vector<int> colIds;
	vector<int> min,max,factors;
	vector<float> probs;
	
	ADD_SERIALIZATION_FUNCS(names,colIds,min,max,factors,probs)
};

MEDSERIALIZE_SUPPORT(priors_info)

class MedAdjustedModel : public MedPredictor {
public:
	// Model
	MedPredictor *predictor;
	vector<string> onlyPriorsFeatures; ///< Features that are used only for prior adjustment;
	unordered_set<string> resolvedOnlyPriorsFeatures;
	priors_info priors;
	int n_preds = 1;
	float odds=-1; ///< over all odds. learn if not given

	/// Parameters
	string predictor_name;
	string predictor_params;
	string priorsFile;



	// Function
	MedAdjustedModel() { classifier_type = MODEL_ADJUSTED; };

	~MedAdjustedModel() {};

	int set_params(map<string, string>& mapper);

	void init_defaults() {};

	// learn - we only have the medmat option
	int learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts);
	int Learn(float *x, float *y, int nsamples, int nftrs) { HMTHROW_AND_ERR("MedAdjustedModel: Learn(float *,...) not implemented, used the MedFeatures API instead\n"); };
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs) { HMTHROW_AND_ERR("MedAdjustedModel: Learn(float *,...) not implemented, used the MedFeatures API instead\n"); };

	// predict - we only have the medmat option
	int predict(MedMat<float> &x, vector<float> &preds) const;
	int Predict(float *x, float *&preds, int nsamples, int nftrs) { HMTHROW_AND_ERR("MedAdjustedModel: Predict(float *,...) not implemented, used the MedFeatures API instead\n"); };
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) { HMTHROW_AND_ERR("MedAdjustedModel: Predict(float *,...) not MedFeatures, used the MedMat API instead\n"); };

	int n_preds_per_sample() { return n_preds; }

	// helpers
	void readPriors();
	void getOdds(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts);
	void getPosteriors(vector<float>& preds, MedMat<float>& priorsX) const;
	void getSubMatrices(MedMat<float> &x, MedMat<float> &priorsX, MedMat<float> &predictorX) const;

	ADD_CLASS_NAME(MedAdjustedModel)
	ADD_SERIALIZATION_FUNCS(classifier_type, predictor, onlyPriorsFeatures, priors, n_preds, priorsFile,odds)
};

MEDSERIALIZE_SUPPORT(MedAdjustedModel)