#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>

//========================================================================================
// TQRF Wrapper
//========================================================================================
class MedTQRF : public MedPredictor {
public:
	TQRF_Forest _tqrf;

	MedTQRF() { classifier_type = MODEL_TQRF; }
	~MedTQRF() {};

	void init_defaults() {};

	// initialize using the init_from_string() method (inherited from SerializableObject)

	virtual int init(map<string, string>& mapper) { return _tqrf.init(mapper); }
	virtual int set_params(map<string, string>& mapper) { return _tqrf.init(mapper); }

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
		HMTHROW_AND_ERR("MedTQRF does not support the Learn(float *x, float *y, float *w, int nsamples, int nftrs). Use Learn(MedFeatures &feats) API instead\n");
	};
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const {
		HMTHROW_AND_ERR("MedTQRF does not support the Predict(float *x, float *&preds, int nsamples, int nftrs). Use Predict(MedMat<float> &x, vector<float> &preds) API instead\n");
	}

	int Learn(const MedFeatures &feats) { return _tqrf.Train(feats); }
	int Predict(MedMat<float> &x, vector<float> &preds) const { return _tqrf.Predict(x, preds); }

	ADD_CLASS_NAME(MedTQRF)
		ADD_SERIALIZATION_FUNCS(classifier_type, _tqrf)

};

MEDSERIALIZE_SUPPORT(MedTQRF)