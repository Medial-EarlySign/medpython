#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>


//======================================================================================
// micNet: Home brewed implimenatation for Neural Nets and Deep Learning
//======================================================================================
struct MedMicNetParams {

	string init_string;
};

class MedMicNet : public MedPredictor {
private:
	vector<micNet> model_per_thread;
	bool is_prepared = false;
public:
	// Model 
	micNet mic;

	/// Parameters
	MedMicNetParams mic_params;

	// Function
	MedMicNet() { classifier_type = MODEL_MIC_NET; mic.params.init_defaults(); }
	MedMicNet(void *params) { mic_params = *(MedMicNetParams *)params; mic.init_from_string(mic_params.init_string); }
	MedMicNet(MedMicNetParams& params) { mic_params = params; mic.init_from_string(mic_params.init_string); }
	int init(void *params) { mic_params = *(MedMicNetParams *)params; return mic.init_from_string(mic_params.init_string); }
	/// The parsed fields from init command.
	/// @snippet micNet.cpp micNetParams::init_from_string
	int init_from_string(string initialization_text) {
		cerr << "MedMicNet init_from_string ! :: " << initialization_text << "\n";
		mic_params.init_string = initialization_text;
		cerr << "calling init_from_string of micNet\n"; fflush(stderr);
		return mic.init_from_string(initialization_text);
	}

	///MedMicNet:: init map :: not supported, only init_from_string supported 
	int init(map<string, string>& mapper) {
		cerr << "MedMicNet:: init map :: not supported, only init_from_string supported....\n";
		return -1;
	}
	int set_params(map<string, string>& mapper) {
		cerr << "MedMicNet:: init map :: not supported, only init_from_string supported....\n";
		return -1;
	}
	//	int init(const string &init_str); // allows init of parameters from a string. Format is: param=val,... , for sampsize: 0 is NULL, a list of values is separated by ; (and not ,)
	void init_defaults() { mic_params.init_string = ""; mic.params.init_defaults(); }

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
		cerr << "MedMicNet:: Learn :: API's with MedMat are preferred....\n";
		MedMat<float> xmat; xmat.load(x, nsamples, nftrs);
		MedMat<float> ymat; ymat.load(y, nsamples, 1);
		MedMat<float> wmat;
		if (w != NULL) wmat.load(w, nsamples, 1);
		return learn(xmat, ymat, wmat.get_vec());
	}

	int Predict(float *x, float *&preds, int nsamples, int nftrs) const {
		cerr << "MedMicNet:: Predict :: API's with MedMat are preferred....\n";
		MedMat<float> xmat; xmat.load(x, nsamples, nftrs);
		vector<float> vpreds;
		int rc = predict(xmat, vpreds);
		if (preds == NULL) preds = new float[nsamples];
		memcpy(preds, &vpreds[0], sizeof(float)*nsamples);
		return rc;
	}

	int learn(MedMat<float> &x, MedMat<float> &y, vector<float> &wgt) { return mic.learn(x, y, wgt); }
	int learn(MedMat<float> &x, MedMat<float> &y) { return mic.learn(x, y); }
	int predict(MedMat<float> &x, vector<float> &preds) const { micNet mutable_net = mic; return mutable_net.predict(x, preds); }

	void prepare_predict_single();
	void predict_single(const vector<float> &x, vector<float> &preds) const;

	// Predictions per sample
	int n_preds_per_sample() const { return mic.n_preds_per_sample(); }

	// (De)Serialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	ADD_CLASS_NAME(MedMicNet)
	ADD_SERIALIZATION_FUNCS(classifier_type, mic_params.init_string, mic)


};

MEDSERIALIZE_SUPPORT(MedMicNet)
