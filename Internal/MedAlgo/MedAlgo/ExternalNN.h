//
// ExternalNN predictor is used in cases we train the model in external code such as tensorflow or keras,
// and then wish to use it as a predictor. We use the ApplyKeras mechanism for that.
// 

//#pragma once
#ifndef __EXTERNAL_NN_H__
#define __EXTERNAL_NN_H__

#include <Logger/Logger/Logger.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedEmbed/MedEmbed/ApplyKeras.h>


class MedExternalNN : public MedPredictor {
public:
	// Model
	ApplyKeras ak;
	int n_preds = 1;

	/// Parameters
	string init_file = "";

	// Function
	MedExternalNN() { classifier_type = MODEL_EXTERNAL_NN; };

	~MedExternalNN() {};

	int set_params(map<string, string>& mapper) {
		for (auto &e : mapper)
			if (e.first == "init_file") {
				init_file = e.second;
				if (ak.init_from_text_file(init_file) < 0)
					HMTHROW_AND_ERR("ERROR: Failed reading layers file %s\n", init_file.c_str());
				n_preds = ak.get_output_dimension();
				//fprintf(stderr, " ===> n_preds is %d\n", n_preds);
			}
		return 0;
	}

	void init_defaults() { init_file = ""; };

	// the following simply initializes 'ak' and 'n_preds' from init_file
	int external_nn_learn() {

		return 0;
	};

	// learn simply calls init from file
	int learn(const MedFeatures& features) { return external_nn_learn(); }
	int learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts) { return external_nn_learn(); }
	int Learn(float *x, float *y, int nsamples, int nftrs) { return external_nn_learn(); }
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs) { return external_nn_learn();  };

	// predict - we only have the medmat option
	int Predict(MedMat<float> &x, vector<float> &preds)	{

		MedMat<float> res;
		ak.apply(x, res);
		res.copy_vec(preds);

		return 0;
	}

	// following are not implemented and will simply HMTHROW_AND_ERR

	int Predict(float *x, float *&preds, int nsamples, int nftrs) { HMTHROW_AND_ERR("ExternalNN: Predict(float *,...) not implemented, used the MedMat API instead\n"); };
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) { HMTHROW_AND_ERR("ExternalNN: Predict(float *,...) not implemented, used the MedMat API instead\n"); }; 

	int n_preds_per_sample() const { return n_preds; } 

	bool predict_single_not_implemented() { return true; }

	ADD_CLASS_NAME(MedExternalNN)
	ADD_SERIALIZATION_FUNCS(classifier_type, n_preds, ak)
};

MEDSERIALIZE_SUPPORT(MedExternalNN)


#endif