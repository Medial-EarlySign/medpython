//
// MedBooster.h
//
// A general wrapper allowing several boosting strategies for predictors.
//

#ifndef __MEDBOOSTER__H__
#define __MEDBOOSTER__H__

#include "MedAlgo.h"

//============================================================================================================
// currently parameters should be initialized from a string
// The format for that is :
// boost#booster_type=...,loss_func=...,nrounds=...,shrinkage=....,internal=...#internal_params#....
// example:
// boost@booster_type=gradient,loss_func=square,nrounds=40,shrinkage=0.05,internal_model=linear_model@internal_params@rfactor=0.8
//
class MedBoosterParams {
	public:

		string init_string;

		int nrounds;
		float shrinkage;

		string booster_type;		// adaboost OR gradient
		string loss_func;			// one of: square , softmax , logistic

		string internal_booster_model;
		string internal_booster_params;

		int print_flag;				// weather to print something in each round in training

		void init_from_string(const string &init_str);

};

//============================================================================================================
class MedBooster : public MedPredictor {

	public:

		MedBoosterParams params;

		vector<MedPredictor *> predictors;		// actual predictors for each round.
		MedMat<float> alpha;					// linear coefficient of each predictor (could be multi channel), last one is always the bias

		MedBooster() { classifier_type = MODEL_BOOSTER; init_defaults(); }
		~MedBooster() {
			for (auto &predictor : predictors) if (predictor) { delete predictor; predictor = NULL; }
		}

		MedBooster(void *_params) { params = *(MedBoosterParams *)_params; }
		MedBooster(MedBoosterParams& _params) { params = _params; }
		int init(void *_params) { params = *(MedBoosterParams *)_params; return 0; }
		
		/// The parsed fields from init command.
		/// @snippet MedBooster.cpp MedBoosterParams::init_from_string
		int init_from_string(string initialization_text) {
			cerr << "MedBooster init_from_string ! :: " << initialization_text << "\n";
			params.init_from_string(initialization_text);
			return 0;
		}
		///MedBooster:: init map :: not supported, only init_from_string supported 
		int init(map<string, string>& mapper) {
			cerr << "MedBooster:: init map :: not supported, only init_from_string supported....\n";
			return -1;
		}

		int set_params(map<string, string>& mapper) {
			cerr << "MedBooster:: init map :: not supported, only init_from_string supported....\n";
			return -1;
		}

		//	int init(const string &init_str); // allows init of parameters from a string. Format is: param=val,... , for sampsize: 0 is NULL, a list of values is separated by ; (and not ,)
		void init_defaults() { params.init_from_string("boost@booster_type=gradient,loss_func=square,nrounds=10,shrinkage=0.1,internal_model=linear_model@internal_params@rfactor=0.8"); }

		int Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
			cerr << "MedBooster:: Learn :: API's with MedMat are preferred....\n";
			MedMat<float> xmat; xmat.load(x, nsamples, nftrs);
			MedMat<float> ymat; ymat.load(y, nsamples, 1);
			return learn(xmat, ymat);
		}

		int Predict(float *x, float *&preds, int nsamples, int nftrs) const {
			cerr << "MedBooster:: Learn :: API's with MedMat are preferred....\n";
			MedMat<float> xmat; xmat.load(x, nsamples, nftrs);
			vector<float> vpreds;
			int rc = predict(xmat, vpreds);
			if (preds == NULL) preds = new float[nsamples];
			memcpy(preds, &vpreds[0], sizeof(float)*nsamples);
			return rc;
		}

		int learn(MedMat<float> &x, MedMat<float> &y);
		int learn(MedMat<float> &x, MedMat<float> &y, vector<float> &wgt) { return learn(x, y); }
		int predict(MedMat<float> &x, vector<float> &preds) const;

		// serializations
		ADD_CLASS_NAME(MedBooster)
		size_t get_size();
		size_t serialize(unsigned char *blob);
		size_t deserialize(unsigned char *blob);


		// internal funcs
		int learn_gradient_booster(MedMat<float> &x, MedMat<float> &y);
		MedPredictor *train_internal_model(MedMat<float> &x, MedMat<float> &y, string model_name, string model_params);
		int gradient_get_round_alphas(MedMat<float> &residual, vector<float> &preds, vector<float> &a);
		int gradient_get_alphas_square_loss(MedMat<float> &residual, vector<float> &preds, vector<float> &a);
		int update_round_residuals(MedMat<float> &y, vector<float> &all_preds, MedMat<float> &residual, vector<float> &preds, vector<float> &a);
		int update_round_residuals_single_linear(MedMat<float> &y, vector<float> &all_preds, MedMat<float> &residual, vector<float> &preds, vector<float> &a);
};

//=================================================================
// Joining the MedSerialize Wagon
//=================================================================

MEDSERIALIZE_SUPPORT(MedBooster)
#endif