#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>

//==============================================================================================
// Linear Models2: Linear regression (with Ridge and/or Lasso), using Gradient Descent variants
//==============================================================================================
struct MedGDLMParams : public SerializableObject {

	// Required params
	int max_iter;
	float stop_at_err; ///< stop criteria
	int max_times_err_grows;
	string method; ///< gd or sgd
	int batch_size;	///< for sgd
	float rate;
	float rate_decay;
	float momentum;

	int last_is_bias;
	int print_model;
	bool verbose_learn;

	// Optional params
	float l_ridge; ///< lambda for ridge
	vector<float> ls_ridge; ///< lambdas for ridge
	float l_lasso; ///< labmda for lasso
	vector<float> ls_lasso;; ///< labmdas for lasso

	int nthreads;  ///< 0 -> auto choose, >0 - user set.
	int err_freq;  ///< the frequency in which the stopping err on loss will be tested, reccomended > 10

	int normalize = 0;

	MedGDLMParams() {
		max_iter = 500; stop_at_err = (float)1e-4; max_times_err_grows = 20; method = "logistic_sgd"; batch_size = 512; rate = (float)0.01; rate_decay = (float)1.0; momentum = (float)0.95; last_is_bias = 0;
		l_ridge = (float)0; l_lasso = (float)0;  nthreads = 0; err_freq = 10; normalize = 0; verbose_learn = true;
	}

	ADD_CLASS_NAME(MedGDLMParams)
		ADD_SERIALIZATION_FUNCS(method, last_is_bias, max_iter, stop_at_err, max_times_err_grows, batch_size, rate, rate_decay, l_ridge, l_lasso, ls_lasso, ls_ridge, nthreads, err_freq)
};

class MedGDLM : public MedPredictor {
public:
	// Model
	int n_ftrs;
	vector<float> b;
	float b0;

	// Parameters
	MedGDLMParams params;

	// Function
	MedGDLM();
	~MedGDLM() {};
	MedGDLM(void *params);
	MedGDLM(MedGDLMParams& params);
	/// The parsed fields from init command.
	/// @snippet MedGDLM.cpp MedGDLM::init
	int set_params(map<string, string>& mapper);
	int init(void *params);
	void init_defaults();

	//int learn(MedMat<float> &x, MedMat<float> &y) {return (MedPredictor::learn(x,y));}; 	// Special case - un-normalized Y

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) const;

	int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std);

	void print(FILE *fp, const string& prefix, int level = 0) const;

	void calc_feature_contribs(MedMat<float> &x, MedMat<float> &contribs);

	ADD_CLASS_NAME(MedGDLM)
	ADD_SERIALIZATION_FUNCS(classifier_type, params, n_ftrs, b, b0, model_features, features_count)


	// actual computation functions
	int Learn_full(float *x, float *y, const float *w, int nsamples, int nftrs); // full non-iterative solution, not supporting lasso
	int Learn_gd(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Learn_sgd(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Learn_logistic_sgd(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Learn_logistic_sgd_threaded(float *x, float *y, const float *w, int nsamples, int nftrs);
private:
	void set_eigen_threads() const;
	void calc_feature_importance(vector<float> &features_importance_scores, const string &general_params, const MedFeatures *features);
};

MEDSERIALIZE_SUPPORT(MedGDLMParams)
MEDSERIALIZE_SUPPORT(MedGDLM)
