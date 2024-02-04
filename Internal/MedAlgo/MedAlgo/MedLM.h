#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>

#define NITER 200
#define EITER 0.00001
#define RFACTOR 0.99

#define XIDX(i,j,ncol) ((i)*(ncol) + (j))

//======================================================================================
// Linear Model: Linear Regression (+ Ridge)
//======================================================================================

#define LM_NITER 200
#define LM_EITER 0.00001

struct MedLMParams {

	// Required params
	float eiter;
	int niter;


	/// A simple way to check a single column , default is -1, but if >=0 the algorithm will simply return this column as prediction
	int get_col = -1;

	// Optional params
	float rfactor;
	float *rfactors;
	float *corrs;
	float *sumxx;

	MedLMParams() { eiter = (float)EITER; niter = NITER; rfactor = 1.0; rfactors = NULL; corrs = NULL; sumxx = NULL; get_col = -1; }

};

class MedLM : public MedPredictor {
public:
	// Model
	int n_ftrs;
	vector<float> b;
	float b0;
	float err;

	/// Parameters
	MedLMParams params;

	// Function
	MedLM();
	~MedLM() {};
	MedLM(void *params);
	MedLM(MedLMParams& params);
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedLM.cpp MedLM::init
	virtual int set_params(map<string, string>& mapper);
	void init_defaults();

	//int learn(MedMat<float> &x, MedMat<float> &y) {return (MedPredictor::learn(x,y));}; 	// Special case - un-normalized Y
	void calc_feature_contribs(MedMat<float> &x, MedMat<float> &contribs);

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) const;

	void normalize_x_and_y(float *x, float *y, const float *w, int nsamples, int nftrs, vector<float>& x_avg, vector<float>& x_std, float& y_avg, float& y_std);
	int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std);
	void print(FILE *fp, const string& prefix, int level = 0) const;

	bool predict_single_not_implemented() { return true; }

	ADD_CLASS_NAME(MedLM)
		ADD_SERIALIZATION_FUNCS(classifier_type, n_ftrs, b0, b, err)
};

// Ancillary function for string analysis
int init_farray(string& in, float **out);
int init_darray(string& in, int **out);

//======================================================================================
// Linear Model: Linear Regression (+ Lasso)
//======================================================================================

#define LASSO_LAMBDA 0
#define LASSO_NITER 1000

struct MedLassoParams {

	// Required params
	double lambda;
	int num_iterations;
	MedLassoParams() { lambda = LASSO_LAMBDA; num_iterations = LASSO_NITER; }

};

class MedLasso : public MedPredictor {

public:
	// Model
	int n_ftrs;
	vector<float> b;
	float b0;

	/// Parameters
	MedLassoParams params;

	// Work variables
	double **trainx;
	double *y1;

	// Function
	void Initb();
	MedLasso();
	~MedLasso() {};
	MedLasso(void *params);
	MedLasso(MedLassoParams& params);
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedLasso.cpp MedLasso::init
	int set_params(map<string, string>& mapper);
	void init_defaults();

	//int learn(MedMat<float> &x, MedMat<float> &y) {return (MedPredictor::learn(x,y));}; 	// Special case - un-normalized Y

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) const;

	void normalize_x_and_y(float *x, float *y, const float *w, int nsamples, int nftrs, vector<float>& x_avg, vector<float>& x_std, float& y_avg, float& y_std);
	int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std);

	void initialize_vars(float *x_in, float *y_in, const float *w, vector<float>& b, int nrow_train, int n_ftrs);
	void lasso_regression(vector<float>& b, int nrow_train, int n_ftrs, double lambda, int num_iterations);
	void print(FILE *fp, const string& prefix, int level = 0) const;

	ADD_CLASS_NAME(MedLasso)
		ADD_SERIALIZATION_FUNCS(classifier_type, n_ftrs, b0, b)
};
/// Least Square direct iterations solution
int learn_lm(float *x, float *_y, const float *w, int nsamples, int nftrs, int niter, float eiter, float *rfactors, float *b, float *err, float *corrs);
/// Least Square direct iterations solution
int learn_lm(float *x, float *_y, const float *w, int nsamples, int nftrs, int niter, float eiter, float *rfactors, float *b, float *err, float *corrs, float *sumxx);

void init_default_lm_params(MedLMParams& _parmas);

MEDSERIALIZE_SUPPORT(MedLM)
MEDSERIALIZE_SUPPORT(MedLasso)
