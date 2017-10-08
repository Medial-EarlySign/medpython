//
// MedAlgo - APIs to different algorithms: Linear Models, RF, GBM, KNN, and more
//

#ifndef __MED_ALGO_H__
#define __MED_ALGO_H__

#include "Logger/Logger/Logger.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "MedStat/MedStat/MedStat.h"
#include "MedFeat/MedFeat/MedFeat.h"
#include "QRF/QRF/QRF.h"
#include "gbm/gbm/gbm_utils.h"
#include "micNet/micNet/micNet.h"
#include "string.h"
#include "limits.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "svm.h"

// #include "MedBooster.h" // this include is at the end of file as it depends on following definitions to come first

// Forward Declaration
class MedFeatures;

#pragma warning(disable: 4297) //disable annoying " function assumed not to throw an exception but does "

// ToDo:: move defalts from classifiers in here (NITERS etc...)
//#include "classifiers/classifiers/classifiers.h"
#define NITER 200
#define EITER 0.00001
#define RFACTOR 0.99

#include <map>
#include <string>


using namespace std;


//Macros area

// Access to float* matrix:

#define XIDX(i,j,ncol) ((i)*(ncol) + (j))

//================================================================================
// MedPredictor - wrapper for classical learn/predict algorithms
//================================================================================


// Holder for 

// Model Types
typedef enum {
	MODEL_LINEAR_MODEL = 0, // Linear Model
	MODEL_QRF = 1, // Q-Random-Forest
	MODEL_GBM = 2, // Gradient Boosting Model
	MODEL_KNN = 3, // K Nearest Neighbour
	MODEL_BP = 4, // Neural Network Back Propagation
	MODEL_MARS = 5, // Multivariate Adaptive Regression Splines
	MODEL_GD_LINEAR = 6, // Gradient Descent/Full solution ridge
	MODEL_MULTI_CLASS = 7, // general one vs. all multi class extention
	MODEL_XGB = 8, // XGBoost
	MODEL_LASSO = 9, //Lasso model
	MODEL_MIC_NET = 10, //Home brew Neural Net implementation (Allows deep learning)
	MODEL_BOOSTER = 11, //general booster (meta algorithm)
	MODEL_DEEP_BIT = 12, //general booster (meta algorithm)
	MODEL_LIGHTGBM = 13, // the celebrated LightGBM algorithm
	MODEL_SPECIFIC_GROUPS_MODELS = 14,
	MODEL_SVM = 15,
	MODEL_LINEAR_SGD = 16,
	MODEL_LAST
} MedPredictorTypes;

MedPredictorTypes predictor_name_to_type(const string& model_name);

class MedPredictor : public SerializableObject {
public:
	MedPredictorTypes classifier_type;

	// General constructor
	MedPredictor() {}
	virtual ~MedPredictor() {};

	bool transpose_for_learn;
	bool normalize_for_learn;
	bool normalize_y_for_learn;

	bool transpose_for_predict;
	bool normalize_for_predict;

	// Each wrapped algorithm needs to implement the following:
	//.........................................................
	// Init
	virtual int init(void *classifier_params) { return 0; };
	int init_from_string(string initialization_text);
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};


	// Learn
	// should be implemented for each model. This API always assumes the data is already normalized/transposed as needed, 
	// and never changes data in x,y,w. method should support calling with w=NULL.
	virtual int Learn(float *x, float *y, float *w, int n_samples, int n_ftrs) { return 0; };

	// Predict
	// should be implemented for each model. This API assumes x is normalized/transposed if needed.
	// preds should either be pre-allocated or NULL - in which case the predictor should allocate it to the right size.
	virtual int Predict(float *x, float *&preds, int n_samples, int n_ftrs) { return 0; }

	virtual size_t get_size() { return 0; }
	virtual size_t serialize(unsigned char *blob) { return 0; }
	virtual size_t deserialize(unsigned char *blob) { return 0; }

	// Print
	virtual void print(FILE *fp, const string& prefix);

	// Number of predictions per sample. typically 1 - but some models return several per sample (for example a probability vector)
	virtual int n_preds_per_sample() { return 1; };

	virtual int denormalize_model(float *f_avg, float *f_std, float label_avg, float label_std) { return 0; };

	// methods relying on virtual methods, and applicable to all predictors: (one can still reimplement in derived class if needed)
	//..............................................................................................................................

	// simple no weights call
	int learn(float *x, float *y, int nsamples, int nftrs) { return Learn(x, y, NULL, nsamples, nftrs); }

	// simple c++ style learn

	// MedMat x,y : will transpose/normalize x,y if needed by algorithm
	// The convention is that untransposed mats are always samples x features, and transposed are features x samples
	int learn(MedMat<float> &x, MedMat<float> &y, vector<float> &wgts);
	int learn(MedMat<float> &x, MedMat<float> &y) { vector<float> w; return(learn(x, y, w)); }

	// MedMat x, vector y: will transpose normalize x if needed (y assumed to be normalized)
	int learn(MedMat<float> &x, vector<float> &y, vector<float> &wgts);
	int learn(MedMat<float> &x, vector<float> &y) { vector<float> w; return(learn(x, y, w)); }

	// vector x,y: transpose/normalizations not done.
	int learn(vector<float> &x, vector<float> &y, vector<float> &wgts, int n_samples, int n_ftrs);
	int learn(vector<float> &x, vector<float> &y, int n_samples, int n_ftrs) { vector<float> w; return learn(x, y, w, n_samples, n_ftrs); }

	// simple c++ style predict
	int predict(MedMat<float> &x, vector<float> &preds);
	int predict(vector<float> &x, vector<float> &preds, int n_samples, int n_ftrs);
	int threaded_predict(MedMat<float> &x, vector<float> &preds, int nthreads);

	// MedFeaturesData related
	int learn(MedFeaturesData& data, int isplit);

	int learn(MedFeaturesData& data);
	int learn(MedFeatures& features);
	int learn(MedFeatures& features, vector<string>& names);
	int cross_validate_splits(MedFeaturesData& data);
	int predict(MedFeaturesData& data, int isplit);
	int predict_on_train(MedFeaturesData& data, int isplit);
	int predict(MedFeaturesData& data);
	int predict(MedFeatures& features);

	//caliberation for probability using training:
	int learn_prob_calibration(MedMat<float> &x, vector<float> &y,
		vector<float> &min_range, vector<float> &max_range, vector<float> &map_prob, int min_bucket_size = 10000);
	int convert_scores_to_prob(const vector<float> &preds, const vector<float> &min_range,
		const vector<float> &max_range, const vector<float> &map_prob, vector<float> &probs);
	int learn_prob_calibration(MedMat<float> &x, vector<float> &y, float &A, float &B, int min_bucket_size = 10000);
	int convert_scores_to_prob(const vector<float> &preds, float A, float B, vector<float> &probs);

	// init
	static MedPredictor *make_predictor(string model_type);
	static MedPredictor *make_predictor(MedPredictorTypes model_type);
	static MedPredictor *make_predictor(string model_type, string params);
	static MedPredictor *make_predictor(MedPredictorTypes model_type, string params);

	// (De)Serialize
	size_t get_predictor_size();
	size_t predictor_serialize(unsigned char *blob);

	// write/read from file
	int read_from_file(const string &fname); // read and deserialize model
	int write_to_file(const string &fname);  // serialize model and write to file

private:
	// some needed helpers
	void prepare_x_mat(MedMat<float> &x, vector<float> &wgts, int &nsamples, int &nftrs, bool transpose_needed);
	void predict_thread(void *p);
	void build_learning_x_mat_for_split(MedFeaturesData & ftrs_data, vector<float>& signal, int isplit, MedMat<float>& x);

};

//======================================================================================
// Linear Model: Linear Regression (+ Ridge)
//======================================================================================

#define LM_NITER 200
#define LM_EITER 0.00001

struct MedLMParams {

	// Required params
	float eiter;
	int niter;


	// A simple way to check a single column , default is -1, but if >=0 the algorithm will simply return this column as prediction
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

	// Parameters
	MedLMParams params;

	// Function
	MedLM();
	MedLM(void *params);
	MedLM(MedLMParams& params);
	int init(void *params);
	virtual int init(map<string, string>& mapper);
	void init_defaults();

	//int learn(MedMat<float> &x, MedMat<float> &y) {return (MedPredictor::learn(x,y));}; 	// Special case - un-normalized Y

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag);

	int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std);

	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void print(FILE *fp, const string& prefix);
};

// Ancillary function for string analysis
int init_farray(string& in, float **out);
int init_darray(string& in, int **out);

//======================================================================================
// Linear Model: Linear Regression (+ Lasso)
//======================================================================================

#define LASSO_LAMBDA 0;
#define LASSO_NITER 1000;

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

	// Parameters
	MedLassoParams params;

	// Work variables
	double **trainx;
	double *y1;

	// Function
	void Initb();
	MedLasso();
	MedLasso(void *params);
	MedLasso(MedLassoParams& params);
	int init(void *params);
	virtual int init(map<string, string>& mapper);
	void init_defaults();

	//int learn(MedMat<float> &x, MedMat<float> &y) {return (MedPredictor::learn(x,y));}; 	// Special case - un-normalized Y

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag);

	int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std);

	void initialize_vars(float *x_in, float *y_in, float *w, vector<float>& b, int nrow_train, int n_ftrs);
	void lasso_regression(vector<float>& b, int nrow_train, int n_ftrs, double lambda, int num_iterations);

	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void print(FILE *fp, const string& prefix);
};
// Least Square direct iterations solution
int learn_lm(float *x, float *_y, float *w, int nsamples, int nftrs, int niter, float eiter, float *rfactors, float *b, float *err, float *corrs);
int learn_lm(float *x, float *_y, float *w, int nsamples, int nftrs, int niter, float eiter, float *rfactors, float *b, float *err, float *corrs, float *sumxx);

//==============================================================================================
// Linear Models2: Linear regression (with Ridge and/or Lasso), using Gradient Descent variants
//==============================================================================================
struct MedGDLMParams : public SerializableObject {

	// Required params
	int max_iter;
	float stop_at_err; // stop criteria
	int max_times_err_grows;
	string method; // gd or sgd
	int batch_size;	// for sgd
	float rate;
	float rate_decay;
	float momentum;

	int last_is_bias;
	int print_model;

	// Optional params
	float l_ridge; // lambda for ridge
	float l_lasso; // labmda for lasso

	int nthreads;  // 0 -> auto choose, >0 - user set.
	int err_freq;  // the frequency in which the stopping err on loss will be tested, reccomended > 10

	MedGDLMParams() {
		max_iter = 500; stop_at_err = (float)1e-4; max_times_err_grows = 20; method = "logistic_sgd"; batch_size = 512; rate = (float)0.01; rate_decay = (float)1.0; momentum = (float)0.95; last_is_bias = 0;
		l_ridge = (float)0; l_lasso = (float)0; nthreads = 0; err_freq = 10;
	}

	ADD_SERIALIZATION_FUNCS(method, last_is_bias, max_iter, stop_at_err, max_times_err_grows, batch_size, rate, rate_decay, l_ridge, l_lasso, nthreads, err_freq);
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
	MedGDLM(void *params);
	MedGDLM(MedGDLMParams& params);
	int init(map<string, string>& mapper); ;
	int init(void *params);
	void init_defaults();

	//int learn(MedMat<float> &x, MedMat<float> &y) {return (MedPredictor::learn(x,y));}; 	// Special case - un-normalized Y

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag);

	ADD_SERIALIZATION_FUNCS(params, n_ftrs, b, b0);

	int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std);

	void print(FILE *fp, const string& prefix);

	// actual computation functions
	int Learn_full(float *x, float *y, float *w, int nsamples, int nftrs); // full non-iterative solution, not supporting lasso
	int Learn_gd(float *x, float *y, float *w, int nsamples, int nftrs);
	int Learn_sgd(float *x, float *y, float *w, int nsamples, int nftrs);
	int Learn_logistic_sgd(float *x, float *y, float *w, int nsamples, int nftrs);
	int Learn_logistic_sgd_threaded(float *x, float *y, float *w, int nsamples, int nftrs);
private:
	void set_eigen_threads();
};

void init_default_lm_params(MedLMParams& _parmas);

//======================================================================================
// QRF: Quantized Regression/Classification random forest
//======================================================================================
#define MED_QRF_DEF_NTREES 100
#define MED_QRF_DEF_MAXQ 200
#define MED_QRF_DEF_MIN_NODE 50
#define MED_QRF_DEF_LEARN_NTHREADS 8
#define MED_QRF_DEF_PREDICT_NTHREADS 8
#define MED_QRF_DEF_SPREAD	0.1

struct MedQRFParams {

	// Required
	int ntrees;
	int maxq;
	int learn_nthreads, predict_nthreads;
	QRF_TreeType type;

	// Optional
	int max_samp; // if > 0 & sampsize is NULL : the maximal sampsize we will take from each category
	float samp_factor; // if > 0 & sampsize if NULL : the maximal factor of samples between the 2 largest categories
	vector<int> samp_vec; // to be used when sampsize is NULL and max_samp,samp_vector > 0
	int *sampsize;
	int ntry;
	int get_only_this_categ;
	int max_depth; //maximial depth of tree branches - if 0 no limit

	// Regression
	float spread;
	bool keep_all_values; // For quantile regression

	// categorical
	int min_node;
	int n_categ;

	int collect_oob;

	// For Prediction
	int get_count;
	vector<float> quantiles; // For quantile regression
};

class MedQRF : public MedPredictor {
public:
	// Model 
	QRF_Forest qf;

	// Parameters
	MedQRFParams params;

	// Function
	MedQRF();
	MedQRF(void *params);
	MedQRF(MedQRFParams& params);
	int init(void *params);
	virtual int init(map<string, string>& mapper);
	//	int init(const string &init_str); // allows init of parameters from a string. Format is: param=val,... , for sampsize: 0 is NULL, a list of values is separated by ; (and not ,)
	void init_defaults();
	QRF_TreeType get_tree_type(string name);

	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int get_count);

	//int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std) {return 0;};

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Print
	void print(FILE *fp, const string& prefix);

	void printTrees(const vector<string> &modelSignalNames, const string &outputPath);
	// Predictions per sample
	int n_preds_per_sample();

private:
	void set_sampsize(float *y, int nsamples); // checking if there's a need to prep sampsize based on max_samp and samp_factor
};

//======================================================================================
// micNet: Home brewed implimenatation for Neural Nets and Deep Learning
//======================================================================================
struct MedMicNetParams {

	string init_string;
};

class MedMicNet : public MedPredictor {
public:
	// Model 
	micNet mic;

	// Parameters
	MedMicNetParams mic_params;

	// Function
	MedMicNet() { classifier_type = MODEL_MIC_NET; mic.params.init_defaults(); }
	MedMicNet(void *params) { mic_params = *(MedMicNetParams *)params; mic.init_from_string(mic_params.init_string); }
	MedMicNet(MedMicNetParams& params) { mic_params = params; mic.init_from_string(mic_params.init_string); }
	int init(void *params) { mic_params = *(MedMicNetParams *)params; return mic.init_from_string(mic_params.init_string); }
	int init_from_string(string initialization_text) {
		cerr << "MedMicNet init_from_string ! :: " << initialization_text << "\n";
		mic_params.init_string = initialization_text;
		cerr << "calling init_from_string of micNet\n"; fflush(stderr);
		return mic.init_net(initialization_text);
	}

	int init(map<string, string>& mapper) {
		cerr << "MedMicNet:: init map :: not supported, only init_from_string supported....\n";
		return -1;
	}
	//	int init(const string &init_str); // allows init of parameters from a string. Format is: param=val,... , for sampsize: 0 is NULL, a list of values is separated by ; (and not ,)
	void init_defaults() { mic_params.init_string = ""; mic.params.init_defaults(); }

	int Learn(float *x, float *y, float *w, int nsamples, int nftrs) {
		cerr << "MedMicNet:: Learn :: API's with MedMat are preferred....\n";
		MedMat<float> xmat; xmat.load(x, nsamples, nftrs);
		MedMat<float> ymat; ymat.load(y, nsamples, 1);
		return learn(xmat, ymat);
	}

	int Predict(float *x, float *&preds, int nsamples, int nftrs) {
		cerr << "MedMicNet:: Learn :: API's with MedMat are preferred....\n";
		MedMat<float> xmat; xmat.load(x, nsamples, nftrs);
		vector<float> vpreds;
		int rc = predict(xmat, vpreds);
		if (preds == NULL) preds = new float[nsamples];
		memcpy(preds, &vpreds[0], sizeof(float)*nsamples);
		return rc;
	}

	int learn(MedMat<float> &x, MedMat<float> &y) { return mic.learn(x, y); }
	int predict(MedMat<float> &x, vector<float> &preds) { return mic.predict(x, preds); }

	// (De)Serialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	size_t get_size() { return mic.get_size(); }
	size_t serialize(unsigned char *blob) { return mic.serialize(blob); }
	size_t deserialize(unsigned char *blob) { return mic.deserialize(blob); }

	// Print
	//void print(FILE *fp, const string& prefix);

	// Predictions per sample
	int n_preds_per_sample() { return mic.n_preds_per_sample(); }

};

//======================================================================================
// Mars: C++ version of Mars
//======================================================================================
#define MED_MARS_DEF_MAXTERMS			50
#define MED_MARS_DEF_DEGREE				2
#define MED_MARS_DEF_PENALTY			3
#define MED_MARS_DEF_THRESH				0.001
#define MED_MARS_DEF_MINSPAN			0
#define MED_MARS_DEF_PRUNE				1
#define MED_MARS_DEF_FASTK				20
#define MED_MARS_DEF_FASTBETA			0
#define MED_MARS_DEF_NEWVAR_PENALTY		0
#define MED_MARS_DEF_USEBETACACHE		1
#define MED_MARS_DEF_TRACE				3

struct MedMarsParams {

	// Required
	int MaxTerms;	// Maximal number of Terms in final model
	int MaxDegree;	// Model max linear degree : 1 - linear, 2 - square, 3 - cubic, etc...
	double Penalty;
	double Thresh;
	int MinSpan;
	bool Prune;
	int FastK;
	double FastBeta;
	double NewVarPenalty;
	bool UseBetaCache;
	double Trace;	// debug prints during algorithm run (recommended): 0: no prints , 3: print all (recommended)
};

class MedMars : public MedPredictor {
public:
	// Model 
	int nMaxTerms;
	int nTerms;
	int nPreds;
	bool *BestSet;			// size: nMaxTerms
	vector<int> Dirs;		// size: nMaxTerms*nPreds
	vector<double> Cuts;	// size: nMaxTerms*nPreds
	vector<double> Betas;	// size: nMaxTerms

	// Model Inner quality measures
	double BestGcv;
	vector<double> bx;
	vector<double> Residuals;

	// Parameters
	MedMarsParams params;

	// Function
	MedMars();
	MedMars(void *params);
	MedMars(MedMarsParams& params);
	int init(void *params);
	virtual int init(map<string, string>& mapper);
	void init_defaults();

	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	//int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std) {return 0;};

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Print
	void print(FILE *fp, const string& prefix);

	// Predictions per sample
	int n_preds_per_sample();
};

// Initialization of parameters
void init_default_mars_params(MedMarsParams& _params);

//======================================================================================
// GBM: C++ version of GBM from R.
//======================================================================================
#define MED_GBM_DEF_SHRINKAGE 0.1
#define MED_GBM_DEF_BAG_P 0.5
#define MED_GBM_DEF_NTREES 100
#define MED_GBM_DEF_DEPTH 6
#define MED_GBM_DEF_TAKE_ALL_POS false
#define MED_GBM_DEF_MIN_OBS_IN_NODE 10

typedef gbm_parameters MedGBMParams;

class MedGBM : public MedPredictor {
public:
	// Loss function
	GBM_LossFunctions loss_function;

	// Alpha for quantile loss function
	double alpha_quantile;

	// Model 
	full_gbm_learn_info_t gbm_model;

	// Parameters
	MedGBMParams params;

	// Predicting on subset of trees
	int predict_ntrees;

	// Function
	MedGBM();
	MedGBM(void *params);
	MedGBM(MedGBMParams& params);
	int init(void *params);
	virtual int init(map<string, string>& mapper); ;
	void init_defaults();
	GBM_LossFunctions get_loss_function(string name);
	~MedGBM();

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Print
	void print(FILE *fp, const string& prefix);
};

// Initialization of parameters
void init_default_gbm_params(MedGBMParams& _params);

//======================================================================================
// KNN
//======================================================================================
typedef enum {
	KNN_DIST_MEAN,
	KNN_1_DIST,
	KNN_WEIGHTEDLS,
	KNN_AVG_LAST
} knnAveraging;

typedef enum {
	KNN_L1,
	KNN_L2,
	KNN_METRIC_LAST
}knnMetric;

struct MedKNNParams {

	int k;
	knnAveraging knnAv;
	knnMetric knnMetr;
};

class MedKNN : public MedPredictor {
public:
	// Model
	int nsamples;
	int nftrs;
	float *x;
	float *y;
	float *w;


	// Parameters
	MedKNNParams params;


	// Function
	MedKNN();
	MedKNN(void *params);
	MedKNN(MedKNNParams& params);
	virtual int init(map<string, string>& mapper); ;
	int init(void *params);
	~MedKNN();
	knnAveraging get_knn_averaging(string name);
	knnMetric get_knn_metric(string name);

	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	//	void print(FILE *fp, const string& prefix) ;
};


//======================================================================================
// BackProp 
//======================================================================================
typedef enum { SIGMOID, RELU, LINEAR }neuronFunT;

typedef struct {
	int layerIndex;
	neuronFunT neuronFunction;
	int x;
	int y;
	int  firstWeight, lastWeight;
	int firstSource, lastSource;
	double value;
	double error;
	double delta;
}neuronStruct;

//================================================================
typedef struct {
	neuronStruct *neuron;


	int *source;
	double *weight;
	int numLayers;
	int numNeurons;
	int numInputs, numOutputs;
	int numWeights, numSource;

}netStruct;

struct MedBPParams {

	int numLayers;

	int numIterations;
	double alpha; // learning rate
	double beta;// parameter of logistic function
};

class MedBP : public MedPredictor {
public:
	// Model

	int nsamples;
	int nftrs;
	/*double **x;
	double **y;
	float *w;
	*/

	// Function
	MedBP();
	MedBP(void *params);
	MedBP(MedBPParams& params);
	int init(void *params);
	virtual int init(map<string, string>& mapper);
	~MedBP();

	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	//	void print(FILE *fp, const string& prefix) ;
		// Parameters
private:
	MedBPParams params;

	netStruct network;
};

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
	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Print
	void print(FILE *fp, const string& prefix);

	// Predictions per sample
	int n_preds_per_sample();
};


//================================================================
// Unsupervised
//================================================================

// K-Means
int KMeans(MedMat<float> &x, int K, MedMat<float> &centers, vector<int> &clusters, MedMat<float> &dists);
int KMeans(MedMat<float> &x, int K, int max_iter, MedMat<float> &centers, vector<int> &clusters, MedMat<float> &dists);
int KMeans(float *x, int nrows, int ncols, int K, float *centers, int *clusters, float *dists);

int KMeans(float *x, int nrows, int ncols, int K, int max_iter, float *centers, int *clusters, float *dists); // actual implemetation routine

// PCA

// given a matrix, returns the base PCA matrix and the cummulative relative variance explained by them.
// it is highly recommended to normalize the input matrix x before calling.
int MedPCA(MedMat<float> &x, MedMat<float> &pca_base, vector<float> &varsum);

// returns the projection of the pca base on the first dim dimensions.
int MedPCA_project(MedMat<float> &x, MedMat<float> &pca_base, int dim, MedMat<float> &projected);

//wrapper for MedPredictor for certian groups - routes the input to correct model group.
//for example may be used to train specific model for each age group
class MedSpecificGroupModels : public MedPredictor {
public:
	// Model

	int nsamples;
	int nftrs;
	/*double **x;
	double **y;
	float *w;
	*/

	// Function
	MedSpecificGroupModels();
	~MedSpecificGroupModels();

	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);
	MedSpecificGroupModels *clone();
	ADD_SERIALIZATION_FUNCS(featNum, feat_ths, predictors)

		//	void print(FILE *fp, const string& prefix) ;
		// Parameters
		void set_predictors(const vector<MedPredictor *> &predictors); //for each group index
	void set_group_selection(int featNum, const vector<float> &feat_ths);
	MedPredictor *get_model(int ind);
	int model_cnt();
private:
	vector<MedPredictor *> predictors;
	int featNum;
	vector<float> feat_ths;
	int selectPredictor(const float *x); //retrieve predictor index
};

class MedSvm : public MedPredictor {
public:
	// Model

	struct svm_parameter params;
	struct svm_model *model;
	/*double **x;
	double **y;
	float *w;
	*/

	// Function
	MedSvm();
	MedSvm(void *params);
	MedSvm(struct svm_parameter &params);
	~MedSvm();

	void init_defaults();
	int init(void *params);
	virtual int init(map<string, string>& mapper);
	int init(struct svm_parameter &params);

	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

private:


};

//================================================================
// dependent includes
//================================================================

#include "MedBooster.h"
//#include "MedLightGBM.h"


//=================================================================
// Joining the MedSerialize Wagon
//=================================================================

MEDSERIALIZE_SUPPORT(MedLM)
MEDSERIALIZE_SUPPORT(MedLasso)
MEDSERIALIZE_SUPPORT(MedGDLMParams)
MEDSERIALIZE_SUPPORT(MedGDLM)
MEDSERIALIZE_SUPPORT(MedQRF)
MEDSERIALIZE_SUPPORT(MedMicNet)
MEDSERIALIZE_SUPPORT(MedBP)
MEDSERIALIZE_SUPPORT(MedMars)
MEDSERIALIZE_SUPPORT(MedKNN)
MEDSERIALIZE_SUPPORT(MedGBM)
MEDSERIALIZE_SUPPORT(MedMultiClass)
MEDSERIALIZE_SUPPORT(MedSpecificGroupModels)


#endif