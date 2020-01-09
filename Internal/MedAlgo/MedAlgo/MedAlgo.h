/// @file
/// MedAlgo - APIs to different algorithms: Linear Models, RF, GBM, KNN, and more
///

#ifndef __MED_ALGO_H__
#define __MED_ALGO_H__

#if __GNUC__  >= 5 || (defined(_MSC_VER) && !defined(_DEBUG))
#define NEW_COMPLIER true
#else
#define NEW_COMPLIER false
#endif

#include <Logger/Logger/Logger.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <MedStat/MedStat/MedStat.h>
#include <MedFeat/MedFeat/MedFeat.h>
#include <QRF/QRF/QRF.h>
#include <gbm/gbm/gbm_utils.h>
#include <micNet/micNet/micNet.h>
#include <string.h>
#include <limits.h>
#include <MedProcessTools/MedProcessTools/MedProcessUtils.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <TQRF/TQRF/TQRF.h>
#include "svm.h"
#include <unordered_map>
#include <random>

// #include "MedBooster.h" // this include is at the end of file as it depends on following definitions to come first

// Forward Declaration
class MedFeatures;

#pragma warning(disable: 4297) //disable annoying " function assumed not to throw an exception but does "

/// @todo move defalts from classifiers in here (NITERS etc...)
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

/// @enum
/// Model Types options
typedef enum {
	MODEL_LINEAR_MODEL = 0, ///< to_use:"linear_model" Linear %Model - creates MedLM
	MODEL_QRF = 1, ///< to_use:"qrf" Q-Random-Forest - creates MedQRF
	MODEL_GBM = 2, ///< to_use:"gbm" Gradient Boosting %Model - creates MedGBM
	MODEL_KNN = 3, ///< to_use:"knn" K Nearest Neighbour - creates MedKNN
	MODEL_BP = 4, ///< to_use:"BP" Neural Network Back Propagation - creates MedBP
	MODEL_MARS = 5, ///< to_use:"mars" Multivariate Adaptive Regression Splines - creates MedMars
	MODEL_GD_LINEAR = 6, ///< to_use:"gdlm" Gradient Descent/Full solution ridge - creates MedGDLM
	MODEL_MULTI_CLASS = 7, ///< to_use:"multi_class" general one vs. all multi class extention - creates MedMultiClass
	MODEL_XGB = 8, ///< to_use:"xgb" XGBoost - creates MedXGB
	MODEL_LASSO = 9, ///< to_use:"lasso" Lasso model - creates MedLasso
	MODEL_MIC_NET = 10, ///< to_use:"micNet" Home brew Neural Net implementation (Allows deep learning) - creates MedMicNet
	MODEL_BOOSTER = 11, ///< to_use:"booster" general booster (meta algorithm) - creates MedBooster
	MODEL_DEEP_BIT = 12, ///< to_use:"deep_bit" Nir\'s DeepBit method - creates MedDeepBit
	MODEL_LIGHTGBM = 13, ///< to_use:"lightgbm" the celebrated LightGBM algorithm - creates MedLightGBM
	MODEL_SPECIFIC_GROUPS_MODELS = 14, ///< to_use:"multi_models" spliting model by specific value (for example age-range) and train diffretn model for each bin - creates MedSpecificGroupModels
	MODEL_SVM = 15, ///< to_use:"svm" Svm model - creates MedSvm 
	MODEL_LINEAR_SGD = 16, ///< to_use:"linear_sgd" linear model using our customized SGD - creates MedLinearModel
	MODEL_VW = 17, ///< to_use:"vw" %VowpalWabbit yahoo reasearch library - creates MedVW
	MODEL_TQRF = 18, ///< to_use:"tqrf" TQRF model
	MODEL_BART = 19, ///< to_use:"bart" MedBART model using BART
	MODEL_LAST
} MedPredictorTypes;

///Maping from predictor enum type ::MedPredictorTypes to model name in string
extern unordered_map<int, string> predictor_type_to_name;
///Maping from model name in string to enum ::MedPredictorTypes 
MedPredictorTypes predictor_name_to_type(const string& model_name);

/**
* Base Interface for predictor
*/
class MedPredictor : public SerializableObject {
public:
	MedPredictorTypes classifier_type; ///<The Predicotr enum type

	// General constructor
	MedPredictor() {}
	virtual ~MedPredictor() {};

	bool transpose_for_learn; ///<True if need to transpose before learn
	bool normalize_for_learn; ///<True if need to normalize before learn
	bool normalize_y_for_learn; ///<True if need to normalize labels before learn

	bool transpose_for_predict; ///<True if need to transpose before predict
	bool normalize_for_predict; ///<True if need to normalize before predict

	vector<string> model_features; ///<The model features used in Learn, to validate when caling predict
	///The model features count used in Learn, to validate when caling predict. 
	///used if model_features is empty because feature names aren't availabe during learn
	int features_count = 0;

	// Each wrapped algorithm needs to implement the following:
	//.........................................................
	// Init
	virtual int init(void *classifier_params) { return 0; };
	int init_from_string(string initialization_text);
	virtual int init(map<string, string>& mapper);
	virtual int set_params(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};


	/// Learn
	/// should be implemented for each model. This API always assumes the data is already normalized/transposed as needed, 
	/// and never changes data in x,y,w. method should support calling with w=NULL.
	virtual int Learn(float *x, float *y, const float *w, int n_samples, int n_ftrs) { return 0; };

	/// Predict
	/// should be implemented for each model. This API assumes x is normalized/transposed if needed.
	/// preds should either be pre-allocated or NULL - in which case the predictor should allocate it to the right size.
	virtual int Predict(float *x, float *&preds, int n_samples, int n_ftrs) const { return 0; }

	// Print
	virtual void print(FILE *fp, const string& prefix, int level = 0) const;

	/// Number of predictions per sample. typically 1 - but some models return several per sample (for example a probability vector)
	virtual int n_preds_per_sample() const { return 1; };

	virtual int denormalize_model(float *f_avg, float *f_std, float label_avg, float label_std) { return 0; };

	// methods relying on virtual methods, and applicable to all predictors: (one can still reimplement in derived class if needed)
	//..............................................................................................................................

	/// simple no weights call
	int learn(float *x, float *y, int nsamples, int nftrs) { return Learn(x, y, NULL, nsamples, nftrs); }

	// simple c++ style learn

	/// MedMat x,y : will transpose/normalize x,y if needed by algorithm
	/// The convention is that untransposed mats are always samples x features, and transposed are features x samples
	int learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts);
	/// MedMat x,y : will transpose/normalize x,y if needed by algorithm
	/// The convention is that untransposed mats are always samples x features, and transposed are features x samples
	int learn(MedMat<float> &x, MedMat<float> &y) { vector<float> w; return(learn(x, y, w)); }

	/// MedMat x, vector y: will transpose normalize x if needed (y assumed to be normalized)
	int learn(MedMat<float> &x, vector<float> &y, const vector<float> &wgts);
	/// MedMat x, vector y: will transpose normalize x if needed (y assumed to be normalized)
	int learn(MedMat<float> &x, vector<float> &y) { vector<float> w; return(learn(x, y, w)); }

	/// vector x,y: transpose/normalizations not done.
	int learn(vector<float> &x, vector<float> &y, const vector<float> &wgts, int n_samples, int n_ftrs);
	/// vector x,y: transpose/normalizations not done.
	int learn(vector<float> &x, vector<float> &y, int n_samples, int n_ftrs) { vector<float> w; return learn(x, y, w, n_samples, n_ftrs); }

	// simple c++ style predict
	int predict(MedMat<float> &x, vector<float> &preds) const;
	int predict(vector<float> &x, vector<float> &preds, int n_samples, int n_ftrs) const;
	int threaded_predict(MedMat<float> &x, vector<float> &preds, int nthreads) const;

	int learn(const MedFeatures& features);
	int learn(const MedFeatures& features, vector<string>& names);
	int predict(MedFeatures& features) const;

	///Feature Importance - assume called after learn
	virtual void calc_feature_importance(vector<float> &features_importance_scores,
		const string &general_params)
	{
		const MedFeatures *features = NULL;
		calc_feature_importance(features_importance_scores,
			general_params, features);
	}
	virtual void calc_feature_importance(vector<float> &features_importance_scores,
		const string &general_params, const MedFeatures *features)  {
		string model_name = "model_id=" + to_string(classifier_type);
		if (predictor_type_to_name.find(classifier_type) != predictor_type_to_name.end())
			model_name = predictor_type_to_name[classifier_type];
		throw logic_error("ERROR:: operation calc_feature_importance "
			"isn't supported for " + model_name + " yet.");
	};

	///Feature contributions explains the prediction on each sample (aka BUT_WHY)
	virtual void calc_feature_contribs(MedMat<float> &x, MedMat<float> &contribs) {
		string model_name = "model_id=" + to_string(classifier_type);
		if (predictor_type_to_name.find(classifier_type) != predictor_type_to_name.end())
			model_name = predictor_type_to_name[classifier_type];
		throw logic_error("ERROR:: operation calc_feature_contribs "
			"isn't supported for " + model_name + " yet.");
	};

	virtual void export_predictor(const string &output_fname) {
		string model_name = "model_id=" + to_string(classifier_type);
		if (predictor_type_to_name.find(classifier_type) != predictor_type_to_name.end())
			model_name = predictor_type_to_name[classifier_type];
		throw logic_error("ERROR:: operation export_predictor "
			"isn't supported for " + model_name + " yet.");
	}

	/// <summary>
	/// calibration for probability using training data
	/// @param x The training matrix
	/// @param y The Labels
	/// @param min_bucket_size The minimal observations to create probability bin
	/// @param min_score_jump The minimal diff in scores to create bin
	/// @param min_prob_jump The minimal diff in probabilties to create bin
	/// @param fix_prob_order If true will unite bins that are sorted in wrong way
	/// </summary>
	/// <returns>
	/// @param min_range - writes a corresponding vector with minimal score range
	/// @param max_range - writes a corresponding vector with maximal score range
	/// @param map_prob - writes a corresponding vector with probability for score range
	/// </returns>
	int learn_prob_calibration(MedMat<float> &x, vector<float> &y,
		vector<float> &min_range, vector<float> &max_range, vector<float> &map_prob, int min_bucket_size = 10000,
		float min_score_jump = 0.001, float min_prob_jump = 0.005, bool fix_prob_order = false);
	/// <summary>
	/// If you have ran learn_prob_calibration before, you have min_range,max_range,map_prob from
	/// This function - that is used to convert preds to probs
	/// </summary>
	int convert_scores_to_prob(const vector<float> &preds, const vector<float> &min_range,
		const vector<float> &max_range, const vector<float> &map_prob, vector<float> &probs) const;
	/// <summary>
	/// Will create probability bins using Platt scale method
	/// @param x The training matrix
	/// @param y The Labels
	/// @param poly_rank the polynom rank for the Platt scale fit
	/// @param min_bucket_size The minimal observations to create probability bin
	/// @param min_score_jump The minimal diff in scores to create bin
	/// </summary>
	/// <returns>
	/// @param params Stores the Platt scale model params for conversion
	/// </returns>
	int learn_prob_calibration(MedMat<float> &x, vector<float> &y, int poly_rank, vector<double> &params, int min_bucket_size = 10000, float min_score_jump = 0.001);
	/// <summary>
	/// Converts probability from Platt scale model
	/// </summary>
	template<class T, class L> int convert_scores_to_prob(const vector<T> &preds, const vector<double> &params, vector<L> &converted) const;

	// init
	static MedPredictor *make_predictor(string model_type);
	static MedPredictor *make_predictor(MedPredictorTypes model_type);
	static MedPredictor *make_predictor(string model_type, string params);
	static MedPredictor *make_predictor(MedPredictorTypes model_type, string params);

	/// Prepartion function for fast prediction on single item each time
	virtual void prepare_predict_single() {};
	virtual void predict_single(const vector<float> &x, vector<float> &preds) const;
	virtual void predict_single(const vector<double> &x, vector<double> &preds) const;
	virtual void calc_feature_importance_shap(vector<float> &features_importance_scores, string &importance_type, const MedFeatures *features);

	// (De)Serialize
	ADD_CLASS_NAME(MedPredictor)
		ADD_SERIALIZATION_FUNCS(classifier_type)
		void *new_polymorphic(string derived_class_name);
	size_t get_predictor_size();
	size_t predictor_serialize(unsigned char *blob);


protected:
	// some needed helpers
	void prepare_x_mat(MedMat<float> &x, const vector<float> &wgts, int &nsamples, int &nftrs, bool transpose_needed) const;
	void predict_thread(void *p) const;

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

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) const;

	void normalize_x_and_y(float *x, float *y, const float *w, int nsamples, int nftrs, vector<float>& x_avg, vector<float>& x_std, float& y_avg, float& y_std);
	int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std);
	void print(FILE *fp, const string& prefix, int level = 0) const;

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
	void print(FILE *fp, const string& prefix, int level=0) const;

	ADD_CLASS_NAME(MedLasso)
	ADD_SERIALIZATION_FUNCS(classifier_type, n_ftrs, b0, b)	
};
/// Least Square direct iterations solution
int learn_lm(float *x, float *_y, const float *w, int nsamples, int nftrs, int niter, float eiter, float *rfactors, float *b, float *err, float *corrs);
/// Least Square direct iterations solution
int learn_lm(float *x, float *_y, const float *w, int nsamples, int nftrs, int niter, float eiter, float *rfactors, float *b, float *err, float *corrs, float *sumxx);

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

struct MedQRFParams : public SerializableObject {

	// Required
	int ntrees;
	int maxq;
	int learn_nthreads, predict_nthreads;
	QRF_TreeType type;

	// Optional
	int max_samp; ///<M if > 0 & sampsize is NULL : the maximal sampsize we will take from each category
	float samp_factor; ///< if > 0 & sampsize if NULL : the maximal factor of samples between the 2 largest categories
	vector<int> samp_vec; ///< to be used when sampsize is NULL and max_samp,samp_vector > 0
	int *sampsize;
	int ntry;
	int get_only_this_categ;
	int max_depth; ///<maximial depth of tree branches - if 0 no limit
	bool take_all_samples; ///<use all samples - no sampling in building tree

	// Regression
	float spread;
	bool keep_all_values; ///< For quantile regression
	bool sparse_values; ///< For keeping all values as a value-index(int):count(char) vector

	// categorical
	int min_node;
	int n_categ;

	int collect_oob;

	// For Prediction
	int get_count;
	vector<float> quantiles; ///< For quantile regression

	ADD_CLASS_NAME(MedQRFParams)
		ADD_SERIALIZATION_FUNCS(ntrees, maxq, learn_nthreads, predict_nthreads, type, max_samp, samp_factor, samp_vec,
			ntry, get_only_this_categ, max_depth, take_all_samples, spread, keep_all_values, sparse_values, min_node, n_categ, collect_oob, get_count, quantiles)
		void post_deserialization() { if (samp_vec.size() == 0) sampsize = NULL;  else sampsize = &samp_vec[0]; }

};

class MedQRF : public MedPredictor {
public:
	/// Model 
	QRF_Forest qf;

	/// Parameters
	MedQRFParams params;

	// Function
	MedQRF();
	~MedQRF() {};
	MedQRF(void *params);
	MedQRF(MedQRFParams& params);
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedQRF.cpp MedQRF::init
	virtual int set_params(map<string, string>& mapper);
	//	int init(const string &init_str); // allows init of parameters from a string. Format is: param=val,... , for sampsize: 0 is NULL, a list of values is separated by ; (and not ,)
	void init_defaults();

	/// @snippet MedQRF.cpp MedQRF_get_types
	QRF_TreeType get_tree_type(string name);

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	//int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std) {return 0;};

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	ADD_CLASS_NAME(MedQRF)
		ADD_SERIALIZATION_FUNCS(classifier_type, qf, params, model_features, features_count)

	// Print
	void print(FILE *fp, const string& prefix, int level=0) const;
	void printTrees(const vector<string> &modelSignalNames, const string &outputPath) const;
	void calc_feature_importance(vector<float> &features_importance_scores, const string &general_params, const MedFeatures *features);

	// Predictions per sample
	int n_preds_per_sample() const;

	void prepare_predict_single();
	void predict_single(const vector<float> &x, vector<float> &preds) const;

private:
	void set_sampsize(float *y, int nsamples); // checking if there's a need to prep sampsize based on max_samp and samp_factor
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int get_count) const;

	vector<pair<float, int>> _indexd_quantiles;
	vector<float> _sorted_quantiles;
	qrf_scoring_thread_params _single_pred_args;
	bool prepared_single;
};

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
	int MaxTerms;	///< Maximal number of Terms in final model
	int MaxDegree;	///< Model max linear degree : 1 - linear, 2 - square, 3 - cubic, etc...
	double Penalty;
	double Thresh;
	int MinSpan;
	bool Prune;
	int FastK;
	double FastBeta;
	double NewVarPenalty;
	bool UseBetaCache;
	double Trace;	///< debug prints during algorithm run (recommended): 0: no prints , 3: print all (recommended)
};

class MedMars : public MedPredictor {
public:
	// Model 
	int nMaxTerms;
	int nTerms;
	int nPreds;
	bool *BestSet;			///< size: nMaxTerms
	vector<int> Dirs;		///< size: nMaxTerms*nPreds
	vector<double> Cuts;	///< size: nMaxTerms*nPreds
	vector<double> Betas;	///< size: nMaxTerms

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
	/// The parsed fields from init command.
	/// @snippet MedMars.cpp MedMars::init
	virtual int set_params(map<string, string>& mapper);
	void init_defaults();

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	//int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std) {return 0;};

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	ADD_CLASS_NAME(MedMars)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Print
	void print(FILE *fp, const string& prefix, int level = 0) const;

	// Predictions per sample
	int n_preds_per_sample() const;
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
	/// Loss function
	GBM_LossFunctions loss_function;

	/// Alpha for quantile loss function
	double alpha_quantile;

	/// Model 
	full_gbm_learn_info_t gbm_model;

	/// Parameters
	MedGBMParams params;

	/// Predicting on subset of trees
	int predict_ntrees;

	// Function
	MedGBM();
	MedGBM(void *params);
	MedGBM(MedGBMParams& params);
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedGBM.cpp MedGBM::init
	virtual int set_params(map<string, string>& mapper);
	void init_defaults();
	GBM_LossFunctions get_loss_function(string name);
	~MedGBM();

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	ADD_CLASS_NAME(MedGBM)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Print
	void print(FILE *fp, const string& prefix, int level = 0) const;
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
	/// The parsed fields from init command.
	/// @snippet MedKNN.cpp MedKNN::init
	virtual int set_params(map<string, string>& mapper);
	int init(void *params);
	~MedKNN();
	knnAveraging get_knn_averaging(string name);
	knnMetric get_knn_metric(string name);

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_CLASS_NAME(MedKNN)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
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
	double alpha; ///< learning rate
	double beta;///< parameter of logistic function
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
	/// The parsed fields from init command.
	/// @snippet MedBP.cpp MedBP::init
	virtual int set_params(map<string, string>& mapper);
	~MedBP();

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_CLASS_NAME(MedBP)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
		
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


//================================================================
// Unsupervised
//================================================================

/// K-Means: x is input matrix(each row is sample N*M). K- number of clusters, centers - output centroids of clusters(K*M)
/// clusters - output for each sample the cluster number from 0 to K-1(N*1). 
/// dists - output of distance for each sample form each cluster(N*K)
int KMeans(MedMat<float> &x, int K, MedMat<float> &centers, vector<int> &clusters, MedMat<float> &dists);
/// K-Means: x is input matrix(each row is sample N*M). K- number of clusters, centers - output centroids of clusters(K*M)
/// clusters - output for each sample the cluster number from 0 to K-1(N*1). 
/// dists - output of distance for each sample form each cluster(N*K)
int KMeans(MedMat<float> &x, int K, int max_iter, MedMat<float> &centers, vector<int> &clusters, MedMat<float> &dists);
/// K-Means: x is input matrix(each row is sample N*M). K- number of clusters, centers - output centroids of clusters(K*M)
/// clusters - output for each sample the cluster number from 0 to K-1(N*1). 
/// dists - output of distance for each sample form each cluster(N*K)
int KMeans(float *x, int nrows, int ncols, int K, float *centers, int *clusters, float *dists);

/// K-Means: x is input matrix(each row is sample N*M). K- number of clusters, centers - output centroids of clusters(K*M)
/// clusters - output for each sample the cluster number from 0 to K-1(N*1). 
/// dists - output of distance for each sample form each cluster(N*K)
int KMeans(float *x, int nrows, int ncols, int K, int max_iter, float *centers, int *clusters, float *dists, bool verbose_print = true); // actual implemetation routine

// PCA

/// given a matrix, returns the base PCA matrix and the cummulative relative variance explained by them.
/// it is highly recommended to normalize the input matrix x before calling.
int MedPCA(MedMat<float> &x, MedMat<float> &pca_base, vector<float> &varsum);

/// returns the projection of the pca base on the first dim dimensions.
int MedPCA_project(MedMat<float> &x, MedMat<float> &pca_base, int dim, MedMat<float> &projected);

///wrapper for MedPredictor for certian groups - routes the input to correct model group.
///for example may be used to train specific model for each age group
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

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;
	MedSpecificGroupModels *clone() const;

	ADD_CLASS_NAME(MedSpecificGroupModels)
		ADD_SERIALIZATION_FUNCS(classifier_type, featNum, feat_ths, predictors, model_features, features_count)

		//	void print(FILE *fp, const string& prefix) ;
		// Parameters
		void set_predictors(const vector<MedPredictor *> &predictors); //for each group index
	void set_group_selection(int featNum, const vector<float> &feat_ths);
	MedPredictor *get_model(int ind);
	int model_cnt() const;
private:
	vector<MedPredictor *> predictors;
	int featNum;
	vector<float> feat_ths;
	int selectPredictor(const float *x) const; //retrieve predictor index
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
	/// The parsed fields from init command.
	/// @snippet MedSvm.cpp MedSvm::init
	virtual int set_params(map<string, string>& mapper);
	int init(struct svm_parameter &params);

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_CLASS_NAME(MedSvm)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

private:


};


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
//=========================================================================================

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	*  \brief models namespace
	*/
	namespace models {
		/// \brief returns string to create model with init_string. void * is MedPredictor
		string getParamsInfraModel(void *model);
		/// \brief returns MedPredictor *, a clone copy of given model (params without learned data). if delete_old is true will free old given model
		void *copyInfraModel(void *model, bool delete_old = true);
		/// \brief initialize model which is MedPredictor by copying it's parameters to new address and freeing old one
		void initInfraModel(void *&model);
		/// \brief run Learn on the MedPredictor - wrapper api
		void learnInfraModel(void *model, const vector<vector<float>> &xTrain, vector<float> &y, vector<float> &weights);
		/// \brief run predict on the MedPredictor - wrapper api
		vector<float> predictInfraModel(void *model, const vector<vector<float>> &xTest);
		/// \brief run cross validation where each pid is in diffrent fold and saves the preds.
		void get_pids_cv(MedPredictor *pred, MedFeatures &matrix, int nFolds,
			mt19937 &generator, vector<float> &preds);
		/// \brief run cross validation where each samples can be in diffrent fold and saves the preds.
		void get_cv(MedPredictor *pred, MedFeatures &matrix, int nFolds,
			mt19937 &generator, vector<float> &preds);
	}
	/*!
	*  \brief process namespace
	*/
	namespace process {
		/// \brief compares two matrixes populations. it's also try to seperate between populations
		/// using the predictor parameters if given
		void compare_populations(const MedFeatures &population1, const MedFeatures &population2,
			const string &name1, const string &name2, const string &output_file,
			const string &predictor_type = "", const string &predictor_init = "", int nfolds = 5, int max_learn = 0);
	}
}

//================================================================
// dependent includes
//================================================================

#include "MedBooster.h"
//#include "MedLightGBM.h"



//=================================================================
// Joining the MedSerialize Wagon
//=================================================================
MEDSERIALIZE_SUPPORT(MedPredictor)
MEDSERIALIZE_SUPPORT(MedLM)
MEDSERIALIZE_SUPPORT(MedLasso)
MEDSERIALIZE_SUPPORT(MedGDLMParams)
MEDSERIALIZE_SUPPORT(MedGDLM)
MEDSERIALIZE_SUPPORT(MedQRFParams)
MEDSERIALIZE_SUPPORT(MedQRF)
MEDSERIALIZE_SUPPORT(MedMicNet)
MEDSERIALIZE_SUPPORT(MedBP)
MEDSERIALIZE_SUPPORT(MedMars)
MEDSERIALIZE_SUPPORT(MedKNN)
MEDSERIALIZE_SUPPORT(MedGBM)
MEDSERIALIZE_SUPPORT(MedMultiClass)
MEDSERIALIZE_SUPPORT(MedSpecificGroupModels)
MEDSERIALIZE_SUPPORT(MedTQRF)


#endif