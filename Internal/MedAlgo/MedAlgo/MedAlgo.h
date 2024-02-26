/// @file
/// MedAlgo - APIs to different algorithms: Linear Models, RF, GBM, KNN, and more
///

#ifndef __MED_ALGO_H__
#define __MED_ALGO_H__

#if __GNUC__  >= 5 || (defined(_MSC_VER) && !defined(_DEBUG))
#define NEW_COMPLIER false
#else
#define NEW_COMPLIER false
#endif

#include <Logger/Logger/Logger.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <MedStat/MedStat/MedStat.h>
#include <MedFeat/MedFeat/MedFeat.h>
#include <QRF/QRF/QRF.h>
#include <micNet/micNet/micNet.h>
#include <string.h>
#include <limits.h>
#include <MedProcessTools/MedProcessTools/MedProcessUtils.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <TQRF/TQRF/TQRF.h>
#include "svm.h"
#include <unordered_map>
#include <random>
#include <map>
#include <string>

// Forward Declaration
class MedFeatures;

#pragma warning(disable: 4297) //disable annoying " function assumed not to throw an exception but does "

using namespace std;

//================================================================================
// MedPredictor - wrapper for classical learn/predict algorithms
//================================================================================

/// @enum
/// Model Types options
typedef enum {
	MODEL_LINEAR_MODEL = 0, ///< to_use:"linear_model" Linear %Model - creates MedLM
	MODEL_QRF = 1, ///< to_use:"qrf" Q-Random-Forest - creates MedQRF
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
	MODEL_TQRF = 18, ///< to_use:"tqrf" TQRF model - creates MedTQRF
	MODEL_BART = 19, ///< to_use:"bart" MedBART model using BART
	MODEL_EXTERNAL_NN = 20, ///< to_use: "external_nn" , initialize a neural net using a layers file. creates MedExternalNN
	MODEL_SIMPLE_ENSEMBLE = 21, ///< to_use: "simple_ensemble" , give 1 or more models to train, and ensemble them with given weights from the user. creates MedSimpleEnsemble
	MODEL_BY_MISSING_VALUES_SUBSET = 22, ///< to_use: "by_missing_value_subset", choosed MedPredictor on subset of the features based on missing values. choose best fit - creates MedPredictorsByMissingValues.
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
	virtual int set_params(map<string, string>& mapper) { fprintf(stderr, "????? Using the base class set_params() ?????\n"); fflush(stderr); return 0; };
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
	virtual int learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts);
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
	virtual int predict(MedMat<float> &x, vector<float> &preds) const;
	int predict(vector<float> &x, vector<float> &preds, int n_samples, int n_ftrs) const;
	int threaded_predict(MedMat<float> &x, vector<float> &preds, int nthreads) const;

	int learn(const MedFeatures& features);
	int learn(const MedFeatures& features, vector<string>& names);
	virtual int predict(MedFeatures& features) const;

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

	virtual void calc_feature_contribs_conditional(MedMat<float> &mat_x_in, unordered_map<string, float> &contiditional_variables, MedMat<float> &mat_x_out, MedMat<float> &mat_contribs)
	{
		string model_name = "model_id=" + to_string(classifier_type);
		if (predictor_type_to_name.find(classifier_type) != predictor_type_to_name.end())
			model_name = predictor_type_to_name[classifier_type];
		throw logic_error("ERROR:: operation calc_feature_contribs_conditional "
			"isn't supported for " + model_name + " yet.");
	}

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
	virtual bool predict_single_not_implemented() { return false; }
	virtual void prepare_predict_single();
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


//=================================================================
// Joining the MedSerialize Wagon
//=================================================================
MEDSERIALIZE_SUPPORT(MedPredictor)

#endif