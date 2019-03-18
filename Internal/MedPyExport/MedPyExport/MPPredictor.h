#ifndef __MED__MPPREDICTOR__H__
#define __MED__MPPREDICTOR__H__

#include "MedPyCommon.h"
#include "MPSerializableObject.h"
#include "MPFeatures.h"

class MedPredictor;

class MPPredictorTypes {
public:
	static const string LINEAR_MODEL;
	static const string QRF;
	static const string GBM;
	static const string KNN;
	static const string BP;
	static const string MARS;
	static const string GD_LINEAR;
	static const string MULTI_CLASS;
	static const string XGB;
	static const string LASSO;
	static const string MIC_NET;
	static const string BOOSTER;
	static const string DEEP_BIT;
	static const string LIGHTGBM;
	static const string SPECIFIC_GROUPS_MODELS;
	static const string SVM;
	static const string LINEAR_SGD;
	static const string VW;
	static const string TQRF;
	static const string BART;
};


class MPPredictor {
public:
	MEDPY_IGNORE(MedPredictor* o);
	MPSerializableObject asSerializable();
	
	MPPredictor();
	MPPredictor(string model_type, string params="");

	bool MEDPY_GET_transpose_for_learn();
	void MEDPY_SET_transpose_for_learn(bool new_val);

	bool MEDPY_GET_normalize_for_learn();
	void MEDPY_SET_normalize_for_learn(bool new_val);

	bool MEDPY_GET_normalize_y_for_learn();
	void MEDPY_SET_normalize_y_for_learn(bool new_val);

	bool MEDPY_GET_transpose_for_predict();
	void MEDPY_SET_transpose_for_predict(bool new_val);
	
	bool MEDPY_GET_normalize_for_predict();
	void MEDPY_SET_normalize_for_predict(bool new_val);

	std::vector<std::string> MEDPY_GET_model_features();
	void MEDPY_SET_model_features(std::vector<std::string> new_val);

	int MEDPY_GET_features_count();
	void MEDPY_SET_features_count(int new_val);

	static MPPredictor make_predictor(string model_type);
	static MPPredictor make_predictor(string model_type, string params);

	int learn(MPFeatures& features);
	int predict(MPFeatures& features);

	// Write Predictor to File
	void write_predictor_to_file(string& outFile);
};
/*



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
virtual void print(FILE *fp, const string& prefix) const;

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
const string &general_params) {
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

/// <summary>
/// calibration for probability using training data
/// @param x The training matrix
/// @param y The Labels
/// @param min_bucket_size The minimal observations to create probabilty bin
/// @param min_score_jump The minimal diff in scores to create bin
/// @param min_prob_jump The minimal diff in probabilties to create bin
/// @param fix_prob_order If true will unite bins that are sorted in wrong way
/// </summary>
/// <returns>
/// @param min_range - writes a corresponding vector with minimal score range
/// @param max_range - writes a corresponding vector with maximal score range
/// @param map_prob - writes a corresponding vector with probabilty for score range
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
/// Will create probabilty bins using Platt scale method
/// @param x The training matrix
/// @param y The Labels
/// @param poly_rank the polynom rank for the Platt scale fit
/// @param min_bucket_size The minimal observations to create probabilty bin
/// @param min_score_jump The minimal diff in scores to create bin
/// </summary>
/// <returns>
/// @param params Stores the Platt scale model params for conversion
/// </returns>
int learn_prob_calibration(MedMat<float> &x, vector<float> &y, int poly_rank, vector<double> &params, int min_bucket_size = 10000, float min_score_jump = 0.001);
/// <summary>
/// Converts probabilty from Platt scale model
/// </summary>
template<class T, class L> int convert_scores_to_prob(const vector<T> &preds, const vector<double> &params, vector<L> &converted) const;

// init
static MedPredictor *make_predictor(string model_type);
static MedPredictor *make_predictor(MedPredictorTypes model_type);
static MedPredictor *make_predictor(string model_type, string params);
static MedPredictor *make_predictor(MedPredictorTypes model_type, string params);

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

*/


#endif //__MED__MPPREDICTOR__H__