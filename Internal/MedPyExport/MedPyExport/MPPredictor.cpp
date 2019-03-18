#include "MPPredictor.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"

const string MPPredictorTypes::LINEAR_MODEL = "linear_model";
const string MPPredictorTypes::QRF = "qrf";
const string MPPredictorTypes::GBM = "gbm";
const string MPPredictorTypes::KNN = "knn";
const string MPPredictorTypes::BP = "BP";
const string MPPredictorTypes::MARS = "";
const string MPPredictorTypes::GD_LINEAR = "gdlm";
const string MPPredictorTypes::MULTI_CLASS = "multi_class";
const string MPPredictorTypes::XGB = "xgb";
const string MPPredictorTypes::LASSO = "lasso";
const string MPPredictorTypes::MIC_NET = "micNet";
const string MPPredictorTypes::BOOSTER = "booster";
const string MPPredictorTypes::DEEP_BIT = "deep_bit";
const string MPPredictorTypes::LIGHTGBM = "lightgbm";
const string MPPredictorTypes::SPECIFIC_GROUPS_MODELS = "multi_models";
const string MPPredictorTypes::SVM = "svm";
const string MPPredictorTypes::LINEAR_SGD = "linear_sgd";
const string MPPredictorTypes::VW = "vw";
const string MPPredictorTypes::TQRF = "tqrf";
const string MPPredictorTypes::BART = "bart";


MPPredictor::MPPredictor() { o = nullptr; };
MPPredictor::MPPredictor(string model_type, string params) {
	if(params!="")
		o = MedPredictor::make_predictor(model_type, params);
	else 
		o = MedPredictor::make_predictor(model_type);
}

MPPredictor MPPredictor::make_predictor(string model_type) { return MPPredictor(model_type); }
MPPredictor MPPredictor::make_predictor(string model_type, string params) { return MPPredictor(model_type, params); }

bool MPPredictor::MEDPY_GET_transpose_for_learn() { return o->transpose_for_learn; };
void MPPredictor::MEDPY_SET_transpose_for_learn(bool new_val) { o->transpose_for_learn = new_val; };

bool MPPredictor::MEDPY_GET_normalize_for_learn() { return o->normalize_for_learn; };
void MPPredictor::MEDPY_SET_normalize_for_learn(bool new_val) { o->normalize_for_learn = new_val; };

bool MPPredictor::MEDPY_GET_normalize_y_for_learn() { return o->normalize_y_for_learn; };
void MPPredictor::MEDPY_SET_normalize_y_for_learn(bool new_val) { o->normalize_y_for_learn = new_val; };

bool MPPredictor::MEDPY_GET_transpose_for_predict() { return o->transpose_for_predict; };
void MPPredictor::MEDPY_SET_transpose_for_predict(bool new_val) { o->transpose_for_predict = new_val; };

bool MPPredictor::MEDPY_GET_normalize_for_predict() { return o->normalize_for_predict; };
void MPPredictor::MEDPY_SET_normalize_for_predict(bool new_val) { o->normalize_for_predict = new_val; };

std::vector<std::string> MPPredictor::MEDPY_GET_model_features() { return o->model_features; };
void MPPredictor::MEDPY_SET_model_features(std::vector<std::string> new_val) { o->model_features = new_val; };

int MPPredictor::MEDPY_GET_features_count() { return o->features_count; };
void MPPredictor::MEDPY_SET_features_count(int new_val) { o->features_count = new_val; };


int MPPredictor::learn(MPFeatures& features) { return o->learn(*(features.o)); };
int MPPredictor::predict(MPFeatures& features) { return o->predict(*(features.o)); };

MPSerializableObject MPPredictor::asSerializable() { return MPSerializableObject(o); }

void MPPredictor::write_predictor_to_file(string& outFile) {
	size_t predictor_size = o->get_predictor_size();
	vector<unsigned char> buffer(predictor_size);
	o->predictor_serialize(&(buffer[0]));

	if (write_binary_data(outFile, &(buffer[0]), predictor_size) != 0)
		throw runtime_error(string("Error writing model to file ")+outFile);
}
