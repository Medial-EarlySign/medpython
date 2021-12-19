#ifndef __LINEAR_MODEL_H__
#define __LINEAR_MODEL_H__
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedAlgo/MedAlgo/PredictiveModel.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>

using namespace std;

/**
* Linear Model with customizable SGD support
*/
class MedLinearModel : public MedPredictor, public PredictiveModel
{
public:
	MedLinearModel();
	/// The parsed fields from init command.
	/// @snippet MedLinearModel.cpp MedLinearModel::init
	int init(map<string, string>& mapper) { return set_params(mapper); }
	int set_params(map<string, string>& mapper);

	subGradientFunction getSubGradients(); ///<Subgradient of RMSE loss function
	subGradientFunction  getSubGradientsAUC(); ///<Subgradient of smooth auc loss function
	subGradientFunction  getSubGradientsSvm(); ///<Subgradient of svm loss function
	double predict(const vector<float> &input) const;
	void predict(const vector<vector<float>> &inputs, vector<double> &preds) const;
	PredictiveModel *clone() const;

	void print(const vector<string> &signalNames) const;
	void set_normalization(const vector<float> &meanShift, const vector<float> &factors); ///<Normalization
	void apply_normalization(vector<vector<float>> &input) const; ///<apply Normalization
	void get_normalization(vector<float> &meanShift, vector<float> &factors) const;
	
	//Set Loss Fucntions to learn:
	double(*loss_function)(const vector<double> &got, const vector<float> &y, const vector<float> *weights);///<The custom loss_function
	///The custom loss_function step for sgd
	double(*loss_function_step)(const vector<double> &, const vector<float> &, const vector<double> &, const vector<float> *);
	int sample_count; ///<The sample count of sgd
	int tot_steps; ///<The total iteration count of sgd
	double learning_rate; ///<The learning rate  of sgd
	float block_num; ///<The blocking norm for parameter search in sgd
	bool norm_l1; ///<The blocking norm should be n1 or n2?
	int print_steps; ///< how many prints for learn to print
	bool print_model; ///< If true will print model coeff in the end
	int poly_degree; ///< add polynom degree
	int min_cat; ///< control minimal samples per categ in categories
	bool print_auc = false; ///< internal var if o print also auc metric
	bool normalize = true; ///< iof true will do internal method for normalizing

	//MedPredictor Api:
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	void calc_feature_importance(vector<float> &features_importance_scores,
		const string &general_params, const MedFeatures *features);

	ADD_CLASS_NAME(MedLinearModel)
	ADD_SERIALIZATION_FUNCS(classifier_type, model_params, _meanShift, _factor, model_features, features_count, normalize, poly_degree)

private:
	vector<float> _meanShift;
	vector<float> _factor;
	bool mark_learn_finish;
};


double _linear_loss_target_auc(const vector<double> &preds, const vector<float> &y, const vector<float> *weights);
double _linear_loss_step_auc(const vector<double> &preds, const vector<float> &y, const vector<double> &params, const vector<float> *weights);
double _linear_loss_step_auc_fast(const vector<double> &preds, const vector<float> &y, const vector<float> *weights);
double _linear_loss_target_work_point(const vector<double> &preds, const vector<float> &y, const vector<float> *weights);
double _linear_loss_step_work_point(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params, const vector<float> *weights);
double _linear_loss_target_rmse(const vector<double> &preds, const vector<float> &y, const vector<float> *weights);
double _linear_loss_step_rmse(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params, const vector<float> *weights);
double _linear_loss_target_svm(const vector<double> &preds, const vector<float> &y, const vector<float> *weights);
double _linear_loss_step_svm(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params, const vector<float> *weights);

MEDSERIALIZE_SUPPORT(MedLinearModel)

#endif

