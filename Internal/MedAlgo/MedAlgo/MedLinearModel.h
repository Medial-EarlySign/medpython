#ifndef __LINEAR_MODEL_H__
#define __LINEAR_MODEL_H__
#include <MedProcessTools/MedProcessTools/SerializableObject.h>
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
	int init(map<string, string>& mapper) { set_params(mapper); }
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
	double(*loss_function)(const vector<double> &got, const vector<float> &y);///<The custom loss_function
	///The custom loss_function step for sgd
	double(*loss_function_step)(const vector<double> &, const vector<float> &, const vector<double> &);
	int sample_count; ///<The sample count of sgd
	int tot_steps; ///<The total iteration count of sgd
	double learning_rate; ///<The learning rate  of sgd
	float block_num; ///<The blocking norm for parameter search in sgd
	bool norm_l1; ///<The blocking norm should be n1 or n2?

	//MedPredictor Api:
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_SERIALIZATION_FUNCS(model_params, _meanShift, _factor, model_features, features_count)

private:
	vector<float> _meanShift;
	vector<float> _factor;
	bool mark_learn_finish;
};


double _linear_loss_target_auc(const vector<double> &preds, const vector<float> &y);
double _linear_loss_step_auc(const vector<double> &preds, const vector<float> &y, const vector<double> &params);
double _linear_loss_step_auc_fast(const vector<double> &preds, const vector<float> &y);
double _linear_loss_target_work_point(const vector<double> &preds, const vector<float> &y);
double _linear_loss_step_work_point(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params);
double _linear_loss_target_rmse(const vector<double> &preds, const vector<float> &y);
double _linear_loss_step_rmse(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params);
double _linear_loss_target_svm(const vector<double> &preds, const vector<float> &y);
double _linear_loss_step_svm(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params);

MEDSERIALIZE_SUPPORT(MedLinearModel)

#endif

