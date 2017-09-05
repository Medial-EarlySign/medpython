#ifndef __LINEAR_MODEL_H__
#define __LINEAR_MODEL_H__
#include "PredictiveModel.h"
#include "MedAlgo.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"

using namespace std;

class MedLinearModel : public MedPredictor, public PredictiveModel
{
public:
	MedLinearModel(int numOdSignals);
	subGradientFunction getSubGradients();
	subGradientFunction  getSubGradientsAUC();
	subGradientFunction  getSubGradientsSvm();
	double predict(const vector<float> &input);
	void predict(const vector<vector<float>> &inputs, vector<double> &preds);
	PredictiveModel *clone();

	void print(const vector<string> &signalNames);
	void set_normalization(const vector<float> &meanShift, const vector<float> &factors);
	void apply_normalization(vector<vector<float>> &input);
	
	//Set Loss Fucntions to learn:
	double(*loss_function)(const vector<double> &got, const vector<float> &y);
	double(*loss_function_step)(const vector<double> &, const vector<float> &, const vector<double> &);
	int sample_count;
	int tot_steps;
	double learning_rate;
	float block_num;
	bool norm_l1;

	//MedPredictor Api:
	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	ADD_SERIALIZATION_FUNCS(model_params, _meanShift, _factor)

private:
	vector<float> _meanShift;
	vector<float> _factor;
	bool mark_learn_finish;
};

MEDSERIALIZE_SUPPORT(MedLinearModel)
#endif

