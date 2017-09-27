#ifndef __PREDICTIVE_MODEL_H__
#define __PREDICTIVE_MODEL_H__

//Predictive Model is abstract class of predictor model which has parameters for GD or SGD uses
// it also has function to retrieve direct sub-gradients for the loss function of the model

#include <vector>
#include <string>

using namespace std;
typedef double(*subGradientFunction)(int, const vector<double> &, const vector<vector<float>> &, const vector<float> &);

class PredictiveModel {

public:
	PredictiveModel(string name);
	virtual double predict(const vector<float> &input) = 0;
	virtual subGradientFunction getSubGradients();
	virtual void predict(const vector<vector<float>> &inputs, vector<double> &preds); //virtual to allow more efficeint implemention
	virtual void print(const vector<string> &signalNames) = 0;
	virtual PredictiveModel *clone() = 0;

	vector<double> model_params;
	string model_name;
};

#endif // !__PREDICTIVE_MODEL_H__