#ifndef __PREDICTIVE_MODEL_H__
#define __PREDICTIVE_MODEL_H__

#include <vector>
#include <string>

using namespace std;
typedef double(*subGradientFunction)(int, const vector<double> &, const vector<vector<float>> &, const vector<float> &);

/** Predictive Model is abstract class of predictor model which has parameters for GD or SGD uses
* it also has function to retrieve direct sub-gradients for the loss function of the model
*/
class PredictiveModel {

public:
	PredictiveModel(string name); ///<The name of the model
	virtual double predict(const vector<float> &input) = 0;
	virtual subGradientFunction getSubGradients();///<Subgradient function to calc directly the gradient descent
	virtual void predict(const vector<vector<float>> &inputs, vector<double> &preds); ///<virtual to allow more efficeint implemention
	virtual void print(const vector<string> &signalNames) = 0; ///<print model to stdout
	virtual PredictiveModel *clone() = 0;///<copy model

	vector<double> model_params;///<model parameters
	string model_name;///<model name
};

#endif // !__PREDICTIVE_MODEL_H__