#include <stdexcept>
#include "PredictiveModel.h"

PredictiveModel::PredictiveModel(string name) {
	this->model_name = name;
}

vector<float> takeSignalsForRow(const vector<vector<float>> &inputs, int ind) {
	if (inputs.size() == 0) {
		throw invalid_argument("must have at least one signal");
	}
	vector<float> res(inputs.size());
	for (size_t i = 0; i < inputs.size(); ++i)
	{
		res[i] = inputs[i][ind];
	}
	return res;
}

 void PredictiveModel::predict(const vector<vector<float>> &inputs, vector<double> &preds) const {
	if (inputs.size() == 0) 
		throw invalid_argument("must have at least one signal");
	if (preds.size() < inputs[0].size()) 
		preds.resize(inputs[0].size());
	
	for (size_t i = 0; i < preds.size(); ++i)
	{
		vector<float> inp = takeSignalsForRow(inputs, (int)i);
		double y = this->predict(inp);
		preds[i] = y;
	}

}

subGradientFunction PredictiveModel::getSubGradients() {
	return NULL;
}