#include <stdexcept>
#include <map>
#include <cmath>
#include <iostream>
#include "MedLinearModel.h"
#include "SGD.h"

#define MED_MISSSING_VALUE -65336

//not in use
double MedLinearModel::predict(const vector<float> &input) {
	if (input.size() != model_params.size() - 1) {
		throw invalid_argument("input has wrong number of signals. expeced" + to_string(model_params.size() - 1) + " got "
			+ to_string(input.size()));
	}
	double res = model_params[0];
	for (size_t i = 0; i < input.size(); ++i)
	{
		res += input[i] * model_params[i + 1];
	}
	//res += model_params[model_params.size() - 1];
	return res;
}

void MedLinearModel::predict(const vector<vector<float>> &inputs, vector<double> &preds) {
	vector<vector<float>> copy_inp;
	const vector<float> *access_arr = inputs.data();
	if (inputs.size() == 0)
		throw invalid_argument("must have at least one signal");
	if (mark_learn_finish && _meanShift.size() > 0) {
		copy_inp = vector<vector<float>>(inputs);
		apply_normalization(copy_inp);
		access_arr = copy_inp.data();
	}
	if (preds.size() < inputs[0].size())
		preds.resize(inputs[0].size());

	for (size_t i = 0; i < preds.size(); ++i)
	{
		preds[i] = model_params[0];
		for (size_t k = 0; k < inputs.size(); ++k)
			preds[i] += access_arr[k][i] * model_params[k + 1];
	}
}

subGradientFunction  MedLinearModel::getSubGradients() {
	//This is subGradient in L2, for other loss function you need to change this function
	subGradientFunction func = [](int ind, const vector<double> &params, const vector<vector<float>> &x, const vector<float> &y) {
		double res = 0;
		for (size_t i = 0; i < y.size(); ++i)
		{
			double productRes = params[0];
			for (size_t k = 0; k < x.size(); ++k)
			{
				productRes += params[k + 1] * x[k][i];
			}
			float x_val = 1;
			if (ind > 0) {
				x_val = x[ind - 1][i];
			}
			res += 2 * (params[ind] * x_val * x_val + (productRes - params[ind] * x_val)*x_val - y[i] * x_val);
		}
		//res /= y.size(); - constant, not needed

		return res;
	};

	return func;
}
#define REG_LAMBDA 1.0
subGradientFunction  MedLinearModel::getSubGradientsAUC() {
	//This is subGradient in L2, for other loss function you need to change this function. (need fix for bias param, W{n+1})
	subGradientFunction func = [](int ind, const vector<double> &params, const vector<vector<float>> &x, const vector<float> &y) {
		if (ind == 0) {
			return (double)0;
		}

		double res = 0;
		map<int, vector<int>> targetToInd;
		for (size_t i = 0; i < y.size(); ++i)
		{
			targetToInd[(int)y[i]].push_back((int)i);
		}
		vector<int> posInds = targetToInd[1];
		vector<int> negInds = targetToInd[0]; //change to -1 if y is given that way
		for (size_t i = 0; i < posInds.size(); ++i)
		{
			int posIndex = posInds[i];
			for (size_t j = 0; j < negInds.size(); ++j)
			{
				int negIndex = negInds[j];
				double sumDiff = 0;
				for (size_t k = 1; k < params.size(); ++k) //param 0 should be zero
				{
					sumDiff += params[k] * (x[k - 1][posIndex] - x[k - 1][negIndex]);
				}
				sumDiff *= 1;
				if (sumDiff > 100) {
					res += 0;
					continue;
				}
				if (sumDiff < -100) {
					res += (x[ind - 1][posIndex] - x[ind - 1][negIndex]);
					continue;
				}
				double divider = 1 + exp(-sumDiff);
				//avoid overflow:
				if (divider < 1e10) {
					divider = divider * divider;
					res += (x[ind - 1][posIndex] - x[ind - 1][negIndex]) * exp(-sumDiff) / divider;
				}

			}

		}
		res /= (posInds.size() * negInds.size());
		res = -res; //because we need to minimize and auc we need to maximize
		res += REG_LAMBDA * 2 * params[ind];  //regularization
		return res;
	};

	return func;
}
subGradientFunction  MedLinearModel::getSubGradientsSvm() {
	//This is subGradient in L2, for other loss function you need to change this function. (need fix for bias param, W{n+1})
	subGradientFunction func = [](int ind, const vector<double> &params, const vector<vector<float>> &x, const vector<float> &y) {
		double res = 0;
		for (size_t i = 0; i < y.size(); ++i)
		{
			double fx = 1; //first param (ind ==0) is bais like x vector has 1 vector;
			if (ind > 0)
				fx = x[ind - 1][i];
			double diff = 1 - ((y[i] > 0) * 2 - 1) * params[ind] * fx;
			if (diff > 0)
				res += 1 - ((y[i] > 0) * 2 - 1) * fx;
		}

		res /= y.size();
		//add regularization:
		double reg = 0;
		if (REG_LAMBDA > 0) {
			double n_params = 0;
			for (size_t i = 0; i < params.size(); ++i)
				n_params += params[i] * params[i];
			n_params = sqrt(n_params);
			reg = -params[ind] / n_params;
		}

		return res + REG_LAMBDA * reg;
	};

	return func;
}
#define SMOOTH false
#define REGULARIZATION_GEOM 0.0
double _linear_loss_target_auc(const vector<double> &preds, const vector<float> &y) {
	map<int, vector<int>> targetToInd;
	for (size_t i = 0; i < y.size(); ++i)
	{
		targetToInd[(int)y[i]].push_back((int)i);
	}
	vector<int> posInds = targetToInd[1];
	vector<int> negInds = targetToInd[0]; //change to -1 if y is given that way
	double res = 0;
	for (size_t i = 0; i < posInds.size(); ++i)
	{
		for (size_t j = 0; j < negInds.size(); ++j)
		{
			double diffScores = preds[posInds[i]] - preds[negInds[j]];
			//res += 1.0 / (1.0 + exp(-10 * diffScores));
			if (SMOOTH) {
				diffScores *= 5;
				if (diffScores > -100 && diffScores < 100)  res += 1.0 / (1.0 + exp(-diffScores)); else if (diffScores >= 100) res += 1;
			}
			else {
				res += diffScores > REGULARIZATION_GEOM;
				res += 0.5*(diffScores == REGULARIZATION_GEOM);
				//res += diffScores;
			}
		}
	}
	res /= (posInds.size() * negInds.size());
	res = -res; //auc needs to be maximize
	return res;
}
double _linear_loss_step_auc(const vector<double> &preds, const vector<float> &y, const vector<double> &params) {
	double res = _linear_loss_target_auc(preds, y);
	double nrm = 0;
	for (size_t i = 0; i < params.size(); ++i)
		nrm += params[i] * params[i];
	nrm /= 2;
	//res += REG_LAMBDA*nrm; //not needed projecting to 1 after each iteration
	return res;
}
double _linear_loss_step_auc_fast(const vector<double> &preds, const vector<float> &y) {
	map<int, vector<int>> targetToInd;
	for (size_t i = 0; i < y.size(); ++i)
	{
		targetToInd[(int)y[i]].push_back((int)i);
	}
	double res = 0;
	//random_device rd;
	//auto gen = mt19937(rd());
	std::default_random_engine gen;
	gen.seed( /* your seed for the RNG goes here */);
	vector<int> posInds = targetToInd[1];
	vector<int> negInds = targetToInd[0]; //change to -1 if y is given that way
	uniform_int_distribution<> pos_rand(0, (int)posInds.size() - 1);
	uniform_int_distribution<> neg_rand(0, (int)negInds.size() - 1);
	for (size_t i = 0; i < 100; ++i)
	{
		int posInd = posInds[pos_rand(gen)];
		int negInd = negInds[neg_rand(gen)];
		double diffScores = preds[posInd] - preds[negInd];
		//res += 1.0 / (1.0 + exp(-diffScores));
		res += diffScores > 0;
		res += (diffScores == 0)* 0.5;
	}
	//res /= 100;
	res = -res; //auc needs to be maximize

	return res;
}
double _linear_loss_target(const vector<double> &preds, const vector<float> &y) {
	double res = 0;
	int totPos = 0;
	float deired_sen = (float)0.75; //Take AUC @ deired_sen (smooth local AUC)
	int local_points_diff = 1;
	for (size_t i = 0; i < y.size(); ++i)
	{
		totPos += int(y[i] > 0);
	}
	if (totPos == 0) {
		//can't be!! 
		res = INFINITY; //exception
		return res; //change to mssing value - take other sample - can't be
	}

	int stopCnt = int(totPos * deired_sen);
	vector<tuple<double, bool>> predSorted((int)preds.size());
	for (size_t i = 0; i < predSorted.size(); ++i)
	{
		predSorted[i] = tuple<double, bool>(preds[i], y[i] > 0);
	}
	sort(predSorted.begin(), predSorted.end());

	//calc AUC on reversed matchY:
	int posCnt = 0;
	int negCnt = 0;
	int bin_size = 0;
	for (int i = (int)predSorted.size() - 1; i >= 0; --i) {
		posCnt += get<1>(predSorted[i]);
		negCnt += !get<1>(predSorted[i]);
		if (posCnt >= stopCnt - local_points_diff) {
			res += float(negCnt) / (y.size() - totPos);
			++bin_size;
		}
		if (posCnt >= stopCnt + local_points_diff) {
			break;
		}
	}
	if (bin_size > 0)
		res /= bin_size;
	else
		res = 0.5; //unknown

	return res;
}
double _linear_loss_step(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params) {
	float loss_val = (float)_linear_loss_target(preds, y);
	double nrm = 0;
	if (REG_LAMBDA > 0) {
		for (size_t i = 0; i < model_params.size(); ++i)
		{
			nrm += model_params[i] * model_params[i];
		}
		nrm /= model_params.size();
	}

	return loss_val + REG_LAMBDA*nrm;
}
double _linear_loss_target_rmse(const vector<double> &preds, const vector<float> &y) {
	double res = 0;
	for (size_t i = 0; i < y.size(); ++i)
	{
		res += (y[i] - preds[i]) * (y[i] - preds[i]);
	}
	res /= y.size();
	res = sqrt(res);
	return res;
}
double _linear_loss_target_svm(const vector<double> &preds, const vector<float> &y) {
	double res = 0;
	for (size_t i = 0; i < y.size(); ++i)
	{
		double diff = 1 - (2 * (y[i] > 0) - 1) * preds[i];
		if (diff > 0)
			res += diff;
	}
	res /= y.size();

	return res; //no reg - maybe count only accourcy beyond 1 and beyond 0
}
double _linear_loss_step_svm(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params) {
	double res = 0;
	for (size_t i = 0; i < y.size(); ++i)
	{
		double diff = 1 - (2 * (y[i] > 0) - 1) * preds[i];
		if (diff > 0)
			res += diff;
	}
	res /= y.size();

	double reg = 0;
	if (REG_LAMBDA > 0) {
		for (size_t i = 0; i < model_params.size(); ++i)
			reg += model_params[i] * model_params[i];
		reg = sqrt(reg);
	}

	return res + REG_LAMBDA * reg;
}
double _linear_loss_step_rmse(const vector<double> &preds, const vector<float> &y, const vector<double> &model_params) {
	double res = 0;
	for (size_t i = 0; i < y.size(); ++i)
	{
		res += (y[i] - preds[i]) * (y[i] - preds[i]);
	}
	res /= y.size();
	res = sqrt(res);

	double reg = 0;
	if (REG_LAMBDA > 0) {
		for (size_t i = 0; i < model_params.size(); ++i)
			reg += model_params[i] * model_params[i];
		reg = sqrt(reg);
	}

	return res + REG_LAMBDA * reg;
}

MedLinearModel::MedLinearModel(int numOdSignals) : PredictiveModel("LinearModel(" + to_string(numOdSignals) + ")") {
	model_params = vector<double>(numOdSignals + 1); //for bias
	mark_learn_finish = false;
	loss_function = _linear_loss_target_auc; //default target function, can be changed programitaclly
	loss_function_step = _linear_loss_step_auc; //default gradient function (with regularzation), can be changed programitaclly

	transpose_for_learn = false;
	normalize_y_for_learn = false;
	transpose_for_predict = false;

	normalize_for_learn = false; //doing internal and not with MedMat to save normalization params for predict
	normalize_for_predict = false;
	classifier_type = MODEL_LINEAR_SGD;
}

void MedLinearModel::print(const vector<string> &signalNames) {
	for (size_t i = 0; i < model_params.size(); ++i)
	{
		if (i == 0) {
			cout << "Param0=" << model_params[i] << endl;
			continue;
		}
		cout << "Param" << i << " " << signalNames[i - 1] << "=" << model_params[i] << endl;
	}
}

void MedLinearModel::set_normalization(const vector<float> &meanShift, const vector<float> &factors) {
	_meanShift = meanShift;
	_factor = factors;
}

void MedLinearModel::apply_normalization(vector<vector<float>> &input) {
	for (size_t i = 0; i < input.size(); ++i)
	{
		for (size_t j = 0; j < input[i].size(); ++j)
		{
			if (input[i][j] != MED_MISSSING_VALUE)
				input[i][j] = (input[i][j] - _meanShift[i]) / _factor[i];
			else
				input[i][j] = 0;
		}
	}
}

PredictiveModel *MedLinearModel::clone() {
	PredictiveModel *copy = new MedLinearModel((int)model_params.size() - 1);
	//dont copy values of params and normalization - not need for now
	return copy;
}

int MedLinearModel::Predict(float *x, float *&preds, int nsamples, int nftrs) {

	for (size_t i = 0; i < nsamples; ++i)
	{
		double p = model_params[0];
		for (size_t k = 0; k < nftrs; ++k) {
			float val = x[i*nftrs + k];
			// has normalization in MedMat - but want to use same from train. when calling this function, it's always need normalizations
			if (val == MED_MISSSING_VALUE)
				val = 0;
			else
				val = (val - _meanShift[k]) / _factor[k];
			p += val * model_params[k + 1];
		}
		preds[i] = (float)p;
	}
	return 0;
}

float _maxDiffVec(const vector<float> &vec) {
	if (vec.size() == 0)
		throw invalid_argument("vector can't  be empty");

	float mmax = vec[0];
	float mmin = vec[0];

	for (size_t i = 0; i < vec.size(); ++i)
	{
		if (mmax < vec[i]) {
			mmax = vec[i];
		}
		if (mmin > vec[i]) {
			mmin = vec[i];
		}
	}

	return mmax - mmin;
}

void _learnModel(SGD &learner, const vector<vector<float>> &xData, const vector<float> &yData, int categoryCnt,
	int T_Steps, int print_steps, double learnRate, int sampleSize) {
	float h_size = (float)0.01;
	int numSteps = T_Steps;
	//float learn_rate = (float)0.0000001;
	float blocking_num = (float)sqrt((int)xData.size()) * 1;
	if (learner.get_blocking() == 0) {
		learner.set_blocking(blocking_num);
	}
	else {
		blocking_num = abs(learner.get_blocking());
	}

	learner.set_gradient_params(sampleSize, h_size, categoryCnt);

	float maxP = _maxDiffVec(xData[0]); //p = the maximal number in data x
	for (size_t i = 1; i < xData.size(); ++i)
	{
		float maxSignal = _maxDiffVec(xData[i]);
		if (maxSignal > maxP) {
			maxP = maxSignal;
		}
	}
	cout << "maxDiffSignal = " << maxP << endl;
	float blockDer = maxP / h_size;
	if (learnRate > 0) {
		learner.set_learing_rate((float)learnRate);
	}
	else {
		learner.set_learing(blocking_num, blockDer, T_Steps);
	}
	//learner.subGradientI = model->getSubGradients();
	//learner.subGradientI = ((LinearModel *)model)->getSubGradientsAUC(); 

	//learner.set_learing_rate(learner.get_learing_rate() / 5);
	if (print_steps > 0) {
		learner.output_num = T_Steps / print_steps;
	}
	cout << "learning_rate = " << learner.get_learing_rate() << ", eppsilon=" << learner.get_learing_eppsilon(blocking_num, blockDer, T_Steps) << endl;

	learner.Learn(xData, yData, numSteps);
}

template<class T> T _avgVec(const vector<T> &vec) {
	T res = 0;
	int sz = 0;
	for (size_t i = 0; i < vec.size(); ++i)
	{
		if (vec[i] == MED_MISSSING_VALUE) {
			continue;
		}
		res += vec[i];
		++sz;
	}
	if (sz == 0) {
		return MED_MISSSING_VALUE;
	}
	return res / sz;
}

float _stdVec(const vector<float> &vec, float avg) {
	float res = 0;
	int sz = 0;
	for (size_t i = 0; i < vec.size(); ++i)
	{
		if (vec[i] == MED_MISSSING_VALUE) {
			continue;
		}
		res += (vec[i] - avg)*(vec[i] - avg);
		++sz;
	}
	if (sz == 0) {
		return MED_MISSSING_VALUE;
	}
	res /= sz;
	res = sqrt(res);

	return res;
}

void _normalizeSignalToAvg(vector<vector<float>> &xData, vector<float> &meanShift, vector<float> &factor) {
	meanShift = vector<float>((int)xData.size());
	factor = vector<float>((int)xData.size());
	for (size_t i = 0; i < xData.size(); ++i)
	{
		//fixOutlyers(xData[i]);
		float avg = _avgVec(xData[i]);
		meanShift[i] = avg;
		float std = _stdVec(xData[i], avg);
		factor[i] = std;
		for (size_t k = 0; k < xData[i].size(); ++k)
		{
			if (xData[i][k] != MED_MISSSING_VALUE && std != 0) {
				xData[i][k] = (xData[i][k] - avg) / std; //z-Score
			}
			else {
				xData[i][k] = 0; //maybe Avg?
			}
		}
	}
}

int MedLinearModel::Learn(float *x, float *y, float *w, int nsamples, int nftrs) {

	vector<float> avg_diff, factors;
	vector<float> yData(y, y + nsamples - 1);
	vector<vector<float>> xData(nftrs);
	for (size_t i = 0; i < nftrs; ++i)
	{
		xData[i].resize(nsamples);
		for (size_t j = 0; j < nsamples; ++j)
			xData[i][j] = x[j* nftrs + i];
	}
	_normalizeSignalToAvg(xData, avg_diff, factors);
	set_normalization(avg_diff, factors);
	SGD learner(this, loss_function);
	learner.subGradientI = NULL; //((LinearModel *)mdl)->getSubGradientsSvm();
	learner.set_blocking((float)1);
	//learner.set_model_precision(1e-5);
	learner.set_special_step_func(loss_function_step); //not in use if learner.subGradientI is not NULL

	int minCat = 1;
	int sampleSize = 500;
	int tot_steps = 10000;
	double learning_rate = 3 * 1e-7;

	mark_learn_finish = false;
	_learnModel(learner, xData, yData, minCat, tot_steps, 10, learning_rate, sampleSize);
	mark_learn_finish = true;

	return 0;
}
