#include "SGD.h"
#include <iostream>
#include "MedUtils/MedUtils/MedPlot.h"
#include <ctime>
#include <cmath>
#include <omp.h>
#include <random>
#include <algorithm>
#include <functional>
#include "MedStat/MedStat/MedStat.h"

#include <fenv.h>
#ifndef  __unix__
#pragma float_control( except, on )
#endif

#define IMPROVE_LINEAR_SPEED

SGD::SGD(PredictiveModel *mdl, double(*loss_funct)(const vector<double> &got, const vector<float> &y, const vector<float> *weights)) {
#if defined(__unix__)
	#if defined(SO_COMPILATION)
		#pragma message ( "You are compiling shared library - no exceptions in floating points" )
	#else
		feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	#endif
#endif
	this->_model = mdl;
	this->loss_function = loss_funct;
	this->step_loss_function = NULL;
	subGradientI = NULL;
	has_learing_rate = false;

	_h = 0;
	_sampleSize = 0;

	random_device rd;
	_gen = mt19937(rd());

	_blocking_val = 0;
	_min_precision = 0;
	output_num = 0;
	norm_l1 = false;
}

void SGD::set_learing_rate(float val) {
	if (val <= 0) {
		throw invalid_argument("learning rate must be positive");
	}
	_learning_rate = val;
	has_learing_rate = true;
}
void SGD::set_learing(float blockVals, float blockDerivate, int T_steps) {
	if (T_steps <= 0) {
		throw invalid_argument("T_Step must be positive");
	}
	_learning_rate = (abs(blockVals) / abs(blockDerivate)) / (float)sqrt(T_steps);
	has_learing_rate = true;
}
float SGD::get_learing_rate() {
	return _learning_rate;
}
void SGD::set_special_step_func(double(*function)(const vector<double> &, const vector<float> &, const vector<double> &, const vector<float> *)) {
	step_loss_function = function;
}

float SGD::get_learing_eppsilon(float blockVals, float blockDerivate, int T_steps) {
	return abs(blockVals)  * abs(blockDerivate) / (float)sqrt(T_steps);
}

double product(const vector<double> &a, const vector<double> &b) {
	double res = 0;
	if (a.size() != b.size()) {
		throw invalid_argument("a and b vectors must have same size");
	}
	for (size_t i = 0; i < a.size(); ++i)
	{
		res += a[i] * b[i];
	}
	return res;
}

void add(vector<double> &a, const vector<double> &b) {
	if (a.size() != b.size()) {
		throw invalid_argument("a and b vectors must have same size");
	}
	for (size_t i = 0; i < a.size(); ++i)
	{
		a[i] = a[i] + b[i];
	}
}

void factor(vector<double> &a, double fact) {
	for (size_t i = 0; i < a.size(); ++i)
	{
		a[i] = a[i] * fact;
	}
}

double vectorNorm(const vector<double> &v, bool norm_l1 = false) {
	double res = 0;
	if (!norm_l1) {
		for (size_t i = 0; i < v.size(); ++i)
		{
			res += v[i] * v[i];
		}
		res = sqrt(res);
	}
	else {
		for (size_t i = 0; i < v.size(); ++i)
		{
			res += abs(v[i]);
		}
	}
	return res;
}

void SGD::set_gradient_params(int samplePointCnt, float h, int minSampForCat) {
	if (h <= 0) {
		throw invalid_argument("h param for derivate must be positive");
	}
	if (samplePointCnt < 0) {
		throw invalid_argument("samplePointCnt param for derivate must be positive");
	}

	_h = h;
	_sampleSize = samplePointCnt;
	_minSampForCat = minSampForCat;
}

void SGD::set_blocking(float val) {
	_blocking_val = val;
}

float SGD::get_blocking() {
	return _blocking_val;
}

vector<double> SGD::_step(const vector<vector<float>> &xData, const vector<float> &yData, const vector<float> *weights) {
	if (!has_learing_rate) {
		throw logic_error("please initialize learing rate using: set_learing or set_learing_rate directly");
	}
	if (xData.size() == 0) {
		throw invalid_argument("xData must be at least in size 1");
	}
	int dataPointsCnt = (int)xData[0].size();

	vector<double> currentGrad(_model->model_params.size());

	uniform_int_distribution<> random_num(1, dataPointsCnt);
	auto gen_full = std::bind(random_num, _gen);
	if (subGradientI != NULL) {
		if (_sampleSize == 0 || _sampleSize >= yData.size()) {
#pragma omp parallel for
			for (int i = 0; i < _model->model_params.size(); ++i)
			{
				double res = subGradientI((int)i, _model->model_params, xData, yData, weights);
#pragma omp critical
				currentGrad[i] = res;
			}
		}
		else {
			int addedCnt = 0;
			if (_minSampForCat > 0) { //statsfy _min category if have
				for (auto it = _categoryIndex.begin(); it != _categoryIndex.end(); ++it)
				{
					uniform_int_distribution<> category_num(0, (int)it->second.size() - 1);
					vector<int> allRandomInds(_minSampForCat);
					auto gen = std::bind(category_num, _gen);
					generate(allRandomInds.begin(), allRandomInds.end(), gen);
					for (size_t m = 0; m < _minSampForCat; ++m)
					{
						//int randomIndex = categoryIndexes[category_num(_gen)];
						int randomIndex = it->second[allRandomInds[m]];

						float expectedVal = yData[randomIndex];
						_sampleY[addedCnt] = expectedVal;
						if (weights != NULL && !weights->empty())
							_sampleW[addedCnt] = weights->at(randomIndex);
						for (size_t k = 0; k < xData.size(); ++k)
						{
							_sampleX[k][addedCnt] = xData[k][randomIndex];
						}
						++addedCnt;
					}
				}
			}

			vector<int> restInds(_sampleSize - addedCnt);
			generate(restInds.begin(), restInds.end(), gen_full);
			for (size_t k = addedCnt; k < _sampleSize; ++k)
			{
				//may choose same data point more than once
				//int num = random_num(_gen) - 1;
				int num = restInds[k - addedCnt] - 1;
				float expectedVal = yData[num];
				_sampleY[k] = expectedVal;
				if (weights != NULL && !weights->empty())
					_sampleW[k] = weights->at(num);
				for (size_t m = 0; m < xData.size(); ++m)
				{
					_sampleX[m][k] = xData[m][num];
				}
			}

#pragma omp parallel for
			for (int i = 0; i < _model->model_params.size(); ++i)
			{
				//clock_t start = clock();
				double res = subGradientI((int)i, _model->model_params, _sampleX, _sampleY, &_sampleW);
#pragma omp critical
				currentGrad[i] = res;
				//float dur = (clock() - start) / (float)CLOCKS_PER_SEC;
				//cout << "took " << dur << endl;
			}

		}

	}
	else {
		if (_h == 0) {
			throw logic_error("please initialize h & sampleSize when running SGD without specifiying subGradientI directly");
		}
		if (step_loss_function == NULL) {
			throw logic_error("please initialize loss_fuction when using SGD witout specifiying subGradientI directly");
		}

		//calc grad by exp & update currentGrad:

		int sampleSize = _sampleSize;

		if (sampleSize == 0) {
			sampleSize = (int)yData.size();
		}
		else {
			int addedCnt = 0;
			//statisfy _min category if have
			if (_minSampForCat > 0) {
				for (auto it = _categoryIndex.begin(); it != _categoryIndex.end(); ++it)
				{
					uniform_int_distribution<> category_num(0, (int)it->second.size() - 1);
					for (size_t m = 0; m < _minSampForCat; ++m)
					{
						int randomIndex = it->second[category_num(_gen)];
						float expectedVal = yData[randomIndex];
						_sampleY[addedCnt] = expectedVal;
						if (weights != NULL && !weights->empty())
							_sampleW[addedCnt] = weights->at(randomIndex);
						for (size_t k = 0; k < xData.size(); ++k)
							_sampleX[k][addedCnt] = xData[k][randomIndex];

						++addedCnt;
					}
				}
			}
			for (size_t k = addedCnt; k < _sampleSize; ++k)
			{
				//may choose same data point more than once
				int num = random_num(_gen) - 1;
				float expectedVal = yData[num];
				_sampleY[k] = expectedVal;
				if (weights != NULL && !weights->empty())
					_sampleW[k] = weights->at(num);
				for (size_t m = 0; m < xData.size(); ++m)
					_sampleX[m][k] = xData[m][num];
			}
		}

#if defined(IMPROVE_LINEAR_SPEED)
		//calc numeric derivate for each variable:
		_model->predict(_sampleX, _preds_base);
		//now just adjust scores by the changed param:
		for (int i = 0; i < _model->model_params.size(); ++i)
		{
			currentGrad[i] = 0;
			//iterate on samples for this variable:
			if (i == 0)  //free param
				for (size_t j = 0; j < _preds_plus.size(); ++j) {
					_preds_plus[j] = _preds_base[j] + _h / 2;
					_preds_minus[j] = _preds_base[j] - _h / 2;
				}
			else
				for (size_t j = 0; j < _preds_plus.size(); ++j) {
					_preds_plus[j] = _preds_base[j] + (_h / 2 * _sampleX[i - 1][j]);
					_preds_minus[j] = _preds_base[j] - (_h / 2 * _sampleX[i - 1][j]);
				}

			_model->model_params[i] = _model->model_params[i] + (_h / 2);
			double lossVal_plus = step_loss_function(_preds_plus, _sampleY, _model->model_params, &_sampleW); //can also improve speed by calcing only diff
			_model->model_params[i] = _model->model_params[i] - _h;
			double lossVal_minus = step_loss_function(_preds_minus, _sampleY, _model->model_params, &_sampleW);
			_model->model_params[i] = _model->model_params[i] + (_h / 2);

			double subgradientSampleVal = (lossVal_plus - lossVal_minus) / _h;
			currentGrad[i] += subgradientSampleVal;
			currentGrad[i] /= sampleSize;
		}
#else
		//calc numeric derivate for each variable:
		for (int i = 0; i < _model->model_params.size(); ++i)
		{
			currentGrad[i] = 0;
			//iterate on sample for this variable:
			_model->model_params[i] = _model->model_params[i] + (_h / 2);
			_model->predict(_sampleX, _preds_plus);
			double lossVal_plus = step_loss_function(_preds_plus, _sampleY, _model->model_params, &_sampleW);
			_model->model_params[i] = _model->model_params[i] - _h;
			_model->predict(_sampleX, _preds_minus);
			double lossVal_minus = step_loss_function(_preds_minus, _sampleY, _model->model_params, &_sampleW);
			_model->model_params[i] = _model->model_params[i] + (_h / 2);

			double subgradientSampleVal = (lossVal_plus - lossVal_minus) / _h;
			currentGrad[i] += subgradientSampleVal;
			currentGrad[i] /= sampleSize;
			//TODO: calc STD and median - maybe it's not stable ? alert or show this data?
		}
#endif

	}


	for (size_t i = 0; i < _model->model_params.size(); ++i)
		_model->model_params[i] = _model->model_params[i] - _learning_rate * currentGrad[i];

	return _model->model_params;
}

void updateVector(vector<double> &target, const vector<double> &src) {
	if (src.size() != target.size()) {
		throw invalid_argument("vectors must be same size");
	}
	for (size_t i = 0; i < target.size(); ++i)
	{
		target[i] = src[i];
	}
}

void SGD::_projection_step(vector<double> &params) {
	if (_blocking_val == 0) {
		return;
	}
	//search for closet solution till _blocking_val norm size in L2
	double currentNorm = vectorNorm(params, norm_l1);
	if (_blocking_val > 0 && currentNorm <= _blocking_val) {
		return;
	}
	factor(params, abs(_blocking_val) / currentNorm); //this is the optimal solution that holds ||W|| <= _blocking_val & ||W-W*|| is minimal. connecting a line to central of 
													  //a ball and the intersection with the ball is the closet point to W*. this is geometric solution
}

void SGD::set_model_precision(double val) {
	if (val < 0) {
		throw invalid_argument("val must be positive");
	}
	_min_precision = val;
}

double SGD::get_model_precision() {
	return _min_precision;
}

void SGD::_round_step(vector<double> &params) {
	if (_min_precision == 0) {
		return;
	}

	for (size_t i = 0; i < params.size(); ++i)
	{
		params[i] = round(params[i] / _min_precision) * _min_precision;
	}
}

vector<int> randomGroup(int grp_size, int max_ind) {
	vector<int> res(grp_size);
	vector<bool> in_groop(max_ind, false);

	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> int_dis(0, max_ind - 1);

	int i = 0;
	while (i < grp_size) {
		int num = int_dis(gen);
		if (in_groop[num])
			continue;
		in_groop[num] = true;
		res[i] = num;
		++i;
	}

	return res;
}

void SGD::Learn(const vector<vector<float>> &xData, const vector<float> &yData, int T_Steps, const vector<float> *weights, bool print_auc) {
	vector<double> Wt;
	vector<double> W_final(_model->model_params.size());
	size_t ii = 0;
	int max_sample_size = 10000000;
	vector<int> inds_selected;
	vector<double> modelRes;
	if (yData.size() >= max_sample_size * 10)
		inds_selected = randomGroup(max_sample_size, (int)yData.size());

	if (_minSampForCat > 0) { //categorical
		for (size_t i = 0; i < yData.size(); ++i)
		{
			_categoryIndex[yData[i]].push_back(int(i));
		}
		if (_sampleSize > 0 && _categoryIndex.size() * _minSampForCat >= _sampleSize) {
			throw invalid_argument("Ilegal use of minSampForCat, sampleSize (please increase sampleSize or decrease minSampForCat)");
		}
	}
	if (_sampleSize > 0) {
		_sampleY.resize(_sampleSize); //take only random samples
		if (weights != NULL && !weights->empty())
			_sampleW.resize(_sampleSize, 1);
		_sampleX.resize((int)xData.size());
		for (size_t k = 0; k < xData.size(); ++k)
			_sampleX[k] = vector<float>(_sampleSize);
		_preds_plus.resize(_sampleSize);
		_preds_minus.resize(_sampleSize);
#if defined(IMPROVE_LINEAR_SPEED)
		_preds_base.resize(_sampleSize);
#endif
	}
	else {
		//preformance imporvment with using pointers aren't supported for this option right now
		_sampleY = yData;
		_sampleX = xData;
		if (weights != NULL && !weights->empty())
			_sampleW = *weights;
		_preds_plus.resize(yData.size());
		_preds_minus.resize(yData.size());
#if defined(IMPROVE_LINEAR_SPEED)
		_preds_base.resize(yData.size());
#endif
	}

	time_t start = time(NULL);
	double duration;
	float avgLoss, prev_loss = 0;
	bool firstTime = true;
	for (size_t i = 0; i < T_Steps; ++i)
	{
		Wt = _step(xData, yData, weights);
		_projection_step(Wt);
		_round_step(Wt);
		add(W_final, Wt);
		updateVector(_model->model_params, Wt);
		++ii;
		if (output_num > 0 && ii % output_num == 0) {
			duration = difftime(time(NULL), start);
			start = time(NULL);
			cout << "Done " << ii << " Iterations in " << float2Str((float)duration) << " seconds" << endl;

			if (yData.size() >= max_sample_size * 10) {
				//inds_selected = randomGroup(max_sample_size, yData.size()); //did once no need to do anymore
				//commit selection:
				vector<float> yf(max_sample_size);
				vector<vector<float>> xf((int)xData.size());
				vector<float> wf;
				if (weights != NULL && !weights->empty())
					wf.resize(max_sample_size, 1);
				for (size_t kk = 0; kk < xData.size(); ++kk)
					xf[kk].resize(max_sample_size);
				for (size_t k = 0; k < inds_selected.size(); ++k)
				{
					yf[k] = yData[inds_selected[k]];
					for (size_t kk = 0; kk < xData.size(); ++kk)
						xf[kk][k] = xData[kk][inds_selected[k]];
					if (weights != NULL && !weights->empty())
						wf[k] = weights->at(inds_selected[k]);
				}

				_model->predict(xf, modelRes);
				avgLoss = (float)loss_function(modelRes, yf, &wf);
				if (step_loss_function != NULL) {
					float avg_loss_step = (float)step_loss_function(modelRes, yf, _model->model_params, &wf);

					if (print_auc) {
						float auc_val = medial::performance::auc_q(modelRes, yf, &wf);
						cout << "Learned Model \"" << _model->model_name
							<< "\" with average loss of " << float2Str(avgLoss) << " step loss " << float2Str(avg_loss_step) << " AUC " << auc_val << endl;
					}
					else
						cout << "Learned Model \"" << _model->model_name
						<< "\" with average loss of " << float2Str(avgLoss) << " step loss " << float2Str(avg_loss_step) << endl;
				}
				else
					cout << "Learned Model \"" << _model->model_name << "\" with average loss of " << float2Str(avgLoss) << endl;
			}
			else {
				_model->predict(xData, modelRes);
				avgLoss = (float)loss_function(modelRes, yData, weights);
				if (step_loss_function != NULL) {
					float avg_loss_step = (float)step_loss_function(modelRes, yData, _model->model_params, weights);
					if (print_auc) {
						float auc_val = medial::performance::auc_q(modelRes, yData, weights);
						cout << "Learned Model \"" << _model->model_name
							<< "\" with average loss of " << float2Str(avgLoss) << " step loss " << float2Str(avg_loss_step) << " AUC " << auc_val << endl;
					}
					else
						cout << "Learned Model \"" << _model->model_name
						<< "\" with average loss of " << float2Str(avgLoss) << " step loss " << float2Str(avg_loss_step) << endl;
				}
				else
					cout << "Learned Model \"" << _model->model_name << "\" with average loss of " << float2Str(avgLoss) << endl;
			}
			if (!firstTime && abs(avgLoss - prev_loss) < 1e-5) {
				cout << "loss is stable - stopping loop" << endl;
				break;
			}
			firstTime = false;
			prev_loss = avgLoss;
		}
	}
	duration = difftime(time(NULL), start);
	if (output_num > 0)
		cout << "Finished in " << float2Str((float)duration) << " seconds" << endl;
	/*factor(W_final, 1 / (float)T_Steps);
	updateVector(_model->model_params, W_final);*/
}

PredictiveModel * SGD::get_model() {
	return _model;
}