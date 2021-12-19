#ifndef __SGD_H__
#define __SGD_H__

#include <MedAlgo/MedAlgo/PredictiveModel.h>
#include <random>
#include <map>

using namespace std;

class SGD
{
public:
	SGD(PredictiveModel *mdl, double(*loss_funct)(const vector<double> &got, const vector<float> &y, const vector<float> *weights));
	void Learn(const vector<vector<float>> &xData, const vector<float> &yData, int T_Steps, const vector<float> *weights = NULL, bool print_auc = false);
	//use one of these techniques to set params for subgradients
	double(*subGradientI)(int param_number, const vector<double> &param_values, const vector<vector<float>> &x, const vector<float> &y, const vector<float> *weights);
	void set_gradient_params(int samplePointCnt, float h, int minSampForCat = 0);

	//use one of them to set learning_rate
	void set_learing_rate(float val);
	void set_learing(float blockVals, float blockDerivate, int T_steps);
	void set_special_step_func(double(*function)(const vector<double> &, const vector<float> &, const vector<double> &, const vector<float> *));

	//use to block parameter search, not mandatory: pass val < 0 to force each step projection to abs(val)
	void set_blocking(float val);

	//use to round model params - act as some regularization 
	void set_model_precision(double val);
	double get_model_precision();

	float get_learing_rate();
	float get_learing_eppsilon(float blockVals, float blockDerivate, int T_steps);
	PredictiveModel *get_model();
	float get_blocking();

	size_t output_num; //The number of lines before printing output, 0= no printing
	bool norm_l1;
private:
	PredictiveModel * _model;
	vector<PredictiveModel *> _models_par;
	double(*loss_function)(const vector<double> &got, const vector<float> &y, const vector<float> *weights);
	double(*step_loss_function)(const vector<double> &got, const vector<float> &y, const vector<double> &model_params, const vector<float> *weights);
	float _learning_rate;
	bool has_learing_rate;
	int _minSampForCat;
	double _min_precision;
	map<float, vector<int>> _categoryIndex;

	float _h;
	int _sampleSize;
	mt19937 _gen;

	float _blocking_val;

	//allocate memory for step function once
	vector<float> _sampleY;
	vector<float> _sampleW;
	vector<vector<float>> _sampleX;
	vector<double> _preds_plus;
	vector<double> _preds_minus;
	vector<double> _preds_base;

	vector<double> _step(const vector<vector<float>> &xData, const vector<float> &yData, const vector<float> *weights);
	void _projection_step(vector<double> &params);
	void _round_step(vector<double> &params);
};

void factor(vector<double> &a, double fact);
void add(vector<double> &a, const vector<double> &b);
double product(const vector<double> &a, const vector<double> &b);

#endif