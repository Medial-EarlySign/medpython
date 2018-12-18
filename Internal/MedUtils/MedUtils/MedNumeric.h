// General numeric macros and functions

#ifndef __MED_NUMERIC_H__
#define __MED_NUMERIC_H__

#include <vector>

using namespace std;

//==========================================================================
// Simple Fit of parameters to general 1d->1d functions
//==========================================================================
class MedFitParams {
	public:
		int n_params;
		vector<pair<float,float>> range;
		vector<float> jump;
};

double square_dist_for_fit(float (*f_to_fit)(vector<float> &p, float xp), vector<float> &p, vector<float> &x, vector<float> &y, vector<float> &w);
int full_enumaration_fit(float (*f_to_fit)(vector<float> &p, float xp), MedFitParams &params, vector<float> &x, vector<float> &y, vector<float> &w, vector<float> &res, double &res_err);

// some usefull fitting functions
float scaled_normal_dist(vector<float> &p, float xp);
float scaled_skewed_normal_dist(vector<float> &p, float xp);
float scaled_log_normal_dist(vector<float> &p, float xp);
float scaled_inv_gauss_dist(vector<float> &p, float xp);

float normal_cdf(float mean, float std, float x);
float log_normal_cdf(float mean, float std, float x);

int get_normal_dist_quantiles(float mean, float sd, vector<float> &q, vector<float> &qv);
int get_skewed_normal_dist_quantiles(float mean, float sd, float shape, vector<float> &q, vector<float> &qv);
int get_log_normal_dist_quantiles(float mean, float sd, vector<float> &q, vector<float> &qv);
int get_inv_gauss_dist_quantiles(float mean, float sd, vector<float> &q, vector<float> &qv);

#endif