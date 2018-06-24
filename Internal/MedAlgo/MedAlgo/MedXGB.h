#ifndef __MEDXGB_H__
#define __MEDXGB_H__
#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedProcessTools/MedProcessTools/MedProcessUtils.h>
#include <xgboost/learner.h>
#include <xgboost/data.h>
#include <xgboost/c_api.h>

struct MedXGBParams {
	string booster; // gbtree or gblinear
	string objective; // binary:logistic is logistic regression loss function for binary classification
	float eta; // step size shrinkage
	float gamma; // minimum loss reduction required to make a further partition
	int min_child_weight; // minimum sum of instance weight(hessian) needed in a child
	int max_depth; // maximum depth of a tree
	int num_round; // the number of rounds to do boosting
	string eval_metric; // when not silent, report this metric
	int silent; // debug mode
	float missing_value; // which value in the input is representing missing
	int num_class; // needed for multi:softmax
	float colsample_bytree;
	float colsample_bylevel;
	float subsample;
	float scale_pos_weight;
	string tree_method;
	float lambda;
	float alpha;
	int seed; // randomization seed
	string split_penalties; // feature-dependent splitting penalty. string format is "number:value,number:value,..."


	MedXGBParams() {
		booster = "gbtree"; objective = "binary:logistic"; eta = 1.0; gamma = 1.0;
		min_child_weight = 1; max_depth = 3; num_round = 500; silent = 1; eval_metric = "auc"; missing_value = MED_DEFAULT_MISSING_VALUE;
		num_class = 2;
		colsample_bytree = 1.0; colsample_bylevel= 1.0; subsample = 1.0; scale_pos_weight = 1.0; tree_method = "auto"; lambda = 1; alpha = 0;
		seed = 0;
	}
};

class MedXGB : public MedPredictor {
public:
	BoosterHandle my_learner = NULL;
	xgboost::DMatrix* dvalidate = NULL;
	MedXGBParams params;
	void init_defaults();
	virtual int init(void *classifier_params) { this->params = *((MedXGBParams*)classifier_params); return 0; };
	/// The parsed fields from init command.
	/// @snippet MedXGB.cpp MedXGB::init
	virtual int init(map<string, string>& initialization_map);

	// Function
	MedXGB() { init_defaults(); };
	~MedXGB() {};

	int validate_me_while_learning(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	virtual void print(FILE *fp, const string& prefix);

	void calc_feature_importance(vector<float> &features_importance_scores,
		const string &general_params);
	virtual void calc_feature_contribs(MedMat<float> &x, MedMat<float> &contribs);

	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

private:
	bool _mark_learn_done;

	void translate_split_penalties(string& split_penalties_s);
};

//=================================================================
// Joining the MedSerialize Wagon
//=================================================================

MEDSERIALIZE_SUPPORT(MedXGB)

//#endif
#endif