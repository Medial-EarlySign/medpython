#ifndef __MEDXGB_H__
#define __MEDXGB_H__
#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedProcessTools/MedProcessTools/MedProcessUtils.h>
#include <xgboost/learner.h>
#include <xgboost/data.h>
#include <xgboost/c_api.h>
#include "MedProcessTools/MedProcessTools/MedSamples.h"

struct MedXGBParams : public SerializableObject {
	string booster; // gbtree or gblinear
	string objective; // binary:logistic is logistic regression loss function for binary classification
	float eta; // step size shrinkage
	float gamma; // minimum loss reduction required to make a further partition
	int min_child_weight; // minimum sum of instance weight(hessian) needed in a child
	int max_depth; // maximum depth of a tree
	int num_round; // the number of rounds to do boosting
	vector<string> eval_metric; // when not silent, report this metric
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
	int verbose_eval;
	float validate_frac; // how much of the training set is used as validation  for evaluation. should be between 0 and 1.
	string split_penalties; // feature-dependent splitting penalty. string format is "number:value,number:value,..."
	string monotone_constraints; // feature-dependent monotonic constraint. string format is "part_of_feature_name:part_of_feature_name,number:value,..."


	MedXGBParams() {
		booster = "gbtree"; objective = "binary:logistic"; eta = 1.0; gamma = 1.0;
		min_child_weight = 1; max_depth = 3; num_round = 500; silent = 1; eval_metric.push_back("auc"); missing_value = MED_MAT_MISSING_VALUE;
		num_class = 1; //only set when multiclass
		colsample_bytree = 1.0; colsample_bylevel = 1.0; subsample = 1.0; scale_pos_weight = 1.0; tree_method = "auto"; lambda = 1; alpha = 0;
		seed = 0;
		verbose_eval = 0;
		validate_frac = 0;
	}

	ADD_CLASS_NAME(MedXGBParams)
		ADD_SERIALIZATION_FUNCS(booster, objective, eta, gamma, min_child_weight, max_depth, num_round, eval_metric, silent, missing_value, num_class,
			colsample_bytree, colsample_bylevel, subsample, scale_pos_weight, tree_method, lambda, alpha, seed, verbose_eval, validate_frac, split_penalties, monotone_constraints)
};

class XGBBooster {
public:
	explicit XGBBooster(const std::vector<std::shared_ptr<xgboost::DMatrix> >& cache_mats)
		: configured_(false),
		initialized_(false),
		learner_(xgboost::Learner::Create(cache_mats)) {}

	inline xgboost::Learner* learner() {
		return learner_.get();
	}

	inline void SetParam(const std::string& name, const std::string& val) {
		auto it = std::find_if(cfg_.begin(), cfg_.end(),
			[&name, &val](decltype(*cfg_.begin()) &x) {
			if (name == "eval_metric") {
				return x.first == name && x.second == val;
			}
			return x.first == name;
		});
		if (it == cfg_.end()) {
			cfg_.push_back(std::make_pair(name, val));
		}
		else {
			(*it).second = val;
		}
		if (configured_) {
			learner_->SetParams(cfg_);
			learner_->Configure();
		}
	}

	inline void LazyInit() {
		if (!configured_) {
			learner_->SetParams(cfg_);
			learner_->Configure();
			configured_ = true;
		}
		if (!initialized_) {
			//learner_->Configure();
			initialized_ = true;
		}
	}

	inline void LoadModel(dmlc::Stream* fi) {
		learner_->Load(fi);
		initialized_ = true;
	}

public:
	bool configured_;
	bool initialized_;
	std::unique_ptr<xgboost::Learner> learner_;
	std::vector<std::pair<std::string, std::string> > cfg_;
};

class MedXGB : public MedPredictor {
public:
	BoosterHandle my_learner = NULL;
	xgboost::DMatrix* dvalidate = NULL;
	MedXGBParams params;
	void init_defaults();
	int feat_contrib_flags = 0;
	virtual int init(void *classifier_params) { this->params = *((MedXGBParams*)classifier_params); return 0; };
	/// The parsed fields from init command.
	/// @snippet MedXGB.cpp MedXGB::init
	virtual int set_params(map<string, string>& initialization_map);

	// Function
	MedXGB() { init_defaults(); };
	~MedXGB();

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;
	void prepare_mat_handle(float *x, float *y, const float *w, int nsamples, int nftrs, DMatrixHandle &matrix_handle);

	virtual void print(FILE *fp, const string& prefix, int level = 0) const;

	void calc_feature_importance(vector<float> &features_importance_scores,
		const string &general_params, const MedFeatures *features);

	void calc_feature_contribs(MedMat<float> &x, MedMat<float> &contribs);

	void calc_feature_contribs_conditional(MedMat<float> &mat_x_in, unordered_map<string, float> &contiditional_variables, MedMat<float> &mat_x_out, MedMat<float> &mat_contribs);


	void export_predictor(const string &output_fname);

	int n_preds_per_sample() const;

	void pre_serialization() {
		const char* out_dptr;
		xgboost::bst_ulong len;
		string cfg_js = "{ \"format\":\"json\" }";
		if (my_learner != NULL) {
			if (XGBoosterSaveModelToBuffer(my_learner, cfg_js.c_str(), &len, &out_dptr) != 0)
				throw runtime_error("failed XGBoosterGetModelRaw\n");
			serial_xgb.resize(len);
			memcpy(&serial_xgb[0], out_dptr, len);
		}
		else
			serial_xgb.clear();
	}

	void post_deserialization() {
		if (this->my_learner != NULL)
			XGBoosterFree(this->my_learner);
		if (!serial_xgb.empty()) {
			DMatrixHandle h_train_empty[1];
			if (XGBoosterCreate(h_train_empty, 0, &my_learner) != 0)
				throw runtime_error("failed XGBoosterCreate\n");
			if (XGBoosterLoadModelFromBuffer(my_learner, &serial_xgb[0], serial_xgb.size()) != 0)
				throw runtime_error("failed XGBoosterLoadModelFromBuffer\n");
			serial_xgb.clear();
		}
	}

	void prepare_predict_single();
	void predict_single(const vector<float> &x, vector<float> &preds) const;

	void get_json(const char ***json, int& len, string type) {
		if (my_learner != NULL) {
			string no_fmap = "";
			xgboost::bst_ulong _len;
			int succ = XGBoosterDumpModelEx(my_learner, no_fmap.c_str(), 1, type.c_str(), &_len, json);
			if (succ < 0)
				HMTHROW_AND_ERR("Error MedXGB::get_json - can't get model\n");
			len = (int)_len;
		}
		else
			len = 0;
	}

	ADD_CLASS_NAME(MedXGB)
		ADD_SERIALIZATION_FUNCS(classifier_type, serial_xgb, params, model_features, features_count, _mark_learn_done)

private:
	bool _mark_learn_done;
	bool prepared_single;
	vector<BoosterHandle> learner_per_thread;

	void translate_split_penalties(string& split_penalties_s);
	void translate_monotone_constraints(string& monotone_constraints_s);
	void calc_feature_importance_local(vector<float> &features_importance_scores, string &importance_type);
	vector<char> serial_xgb;
};

//=================================================================
// Joining the MedSerialize Wagon
//=================================================================
MEDSERIALIZE_SUPPORT(MedXGBParams)
MEDSERIALIZE_SUPPORT(MedXGB)

//#endif
#endif