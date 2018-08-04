#pragma once

#ifndef __MED_LIGHT_GBM__
#define __MED_LIGHT_GBM__

//========================================================================================================
// MedLightGBM
//
// Wrapping the LightGBM package into a MedPredictor.
// The package code is in Libs/External/LightGBM.
// 
//========================================================================================================

#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <LightGBM/application.h>
#include <LightGBM/dataset.h>
#include <LightGBM/boosting.h>
#include <LightGBM/objective_function.h>
#include <LightGBM/metric.h>
#include <LightGBM/config.h>
#include <LightGBM/LightGBM/src/boosting/gbdt.h>
#include <LightGBM/LightGBM/src/boosting/dart.hpp>
#include <LightGBM/LightGBM/src/boosting/goss.hpp>


//==================================================================
// Wrapper for LightGBM::Application to handle our special cases of
// loading data from memory, etc...
//==================================================================
namespace LightGBM {
	class DatasetLoader;
	class Dataset;
	class Boosting;
	class ObjectiveFunction;
	class Metric;

	class MemApp : public Application {
	public:
		MemApp(int argc, char **argv) : Application::Application(argc, argv) {};
		MemApp() : Application::Application(0, NULL) {}
		//~MemApp() { Application::~Application(); };

		int init(map<string, string>& initialization_map);

		// train
		int InitTrain(float *xdata, float *ydata, const float *weight, int nrows, int ncols);
		void Train();

		// predict
		void Predict(float *x, int nrows, int ncols, float *&preds) const;

		// initializing the train_data_ object from a float c matrix
		int InitTrainData(float *xdata, float *ydata, const float *weight, int nrows, int ncols);

		// string serializations
		int serialize_to_string(string &str) { str = boosting_->SaveModelToString(-1);	return 0; }
		int deserialize_from_string(string &str) {
			std::unique_ptr<Boosting> ret;
			string type = config_.boosting_type;
			if (type == std::string("gbdt")) ret.reset(new GBDT());
			else if (type == std::string("dart")) ret.reset(new DART());
			else if (type == std::string("goss")) ret.reset(new GOSS());
			else { fprintf(stderr, "deserialize MedLightGBM ERROR: unknown boosting type %s\n", type.c_str()); return -1; }
			if (!ret.get()->LoadModelFromString(str)) return -1;
			boosting_.reset(ret.release());
			return 0;
		}

		// n_preds 
		int n_preds_per_sample() const {

			int num_preb_in_one_row = config_.boosting_config.num_class;
			int is_pred_leaf = config_.io_config.is_predict_leaf_index ? 1 : 0;
			int num_iteration = config_.boosting_config.num_iterations;
			if (is_pred_leaf) {
				int max_iteration = num_iteration;
				if (num_iteration > 0) {
					num_preb_in_one_row *= static_cast<int>(std::min(max_iteration, num_iteration));
				}
				else {
					num_preb_in_one_row *= max_iteration;
				}
			}
			return num_preb_in_one_row;
		}

		void calc_feature_importance(vector<float> &features_importance_scores,
			const string &general_params, int max_feature_idx_);

	};

	class GBDT_Accessor : public GBDT {
	public:
		GBDT_Accessor(LightGBM::Boosting *booster) {
			string mdl = booster->SaveModelToString(-1);
			LoadModelFromString(mdl);
		}

		vector<float> FeatureImportanceTrick(const string &method = "frequency") {
			vector<pair<size_t, string>> res = FeatureImportance(method);

			vector<float> final_res(MaxFeatureIdx() + 1);
			for (size_t i = 0; i < res.size(); ++i) {
				int index = stoi(boost::replace_all_copy(res[i].second, "Column_", ""));
				if (index >= final_res.size() || index < 0)
					throw out_of_range("index is out of range: " + to_string(index) + " max=" +
						to_string(final_res.size()));
				final_res[index] = (float)res[i].first;
			}

			return final_res;
		}
	};

};



//=============================================
// MedLightGBM
//=============================================

struct MedLightGBMParams : public SerializableObject {

	MedLightGBMParams() {

		defaults = "";
		defaults += "boosting_type=gbdt;";
		defaults += "objective=binary;";
		defaults += "metric=binary_logloss,auc;";
		defaults += "metric_freq=1;";
		defaults += "is_training_metric=true;";
		defaults += "max_bin=255;";
		defaults += "num_trees=200;";
		defaults += "learning_rate=0.05;";
		defaults += "tree_learner=serial;";
		defaults += "num_threads=12;";
		defaults += "min_data_in_leaf=50;";
		defaults += "min_sum_hessian_in_leaf=5.0;";
		defaults += "is_enable_sparse=false;";
		defaults += "num_machines=1;";
		defaults += "feature_fraction=0.8;";
		defaults += "bagging_fraction=0.25;";
		defaults += "bagging_freq=4;";
		defaults += "is_unbalance=true;";
		defaults += "num_leaves=80"; // keep last param without a ; at the end 

	}

	string defaults = "";
	string user_params = "";

	ADD_SERIALIZATION_FUNCS(defaults, user_params);
};

class MedLightGBM : public MedPredictor {
public:
	MedLightGBMParams params;

	LightGBM::MemApp mem_app;

	/// please reffer to KeyAliasTransform in LightGBM::ParameterAlias
	int init(map<string, string>& initialization_map) { return mem_app.init(initialization_map); }

	int init_from_string(string init_str) {
		params.user_params += init_str;
		string init = params.defaults + ";" + params.user_params;
		//fprintf(stderr, "Calling MedLightGBM init with :\ninit_str %s\n user_params %s\n all %s\n", init_str.c_str(), params.user_params.c_str(), init.c_str());
		map<string, string> init_map;
		initialization_text_to_map(init, init_map);
		return MedLightGBM::init(init_map);
	}


	// Function
	MedLightGBM() {
		classifier_type = MODEL_LIGHTGBM;
		normalize_for_learn = false; //true;
		normalize_for_predict = false; //true;
		normalize_y_for_learn = false;
		transpose_for_learn = false;
		transpose_for_predict = false;
		init_from_string("");
		_mark_learn_done = false;
	};
	~MedLightGBM() {};

	// learn predict
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
		fprintf(stderr, "Starting a LightGBM train session...\n");
		mem_app.InitTrain(x, y, w, nsamples, nftrs);
		mem_app.Train();
		_mark_learn_done = true;
		return 0;
	}
	int Learn(float *x, float *y, int nsamples, int nftrs) { return Learn(x, y, NULL, nsamples, nftrs); }
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const {
		//mem_app.InitPredict(x, nsamples, nftrs);
		mem_app.Predict(x, nsamples, nftrs, preds);
		return 0;
	}

	//virtual void print(FILE *fp, const string& prefix);

	void calc_feature_importance(vector<float> &features_importance_scores,
		const string &general_params);


	int n_preds_per_sample() const { return mem_app.n_preds_per_sample(); }


	// serializations
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
private:
	bool _mark_learn_done;
};


MEDSERIALIZE_SUPPORT(MedLightGBMParams);
MEDSERIALIZE_SUPPORT(MedLightGBM);


#endif
