#include "PredictorImputer.h"
#include <algorithm>
#include <random>

using namespace std;

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

Predictor_Imputer_Params::Predictor_Imputer_Params() {
	predictor_type = "lightgbm";
	predictor_args = "objective=multiclass;metric=multi_logloss;verbose=0;num_threads=0;num_trees=80;"
		"learning_rate=0.05;lambda_l2=0;metric_freq=50;is_training_metric=false;max_bin=255;min_data_in_leaf=30"
		";feature_fraction=0.8;bagging_fraction=0.25;bagging_freq=4;is_unbalance=true;num_leaves=80;silent=2";
	num_class_setup = "num_class";
	calibration_string = "calibration_type=isotonic_regression;verbose=0";
	calibration_save_ratio = (float)0.2;

	bin_settings.init_from_string("split_method=iterative_merge;min_bin_count=100;binCnt=100");
}

int Predictor_Imputer_Params::init(map<string, string>& map) {
	for (auto &it : map)
	{
		if (it.first == "predictor_type")
			predictor_type = it.second;
		else if (it.first == "predictor_args")
			predictor_args = it.second;
		else if (it.first == "num_class_setup")
			num_class_setup = it.second;
		else if (it.first == "calibration_string")
			calibration_string = it.second;
		else if (it.first == "calibration_save_ratio")
			calibration_save_ratio = med_stof(it.second);
		else if (it.first == "bin_settings")
			bin_settings.init_from_string(it.second);
		else
			MTHROW_AND_ERR("Error in PredictorImputer::init - unsupported argument %s\n", it.first.c_str());
	}

	if (calibration_save_ratio >= 1)
		MTHROW_AND_ERR("Error Predictor_Imputer_Params::init - calibration_save_ratio should be in [0,1) range\n");
	return 0;
}

PredictorImputer::~PredictorImputer() {
	if (predictor != NULL)
		delete predictor;
	predictor = NULL;
}

void PredictorImputer::init_defaults() {
	processor_type = FTR_PROCESS_PREDICTOR_IMPUTER;
	missing_value = MED_MAT_MISSING_VALUE;
	predictor = NULL;
	feature_name = "";
	resolved_feature_name = "";
	verbose_learn = true;
	find_real_value = true;
	debug = false;

	random_device rd;
	gen = mt19937(rd());
}

int PredictorImputer::init(map<string, string>& mapper) {

	for (auto &it : mapper)
	{
		//! [PredictorImputer::init]
		if (it.first == "feature_name")
			feature_name = it.second;
		else if (it.first == "missing_value")
			missing_value = med_stof(it.second);
		else if (it.first == "params")
			params.init_from_string(it.second);
		else if (it.first == "verbose_learn")
			verbose_learn = med_stoi(it.second) > 0;
		else if (it.first == "find_real_value")
			find_real_value = med_stoi(it.second) > 0;
		else if (it.first == "debug")
			debug = med_stoi(it.second) > 0;
		else if (it.first == "fp_type" || it.first == "tag") {}
		else
			MTHROW_AND_ERR("Error in PredictorImputer::init - unsupported argument %s\n", it.first.c_str());
		//! [PredictorImputer::init]
	}
	if (params.predictor_type.empty())
		MTHROW_AND_ERR("Error in PredictorImputer::init - must provide params with predictor_type\n");

	return 0;
}

int PredictorImputer::Learn(MedFeatures& features, unordered_set<int>& ids) {
	if (feature_name.empty())
		MTHROW_AND_ERR("Error in PredictorImputer::Learn - must provide feature_name\n");
	resolved_feature_name = resolve_feature_name(features, feature_name);
	int feature_ind = -1;

	vector<string> all_names;
	features.get_feature_names(all_names);
	for (int i = 0; i < all_names.size() && feature_ind < 0; ++i)
		if (all_names[i] == resolved_feature_name)
			feature_ind = i;

	if (feature_ind < 0)
		MTHROW_AND_ERR("Error in PredictorImputer::Learn - can't find feature %s index\n",
			resolved_feature_name.c_str());

	const vector<float> &feat_vals_o = features.data.at(resolved_feature_name);
	//clear missing_value id's:
	vector<float> feat_vals;
	feat_vals.reserve(feat_vals_o.size());
	vector<int> non_miss_idx;
	for (size_t i = 0; i < feat_vals_o.size(); ++i)
		if (feat_vals_o[i] != missing_value) {
			non_miss_idx.push_back((int)i);
			feat_vals.push_back(feat_vals_o[i]);
		}

	int train_sz = (int)feat_vals.size();
	int pred_num_feats = (int)features.data.size() - 1; //use all other feature execept this one to predict value
	predictor_features.resize(pred_num_feats);
	for (size_t i = 0; i < predictor_features.size(); ++i)
	{
		int fix_i = (int)i + int(i >= feature_ind);
		predictor_features[i] = all_names[fix_i];
	}
	vector<float> train_vec(train_sz * pred_num_feats), label_vec(train_sz);

	unordered_set<float> uniq_vals(feat_vals.begin(), feat_vals.end());
	sorted_uniq_vals.clear();
	sorted_uniq_vals.insert(sorted_uniq_vals.end(), uniq_vals.begin(), uniq_vals.end());
	sort(sorted_uniq_vals.begin(), sorted_uniq_vals.end());

	vector<int> empt;
	medial::process::split_feature_to_bins(params.bin_settings, feat_vals, empt, feat_vals);
	unordered_set<float> seen_val(feat_vals.begin(), feat_vals.end());
	int class_num = (int)seen_val.size();
	num_classes = class_num;
	string predictor_init = params.predictor_args;
	//set num classes if needed:
	if (!params.num_class_setup.empty()) {
		//std::regex rgx(params.num_class_setup + "=[^;]+");
		//predictor_init = std::regex_replace(predictor_init, rgx, empty_str);
		//boost::replace_all(predictor_init, ";;", ";");
		predictor_init += ";" + params.num_class_setup + "=" + to_string(class_num);
		//change predictor_init
	}
	sorted_bin_vals.clear();
	sorted_bin_vals.insert(sorted_bin_vals.end(), seen_val.begin(), seen_val.end());
	sort(sorted_bin_vals.begin(), sorted_bin_vals.end());

	predictor = MedPredictor::make_predictor(params.predictor_type, predictor_init);
	if (verbose_learn)
		MLOG("Feature %s has %d categories\n", resolved_feature_name.c_str(), class_num);
	//change labels to be 0 to K-1:
	unordered_map<float, int> map_categ;
	//calc by order:
	for (size_t ii = 0; ii < sorted_bin_vals.size(); ++ii)
		map_categ[sorted_bin_vals[ii]] = (int)ii;
	//commit to create labels:
	for (size_t ii = 0; ii < feat_vals.size(); ++ii)
		label_vec[ii] = (float)map_categ.at(feat_vals[ii]);
	//prepate train:
	for (size_t ii = 0; ii < train_sz; ++ii) {
		int real_idx = non_miss_idx[ii];
		for (size_t jj = 0; jj < pred_num_feats; ++jj) {
			int fixed_idx = (int)jj + int(jj >= feature_ind); //skip current feature
			train_vec[ii* pred_num_feats + jj] = float(features.data.at(all_names[fixed_idx])[real_idx]);
		}
	}

	//learn:
	if (params.calibration_save_ratio > 0) {
		//use Calibrator
		calibrators.resize(class_num);

		for (size_t kk = 0; kk < calibrators.size(); ++kk)
			calibrators[kk].init_from_string(params.calibration_string);

		int calib_ratio = params.calibration_save_ratio * (int)label_vec.size();
		int train_ratio = (int)label_vec.size() - calib_ratio;
		uniform_int_distribution<> sel_rnd(0, (int)label_vec.size() - 1);
		vector<bool> seen_sel(label_vec.size());
		vector<float> pred_train_vec(train_ratio * pred_num_feats), pred_label_vec(train_ratio);
		vector<vector<MedSample>> pred_calib_train(class_num);
		MedFeatures pred_calib_mat;
		for (const string &name_feat : all_names)
			pred_calib_mat.attributes[name_feat].denorm_mean = 0;
		pred_calib_mat.samples.resize(calib_ratio);

		for (size_t j = 0; j < calib_ratio; ++j)
		{
			int sel_idx = sel_rnd(gen);
			while (seen_sel[sel_idx])
				sel_idx = sel_rnd(gen);
			seen_sel[sel_idx] = true;

			int categ = (int)label_vec[sel_idx];
			MedSample smp; smp.id = (int)j;

			for (size_t k = 0; k < pred_calib_train.size(); ++k)
			{
				smp.outcome = k == categ; // set outcome := 1 (as case) only for categ
				pred_calib_train[k].push_back(smp);
			}
			pred_calib_mat.samples[j].id = (int)j;
			pred_calib_mat.samples[j].outcome = 0; //doesn't matter for prediction only
			for (size_t k = 0; k < pred_num_feats; ++k) {
				int fixed_idx = (int)k + int(k >= feature_ind); //skip current
				pred_calib_mat.data[all_names[fixed_idx]].push_back(train_vec[sel_idx * pred_num_feats + k]);
			}
		}

		//build pred:
		int idx_train = 0;
		for (size_t j = 0; j < seen_sel.size(); ++j)
		{
			if (seen_sel[j])
				continue;
			for (size_t k = 0; k < pred_num_feats; ++k)
				pred_train_vec[idx_train * pred_num_feats + k] = train_vec[j * pred_num_feats + k];
			pred_label_vec[idx_train] = label_vec[j];
			++idx_train;
		}

		predictor->learn(pred_train_vec, pred_label_vec, train_ratio, pred_num_feats);
		//get predictions for pred_calib_train to learn calibrator:
		predictor->predict(pred_calib_mat);
		//get predictions into pred_calib_train:
		for (size_t j = 0; j < calibrators.size(); ++j)
			for (size_t k = 0; k < pred_calib_train[j].size(); ++k)
				pred_calib_train[j][k].prediction = { pred_calib_mat.samples[k].prediction[j] };

		//Learn calibrators - for each pred bin:
		for (size_t k = 0; k < calibrators.size(); ++k)
			calibrators[k].Learn(pred_calib_train[k]);
	}
	else
		predictor->learn(train_vec, label_vec, (int)label_vec.size(), pred_num_feats);

	return 0;
}

int PredictorImputer::_apply(MedFeatures& features, unordered_set<int>& ids) {
	resolved_feature_name = resolve_feature_name(features, feature_name);
	vector<float> &feat_vals = features.data.at(resolved_feature_name);

	vector<int> all_idx_to_impute;
	for (size_t i = 0; i < feat_vals.size(); ++i)
		if (feat_vals[i] == missing_value)
			all_idx_to_impute.push_back((int)i);
	if (all_idx_to_impute.empty())
		return 0;
	int model_feat_use = (int)predictor_features.size();
	vector<float> flat_x(all_idx_to_impute.size() * model_feat_use);
	vector<const vector<float> *> fast_access(predictor_features.size());
	//need same order of predcitor features
	for (size_t i = 0; i < fast_access.size(); ++i) {
		if (features.data.find(predictor_features[i]) == features.data.end())
			MTHROW_AND_ERR("Error in PredictorImputer::_apply - can't find input feature %s\n",
				predictor_features[i].c_str());
		fast_access[i] = &features.data.at(predictor_features[i]);
	}
	//fill:
	for (size_t i = 0; i < all_idx_to_impute.size(); ++i)
	{
		int idx = all_idx_to_impute[i];
		for (size_t j = 0; j < model_feat_use; ++j)
			flat_x[i* model_feat_use + j] = fast_access[j]->at(idx);
	}
	//full predict:
	vector<float> preds;
	predictor->predict(flat_x, preds, (int)all_idx_to_impute.size(), model_feat_use);
	if (preds.size() != all_idx_to_impute.size() * num_classes)
		MTHROW_AND_ERR("Error in PredictorImputer::_apply - after apply predictor, should contain %zu*%d preds and has %zu\n",
			all_idx_to_impute.size(), num_classes, preds.size());
	//do calibration for predictions if needed
	if (!calibrators.empty()) {
		vector<MedSample> cal_input(all_idx_to_impute.size());
		for (size_t i = 0; i < cal_input.size(); ++i)
		{
			cal_input[i].id = (int)i;
			cal_input[i].prediction.resize(1);
		}
		//apply for each calibrator class:
		for (size_t k = 0; k < num_classes; ++k)
		{
			//set predictions for class k
			for (size_t i = 0; i < cal_input.size(); ++i)
				cal_input[i].prediction[0] = preds[i * num_classes + k];
			//get results:
			calibrators[k].Apply(cal_input);
			//fetch again to preds:
			for (size_t i = 0; i < cal_input.size(); ++i)
				preds[i * num_classes + k] = cal_input[i].prediction[0];
		}
	}

	//sample randomliy by preds (which are calibrated)
	int print_ex = 0;
	int max_exm = 10;
	for (size_t i = 0; i < all_idx_to_impute.size(); ++i) {
		//for each needed to impute:
		double tot_num = 0;
		int curr_idx = (int)i * num_classes;
		for (size_t k = 0; k < num_classes; ++k)
			tot_num += preds[curr_idx + k];
		uniform_real_distribution<> real_dist(0, tot_num);
		double sel = real_dist(gen);

		//now select correspond bin value:
		double tot_num2 = 0;
		int sel_idx = 0;
		while (sel_idx < num_classes && tot_num2 + preds[curr_idx + sel_idx] < sel) {
			tot_num2 += preds[curr_idx + sel_idx];
			++sel_idx;
		}
		if (print_ex < max_exm) {
			++print_ex;
			stringstream input;
			bool first = true;
			for (const string &feat : predictor_features)
			{
				if (!first) {
					input << ", ";
					first = false;
				}
				input << feat << "=" << features.data.at(feat)[all_idx_to_impute[i]];
			}
			MLOG("DEBUG: feature %s :: [%s] => total_sum=%f, got = %f(%d), range= [%f, %f], parsed= %f\n",
				feature_name.c_str(), input.str().c_str(), tot_num, sel, sel_idx, sorted_bin_vals[0], sorted_bin_vals.back(),
				sorted_bin_vals[sel_idx]);
		}

		//now sel_idx is the index of the category:
		if (sel_idx >= sorted_bin_vals.size())
			sel_idx = (int)sorted_bin_vals.size() - 1;
		float impute_val = sorted_bin_vals[sel_idx];
		//if need to round to most close value in sorted_uniq_vals:
		if (find_real_value) {
			int pos = medial::process::binary_search_position(sorted_uniq_vals.data(),
				sorted_uniq_vals.data() + sorted_uniq_vals.size() - 1, impute_val);
			if (pos == 0)
				impute_val = sorted_uniq_vals[0];
			else {
				if (pos >= sorted_uniq_vals.size())
					impute_val = sorted_uniq_vals.back();
				else {
					float diff_next = abs(impute_val - sorted_uniq_vals[pos]);
					float diff_prev = abs(impute_val - sorted_uniq_vals[pos - 1]);
					if (diff_prev < diff_next)
						impute_val = sorted_uniq_vals[pos - 1];
					else
						impute_val = sorted_uniq_vals[pos];
				}
			}

		}

		//store value in matrix:
		feat_vals[all_idx_to_impute[i]] = impute_val;
	}
	return 0;
}