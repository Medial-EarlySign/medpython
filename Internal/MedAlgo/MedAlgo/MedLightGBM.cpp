#include "MedLightGBM.h"

#include <LightGBM/application.h>

#include <LightGBM/utils/common.h>
#include <LightGBM/utils/text_reader.h>

#include <LightGBM/network.h>
#include <LightGBM/dataset.h>
#include <LightGBM/dataset_loader.h>
#include <LightGBM/boosting.h>
#include <LightGBM/objective_function.h>
#include <LightGBM/prediction_early_stop.h>
#include <LightGBM/metric.h>
#include <LightGBM/c_api.h>

//#include "predictor.hpp"

#include <LightGBM/utils/openmp_wrapper.h>
#include <LightGBM/../../src/application/predictor.hpp>

#include <cstdio>
#include <ctime>

#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;



namespace LightGBM {

	std::function<std::vector<double>(int row_idx)>	RowFunctionFromDenseMatric(const void* data, int num_row, int num_col, int data_type, int is_row_major);
	std::function<std::vector<std::pair<int, double>>(int row_idx)> RowPairFunctionFromDenseMatric(const void* data, int num_row, int num_col, int data_type, int is_row_major);
	//-------------------------------------------------------------------------------------------------
	int MemApp::set_params(map<string, string>& init_params)
	{
		bool prev_silent = is_silent;
		is_silent = false;
		only_fatal = false;
		is_silent = init_params.find("verbose") != init_params.end()
			&& stoi(init_params.at("verbose")) <= 0;
		is_silent |= (init_params.find("verbosity") != init_params.end()
			&& stoi(init_params.at("verbosity")) <= 0);
		is_silent |= (init_params.find("silent") != init_params.end()
			&& stoi(init_params.at("silent")) > 0);
		if (global_logger.levels[LOG_MEDALGO] > LOG_DEF_LEVEL || is_silent)
			Log::ResetLogLevel(LogLevel::Warning);
		if (init_params.find("silent") != init_params.end() && stoi(init_params.at("silent")) > 1) {
			Log::ResetLogLevel(LogLevel::Fatal);
			only_fatal = true;
		}

		unordered_map<string, string> params;
		for (auto &e : init_params) params[e.first] = e.second;
		ParameterAlias::KeyAliasTransform(&params);
		// load configs
		config_.Set(params);

		if (!prev_silent)
			Log::Info("Finished loading parameters");
		return 0;
	}

	//-------------------------------------------------------------------------------------------------
	int MemApp::InitTrainData(float *xdata, float *ydata, const float *weight, int nrows, int ncols)
	{
		if (global_logger.levels[LOG_MEDALGO] > LOG_DEF_LEVEL || is_silent)
			Log::ResetLogLevel(LogLevel::Warning);
		if (only_fatal)
			Log::ResetLogLevel(LogLevel::Fatal);
		Log::Info("init train data %d x %d", nrows, ncols);
		if (config_.num_threads > 0) omp_set_num_threads(config_.num_threads);

		std::unique_ptr<Dataset> ret;
		auto get_row_fun = RowFunctionFromDenseMatric(xdata, nrows, ncols, C_API_DTYPE_FLOAT32, 1);

		// sample data first
		Random rand(config_.seed);
		int sample_cnt = static_cast<int>(nrows < config_.bin_construct_sample_cnt ? nrows : config_.bin_construct_sample_cnt);
		auto sample_indices = rand.Sample(nrows, sample_cnt);
		sample_cnt = static_cast<int>(sample_indices.size());
		std::vector<std::vector<double>> sample_values(ncols);
		std::vector<std::vector<int>> sample_idx(ncols);
		for (size_t i = 0; i < sample_indices.size(); ++i) {
			auto idx = sample_indices[i];
			auto row = get_row_fun(static_cast<int>(idx));
			for (size_t j = 0; j < row.size(); ++j) {
				if (std::fabs(row[j]) > kEpsilon) {
					sample_values[j].emplace_back(row[j]);
					sample_idx[j].emplace_back(static_cast<int>(i));
				}
			}
		}
		DatasetLoader loader(config_, nullptr, 1, nullptr);
		train_data_.reset(loader.CostructFromSampleData(Common::Vector2Ptr<double>(sample_values).data(),
			Common::Vector2Ptr<int>(sample_idx).data(),
			static_cast<int>(sample_values.size()),
			Common::VectorSize<double>(sample_values).data(),
			sample_cnt, nrows));

		OMP_INIT_EX();
#pragma omp parallel for schedule(static)
		for (int i = 0; i < nrows; ++i) {
			OMP_LOOP_EX_BEGIN();
			const int tid = omp_get_thread_num();
			auto one_row = get_row_fun(i);
			train_data_->PushOneRow(tid, i, one_row);
			OMP_LOOP_EX_END();
		}
		OMP_THROW_EX();
		train_data_->FinishLoad();

		// load label
		train_data_->SetFloatField("label", ydata, nrows);

		// load weight
		if (weight != NULL)
			train_data_->SetFloatField("weight", weight, nrows);

		// create training metric
		Log::Info("training eval bit %d", config_.is_provide_training_metric);
		if (config_.is_provide_training_metric) {
			Log::Info("Creating training metrics: types %d [%s]", 
				config_.metric.size(), medial::io::get_list(config_.metric, ", ").c_str());
			for (auto metric_type : config_.metric) {
				auto metric = std::unique_ptr<Metric>(Metric::CreateMetric(metric_type, config_));
				if (metric == nullptr) { continue; }
				metric->Init(train_data_->metadata(), train_data_->num_data());
				train_metric_.push_back(std::move(metric));
			}
		}
		train_metric_.shrink_to_fit();

		Log::Info("finished loading train mat");
		return 0;
	}

	//-------------------------------------------------------------------------------------------------
	int MemApp::InitTrain(float *xdata, float *ydata, const float *weight, int nrows, int ncols) {
		if (global_logger.levels[LOG_MEDALGO] > LOG_DEF_LEVEL || is_silent)
			Log::ResetLogLevel(LogLevel::Warning);
		if (only_fatal)
			Log::ResetLogLevel(LogLevel::Fatal);
		if (config_.is_parallel) { Log::Info("parallel mode not supported yet for MedLightGBM !!"); return -1; }

		// create boosting
		boosting_.reset(Boosting::CreateBoosting(config_.boosting, config_.input_model.c_str()));

		// create objective function
		objective_fun_.reset(ObjectiveFunction::CreateObjectiveFunction(config_.objective, config_));

		// load training data
		InitTrainData(xdata, ydata, weight, nrows, ncols);

		// initialize the objective function
		objective_fun_->Init(train_data_->metadata(), train_data_->num_data());

		// initialize the boosting
		boosting_->Init(&config_, train_data_.get(), objective_fun_.get(), Common::ConstPtrInVectorWrapper<Metric>(train_metric_));

		// add validation data into boosting ==> Currently not used, as we do not allow loading validation data at this stage in MedLightGBM
		//for (size_t i = 0; i < valid_datas_.size(); ++i)
		//	boosting_->AddValidDataset(valid_datas_[i].get(), Common::ConstPtrInVectorWrapper<Metric>(valid_metrics_[i]));

		Log::Info("Finished initializing training");
		return 0;
	}

	//-------------------------------------------------------------------------------------------------
	void MemApp::Train() {
		if (global_logger.levels[LOG_MEDALGO] > LOG_DEF_LEVEL || is_silent)
			Log::ResetLogLevel(LogLevel::Warning);
		if (only_fatal)
			Log::ResetLogLevel(LogLevel::Fatal);
		Log::Info("Started training...");
		int total_iter = config_.num_iterations;
		bool is_finished = false;
		bool need_eval = true;
		auto start_time = std::chrono::steady_clock::now();
		Log::Info("total_iter %d is_finished %d need_eval %d", total_iter, (int)is_finished, (int)need_eval);
		for (int iter = 0; iter < total_iter && !is_finished; ++iter) {
			is_finished = boosting_->TrainOneIter(nullptr, nullptr);
			auto end_time = std::chrono::steady_clock::now();
			// output used time per iteration
			if ((((iter + 1) % config_.metric_freq) == 0) || (iter == total_iter - 1)) {

				if (!config_.metric.empty()) {
					vector<double> m_res = boosting_->GetEvalAt(0);
					stringstream eval_str;


					eval_str << config_.metric[0] << "=";
					if (!m_res.empty())
						eval_str << m_res[0];
					for (size_t i = 1; i < config_.metric.size(); ++i) {
						eval_str << ", " << config_.metric[i] << "=";
						if (i < m_res.size())
							eval_str << m_res[i];
					}

					Log::Info("%f seconds elapsed, finished iteration %d. [%s]",
						std::chrono::duration<double, std::milli>(end_time - start_time) * 1e-3, iter + 1,
						eval_str.str().c_str());
				}
				else
					Log::Info("%f seconds elapsed, finished iteration %d", std::chrono::duration<double, std::milli>(end_time - start_time) * 1e-3, iter + 1);
			}
		}
		Log::Info("Finished training");
	}

	void MemApp::fetch_boosting(LightGBM::Boosting *&res) {
		res = boosting_.get();
	}

	void MemApp::fetch_early_stop(LightGBM::PredictionEarlyStopInstance &early_stop_) {
		LightGBM::Boosting *_boosting = boosting_.get();

		early_stop_ = CreatePredictionEarlyStopInstance("none", LightGBM::PredictionEarlyStopConfig());
		if (config_.pred_early_stop && !_boosting->NeedAccuratePrediction()) {
			LightGBM::PredictionEarlyStopConfig pred_early_stop_config;
			pred_early_stop_config.margin_threshold = config_.pred_early_stop_margin;
			pred_early_stop_config.round_period = config_.pred_early_stop_freq;
			if (_boosting->NumberOfClasses() == 1) {
				early_stop_ = CreatePredictionEarlyStopInstance("binary", pred_early_stop_config);
			}
			else {
				early_stop_ = CreatePredictionEarlyStopInstance("multiclass", pred_early_stop_config);
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	void MemApp::Predict(float *x, int nrows, int ncols, float *&preds) const
	{
		auto get_row_fun = RowPairFunctionFromDenseMatric(x, nrows, ncols, C_API_DTYPE_FLOAT32, 1);

		// create boosting
		Predictor predictor(boosting_.get(), config_.num_iteration_predict, config_.predict_raw_score, config_.predict_leaf_index, config_.predict_contrib,
			config_.pred_early_stop, config_.pred_early_stop_freq, config_.pred_early_stop_margin);
		int64_t num_pred_in_one_row = boosting_->NumPredictOneRow(config_.num_iteration_predict, config_.predict_leaf_index, config_.predict_contrib);
		auto pred_fun = predictor.GetPredictFunction();

		//string str;
		//serialize_to_string(str);

		int64_t len_res = nrows * num_pred_in_one_row;
		//MLOG("[MedLightGBM] predict: nrows %d , num_pred %d , len_res %d\n", nrows, num_pred_in_one_row, len_res);
		//MLOG("[MedLightGBM] predict: num_iter %d , is_raw %d , is_leaf %d\n", 
		//	config_.io_config.num_iteration_predict, config_.io_config.is_predict_raw_score ? 1 : 0, config_.io_config.is_predict_leaf_index ? 1:0);
		vector<double> out_result_vec(len_res);
		double *out_result = &out_result_vec[0];
		if (preds == NULL) preds = new float[len_res];

		OMP_INIT_EX();
#pragma omp parallel for schedule(static)
		for (int i = 0; i < nrows; ++i) {
			OMP_LOOP_EX_BEGIN();
			auto one_row = get_row_fun(i);
			auto pred_wrt_ptr = out_result + static_cast<size_t>(num_pred_in_one_row) * i;
			pred_fun(one_row, pred_wrt_ptr);
			OMP_LOOP_EX_END();
		}
		OMP_THROW_EX();

		for (int64_t i = 0; i < len_res; i++) preds[i] = (float)out_result[i];
	}

	//-------------------------------------------------------------------------------------------------
	void MemApp::PredictShap(float *x, int nrows, int ncols, float *&shap_vals) const
	{
		auto get_row_fun = RowPairFunctionFromDenseMatric(x, nrows, ncols, C_API_DTYPE_FLOAT32, 1);

		// create boosting
		Predictor predictor(boosting_.get(), config_.num_iteration_predict, config_.predict_raw_score, config_.predict_leaf_index, true,
			config_.pred_early_stop, config_.pred_early_stop_freq, config_.pred_early_stop_margin);
		int64_t num_pred_in_one_row = boosting_->NumPredictOneRow(config_.num_iteration_predict, config_.predict_leaf_index, true);
		auto pred_fun = predictor.GetPredictFunction();

		int64_t len_res = nrows * num_pred_in_one_row;

		vector<double> out_result_vec(len_res);
		double *out_result = &out_result_vec[0];
		if (shap_vals == NULL)
			MTHROW_AND_ERR("Error MemApp::PredictShap - shap_vals was NULL\n");

		OMP_INIT_EX();
#pragma omp parallel for schedule(static)
		for (int i = 0; i < nrows; ++i) {
			OMP_LOOP_EX_BEGIN();
			auto one_row = get_row_fun(i);
			auto pred_wrt_ptr = out_result + static_cast<size_t>(num_pred_in_one_row) * i;
			pred_fun(one_row, pred_wrt_ptr);
			for (int j = 0; j < num_pred_in_one_row; j++)
				shap_vals[i*num_pred_in_one_row + j] = (float)pred_wrt_ptr[j];
			OMP_LOOP_EX_END();
		}
		OMP_THROW_EX();

	}

	//-----------------------------------------------------------------------------------------------------------
	//----- start of some help functions
	//-----------------------------------------------------------------------------------------------------------
	std::function<std::vector<std::pair<int, double>>(int row_idx)>
		RowPairFunctionFromDenseMatric(const void* data, int num_row, int num_col, int data_type, int is_row_major) {
		auto inner_function = RowFunctionFromDenseMatric(data, num_row, num_col, data_type, is_row_major);
		if (inner_function != nullptr) {
			return [inner_function](int row_idx) {
				auto raw_values = inner_function(row_idx);
				std::vector<std::pair<int, double>> ret;
				for (int i = 0; i < static_cast<int>(raw_values.size()); ++i) {
					if (std::fabs(raw_values[i]) > 1e-15) {
						ret.emplace_back(i, raw_values[i]);
					}
				}
				return ret;
			};
		}
		return nullptr;
	}

	std::function<std::vector<double>(int row_idx)>
		RowFunctionFromDenseMatric(const void* data, int num_row, int num_col, int data_type, int is_row_major) {
		if (data_type == C_API_DTYPE_FLOAT32) {
			const float* data_ptr = reinterpret_cast<const float*>(data);
			if (is_row_major) {
				return [data_ptr, num_col, num_row](int row_idx) {
					std::vector<double> ret(num_col);
					auto tmp_ptr = data_ptr + static_cast<size_t>(num_col) * row_idx;
					for (int i = 0; i < num_col; ++i) {
						ret[i] = static_cast<double>(*(tmp_ptr + i));
						if (std::isnan(ret[i])) {
							ret[i] = 0.0f;
						}
					}
					return ret;
				};
			}
			else {
				return [data_ptr, num_col, num_row](int row_idx) {
					std::vector<double> ret(num_col);
					for (int i = 0; i < num_col; ++i) {
						ret[i] = static_cast<double>(*(data_ptr + static_cast<size_t>(num_row) * i + row_idx));
						if (std::isnan(ret[i])) {
							ret[i] = 0.0f;
						}
					}
					return ret;
				};
			}
		}
		else if (data_type == C_API_DTYPE_FLOAT64) {
			const double* data_ptr = reinterpret_cast<const double*>(data);
			if (is_row_major) {
				return [data_ptr, num_col, num_row](int row_idx) {
					std::vector<double> ret(num_col);
					auto tmp_ptr = data_ptr + static_cast<size_t>(num_col) * row_idx;
					for (int i = 0; i < num_col; ++i) {
						ret[i] = static_cast<double>(*(tmp_ptr + i));
						if (std::isnan(ret[i])) {
							ret[i] = 0.0f;
						}
					}
					return ret;
				};
			}
			else {
				return [data_ptr, num_col, num_row](int row_idx) {
					std::vector<double> ret(num_col);
					for (int i = 0; i < num_col; ++i) {
						ret[i] = static_cast<double>(*(data_ptr + static_cast<size_t>(num_row) * i + row_idx));
						if (std::isnan(ret[i])) {
							ret[i] = 0.0f;
						}
					}
					return ret;
				};
			}
		}
		throw std::runtime_error("unknown data type in RowFunctionFromDenseMatric");
	}

	void  MemApp::calc_feature_importance(vector<float> &features_importance_scores,
		const string &general_params, int max_feature_idx_) {

		map<string, string> params;
		MedSerialize::init_map_from_string(general_params, params);
		string importance_type = "gain"; //"frequency"; //"gain";
		if (params.find("importance_type") != params.end())
			importance_type = params.at("importance_type");

		GBDT_Accessor booster_access(boosting_.get());
		features_importance_scores = booster_access.FeatureImportanceTrick(importance_type);
	}

}

//===============================================================================================
// MedLightGBM
//===============================================================================================

void MedLightGBM::calc_feature_importance(vector<float> &features_importance_scores,
	const string &general_params, const MedFeatures *features) {
	if (!_mark_learn_done)
		MTHROW_AND_ERR("ERROR:: Requested calc_feature_importance before running learn\n");

	map<string, string> params;

	unordered_set<string> local_types = { "gain", "split" };
	unordered_set<string> legal_types = { "gain", "split", "shap" };

	MedSerialize::initialization_text_to_map(general_params, params);
	string importance_type = "gain"; // default
	for (auto it = params.begin(); it != params.end(); ++it)
		if (it->first == "importance_type")
			importance_type = it->second;
		else
			MTHROW_AND_ERR("Unsupported calc_feature_importance param \"%s\"\n", it->first.c_str());

	if (legal_types.find(importance_type) == legal_types.end())
		MTHROW_AND_ERR("Ilegal importance_type value \"%s\" "
			"- should by one of [%s]\n",
			importance_type.c_str(), medial::io::get_list(legal_types, ", ").c_str());

	features_importance_scores.resize(model_features.empty() ? features_count : (int)model_features.size());

	if (local_types.count(importance_type) > 0) {
		mem_app.calc_feature_importance(features_importance_scores, "importance_type=" + importance_type,
			(model_features.empty() ? features_count : (int)model_features.size()));
		return;
	}

	// shap option
	MedMat<float> feat_mat, contribs_mat;
	if (features == NULL)
		MTHROW_AND_ERR("SHAP values feature importance requires features \n");

	features->get_as_matrix(feat_mat);
	calc_feature_contribs(feat_mat, contribs_mat);
#pragma omp parallel for
	for (int j = 0; j < contribs_mat.ncols; ++j)
	{
		float col_sum = 0;

		for (int i = 0; i < contribs_mat.nrows; ++i)
		{
			col_sum += abs(contribs_mat.get(i, j));
		}
		features_importance_scores[j] = col_sum / (float)contribs_mat.nrows;
	}

}

void MedLightGBM::calc_feature_contribs(MedMat<float> &x, MedMat<float> &contribs)
{
	int nrows = x.nrows;
	int ncols = x.ncols;

	contribs.resize(nrows, ncols + 1);
	// copy metadata
	contribs.signals.insert(contribs.signals.end(), x.signals.begin(), x.signals.end());
	contribs.signals.push_back("b0");
	contribs.recordsMetadata.insert(contribs.recordsMetadata.end(), x.recordsMetadata.begin(), x.recordsMetadata.end());

	float *contribs_ptr = contribs.data_ptr();
	float *x_ptr = x.data_ptr();
	mem_app.PredictShap(x_ptr, nrows, ncols, contribs_ptr);
}

void MedLightGBM::prepare_predict_single() {
	if (!prepared_single) {
		num_preds = n_preds_per_sample();

		mem_app.fetch_boosting(_boosting);
		mem_app.fetch_early_stop(early_stop_);
		prepared_single = true;
	}
}

void MedLightGBM::predict_single(const vector<float> &x, vector<float> &preds) const {
	int n_ftrs = (int)x.size();
	vector<double> one_row(n_ftrs);
	for (int i = 0; i < n_ftrs; ++i)
		one_row[i] = static_cast<double>(x[i]);

	vector<double> out_result_vec(num_preds);
	predict_single(one_row, out_result_vec);

	preds.resize(num_preds);
	for (int64_t i = 0; i < num_preds; i++) preds[i] = (float)out_result_vec[i];
}

void MedLightGBM::predict_single(const vector<double> &x, vector<double> &preds) const {
	preds.resize(num_preds);
	double *out_result = &preds[0];
	_boosting->Predict(x.data(), out_result, &early_stop_);
}

void MedLightGBM::export_predictor(const string &output_fname)
{
	string predictor_str;
	mem_app.serialize_to_string(predictor_str);
	ofstream ofs(output_fname, std::ios::binary);
	if (!ofs)
	{
		MTHROW_AND_ERR("MedLightGBM::export_predictor failed. couldn't write %s \n", output_fname.c_str());
	}
	ofs << predictor_str;
	ofs.close();
}


void MedLightGBM::print(FILE *fp, const string& prefix, int level) const {

	if (level == 0)
		fprintf(fp, "%s: MedLightGBM ()\n", prefix.c_str());
	else {
		string predictor_str;
		mem_app.serialize_to_string(predictor_str);
		fprintf(fp, "%s: MedLightGBM ()\n%s\n", prefix.c_str(), predictor_str.c_str());
	}

}