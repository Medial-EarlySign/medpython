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
#include <LightGBM/LightGBM/src/application/predictor.hpp>

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
	int MemApp::init(map<string, string>& init_params)
	{
		unordered_map<string, string> params;
		for (auto &e : init_params) params[e.first] = e.second;
		ParameterAlias::KeyAliasTransform(&params);
		// load configs
		config_.Set(params);
		Log::Info("Finished loading parameters");
		return 0;
	}

	//-------------------------------------------------------------------------------------------------
	int MemApp::InitTrainData(float *xdata, float *ydata, float *weight, int nrows, int ncols)
	{
		MLOG("MedLightGBM:: init train data %d x %d\n", nrows, ncols);
		if (config_.num_threads > 0) omp_set_num_threads(config_.num_threads);

		std::unique_ptr<Dataset> ret;
		auto get_row_fun = RowFunctionFromDenseMatric(xdata, nrows, ncols, C_API_DTYPE_FLOAT32, 1);

		// sample data first
		Random rand(config_.io_config.data_random_seed);
		int sample_cnt = static_cast<int>(nrows < config_.io_config.bin_construct_sample_cnt ? nrows : config_.io_config.bin_construct_sample_cnt);
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
		DatasetLoader loader(config_.io_config, nullptr, 1, nullptr);
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
		MLOG("training eval bit %d\n", config_.boosting_config.is_provide_training_metric);
		if (config_.boosting_config.is_provide_training_metric) {
			MLOG("Creating training metrics: types %d\n", config_.metric_types.size());
			for (auto metric_type : config_.metric_types) {
				auto metric = std::unique_ptr<Metric>(Metric::CreateMetric(metric_type, config_.metric_config));
				if (metric == nullptr) { continue; }
				metric->Init(train_data_->metadata(), train_data_->num_data());
				train_metric_.push_back(std::move(metric));
			}
		}
		train_metric_.shrink_to_fit();

		MLOG("MedLightGBM:: finished loading train mat\n");
		return 0;
	}

	//-------------------------------------------------------------------------------------------------
	int MemApp::InitTrain(float *xdata, float *ydata, float *weight, int nrows, int ncols) {
		if (config_.is_parallel) { Log::Info("parallel mode not supported yet for MedLightGBM !!"); return -1; }

		// create boosting
		boosting_.reset(Boosting::CreateBoosting(config_.boosting_type, config_.io_config.input_model.c_str()));

		// create objective function
		objective_fun_.reset(ObjectiveFunction::CreateObjectiveFunction(config_.objective_type, config_.objective_config));

		// load training data
		InitTrainData(xdata, ydata, weight, nrows, ncols);

		// initialize the objective function
		objective_fun_->Init(train_data_->metadata(), train_data_->num_data());

		// initialize the boosting
		boosting_->Init(&config_.boosting_config, train_data_.get(), objective_fun_.get(), Common::ConstPtrInVectorWrapper<Metric>(train_metric_));

		// add validation data into boosting ==> Currently not used, as we do not allow loading validation data at this stage in MedLightGBM
		//for (size_t i = 0; i < valid_datas_.size(); ++i)
		//	boosting_->AddValidDataset(valid_datas_[i].get(), Common::ConstPtrInVectorWrapper<Metric>(valid_metrics_[i]));

		Log::Info("Finished initializing training");
		return 0;
	}

	//-------------------------------------------------------------------------------------------------
	void MemApp::Train() {
		Log::Info("Started training...");
		int total_iter = config_.boosting_config.num_iterations;
		bool is_finished = false;
		bool need_eval = true;
		auto start_time = std::chrono::steady_clock::now();
		Log::Info("total_iter %d is_finished %d need_eval %d\n", total_iter, (int)is_finished, (int)need_eval);
		for (int iter = 0; iter < total_iter && !is_finished; ++iter) {
			is_finished = boosting_->TrainOneIter(nullptr, nullptr, need_eval);
			auto end_time = std::chrono::steady_clock::now();
			// output used time per iteration
			if ((((iter + 1) % config_.boosting_config.output_freq) == 0) || (iter == total_iter - 1))
				Log::Info("%f seconds elapsed, finished iteration %d", std::chrono::duration<double, std::milli>(end_time - start_time) * 1e-3, iter + 1);
		}
		Log::Info("Finished training");
	}


//-------------------------------------------------------------------------------------------------
void MemApp::Predict(float *x, int nrows, int ncols, float *&preds)
{
	auto get_row_fun = RowPairFunctionFromDenseMatric(x, nrows, ncols, C_API_DTYPE_FLOAT32, 1);

// create boosting
	Predictor predictor(boosting_.get(), config_.io_config.num_iteration_predict, config_.io_config.is_predict_raw_score, config_.io_config.is_predict_leaf_index,
		config_.io_config.pred_early_stop, config_.io_config.pred_early_stop_freq, config_.io_config.pred_early_stop_margin);
	int64_t num_preb_in_one_row = boosting_->NumPredictOneRow(config_.io_config.num_iteration_predict, config_.io_config.is_predict_leaf_index);
	auto pred_fun = predictor.GetPredictFunction();

	//string str;
	//serialize_to_string(str);

	int64_t len_res = nrows * num_preb_in_one_row;
	//MLOG("[MedLightGBM] predict: nrows %d , num_pred %d , len_res %d\n", nrows, num_preb_in_one_row, len_res);
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
		auto pred_wrt_ptr = out_result + static_cast<size_t>(num_preb_in_one_row) * i;
		pred_fun(one_row, pred_wrt_ptr);
		OMP_LOOP_EX_END();
	}
	OMP_THROW_EX();

	for (int64_t i=0; i<len_res; i++) preds[i] = (float)out_result[i];
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



}

//===============================================================================================
// MedLightGBM
//===============================================================================================

// serializations

//------------------------------------------------------------------------------------------------------------------------------------------
size_t MedLightGBM::get_size()
{
	size_t size = 0;
	size += MedSerialize::get_size(params);
	string str;
	if (mem_app.serialize_to_string(str) < 0)
		MERR("MedLightGBM::get_size() failed moving model to string\n");
	size += MedSerialize::get_size(str);

	size += MedSerialize::get_size(model_features);
	size += MedSerialize::get_size(features_count);
	return size;
}
//------------------------------------------------------------------------------------------------------------------------------------------
size_t MedLightGBM::serialize(unsigned char *blob)
{
	size_t size = 0;

	size += MedSerialize::serialize(blob, params);
	string str;
	if (mem_app.serialize_to_string(str) < 0)
		MERR("MedLightGBM::serialize() failed moving model to string\n");
	size += MedSerialize::serialize(blob + size, str);
	size += MedSerialize::serialize(blob + size, model_features);
	size += MedSerialize::serialize(blob + size, features_count);
	return size;
}
//------------------------------------------------------------------------------------------------------------------------------------------
size_t MedLightGBM::deserialize(unsigned char *blob)
{
	size_t size = MedSerialize::deserialize(blob, params);
	init_from_string(""); //loading the params as they were saved
	string str;
	size += MedSerialize::deserialize(blob + size, str);
	if (mem_app.deserialize_from_string(str) < 0)
		MERR("MedLightGBM::deserialize() failed moving model to string\n");
	size += MedSerialize::deserialize(blob + size, model_features);
	size += MedSerialize::deserialize(blob + size, features_count);
	return size;
}

