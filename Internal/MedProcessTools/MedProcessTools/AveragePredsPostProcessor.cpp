#include "AveragePredsPostProcessor.h"
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL LOG_DEF_LEVEL

//defaults:
AveragePredsPostProcessor::AveragePredsPostProcessor() {
	feature_processor_type = "";
	feature_processor_args = "";
	use_median = false;
	resample_cnt = 50;
	batch_size = 10000;
	force_cancel_imputations = false;
}

int AveragePredsPostProcessor::init(map<string, string> &mapper) {
	for (const auto &it : mapper)
	{
		if (it.first == "feature_processor_type")
			feature_processor_type = it.second;
		else if (it.first == "feature_processor_args")
			feature_processor_args = it.second;
		else if (it.first == "use_median")
			use_median = med_stoi(it.second) > 0;
		else if (it.first == "resample_cnt")
			resample_cnt = med_stoi(it.second);
		else if (it.first == "batch_size")
			batch_size = med_stoi(it.second);
		else if (it.first == "force_cancel_imputations")
			force_cancel_imputations = med_stoi(it.second) > 0;
		else if (it.first == "pp_type") {} //ignore
		else
			MTHROW_AND_ERR("Error AveragePredsPostProcessor::init - unknown argument %s\n",
				it.first.c_str());
	}

	if (feature_processor_type.empty())
		MTHROW_AND_ERR("Error AveragePredsPostProcessor::init - must provide feature_processor_type\n");
	if (resample_cnt <= 0)
		MTHROW_AND_ERR("Error AveragePredsPostProcessor::init - resample_cnt > 0\n");
	if (batch_size <= 0)
		MTHROW_AND_ERR("Error AveragePredsPostProcessor::init - batch_size > 0\n");

	if (!boost::starts_with(feature_processor_type, "MODEL::")) {
		feature_processor = FeatureProcessor::make_processor(feature_processor_type, feature_processor_args);
	}

	return 0;
}

void AveragePredsPostProcessor::init_post_processor(MedModel& model)
{
	model_predictor = model.predictor;
	if (boost::starts_with(feature_processor_type, "MODEL::")) {
		feature_processor = NULL;
		//store feature processor if has or need from model.
		string feature_type = feature_processor_type.substr(7);
		vector<FeatureProcessor *> fps_flat;
		for (size_t i = 0; i < model.feature_processors.size(); ++i)
		{
			if (feature_processor->processor_type == FTR_PROCESS_MULTI) {
				vector<FeatureProcessor *> &to_add = static_cast<MultiFeatureProcessor *>(model.feature_processors[i])->processors;
				fps_flat.insert(fps_flat.end(), to_add.begin(), to_add.end());
			}
			else
				fps_flat.push_back(model.feature_processors[i]);
		}

		int has_before = 0;
		for (FeatureProcessor *f : fps_flat)
		{
			if (f->my_class_name() == feature_type) {
				if (feature_processor != NULL)
					MTHROW_AND_ERR("Error AveragePredsPostProcessor::init_post_processor - found multiple feature processors of type %s\n",
						feature_type.c_str());
				feature_processor = f;
			}
			else {
				if (feature_processor != NULL)
					after_processors.push_back(f);
				else
					++has_before;
			}
		}
		if (feature_processor == NULL)
			MTHROW_AND_ERR("Error AveragePredsPostProcessor::init_post_processor - can't find feature processors of type %s\n",
				feature_type.c_str());
		if (has_before > 0)
			MWARN("WARN:: AveragePredsPostProcessor :: found %d processors before\n", has_before);
		if (!after_processors.empty())
			MLOG("INFO:: AveragePredsPostProcessor :: found %zu processors after\n", after_processors.size());
	}
	else if (force_cancel_imputations) {
		//find imputer and store all what happens after him:
		vector<bool> processors_tp(FTR_PROCESS_LAST);
		processors_tp[FTR_PROCESS_IMPUTER] = true;
		processors_tp[FTR_PROCESS_ITERATIVE_IMPUTER] = true;
		processors_tp[FTR_PROCESS_PREDICTOR_IMPUTER] = true;

		vector<FeatureProcessor *> fps_flat;
		for (size_t i = 0; i < model.feature_processors.size(); ++i)
		{
			if (feature_processor->processor_type == FTR_PROCESS_MULTI) {
				vector<FeatureProcessor *> &to_add = static_cast<MultiFeatureProcessor *>(model.feature_processors[i])->processors;
				fps_flat.insert(fps_flat.end(), to_add.begin(), to_add.end());
			}
			else
				fps_flat.push_back(model.feature_processors[i]);
		}

		bool found = false;
		int has_before = 0;
		for (FeatureProcessor *f : fps_flat)
		{
			if (found)
				after_processors.push_back(f);
			else {
				if (processors_tp[f->processor_type])
					found = true;
				else
					++has_before;
			}
		}

		if (has_before > 0)
			MWARN("WARN:: AveragePredsPostProcessor :: found %d processors before\n", has_before);
		if (!after_processors.empty())
			MLOG("INFO:: AveragePredsPostProcessor :: found %zu processors after\n", after_processors.size());
	}
}

void AveragePredsPostProcessor::Learn(const MedFeatures &train_mat) {
	if (!boost::starts_with(feature_processor_type, "MODEL::")) {
		unordered_set<int> empt;
		MedFeatures mat = train_mat;
		feature_processor->Learn(mat, empt);
	}
}

void clear_map(map<string, vector<float>> &data) {
	for (auto &it : data)
		data[it.first].clear();
}

void AveragePredsPostProcessor::Apply(MedFeatures &matrix) const {
	//Apply plan, do in batches:
	//1. resample input - apply feature_processor multiple times for each sample (if imputer and using existing in model. need to know where are the missing values/get matrix without this processor in more general)
	// Currently will use imputation mask. TODO: later implement in a cleaner way
	//2. predict with model_predictor
	//3. aggregate predictions - fetch mean,median,std,ci - the rest in attributes
	vector<float> prctile_list = { (float)0.05, (float)0.5, (float)0.95 };

	if (boost::starts_with(feature_processor_type, "MODEL::") || force_cancel_imputations) {
		//matrix.masks validate generate_masks_for_features is on
		vector<bool> processors_tp(FTR_PROCESS_LAST);
		processors_tp[FTR_PROCESS_IMPUTER] = true;
		processors_tp[FTR_PROCESS_ITERATIVE_IMPUTER] = true;
		processors_tp[FTR_PROCESS_PREDICTOR_IMPUTER] = true;
		if (processors_tp[feature_processor->processor_type] || force_cancel_imputations) {
			//need to use masks
			if (matrix.masks.empty())
				MTHROW_AND_ERR("Error AveragePredsPostProcessor::Apply - model should run with generate_masks_for_features ON. no masks for imputations\n");
			//Erase imputed values - prepare for imputataions:
			vector<string> names;
			matrix.get_feature_names(names);
			for (size_t i = 0; i < names.size(); ++i)
			{
				const vector<unsigned char> &m = matrix.masks.at(names[i]);
				vector<float> &to_change = matrix.data.at(names[i]);
				if (m.empty())
					MTHROW_AND_ERR("Error AveragePredsPostProcessor::Apply - model should run with generate_masks_for_features ON. no masks for imputations\n");
				for (size_t j = 0; j < m.size(); ++j)
				{
					if (m[j] & matrix.imputed_mask)
						to_change[j] = matrix.medf_missing_value;
				}
			}
		}
	}

	//1. resample input - apply feature_processor multiple times for each sample	
	MedFeatures batch;
	batch.attributes = matrix.attributes;
	batch.tags = matrix.tags;
	batch.time_unit = matrix.time_unit;
	batch.medf_missing_value = matrix.medf_missing_value;
	batch.samples.resize(batch_size * resample_cnt);
	for (int i = 0; i < batch.samples.size(); ++i)
		batch.samples[i].id = i;
	batch.init_pid_pos_len();
	//data, samples
	int i = 0;
	vector<unordered_map<string, float>> samples_res(matrix.samples.size());
	while (i < matrix.samples.size())
	{
		//prepate batch
		int curr_sz = 0;
		MedFeatures apply_batch = batch;
		for (auto &it : matrix.data)
			apply_batch.data[it.first].resize(resample_cnt * batch_size);
		while (curr_sz < batch_size && i < matrix.samples.size()) {
			//add data from matrix
			for (auto &it : matrix.data) {
				for (size_t j = 0; j < resample_cnt; ++j)
					apply_batch.data[it.first][curr_sz*resample_cnt + j] = it.second[i];
			}
			++curr_sz;
			++i;
		}
		//by curr_sz:
		if (curr_sz < batch_size) {//last batch - remove last samples
			apply_batch.samples.resize(curr_sz*resample_cnt);
			apply_batch.init_pid_pos_len();
			for (auto &it : matrix.data)
				apply_batch.data[it.first].resize(resample_cnt*curr_sz);
		}
		//apply feature processor on all duplicated batch:
		feature_processor->apply(apply_batch);
		//Apply after batch
		for (FeatureProcessor *fp : after_processors)
			fp->apply(apply_batch);
		//apply batch with MedPredictor:
		model_predictor->predict(apply_batch);
		//collect preds from samples: each row was duplicated resample_cnt times
		vector<vector<float>> collected_preds(curr_sz); //for each sample
		for (size_t j = 0; j < curr_sz; ++j)
		{
			vector<float> &v = collected_preds[j];
			v.resize(resample_cnt);
			//add preds:
			for (size_t k = 0; k < resample_cnt; ++k)
				v[k] = apply_batch.samples[j*resample_cnt + k].prediction[0];
		}
		//aggregate results using collected_preds for each original pred:

		for (size_t j = 0; j < curr_sz; ++j) {
			unordered_map<string, float> &dict = samples_res[j];
			vector<float> &dt = collected_preds[j];
			float mean, std;
			medial::stats::get_mean_and_std_without_cleaning(dt, mean, std);
			dict["Mean"] = mean;
			dict["Std"] = std;
			vector<float> res;
			medial::stats::get_percentiles(dt, prctile_list, res);
			dict["CI_Lower"] = res[0];
			dict["Median"] = res[1];
			dict["CI_Upper"] = res[2];
		}
	}

	//store in final matrix:
	for (size_t j = 0; j < matrix.samples.size(); ++j)
	{
		const unordered_map<string, float> &dict = samples_res[j];
		MedSample &s = matrix.samples[j];

		if (use_median)
			s.prediction[0] = dict.at("Median");
		else
			s.prediction[0] = dict.at("Mean");

		//store in attributes the rest:
		s.attributes["pred.std"] = dict.at("Std");
		s.attributes["pred.ci_lower"] = dict.at("CI_Lower");
		s.attributes["pred.ci_upper"] = dict.at("CI_Upper");
		s.attributes["pred.mean"] = dict.at("Mean");
		s.attributes["pred.median"] = dict.at("Median");
	}

}

AveragePredsPostProcessor::~AveragePredsPostProcessor() {
	if (boost::starts_with(feature_processor_type, "MODEL::")) {
		if (feature_processor != NULL) {
			delete feature_processor;
			feature_processor = NULL;
		}
	}
}

void AveragePredsPostProcessor::dprint(const string &pref) const {
	MLOG("%s using %s preidctor, feature_processor of %s\n", pref.c_str(), model_predictor->my_class_name().c_str(),
		feature_processor->my_class_name().c_str());
}