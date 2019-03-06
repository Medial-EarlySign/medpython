#include "GibbsSampler.h"
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedMat/MedMat/MedMatConstants.h>

#define LOCAL_SECTION LOG_MEDSTAT

Gibbs_Params::Gibbs_Params() {
	burn_in_count = 1000;
	jump_between_samples = 10;
	samples_count = 1;
	predictor_type = "linear_model";
	predictor_args = "";
	predictors_counts = 10;
	selection_ratio = (float)0.7;
}

int Gibbs_Params::init(map<string, string>& map) {

	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "burn_in_count")
			burn_in_count = med_stoi(it->second);
		else if (it->first == "jump_between_samples")
			jump_between_samples = med_stoi(it->second);
		else if (it->first == "samples_count")
			samples_count = med_stoi(it->second);
		else if (it->first == "predictor_type")
			predictor_type = it->second;
		else if (it->first == "predictor_args")
			predictor_args = it->second;
		else if (it->first == "predictors_counts")
			predictors_counts = med_stoi(it->second);
		else if (it->first == "selection_ratio")
			selection_ratio = med_stof(it->second);
		else
			MTHROW_AND_ERR("Error in Gibbs_Params::init - no parameter \"%s\"\n", it->first.c_str());
	}

	return 0;
}

PredictorOrEmpty::PredictorOrEmpty() {
	random_device rd;
	gen = mt19937(rd());
}

float PredictorOrEmpty::get_sample(vector<float> &x) {
	if (!sample_cohort.empty()) {
		uniform_int_distribution<> rnd_gen(0, (int)sample_cohort.size() - 1);
		int sel = rnd_gen(gen);
		return sample_cohort[sel];
	}
	else if (!predictors.empty()) {
		uniform_int_distribution<> rnd_gen(0, (int)predictors.size() - 1);
		int sel = rnd_gen(gen);
		MedPredictor *sel_pred = predictors[sel];
		//get prediction and return 

		vector<float> pred_res;
		sel_pred->predict(x, pred_res, 1, (int)x.size());

		return pred_res[0];
	}

	MTHROW_AND_ERR("Error PredictorOrEmpty - not initialized");
}

//TODO: seperate the to prepare and sample - to learn variables predictors seperatlly and than just activate
void medial::stats::gibbs_sampling(const map<string, vector<float>> &cohort_data, const Gibbs_Params &params,
	map<string, vector<float>> &results, const vector<bool> *mask, const vector<float> *mask_values) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> rnd_num(0, 1);

	vector<bool> mask_f(cohort_data.size());
	vector<float> mask_values_f(cohort_data.size());
	if (mask == NULL)
		mask = &mask_f;
	if (mask_values == NULL) //and with init values
		mask_values = &mask_values_f;
	if (cohort_data.empty())
		MTHROW_AND_ERR("Error in medial::stats::gibbs_sampling - cohort_data can't be empty\n");
	//fix mask values and sample gibbs for the rest by cohort_data as statistical cohort for univariate marginal dist
	int sample_loop = params.burn_in_count + (params.samples_count - 1) * (params.jump_between_samples + 1) + 1;

	vector<string> all_names; all_names.reserve(cohort_data.size());
	for (auto it = cohort_data.begin(); it != cohort_data.end(); ++it)
		all_names.push_back(it->first);

	vector<float> current_sample(cohort_data.size());
	for (size_t i = 0; i < mask->size(); ++i)
	{
		if (mask->at(i))
			current_sample[i] = mask_values->at(i);
		else
			current_sample[i] = mask_values->at(i); //init value - not fixed to be this value
	}
	vector<int> idx_iter; idx_iter.reserve(mask->size());
	for (int i = 0; i < mask->size(); ++i)
		if (mask->at(i))
			idx_iter.push_back(i);
	int cohort_size = (int)cohort_data.begin()->second.size(); //assume not empty
	//build predictor for each unknown parameter in the mask based on known parameters + rest_unknown in the mask:
	vector<PredictorOrEmpty> feats_predictors(idx_iter.size());
	int pred_num_feats = (int)cohort_data.size() - 1;
	if (pred_num_feats == 0) {
		for (size_t i = 0; i < idx_iter.size(); ++i) {
			//just test for values as distribution mean, variance
			feats_predictors[i].sample_cohort = cohort_data.at(all_names[idx_iter[i]]);
		}
	}
	else {
		for (size_t i = 0; i < idx_iter.size(); ++i)
		{

			feats_predictors[i].predictors.resize(params.predictors_counts);
			for (size_t k = 0; k < params.predictors_counts; ++k)
			{
				//create predictors_count predictors on random selected samples
				feats_predictors[i].predictors[k] = MedPredictor::make_predictor(params.predictor_type, params.predictor_args);
				vector<float> train_vec, label_vec;
				train_vec.reserve(cohort_size * pred_num_feats); label_vec.reserve(cohort_size);
				for (size_t ii = 0; ii < cohort_size; ++ii) {
					float random_res = rnd_num(gen);
					if (random_res >= 1 - params.selection_ratio) {
						for (size_t jj = 0; jj < pred_num_feats; ++jj) {
							int fixed_idx = (int)jj + int(jj >= i); //skip current
							train_vec.push_back(cohort_data.at(all_names[fixed_idx])[ii]);
						}
						label_vec.push_back(cohort_data.at(all_names[idx_iter[i]])[ii]);
					}
				}
				//Learn Predictor:
				feats_predictors[i].predictors[k]->learn(train_vec, label_vec, cohort_size, pred_num_feats);
			}
		}
	}

	//can parallel for random init of initiale values (just burn in)
	//TODO: test each feature for distribution, normal, log_normal test without conditioning - assum it remains like that
	//For now all have normal dist that depends on all the other variables

	for (size_t i = 0; i < sample_loop; ++i)
	{
		//create sample - iterate over all variables not in mask:
		for (int idx = 0; idx < idx_iter.size(); ++idx)
		{
			vector<float> curr_x(pred_num_feats);
			for (size_t k = 0; k < curr_x.size(); ++k)
			{
				int fixxed_idx = (int)k + int(k >= i);
				curr_x[k] = current_sample[fixxed_idx];
			}
			float val = feats_predictors[idx].get_sample(curr_x); //based on dist (or predictor - value bin dist)
#pragma omp critical
			current_sample[idx_iter[idx]] = val; //update current pos variable
		}

		if (i >= params.burn_in_count && ((i - params.burn_in_count) % params.jump_between_samples) == 0) {
			//collect sample to result:
			for (size_t k = 0; k < all_names.size(); ++k)
				results[all_names[k]].push_back(current_sample[k]);
		}
	}


}