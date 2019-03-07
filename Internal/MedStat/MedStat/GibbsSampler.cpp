#include "GibbsSampler.h"
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedAlgo/MedAlgo/BinSplitOptimizer.h>
#include <MedMat/MedMat/MedMatConstants.h>
#include <omp.h>

#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL LOG_DEF_LEVEL

Gibbs_Params::Gibbs_Params() {
	burn_in_count = 1000;
	jump_between_samples = 10;
	samples_count = 1;
	predictor_type = "linear_model";
	predictor_args = "";
	predictors_counts = 10;
	selection_ratio = (float)0.7;
	select_with_repeats = false;
	kmeans = 0;
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
		else if (it->first == "select_with_repeats")
			select_with_repeats = med_stoi(it->second) > 0;
		else if (it->first == "kmeans")
			kmeans = med_stoi(it->second);
		else
			MTHROW_AND_ERR("Error in Gibbs_Params::init - no parameter \"%s\"\n", it->first.c_str());
	}

	return 0;
}

PredictorOrEmpty::PredictorOrEmpty() {
	random_device rd;
	gen = mt19937(rd());
}

int GibbsSampler::init(map<string, string>& map) {
	return params.init(map);
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

void GibbsSampler::learn_gibbs(const map<string, vector<float>> &cohort_data) {
	random_device rd;
	mt19937 gen(rd());
	if (params.selection_ratio > 1)
		MTHROW_AND_ERR("ERROR in GibbsSampler::learn_gibbs - params.selection_ratio can't be bigger than 1\n");

	vector<string> all_names; all_names.reserve(cohort_data.size());
	for (auto it = cohort_data.begin(); it != cohort_data.end(); ++it)
		all_names.push_back(it->first);
	all_feat_names = all_names;
	int cohort_size = (int)cohort_data.begin()->second.size(); //assume not empty
	uniform_int_distribution<> rnd_num(0, cohort_size - 1);

	feats_predictors.resize(all_names.size());
	int pred_num_feats = (int)cohort_data.size() - 1;
	if (pred_num_feats == 0) {
		for (size_t i = 0; i < all_names.size(); ++i) {
			//just test for values as distribution mean, variance
			feats_predictors[i].sample_cohort = cohort_data.at(all_names[i]);
		}
	}
	else {
		vector<int> clusters;
		if (params.kmeans > 0) {
			//seperate the X space to k clusters:
			vector<float> full_vec(cohort_size * all_names.size());
			for (size_t i = 0; i < cohort_size; ++i)
				for (size_t j = 0; j < all_names.size(); ++j)
					full_vec[i * pred_num_feats + j] = cohort_data.at(all_names[j])[i];
			int k = params.kmeans;
			vector<float> centers(k * all_names.size()), dists(k * cohort_size);
			clusters.resize(cohort_size);
			KMeans(full_vec.data(), cohort_size, (int)all_names.size(), k, 1000, centers.data(), clusters.data(), dists.data());
		}

		MedTimer tm;
		tm.start();
		chrono::high_resolution_clock::time_point tm_prog = chrono::high_resolution_clock::now();
		int progress = 0;
		int max_loop = (int)all_names.size() * params.predictors_counts;

		int train_sz = int(cohort_size * params.selection_ratio);
		for (size_t i = 0; i < all_names.size(); ++i)
		{
			feats_predictors[i].predictors.resize(params.predictors_counts);
			for (size_t k = 0; k < params.predictors_counts; ++k)
			{
				//create predictors_count predictors on random selected samples
				feats_predictors[i].predictors[k] = MedPredictor::make_predictor(params.predictor_type, params.predictor_args);
				vector<float> train_vec(train_sz * pred_num_feats), label_vec(train_sz);
				vector<bool> seen;
				if (!params.select_with_repeats)
					seen.resize(cohort_size);
				vector<int> sel_ls;
				if (params.kmeans > 0)
					sel_ls.resize(train_sz);
				for (size_t ii = 0; ii < train_sz; ++ii) {
					int random_idx = rnd_num(gen);
					if (!params.select_with_repeats) { //if need to validate no repeats - do it
						while (seen[random_idx])
							random_idx = rnd_num(gen);
						seen[random_idx] = true;
					}
					for (size_t jj = 0; jj < pred_num_feats; ++jj) {
						int fixed_idx = (int)jj + int(jj >= i); //skip current
						train_vec[ii* pred_num_feats + jj] = (cohort_data.at(all_names[fixed_idx])[random_idx]);
					}
					label_vec[ii] = cohort_data.at(all_names[i])[random_idx];
					if (params.kmeans > 0)
						sel_ls[ii] = random_idx;
				}
				//Learn Predictor:
				if (params.kmeans > 0) {
					//randomize y for each cluster - select random y from cluster to learn
					vector<vector<float>> cluster_labels(params.kmeans);
					for (size_t kk = 0; kk < label_vec.size(); ++kk)
					{
						int cluster_id = clusters[sel_ls[kk]];
						cluster_labels[cluster_id].push_back(label_vec[kk]);
					}
					vector<float> cluster_sel(params.kmeans);
					for (size_t kk = 0; kk < params.kmeans; ++kk)
					{
						if (cluster_labels[kk].empty())
							continue;
						uniform_int_distribution<> sel_rnd(0, (int)cluster_labels[kk].size() - 1);
						cluster_sel[kk] = cluster_labels[kk][sel_rnd(gen)];
					}
					//override labels by selection:
					for (size_t kk = 0; kk < label_vec.size(); ++kk)
					{
						int cluster_id = clusters[sel_ls[kk]];
						label_vec[kk] = cluster_sel[cluster_id];
					}
				}
				feats_predictors[i].predictors[k]->learn(train_vec, label_vec, train_sz, pred_num_feats);

				++progress;
				double duration = (unsigned long long)(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
					- tm_prog).count()) / 1000000.0;
				if (duration > 30) {
#pragma omp critical
					tm_prog = chrono::high_resolution_clock::now();
					double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
						- tm.t[0]).count()) / 1000000.0;
					double estimate_time = int(double(max_loop - progress) / double(progress) * double(time_elapsed));
					MLOG("Processed %d out of %d(%2.2f%%) time elapsed: %2.1f Minutes, "
						"estimate time to finish %2.1f Minutes\n",
						progress, (int)max_loop, 100.0*(progress / float(max_loop)), time_elapsed / 60,
						estimate_time / 60.0);
				}
			}
		}
	}
}

void GibbsSampler::get_samples(map<string, vector<float>> &results, const vector<bool> *mask, const vector<float> *mask_values) {

	vector<bool> mask_f(all_feat_names.size());
	vector<float> mask_values_f(all_feat_names.size());
	if (mask == NULL)
		mask = &mask_f;
	if (mask_values == NULL) //and with init values
		mask_values = &mask_values_f;
	if (all_feat_names.empty())
		MTHROW_AND_ERR("Error in medial::stats::gibbs_sampling - cohort_data can't be empty\n");
	//fix mask values and sample gibbs for the rest by cohort_data as statistical cohort for univariate marginal dist
	int sample_loop = params.burn_in_count + (params.samples_count - 1) * (params.jump_between_samples + 1) + 1;

	vector<string> &all_names = all_feat_names;

	vector<float> current_sample(all_feat_names.size());
	for (size_t i = 0; i < mask->size(); ++i)
	{
		if (mask->at(i))
			current_sample[i] = mask_values->at(i);
		else
			current_sample[i] = mask_values->at(i); //init value - not fixed to be this value
	}
	vector<int> idx_iter; idx_iter.reserve(mask->size());
	for (int i = 0; i < mask->size(); ++i)
		if (!mask->at(i))
			idx_iter.push_back(i);
	int pred_num_feats = (int)all_feat_names.size() - 1;

	//can parallel for random init of initiale values (just burn in)
	for (size_t i = 0; i < sample_loop; ++i)
	{
		//create sample - iterate over all variables not in mask:
		for (int idx = 0; idx < idx_iter.size(); ++idx)
		{
			int f_idx = idx_iter[idx];
			vector<float> curr_x(pred_num_feats);
			for (size_t k = 0; k < curr_x.size(); ++k)
			{
				int fixxed_idx = (int)k + int(k >= i);
				curr_x[k] = current_sample[fixxed_idx];
			}
			float val = feats_predictors[f_idx].get_sample(curr_x); //based on dist (or predictor - value bin dist)
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

void GibbsSampler::get_parallel_samples(map<string, vector<float>> &results, uniform_real_distribution<> &real_dist,
	const vector<bool> *mask) {
	random_device rd;
	mt19937 gen;
	vector<bool> mask_f(all_feat_names.size());

	if (mask == NULL)
		mask = &mask_f;
	if (all_feat_names.empty())
		MTHROW_AND_ERR("Error in medial::stats::gibbs_sampling - cohort_data can't be empty\n");
	int N_tot_threads = omp_get_max_threads();
	vector<GibbsSampler> copy_gibbs(N_tot_threads);
	for (size_t i = 0; i < copy_gibbs.size(); ++i) {
		copy_gibbs[i] = *this;
		copy_gibbs[i].params.samples_count = 1;
	}

	MedTimer tm;
	tm.start();
	chrono::high_resolution_clock::time_point tm_prog = chrono::high_resolution_clock::now();
	int progress = 0;
	int max_loop = params.samples_count;
#pragma omp parallel for
	for (int i = 0; i < params.samples_count; ++i)
	{
		int n_th = omp_get_thread_num();
		GibbsSampler &g = copy_gibbs[n_th];

		vector<float> mask_vals(all_feat_names.size());
		for (size_t i = 0; i < mask_vals.size(); ++i)
			mask_vals[i] = real_dist(gen);
		map<string, vector<float>> res;

		g.get_samples(res, mask, &mask_vals);

#pragma omp critical
		for (auto it = res.begin(); it != res.end(); ++it)
			results[it->first].push_back(it->second[0]);

#pragma omp atomic
		++progress;
		double duration = (unsigned long long)(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
			- tm_prog).count()) / 1000000.0;
		if (duration > 30 && progress % 50 == 0) {
#pragma omp critical
			tm_prog = chrono::high_resolution_clock::now();
			double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
				- tm.t[0]).count()) / 1000000.0;
			double estimate_time = int(double(max_loop - progress) / double(progress) * double(time_elapsed));
			MLOG("Processed %d out of %d(%2.2f%%) time elapsed: %2.1f Minutes, "
				"estimate time to finish %2.1f Minutes\n",
				progress, (int)max_loop, 100.0*(progress / float(max_loop)), time_elapsed / 60,
				estimate_time / 60.0);
		}
	}

}