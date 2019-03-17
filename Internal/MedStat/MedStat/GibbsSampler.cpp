#include "GibbsSampler.h"
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedAlgo/MedAlgo/BinSplitOptimizer.h>
#include <MedMat/MedMat/MedMatConstants.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <regex>

#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL LOG_DEF_LEVEL

Gibbs_Params::Gibbs_Params() {
	burn_in_count = 1000;
	jump_between_samples = 10;
	samples_count = 1;
	kmeans = 0;
	max_iters = 500;
	selection_ratio = (float)0.7;
	find_real_value_bin = true;
	select_with_repeats = false;

	predictor_type = "lightgbm";
	predictor_args = "objective=multiclassova;metric=multi_logloss;verbose=0;num_threads=15;"
		"num_trees=250;learning_rate=0.05;lambda_l2=0;metric_freq=50;num_class=50";
	num_class_setup = "num_class";
	bin_settings.init_from_string("split_method=iterative_merge;binCnt=50");
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
		else if (it->first == "kmeans")
			kmeans = med_stoi(it->second);
		else if (it->first == "max_iters")
			max_iters = med_stoi(it->second);
		else if (it->first == "selection_ratio")
			selection_ratio = med_stof(it->second);
		else if (it->first == "find_real_value_bin")
			find_real_value_bin = med_stoi(it->second) > 0;
		else if (it->first == "select_with_repeats")
			select_with_repeats = med_stoi(it->second) > 0;
		else if (it->first == "predictor_type")
			predictor_type = it->second;
		else if (it->first == "predictor_args")
			predictor_args = it->second;
		else if (it->first == "num_class_setup")
			num_class_setup = it->second;
		else if (it->first == "bin_settings")
			bin_settings.init_from_string(it->second);
		else
			MTHROW_AND_ERR("Error in Gibbs_Params::init - no parameter \"%s\"\n", it->first.c_str());
	}

	return 0;
}

PredictorOrEmpty::PredictorOrEmpty() {
	predictor = NULL;
}

PredictorOrEmpty::~PredictorOrEmpty() {
	if (predictor != NULL)
		delete predictor;
	predictor = NULL;
}

int GibbsSampler::init(map<string, string>& map) {
	return params.init(map);
}

float PredictorOrEmpty::get_sample(vector<float> &x, mt19937 &gen) const {
	if (!sample_cohort.empty()) {
		uniform_int_distribution<> rnd_gen(0, (int)sample_cohort.size() - 1);
		int sel = rnd_gen(gen);
		return sample_cohort[sel];
	}
	else if (!cluster_centers.empty()) {
		//find closet cluster:
		int close_idx = -1;
		float min_diff = -1;
		int k = (int)cluster_centers.size() / input_size;

		for (size_t i = 0; i < k; ++i)
		{
			if (clusters_y[i].empty())
				continue; // skip empty clusters
			float curr_dif = 0;
			for (size_t j = 0; j < input_size; ++j)
				curr_dif += pow(cluster_centers[i * input_size + j] - x[j], 2);
			if (close_idx == -1 || min_diff > curr_dif) {
				min_diff = curr_dif;
				close_idx = (int)i;
			}
		}
		//the relevant cluster is close_idx - let's randomize y value from it:
		const vector<float> &sample_from = clusters_y[close_idx];
		uniform_int_distribution<> rnd_gen(0, (int)sample_from.size() - 1);
		int sel = rnd_gen(gen);
		return sample_from[sel];
	}
	else if (predictor != NULL) {
		vector<float> prd; //for each class:
		predictor->predict(x, prd, 1, (int)x.size());
		float tot_num = 0;
		for (size_t i = 0; i < prd.size(); ++i)
			tot_num += prd[i];
		uniform_real_distribution<> real_dist(0, tot_num);
		float sel = real_dist(gen);

		//now select correspond bin value:
		tot_num = 0;
		int sel_idx = 0;
		while (sel_idx < prd.size() && tot_num + prd[sel_idx] < sel) {
			tot_num += prd[sel_idx];
			++sel_idx;
		}

		return bin_vals[sel_idx];
	}

	MTHROW_AND_ERR("Error PredictorOrEmpty - not initialized");
}

GibbsSampler::GibbsSampler() {
	random_device rd;
	_gen = mt19937(rd());
}

void GibbsSampler::learn_gibbs(const map<string, vector<float>> &cohort_data) {
	random_device rd;
	mt19937 gen(rd());
	if (params.selection_ratio > 1) {
		MWARN("Warning - GibbsSampler::learn_gibbs - params.selection_ratio is bigger than 1 - setting to 1");
		params.selection_ratio = 1;
	}

	vector<string> all_names; all_names.reserve(cohort_data.size());
	for (auto it = cohort_data.begin(); it != cohort_data.end(); ++it)
		all_names.push_back(it->first);
	all_feat_names = all_names;
	int cohort_size = (int)cohort_data.begin()->second.size(); //assume not empty
	uniform_int_distribution<> rnd_num(0, cohort_size - 1);

	feats_predictors.resize(all_names.size());
	uniqu_value_bins.resize(all_names.size());
	int pred_num_feats = (int)cohort_data.size() - 1;
	for (size_t i = 0; i < all_names.size(); ++i)
		feats_predictors[i].input_size = pred_num_feats;
	if (pred_num_feats == 0) {
		for (size_t i = 0; i < all_names.size(); ++i) {
			//just test for values as distribution mean, variance
			feats_predictors[i].sample_cohort = cohort_data.at(all_names[i]);
		}
	}
	else {
		MedTimer tm;
		tm.start();
		chrono::high_resolution_clock::time_point tm_prog = chrono::high_resolution_clock::now();
		int progress = 0;

		int max_loop = (int)all_names.size();
		int train_sz = int(cohort_size * params.selection_ratio);

#pragma omp parallel for
		for (int i = 0; i < all_names.size(); ++i)
		{
			unordered_set<float> uniq_vals;
			for (size_t k = 0; k < cohort_data.at(all_names[i]).size(); ++k)
				uniq_vals.insert(cohort_data.at(all_names[i])[k]);
#pragma omp critical
			{
				uniqu_value_bins[i].insert(uniqu_value_bins[i].end(), uniq_vals.begin(), uniq_vals.end());
				sort(uniqu_value_bins[i].begin(), uniqu_value_bins[i].end());
			}

			vector<int> clusters;
			vector<float> train_vec(train_sz * pred_num_feats), label_vec(train_sz);
			vector<bool> seen;
			if (!params.select_with_repeats)
				seen.resize(cohort_size);
			vector<int> sel_ls(train_sz);
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
				sel_ls[ii] = random_idx;
			}
			if (params.kmeans > 0) {
				//seperate the X space to k clusters:
				int k = params.kmeans;
				if (INT_MAX / train_sz < k) {
					k = INT_MAX / train_sz - 1;
					MWARN("Warning: k=%d for kMeans is too large for that sample size of %d shrinking k to %d\n",
						params.kmeans, train_sz, k);
				}
				//vector<float> centers(k * pred_num_feats);
#pragma omp critical 
				{
					feats_predictors[i].cluster_centers.resize(k * pred_num_feats);
					feats_predictors[i].clusters_y.resize(k);
				}
				vector<float> dists(k * train_sz);
				clusters.resize(train_sz);
				//MLOG("Running kMeans for %s (%zu / %zu)\n", all_names[i].c_str(), i + 1, all_names.size());
				KMeans(train_vec.data(), train_sz, pred_num_feats, k, params.max_iters,
					feats_predictors[i].cluster_centers.data(), clusters.data(), dists.data(), false);
				//calc feats_predictors[i].clusters_y:
#pragma omp critical 
				for (size_t j = 0; j < train_sz; ++j)
					feats_predictors[i].clusters_y[clusters[j]].push_back(cohort_data.at(all_names[i])[j]);
			}
			else {

				//use predictors to train on train_vec and predcit on label_vec:
				//do binning for label_vec:
				vector<int> empt;
				medial::process::split_feature_to_bins(params.bin_settings, label_vec, empt, label_vec);
				//count num of classes:
				unordered_set<float> seen_val;
				for (size_t ii = 0; ii < label_vec.size(); ++ii)
					seen_val.insert(label_vec[ii]);
				int class_num = (int)seen_val.size();
				string predictor_init = params.predictor_args;
				//set num classes if needed:
				string empty_str = "";
				if (!params.num_class_setup.empty()) {
					//std::regex rgx(params.num_class_setup + "=[^;]+");
					//predictor_init = std::regex_replace(predictor_init, rgx, empty_str);
					//boost::replace_all(predictor_init, ";;", ";");
					predictor_init += ";" + params.num_class_setup + "=" + to_string(class_num);
					//change predictor_init
				}
				//init predictor
				MedPredictor *train_pred = MedPredictor::make_predictor(params.predictor_type, predictor_init);

				vector<float> sorted_vals(seen_val.begin(), seen_val.end());
				sort(sorted_vals.begin(), sorted_vals.end());
				MLOG("Feature %s has %d categories\n", all_names[i].c_str(), class_num);

#pragma omp critical
				feats_predictors[i].bin_vals.insert(feats_predictors[i].bin_vals.end(), sorted_vals.begin(), sorted_vals.end());
				//learn predictor
				//change labels to be 0 to K-1:
				unordered_map<float, int> map_categ;
				//calc by order:
				for (size_t ii = 0; ii < sorted_vals.size(); ++ii)
					map_categ[sorted_vals[ii]] = (int)ii;
				//commit:
				for (size_t ii = 0; ii < label_vec.size(); ++ii)
					label_vec[ii] = (float)map_categ.at(label_vec[ii]);

				train_pred->learn(train_vec, label_vec, (int)label_vec.size(), pred_num_feats);
#pragma omp critical
				feats_predictors[i].predictor = train_pred;
			}

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
			float val = feats_predictors[f_idx].get_sample(curr_x, _gen); //based on dist (or predictor - value bin dist)
			//find best bin:
			if (params.find_real_value_bin) {
				int pos = medial::process::binary_search_position(uniqu_value_bins[f_idx].data(), uniqu_value_bins[f_idx].data() + uniqu_value_bins[f_idx].size() - 1, val);
				if (pos == 0)
					val = uniqu_value_bins[f_idx][0];
				else {
					float diff_next = abs(val - uniqu_value_bins[f_idx][pos]);
					float diff_prev = abs(val - uniqu_value_bins[f_idx][pos - 1]);
					if (diff_prev < diff_next)
						val = uniqu_value_bins[f_idx][pos - 1];
					else
						val = uniqu_value_bins[f_idx][pos];
				}
			}
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
		copy_gibbs[i]._gen = mt19937(rd());
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
			mask_vals[i] = real_dist(g._gen);
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

	for (size_t i = 0; i < copy_gibbs.size(); ++i)
		for (size_t j = 0; j < copy_gibbs[i].feats_predictors.size(); ++j)
			copy_gibbs[i].feats_predictors[j].predictor = NULL; //that won't be cleaned from memory here - just a copy

}

void GibbsSampler::filter_samples(const map<string, vector<float>> &cohort_data,
	map<string, vector<float>> &results, const string &predictor_type, const string &predictor_args, float filter_sens) {
	random_device rd;
	mt19937 gen(rd());

	MedFeatures new_data;
	for (auto it = cohort_data.begin(); it != cohort_data.end(); ++it)
		new_data.attributes[it->first].normalized = false;

	int cohort_size = (int)cohort_data.begin()->second.size();
	new_data.samples.resize(cohort_size + results.begin()->second.size());
	for (size_t i = 0; i < new_data.samples.size(); ++i) {
		new_data.samples[i].id = (int)i;
		new_data.samples[i].outcome = (float)int(i < cohort_size);
	}
	//change outcome to be population label: is population 1?
	vector<float> labels(new_data.samples.size());
	for (size_t i = 0; i < new_data.samples.size(); ++i)
		labels[i] = new_data.samples[i].outcome;
	new_data.init_pid_pos_len();
	for (auto it = cohort_data.begin(); it != cohort_data.end(); ++it)
	{
		new_data.data[it->first] = it->second;
		new_data.data[it->first].insert(new_data.data[it->first].end(),
			results.at(it->first).begin(), results.at(it->first).end());
	}

	//lets get auc on this problem:
	MedPredictor *predictor = MedPredictor::make_predictor(predictor_type, predictor_args);
	//lets fix labels weight that cases will be less common
	vector<float> preds;
	medial::models::get_pids_cv(predictor, new_data, 5, gen, preds);

	float auc = medial::performance::auc_q(preds, labels, &new_data.weights);
	MLOG("predictor AUC with CV to diffrentiate between populations is %2.3f\n", auc);

	//do filter: take FPR on SENS
	unordered_map<float, vector<int>> pred_idx;
	vector<float> sorted_preds;
	for (size_t i = 0; i < preds.size(); ++i) {
		if (pred_idx.find(preds[i]) == pred_idx.end())
			sorted_preds.push_back(preds[i]);
		pred_idx[preds[i]].push_back((int)i);
	}
	sort(sorted_preds.begin(), sorted_preds.end());

	double t_cnt = 0;
	double f_cnt = 0;
	double tot_true_labels = cohort_size;
	double tot_false_labels = results.begin()->second.size();
	vector<float> true_rate = vector<float>((int)sorted_preds.size());
	vector<float> false_rate = vector<float>((int)sorted_preds.size());
	int st_size = (int)sorted_preds.size() - 1;
	for (int i = st_size; i >= 0; --i)
	{
		const vector<int> &indexes = pred_idx[sorted_preds[i]];
		//calc AUC status for this step:
		for (int ind : indexes)
		{
			bool true_label = labels[ind] > 0;
			t_cnt += int(true_label);
			f_cnt += int(!true_label);
		}
		true_rate[st_size - i] = float(t_cnt / tot_true_labels);
		false_rate[st_size - i] = float(f_cnt / tot_false_labels);
	}

	//stop on SENS point:
	int stop_idx = 0;
	while (stop_idx < true_rate.size() && true_rate[stop_idx] < filter_sens)
		++stop_idx;
	if (stop_idx >= true_rate.size())
		--stop_idx;
	stop_idx = st_size - stop_idx;
	//collect all indexes above that score
	vector<int> filter_sel;
	for (int i = st_size; i >= stop_idx; --i)
	{
		const vector<int> &indexes = pred_idx[sorted_preds[i]];
		for (int ind : indexes)
			if (!labels[ind])
				filter_sel.push_back(ind - cohort_size);
	}

	//commit selection:
	map<string, vector<float>> filterd;
	for (auto it = results.begin(); it != results.end(); ++it) {
		filterd[it->first].resize(filter_sel.size());
		for (size_t i = 0; i < filter_sel.size(); ++i)
			filterd[it->first][i] = (it->second[filter_sel[i]]);
	}

	results.swap(filterd);

}