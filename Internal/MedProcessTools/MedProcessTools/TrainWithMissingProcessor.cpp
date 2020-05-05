#include "TrainWithMissingProcessor.h"
#include "ExplainWrapper.h"

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void TrainMissingProcessor::init_defaults() {
	missing_value = MED_MAT_MISSING_VALUE;
	add_new_data = 0;
	sample_masks_with_repeats = true;
	uniform_rand = false;
	use_shuffle = true;
	subsample_train = 0;
	limit_mask_size = 0;
	verbose = true;
	processor_type = FTR_PROCESS_ADD_MISSING_TO_LEARN;
	grouping = "";
}

int TrainMissingProcessor::init(map<string, string>& mapper) {
	for (const auto &it : mapper)
	{
		//! [TrainMissingProcessor::init]
		if (it.first == "missing_value")
			missing_value = med_stof(it.second);
		else if (it.first == "add_new_data")
			add_new_data = med_stoi(it.second);
		else if (it.first == "sample_masks_with_repeats")
			sample_masks_with_repeats = med_stoi(it.second) > 0;
		else if (it.first == "uniform_rand")
			uniform_rand = med_stoi(it.second) > 0;
		else if (it.first == "use_shuffle")
			use_shuffle = med_stoi(it.second) > 0;
		else if (it.first == "subsample_train")
			subsample_train = med_stoi(it.second);
		else if (it.first == "limit_mask_size")
			limit_mask_size = med_stoi(it.second);
		else if (it.first == "verbose")
			verbose = med_stoi(it.second) > 0;
		else if (it.first == "grouping")
			grouping = it.second;
		else if (it.first == "fp_type" || it.first == "use_parallel_learn" || it.first == "use_parallel_apply") {}
		else
			MTHROW_AND_ERR("Error in TrainMissingProcessor::init - unsupported argument %s\n", it.first.c_str());
		//! [TrainMissingProcessor::init]
	}

	return 0;
}

void TrainMissingProcessor::dprint(const string &pref, int fp_flag) {
	string res = this->object_json();
	boost::replace_all(res, "\n", " ");
	MLOG("%s :: %s\n", pref.c_str(), res.c_str());
}

int _count_msn(const float *vals, int sz, float val) {
	int res = 0;
	for (size_t i = 0; i < sz; ++i)
		res += int(vals[i] == val);
	return res;
}

int TrainMissingProcessor::Learn(MedFeatures& features, unordered_set<int>& ids) {
	vector<string> features_nms;
	features.get_feature_names(features_nms);
	if (!grouping.empty())
		ExplainProcessings::read_feature_grouping(grouping, features_nms, group2Inds, groupNames);
	//use group2Inds, groupNames

	if (limit_mask_size >= group2Inds.size()) {
		MWARN("WARNING: limit_mask_size=%d which is bigger than number of groups/features(%zu)\n",
			limit_mask_size, group2Inds.size());
		limit_mask_size = (int)group2Inds.size(); //problem with arguments
	}

	mt19937 gen(globalRNG::rand());
	int nftrs_grp = (int)group2Inds.size();
	int nftrs = (int)features.data.size();
	int train_mat_size = (int)features.samples.size();

	vector<int> missing_hist(nftrs + 1), added_missing_hist(nftrs + 1), added_grp_hist(nftrs_grp + 1);
	MedMat<float> x_mat;
	vector<int> original_samples_id(features.samples.size());
	features.get_as_matrix(x_mat);
	for (int i = 0; i < original_samples_id.size(); ++i)
		original_samples_id[i] = i;

	vector<int> miss_cnts(train_mat_size + add_new_data);
	vector<int> mask_group_sizes(train_mat_size + add_new_data); //stores for each sample - how many missings in groups manner:
	for (size_t i = 0; i < train_mat_size; ++i)
	{
		//check how many groups missings:
		int grp_misses = 0;
		for (int j = 0; j < nftrs_grp; ++j) {
			bool has_missing = false;
			for (size_t k = 0; k < group2Inds[j].size() && !has_missing; ++k)
				has_missing = x_mat(i, group2Inds[j][k]) == missing_value;
			grp_misses += int(has_missing);
		}
		mask_group_sizes[i] = grp_misses;
	}

	if (add_new_data > 0) {
		original_samples_id.reserve(original_samples_id.size() + add_new_data);
		vector<float> rows_m(add_new_data * nftrs);
		unordered_set<vector<bool>> seen_mask;
		uniform_int_distribution<> rnd_row(0, train_mat_size - 1);
		double log_max_opts = log(add_new_data) / log(2.0);
		if (log_max_opts >= nftrs_grp) {
			if (!sample_masks_with_repeats)
				MWARN("Warning: you have request to sample masks without repeats, but it can't be done. setting sample with repeats\n");
			sample_masks_with_repeats = true;
		}
		if (verbose)
			MLOG("Adding %d Data points (has %d features with %d groups)\n", add_new_data, nftrs, nftrs_grp);
		MedProgress add_progress("Add_Train_Data", add_new_data, 30, 1);
		for (size_t i = 0; i < add_new_data; ++i)
		{
			float *curr_row = &rows_m[i *  nftrs];
			//select row:
			int row_sel = rnd_row(gen);

			vector<bool> curr_mask; curr_mask.resize(nftrs_grp);
			for (int j = 0; j < nftrs_grp; ++j) {
				bool has_missing = false;
				for (size_t k = 0; k < group2Inds[j].size() && !has_missing; ++k)
					has_missing = x_mat(row_sel, group2Inds[j][k]) == missing_value;
				curr_mask[j] = !has_missing;
			}

			medial::shapley::generate_mask_(curr_mask, nftrs_grp, gen, uniform_rand, use_shuffle, limit_mask_size);
			while (!sample_masks_with_repeats && seen_mask.find(curr_mask) != seen_mask.end())
				medial::shapley::generate_mask_(curr_mask, nftrs_grp, gen, uniform_rand, use_shuffle, limit_mask_size);
			if (!sample_masks_with_repeats)
				seen_mask.insert(curr_mask);

			//commit mask to curr_row
			int msn_cnt = 0;
			for (int j = 0; j < nftrs_grp; ++j)
			{
				if (curr_mask[j]) {
					for (size_t k = 0; k < group2Inds[j].size(); ++k)
						curr_row[group2Inds[j][k]] = x_mat(row_sel, group2Inds[j][k]);
				}
				else {
					for (size_t k = 0; k < group2Inds[j].size(); ++k)
						curr_row[group2Inds[j][k]] = missing_value;
				}
				msn_cnt += int(!curr_mask[j]); //how many missings
			}
			original_samples_id.push_back(row_sel);
			++added_grp_hist[msn_cnt];
			add_progress.update();
			mask_group_sizes[train_mat_size + i] = msn_cnt;
		}
		x_mat.add_rows(rows_m);
	}

	// Add data with missing values according to sample masks
	vector<int> grp_missing_hist_all(nftrs_grp + 1);

	for (int i = 0; i < x_mat.nrows; ++i) {
		miss_cnts[i] = _count_msn(x_mat.data_ptr(i, 0), nftrs, missing_value);
		++missing_hist[miss_cnts[i]];
		if (i >= train_mat_size)
			++added_missing_hist[miss_cnts[i]];
		++grp_missing_hist_all[mask_group_sizes[i]];
	}
	if (verbose) {
		medial::print::print_hist_vec(miss_cnts, "missing_values_cnt percentiles [0 - " + to_string(nftrs) + "] (with added samples - no groups)", "%d");
		medial::print::print_hist_vec(mask_group_sizes, "mask_group_sizes percentiles [0 - " + to_string(nftrs_grp) + "] (with added samples - for groups)", "%d");
		medial::print::print_hist_vec(added_missing_hist, "selected counts in hist of missing_values_cnt (only for added - no groups)", "%d");
		if (added_grp_hist.size() < 50)
			medial::print::print_vec(added_grp_hist, "grp hist (only for added - on groups)", "%d");
		else
			medial::print::print_hist_vec(added_grp_hist, "hist of added_grp_hist (only for added - on groups)", "%d");
	}

	if (subsample_train > 0 && subsample_train < train_mat_size) {
		//do subsampling:
		MLOG("INFO:: TrainMissingProcessor::Learn - subsampling original train matrix");
		unordered_set<int> selected_idx;

		uniform_int_distribution<> rnd_opts(0, train_mat_size - 1);
		for (size_t i = 0; i < subsample_train; ++i)
		{
			int sel_idx = rnd_opts(gen);
			while (selected_idx.find(sel_idx) != selected_idx.end())
				sel_idx = rnd_opts(gen);
			selected_idx.insert(sel_idx);
		}
		//add all rest:
		vector<int> empty_is_all;
		vector<int> selected_idx_vec(selected_idx.begin(), selected_idx.end());
		for (int i = train_mat_size; i < x_mat.nrows; ++i)
			selected_idx_vec.push_back(i);
		sort(selected_idx_vec.begin(), selected_idx_vec.end());

		//commit selection xmat and labels, weights:
		vector<int> new_samples_ids(selected_idx_vec.size());
		x_mat.get_sub_mat(selected_idx_vec, empty_is_all);
		for (size_t i = 0; i < selected_idx_vec.size(); ++i)
			new_samples_ids[i] = original_samples_id[selected_idx_vec[i]];
		original_samples_id = move(new_samples_ids);
	}
	//use x_mat, original_samples_id to manipulate features => data and samples:
	vector<MedSample> new_smps(original_samples_id.size());
	map<string, vector<float>> new_data;
	for (size_t i = 0; i < original_samples_id.size(); ++i)
		new_smps[i] = features.samples[original_samples_id[i]];
	for (size_t i = 0; i < features_nms.size(); ++i)
	{
		vector<float> &vec = new_data[features_nms[i]];
		vec.resize(x_mat.nrows);
		for (size_t j = 0; j < x_mat.nrows; ++j)
			vec[j] = x_mat(j, i);
	}

	features.samples = move(new_smps);
	features.data = move(new_data);
	features.init_pid_pos_len();

	return 0;
}