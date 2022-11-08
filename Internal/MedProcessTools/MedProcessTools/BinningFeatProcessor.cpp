#include "FeatureProcess.h"

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL


int Binning_Wrapper::init(map<string, string>& mapper) {
	for (const auto &it : mapper)
	{
		//! [Binning_Wrapper::init]
		if (it.first == "bin_cutoffs") {
			vector<string> num_str;
			boost::split(num_str, it.second, boost::is_any_of(","));
			bin_cutoffs.resize(num_str.size());
			for (size_t i = 0; i < num_str.size(); ++i)
				bin_cutoffs[i] = stod(num_str[i]);
		}
		else if (it.first == "bin_repr_vals") {
			vector<string> num_str;
			boost::split(num_str, it.second, boost::is_any_of(","));
			bin_repr_vals.resize(num_str.size());
			for (size_t i = 0; i < num_str.size(); ++i)
				bin_repr_vals[i] = stod(num_str[i]);
		}
		else if (it.first == "use_bin_settings")
			use_bin_settings = it.second;
		//! [Binning_Wrapper::init]
		else
			MTHROW_AND_ERR("Error Binning::init - unknown arg \"%s\"\n", it.first.c_str());
	}

	if (bin_repr_vals.empty() && !bin_cutoffs.empty()) {
		bin_repr_vals.resize(bin_cutoffs.size() + 1);
		for (size_t i = 0; i < bin_cutoffs.size(); ++i)
			bin_repr_vals[i] = bin_cutoffs[i];
		bin_repr_vals[bin_cutoffs.size()] = bin_cutoffs.back() + 1e-5;
	}

	return 0;
}

void Binning_Wrapper::load_bin_settings(const vector<float> &nums, vector<float> &y) {
	if (use_bin_settings.empty())
		return;

	BinSettings bs;
	bs.init_from_string(use_bin_settings);
	vector<float> n = nums;
	vector<int> all;
	medial::process::split_feature_to_bins(bs, n, all, y);

	set<float> uniq_nums(n.begin(), n.end());

	bin_cutoffs.resize(uniq_nums.size() + 1);
	int ii = 1;
	bin_cutoffs[0] = *uniq_nums.begin() - 1e-5;
	for (float nm : uniq_nums)
	{
		bin_cutoffs[ii] = nm;
		++ii;
	}

	bin_repr_vals.resize(bin_cutoffs.size() + 1);
	for (size_t i = 0; i < bin_cutoffs.size(); ++i)
		bin_repr_vals[i] = bin_cutoffs[i];
	bin_repr_vals[bin_cutoffs.size()] = bin_cutoffs.back() + 1e-5;
}

int Binning_Wrapper::get_idx(float v) const {
	int res = medial::process::binary_search_position(bin_cutoffs, (double)v);
	return res;
}

int Binning_Wrapper::num_of_bins() const {
	return (int)bin_cutoffs.size() + 1;
}

float Binning_Wrapper::normalize(float v) const {
	int idx = get_idx(v);
	return bin_repr_vals[idx];
}

int BinningFeatProcessor::init(map<string, string>& mapper) {
	for (auto e : mapper) {
		//! [BinningFeatProcessor::init]
		if (e.first == "name") feature_name = e.second;
		else if (e.first == "missing_value") missing_value = med_stof(e.second);
		else if (e.first == "missing_target_val") missing_target_val = med_stof(e.second);
		else if (e.first == "bin_format") bin_format = e.second;
		else if (e.first == "one_hot") one_hot = med_stoi(e.second) > 0;
		else if (e.first == "remove_origin") remove_origin = med_stoi(e.second) > 0;
		else if (e.first == "keep_original_val") keep_original_val = med_stoi(e.second) > 0;
		else if (e.first == "bin_sett") bin_sett.init_from_string(e.second);
		else if (e.first != "names" && e.first != "fp_type" && e.first != "tag")
			MTHROW_AND_ERR("Unknown parameter \'%s\' for BinningFeatProcessor\n", e.first.c_str());
		//! [BinningFeatProcessor::init]
	}

	if (keep_original_val && !one_hot)
		MTHROW_AND_ERR("ERROR - BinningFeatProcessor::init Error Can't turn on keep_original_val when one_hot mode is off\n");

	return 0;
}

int BinningFeatProcessor::Learn(MedFeatures& features, unordered_set<int>& ids) {
	resolved_feature_name = resolve_feature_name(features, feature_name);
	vector<float> y(features.samples.size());
	for (size_t i = 0; i < y.size(); ++i)
		y[i] = features.samples[i].outcome;
	bin_sett.load_bin_settings(features.data.at(resolved_feature_name), y);
	return 0;
}

char buffer[1000];
string BinningFeatProcessor::get_bin_name(float num) const {
	snprintf(buffer, sizeof(buffer), bin_format.c_str(), num);
	return resolved_feature_name + ".BINNED_" + string(buffer);
}

int BinningFeatProcessor::_apply(MedFeatures& features, unordered_set<int>& ids) {
	resolved_feature_name = resolve_feature_name(features, feature_name);

	vector<float>& data = features.data[resolved_feature_name];
	vector<float> *data_p = &data;
	vector<float> data_cp;
	if (keep_original_val && one_hot) { //only optional when using one hot
		data_cp.resize(features.samples.size());
		data_p = &data_cp;
	}
	for (unsigned int i = 0; i < features.samples.size(); i++) {
		if ((*data_p)[i] == missing_value)
			(*data_p)[i] = missing_target_val;
		else
			(*data_p)[i] = bin_sett.normalize(data[i]);
	}

	//Split to one hot based on data uniq vals:
	if (one_hot) {
		unordered_set<string> new_features;
		for (size_t i = 0; i < bin_sett.bin_repr_vals.size(); ++i)
			new_features.insert(get_bin_name(bin_sett.bin_repr_vals[i]));
		//Add all features:
#pragma omp critical 
		{
			for (const string &full_name : new_features)
			{
				features.data[full_name].resize(features.samples.size());
				features.attributes[full_name].normalized = true;
			}
		}

		//calculate feature group for each row
		vector<string> parsed_names(data_p->size());
		for (size_t i = 0; i < data_p->size(); ++i)
			parsed_names[i] = get_bin_name((*data_p)[i]);

		//Write values in binary bits:
		for (size_t i = 0; i < data_p->size(); ++i)
			features.data.at(parsed_names[i])[i] = 1;
	}

#pragma omp critical
	if (remove_origin) {
		features.data.erase(resolved_feature_name);
		features.attributes.erase(resolved_feature_name);
	}
	return 0;
}