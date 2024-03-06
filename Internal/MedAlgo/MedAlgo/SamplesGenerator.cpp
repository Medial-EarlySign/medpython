#include "SamplesGenerator.h"
#include "medial_utilities/medial_utilities/globalRNG.h"
#include <random>

#include <omp.h>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL LOG_DEF_LEVEL

string GeneratorType_toStr(GeneratorType type) {
	switch (type)
	{
	case GIBBS:
		return "GIBBS";
	case GAN:
		return "GAN";
	case MISSING:
		return "MISSING";
	case RANDOM_DIST:
		return "RANDOM_DIST";
	case UNIVARIATE_DIST:
		return "UNIVARIATE_DIST";
	default:
		MTHROW_AND_ERR("Unknown type %d\n", type);
	}
}
GeneratorType GeneratorType_fromStr(const string &type) {
	string tp = boost::to_upper_copy(type);
	if (tp == "GAN")
		return GeneratorType::GAN;
	else if (tp == "GIBBS")
		return GeneratorType::GIBBS;
	else if (tp == "MISSING")
		return GeneratorType::MISSING;
	else if (tp == "RANDOM_DIST")
		return GeneratorType::RANDOM_DIST;
	else if (tp == "UNIVARIATE_DIST")
		return GeneratorType::UNIVARIATE_DIST;
	else
		MTHROW_AND_ERR("Unknown type %s\n", type.c_str());
}

template<typename T> void SamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values) {
	MTHROW_AND_ERR("SamplesGenerator<T>::Not Implemented\n");
}
template<typename T> void SamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values) {
	MTHROW_AND_ERR("SamplesGenerator<T>::Not Implemented\n");
}
template<typename T> void SamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const {
	MTHROW_AND_ERR("SamplesGenerator<T>::Not Implemented\n");
}
template<typename T> void SamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const {
	MTHROW_AND_ERR("SamplesGenerator<T>::Not Implemented\n");
}
template<typename T> void SamplesGenerator<T>::learn(const map<string, vector<T>> &data) {
	vector<string> all(data.size());
	int ii = 0;
	for (const auto &it : data)
	{
		all[ii] = it.first;
		++ii;
	}
	learn(data, all, false);
}

template<typename T> SamplesGenerator<T>::SamplesGenerator() {
	use_vector_api = true;
}

template<typename T> void SamplesGenerator<T>::pre_serialization() {}
template<typename T> void SamplesGenerator<T>::post_deserialization() {}
template<typename T> void GibbsSamplesGenerator<T>::pre_serialization() {}
template<typename T> void GibbsSamplesGenerator<T>::post_deserialization() {}
template<typename T> void MaskedGAN<T>::pre_serialization() {}
template<typename T> void MaskedGAN<T>::post_deserialization() {}
template<typename T> void MissingsSamplesGenerator<T>::pre_serialization() {}
template<typename T> void MissingsSamplesGenerator<T>::post_deserialization() {}
template<typename T> void RandomSamplesGenerator<T>::pre_serialization() {}
template<typename T> void RandomSamplesGenerator<T>::post_deserialization() {}
template<typename T> void UnivariateSamplesGenerator<T>::pre_serialization() {}
template<typename T> void UnivariateSamplesGenerator<T>::post_deserialization() {}

template<typename T> GibbsSamplesGenerator<T>::GibbsSamplesGenerator() : SamplesGenerator<T>(false) {
	_gibbs = NULL;
	_do_parallel = true;
	no_need_to_clear_mem = false;
}

template<typename T> GibbsSamplesGenerator<T>::~GibbsSamplesGenerator() {
	if (_gibbs != NULL && !no_need_to_clear_mem) {
		delete _gibbs;
		_gibbs = NULL;
	}
}

template<typename T> SamplesGenerator<T>::SamplesGenerator(bool _use_vector_api) {
	use_vector_api = _use_vector_api;
}

// Gibbs samples generator
template<typename T> GibbsSamplesGenerator<T>::GibbsSamplesGenerator(GibbsSampler<T> &gibbs, bool do_parallel, bool no_need_clear_mem)
	: SamplesGenerator<T>(false) {
	_gibbs = &gibbs;
	_do_parallel = do_parallel;
	no_need_to_clear_mem = no_need_clear_mem;
}

template<typename T> void GibbsSamplesGenerator<T>::learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing) {
	_gibbs->learn_gibbs(data, learn_features, skip_missing);
}

template<typename T> void *SamplesGenerator<T>::new_polymorphic(string derived_name) {
	CONDITIONAL_NEW_CLASS(derived_name, GibbsSamplesGenerator<T>);
	CONDITIONAL_NEW_CLASS(derived_name, MaskedGAN<float>);
	CONDITIONAL_NEW_CLASS(derived_name, MissingsSamplesGenerator<float>);
	CONDITIONAL_NEW_CLASS(derived_name, RandomSamplesGenerator<float>);
	if (boost::starts_with(derived_name, "GibbsSamplesGenerator")) return new GibbsSamplesGenerator<T>;
	if (boost::starts_with(derived_name, "MaskedGAN")) return new MaskedGAN<float>;
	if (boost::starts_with(derived_name, "MissingsSamplesGenerator")) return new MissingsSamplesGenerator<float>;
	if (boost::starts_with(derived_name, "RandomSamplesGenerator")) return new RandomSamplesGenerator<float>;
	if (boost::starts_with(derived_name, "UnivariateSamplesGenerator")) return new UnivariateSamplesGenerator<float>;
	MTHROW_AND_ERR("SamplesGenerator<T>::new_polymorphic:: Unsupported object %s\n", derived_name.c_str());
	return NULL;
}

template<typename T> void GibbsSamplesGenerator<T>::prepare(void *params) {
	_gibbs->prepare_predictors();
}

template<typename T> void GibbsSamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const {
	GibbsSamplingParams *sampling_params = (GibbsSamplingParams *)params;
	_gibbs->get_samples(data, *sampling_params, rnd_gen, &mask, &mask_values);
}

template<typename T> void GibbsSamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values) {
	GibbsSamplingParams *sampling_params = (GibbsSamplingParams *)params;
	if (_do_parallel && sampling_params->samples_count >= 10)
		_gibbs->get_parallel_samples(data, *sampling_params, &mask, &mask_values);
	else
		_gibbs->get_samples(data, *sampling_params, &mask, &mask_values);
}

template<typename T> void GibbsSamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values) {
	MTHROW_AND_ERR("Error no supported in Gibbs");
}
template<typename T> void GibbsSamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const {
	MTHROW_AND_ERR("Error no supported in Gibbs");
}

template class SamplesGenerator<float>;
template class SamplesGenerator<double>;
template class GibbsSamplesGenerator<float>;
template class GibbsSamplesGenerator<double>;

int MaskedGANParams::init(map<string, string> &mapper) {
	for (const auto &it : mapper)
	{
		if (it.first == "keep_original_values")
			keep_original_values = med_stoi(it.second) > 0;
		else
			MTHROW_AND_ERR("Error in MaskedGANParams::init Unsupported argument %s\n", it.first.c_str());
	}
	return 0;
}

// Masked GAN
template<typename T> MaskedGAN<T>::MaskedGAN()
	: SamplesGenerator<T>(true) {
	_gen = mt19937(globalRNG::rand());
	norm_by_by_file = false;
}

template<typename T> void MaskedGAN<T>::get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values) {
	MTHROW_AND_ERR("Error: Mode not supported for MaskedGAN");
}
template<typename T> void MaskedGAN<T>::get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const {
	MTHROW_AND_ERR("Error: Mode not supported for MaskedGAN");
}

template<> void MaskedGAN<float>::prepare(void *params) {
	set_params(params);
}

template<> void MaskedGAN<float>::get_samples(MedMat<float> &data, int sample_per_row, void *params, const vector<vector<bool>> &masks, const MedMat<float> &mask_values, mt19937 &rnd_gen) const {

	// Sanity
	if ((int)masks.size() != mask_values.nrows)
		MTHROW_AND_ERR("size mismatch between mask (%d samples) and mask_values (%lld samples)\n", (int)masks.size(), mask_values.nrows);

	if (mask_values.nrows == 0)
		return;

	// Sample
	int nSamples = mask_values.nrows;
	int nFtrs = mask_values.ncols;

	data.resize(sample_per_row * nSamples, nFtrs);
	MedMat<float> input(sample_per_row * nSamples, 3 * nFtrs);

	normal_distribution<> norm_dist(0, 1);
	int index = 0;
	for (int i = 0; i < nSamples; i++) {
		for (int j = 0; j < sample_per_row; j++) {
			// Generate input Z + ZX + I
			for (int k = 0; k < nFtrs; k++) {
				if (masks[i][k]) {
					input(index, k) = (float)norm_dist(rnd_gen);
					input(index, k + nFtrs) = mask_values(i, k);
					input(index, k + 2 * nFtrs) = 1.0;
					if (norm_by_by_file) {
						//normalize values
						input(index, k + nFtrs) = input(index, k + nFtrs) - mean_feature_vals[k];
						if (std_feature_vals[k] > 0)
							input(index, k + nFtrs) /= std_feature_vals[k];
					}
				}
				else {
					input(index, k) = (float)norm_dist(rnd_gen);
					input(index, k + nFtrs) = 0.0;
					input(index, k + 2 * nFtrs) = 0.0;
				}
			}

			index++;
		}
	}

	// Apply generator
	generator.apply(input, data);

	// Mask
	index = 0;
	for (int i = 0; i < nSamples; i++) {
		for (int j = 0; j < sample_per_row; j++) {
			for (int k = 0; k < nFtrs; k++) {
				if (masks[i][k])
					data(index, k) = mask_values(i, k);
				else {
					//unorm if needed:
					if (norm_by_by_file) {
						//unorm outData:
						if (std_feature_vals[k] > 0)
							data(index, k) = data(index, k) * std_feature_vals[k] + mean_feature_vals[k];
					}
					if (!mg_params.keep_original_values)
						data(index, k) = round_to_allowed_values(data(i, k), allowed_values[k]);
				}
			}
			index++;
		}
	}
}

template<> void MaskedGAN<float>::get_samples(MedMat<float> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<float> &mask_values) {
	set_params(params);
	get_samples(data, sample_per_row, params, mask, mask_values, _gen);
}

template<> void MaskedGAN<float>::get_samples_from_Z(MedMat<float> &data, void *params, const vector<vector<bool>> &masks, const MedMat<float> &mask_values, const MedMat<float> &Z) {

	set_params(params);

	// Sanity
	if ((int)masks.size() != Z.nrows || Z.nrows != mask_values.nrows)
		MTHROW_AND_ERR("size mismatch between mask (%d samples), Z (%lld samples) and mask_values (%lld samples)\n", (int)masks.size(), Z.nrows, mask_values.nrows);

	if (mask_values.nrows == 0)
		return;

	// Sample
	int nSamples = mask_values.nrows;
	int nFtrs = mask_values.ncols;

	data.resize(nSamples, nFtrs);
	MedMat<float> input(nSamples, 3 * nFtrs);

	normal_distribution<> norm_dist;
	for (int i = 0; i < nSamples; i++) {
		// Generate input Z + ZX + I
		for (int k = 0; k < nFtrs; k++) {
			if (masks[i][k]) {
				input(i, k) = Z(i, k);
				input(i, k + nFtrs) = mask_values(i, k);
				if (norm_by_by_file) {
					//normalize values
					input(i, k + nFtrs) = input(i, k + nFtrs) - mean_feature_vals[k];
					if (std_feature_vals[k] > 0)
						input(i, k + nFtrs) /= std_feature_vals[k];
				}
				input(i, k + 2 * nFtrs) = 1.0;
			}
			else {
				input(i, k) = Z(i, k);
				input(i, k + nFtrs) = 0.0;
				input(i, k + 2 * nFtrs) = 0.0;
			}
		}
	}

	// Apply generator
	generator.apply(input, data);

	// Mask
	for (int i = 0; i < nSamples; i++) {
		for (int k = 0; k < nFtrs; k++) {
			if (masks[i][k])
				data(i, k) = mask_values(i, k);
			else {
				//unorm if needed:
				if (norm_by_by_file) {
					//unorm outData:
					if (std_feature_vals[k] > 0)
						data(i, k) = data(i, k) * std_feature_vals[k] + mean_feature_vals[k];
				}
				if (!mg_params.keep_original_values)
					data(i, k) = round_to_allowed_values(data(i, k), allowed_values[k]);
			}
		}
	}

}

template<> void MaskedGAN<float>::read_from_text_file(const string& file_name) {

	// Read geneator (ApplyKeras) object
	generator.init_from_text_file(file_name);

	// Read allowed values
	string allowed_values_file_name = file_name + ".allowed_values";

	ifstream inf(allowed_values_file_name);
	if (!inf.is_open())
		MTHROW_AND_ERR("Cannot opend allowed-values file \'%s\' for reading\n", allowed_values_file_name.c_str());

	string curr_line;
	vector<string> fields;
	while (getline(inf, curr_line)) {
		boost::split(fields, curr_line, boost::is_any_of(","));
		allowed_values.push_back(vector<float>(fields.size()));
		for (unsigned int i = 0; i < fields.size(); i++)
			allowed_values.back()[i] = stof(fields[i]);
	}

	//check for norm file:
	string norm_file_name = file_name + ".norm_vector";
	if (file_exists(norm_file_name)) {
		ifstream f_norm(norm_file_name);
		norm_by_by_file = true;
		string line;
		getline(f_norm, line);
		vector<string> tokens;
		boost::split(tokens, line, boost::is_any_of(","));
		mean_feature_vals.resize(tokens.size());
		for (size_t i = 0; i < tokens.size(); ++i)
			mean_feature_vals[i] = med_stof(tokens[i]);
		getline(f_norm, line);
		tokens.clear();
		boost::split(tokens, line, boost::is_any_of(","));
		std_feature_vals.resize(tokens.size());
		for (size_t i = 0; i < tokens.size(); ++i)
			std_feature_vals[i] = med_stof(tokens[i]);
		f_norm.close();
	}
}

template<typename T> T MaskedGAN<T>::round_to_allowed_values(T in_value, const vector<T>& curr_allowed_values) const {

	// Perform binary search
	unsigned int start = 0;
	unsigned int end = (unsigned int)(curr_allowed_values.size() - 1);
	while (end > start + 1) {
		int mid = (start + end) / 2;
		if (in_value > curr_allowed_values[mid])
			start = mid;
		else
			end = mid;
	}

	if (abs(in_value - curr_allowed_values[end]) < abs(in_value - curr_allowed_values[start]))
		return curr_allowed_values[end];
	else
		return curr_allowed_values[start];

}

template<typename T> void MaskedGAN<T>::set_params(void *params) {

	if (params != NULL) {
		MaskedGANParams *_params = (MaskedGANParams *)params;
		mg_params.keep_original_values = _params->keep_original_values;
	}
}

template class MaskedGAN<float>;

template<typename T> void MissingsSamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params,
	const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const {
	for (size_t i = 0; i < names.size(); ++i)
	{
		if (mask[i])
			data[names[i]].push_back(mask_values[i]);
		else
			data[names[i]].push_back(missing_value);
	}
}
template<typename T> void MissingsSamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params,
	const vector<bool> &mask, const vector<T> &mask_values) {
	mt19937 rnd_not_used;
	get_samples(data, params, mask, mask_values, rnd_not_used);
}

template<typename T> void MissingsSamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params,
	const vector<vector<bool>> &masks, const MedMat<T> &mask_values, mt19937 &rnd_gen) const {

	if (!masks.empty()) {
		data.resize((int)masks.size(), (int)masks[0].size());
		for (int i = 0; i < masks.size(); ++i)
		{
			for (int j = 0; j < names.size(); ++j)
			{
				if (masks[i][j])
					data(i, j) = mask_values(i, j);
				else
					data(i, j) = missing_value;
			}
		}
	}
	else
		data.clear();
}

template<typename T> void MissingsSamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params,
	const vector<vector<bool>> &mask, const MedMat<T> &mask_values) {
	mt19937 rnd_not_used;
	get_samples(data, sample_per_row, params, mask, mask_values, rnd_not_used);
}

template<typename T> MissingsSamplesGenerator<T>::MissingsSamplesGenerator(float miss_valu)
	: SamplesGenerator<T>(false) {
	missing_value = miss_valu;
}
template<typename T> MissingsSamplesGenerator<T>::MissingsSamplesGenerator()
	: SamplesGenerator<T>(false) {}

template<typename T> void MissingsSamplesGenerator<T>::learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing) {
	names.reserve(data.size());
	for (auto it = data.begin(); it != data.end(); ++it)
		names.push_back(it->first);
}

template<typename T> RandomSamplesGenerator<T>::RandomSamplesGenerator() : SamplesGenerator<T>(false) {
	mean_value = (T)0;
	std_value = (T)5;
}

template<typename T> RandomSamplesGenerator<T>::RandomSamplesGenerator(T mean_val, T std_val) : SamplesGenerator<T>(false) {
	mean_value = mean_val;
	std_value = std_val;
}

template<typename T> void RandomSamplesGenerator<T>::learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing) {
	names.reserve(data.size());
	for (auto it = data.begin(); it != data.end(); ++it)
		names.push_back(it->first);
}

template<typename T> void RandomSamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params,
	const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const {
	normal_distribution<> norm_gen(mean_value, std_value);
	int smp_count = *(int *)params;
	for (size_t s = 0; s < smp_count; ++s)
	{
		for (size_t i = 0; i < names.size(); ++i)
		{
			if (mask[i])
				data[names[i]].push_back(mask_values[i]);
			else
				data[names[i]].push_back(norm_gen(rnd_gen));
		}
	}
}
template<typename T> void RandomSamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params,
	const vector<bool> &mask, const vector<T> &mask_values) {
	random_device rd;
	mt19937 rnd_gen(rd());
	get_samples(data, params, mask, mask_values, rnd_gen);
}

template<typename T> void RandomSamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params,
	const vector<vector<bool>> &masks, const MedMat<T> &mask_values, mt19937 &rnd_gen) const {
	normal_distribution<> norm_gen(mean_value, std_value);
	if (!masks.empty()) {
		data.resize((int)masks.size(), (int)masks[0].size());
		for (int i = 0; i < masks.size(); ++i)
		{
			for (int j = 0; j < names.size(); ++j)
			{
				if (masks[i][j])
					data(i, j) = mask_values(i, j);
				else
					data(i, j) = norm_gen(rnd_gen);
			}
		}
	}
	else
		data.clear();
}

template<typename T> void RandomSamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params,
	const vector<vector<bool>> &mask, const MedMat<T> &mask_values) {
	random_device rd;
	mt19937 rnd_gen(rd());
	get_samples(data, sample_per_row, params, mask, mask_values, rnd_gen);
}

template class RandomSamplesGenerator<float>;

template class MissingsSamplesGenerator<float>;

template<typename T> UnivariateSamplesGenerator<T>::UnivariateSamplesGenerator() : SamplesGenerator<T>(false) {
	min_samples = 50;
}

template<typename T> void UnivariateSamplesGenerator<T>::learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing) {
	if (data.empty())
		return;
	for (const auto &it : data)
		names.push_back(it.first);
	vector<string> final_feats = learn_features;
	if (learn_features.empty())
		final_feats = names;
	MLOG("INFO :: UnivariateSamplesGenerator Learn started (%zu)\n", final_feats.size());

	strata_settings.getFactors();
	vector<const vector<float> *> strataData(strata_settings.nStratas());
	for (int j = 0; j < strata_settings.nStratas(); j++) {
		string resolved_strata_name = names[find_in_feature_names(names, strata_settings.stratas[j].name)];
		strataData[j] = &data.at(resolved_strata_name);
	}
	if (strata_settings.nStratas() > 0)
		MLOG("INFO:: UnivariateSamplesGenerator::learn - has %zu stratas with %d bins\n",
			strata_settings.nStratas(), strata_settings.nValues());
	//already sorted because map
	feature_values.resize(names.size());
	feature_val_probs.resize(names.size());
	unordered_set<string> learn_set(final_feats.begin(), final_feats.end());
	for (int i = 0; i < names.size(); ++i)
	{
		if (learn_set.find(names[i]) == learn_set.end())
			continue;
		const vector<T> &v = data.at(names[i]);
		map<T, double> val_to_prob;

		double tot_cnt = 0; //count non missings or all
		for (T val : v)
		{
			if (skip_missing && val == missing_value)
				continue;
			//find val index in val_to_val
			++val_to_prob[val];
			++tot_cnt;
		}
		//process v into val_to_prob:
		if (tot_cnt > 0)
			for (auto &it : val_to_prob)
				val_to_prob[it.first] /= tot_cnt;

		vector<T> &val_to_val = feature_values[i];
		vector<double> &val_to_p = feature_val_probs[i];
		//update val_to_val,val_to_prob from val_to_prob which is already ordered (map)
		double cumsum = 0;
		for (auto &it : val_to_prob) {
			val_to_val.push_back(it.first);
			cumsum += it.second;
			val_to_p.push_back(cumsum);
		}

	}

	//calc for stratas: strata_sizes, strata_feature_val_agg:
	int num_of_rows = (int)data.begin()->second.size();
	vector<int> strata_ind(num_of_rows);
	strata_sizes.resize(strata_settings.nValues());
	strata_feature_val_agg_val.resize(strata_sizes.size());
	strata_feature_val_agg_prob.resize(strata_sizes.size());
	vector<vector<int>> strata_to_indexes(strata_sizes.size());
	for (int i = 0; i < num_of_rows; ++i)
	{
		strata_ind[i] = strata_settings.getIndex(missing_value, strataData, i);
		++strata_sizes[strata_ind[i]];
		strata_to_indexes[strata_ind[i]].push_back(i);
	}
	//update strata_feature_val_agg:
	int skip_cnt = 0;
	for (int i = 0; i < strata_sizes.size(); ++i)
	{
		vector<vector<T>> &feat_val_val = strata_feature_val_agg_val[i];
		vector<vector<double>> &feat_val_prb = strata_feature_val_agg_prob[i];
		if (strata_sizes[i] >= min_samples) {
			feat_val_val.resize(names.size());
			feat_val_prb.resize(names.size());
			for (int j = 0; j < names.size(); ++j)
			{
				const string &feat = names[j];
				if (learn_set.find(feat) == learn_set.end())
					continue;
				const vector<T> &v = data.at(feat);
				map<T, double> val_to_prob;
				//map<float, double> &val_to_prob = strata_feature_val_agg[i][feat];
				double tot_cnt = 0; //count non missings or all
				const vector<int> &relevant_idx = strata_to_indexes[i];
				for (int idx : relevant_idx) //go over relevant indexes only
				{
					T val = v[idx];
					if (skip_missing && val == missing_value)
						continue;
					++val_to_prob[val];
					++tot_cnt;
				}
				//process v into val_to_prob:
				if (tot_cnt > 0)
					for (auto &it : val_to_prob)
						val_to_prob[it.first] /= tot_cnt;

				vector<T> &val_vec = feat_val_val[j];
				vector<double> &prob_vec = feat_val_prb[j];
				//update val_vec,prob_vec from val_to_prob 
				double cumsum = 0;
				for (auto &it : val_to_prob) {
					val_vec.push_back(it.first);
					cumsum += it.second;
					prob_vec.push_back(cumsum);
				}
			}

		}
		else
			++skip_cnt;
	}
	if (skip_cnt > 0)
		MWARN("Has %d skipped strats with few samples\n", skip_cnt);
}

template<typename T> T UnivariateSamplesGenerator<T>::find_pos(const vector<T> &v, const vector<double> &cumsum, double p) const {
	int pos = medial::process::binary_search_position(cumsum, p);
	//int pos = medial::process::binary_search_position(cumsum.data(), cumsum.data() + cumsum.size() - 1, p);
	if (pos >= v.size())
		pos = (int)v.size() - 1;
	return v[pos];
	
}

template<typename T> void UnivariateSamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params,
	const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const {
	int smp_count = 1;
	if (params != NULL)
		smp_count = *(int *)params;
	//strata_settings.getFactors();
	vector<const vector<float> *> strataData(strata_settings.nStratas());
	vector<vector<float>> data_s(strataData.size());
	for (int j = 0; j < strata_settings.nStratas(); j++) {
		int starta_feat_idx = find_in_feature_names(names, strata_settings.stratas[j].name);
		data_s[j].push_back(mask_values[starta_feat_idx]);
		strataData[j] = &data_s[j];
	}

	uniform_real_distribution<> rnd_prob(0, 1);
	for (size_t s = 0; s < smp_count; ++s)
	{
		//generate sample with all features:
		for (size_t i = 0; i < names.size(); ++i)
		{
			if (mask[i])
				data[names[i]].push_back(mask_values[i]);
			else {
				//get strata index for this one sample:
				int strata_index = strata_settings.getIndex(missing_value, strataData, 0);

				//const map<T, double> *feat_prob = &feature_val_agg.at(names[i]);
				const vector<T> *feat_vals = &feature_values[i];
				const vector<double> *feat_prbs = &feature_val_probs[i];

				if (strata_sizes[strata_index] >= min_samples && !strata_feature_val_agg_val[strata_index].empty()) {
					feat_vals = &strata_feature_val_agg_val[strata_index][i];
					feat_prbs = &strata_feature_val_agg_prob[strata_index][i];
					//feat_prob = &strata_feature_val_agg[strata_index].at(names[i]);
				}

				double p = rnd_prob(rnd_gen);
				T val = find_pos(*feat_vals, *feat_prbs, p);

				data[names[i]].push_back(val);
			}
		}
	}
}
template<typename T> void UnivariateSamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params,
	const vector<bool> &mask, const vector<T> &mask_values) {
	random_device rd;
	mt19937 rnd_gen(rd());
	get_samples(data, params, mask, mask_values, rnd_gen);
}

template<typename T> void UnivariateSamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row,
	void *params, const vector<vector<bool>> &masks, const MedMat<T> &mask_values, mt19937 &rnd_gen) const {
	if (!masks.empty()) {
		//strata_settings.getFactors();
		vector<const vector<float> *> strataData(strata_settings.nStratas());
		vector<vector<float>> data_s(strataData.size());
		for (int j = 0; j < strata_settings.nStratas(); j++) {
			int starta_feat_idx = find_in_feature_names(names, strata_settings.stratas[j].name);
			data_s[j].resize(mask_values.nrows);
			for (int i = 0; i < mask_values.nrows; ++i)
				data_s[j][i] = mask_values(i, starta_feat_idx);
			strataData[j] = &data_s[j];
		}

		data.resize((int)masks.size(), (int)masks[0].size());
		uniform_real_distribution<> rnd_prob(0, 1);
		for (int i = 0; i < masks.size(); ++i)
		{
			int strata_index = strata_settings.getIndex(missing_value, strataData, i);

			for (int j = 0; j < names.size(); ++j)
			{
				if (masks[i][j])
					data(i, j) = mask_values(i, j);
				else {


					//const map<T, double> *feat_prob = &feature_val_agg.at(names[j]);
					const vector<T> *feat_vals = &feature_values[j];
					const vector<double> *feat_prbs = &feature_val_probs[j];
					if (strata_sizes[strata_index] >= min_samples && !strata_feature_val_agg_val[strata_index].empty()) {
						feat_vals = &strata_feature_val_agg_val[strata_index][j];
						feat_prbs = &strata_feature_val_agg_prob[strata_index][j];
						//feat_prob = &strata_feature_val_agg[strata_index].at(names[i]);
					}

					double p = rnd_prob(rnd_gen);
					T val = find_pos(*feat_vals, *feat_prbs, p);
					data(i, j) = val;
				}
			}
		}
	}
	else
		data.clear();
}
template<typename T> void UnivariateSamplesGenerator<T>::get_samples(MedMat<T> &data, int sample_per_row, void *params,
	const vector<vector<bool>> &mask, const MedMat<T> &mask_values) {
	random_device rd;
	mt19937 rnd_gen(rd());
	get_samples(data, sample_per_row, params, mask, mask_values, rnd_gen);
}

void addStrata(string& init_string, featureSetStrata &imputerStrata) {

	vector<string> fields;
	boost::split(fields, init_string, boost::is_any_of(","));

	if (fields.size() != 4)
		MLOG("Cannot initialize strata from \'%s\'. Ignoring\n", init_string.c_str());
	else
		imputerStrata.stratas.push_back(featureStrata(fields[0], stof(fields[3]), stof(fields[1]), stof(fields[2])));
}

template<typename T> int UnivariateSamplesGenerator<T>::init(map<string, string>& mapper) {
	vector<string> strata;
	for (auto &it : mapper)
	{
		if (it.first == "strata" || it.first == "strata_settings") {
			boost::split(strata, it.second, boost::is_any_of(":"));
			for (string& stratum : strata) addStrata(stratum, strata_settings);
		}
		else if (it.first == "min_samples")
			min_samples = med_stoi(it.second);
		else if (it.first == "missing_value")
			missing_value = med_stof(it.second);
		else
			MTHROW_AND_ERR("Error UnivariateSamplesGenerator::init - not found arg %s\n",
				it.first.c_str());
	}
	return 0;
}

template class UnivariateSamplesGenerator<float>;