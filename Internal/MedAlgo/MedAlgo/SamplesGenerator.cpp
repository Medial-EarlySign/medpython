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
	for (const auto it : mapper)
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