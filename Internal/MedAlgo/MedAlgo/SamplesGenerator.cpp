#include "SamplesGenerator.h"

#include <omp.h>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL LOG_DEF_LEVEL

template<typename T> void SamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values) {
	MTHROW_AND_ERR("SamplesGenerator<T>::Not Implemented\n");
}
template<typename T> void SamplesGenerator<T>::get_samples(vector<vector<T>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values) {
	MTHROW_AND_ERR("SamplesGenerator<T>::Not Implemented\n");
}
template<typename T> void SamplesGenerator<T>::get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const {
	MTHROW_AND_ERR("SamplesGenerator<T>::Not Implemented\n");
}
template<typename T> void SamplesGenerator<T>::get_samples(vector<vector<T>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values, mt19937 &rnd_gen) const {
	MTHROW_AND_ERR("SamplesGenerator<T>::Not Implemented\n");
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

template<typename T> void GibbsSamplesGenerator<T>::learn(const map<string, vector<T>> &data) {
	_gibbs->learn_gibbs(data);
}

template<typename T> void *SamplesGenerator<T>::new_polymorphic(string derived_name) {
	CONDITIONAL_NEW_CLASS(derived_name, GibbsSamplesGenerator<T>);
	CONDITIONAL_NEW_CLASS(derived_name, MaskedGAN<float>);
	CONDITIONAL_NEW_CLASS(derived_name, MissingsSamplesGenerator<float>);
	if (boost::starts_with(derived_name, "GibbsSamplesGenerator")) return new GibbsSamplesGenerator<T>;
	if (boost::starts_with(derived_name, "MaskedGAN")) return new MaskedGAN<float>;
	if (boost::starts_with(derived_name, "MissingsSamplesGenerator")) return new MissingsSamplesGenerator<float>;
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

template<typename T> void GibbsSamplesGenerator<T>::get_samples(vector<vector<T>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values) {
	MTHROW_AND_ERR("Error no supported in Gibbs");
}
template<typename T> void GibbsSamplesGenerator<T>::get_samples(vector<vector<T>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values, mt19937 &rnd_gen) const {
	MTHROW_AND_ERR("Error no supported in Gibbs");
}

template class SamplesGenerator<float>;
template class SamplesGenerator<double>;
template class GibbsSamplesGenerator<float>;
template class GibbsSamplesGenerator<double>;

// Masked GAN
template<typename T> MaskedGAN<T>::MaskedGAN()
	: SamplesGenerator<T>(true) {
	random_device rd;
	_gen = mt19937(rd());
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

template<> void MaskedGAN<float>::get_samples(vector<vector<float>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<float>> &mask_values, mt19937 &rnd_gen) const {
	// Sanity
	if (mask.size() != mask_values.size())
		MTHROW_AND_ERR("size mismatch between mask (%d samples) and mask_values (%d samples)\n", (int)mask.size(), (int)mask_values.size());

	if (mask.empty())
		return;

	// Sample
	int nFtrs = (int)mask[0].size();

	data.resize(sample_per_row * mask.size());
	vector<float> input(2 * nFtrs);

	uniform_real_distribution<> real_dist(-1.0, 1.0);
	int index = 0;
	for (int i = 0; i < mask.size(); i++) {
		for (int j = 0; j < sample_per_row; j++) {
			// Generate input Z
			for (int k = 0; k < nFtrs; k++) {
				if (mask[i][k]) {
					input[k] = mask_values[i][k];
					input[k + nFtrs] = 1.0;
				}
				else {
					input[k] = (float)real_dist(rnd_gen);
					input[k + nFtrs] = 0.0;
				}
			}

			// Apply generator
			generator.apply(input, data[index]);

			// Mask
			for (int k = 0; k < nFtrs; k++) {
				if (mask[i][k])
					data[index][k] = mask_values[i][k];
				else if (!mg_params.keep_original_values)
					data[i][k] = round_to_allowed_values(data[i][k], allowed_values[k]);
			}

			index++;
		}
	}
}

template<> void MaskedGAN<float>::get_samples(vector<vector<float>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<float>> &mask_values) {
	set_params(params);
	get_samples(data, sample_per_row, params, mask, mask_values, _gen);
}

template<> void MaskedGAN<float>::get_samples_from_Z(vector<vector<float>> &data, void *params, const vector<vector<bool>> &mask, const vector<vector<float>> &Z) {

	set_params(params);

	// Sanity
	if (mask.size() != Z.size())
		MTHROW_AND_ERR("size mismatch between mask (%d samples) and Z (%d samples)\n", (int)mask.size(), (int)Z.size());

	if (mask.empty())
		return;

	// Sample
	int nFtrs = (int)mask[0].size();

	data.resize(mask.size());
	vector<float> input(2 * nFtrs);

	for (int i = 0; i < mask.size(); i++) {
		// Generate input Z
		for (int k = 0; k < nFtrs; k++) {
			input[k] = Z[i][k];
			if (mask[i][k])
				input[k + nFtrs] = 1.0;
			else
				input[k + nFtrs] = 0.0;
		}

		// Apply generator
		generator.apply(input, data[i]);

		// Mask + round to allowed values
		for (int k = 0; k < nFtrs; k++) {
			if (mask[i][k])
				data[i][k] = Z[i][k];
			else if (!mg_params.keep_original_values)
				data[i][k] = round_to_allowed_values(data[i][k], allowed_values[k]);
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

template<typename T> void MissingsSamplesGenerator<T>::get_samples(vector<vector<T>> &data, int sample_per_row, void *params,
	const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values, mt19937 &rnd_gen) const {

	data.resize(mask.size());
	for (size_t i = 0; i < mask.size(); ++i)
	{
		data[i].resize(mask[i].size());
		for (size_t j = 0; j < names.size(); ++j)
		{
			if (mask[i][j])
				data[i][j] = mask_values[i][j];
			else
				data[i][j] = missing_value;
		}
	}
}

template<typename T> void MissingsSamplesGenerator<T>::get_samples(vector<vector<T>> &data, int sample_per_row, void *params,
	const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values) {
	mt19937 rnd_not_used;
	get_samples(data, sample_per_row, params, mask, mask_values, rnd_not_used);
}

template<typename T> MissingsSamplesGenerator<T>::MissingsSamplesGenerator(float miss_valu)
	: SamplesGenerator<T>(false) {
	missing_value = miss_valu;
}
template<typename T> MissingsSamplesGenerator<T>::MissingsSamplesGenerator()
	: SamplesGenerator<T>(false) {}

template<typename T> void MissingsSamplesGenerator<T>::learn(const map<string, vector<T>> &data) {
	names.reserve(data.size());
	for (auto it = data.begin(); it != data.end(); ++it)
		names.push_back(it->first);
}

template class MissingsSamplesGenerator<float>;