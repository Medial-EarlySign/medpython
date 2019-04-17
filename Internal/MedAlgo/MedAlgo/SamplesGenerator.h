// Classes for generating samples from partial information

#ifndef __SAMPLES_GENERATOR_H__
#define __SAMPLES_GENERATOR_H__

#include <vector>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedStat/MedStat/GibbsSampler.h>
#include <MedEmbed/MedEmbed/ApplyKeras.h>
#include <random>

using namespace std;

/**
* Abstract Random Samples generator
*/
template<typename T> class SamplesGenerator : public SerializableObject {
protected:
	SamplesGenerator(bool _use_vector_api);
public:
	bool use_vector_api = true; ///< In gibbs it's faster to use map<string, float> api

	SamplesGenerator();

	/// <summary>
	/// learn of sample generator
	/// </summary>
	virtual void learn(const map<string, vector<T>> &data) {};

	/// <summary>
	/// apply of sample generator - deafult arguments without mask, and mask values to generate all values
	/// </summary>
	virtual void get_samples(map<string, vector<T>> &data, void *params);

	/// <summary>
	/// apply of sample generator - deafult arguments with mask, and mask values to generate values in mask, where mask[i]==false. 
	/// when mask[i]==true fix values from mask_values
	/// </summary>
	virtual void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);

	/// <summary>
	/// vector api from generating samples
	/// </summary>
	virtual void get_samples(vector<vector<T>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values);

	void *new_polymorphic(string derived_name);

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(SamplesGenerator<T>)
		ADD_SERIALIZATION_FUNCS(use_vector_api)

};

/**
* Samples generator using GibbsSampler object to sample from data dist
*/
template<typename T> class GibbsSamplesGenerator : public SamplesGenerator<T> {
private:
	GibbsSampler<T> * _gibbs;
	bool _do_parallel;
public:
	GibbsSamplesGenerator();

	GibbsSamplesGenerator(GibbsSampler<T> &gibbs, bool do_parallel = true);

	void learn(const map<string, vector<T>> &data);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);

	void get_samples(vector<vector<T>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values);

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(GibbsSampler<T>)
		ADD_SERIALIZATION_FUNCS(_gibbs, _do_parallel)
};

/**
* MaskedGAN parameters
*/
class MaskedGANParams : public SerializableObject {
public:
	bool keep_original_values = false;
	ADD_CLASS_NAME(MaskedGANParams)
		ADD_SERIALIZATION_FUNCS(keep_original_values)
};

/**
* Masked GAN object
*/
template<typename T> class MaskedGAN : public SamplesGenerator<T> {
private:
	ApplyKeras generator;
	vector<vector<T>> allowed_values;
	mt19937 _gen;
	MaskedGANParams mg_params;

	T round_to_allowed_values(T in_value, vector<T>& curr_allowed_values);
	void set_params(void *params);

public:
	MaskedGAN();

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);
	void get_samples(vector<vector<T>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values);
	void get_samples_from_Z(vector<vector<T>> &data, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &Z);

	void read_from_text_file(const string& file_name);

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(MaskedGAN<T>)
		ADD_SERIALIZATION_FUNCS(generator, allowed_values, mg_params)
};

/**
* simple - just puts missing value by mask
*/
template<typename T> class MissingsSamplesGenerator : public SamplesGenerator<T> {
public:
	T missing_value = 0;
	vector<string> names;

	MissingsSamplesGenerator();

	MissingsSamplesGenerator(float miss_valu);

	void learn(const map<string, vector<T>> &data);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);

	void get_samples(vector<vector<T>> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const vector<vector<T>> &mask_values);

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(MissingsSamplesGenerator<T>)
		ADD_SERIALIZATION_FUNCS(missing_value, names)
};

MEDSERIALIZE_SUPPORT(MaskedGANParams)
MEDSERIALIZE_SUPPORT(SamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(SamplesGenerator<double>)
MEDSERIALIZE_SUPPORT(MaskedGAN<float>)
MEDSERIALIZE_SUPPORT(GibbsSamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(GibbsSamplesGenerator<double>)
MEDSERIALIZE_SUPPORT(MissingsSamplesGenerator<float>)

#endif