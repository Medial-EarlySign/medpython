// Classes for generating samples from partial information

#ifndef __SAMPLES_GENERATOR_H__
#define __SAMPLES_GENERATOR_H__

#include <vector>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedStat/MedStat/GibbsSampler.h>
#include <MedEmbed/MedEmbed/ApplyKeras.h>
#include <MedMat/MedMat/MedMat.h>
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
	/// prepare to generate
	/// </summary>
	virtual void prepare(void *params) {};

	/// <summary>
	/// learn of sample generator
	/// </summary>
	virtual void learn(const map<string, vector<T>> &data) {};

	/// <summary>
	/// apply of sample generator - deafult arguments with mask, and mask values to generate values in mask, where mask[i]==false. 
	/// when mask[i]==true fix values from mask_values
	/// </summary>
	virtual void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);

	/// <summary>
	/// vector api from generating samples
	/// </summary>
	virtual void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values);

	/// <summary>
	/// apply of sample generator - deafult arguments with mask, and mask values to generate values in mask, where mask[i]==false. 
	/// when mask[i]==true fix values from mask_values
	/// </summary>
	virtual void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const;

	/// <summary>
	/// vector api from generating samples
	/// </summary>
	virtual void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const;

	void *new_polymorphic(string derived_name);

	void pre_serialization();
	void post_deserialization();

	virtual ~SamplesGenerator() {};

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
	bool no_need_to_clear_mem;
public:
	GibbsSamplesGenerator();

	GibbsSamplesGenerator(GibbsSampler<T> &gibbs, bool do_parallel = true, bool no_need_clear_mem = true);

	void prepare(void *params);

	void learn(const map<string, vector<T>> &data);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const;
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const;

	void pre_serialization();
	void post_deserialization();

	~GibbsSamplesGenerator();

	ADD_CLASS_NAME(GibbsSamplesGenerator<T>)
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
	vector<float> mean_feature_vals;
	vector<float> std_feature_vals;
	bool norm_by_by_file;

	T round_to_allowed_values(T in_value, const vector<T>& curr_allowed_values) const;
	void set_params(void *params);

public:
	MaskedGAN();

	void prepare(void *params);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values);
	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const;
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const;
	void get_samples_from_Z(MedMat<T> &data, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, const MedMat<T> &Z);

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
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values);
	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const;
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const;

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(MissingsSamplesGenerator<T>)
		ADD_SERIALIZATION_FUNCS(missing_value, names)
};

template<typename T> class RandomSamplesGenerator : public SamplesGenerator<T> {
public:
	T mean_value;
	T std_value;
	vector<string> names;

	RandomSamplesGenerator();

	RandomSamplesGenerator(T mean_val, T std_val);

	void learn(const map<string, vector<T>> &data);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values);
	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const;
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const;

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(RandomSamplesGenerator<T>)
		ADD_SERIALIZATION_FUNCS(mean_value, std_value, names)
};


MEDSERIALIZE_SUPPORT(MaskedGANParams)
MEDSERIALIZE_SUPPORT(SamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(SamplesGenerator<double>)
MEDSERIALIZE_SUPPORT(MaskedGAN<float>)
MEDSERIALIZE_SUPPORT(GibbsSamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(GibbsSamplesGenerator<double>)
MEDSERIALIZE_SUPPORT(MissingsSamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(RandomSamplesGenerator<float>)

#endif