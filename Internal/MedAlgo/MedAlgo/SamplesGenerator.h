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

/// @enum
/// Samples Generator Types options
enum GeneratorType
{
	GIBBS = 0, ///< "GIBBS" - to use GibbsSampler
	GAN = 1, ///< "GAN" to use GAN generator, accepts GAN path
	MISSING = 2, ///< "MISSING" to use no generator, just puts missing values where mask[i]==0
	RANDOM_DIST = 3, ///< "RANDOM_DIST" to use random normal distributaion on missing values
	UNIVARIATE_DIST = 4 ///< "UNIVARIATE_DIST" to use sampling from each feature independently
};

/// convert function for generator type to string
string GeneratorType_toStr(GeneratorType type);
/// convert function for generator
GeneratorType GeneratorType_fromStr(const string &type);

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
	void learn(const map<string, vector<T>> &data);

	/// <summary>
	/// learn of sample generator
	/// </summary>
	virtual void learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing) {};

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

	void learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing);

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
	int init(map<string, string> &mapper);

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
	
	vector<float> mean_feature_vals;
	vector<float> std_feature_vals;
	bool norm_by_by_file;

	T round_to_allowed_values(T in_value, const vector<T>& curr_allowed_values) const;
	void set_params(void *params);

public:
	MaskedGANParams mg_params;

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
		ADD_SERIALIZATION_FUNCS(generator, allowed_values, mg_params, mean_feature_vals, std_feature_vals, norm_by_by_file)
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

	void learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values);
	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const;
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const;

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(MissingsSamplesGenerator<T>)
		ADD_SERIALIZATION_FUNCS(missing_value, names)
};

/**
* puts random values from normal distribution in missing values
*/
template<typename T> class RandomSamplesGenerator : public SamplesGenerator<T> {
public:
	T mean_value;
	T std_value;
	vector<string> names;

	RandomSamplesGenerator();

	RandomSamplesGenerator(T mean_val, T std_val);

	void learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values);
	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const;
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const;

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(RandomSamplesGenerator<T>)
		ADD_SERIALIZATION_FUNCS(mean_value, std_value, names)
};

/**
* puts values in each feature selected randomly from it's distribution
*/
template<typename T> class UnivariateSamplesGenerator : public SamplesGenerator<T> {
public:
	T missing_value = MED_MAT_MISSING_VALUE;

	int min_samples; ///< minimal count of samples in strata size to use strata
	featureSetStrata strata_settings; ///< strata settings

	int init(map<string, string>& mapper);

	UnivariateSamplesGenerator();

	void learn(const map<string, vector<T>> &data, const vector<string> &learn_features, bool skip_missing);

	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values);
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values);
	void get_samples(map<string, vector<T>> &data, void *params, const vector<bool> &mask, const vector<T> &mask_values, mt19937 &rnd_gen) const;
	void get_samples(MedMat<T> &data, int sample_per_row, void *params, const vector<vector<bool>> &mask, const MedMat<T> &mask_values, mt19937 &rnd_gen) const;

	void pre_serialization();
	void post_deserialization();

	ADD_CLASS_NAME(UnivariateSamplesGenerator<T>)
		ADD_SERIALIZATION_FUNCS(feature_values, feature_val_probs, strata_feature_val_agg_prob, names, missing_value, strata_settings,
			strata_sizes, min_samples, strata_feature_val_agg_val)
private:
	//new data
	//global
	vector<vector<T>> feature_values; ///< first index is feature name, second is order index
	vector<vector<double>> feature_val_probs; ///< feature name to map of value and prob - first index is feature name, second is order index
	//by strata
	vector<vector<vector<T>>> strata_feature_val_agg_val; ///< indexed by strata, feature_name, index of sorted value
	vector<vector<vector<double>>> strata_feature_val_agg_prob; ///< indexed by strata, feature_name, index of sorted value
	

	//unordered_map<string, map<T, double>> feature_val_agg; ///< feature name to map of value and prob
	//vector<unordered_map<string, map<T, double>>> strata_feature_val_agg; ///<indexed by strata
	vector<int> strata_sizes; ///<the strata size
	vector<string> names; ///< names for all features

	T find_pos(const vector<T> &v, const vector<double> &cumsum, double p) const;
};

MEDSERIALIZE_SUPPORT(MaskedGANParams)
MEDSERIALIZE_SUPPORT(SamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(SamplesGenerator<double>)
MEDSERIALIZE_SUPPORT(MaskedGAN<float>)
MEDSERIALIZE_SUPPORT(GibbsSamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(GibbsSamplesGenerator<double>)
MEDSERIALIZE_SUPPORT(MissingsSamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(RandomSamplesGenerator<float>)
MEDSERIALIZE_SUPPORT(UnivariateSamplesGenerator<float>)

#endif