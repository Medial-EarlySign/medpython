#pragma once
#include "FeatureGenerator.h"

//.......................................................................................
//.......................................................................................
// User defined calculations on other features.
// NOTE: it is the user responsibility to put these generators after their source features generators
//.......................................................................................
//.......................................................................................

class DoCalcFeatGenerator : public FeatureGenerator {
public:
	// target_feature_name as specified by the user, will be decorated for uniqueness
	string raw_target_feature_name = "";

	// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_source_feature_names;

	// source_feature_names after resolving
	vector<string> source_feature_names;

	// user function selector (e.g. sum, ratio)
	string calc_type;

	// when a source_feature == missing_value, the calculation would also be missing_value
	float missing_value;

	// for sum
	vector<float> weights;

	DoCalcFeatGenerator() : FeatureGenerator() { init_defaults();  }
	~DoCalcFeatGenerator() {};

	void init_defaults();

	// init_from_string
	int init(map<string, string>& mapper);

	// resolve raw_source_feature_names
	virtual void init(MedFeatures &features);

	// decorate raw_target_feature_name and make it unique 
	virtual void set_names();

	// Learn nothing, doCalc just uses known formulas
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	float sum(vector<float*> p_sources, int offset);

	// Serialization
	size_t get_size() { return MedSerialize::get_size(generator_type, names, raw_target_feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, names, raw_target_feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, names, raw_target_feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights); }
};