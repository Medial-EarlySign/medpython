#pragma once
#include "FeatureGenerator.h"

class DoCalcFeatGenerator : public FeatureGenerator {
public:
	// Signal Id
	int signalId;
	string signalName;

	void init_defaults();
	int init(map<string, string>& mapper);

	// Constructor/Destructor
	DoCalcFeatGenerator() : FeatureGenerator() { init_defaults();  }
	~DoCalcFeatGenerator() {};

	virtual void set_names();
	// Learn nothing, doCalc just uses known formulas
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Serialization
	size_t get_size() { return MedSerialize::get_size(generator_type, names); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, names); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, names); }
	DEC_FEATURE_GENERATOR(DoCalcFeatGenerator);

};