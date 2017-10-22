#pragma once
#include "FeatureGenerator.h"

class AlcoholGenerator : public FeatureGenerator {
public:

	// Constructor/Destructor
	AlcoholGenerator() : FeatureGenerator() { generator_type = FTR_GEN_ALCOHOL; req_signals.assign(1, "ALCOHOL"); }
	~AlcoholGenerator() {};

	virtual int init(map<string, string>& mapper);

	// Name
	void set_names();

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<AlcoholGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Signal Ids
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.assign(1, dict.id("ALCOHOL")); }

	// Serialization
	size_t get_size() { return MedSerialize::get_size(generator_type, names, tags); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, names, tags); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, names, tags); }
};

MEDSERIALIZE_SUPPORT(AlcoholGenerator);
