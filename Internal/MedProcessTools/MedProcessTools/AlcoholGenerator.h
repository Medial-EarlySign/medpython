#pragma once
#include "FeatureGenerator.h"

void generateAlcoholRangeSignal(SDateVal2* rawSignal, SDateRangeVal *outRangeSignal);

class AlcoholGenerator : public FeatureGenerator {
public:
	// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_feature_names;
	string future_ind = "0";

	// Constructor/Destructor
	AlcoholGenerator() : FeatureGenerator() { generator_type = FTR_GEN_ALCOHOL; req_signals.push_back("Alcohol_quantity"); req_signals.push_back("BYEAR");	}
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
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.push_back(dict.id("Alcohol_quantity")); req_signal_ids.push_back(dict.id("BYEAR"));
	}

	// Serialization
	size_t get_size() { return MedSerialize::get_size(generator_type, names, tags, future_ind); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, names, tags, future_ind); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, names, tags, future_ind); }
};

MEDSERIALIZE_SUPPORT(AlcoholGenerator);
