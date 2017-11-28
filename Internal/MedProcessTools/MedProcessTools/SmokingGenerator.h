#pragma once
#include "FeatureGenerator.h"

void generateSmokingRangeSignal(SDateVal2* rawSignal, SDateRangeVal *outRangeSignal);

class SmokingGenerator : public FeatureGenerator {
public:

	// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_feature_names;
	string smoking_method = "SMOKING_ENRICHED";
	string future_ind = "0";

	// Constructor/Destructor
	//SmokingGenerator() : FeatureGenerator() { generator_type = FTR_GEN_SMOKING; req_signals.assign(1, "SMOKING_ENRICHED"); }
	SmokingGenerator() : FeatureGenerator() { generator_type = FTR_GEN_SMOKING; req_signals.push_back("Smoking_quantity");  req_signals.push_back("SMOKING_ENRICHED"); req_signals.push_back("BYEAR"); }
	
	~SmokingGenerator() {};

	virtual int init(map<string, string>& mapper);

	// Name
	void set_names();

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<SmokingGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Signal Ids
	//void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.assign(1, dict.id("SMOKING_ENRICHED")); }
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.push_back(dict.id("SMOKING_ENRICHED")); req_signal_ids.push_back(dict.id("Smoking_quantity")); req_signals.push_back("BYEAR");
	}

	// Serialization
	virtual int version() { return  1; } // ihadanny 20171128 - added the smoking_method and future_ind to the serialization 
	size_t get_size() { return MedSerialize::get_size(generator_type, raw_feature_names, names, tags, iGenerateWeights, smoking_method, future_ind); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, raw_feature_names, names, tags, iGenerateWeights, smoking_method, future_ind); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, raw_feature_names, names, tags, iGenerateWeights, smoking_method, future_ind); }
};

MEDSERIALIZE_SUPPORT(SmokingGenerator);
