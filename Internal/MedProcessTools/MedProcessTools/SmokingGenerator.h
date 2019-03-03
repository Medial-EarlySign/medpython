#pragma once
#include "FeatureGenerator.h"

typedef enum {
	SMX_CURRENT_SMOKER,
	SMX_EX_SMOKER,
	SMX_YEARS_SINCE_QUITTING,
	SMX_SMOKING_YEARS,
	SMX_SMOK_PACK_YEARS,
	SMX_PLM_SMOKING_LEVEL,
	SMX_NEVER_SMOKER,
	SMX_UNKNOWN_SMOKER,
	SMX_SMOKING_QUANTITY,
	SMX_LAST
} SmokingGeneratorFields ;

/** @file
* Generation of Smoking use
*/
class SmokingGenerator : public FeatureGenerator {
public:

	// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_feature_names;
	string smoking_method = "SMOKING_ENRICHED";
	string future_ind = "1";

	// Constructor/Destructor
	SmokingGenerator() : FeatureGenerator() { generator_type = FTR_GEN_SMOKING; }
	
	~SmokingGenerator() {};

	/// The parsed fields from init command.
	/// @snippet SmokingGenerator.cpp SmokingGenerator::init
	virtual int init(map<string, string>& mapper);

	// Name
	void set_names();

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<SmokingGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// get pointers to data
	void get_p_data(MedFeatures& features, vector<float *> &_p_data);

	// Serialization
	ADD_CLASS_NAME(SmokingGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, raw_feature_names, names, tags, iGenerateWeights, smoking_method, future_ind, req_signals)
};

MEDSERIALIZE_SUPPORT(SmokingGenerator)
