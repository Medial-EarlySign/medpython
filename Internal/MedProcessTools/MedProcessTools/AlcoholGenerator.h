#pragma once
#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>

void generateAlcoholRangeSignal(SDateVal2* rawSignal, SDateRangeVal *outRangeSignal);

typedef enum {
	ALC_CURRENT_DRINKER,
	ALC_EX_DRINKER,
	ALC_DRINKER_YEARS_SINCE_QUITTING,
	ALC_DRINKING_YEARS,
	ALC_DRINKING_UNIT_YEARS,
	ALC_PLM_DRINKING_LEVEL,
	ALC_NEVER_DRINKER,
	ALC_UNKNOWN_DRINKER,
	ALC_DRINKER_QUANTITY,
	ALC_CURRENT_ALCOHOLIC,
	ALC_EX_ALCOHOLIC,
	ALC_LAST
} AlcoholGeneratorFields;

/** @file
* Generation of Alcohol use
*/
class AlcoholGenerator : public FeatureGenerator {
public:
	// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_feature_names;
	string future_ind = "1";

	// Constructor/Destructor
	AlcoholGenerator() : FeatureGenerator() { generator_type = FTR_GEN_ALCOHOL; req_signals.push_back("Alcohol_quantity"); req_signals.push_back("BDATE");	}
	~AlcoholGenerator() {};


	/// The parsed fields from init command.
	/// @snippet AlcoholGenerator.cpp AlcoholGenerator::init
	virtual int init(map<string, string>& mapper);
	

	// Name
	void set_names();

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<AlcoholGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// get pointers to data
	void get_p_data(MedFeatures& features, vector<float *> &_p_data);

	// Signal Ids
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.push_back(dict.id("Alcohol_quantity")); req_signal_ids.push_back(dict.id("BDATE"));
	}

	// Serialization
	ADD_CLASS_NAME(AlcoholGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, names, tags, future_ind, req_signals)
};

MEDSERIALIZE_SUPPORT(AlcoholGenerator)
