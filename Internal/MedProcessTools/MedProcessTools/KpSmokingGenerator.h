#pragma once
#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>

typedef enum {
	SMX_KP_CURRENT_SMOKER,
	SMX_KP_EX_SMOKER,
	SMX_KP_UNKNOWN_SMOKER,
	SMX_KP_NEVER_SMOKER,
	SMX_KP_PASSIVE_SMOKER,
	SMX_KP_DAYS_SINCE_QUITTING,
	SMX_KP_SMOK_PACK_YEARS_MAX,
	SMX_KP_SMOK_PACK_YEARS_LAST,
	NLST_CRITERION,
	SMX_KP_LAST
} KpSmokingGeneratorFields;

#define KP_NEVER_SMOKER_QUIT_TIME (19000101)
class KpSmokingGenerator : public FeatureGenerator {
public:
	float nlstPackYears, nlstQuitTimeYears, nlstMinAge, nlstMaxAge;
	bool nonDefaultNlstCriterion;

	// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_feature_names;
	// Constructor/Destructor
	KpSmokingGenerator() : FeatureGenerator() { missing_val = MED_MAT_MISSING_VALUE, generator_type = FTR_GEN_KP_SMOKING; }

	~KpSmokingGenerator() {};

	/// The parsed fields from init command.
	/// @snippet KpSmokingGenerator.cpp KpSmokingGenerator::init
	virtual int init(map<string, string>& mapper);

	// Name
	void set_names();

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<KpSmokingGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// get pointers to data
	void get_p_data(MedFeatures& features, vector<float *> &_p_data);

	int calcNlst(int age, int unknownSmoker, int daysSinceQuitting, float lastPackYears);

	// Serialization
	ADD_CLASS_NAME(KpSmokingGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, raw_feature_names, names, tags, iGenerateWeights, req_signals)
};

MEDSERIALIZE_SUPPORT(KpSmokingGenerator)
