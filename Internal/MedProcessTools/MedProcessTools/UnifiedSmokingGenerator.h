#pragma once
#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>

typedef enum {
	SMX_UNIFIED_CURRENT_SMOKER,
	SMX_UNIFIED_EX_SMOKER,
	SMX_UNIFIED_UNKNOWN_SMOKER,
	SMX_UNIFIED_NEVER_SMOKER,
	SMX_UNIFIED_PASSIVE_SMOKER,
	SMX_UNIFIED_DAYS_SINCE_QUITTING,
	SMX_UNIFIED_SMOK_PACK_YEARS_MAX,
	SMX_UNIFIED_SMOK_PACK_YEARS_LAST,
	SMX_UNIFIED_SMOK_PACK_YEARS,
	UNIFIED_NLST_CRITERION,
	SMX_UNIFIED_YEARS_SINCE_QUITTING,
	SMX_UNIFIED_SMOKING_INTENSITY,
	SMX_UNIFIED_SMOKING_YEARS,
	SMX_UNIFIED_LAST
} UnifiedSmokingGeneratorFields;

typedef enum {
	UNKNOWN_SMOKER,
	NEVER_SMOKER,
	PASSIVE_SMOKER,
	EX_SMOKER,
	CURRENT_SMOKER,
	NEVER_OR_EX_SMOKER
} SMOKING_STATUS;

struct RangeStatus
{
	int startDate, endDate;
	SMOKING_STATUS smokingStatus;
};

#define NA_SMOKING_DATE -1
#define MAX_PACK_YEARS 200
#define MAX_INTENSITY_TO_TRIM 140
#define MAX_INTENSITY_TO_REMOVE 200
#define AGE_AT_START_SMOKING 20
#define PACK_SIZE 20

class UnifiedSmokingGenerator : public FeatureGenerator {
public:
	float nlstPackYears, nlstQuitTimeYears, nlstMinAge, nlstMaxAge;
	bool nonDefaultNlstCriterion, useDataComplition;
	FILE *fp = stderr;
	string debug_file = "";
	set<vector<SMOKING_STATUS>> possibleCombinations;
	float timeSinceQuittingModelSlope = missing_val;
	float timeSinceQuittingModelConst = missing_val;

	// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_feature_names;
	// Constructor/Destructor
	UnifiedSmokingGenerator() : FeatureGenerator() { 
		missing_val = MED_MAT_MISSING_VALUE, generator_type = FTR_GEN_UNIFIED_SMOKING;
		// Set NLST default values:
		nlstMinAge = 55;
		nlstMaxAge = 80;
		nlstPackYears = 30;
		nlstQuitTimeYears = 15;
		nonDefaultNlstCriterion = false;
		useDataComplition = false;
	}

	~UnifiedSmokingGenerator() {
		if (fp != stderr)
		{
			fclose(fp);
		}  
		char *smokingStatusDesc[] = { "UNKNOWN_SMOKER", "NEVER_SMOKER", "PASSIVE_SMOKER", "EX_SMOKER", "CURRENT_SMOKER" };
		for (auto &comb : possibleCombinations)
		{
			for (auto stat : comb)
				cout << smokingStatusDesc[stat] << " ";
			cout << endl;
		}

	};

	void init_tables(MedDictionarySections& dict);

	/// The parsed fields from init command.
	/// @snippet UnifiedSmokingGenerator.cpp UnifiedSmokingGenerator::init
	virtual int init(map<string, string>& mapper);

	virtual int update(map<string, string>& mapper);

	// Name
	void set_names();

	void set_signal_ids(MedSignals& sigs);

	void fit_for_repository(MedPidRepository &rep);

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<UnifiedSmokingGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository & rep, const MedSamples & samples, vector<RepProcessor*> processors);

	void calcSmokingDuration(int neverSmoker, int unknownSmoker, vector<RangeStatus>& smokeRanges, int birthDate, int lastPackYearsDate, UniversalSigVec & SmokingDurationUsv, int testDate, float & smokingDurationBeforeLastPackYears, float & smokingDuration);

	void getLastSmokingDuration(int birthDate, UniversalSigVec & SmokingDurationUsv, int testDate, int &lastDurationDate, float &lastDurationValue);

	void calcPackYears(UniversalSigVec & SmokingPackYearsUsv, int testDate, int & neverSmoker, int & currentSmoker, int & formerSmoker, int & lastPackYearsDate, float & lastPackYears, float & maxPackYears);

	void calcQuitTimeOriginalData(PidDynamicRec& rec, UniversalSigVec & smokingStatusUsv, UniversalSigVec & quitTimeUsv, int testDate, int formerSmoker, int neverSmoker, int currentSmoker, float & daysSinceQuittingOriginal);
	
	void calcPackYearsOriginalData(int testDate, int lastPackYearsDate, float lastPackYears, float & lastPackYearsOriginal, UniversalSigVec SmokingIntensityUsv, UniversalSigVec SmokingDurationUsv);
	
	void fixPackYearsSmokingIntensity(float smokingDurationSinceLastPackYears, float & smokingIntensity, float smokingDuration, float & lastPackYears, float & maxPackYears);

	void printDebug(vector<RangeStatus>& smokeRanges, int qa_print, UniversalSigVec & smokingStatusUsv, UniversalSigVec & SmokingIntensityUsv, int birthDate, int testDate, vector<pair<SMOKING_STATUS, int>>& smokingStatusVec, PidDynamicRec & rec, UniversalSigVec & quitTimeUsv, UniversalSigVec & SmokingPackYearsUsv, float smokingIntensity, float smokingDuration, float yearsSinceQuitting, float maxPackYears);

	void addDataToMat(vector<float*>& _p_data, int index, int i, int age, int currentSmoker, int formerSmoker, float daysSinceQuitting, float daysSinceQuittingOriginal, float maxPackYears, float lastPackYears, float lastPackYearsOriginal, int neverSmoker, int unknownSmoker, int passiveSmoker, float yearsSinceQuitting, float smokingIntensity, float smokingDuration);

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// get pointers to data
	void get_p_data(MedFeatures& features, vector<float *> &_p_data);

	int calcNlst(int age, int unknownSmoker, int daysSinceQuitting, float lastPackYears);

	SMOKING_STATUS val2SmokingStatus(int sigVal, int smokingStatusSid, PidRec &rec);

	void genSmokingVec(PidDynamicRec & rec, UniversalSigVec & smokingStatusUsv, vector<pair<SMOKING_STATUS, int>>& smokingStatusVec, int testDate, int & unknownSmoker, int & neverSmoker, int & passiveSmoker, int & formerSmoker, int & currentSmoker);

	void genFirstLastSmokingDates(PidDynamicRec & rec, UniversalSigVec & smokingStatusUsv, UniversalSigVec & quitTimeUsv, int testDate, map<SMOKING_STATUS, pair<int, int>>& smokingStatusDates, vector<int>& dates, int birth_date);

	void genSmokingStatus(map<SMOKING_STATUS, pair<int, int>>& smokingStatusDates, vector<int>& dates, int testDate, int birthDate, vector<pair<SMOKING_STATUS, int>>& smokingStatusVec);

	void genSmokingRanges(vector<pair<SMOKING_STATUS, int>>& smokingStatusVec, int testDate, int birthDate, vector<RangeStatus>& smokeRanges);

	void genLastStatus(vector<RangeStatus>& smokeRanges, int & unknownSmoker, int & neverSmoker, int & formerSmoker, int & currentSmoker, int & passiveSmoker);

	void calcQuitTime(vector<RangeStatus> &smokeRange, int formerSmoker, int neverSmoker, int currentSmoker, int testDate, int birthDate, float & daysSinceQuitting, float & yearsSinceQuitting);
	
	void calcSmokingIntensity(UniversalSigVec & SmokingIntensityUsv, int testDate, int neverSmoker, float & smokingIntensity);

	void getQuitAge(PidDynamicRec& rec, int lastDate, float &ageAtEx, float &deltaTime);
		
	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;
	// Serialization
	ADD_CLASS_NAME(UnifiedSmokingGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, raw_feature_names, names, tags, iGenerateWeights, req_signals, timeSinceQuittingModelSlope, timeSinceQuittingModelConst, useDataComplition)
private:
	vector<vector<char>> smoke_status_luts;
	int smoke_status_sec_id;

	int smoking_quit_date_id, smoking_status_id, smoking_intensity_id, smoking_duration_id, smoking_pack_years_id;
	int bdate_sid;
};

MEDSERIALIZE_SUPPORT(UnifiedSmokingGenerator)
