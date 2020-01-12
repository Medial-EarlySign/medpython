#pragma once
#include "FeatureGenerator.h"

/** @file
* Generation of Elixhauser score
*/
class ElixhauserGenerator : public FeatureGenerator {
private:
	vector<string> types = { "CONGESTIVE_HEART_FAILURE","CARDIAC_ARRHYTHMIAS","VALVULAR_DISEASE","PULMONARY_CIRCULATION","PERIPHERAL_VASCULAR","HYPERTENSION","PARALYSIS","OTHER_NEUROLOGICAL",
		"CHRONIC_PULMONARY","DIABETES_UNCOMPLICATED","DIABETES_COMPLICATED","HYPOTHYROIDISM","RENAL_FAILURE","LIVER_DISEASE","PEPTIC_ULCER","AIDS","LYMPHOMA","METASTATIC_CANCER","SOLID_TUMOR",
		"RHEUMATOID_ARTHRITIS","COAGULOPATHY","OBESITY","WEIGHT_LOSS","FLUID_ELECTROLYTE","BLOOD_LOSS_ANEMIA","DEFICIENCY_ANEMIAS","ALCOHOL_ABUSE","DRUG_ABUSE","PSYCHOSES","DEPRESSION" };
	vector<int> weights = {};
	vector<vector<string>> drgSets, diagSets;
	vector<vector<unsigned char>> drgLuts,diagLuts;

	void parseSets(string& init_string, vector<vector<string>>& sets);
	void initSets();

public:
	// Feature Descrption
	string drgSignalName = "DRG_IP", diagSignalName = "DIAGNOSIS_IP";
	int drgSignalId, diagSignalId;

	// parameters (should be serialized)
	int win_from = 0;///< time window for feature: win_from is the minimal time before from the prediction time
	int win_to = 360000;///< time window for feature: win_to is the maximal time before the prediction time			 
	int time_unit_win = MedTime::Undefined;			///< the time unit in which the windows are given. Default: Undefined
	int drg_time_unit_sig = MedTime::Undefined;		///< the time init in which the signal is given. (set correctly from Repository in learn and _generate)
	int diag_time_unit_sig = MedTime::Undefined;		///< the time init in which the signal is given. (set correctly from Repository in learn and _generate)

	// Constructor/Destructor
	ElixhauserGenerator() : FeatureGenerator() {
		generator_type = FTR_GEN_ELIXHAUSER; req_signals = { drgSignalName,diagSignalName }; initSets();
	}

	~ElixhauserGenerator() {};

	/// The parsed fields from init command.
	/// @snippet SmokingGenerator.cpp SmokingGenerator::init
	virtual int init(map<string, string>& mapper);

	// Name
	void set_names() { if (names.empty()) names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".Elixhauser"); tags.push_back("Elixhauser"); }

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<ElixhauserGenerator *>(generator)); }

	/// Init required tables
	void init_tables(MedDictionarySections& dict);

	// Signal Ids
	void set_signal_ids(MedSignals& sigs) { drgSignalId = sigs.sid(drgSignalName); diagSignalId = sigs.sid(diagSignalName); }

	// Learn a generator
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors);

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// Serialization
	ADD_CLASS_NAME(ElixhauserGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, types, drgSignalName, diagSignalName, drgSets, diagSets, req_signals);
};