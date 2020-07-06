#pragma once

#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <json/json.hpp>
/**
* calculate drug coverage of prescription time in defined the time window. a value between 0 to 1.
*/
class DiabetesFinderGenerator : public FeatureGenerator {
	// dm related privates
	int df_drug_sid = -1; // idx for drug signal in usvs, sig_ids, etc...
	int df_diagnosis_sid = -1;
	int df_coded_sid = -1;
	int df_glucose_sid = -1;
	int df_hba1c_sid = -1;
	int df_drug_section_id = -1;
	int df_diagnosis_section_id = -1;
	int df_coded_section_id = -1;

	enum {
		REASON_RECENT_LABS = 1,
		REASON_RECENT_DRUGS = 2,
		REASON_RECENT_DIAGNOSTIC = 4,
		REASON_PAST_LABS = 8,
		REASON_PAST_DRUGS = 16,
		REASON_PAST_DIAGNOSTIC = 32,
	};

	enum {
		DFG_DIABETES_EVENT_GLUCOSE,
		DFG_DIABETES_EVENT_HBA1C,
		DFG_DIABETES_EVENT_DRUG,
		DFG_DIABETES_EVENT_DIAGNOSIS,
		DFG_DIABETES_EVENT_PG_DURING_OGTT,  //Plasma Glucose during Oral Glucose Tolerance Test
	};

	class DiabetesEvent {
	public:
		int time = -1;
		int de_type;
		float val;
		bool is_second = false;
		bool is_first = false;
		bool is_non_dm = false;
		string reason = "";

		DiabetesEvent() {};
		DiabetesEvent(int _type, int _time,  float _val) { time = _time; de_type = _type; val = _val; }
	};
	vector<unsigned char> df_drug_lut;
	vector<unsigned char> df_diagnosis_lut;
	vector<unsigned char> df_coded_lut;

	int _resolve(PidDynamicRec& rec, vector<DiabetesEvent>& df_events, int coded_date, int coded_val, int calc_time, json& json_out);
public:

	bool df_score_is_flag = true;
	bool df_score_is_bitmask = false;

	// dm registry related parameters
	vector<string> df_drug_sets = { "ATC_A10_____" };
	//TODO - Diabetes diagnosis sets?
	vector<string> df_coded_sets; // if not given explicitly will be defaulted to df_diagnosis sets

	vector<string> df_diagnosis_sets;
	string df_diagnosis_sig = "RC"; // These are optional sig + diagnosis codes (in df_diagnosis_sets) that can point to Diabetes but are not in the coded sig and sets.
	string df_coded_sig = "RC"; // This is the signal that will define who's coded (along with the df_coded_sets)
	string df_glucose_sig = "Glucose";
	string df_hba1c_sig = "HbA1C";
	string df_drug_sig = "Drug";
	int df_diagnoses_severity = 4; // 3: need supporting evidence as well, 4: single code is enough
	int df_bio_mode = 0; // bio mode - takes the FIRST suggestive test for a condition 
	int df_output_verbosity = 2; // 1 - "score" only, 2 - add reason for detected date, 3 - add all supporting evidences in history
	int df_output_non_dm_period = 0; // 1 - report also the period of non dm tests prior to the predition time
	
	int df_past_event_days = (365)*3;
	float df_by_single_glucose = 200.0f;
	float df_by_single_hba1c = 7.0f;
	float df_by_second_glucose = 126.0f;
	float df_by_second_hba1c = 6.5f;
	float df_pre_d_hba1c = 5.8f;
	float df_pre_d_glucose = 101.0f;
	int df_by_second_time_delta_days = (365) * 2;

	// Constructor/Destructor
	DiabetesFinderGenerator() : FeatureGenerator() { 
		generator_type = FTR_GEN_DIABETES_FINDER; 
		//names.push_back("df"); 
		req_signals.push_back(df_glucose_sig);
		req_signals.push_back(df_hba1c_sig);
		req_signals.push_back(df_drug_sig);
		req_signals.push_back(df_diagnosis_sig);
		if(df_coded_sig != df_diagnosis_sig)
			req_signals.push_back(df_coded_sig);
		init_defaults();
	};
	~DiabetesFinderGenerator() {};

	/// The parsed fields from init command.
	int init(map<string, string>& mapper);
	
	void init_tables(MedDictionarySections& dict);
	void set_signal_ids(MedSignals& sigs);

	void init_defaults();

	// Naming
	void set_names() { if (names.empty()) names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".DiabetesFinder"); tags.push_back("Diabetes"); }

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<DiabetesFinderGenerator *>(generator)); }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);
	//float get_value(PidDynamicRec &rec, int idx, int time, int sig_outcomeTime);

	// Serialization
	ADD_CLASS_NAME(DiabetesFinderGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, names, tags, iGenerateWeights, req_signals, df_drug_sets, df_coded_sets, df_diagnosis_sets, df_diagnosis_sig, df_coded_sig, df_glucose_sig, 
			df_hba1c_sig, df_drug_sig, df_past_event_days, df_by_single_glucose, df_by_second_glucose, df_by_second_hba1c, df_by_single_hba1c, df_by_second_time_delta_days, df_pre_d_hba1c, df_pre_d_glucose,
			df_diagnoses_severity, df_output_verbosity, df_output_non_dm_period)
};

MEDSERIALIZE_SUPPORT(DiabetesFinderGenerator);
