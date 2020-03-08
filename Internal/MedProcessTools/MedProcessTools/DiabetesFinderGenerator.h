#pragma once

#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <json/json.hpp>
/**
* calculate drug coverage of prescription time in defined the time window. a value between 0 to 1.
*/
class DiabetesFinderGenerator : public FeatureGenerator {
	// dm related privates
	int dm_drug_idx = -1; // idx for drug signal in usvs, sig_ids, etc...
	int dm_diagnoses_idx = -1;
	int dm_glucose_idx = -1;
	int dm_hba1c_idx = -1;

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

		DiabetesEvent() {};
		DiabetesEvent(int _type, int _time,  float _val) { time = _time; de_type = _type; val = _val; }
	};
	vector<unsigned char> dm_drug_lut;
	vector<unsigned char> dm_diagnosis_lut;
	vector<unsigned char> dm_coded_lut;

	int _resolve(vector<DiabetesEvent>& dm_events, int calc_time, json& json_out);
public:

	// dm registry related parameters
	vector<string> dm_drug_sets = { "ATC_A10_____" };
	//TODO - Diabetes diagnosis sets?
	vector<string> dm_coded_sets;

	vector<string> dm_diagnosis_sets;
	string dm_diagnosis_sig = "RC";
	string dm_coded_sig = "RC";
	string dm_glucose_sig = "Glucose";
	string dm_hba1c_sig = "HbA1C";
	string dm_drug_sig = "Drug";
	int dm_diagnoses_severity = 4; // 3: need supporting evidence as well, 4: single code is enough
	int dm_bio_mode = 0; // bio mode - takes the FIRST suggestive test for a condition 
	
	int dm_past_event_days = (365)*3;
	float dm_by_single_glucose = 200.0f;
	float dm_by_second_glucose = 126.0f;
	float dm_by_second_hba1c = 6.5f;
	float dm_by_second_time_delta_days = (365) * 2;
	float dm_by_second_time_delta = -1;

	// Constructor/Destructor
	DiabetesFinderGenerator() : FeatureGenerator() { 
		generator_type = FTR_GEN_DIABETES_FINDER; 
		//names.push_back("df"); 
		req_signals.push_back(dm_glucose_sig);
		req_signals.push_back(dm_hba1c_sig);
		req_signals.push_back(dm_drug_sig);
		req_signals.push_back(dm_diagnosis_sig);
		if(dm_coded_sig != dm_diagnosis_sig)
			req_signals.push_back(dm_coded_sig);
		init_defaults();
	};
	~DiabetesFinderGenerator() {};

	/// The parsed fields from init command.
	int init(map<string, string>& mapper);
	
	void init_tables(MedDictionarySections& dict);

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
		ADD_SERIALIZATION_FUNCS(generator_type, names, tags, iGenerateWeights, req_signals, dm_drug_sets, dm_coded_sets, dm_diagnosis_sets, dm_diagnosis_sig, dm_coded_sig, dm_glucose_sig, 
			dm_hba1c_sig, dm_drug_sig, dm_past_event_days, dm_by_single_glucose, dm_by_second_glucose, dm_by_second_hba1c, dm_by_second_time_delta_days, dm_by_second_time_delta)
};

MEDSERIALIZE_SUPPORT(DiabetesFinderGenerator);
