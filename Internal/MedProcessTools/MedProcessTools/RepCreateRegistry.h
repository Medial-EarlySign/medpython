#ifndef __REP_CREATE_REGISTRY_H__
#define __REP_CREATE_REGISTRY_H__

#include <string>
#include <vector>
#include <map>
#include "RepProcess.h"
#include <MedUtils/MedUtils/MedRegistry.h>

using namespace std;

typedef enum {
	REP_REGISTRY_DM,
	REP_REGISTRY_HT,
	REP_REGISTRY_PROTEINURIA,
	REP_REGISTRY_CKD,
	REP_REGISTRY_CUSTOM,
	REP_REGISTRY_LAST
} RegistryTypes;

// helper class 
typedef enum {
	REG_EVENT_DM_GLUCOSE,
	REG_EVENT_DM_HBA1C,
	REG_EVENT_DM_DRUG,
	REG_EVENT_DM_DIAGNOSES
} DMEventsTypes;

class RegistryEvent {
public:
	int time = -1;
	int event_type;
	float event_val;
	int event_severity;

	RegistryEvent() {};
	RegistryEvent(int _time, int _type, float _val, int _severity) { time = _time; event_type = _type; event_val = _val; event_severity = _severity; }

};

class RegistryDecisionRanges {
public:
	string sig_name;
	int sig_id;
	int usv_idx;
	int is_numeric;
	vector<vector<string>> categories;
	vector<vector<int>> categories_i;
	vector<pair<float, float>> ranges;
};

class RepCreateRegistry : public RepProcessor {
public:

	RegistryTypes registry; ///< type of registry to create
	vector<string> names; ///< name(s) of registry signal(s) to create

	vector<string> signals; ///< Vector of required signals, to override default ones.
	vector<string> registry_values; ///< values of registry (to appear in relevant section of dictionary)
	int time_unit = -1; ///< time-unit of registry

	bool verbose = false; ///< verbosity
	bool ht_systolic_first = false; //TODO: change to true after transferring exising models

						  // Registry specific parameters
						  // Hypertension
	vector<string> ht_identifiers; ///< identifiers (ReadCodes) of HT
	vector<string> chf_identifiers; ///< identifiers (ReadCodes) of CHF
	vector<string> mi_identifiers; ///< identifiers (ReadCodes) of MI
	vector<string> af_identifiers; ///< identifiers (ReadCodes) of AF
	bool ht_identifiers_given = false, chf_identifiers_given = false, mi_identifiers_given = false, af_identifiers_given = false;

	vector<string> ht_drugs = { "ATC_C08C____","ATC_C07B____","ATC_C07C____","ATC_C07D____","ATC_C07F____","ATC_C07A_G__","ATC_C09B____","ATC_C09D____", "ATC_C02D_A01" }; ///< drugs indicative of HT
	vector<string> ht_chf_drugs = { "ATC_C03_____" }; ///< drugs indicative of HT, unless CHF
	vector<string> ht_dm_drugs = { "ATC_C09A____", "ATC_C09C____" }; ///< drugs indicative of HT, unless DM
	vector<string> ht_extra_drugs = { "ATC_C07A_A__", "ATC_C07A_B__" }; ///< drugs indicative of HT, unless CHF/MI/AF.

	int ht_drugs_gap = 120; ///< Gap (in days) from drug input to following indication

	RepCreateRegistry() { processor_type = REP_PROCESS_CREATE_REGISTRY; }
	~RepCreateRegistry() {
		if (custom_registry != NULL) {
			delete custom_registry;
			custom_registry = NULL;
		}
	};

	// dm registry related parameters
	string dm_drug_sig = "Drug";
	vector<string> dm_drug_sets = { "ATC_A10_____" };
	string dm_diagnoses_sig = "RC";
	vector<string> dm_diagnoses_sets;
	string dm_glucose_sig = "Glucose";
	string dm_hba1c_sig = "HbA1C";
	int dm_diagnoses_severity = 4; // 3: need supporting evidence as well, 4: single code is enough
	int dm_bio_mode = 0; // bio mode - takes the FIRST suggestive test for a condition 


						 // proteinuria related parameters
						 // <name>:<0/1: is_numeric (numerics are 1)>:<categs or ranges for normal>:<categs or ranges for medium>:<categs or ranges for severe>
						 // / is the separator between signals in a real input
	vector<string> urine_tests_categories = {
		"Urine_Microalbumin:1:0,30:30,300:300,1000000",
		"UrineTotalProtein:1:0,0.15:0.15,0.60:0.60,1000000",
		"UrineAlbumin:1:0,30:30,300:300,1000000",
		"Urine_dipstick_for_protein:0:Urine_dipstick_for_protein_normal:Urine_dipstick_for_protein_medium:Urine_dipstick_for_protein_severe",
		"Urinalysis_Protein:0:Urinalysis_Protein_normal:Urinalysis_Protein_medium:Urinalysis_Protein_severe",
		"Urine_Protein_Creatinine:1:0,15:15,100:100,1000000",
		"UrineAlbumin_over_Creatinine:1:0,3.5:3.5,27:27,1000000" };


	// ckd related
	string ckd_egfr_sig = "eGFR_CKD_EPI";
	string ckd_proteinuria_sig = "Proteinuria_State";

	//custom args:
	string registry_custom_type = "";
	string registry_custom_args = "";

	/// @snippet RepCreateRegistry.cpp RepCreateRegistry::init
	int init(map<string, string>& mapper);
	void init_lists();

	void post_deserialization();


	// making sure V_ids and sigs_ids are initialized
	void init_tables(MedDictionarySections& dict, MedSignals& sigs);

	// Learning
	/// <summary> In this class there's never learning - we return 0 immediately </summary>
	int _learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors) { return 0; };

	// Applying
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points : Should be implemented for all inheriting classes </summary>
	int _apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat);

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

	void make_summary();

	// serialization
	ADD_CLASS_NAME(RepCreateRegistry)
		ADD_SERIALIZATION_FUNCS(processor_type, registry, names, signals, time_unit, req_signals, aff_signals, virtual_signals, virtual_signals_generic, registry_values,
			ht_identifiers, chf_identifiers, mi_identifiers, af_identifiers, ht_identifiers_given, chf_identifiers_given, mi_identifiers_given, af_identifiers_given,
			ht_drugs, ht_chf_drugs, ht_dm_drugs, ht_extra_drugs, ht_drugs_gap,
			dm_drug_sig, dm_drug_sets, dm_diagnoses_sig, dm_diagnoses_sets, dm_glucose_sig, dm_hba1c_sig, dm_diagnoses_severity,
			urine_tests_categories,
			ckd_egfr_sig, ckd_proteinuria_sig, registry_custom_type, registry_custom_args, custom_registry, ht_systolic_first)

private:
	MedRegistry *custom_registry = NULL;
	int bDateCode = -1;
	vector<int> create_custom_input_sigs;
	string registry_name;
	int found_op_ht = 0;
	int num_ht_measurements = 0;

	/// registry name to type
	const map<string, RegistryTypes> name2type = { { "dm" , REP_REGISTRY_DM },{ "ht", REP_REGISTRY_HT },{ "proteinuria", REP_REGISTRY_PROTEINURIA },{ "ckd", REP_REGISTRY_CKD },{ "custom", REP_REGISTRY_CUSTOM } };

	// output signal name + type
	const map<RegistryTypes, vector<pair<string, string>>> type2Virtuals = { { REP_REGISTRY_DM,{ { "DM_Registry", "T(l,l),V(f)" } } },
	{ REP_REGISTRY_HT,{ { "HT_Registry", "T(l,l),V(f)" } } },
	{ REP_REGISTRY_PROTEINURIA,{ { "Proteinuria_State", "T(i),V(f)" } } } ,
	{ REP_REGISTRY_CKD,{ { "CKD_Registry", "T(i),V(f)" } } },
	{ REP_REGISTRY_CUSTOM,{ { "CUSTOM_Registry", "T(l,l),V(f)" } } } };

	// required signals
	const map<RegistryTypes, vector<string>> type2reqSigs = { { REP_REGISTRY_DM,{ "Glucose","HbA1C","Drug","RC" } },
	{ REP_REGISTRY_HT,{ "BP","RC","Drug","BDATE","DM_Registry" } },
	{ REP_REGISTRY_PROTEINURIA ,{ "Urine_Microalbumin", "UrineTotalProtein" , "UrineAlbumin" , "Urine_dipstick_for_protein" , "Urinalysis_Protein" , "Urine_Protein_Creatinine" , "UrineAlbumin_over_Creatinine" } },
	{ REP_REGISTRY_CKD,{ "Proteinuria_State", "eGFR_CKD_EPI" } } };

	set<int> sig_ids_s;
	vector<int> sig_ids; ///< ids of signals used as input by the calculator (for faster usage at run time: save name conversions)
	vector<int> virtual_ids; ///< ids of signals created by the calculator (for faster usage at run time: save name conversions)

							 // dm related privates
	int dm_drug_idx = -1; // idx for drug signal in usvs, sig_ids, etc...
	int dm_diagnoses_idx = -1;
	vector<char> dm_drug_lut;
	vector<char> dm_diagnoses_lut;
	int dm_glucose_idx = -1;
	int dm_hba1c_idx = -1;
	//DEF     0       DM_Registry_Non_diabetic
	//DEF     1       DM_Registry_Pre_diabetic
	//DEF     2       DM_Registry_Diabetic
	vector<string> dm_reg_values = { "DM_Registry_Non_diabetic", "DM_Registry_Pre_diabetic", "DM_Registry_Diabetic" };


	// proteinuria related privates
	vector<string> proteinuria_reg_values = { "Proteinuria_Normal", "Proteinuria_Medium" , "Proteinuria_Severe" };
	vector<RegistryDecisionRanges> proteinuria_ranges;

	// CKD related states
	vector<string> ckd_reg_values = { "CKD_State_Normal" , "CKD_State_Level_1", "CKD_State_Level_2", "CKD_State_Level_3", "CKD_State_Level_4" };
	int ckd_egfr_idx = -1;
	int ckd_proteinuria_idx = -1;


	vector<int> signal_time_units; ///< time-units of all signals

								   // Registry specific functions and parameters
								   // HT
	void init_ht_registry_tables(MedDictionarySections& dict, MedSignals& sigs);
	void ht_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes);
	void ht_init_defaults();
	void fillLookupTableForHTDrugs(MedDictionary& dict, vector<char>& lut, vector<string>& sets, char val);
	void buildLookupTableForHTDrugs(MedDictionary& dict, vector<char>& lut);

	vector<string> ht_def_values = { "HT_Registry_Non_Hypertensive","HT_Registry_Hypertensive" };
	int bp_idx = 0, rc_idx = 1, drug_idx = 2, bdate_idx = 3, dm_registry_idx = 4;
	vector<char> htLut, chfLut, miLut, afLut;
	vector<char> htDrugLut;

	// DM
	void init_dm_registry_tables(MedDictionarySections& dict, MedSignals& sigs);
	void dm_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes);

	// Proteinuria
	void init_proteinuria_registry_tables(MedDictionarySections& dict, MedSignals& sigs);
	void proteinuria_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes);

	// ckd
	void init_ckd_registry_tables(MedDictionarySections& dict, MedSignals& sigs);
	void ckd_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes);

	void custom_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, 
		const vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes);
};

MEDSERIALIZE_SUPPORT(RepCreateRegistry)

#endif