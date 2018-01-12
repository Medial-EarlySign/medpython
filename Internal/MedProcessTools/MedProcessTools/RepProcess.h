#ifndef _REP_PROCESS_H_
#define _REP_PROCESS_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedProcessTools/MedProcessTools/MedValueCleaner.h"


#define DEFAULT_REP_CLNR_NTHREADS 8

//.......................................................................................
/** Define types of repository processors
*/
//.......................................................................................
/// Define types of repository processors
typedef enum {
	REP_PROCESS_MULTI, ///<"multi_processor" or "multi" to activate RepMultiProcessor
	REP_PROCESS_BASIC_OUTLIER_CLEANER,///<"basic_outlier_cleaner" or "basic_cln" to activate RepBasicOutlierCleaner
	REP_PROCESS_NBRS_OUTLIER_CLEANER,///<"nbrs_outlier_cleaner" or "nbrs_cln" to activate RepNbrsOutlierCleaner
	REP_PROCESS_CONFIGURED_OUTLIER_CLEANER,///<"configured_outlier_cleaner" or "conf_cln" to activate RepConfiguredOutlierCleaner
	REP_PROCESS_RULEBASED_OUTLIER_CLEANER,///<"rulebased_outlier_cleaner" or "rule_cln" to activate RepRuleBasedOutlierCleaner
	REP_PROCESS_LAST
} RepProcessorTypes;

/** @file
* RepProcessor is the parent class for processing a MedRepository or PidDynamicRec\n
* Basic functionalities:\n
*		learn : learn the processoring parameters from a given list of ids and a rpository \n
*		apply : process a dynamic PidDynamicRec\n
*/
class RepProcessor : public SerializableObject {
public:


	RepProcessorTypes processor_type; ///< type of repository processor

	unordered_set<string> req_signals; ///< names of signals required for processsing
	unordered_set<int> req_signal_ids; ///< ids of signals required for processing

	unordered_set<string> aff_signals; ///< names of signals affected by processing
	unordered_set<int> aff_signal_ids; ///< ids of signals affected by processing

	// create a new rep_processor
	/// <summary> create a new repository processor from name </summary>
	static RepProcessor *make_processor(string name);
	/// <summary> create a new repository processor from name and a parameters string</summary>
	static RepProcessor *make_processor(string type, string params);
	/// <summary> create a new repository processor from type </summary>
	static RepProcessor *make_processor(RepProcessorTypes type);
	/// <summary> create a new repository processor from type and a parameters string</summary>
	static RepProcessor *make_processor(RepProcessorTypes type, string params);

	/// <summary> create a new repository processor from parameters string which contains rp_type </summary>
	static RepProcessor *create_processor(string &params);

	/// <summary> initialize from a params object :  Should be implemented for inheriting classes that have parameters </summary>
	virtual int init(void *params) { return 0; };
	/// <summary> initialize from a map :  Should be implemented for inheriting classes that have parameters </summary>
	virtual int init(map<string, string>& mapper) { return 0; };
	/// <summary> initialize to default values :  Should be implemented for inheriting classes that have parameters </summary>
	virtual void init_defaults() {};

	/// <summary> set signal-name :  Should be implemented for inheriting classes that have signalName </summary>
	virtual void set_signal(const string& _signalName) { return; };

	/// <summary> set signal-name :  Should be implemented for inheriting classes that have signalId  </summary>
	virtual void set_signal_ids(MedDictionarySections& dict) { return; }

	// Required Signals functions : get all signals that are required by the processor
	/// <summary> Append required signal names to set : parent function just uses req_signals  </summary>
	virtual void get_required_signal_names(unordered_set<string>& signalNames);
	// Required Signals functions : get all signals that are required by the processor
	/// <summary> Append required signal names to set only if processor is actually required to produce any of preReqSignals : parent function just uses req_signals  </summary>
	virtual void get_required_signal_names(unordered_set<string>& signalNames, unordered_set<string> preReqSignals);
	
	/// <summary> Fill req_signal_ids : parent function just fills from req_signals </summary>
	virtual void set_required_signal_ids(MedDictionarySections& dict);
	
	// Required Signals functions : get all signals that are required by the processor
	/// <summary> Append required signal names to set : parent function just uses req_signals  </summary>
	virtual void get_required_signal_ids(unordered_set<int>& signalIds);
	// Required Signals functions : get all signals that are required by the processor
	/// <summary> Append required signal names to set only if processor is actually required to produce any of preReqSignals : parent function just uses req_signals  </summary>
	virtual void get_required_signal_ids(unordered_set<int>& signalIds, unordered_set<int> preReqSignals);

	// Affected Signals functions;
	/// <summary> Fill aff_signal_ids : parent function just fills from aff_signals </summary>
	virtual void set_affected_signal_ids(MedDictionarySections& dict);
	/// <summary>  Check if a signal is affected by processor </summray>
	/// <returns> true if affected, false if not </returns>
	bool is_signal_affected(int signalId) {return (aff_signal_ids.find(signalId) != aff_signal_ids.end());}
	bool is_signal_affected(string& signalName) { return (aff_signals.find(signalName) != aff_signals.end()); }

	/// <summary> Init required tables : Should be implemented for inheriting classes that have such tables </summary>
	virtual void init_tables(MedDictionarySections& dict) { return; }

	// Learning
	/// <summary> learn processing model on a subset of ids. Apply set of preceeding processors on DynamicPidRec before learning : 
	// Should be implemented for inheriting classes that require learning </summary>
	virtual int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) {  return 0; };
	/// <summary> learn processing model on a subset of ids only if required. Apply set of preceeding processors on DynamicPidRec before learning : 
	// May be implemented for inheriting classes that require learning </summary>
	virtual int _conditional_learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors, unordered_set<int>& neededSignalIds) ;

	/// <summary> learn processing model on a subset of ids. Apply set of preceeding processors on DynamicPidRec before learning : 
	// Should be implemented for inheriting classes that require learning </summary>
//	virtual int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) { return 0; };

	// Learning envelopes - Here because of issues with overloading and inheritance
	/// <summary> learn processing model on a subset of ids. Apply set of preceeding processors on DynamicPidRec before learning </summary>
	int learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) { return _learn(rep,ids,prev_processors); };
	/// <summary> learn on all ids in repository </summary>
	int learn(MedPidRepository& rep, vector<RepProcessor *>& prev_processors) { return _learn(rep, rep.pids, prev_processors); }
	/// <summary> learn on subset of ids without preceesing processors  </summary>
	int learn(MedPidRepository& rep, vector<int>& ids) { vector<RepProcessor *> temp;  return _learn(rep, ids, temp); }
	/// <summary> learn on all ids in repository without preceesing processors  </summary>
	int learn(MedPidRepository& rep) { vector<RepProcessor *> temp; return learn(rep, temp); }
	/// <summary> learn processing model on a subset of ids only if required. Apply set of preceeding processors on DynamicPidRec before learning </summary>
	int conditional_learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors, unordered_set<int>& neededSignalIds) 
		{ return _conditional_learn(rep, ids, prev_processors, neededSignalIds); }

	// Applying
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points : Should be implemented for all inheriting classes </summary>
	virtual int _apply(PidDynamicRec& rec, vector<int>& time_points) = 0;
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points only if required : May be implemented for inheriting classes </summary>
	virtual int _conditional_apply(PidDynamicRec& rec, vector<int>& time_points, unordered_set<int>& neededSignalIds);

	// Applying envelopes - Here because of issues with overloading and inheritance
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points</summary>
	int apply(PidDynamicRec& rec, vector<int>& time_points) {return _apply(rec, time_points);}
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points only if required : if any of the signals in neededSignalIds is actually affected by processor </summary>
	int conditional_apply(PidDynamicRec& rec, vector<int>& time_points, unordered_set<int>& neededSignalIds) { return _conditional_apply(rec, time_points, neededSignalIds); }
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points given by samples </summary>
	int apply(PidDynamicRec& rec, MedIdSamples& samples);
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points given by samples only if required </summary>
	int conditional_apply(PidDynamicRec& rec, MedIdSamples& samples, unordered_set<int>& neededSignalIds);

	// Serialization (including type)
	/// <summary> get size of processor + processor_type </summary>
	size_t get_processor_size();
	/// <summary> seialize processor + processor_type </summary>
	size_t processor_serialize(unsigned char *blob);

	/// <summary> optional printing of processor </summary>
	virtual void print() { fprintf(stderr, "No implementation for print()\n"); }
};

// Utilities
/// <summary> get RepProcessorTypes from name </summary>
RepProcessorTypes rep_processor_name_to_type(const string& procesor_name);

//.......................................................................................
/** RepMultiProcessor is a repository processor which contains a vector of simpler processors that can be 
* learned/applied  in parallel. Useful for applying same cleaners on a set of signals, for example
*/
//.......................................................................................
class RepMultiProcessor : public RepProcessor {
public:
	vector<RepProcessor *> processors; ///< Set of processors

	/// <summary> Constructor </summary>
	RepMultiProcessor() { processor_type = REP_PROCESS_MULTI; };

	/// <summary> Add processors to set  </summary>
	void add_processors_set(RepProcessorTypes type, vector<string>& signals);
	/// <summary> Add processors to set with initialization string  </summary>
	void add_processors_set(RepProcessorTypes type, vector<string>& signals, string init_string);
	/// <summary> Required Signals ids : Fill the member vector - req_signal_ids </summary>
	void set_required_signal_ids(MedDictionarySections& dict); 

	/// <summary> Required Signals names : Fill the unordered set signalNames </summary>
	void get_required_signal_names(unordered_set<string>& signalNames);
	/// <summary> Append required signal names to set only if processor is actually required to produce any of preReqSignals </summary>
	virtual void get_required_signal_names(unordered_set<string>& signalNames, unordered_set<string> preReqSignals);

	/// <summary> Required Signals ids : Fill the unordered set signalNames </summary>
	void get_required_signal_ids(unordered_set<int>& signalIds);
	/// <summary> Append required signal names to set only if processor is actually required to produce any of preReqSignals </summary>
	virtual void get_required_signal_ids(unordered_set<int>& signalIds, unordered_set<int> preReqSignals);

	/// <summary> Affected Signals : Fill the member set aff_signal_ids </summary>
	void set_affected_signal_ids(MedDictionarySections& dict); 

	/// <summary> Set signal-ids for all linked signals </summary>
	void set_signal_ids(MedDictionarySections& dict); 

	/// <summary> learn processors </summary>
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors);
	int _conditional_learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors, unordered_set<int>& neededSignalIds);
	
	/// <summary> Apply processors </summary>
	int _apply(PidDynamicRec& rec, vector<int>& time_points);
	/// <summary> Apply processors that affect any of the needed signals </summary>
	int _conditional_apply(PidDynamicRec& rec, vector<int>& time_points, unordered_set<int>& neededSignals);

	/// serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	/// <summary> Print processors information </summary>
	void print() { for (auto& processor : processors) processor->print(); }
};

#define DEF_REP_TRIMMING_SD_NUM 7
#define DEF_REP_REMOVING_SD_NUM 14

//.......................................................................................
/**
* A simple cleaner considering each value of a certain signal separatley
*/
//.......................................................................................
class RepBasicOutlierCleaner : public RepProcessor, public MedValueCleaner {
public:

	string signalName; 	///< name of signal to clean
	int signalId;	///< id of signal to clean
	int time_channel = 0; ///< time channel to consider in cleaning
	int val_channel = 0; ///< value cahnnel to consider in cleaning

	/// <summary> default constructor </summary>
	RepBasicOutlierCleaner() { init_defaults(); }
	/// <summary> default constructor + setting signal name </summary>
	RepBasicOutlierCleaner(const string& _signalName) { init_defaults(); signalId = -1; signalName = _signalName; init_lists(); }
	/// <summary> default constructor + setting signal name + initialize from string </summary>
	RepBasicOutlierCleaner(const string& _signalName, string init_string) { init_defaults(); signalId = -1; signalName = _signalName; init_from_string(init_string); }
	/// <summary> default constructor + setting signal name + initialize from parameters </summary>
	RepBasicOutlierCleaner(const string& _signalName, ValueCleanerParams *_params) {signalId = -1; signalName = _signalName; init_lists() ; MedValueCleaner::init(_params);}

	/// <summary> Initialize to default values </summary>
	void init_defaults() {
		processor_type = REP_PROCESS_BASIC_OUTLIER_CLEANER;
		params.trimming_sd_num = DEF_REP_TRIMMING_SD_NUM; params.removing_sd_num = DEF_REP_REMOVING_SD_NUM; params.nbrs_sd_num = 0;
		params.take_log = 0;
		params.doTrim = params.doRemove = true;
		signalId = -1;
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};

	/// <summary> Set signal name and fill affected and required signals sets </summary> 
	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName; init_lists(); }

	/// <summary> Set signal id </summary>
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }

	/// <summary> Fill required- and affected-signals sets </summary>
	int init(void *processor_params) { return MedValueCleaner::init(processor_params); };
	/// The parsed fields from init command.
	/// @snippet RepProcess.cpp RepBasicOutlierCleaner::init
	virtual int init(map<string, string>& mapper);
	/// Fill req- and aff-signals vectors
	void init_lists();

	/// <summary> learn cleaning boundaries </summary>
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	/// <summary> Learning : learn cleaning boundaries using MedValueCleaner's iterative approximation of moments </summary>
	int iterativeLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	/// <summary> Learning : learn cleaning boundaries using MedValueCleaner's quantile approximation of moments </summary>
	int quantileLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);

	/// <summary> Apply cleaning model </summary>
	int _apply(PidDynamicRec& rec, vector<int>& time_points);

	/// Serialization
	int version() { return 1; }
	ADD_SERIALIZATION_FUNCS(processor_type, signalName, time_channel, val_channel, req_signals, aff_signals, params.take_log, params.missing_value, params.doTrim, params.doRemove, 
		trimMax, trimMin, removeMax, removeMin)

	/// <summary> Print processors information </summary>
	void print();
};

/** Parameters for configured outliers cleaner
*/
//.......................................................................................
typedef struct {
	float logicalLow, logicalHigh, confirmedLow, confirmedHigh;
    string distLow,distHigh; //"none" "norm" or "log" 
} confRecord;

//.......................................................................................
/** RepConfiguredOutlierCleaner is a simple cleaner considering each value of a certain signal separatley,
* but this time use configuration file that holds for each signal the logical values, statistically confirmed 
* values  and distribution for relearning statistical values
*/
//.......................................................................................
class RepConfiguredOutlierCleaner : public  RepBasicOutlierCleaner {
public:

	string confFileName; ///< configuration file and mapping
	string cleanMethod; ///< cleaning method :  "logical" "confirmed" or "learned"
	map<string,confRecord> outlierParams; ///< a map from signal name to outliers parameters

	/// <summary> Initialize to default values </summary>
	void init_defaults() {
		processor_type = REP_PROCESS_CONFIGURED_OUTLIER_CLEANER;
		params.trimming_sd_num = DEF_REP_TRIMMING_SD_NUM; params.removing_sd_num = DEF_REP_REMOVING_SD_NUM; params.nbrs_sd_num = 0;
		params.take_log = 0;
		params.doTrim = params.doRemove = true;
		signalId = -1;
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};

	/// <summary> learn cleaning boundaries </summary>
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
		
	/// The parsed fields from init command.
	/// @snippet RepProcess.cpp RepConfiguredOutlierCleaner::init
	int init(map<string, string>& mapper);
	
	// Apply cleaning model -inheritted
	

	/// Serialization
	int version() { return 1; }
	ADD_SERIALIZATION_FUNCS(processor_type, signalName, time_channel, val_channel, req_signals, aff_signals, params.take_log, params.missing_value, params.doTrim, params.doRemove,
		trimMax, trimMin, removeMax, removeMin, confFileName, cleanMethod, outlierParams)

	void print();
};

void learnDistributionBorders(float& borderHi, float& borderLo,vector<float> filteredValues);
// a function that takes sorted vector of filtered values and estimates the +- 7 sd borders based on the center of distribution
// predefined calibration constants are used for estimation of the borders. 


/**
* A cleaner that is based on rules that describe relations of signal values to each other.\n
* This is a static cleaner ( no learning involved).\n
\n
Rules:\n
Rule1: BMI = Weight / Height ^ 2 * 1e4\n
Rule2:MCH = Hemoglobin / RBC\n
Rule3:MCV = Hematocrit / RBC\n
Rule4:MCHC - M = MCH / MCV\n
Rule5:Eosinophils# + Monocytes# + Basophils# + Lymphocytes# + Neutrophils# <= WBC\n
Rule6:MPV = Platelets_Hematocrit / Platelets\n
Rule7:UrineAlbumin <= UrineTotalProtein\n
Rule8:UrineAlbumin_over_Creatinine = UrineAlbumin / UrineCreatinine\n
Rule9:LDL + HDL <= Cholesterol\n
Rule10:NonHDLCholesterol + HDL = Cholesterol\n
Rule11:HDL_over_nonHDL = HDL / NonHDLCholesterol\n
Rule12:HDL_over_Cholesterol = HDL / Cholesterol\n
Rule13:HDL_over_LDL = HDL / LDL\n
Rule14:HDL_over_LDL = 1 / LDL_over_HDL\n
Rule15:Cholesterol_over_HDL = Cholesterol / HDL\n
Rule16:Cholesterol_over_HDL = 1 / HDL / Cholesterol\n
Rule17:Cholesterol_over_HDL = 1 / HDL_over_Cholestrol\n
Rule18:LDL_over_HDL = LDL / HDL\n
Rule19:Albumin <= Protein_Total\n
Rule20:FreeT4 <= T4\n
Rule21:NRBC <= RBC\n
Rule22:CHADS2_VASC >= CHADS2_VASC\n
*/
class RepRuleBasedOutlierCleaner : public RepProcessor, public MedValueCleaner {
	// get multiple signals and clean them according to consistency  with other signals from same date
public:

	
	/// Signals to clean
	vector <string> signalNames;
	vector <int> signalIds;
	int time_channel = 0;
	int val_channel = 0;
	MedDictionarySections myDict; ///< keeping it will enable us to get ids at apply stage
	bool addRequiredSignals=false; ///< a flag stating if we want to load signals that are not in the cleaned signal list 
								   /// because they share a rule with the cleaned signals (set it in jason)
	vector<int> consideredRules;///< only rules in this list will be considered in this cleaner (read list from jason)
	                            /// rule number 0 means apply all rules. Empty vector: do nothing in this cleaner.
	
	/// static map from rule to participating signals
	map <int, vector<string>>rules2Signals = { 
	{1,{"BMI","Weight","Height"}},
	{2,{"MCH", "Hemoglobin","RBC"}},
	{3,{"MCV","Hematocrit","RBC"} },
	{4,{"MCHC-M","MCH","MCV"}},
	{5,{"Eosinophils#","Monocytes#","Basophils#","Lymphocytes#","Neutrophils#","WBC" }},
	{6,{ "MPV","Platelets_Hematocrit","Platelets" }},
	{7,{"UrineAlbumin","UrineTotalProtein" }},
	{8,{"UrineAlbumin_over_Creatinine","UrineAlbumin","UrineCreatinine" }},
	{9,{"LDL","HDL","Cholesterol"}},
	{10,{"NonHDLCholesterol","HDL","Cholesterol"}},
	{11,{"HDL_over_nonHDL","HDL","NonHDLCholesterol"}},
	{12,{"HDL_over_Cholesterol","HDL","Cholesterol"}},
	{13,{"HDL_over_LDL","HDL","LDL"}},
	{14,{"HDL_over_LDL","LDL_over_HDL"}},
	{15,{"Cholesterol_over_HDL","Cholesterol","HDL"}},
	{17,{"Cholesterol_over_HDL","HDL_over_Cholestrol"}}, //rule 16 canceled
	{18,{"LDL_over_HDL","LDL","HDL"}},
	{19,{"Albumin","Protein_Total"}},
	{20,{"FreeT4","T4"}},
	{21,{"NRBC","RBC"}},
	{22,{"CHADS2","CHADS2_VASC"}}
	};

	vector <int> rulesToApply;

	// Constructors 
	RepRuleBasedOutlierCleaner() : RepProcessor() {init_defaults(); }
	
	void init_defaults() {
		processor_type = REP_PROCESS_RULEBASED_OUTLIER_CLEANER;
		
		params.take_log = 0;
		params.doRemove = true;
		
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};
	



	// Init
	int init(void *processor_params) { return MedValueCleaner::init(processor_params); };

	/// The parsed fields from init command.
	/// @snippet RepProcess.cpp RepRuleBasedOutlierCleaner::init
	int init(map<string, string>& mapper);
	
	
	// Learn cleaning model  : no learning for this cleaner. only apply
	
	///set signals
	void set_signal_ids(MedDictionarySections& dict) { myDict = dict; } // keep the dict. We will set ids later.

	/// Apply cleaning model 
	int _apply(PidDynamicRec& rec, vector<int>& time_points);

	// Serialization  -static not needed
	//print 
private:
	///ruleUsvs hold the signals in the order they appear in the rule in the rules2Signals above
	/// apply the rule and return true if data is consistent with the rule
	bool  applyRule(int rule, vector <UniversalSigVec> ruleUsvs,vector <int >sPointer); 
	
};

#define DEF_REP_NBRS_NBRS_SD_NUM 5
#define DEF_REP_NBRS_TRIM_SD_NUM 7
#define DEF_REP_NBRS_REMOVING_SD_NUM 14
//.......................................................................................
/** RepNbrsOutlierCleaner is cleaner that looks at the neighbourhood of a certain signal value
*/
//.......................................................................................
class RepNbrsOutlierCleaner : public RepProcessor, public MedValueCleaner {
public:

	string signalName; ///< name of signal to clean
	int signalId; ///< id of signal to clean
	int time_channel = 0; ///< time channel to consider in cleaning
	int val_channel = 0; ///< value cahnnel to consider in cleaning

	int nbr_time_width = 7; ///< size of neighborhood for defining neighboring values
	int nbr_time_unit = MedTime::Days; ///< time unit for defining neighboring values


	/// <summary> default constructor </summary>
	RepNbrsOutlierCleaner() { init_defaults(); }
	/// <summary> default constructor + setting signal name </summary>
	RepNbrsOutlierCleaner(const string& _signalName) { init_defaults(); signalId = -1; signalName = _signalName; init_lists(); }
	/// <summary> default constructor + setting signal name + initialize from string </summary>
	RepNbrsOutlierCleaner(const string& _signalName, string init_string) { init_defaults(); signalId = -1; signalName = _signalName; init_from_string(init_string); }
	/// <summary> default constructor + setting signal name + initialize from parameters </summary>
	RepNbrsOutlierCleaner(const string& _signalName, ValueCleanerParams *_params) { signalId = -1; signalName = _signalName; init_lists(); MedValueCleaner::init(_params); }

	/// <summary> Initialize to default values </summary>
	void init_defaults() {
		processor_type = REP_PROCESS_NBRS_OUTLIER_CLEANER;
		params.trimming_sd_num = DEF_REP_NBRS_TRIM_SD_NUM; params.removing_sd_num = DEF_REP_NBRS_REMOVING_SD_NUM; params.nbrs_sd_num = DEF_REP_NBRS_NBRS_SD_NUM;
		params.take_log = 0;
		params.doTrim = params.doRemove = true;
		signalId = -1;
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};

	/// <summary> Set signal name and fill affected and required signals sets </summary> 
	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName; init_lists(); }

	/// <summary> Set signal id </summary>
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }


	/// <summary> Fill required- and affected-signals sets </summary>
	int init(void *processor_params) { return MedValueCleaner::init(processor_params); };
	/// The parsed fields from init command.
	/// @snippet RepProcess.cpp RepNbrsOutlierCleaner::init
	int init(map<string, string>& mapper); 
	// { init_defaults();  return MedValueCleaner::init(mapper); };
	void init_lists();

	/// <summary> learn cleaning boundaries </summary>
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	/// <summary> Learning : learn cleaning boundaries using MedValueCleaner's iterative approximation of moments </summary>
	int iterativeLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	/// <summary> Learning : learn cleaning boundaries using MedValueCleaner's quantile approximation of moments </summary>
	int quantileLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);

	/// <summary> Apply cleaning model </summary>
	int _apply(PidDynamicRec& rec, vector<int>& time_points);

	// Serialization
	int version() { return 1; }
	ADD_SERIALIZATION_FUNCS(processor_type, signalName, time_channel, val_channel, req_signals, aff_signals, params.take_log, params.missing_value, params.doTrim, params.doRemove,
		trimMax, trimMin, removeMax, removeMin, nbr_time_unit, nbr_time_width, nbrsMax, nbrsMin)

	/// <summary> Print processors information </summary>
	void print();
};

//.......................................................................................
//.......................................................................................
// Utility Functions
//.......................................................................................
//.......................................................................................

/// <summary> Get values of a signal from a set of ids applying a set of preceeding cleaners </summary>
int get_values(MedRepository& rep, vector<int>& ids, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values, vector<RepProcessor *>& prev_cleaners);
/// <summary> Get values of a signal from a set of ids </summary>
int get_values(MedRepository& rep, vector<int>& ids, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values) ;

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(RepMultiProcessor)
MEDSERIALIZE_SUPPORT(RepBasicOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepRuleBasedOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepConfiguredOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepNbrsOutlierCleaner)

#endif
