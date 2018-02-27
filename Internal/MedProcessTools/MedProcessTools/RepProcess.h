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
	REP_PROCESS_CALC_SIGNALS,///<"calc_signals" or "calculator" to activate RepCalcSimpleSignals
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

	/// <summary> 
	/// virtual signals are created only in rep processors but can be used by any rep processor that comes after
	/// or any feture generator as a regular signal.
	/// a virtual signal can be defined as one only once - in the rep processor that actually creates it.
	/// The other rep processors and feature generators using it are simply blind to the fact that it is virtual 
	/// and use it by name as a regular signal.
	/// </summary>
	vector<pair<string, int>> virtual_signals; ///< list of all virtual signals CREATED by this rep processor

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

	/// <summary> rep processors CREATING virtual signals need to implement this: adding their signals to the pile </summary>
	virtual void add_virtual_signals(map<string, int> &_virtual_signals) { return; };
	
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
	inline bool is_signal_affected(int signalId) {return (aff_signal_ids.find(signalId) != aff_signal_ids.end());}
	inline bool is_signal_affected(string& signalName) { return (aff_signals.find(signalName) != aff_signals.end()); }

	// check filtering
	/// <summary> Check if processor (and 'sub'-processors within) should be applied according to set of required signals  </summray>
	/// <returns> true if processor is not required and can be filtered, false otherwise </returns>
	virtual bool filter(unordered_set<string>& reqSignals);

	/// <summary> Init required tables : Should be implemented for inheriting classes that have such tables </summary>
	virtual void init_tables(MedDictionarySections& dict) { return; }

	// Learning
	/// <summary> learn processing model on a subset of samples. Apply set of preceeding processors on DynamicPidRec before learning : 
	// Should be implemented for inheriting classes that require learning </summary>
	virtual int _learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors) {  return 0; };
	/// <summary> learn processing model on a subset of samples only if required. Apply set of preceeding processors on DynamicPidRec before learning : 
	// May be implemented for inheriting classes that require learning </summary>
	virtual int _conditional_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors, unordered_set<int>& neededSignalIds) ;

	// Learning envelopes - Here because of issues with overloading and inheritance
	/// <summary> learn processing model on a subset of ids. Apply set of preceeding processors on DynamicPidRec before learning </summary>
	int learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors) { return _learn(rep,samples,prev_processors); };
	/// <summary> learn on all pids in repository, using fake samples - works only for repProcessors that ignore sample dates</summary>
	int learn(MedPidRepository& rep);
	/// <summary> learn on subset of samples without preceesing processors  </summary>
	int learn(MedPidRepository& rep, MedSamples& samples) { vector<RepProcessor *> temp;  return _learn(rep, samples, temp); }
	/// <summary> learn processing model on a subset of samples only if required. Apply set of preceeding processors on DynamicPidRec before learning : 
	virtual int conditional_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors, unordered_set<int>& neededSignalIds) {
		return _conditional_learn(rep, samples, prev_processors, neededSignalIds); 
	}
	/// <summary> learn processing model on a subset of ids only if required without preceesing processors </summary>
	int conditional_learn(MedPidRepository& rep, MedSamples& samples, unordered_set<int>& neededSignalIds) { vector<RepProcessor *> temp;  return _conditional_learn(rep, samples, temp, neededSignalIds); }

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


	// debug prints
	/// used for debug prints, each inheriting class can overload this one to get a more precise debug print. rp_flag can be used to transfer verbosity levels.
	/// the default print just prints the basic type, etc.
	virtual void dprint(const string &pref, int rp_flag);


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
	/// <summary> Reporting back virtual signals if there are any </summary>
	void add_virtual_signals(map<string, int> &_virtual_signals);

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

	// check filtering
	/// <summary> Check if processor (and 'sub'-processors within) should be applied according to set of required signals  </summray>
	/// <returns> true if processor is not required and can be filtered, false otherwise </returns>
	bool filter(unordered_set<string>& reqSignals);

	/// <summary> Set signal-ids for all linked signals </summary>
	void set_signal_ids(MedDictionarySections& dict); 

	/// <summary> Init required tables : Should be implemented for inheriting classes that have such tables </summary>
	void init_tables(MedDictionarySections& dict) { for (RepProcessor * proc : processors) { proc->init_tables(dict); } }


	/// <summary> learn processors </summary>
	int _learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors);
	int _conditional_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors, unordered_set<int>& neededSignalIds);
	
	/// <summary> Apply processors </summary>
	int _apply(PidDynamicRec& rec, vector<int>& time_points);
	/// <summary> Apply processors that affect any of the needed signals </summary>
	int _conditional_apply(PidDynamicRec& rec, vector<int>& time_points, unordered_set<int>& neededSignals);

	/// debug prints
	void dprint(const string &pref, int rp_flag);

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
	int _learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processor);
	/// <summary> Learning : learn cleaning boundaries using MedValueCleaner's iterative approximation of moments </summary>
	int iterativeLearn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processor);
	/// <summary> Learning : learn cleaning boundaries using MedValueCleaner's quantile approximation of moments </summary>
	int quantileLearn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processor);

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
	int _learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processor);
		
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
Rule16:---------------------\n
Rule17:Cholesterol_over_HDL = 1 / HDL_over_Cholestrol\n
Rule18:LDL_over_HDL = LDL / HDL\n
Rule19:Albumin <= Protein_Total\n
Rule20:FreeT4 <= T4\n
Rule21:NRBC <= RBC\n
Rule22:CHADS2_VASC >= CHADS2\n
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
	int _learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processor);
	/// <summary> Learning : learn cleaning boundaries using MedValueCleaner's iterative approximation of moments </summary>
	int iterativeLearn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processor);
	/// <summary> Learning : learn cleaning boundaries using MedValueCleaner's quantile approximation of moments </summary>
	int quantileLearn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processor);

	/// <summary> Apply cleaning model </summary>
	int _apply(PidDynamicRec& rec, vector<int>& time_points);

	// Serialization
	int version() { return 1; }
	ADD_SERIALIZATION_FUNCS(processor_type, signalName, time_channel, val_channel, req_signals, aff_signals, params.take_log, params.missing_value, params.doTrim, params.doRemove,
		trimMax, trimMin, removeMax, removeMin, nbr_time_unit, nbr_time_width, nbrsMax, nbrsMin)

	/// <summary> Print processors information </summary>
	void print();
};


typedef enum {
	CALC_TYPE_UNDEF, ///< undefined signal, nothing to do or error
	CALC_TYPE_EGFR, ///< calculates eGFR CKD_EPI signal in each point we have Creatinine , depends on Creatinine, GENDER, BYEAR
	CALC_TYPE_DEBUG, ///< just here to help when debugging: currently calculates delta Hemoglobin (relative to last test for each test above second)
	CALC_TYPE_HOSP_MELD, ///< calculates MELD signal. components: CHARTLAB_Bilirubin, CHARTLAB_INR(PT), CHARTLAB_Creatinine, CHARTLAB_Sodium
	CALC_TYPE_HOSP_BP_SYS,
	CALC_TYPE_HOSP_BP_DIA,
	CALC_TYPE_HOSP_BMI,
	CALC_TYPE_HOSP_APRI,
	CALC_TYPE_HOSP_SIDA,
	CALC_TYPE_HOSP_PaO2_FiO2_RATIO,
	CALC_TYPE_HOSP_IS_AFRICAN_AMERICAN,
	CALC_TYPE_HOSP_IS_MECHANICALLY_VENTILATED,
	CALC_TYPE_HOSP_SOFA_NERVOUS,
	CALC_TYPE_HOSP_SOFA_LIVER,
	CALC_TYPE_HOSP_24H_URINE_OUTPUT,
	CALC_TYPE_HOSP_SOFA_COAGULATION,
	CALC_TYPE_HOSP_DOPAMINE_PER_KG,
	CALC_TYPE_HOSP_EPINEPHRINE_PER_KG,
	CALC_TYPE_HOSP_NOREPINEPHRINE_PER_KG,
	CALC_TYPE_HOSP_DOBUTAMINE_PER_KG,
	CALC_TYPE_HOSP_QSOFA,
	CALC_TYPE_HOSP_SIRS,
	CALC_TYPE_HOSP_PRESSURE_ADJUSTED_HR,
	CALC_TYPE_HOSP_MODS,
	CALC_TYPE_HOSP_SHOCK_INDEX,
	CALC_TYPE_HOSP_PULSE_PRESSURE,
	CALC_TYPE_HOSP_EFGR,
	CALC_TYPE_HOSP_SOFA_RESPIRATORY,
	CALC_TYPE_HOSP_SOFA_RENAL,
	CALC_TYPE_HOSP_SOFA_CARDIO,
	CALC_TYPE_HOSP_SOFA
} RepCalcSimpleSignalsType;

//.......................................................................................
/** RepCalcSimpleSignals is a rep processor containing several calculators to calculate new
    signals. It supports an internal list of calculators, and the user can select one of them.

	Each calculator can create one or more new virtual signals. Each of these has a specific type.
	The user can configure the virtual signal names, and can also pass a list of float parameters to
	the calculator, making it parametric.

	When adding a new calculator make sure to :
	(1) fill in the calc2defs map (as explained below), to define the calculator name, type, req signals, 
		default virual names and their types, and default parameters.
	(2) add the new calculator type to the enum
	(3) write the matching _apply function for the specific calculator and make sure the general _apply calls it .

	NOTE: in RepCalcSimpleSignals there is NO learning. Calculated signals that have a learning stage should
	      be implmented in a separate class that will also have space to keep the learning process results.

	supported signals calculated (this will be also the virtual signal name):

	calc_eGFR : calculating EGFR (CKD_EPI formula) at each point in which Creatinine is available

*/
//.......................................................................................
class RepCalcSimpleSignals : public RepProcessor {

	public:
		vector<string> V_names; ///< names of signals created by the calculator (a calculator can create more than a single signal at a time)
		vector<int> V_types; ///< matching signal types for V_names

		string calculator; ///< calculator asked for by user
		int calc_type;	///> type of calculator (one of the allowed enum values)

		float missing_value = (float)MED_MAT_MISSING_VALUE;

		vector<float> coeff; ///< it is possible to transfer a vector of params to the calculator, to enable parametric calculators.

		RepCalcSimpleSignals() { processor_type = REP_PROCESS_CALC_SIGNALS; }

		/// <summary> initialize from a map :  Should be implemented for inheriting classes that have parameters </summary>
		int init(map<string, string>& mapper);

		// making sure V_ids and sigs_ids are initialized
		void init_tables(MedDictionarySections& dict);

		void add_virtual_signals(map<string, int> &_virtual_signals);


		// Learning
		/// <summary> In this class there's never learning - we return 0 immediately </summary>
		int _learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors) { init_tables(rep.dict); return 0; };

		// Applying
		/// <summary> apply processing on a single PidDynamicRec at a set of time-points : Should be implemented for all inheriting classes </summary>
		int _apply(PidDynamicRec& rec, vector<int>& time_points);

		// Utils
		/// <summary> translate a calculator name to RepCalcSignalsType </summary>
		int get_calculator_type(const string &calc_name);

		// calculators implementations + helpers
		static float get_age(int byear, int date);
		int _apply_calc_eGFR(PidDynamicRec& rec, vector<int>& time_points);
		int _apply_calc_debug(PidDynamicRec& rec, vector<int>& time_points);
		static float calc_egfr_ckd_epi(float creatinine, int gender, float age, int ethnicity = 0);

		//process io fields in MIMIC style, depending on their type (understood from the name and the type
		static int process_hosp_signal(const string& name, UniversalSigVec& usv, vector<int>& times, vector<float>& vals);
				
		int _apply_calc_hosp_pointwise(PidDynamicRec& rec, vector<int>& time_points, float (*calcFunc)(const vector<float>&, const vector<float>&));
		int _apply_calc_hosp_time_dependent_pointwise(PidDynamicRec& rec, vector<int>& time_points,
			float(*calcFunc)(const vector<pair<int, float> >&, int, const vector<float>&));
		int _apply_calc_24h_urine_output(PidDynamicRec& rec, vector<int>& time_points);

		static float calc_hosp_MELD(float bilirubin, float inr, float creatinine, float na, int mode);	
		static float calc_hosp_MELD(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_MELD(args.at(0), args.at(1), args.at(2), args.at(3), (int)(params.at(0)));
		}
		//no need for a function for combined blood pressure
		static float calc_hosp_BMI(float weight, float height);
		static float calc_hosp_BMI(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_BMI(args.at(0), args.at(1));
		}
		static float calc_hosp_APRI(float ast, float plt);
		static float calc_hosp_APRI(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_APRI(args.at(0), args.at(1));
		}
		static float calc_hosp_SIDA(float na, float k, float cl, int mode);
		static float calc_hosp_SIDA(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SIDA(args.at(0), args.at(1), args.at(2), (int)(params.at(0)));
		}
		static float calc_hosp_PaO2_FiO2_ratio(float paO2, float fiO2);
		static float calc_hosp_PaO2_FiO2_ratio(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_PaO2_FiO2_ratio(args.at(0), args.at(1));
		}
		static float calc_hosp_is_african_american(float ethnicity);
		static float calc_hosp_is_african_american(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_is_african_american(args.at(0));
		}
		//no need for a function for mechanical ventilation
		static float calc_hosp_SOFA_nervous(float x, bool optimistic); 
		static float calc_hosp_SOFA_nervous(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SOFA_nervous(args.at(0), params.at(0) != 0);
		}
		static float calc_hosp_SOFA_liver(float x, bool optimistic);
		static float calc_hosp_SOFA_liver(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SOFA_liver(args.at(0), params.at(0) != 0);
		}
		//no need for a function for 24h urine output
		static float calc_hosp_SOFA_coagulation(float x, bool optimistic);
		static float calc_hosp_SOFA_coagulation(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SOFA_coagulation(args.at(0), params.at(0) != 0);
		}
		static float calc_hosp_dopamine_per_kg(float medPerKg);
		static float calc_hosp_dopamine_per_kg(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_dopamine_per_kg(args.at(0));
		}
		static float calc_hosp_epinephrine_per_kg(float medPerKg, float med, float kg);
		static float calc_hosp_epinephrine_per_kg(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_epinephrine_per_kg(args.at(0), args.at(1), args.at(2));
		}
		static float calc_hosp_norepinephrine_per_kg(float medPerKg, float med, float kg);
		static float calc_hosp_norepinephrine_per_kg(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_norepinephrine_per_kg(args.at(0), args.at(1), args.at(2));
		}
		static float calc_hosp_dobutamine_per_kg(float medPerKg);
		static float calc_hosp_dobutamine_per_kg(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_dobutamine_per_kg(args.at(0));
		}
		static float calc_hosp_SIRS(float temp, float hr, float resp, float paco2, float wbc, float bands, int mode);
		static float calc_hosp_SIRS(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SIRS(args.at(0), args.at(1), args.at(2), args.at(3), args.at(4), args.at(5), (int)(params.at(0)));
		}
		static float calc_hosp_pressure_adjusted_hr(float hr, float cvp, float map);
		static float calc_hosp_pressure_adjusted_hr(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_pressure_adjusted_hr(args.at(0), args.at(1), args.at(2));
		}
		static float calc_hosp_MODS(float paO2, float fiO2, float plt, float bili, float hr, float cvp, float map, float gcs, float cre);
		static float calc_hosp_MODS(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_MODS(args.at(0), args.at(1), args.at(2), args.at(3), args.at(4), args.at(5), args.at(6), args.at(7), args.at(8));
		}
		static float calc_hosp_shock_index(float hr, float sysBp, float diaBp, int mode);
		static float calc_hosp_shock_index(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_shock_index(args.at(0), args.at(1), args.at(2), (int)(params.at(0)));
		}
		static float calc_hosp_pulse_pressure(float sysBp, float diaBp);
		static float calc_hosp_pulse_pressure(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_pulse_pressure(args.at(0), args.at(1));
		}
		static float calc_hosp_eGFR(float cr, float age, float gender, float isAfricanAmerican, bool useEthnicity, int mode);
		static float calc_hosp_eGFR(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_eGFR(args.at(0), args.at(1), args.at(2), args.at(3), params.at(0) != 0, (int)(params.at(0)));
		}
		static float calc_hosp_SOFA_respiratory(float x, float y, bool optimistic);
		static float calc_hosp_SOFA_respiratory(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SOFA_respiratory(args.at(0), args.at(1), params.at(0) != 0);
		}
		static float calc_hosp_SOFA_renal(float x, float y, bool optimistic);
		static float calc_hosp_SOFA_renal(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SOFA_renal(args.at(0), args.at(1), params.at(0) != 0);
		}
		static float calc_hosp_SOFA_cardio(float x, float y, float z, float u, float v, bool optimistic);
		static float calc_hosp_SOFA_cardio(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SOFA_cardio(args.at(0), args.at(1), args.at(2), args.at(3), args.at(4), params.at(0) != 0);
		}

		static float calc_hosp_SOFA(float snerve, float sliver, float scoag, float sresp, float srenal, float scardio, int mode);
		static float calc_hosp_SOFA(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_SOFA(args.at(0), args.at(1), args.at(2), args.at(3), args.at(4), args.at(5), (int)(params.at(0)));
		}

		static float calc_hosp_qSOFA(float gcs, float sBp, float resp, int mode);
		static float calc_hosp_qSOFA(const vector<float>& args, const vector<float>& params) {
			return calc_hosp_qSOFA(args.at(0), args.at(1), args.at(2), (int)(params.at(0)));
		}

		static float med2KgPlusMed(float med2kg, float med, float kg);
		static float anySeenRecently(const vector<pair<int, float> >& data, int time, const vector<float>& params);
		static float interleave(const vector<pair<int, float> >& data, int time, const vector<float>& params);
				
		static bool isNegative(float val) { return (val < 0.0F ? true : false); }
		static bool isNonPositive(float val) { return (val <= 0.0F ? true : false); }
		static bool isAny(const vector<float>& vals, bool (*pred)(float)) {
			for (const auto& e : vals) {
				if (pred(e))
					return true;
			}
			return false;
		}
	
		static float badValue(void) { return (float)MED_MAT_MISSING_VALUE;}

		static bool isMissingValue(float val) { return val == (float)MED_MAT_MISSING_VALUE; }
		static bool isMissingOrNegative(float val) { return ((val < 0.0F || isMissingValue(val)) ? true : false); }
		static bool isMissingOrNonPositive(float val) { return ((val <= 0.0F || isMissingValue(val)) ? true : false); }
		static int beforeEverything() { return numeric_limits<int>::min(); }
		static int afterEverything() { return numeric_limits<int>::max(); }

		// serialization
		ADD_SERIALIZATION_FUNCS(calculator, calc_type, coeff, V_names, V_types, req_signals, aff_signals, virtual_signals)


	private:

		// definitions and defaults for each calculator - all must be filled in for a new calculator

		/// from a calculator name to a calculator enum type
		const map<string, int> calc2type = { 
			{"calc_eGFR", CALC_TYPE_EGFR}, 
			{"calc_debug", CALC_TYPE_DEBUG},
			{"calc_hosp_MELD", CALC_TYPE_HOSP_MELD},
			{ "calc_hosp_BP_sys", CALC_TYPE_HOSP_BP_SYS},
			{ "calc_hosp_BP_dia", CALC_TYPE_HOSP_BP_DIA},
			{ "calc_hosp_BMI", CALC_TYPE_HOSP_BMI },
			{ "calc_hosp_APRI", CALC_TYPE_HOSP_APRI },
			{ "calc_hosp_SIDA", CALC_TYPE_HOSP_SIDA },
			{ "calc_hosp_PaO2_FiO2_ratio", CALC_TYPE_HOSP_PaO2_FiO2_RATIO },	
			{ "calc_hosp_is_african_american", CALC_TYPE_HOSP_IS_AFRICAN_AMERICAN },
			{ "calc_hosp_is_mechanically_ventilated", CALC_TYPE_HOSP_IS_MECHANICALLY_VENTILATED },
			{ "calc_hosp_SOFA_nervous", CALC_TYPE_HOSP_SOFA_NERVOUS },
			{ "calc_hosp_SOFA_liver", CALC_TYPE_HOSP_SOFA_LIVER },
			{ "calc_hosp_24h_urine_output",CALC_TYPE_HOSP_24H_URINE_OUTPUT },
			{ "calc_hosp_SOFA_coagulation", CALC_TYPE_HOSP_SOFA_COAGULATION },
			{ "calc_hosp_dopamine_per_kg", CALC_TYPE_HOSP_DOPAMINE_PER_KG },
			{ "calc_hosp_epinephrine_per_kg",CALC_TYPE_HOSP_EPINEPHRINE_PER_KG },
			{ "calc_hosp_norepinephrine_per_kg", CALC_TYPE_HOSP_NOREPINEPHRINE_PER_KG },
			{ "calc_hosp_dobutamine_per_kg", CALC_TYPE_HOSP_DOBUTAMINE_PER_KG },
			{ "calc_hosp_qSOFA", CALC_TYPE_HOSP_QSOFA },
			{ "calc_hosp_SIRS", CALC_TYPE_HOSP_SIRS },
			{ "calc_hosp_pressure_adjusted_hr", CALC_TYPE_HOSP_PRESSURE_ADJUSTED_HR },
			{ "calc_hosp_MODS", CALC_TYPE_HOSP_MODS },
			{ "calc_hosp_shock_index", CALC_TYPE_HOSP_SHOCK_INDEX },
			{ "calc_hosp_pulse_pressure", CALC_TYPE_HOSP_PULSE_PRESSURE },
			{ "calc_hosp_eGFR", CALC_TYPE_HOSP_EFGR },
			{ "calc_hosp_SOFA_respiratory", CALC_TYPE_HOSP_SOFA_RESPIRATORY },
			{ "calc_hosp_SOFA_renal", CALC_TYPE_HOSP_SOFA_RENAL },
			{ "calc_hosp_SOFA_cardio", CALC_TYPE_HOSP_SOFA_CARDIO },	
			{ "calc_hosp_SOFA", CALC_TYPE_HOSP_SOFA }
		};

		/// from a calculator name to the list of required signals
		const map<string, vector<string>> calc2req_sigs = { 
			//--------- level 1 - calculated from raw signals (level0)
			{"calc_eGFR", {"Creatinine", "GENDER", "BYEAR"}}, 
			{ "calc_debug",{ "Hemoglobin" }},
			{"calc_hosp_MELD",{ "CHARTLAB_Bilirubin", "CHARTLAB_INR(PT)", "CHARTLAB_Creatinine", "CHARTLAB_Sodium" }},
			{"calc_hosp_BP_sys", {"CHART_Arterial_BP_Sys", "CHART_NonInvasive_BP_Sys"}},
			{"calc_hosp_BP_dia", {"CHART_Arterial_BP_Dia", "CHART_NonInvasive_BP_Dia"}},
			{ "calc_hosp_BMI",{ "CHART_Weight", "CHART_Height" } },
			{ "calc_hosp_APRI",{ "LAB_AST", "CHARTLAB_Platelets" } },
			{ "calc_hosp_SIDA",{ "CHARTLAB_Sodium","CHARTLAB_Potassium","CHARTLAB_Chloride" } },
			{ "calc_hosp_PaO2_FiO2_ratio",{ "CHART_ART_PaO2","CHART_FiO2" } },
			{ "calc_hosp_is_african_american",{ "ETHNICITY" } },
			{ "calc_hosp_is_mechanically_ventilated",{ "CHART_Airway_Size","CHART_Peak_Insp_Pressure","CHART_Plateau_Pressure",
														"CHART_Plateau_Off", "CHART_Tidal_Volume_Set", "CHART_Vent_Mode", "CHART_Vent_No",
														"CHART_Vent_Type", "CHART_Resp_Rate_Set" } },
			{ "calc_hosp_SOFA_nervous",{ "CHART_GCS" } },
			{ "calc_hosp_SOFA_liver",{ "CHARTLAB_Bilirubin" } },
			//need to implement "calc_24h_mean_urine_output" with special care - see my module
			{ "calc_hosp_24h_urine_output",{ "OUTPUT_UrineProduction" } },
			{ "calc_hosp_SOFA_coagulation",{ "CHARTLAB_Platelets" } },
			{ "calc_hosp_dopamine_per_kg",{"INPUT_Dopamine-k_Rate"} },
			{ "calc_hosp_epinephrine_per_kg",{"INPUT_Epinephrine-k_Rate","INPUT_Epinephrine_Rate","CHART_Weight" } },
			{ "calc_hosp_norepinephrine_per_kg",{ "INPUT_Norepinephrine-k_Rate","INPUT_Norepinephrine_Rate","CHART_Weight" } },
			{ "calc_hosp_dobutamine_per_kg",{ "INPUT_Dobutamine-k_Rate" } },
			{ "calc_hosp_SIRS",{ "CHART_Temperature","CHART_Heart_Rate","CHART_Resp_Rate","CHART_Art_PaCO2",
									"CHARTLAB_WBC","LAB_Neutrophils_Bands%" } },
			{ "calc_hosp_pressure_adjusted_hr",{ "CHART_Heart_Rate","CHART_CVP","CHART_Arterial_BP_Mean" } },
			{ "calc_hosp_MODS",{ "CHART_ART_PaO2","CHART_FiO2","CHARTLAB_Platelets","CHARTLAB_Bilirubin","CHART_Heart_Rate","CHART_CVP",
									"CHART_Arterial_BP_Mean","CHART_GCS","CHARTLAB_Creatinine" } },			
			//---------- level 2 - calculated from level < 2
			{"calc_hosp_shock_index", { "CHART_Heart_Rate", "calc_hosp_BP_sys", "calc_hosp_BP_dia" }},							
			{ "calc_hosp_pulse_pressure", { "calc_hosp_BP_sys", "calc_hosp_BP_dia" } },			
			{ "calc_hosp_eGFR", { "CHARTLAB_Creatinine", "Age", "Gender", "calc_hosp_is_african_american" } },
			{ "calc_hosp_SOFA_respiratory", { "calc_hosp_PaO2_FiO2_ratio", "calc_hosp_is_mechanically_ventilated" } },
			{ "calc_hosp_SOFA_renal", {"CHARTLAB_Creatinine", "calc_hosp_24h_urine_output" } },
			{ "calc_hosp_SOFA_cardio", {"CHART_Arterial_BP_Mean", "calc_hosp_dopamine_per_kg", "calc_hosp_epinephrine_per_kg", 
										"calc_hosp_norepinephrine_per_kg", "calc_hosp_dobutamine_per_kg"} },
			{ "calc_hosp_qSOFA", {"CHART_GCS", "calc_hosp_BP_sys", "CHART_Resp_Rate"} },	
			//---------- level 3 - calculated from level < 3
		    { "calc_hosp_SOFA",{ "calc_hosp_SOFA_nervous", "calc_hosp_SOFA_liver", "calc_hosp_SOFA_coagulation", "calc_hosp_SOFA_respiratory",
								 "calc_hosp_SOFA_renal", "calc_hosp_SOFA_cardio" } }
		};

		/// from a calculator name to a list of pairs of virtual names and their types created by the calculator
		/// the virtual names can be changed by the user (but have to be given in the SAME order as here)
		const map<string, vector<pair<string, int>>> calc2virtual = {
			{ "calc_eGFR" ,{ { "calc_eGFR", T_DateVal } } },
			{ "calc_debug" ,{ { "calc_debug", T_DateVal } } },
			{ "calc_hosp_MELD" ,{ { "calc_hosp_MELD", T_TimeVal } } },
			{ "calc_hosp_BP_sys",{ { "calc_hosp_BP_sys", T_TimeVal } } },
			{ "calc_hosp_BP_dia",{ { "calc_hosp_BP_dia", T_TimeVal } } },		
			{ "calc_hosp_BMI",{ { "calc_hosp_BMI", T_TimeVal } } },
			{ "calc_hosp_APRI",{ { "calc_hosp_APRI", T_TimeVal } } },
			{ "calc_hosp_SIDA",{ { "calc_hosp_SIDA", T_TimeVal } } },
			{ "calc_hosp_PaO2_FiO2_ratio",{ { "calc_hosp_PaO2_FiO2_ratio", T_TimeVal } } },		
			{ "calc_hosp_is_african_american",{ { "calc_hosp_is_african_american", T_TimeVal } } },		
			{ "calc_hosp_is_mechanically_ventilated", { { "calc_hosp_is_mechanically_ventilated", T_TimeVal } } },		
			{ "calc_hosp_SOFA_nervous",{ { "calc_hosp_SOFA_nervous", T_TimeVal } } },
			{ "calc_hosp_SOFA_liver",{ { "calc_hosp_SOFA_liver", T_TimeVal } } },
			{ "calc_hosp_24h_urine_output",{ { "calc_hosp_24h_urine_output", T_TimeVal } } },		
			{ "calc_hosp_SOFA_coagulation",{ { "calc_hosp_SOFA_coagulation", T_TimeVal } } },
			{ "calc_hosp_dopamine_per_kg",{ { "calc_hosp_dopamine_per_kg", T_TimeVal } } },
			{ "calc_hosp_epinephrine_per_kg",{ { "calc_hosp_epinephrine_per_kg", T_TimeVal } } },
			{ "calc_hosp_norepinephrine_per_kg",{ { "calc_hosp_norepinephrine_per_kg", T_TimeVal } } },
			{ "calc_hosp_dobutamine_per_kg",{ { "calc_hosp_dobutamine_per_kg", T_TimeVal } } },	
			{ "calc_hosp_qSOFA",{ { "calc_hosp_qSOFA", T_TimeVal } } },
			{ "calc_hosp_SIRS",{ { "calc_hosp_SIRS", T_TimeVal } } },
			{ "calc_hosp_pressure_adjusted_hr",{ { "calc_hosp_pressure_adjusted_hr", T_TimeVal } } },
			{ "calc_hosp_MODS",{ { "calc_hosp_MODS", T_TimeVal } } },
			{ "calc_hosp_shock_index",{ { "calc_hosp_shock_index", T_TimeVal } } },
			{ "calc_hosp_pulse_pressure",{ { "calc_hosp_pulse_pressure", T_TimeVal } } },
			{ "calc_hosp_eGFR",{ { "calc_hosp_eGFR", T_TimeVal } } },
			{ "calc_hosp_SOFA_respiratory",{ { "calc_hosp_SOFA_respiratory", T_TimeVal } } },
			{ "calc_hosp_SOFA_renal",{ { "calc_hosp_SOFA_renal", T_TimeVal } } },
			{ "calc_hosp_SOFA_cardio",{ { "calc_hosp_SOFA_cardio", T_TimeVal } } },
		    { "calc_hosp_SOFA",{ { "calc_hosp_SOFA", T_TimeVal } } }
		};

		/// from a calculator name to the default coefficients (parameters) of it. Can of course be empty for a non parametric calculator.
		const map<string, vector<float>> calc2coeffs = {
			{"calc_eGFR" , {}},
			{ "calc_debug" ,{}},
			{"calc_hosp_MELD" ,{1.0F}},
			{ "calc_hosp_BP_sys",{} },
			{ "calc_hosp_BP_dia",{} },		
			{ "calc_hosp_BMI",{} },
			{ "calc_hosp_APRI",{} },
			{ "calc_hosp_SIDA",{1.0F} },
			{ "calc_hosp_PaO2_FiO2_ratio",{} },
			{ "calc_hosp_is_african_american", {} },		
			{ "calc_hosp_is_mechanically_ventilated", {300.0F} },
			{ "calc_hosp_SOFA_nervous",{ 1.0F } },
			{ "calc_hosp_SOFA_liver",{ 1.0F } },
			{ "calc_hosp_24h_urine_output", {} },
			{ "calc_hosp_SOFA_coagulation",{ 1.0F } },
			{ "calc_hosp_dopamine_per_kg",{} },
			{ "calc_hosp_epinephrine_per_kg",{} },
			{ "calc_hosp_norepinephrine_per_kg",{} },
			{ "calc_hosp_dobutamine_per_kg",{} },
			{ "calc_hosp_qSOFA",{ 1.0F } },
			{ "calc_hosp_SIRS",{ 2.0F } },
			{ "calc_hosp_pressure_adjusted_hr",{} },
			{ "calc_hosp_MODS",{} },
			{ "calc_hosp_shock_index",{ 1.0F } },
			{ "calc_hosp_pulse_pressure",{} },
			{ "calc_hosp_eGFR",{ 1.0F } },
			{ "calc_hosp_SOFA_respiratory",{ 1.0F } },
			{ "calc_hosp_SOFA_renal",{ 1.0F } },
			{ "calc_hosp_SOFA_cardio",{ 1.0F } },
		    { "calc_hosp_SOFA",{ 1.0F } }
		};

		vector<int> V_ids; ///< ids of signals created by the calculator (for faster usage at run time: save name conversions)
		vector<int> sigs_ids; /// <ids of signals used as input by the calculator (for faster usage at run time: save name conversions)

};


//.......................................................................................
//.......................................................................................
// Utility Functions
//.......................................................................................
//.......................................................................................

/// <summary> Get values of a signal from a set of samples applying a set of preceeding cleaners </summary>
int get_values(MedRepository& rep, MedSamples& samples, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values,
	vector<RepProcessor *>& prev_cleaners);
/// <summary> Get values of a signal from a set of samples </summary>
int get_values(MedRepository& rep, MedSamples& samples, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values) ;

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(RepMultiProcessor)
MEDSERIALIZE_SUPPORT(RepBasicOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepRuleBasedOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepConfiguredOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepNbrsOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepCalcSimpleSignals)

#endif
