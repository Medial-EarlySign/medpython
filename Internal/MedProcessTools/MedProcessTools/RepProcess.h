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

	/// Type
	RepProcessorTypes processor_type;

	/// Threading (learn/apply) - 
	int learn_nthreads, apply_nthreads;

	/// Required Signals (and their ids) for processing
	vector<string> req_signals;
	vector<int> req_signal_ids;

	/// Signals (and their ids) affected by processing 
	unordered_set<string> aff_signals;
	unordered_set<int> aff_signal_ids;

	// Constructor/Destructor
	RepProcessor() { learn_nthreads = DEFAULT_REP_CLNR_NTHREADS; apply_nthreads = DEFAULT_REP_CLNR_NTHREADS; };
	~RepProcessor() {};

	/// Init from name or type. optinally with a parameters string
	static RepProcessor *make_processor(string name);
	static RepProcessor *make_processor(string type, string params);
	static RepProcessor *make_processor(RepProcessorTypes type);
	static RepProcessor *make_processor(RepProcessorTypes type, string params);

	/// Create processor from params string (type must be given within string)
	static RepProcessor *create_processor(string &params);

	/// initialize : from object/string/defaults.
	/// Should be implemented for inheriting classes that have parameters
	virtual int init(void *params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	/// Set signalName
	/// Should be implemented for inheriting classes that have signalName
	virtual void set_signal(const string& _signalName) { return; };

	/// Set signalId
	/// Should be implemented for inheriting classes that have signalId
	virtual void set_signal_ids(MedDictionarySections& dict) { return; }

	/// Required Signals functions : get all signals that are required by the processor
	/// Append required signal names to vector : parent function just uses req_signals
	virtual void get_required_signal_names(unordered_set<string>& signalNames);
	/// Fill req_signal_ids : parent function just fills from req_signals
	virtual void set_required_signal_ids(MedDictionarySections& dict);
	// Append required signal ids to vector. call set_required_signal_ids if req_signal_ids empty
	void get_required_signal_ids(unordered_set<int>& signalIds, MedDictionarySections& dict);
	
	/// Affected Signals functions;
	/// Fill aff_signal_ids : parent function just fills from aff_signals
	virtual void set_affected_signal_ids(MedDictionarySections& dict);
	/// Check if a signal is affected by processor
	bool is_signal_affected(int signalId) {return (aff_signal_ids.find(signalId) != aff_signal_ids.end());}

	/// Init required tables
	/// Should be implemented for inheriting classes that have such tables
	virtual void init_tables(MedDictionarySections& dict) { return; }

	/// Learn processing model on a subset of ids. Apply set of preceesing processors on DynamicPidRec before learning
	/// Should be implemented for inheriting classes that require learning
	virtual int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) {  return 0; };

	/// Envelope learning functions
	/// Learn on subset of ids
	int learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) { return Learn(rep, ids, prev_processors); }
	/// Learn on all ids in repository
	int learn(MedPidRepository& rep, vector<RepProcessor *>& prev_processors) { return Learn(rep, rep.pids, prev_processors); }
	/// Learn on subset of ids without preceesing processors
	int learn(MedPidRepository& rep, vector<int>& ids) { vector<RepProcessor *> temp;  return Learn(rep, ids, temp); }
	/// Learn on all ids in repository without preceesing processors
	int learn(MedPidRepository& rep) { vector<RepProcessor *> temp; return learn(rep, temp); }

	/// Apply processing on a single PidDynamicRec at a set of time-points
	/// Should be implemented for all inheriting classes
	virtual int apply(PidDynamicRec& rec, vector<int>& time_points) {return 0; }
	/// Apply processing on a single PidDynamicRec at a set of time-points only if required : if any of the signals in neededSignalIds is actually affected by processor
	int apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignalIds);
	/// Apply processing on a single PidDynamicRec at a set of time-points given by Samples
	int apply(PidDynamicRec& rec, MedIdSamples& samples);
	/// Apply processing on a single PidDynamicRec at a set of time-points given by Samples only if required
	int apply(PidDynamicRec& rec, MedIdSamples& samples, vector<int>& neededSignalIds);	

	/// Serialization (including type)
	size_t get_processor_size();
	size_t processor_serialize(unsigned char *blob);

	/// optional printing of processor
	virtual void print() { fprintf(stderr, "No implementation for print()\n"); }
};

// Utilities
RepProcessorTypes rep_processor_name_to_type(const string& procesor_name);

/**
* A Processor which contains a vector of simpler processors that can be learned/applied
* in parallel. Useful for applying same cleaners on a set of signals, for example
*/
class RepMultiProcessor : public RepProcessor {
public:
	/// Set of processors
	vector<RepProcessor *> processors;

	// Constructor/Destructor
	RepMultiProcessor() { processor_type = REP_PROCESS_MULTI; };
	~RepMultiProcessor() {};

	/// Add processors (with initialization string)
	void add_processors_set(RepProcessorTypes type, vector<string>& signals);
	void add_processors_set(RepProcessorTypes type, vector<string>& signals, string init_string);

	/// Required Signals ids : Fill the member vector - req_signal_ids
	void set_required_signal_ids(MedDictionarySections& dict);

	/// Required Signals names : Fill the unordered set signalNames
	void get_required_signal_names(unordered_set<string>& signalNames);

	/// Affected Signals : Fill the member vector aff_signal_ids
	void set_affected_signal_ids(MedDictionarySections& dict);

	/// Set signal-ids for all linked signals
	void set_signal_ids(MedDictionarySections& dict);

	/// Learn processors
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors);

	/// Apply processors
	int apply(PidDynamicRec& rec, vector<int>& time_points);
	/// Apply processors that affect any of the needed signals
	int apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignals);

	/// serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	/// Print processors information
	void print() { for (auto& processor : processors) processor->print(); }
};

#define DEF_REP_TRIMMING_SD_NUM 7
#define DEF_REP_REMOVING_SD_NUM 14

/**
* A simple cleaner considering each value of a certain signal separatley
*/
class RepBasicOutlierCleaner : public RepProcessor, public MedValueCleaner {
public:

	/// Signal to clean
	string signalName;
	int signalId;
	int time_channel = 0;
	int val_channel = 0;


	// Constructors
	RepBasicOutlierCleaner() : RepProcessor() { init_defaults(); }
	RepBasicOutlierCleaner(const string& _signalName) : RepProcessor() { init_defaults(); signalId = -1; signalName = _signalName; init_lists(); }
	RepBasicOutlierCleaner(const string& _signalName, string init_string) : RepProcessor() { init_defaults(); signalId = -1; signalName = _signalName; init_from_string(init_string); }
	RepBasicOutlierCleaner(const string& _signalName, ValueCleanerParams *_params) : RepProcessor() {signalId = -1; signalName = _signalName; init_lists() ; MedValueCleaner::init(_params);}

	/// Init to default values
	void init_defaults() {
		processor_type = REP_PROCESS_BASIC_OUTLIER_CLEANER;
		params.trimming_sd_num = DEF_REP_TRIMMING_SD_NUM; params.removing_sd_num = DEF_REP_REMOVING_SD_NUM; params.nbrs_sd_num = 0;
		params.take_log = 0;
		params.doTrim = params.doRemove = true;
		signalId = -1;
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};

	/// Set signal name and fill aff/req signals 
	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName; init_lists(); }

	/// set signal id
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }

	/// Init from params
	int init(void *processor_params) { return MedValueCleaner::init(processor_params); };
	/// Init from map
	virtual int init(map<string, string>& mapper);
	/// Fill req- and aff-signals vectors
	void init_lists();

	/// Learn cleaning boundaries
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	/// Learning : learn cleaning boundaries using MedValueCleaner's iterative approximation of moments
	int iterativeLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	/// Learning : learn cleaning boundaries using MedValueCleaner's quantile approximation of moments
	int quantileLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);

	/// Apply cleaning model
	int apply(PidDynamicRec& rec, vector<int>& time_points);

	/// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void print();
};

typedef struct {
	float logicalLow, logicalHigh, confirmedLow, confirmedHigh;
    string distLow,distHigh; //"none" "norm" or "log" 
}confRecord;

/**
* A simple cleaner considering each value of a certain signal separatley, but this time use
* configuration file that holds for each signal the logical values, statistically confirmed values
* and distribution for relearning statistical values
*/
class RepConfiguredOutlierCleaner : public  RepBasicOutlierCleaner {
public:

	// Signal to clean -- inheritted
	
	/// configuration file and mapping
	string confFileName;
	string cleanMethod; // "logical" "confirmed" or "learned"
	map<string,confRecord> outlierParams;



	// Constructors 
	// Set Signal -- inheritted
	// Signal Id -- inheritted
	
	void init_defaults() {
		processor_type = REP_PROCESS_CONFIGURED_OUTLIER_CLEANER;
		params.trimming_sd_num = DEF_REP_TRIMMING_SD_NUM; params.removing_sd_num = DEF_REP_REMOVING_SD_NUM; params.nbrs_sd_num = 0;
		params.take_log = 0;
		params.doTrim = params.doRemove = true;
		signalId = -1;
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};
	// Init
	int init(map<string, string>& mapper);
	
	/// Learn cleaning model
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	

	// Apply cleaning model -inheritted
	

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void print();
};

/// a function that takes sorted vector of filtered values and estimates the +- 7 sd borders based on the center of distribution
/// predefined calibration constants are used for estimation of the borders. 
void learnDistributionBorders(float& borderHi, float& borderLo,vector<float> filteredValues);

/*
// A cleaner that is based on rules that describe relations of signal values to each other.\n
// This is a static cleaner ( no learning involved).\n
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
	
	map <int, vector<string>>rules2Signals = { ///< static map from rule to participating signals
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
	virtual int init(map<string, string>& mapper);



	
	
	
	// Learn cleaning model  : no learning for this cleaner. only apply
	
	///set signals
	void set_signal_ids(MedDictionarySections& dict) { myDict = dict; } // keep the dict. We will set ids later.

	/// Apply cleaning model 
	int apply(PidDynamicRec& rec, vector<int>& time_points);

	// Serialization  -static not needed
	//print 
private:
	bool  applyRule(int rule, vector <UniversalSigVec> ruleUsvs,vector <int >sPointer); // apply the rule and return true if data is consistent with the rule
	///<ruleUsvs hold the signals in the order they appear in the rule in the rules2Signals above
};

#define DEF_REP_NBRS_NBRS_SD_NUM 5
#define DEF_REP_NBRS_TRIM_SD_NUM 7
#define DEF_REP_NBRS_REMOVING_SD_NUM 14

/**
* A cleaner that looks at the neighbourhood of a certain signal value
*/
class RepNbrsOutlierCleaner : public RepProcessor, public MedValueCleaner {
public:

	/// Signal to clean
	string signalName;
	int signalId;
	int time_channel = 0;
	int val_channel = 0;
	int nbr_time_unit = MedTime::Days;
	int nbr_time_width = 7;


	// Constructors
	RepNbrsOutlierCleaner() : RepProcessor() { init_defaults(); }
	RepNbrsOutlierCleaner(const string& _signalName) : RepProcessor() { init_defaults(); signalId = -1; signalName = _signalName; init_lists(); }
	RepNbrsOutlierCleaner(const string& _signalName, string init_string) : RepProcessor() { init_defaults(); signalId = -1; signalName = _signalName; init_lists(); }
	RepNbrsOutlierCleaner(const string& _signalName, ValueCleanerParams *_params) : RepProcessor() { signalId = -1; signalName = _signalName; init_lists(); MedValueCleaner::init(_params); }

	void init_defaults() {
		processor_type = REP_PROCESS_NBRS_OUTLIER_CLEANER;
		params.trimming_sd_num = DEF_REP_NBRS_TRIM_SD_NUM; params.removing_sd_num = DEF_REP_NBRS_REMOVING_SD_NUM; params.nbrs_sd_num = DEF_REP_NBRS_NBRS_SD_NUM;
		params.take_log = 0;
		params.doTrim = params.doRemove = true;
		signalId = -1;
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};

	/// Set Signal
	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName; init_lists(); }

	/// Signal Id
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }

	// Init
	int init(void *processor_params) { return MedValueCleaner::init(processor_params); };
	int init(map<string, string>& mapper); // { init_defaults();  return MedValueCleaner::init(mapper); };
	void init_lists();

	/// Learn cleaning model
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	int iterativeLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	int quantileLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);

	/// Apply cleaning model
	int apply(PidDynamicRec& rec, vector<int>& time_points);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void print();
};
//.......................................................................................
//.......................................................................................
// Utility Functions
//.......................................................................................
//.......................................................................................

/// Get values of a signal from a set of ids
int get_values(MedRepository& rep, vector<int>& ids, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values, vector<RepProcessor *>& prev_cleaners);
int get_values(MedRepository& rep, vector<int>& ids, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values) ;

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(RepMultiProcessor)
MEDSERIALIZE_SUPPORT(RepBasicOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepNbrsOutlierCleaner)

#endif