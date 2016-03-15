// RepProcess : Apply processing on repository before generatring features
// E.g. - Cleaning.

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
//.......................................................................................
// A virtual class of processing repository data
//.......................................................................................
//.......................................................................................

typedef enum {
	REP_PROCESS_MULTI,
	REP_PROCESS_BASIC_OUTLIER_CLEANER,
	REP_PROCESS_NBRS_OUTLIER_CLEANER,
	REP_PROCESS_CONFIGURED_OUTLIER_CLEANER,
	REP_PROCESS_RULEBASED_OUTLIER_CLEANER,
	REP_PROCESS_LAST
} RepProcessorTypes;

class RepProcessor : public SerializableObject {
public:

	RepProcessorTypes processor_type;

	// Threading
	int learn_nthreads, apply_nthreads;

	// Required Signals
	vector<string> req_signals;
	vector<int> req_signal_ids;

	// Affected Signals
	unordered_set<string> aff_signals;
	unordered_set<int> aff_signal_ids;

	// Constructor/Destructor
	RepProcessor() { learn_nthreads = DEFAULT_REP_CLNR_NTHREADS; apply_nthreads = DEFAULT_REP_CLNR_NTHREADS; };
	~RepProcessor() {};

	// Virtual set-id function
	virtual void set_signal(const string& _signalName) { return; };

	// Required Signals functions
	void get_required_signal_ids(unordered_set<int>& signalIds, MedDictionarySections& dict);
	virtual void get_required_signal_names(unordered_set<string>& signalNames);
	virtual void set_required_signal_ids(MedDictionarySections& dict);

	// Affected Signals functions;
	virtual void set_affected_signal_ids(MedDictionarySections& dict);
	bool is_signal_affected(int signalId) {return (aff_signal_ids.find(signalId) != aff_signal_ids.end());}

	// Other signal ids
	virtual void set_signal_ids(MedDictionarySections& dict) { return; }

	// Init required tables
	virtual void init_tables(MedDictionarySections& dict) { return; }

	// Learn cleaning model
	virtual int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) { fprintf(stderr, "Not Learning Anything\n");  return 0; };

	int learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) { return Learn(rep, ids, prev_processors); }
	int learn(MedPidRepository& rep, vector<RepProcessor *>& prev_processors) { return Learn(rep, rep.pids, prev_processors); }

	int learn(MedPidRepository& rep, vector<int>& ids) { vector<RepProcessor *> temp;  return Learn(rep, ids, temp); }
	int learn(MedPidRepository& rep) { vector<RepProcessor *> temp; return learn(rep, temp); }

	// Apply cleaning model - ASSUME OUT_REC is preallocated !
	virtual int apply(PidDynamicRec& rec, vector<int>& time_points) {return 0; }
	virtual int apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignalIds);

	int apply(PidDynamicRec& rec, MedIdSamples& samples);
	int apply(PidDynamicRec& rec, MedIdSamples& samples, vector<int>& neededSignalIds);

	// Init
	static RepProcessor *make_processor(string name);
	static RepProcessor *make_processor(string type, string params);
	static RepProcessor *make_processor(RepProcessorTypes type);
	static RepProcessor *make_processor(RepProcessorTypes type, string params);

	// Create processor with just params (type is a must)
	static RepProcessor *create_processor(string &params);
	

	virtual int init(void *params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	// Serialization
	size_t get_processor_size();
	size_t processor_serialize(unsigned char *blob);

	// optional printing of cleaner
	virtual void print() { fprintf(stderr, "No implementation for print()\n"); }
};

// Utilities
RepProcessorTypes rep_processor_name_to_type(const string& procesor_name);

//.......................................................................................
//.......................................................................................
// A Processor which contains a vector of simpler processors
// Useful for applying same cleaners on a set of signals, for example
//.......................................................................................
//.......................................................................................

class RepMultiProcessor : public RepProcessor {
public:
	// Cleaners
	vector<RepProcessor *> processors;

	// Constructor/Destructor
	RepMultiProcessor() { processor_type = REP_PROCESS_MULTI; };
	~RepMultiProcessor() {};

	// Add processors
	void add_processors_set(RepProcessorTypes type, vector<string>& signals);
	void add_processors_set(RepProcessorTypes type, vector<string>& signals, string init_string);

	// Required Signals
	void set_required_signal_ids(MedDictionarySections& dict);

	void get_required_signal_names(unordered_set<string>& signalNames);

	// Affected Signals
	void set_affected_signal_ids(MedDictionarySections& dict);

	// Other signal ids
	void set_signal_ids(MedDictionarySections& dict);

	// Learn cleaning model
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors);

	// Apply cleaning model
	int apply(PidDynamicRec& rec, vector<int>& time_points);
	int apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignals);

	// serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void print() { for (auto& processor : processors) processor->print(); }
};

//.......................................................................................
//.......................................................................................
// A simple cleaner considering each value of a certain signal separatley
//.......................................................................................
//.......................................................................................

#define DEF_REP_TRIMMING_SD_NUM 7
#define DEF_REP_REMOVING_SD_NUM 14

class RepBasicOutlierCleaner : public RepProcessor, public MedValueCleaner {
public:

	// Signal to clean
	string signalName;
	int signalId;
	int time_channel = 0;
	int val_channel = 0;


	// Constructors
	RepBasicOutlierCleaner() : RepProcessor() { init_defaults(); }
	RepBasicOutlierCleaner(const string& _signalName) : RepProcessor() { init_defaults(); signalId = -1; signalName = _signalName; init_lists(); }
	RepBasicOutlierCleaner(const string& _signalName, string init_string) : RepProcessor() { init_defaults(); signalId = -1; signalName = _signalName; init_from_string(init_string); }
	RepBasicOutlierCleaner(const string& _signalName, ValueCleanerParams *_params) : RepProcessor() {signalId = -1; signalName = _signalName; init_lists() ; MedValueCleaner::init(_params);}

	void init_defaults() {
		processor_type = REP_PROCESS_BASIC_OUTLIER_CLEANER;
		params.trimming_sd_num = DEF_REP_TRIMMING_SD_NUM; params.removing_sd_num = DEF_REP_REMOVING_SD_NUM; params.nbrs_sd_num = 0;
		params.take_log = 0;
		params.doTrim = params.doRemove = true;
		signalId = -1;
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};

	// Set Signal
	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName;  }

	// Signal Id
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }

	// Init
	int init(void *processor_params) { return MedValueCleaner::init(processor_params); };
	virtual int init(map<string, string>& mapper); 
	void init_lists();

	// Learn cleaning model
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	int iterativeLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	int quantileLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);

	// Apply cleaning model
	int apply(PidDynamicRec& rec, vector<int>& time_points);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void print();
};

//.......................................................................................
//.......................................................................................
// A simple cleaner considering each value of a certain signal separatley, but this time use 
// configuration file that holds for each signal the logical values, statistically confirmed values 
// and distribution for relearning statistical values
//.......................................................................................
//.......................................................................................

typedef struct {
	double logicalLow, logicalHigh, confirmedLow, confirmedHigh;
    string distLow,distHigh; //"none" "norm" or "log" 
}confRecord;

class RepConfiguredOutlierCleaner : public  RepBasicOutlierCleaner {
public:

	// Signal to clean -- inheritted
	
	// configuration file and mapping
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
	
	// Learn cleaning model
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	

	// Apply cleaning model -inheritted
	

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void print();
};

void learnDistributionBorders(double& borderHi, double& borderLo,vector<float> filteredValues);
// a function that takes sorted vector of filtered values and estimates the +- 7 sd borders based on the center of distribution
// predefined calibration constants are used for estimation of the borders. 


//.......................................................................................
//.......................................................................................
// A cleaner that is based on rules that describe relations of signal values to each other.
// This is a static cleaner ( no learning involved).
//.......................................................................................
//.......................................................................................
/****************************************************************************
Rules:
Rule1: BMI = Weight / Height ^ 2 * 1e4
Rule2:MCH = Hemoglobin / RBC
Rule3:MCV = Hematocrit / RBC
Rule4:MCHC - M = MCH / MCV
Rule5:Eosinophils# + Monocytes# + Basophils# + Lymphocytes# + Neutrophils# <= WBC
Rule6:MPV = Platelets_Hematocrit / Platelets
Rule7:UrineAlbumin <= UrineTotalProtein
Rule8:UrineAlbumin_over_Creatinine = UrineAlbumin / UrineCreatinine
Rule9:LDL + HDL <= Cholesterol
Rule10:NonHDLCholesterol + HDL = Cholesterol
Rule11:HDL_over_nonHDL = HDL / NonHDLCholesterol
Rule12:HDL_over_Cholesterol = HDL / Cholesterol
Rule13:HDL_over_LDL = HDL / LDL
Rule14:HDL_over_LDL = 1 / LDL_over_HDL
Rule15:Cholesterol_over_HDL = Cholesterol / HDL
Rule16:Cholesterol_over_HDL = 1 / HDL / Cholesterol
Rule17:Cholesterol_over_HDL = 1 / HDL_over_Cholestrol
Rule18:LDL_over_HDL = LDL / HDL
Rule19:Albumin <= Protein_Total
Rule20:FreeT4 <= T4
Rule21:NRBC <= RBC
Rule22:CHADS2_VASC >= CHADS2_VASC

******************************************************************************/

class RepRuleBasedOutlierCleaner : public RepProcessor, public MedValueCleaner {
	// get multiple signals and clean them according to consistency  with other signals from same date
public:

	
	// Signals to clean
	vector <string> signalNames;
	vector <int> signalIds;
	int time_channel = 0;
	int val_channel = 0;
	MedDictionarySections myDict; // keeping it will enable us to get ids at apply stage
	
	
	map <int, vector<string>>rules2Signals = { // static map from rule to participating signals
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
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor) {
		return 0;
	};
	//set signals
	void set_signal_ids(MedDictionarySections& dict) { myDict = dict; } // keep the dict. We will set ids later.

	// Apply cleaning model 
	int apply(PidDynamicRec& rec, vector<int>& time_points);

	// Serialization  -static not needed
	//print 
private:
	bool  applyRule(int rule, vector <UniversalSigVec> ruleUsvs,vector <int >sPointer); // apply the rule and return true if data is consistent with the rule
	//ruleUsvs hold the signals in the order they appear in the rule in the rules2Signals above
};



//.......................................................................................
//.......................................................................................
// A cleaner that looks at the neighbourhood of a certain signal value
//.......................................................................................
//.......................................................................................

#define DEF_REP_NBRS_NBRS_SD_NUM 5
#define DEF_REP_NBRS_TRIM_SD_NUM 7
#define DEF_REP_NBRS_REMOVING_SD_NUM 14

class RepNbrsOutlierCleaner : public RepProcessor, public MedValueCleaner {
public:

	// Signal to clean
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

	// Set Signal
	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName; init_lists(); }

	// Signal Id
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }

	// Init
	int init(void *processor_params) { return MedValueCleaner::init(processor_params); };
	int init(map<string, string>& mapper); // { init_defaults();  return MedValueCleaner::init(mapper); };
	void init_lists();

	// Learn cleaning model
	int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	int iterativeLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);
	int quantileLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processor);

	// Apply cleaning model
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

// Get values of a signal from a set of ids
int get_values(MedRepository& rep, vector<int>& ids, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values, vector<RepProcessor *>& prev_cleaners);
int get_values(MedRepository& rep, vector<int>& ids, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values) ;

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(RepMultiProcessor)
MEDSERIALIZE_SUPPORT(RepBasicOutlierCleaner)
MEDSERIALIZE_SUPPORT(RepNbrsOutlierCleaner)

#endif