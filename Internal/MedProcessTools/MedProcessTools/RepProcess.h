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
// RepProcessor is the parent class for processing a MedRepository or PidDynamicRec
// Basic functionalities:
//		learn : learn the processoring parameters from a given list of ids and a rpository 
//		apply : process a dynamic PidDynamicRec
//.......................................................................................
//.......................................................................................
// Define types of repository processors
typedef enum {
	REP_PROCESS_MULTI,
	REP_PROCESS_BASIC_OUTLIER_CLEANER,
	REP_PROCESS_NBRS_OUTLIER_CLEANER,
	REP_PROCESS_CONFIGURED_OUTLIER_CLEANER,
	REP_PROCESS_LAST
} RepProcessorTypes;

class RepProcessor : public SerializableObject {
public:

	// Type
	RepProcessorTypes processor_type;

	// Threading (learn/apply) - 
	int learn_nthreads, apply_nthreads;

	// Required Signals (and their ids) for processing
	vector<string> req_signals;
	vector<int> req_signal_ids;

	// Signals (and their ids) affected by processing 
	unordered_set<string> aff_signals;
	unordered_set<int> aff_signal_ids;

	// Constructor/Destructor
	RepProcessor() { learn_nthreads = DEFAULT_REP_CLNR_NTHREADS; apply_nthreads = DEFAULT_REP_CLNR_NTHREADS; };
	~RepProcessor() {};

	// Init from name or type. optinally with a parameters string
	static RepProcessor *make_processor(string name);
	static RepProcessor *make_processor(string type, string params);
	static RepProcessor *make_processor(RepProcessorTypes type);
	static RepProcessor *make_processor(RepProcessorTypes type, string params);

	// Create processor from params string (type must be given within string)
	static RepProcessor *create_processor(string &params);

	// initialize : from object/string/defaults.
	// Should be implemented for inheriting classes that have parameters
	virtual int init(void *params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	// Set signalName
	// Should be implemented for inheriting classes that have signalName
	virtual void set_signal(const string& _signalName) { return; };

	// Set signalId
	// Should be implemented for inheriting classes that have signalId
	virtual void set_signal_ids(MedDictionarySections& dict) { return; }

	// Required Signals functions : get all signals that are required by the processor
	// Append required signal names to vector : parent function just uses req_signals
	virtual void get_required_signal_names(unordered_set<string>& signalNames);
	// Fill req_signal_ids : parent function just fills from req_signals
	virtual void set_required_signal_ids(MedDictionarySections& dict);
	// Append required signal ids to vector. call set_required_signal_ids if req_signal_ids empty
	void get_required_signal_ids(unordered_set<int>& signalIds, MedDictionarySections& dict);
	
	// Affected Signals functions;
	// Fill aff_signal_ids : parent function just fills from aff_signals
	virtual void set_affected_signal_ids(MedDictionarySections& dict);
	// Check if a signal is affected by processor
	bool is_signal_affected(int signalId) {return (aff_signal_ids.find(signalId) != aff_signal_ids.end());}

	// Init required tables
	// Should be implemented for inheriting classes that have such tables
	virtual void init_tables(MedDictionarySections& dict) { return; }

	// Learn processing model on a subset of ids. Apply set of preceesing processors on DynamicPidRec before learning
	// Should be implemented for inheriting classes that require learning
	virtual int Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) { fprintf(stderr, "Not Learning Anything\n");  return 0; };

	// Envelope learning functions
	// Learn on subset of ids
	int learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) { return Learn(rep, ids, prev_processors); }
	// Learn on all ids in repository
	int learn(MedPidRepository& rep, vector<RepProcessor *>& prev_processors) { return Learn(rep, rep.pids, prev_processors); }
	// Learn on subset of ids without preceesing processors
	int learn(MedPidRepository& rep, vector<int>& ids) { vector<RepProcessor *> temp;  return Learn(rep, ids, temp); }
	// Learn on all ids in repository without preceesing processors
	int learn(MedPidRepository& rep) { vector<RepProcessor *> temp; return learn(rep, temp); }

	// Apply processing on a single PidDynamicRec at a set of time-points
	// Should be implemented for all inheriting classes
	virtual int apply(PidDynamicRec& rec, vector<int>& time_points) {return 0; }
	// Apply processing on a single PidDynamicRec at a set of time-points only if required : if any of the signals in neededSignalIds is actually affected by processor
	int apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignalIds);
	// Apply processing on a single PidDynamicRec at a set of time-points given by Samples
	int apply(PidDynamicRec& rec, MedIdSamples& samples);
	// Apply processing on a single PidDynamicRec at a set of time-points given by Samples only if required
	int apply(PidDynamicRec& rec, MedIdSamples& samples, vector<int>& neededSignalIds);	

	// Serialization (including type)
	size_t get_processor_size();
	size_t processor_serialize(unsigned char *blob);

	// optional printing of processor
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
	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName; init_lists(); }

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
	RepNbrsOutlierCleaner(const string& _signalName, string init_string) : RepProcessor() {init_defaults(); signalId = -1; signalName = _signalName; init_lists();}
	RepNbrsOutlierCleaner(const string& _signalName, ValueCleanerParams *_params) : RepProcessor() {signalId = -1; signalName = _signalName; init_lists(); MedValueCleaner::init(_params);}

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
	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName; init_lists();}

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