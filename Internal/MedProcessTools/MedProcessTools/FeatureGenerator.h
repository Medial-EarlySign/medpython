// FeatureProcess : Apply processing on features
// E.g. - Cleaning, Normalizing

#ifndef _FTR_GENERATOR_H_
#define _FTR_GENERATOR_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "Logger/Logger/Logger.h"
#include "MedProcessTools/MedProcessTools/RepProcess.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"

#define DEFAULT_FEAT_GNRTR_NTHREADS 8

// TBD : Add wrapper for management of features list (read/write to file, etc.)

//.......................................................................................
//.......................................................................................
// A virtual feature generation class
//.......................................................................................
//.......................................................................................
// Types of feature generators
typedef enum {
	FTR_GEN_BASIC,
	FTR_GEN_AGE,
	FTR_GEN_GENDER,
	FTR_GEN_BINNED_LM,
	FTR_GEN_LAST
} FeatureGeneratorTypes;

class FeatureGenerator : public SerializableObject {
public:

	// Type
	FeatureGeneratorTypes generator_type;

	// Feature name
	vector<string> names;

	// Threading
	int learn_nthreads, pred_nthreads;

	// Missing value
	float missing_val;

	// Naming
	virtual void set_names() {names.clear(); }

	// Initialize
	virtual void init(MedFeatures &features);

	// Constructor/Destructor
	FeatureGenerator() { learn_nthreads = DEFAULT_FEAT_GNRTR_NTHREADS; pred_nthreads = DEFAULT_FEAT_GNRTR_NTHREADS;  missing_val = MED_MAT_MISSING_VALUE; };
	~FeatureGenerator() {};

	// Required Signals
	vector<string> req_signals;
	vector<int> req_signal_ids; 
	void get_required_signal_ids(unordered_set<int>& signalIds, MedDictionarySections& dict);
	void get_required_signal_ids(MedDictionarySections& dict);

	// Learn a generator
	virtual int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) {return 0;}
	int learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { set_names(); return _learn(rep, ids, processors); }
	int learn(MedPidRepository& rep, vector<RepProcessor *> processors) { set_names(); return _learn(rep, rep.pids, processors); }
	int learn(MedPidRepository& rep, vector<int>& ids) { set_names(); return _learn(rep, ids, vector<RepProcessor *>()); }
	int learn(MedPidRepository& rep) { set_names(); return _learn(rep, rep.pids, vector<RepProcessor *>()); }

	// generate feature data from repository
	// We assume the corresponding MedSamples have been inserted to MedFeatures : either at the end or at position index
	virtual int Generate(PidDynamicRec& in_rep, MedFeatures& features, int index, int num) { return 0; }
	int generate(PidDynamicRec& in_rep, MedFeatures& features, int index, int num) { return Generate(in_rep, features, index, num); }
	int generate(PidDynamicRec& in_rep, MedFeatures& features);
	int generate(MedPidRepository& rep, int id, MedFeatures& features) ;
	int generate(MedPidRepository& rep, int id,  MedFeatures& features, int index, int num) ;

	// generate feature data from other features
	virtual int Generate(MedFeatures& features) { return 0; }
	int generate(MedFeatures& features) { return Generate(features); }

	// Init
	// create a generator
	static FeatureGenerator *make_generator(string name);
	static FeatureGenerator *make_generator(string name, string params);
	static FeatureGenerator *make_generator(FeatureGeneratorTypes type);
	static FeatureGenerator *make_generator(FeatureGeneratorTypes type, string params);

	virtual int init(void *generator_params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	// Number of features generated
	virtual int nfeatures() { return 1; }

	// Serialization
	size_t get_generator_size();
	size_t generator_serialize(unsigned char *blob);
};

FeatureGeneratorTypes ftr_generator_name_to_type(const string& generator_name);

//.......................................................................................
//.......................................................................................
// Single signal features that do not require learning(e.g. last hemoglobin)
//.......................................................................................
//.......................................................................................

typedef enum {
	FTR_LAST_VALUE = 0,
	FTR_FIRST_VALUE = 1,
	FTR_LAST2_VALUE = 2,
	FTR_AVG_VALUE = 3,
	FTR_MAX_VALUE = 4,
	FTR_MIN_VALUE = 5,
	FTR_STD_VALUE = 6,
	FTR_LAST_DELTA_VALUE = 7,
	FTR_LAST_DAYS = 8,
	FTR_LAST2_DAYS = 9,

	FTR_LAST
} BasicFeatureTypes;

class BasicFeatGenerator : public FeatureGenerator {
public:
	// Feature Descrption
	string signalName;
	BasicFeatureTypes type;
	int win_from, win_to;		// time window for feature: date-win_to <= t < date-win_from

	int signalId;

	// Naming 
	void set_names();

	// Constructor/Destructor
	BasicFeatGenerator() : FeatureGenerator() { init_defaults(); };
	~BasicFeatGenerator() {};
	void set(string& _signalName, BasicFeatureTypes _type) { set(_signalName, _type, 0, 360000); req_signals.push_back(signalName);}
	void set(string& _signalName, BasicFeatureTypes _type, int _time_win_from, int _time_win_to) { 
		signalName = _signalName;type = _type; win_from = _time_win_from; win_to = _time_win_to;
		set_names(); req_signals.push_back(signalName);
	}

	// Init
	int init(map<string, string>& mapper);
	void init_defaults() { generator_type = FTR_GEN_BASIC; signalId = -1; string _signalName = ""; set(_signalName, FTR_LAST, 0, 360000); ; };

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);
	float get_value(PidDynamicRec& rec, int index, int date);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

//.......................................................................................
//.......................................................................................
// Age
//.......................................................................................
//.......................................................................................
class AgeGenerator : public FeatureGenerator {
public:

	// Age Id
	int byearId;

	// Constructor/Destructor
	AgeGenerator() : FeatureGenerator() { generator_type = FTR_GEN_AGE; names.push_back("Age"); byearId = -1; req_signals.push_back("BYEAR");}
	AgeGenerator(int _byearId) : FeatureGenerator() { generator_type = FTR_GEN_AGE; names.push_back("Age"); byearId = _byearId; req_signals.push_back("BYEAR"); }
	~AgeGenerator() {};

	// Name
	void set_names() { if (names.empty()) names.push_back("Age"); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Serialization
	size_t get_size() { return 0; };
	size_t serialize(unsigned char *blob) { return 0; };
	size_t deserialize(unsigned char *blob) { set_names();  return 0; };
};

//.......................................................................................
//.......................................................................................
// Gender
//.......................................................................................
//.......................................................................................
class GenderGenerator : public FeatureGenerator {
public:

	// Gender Id
	int genderId;

	// Constructor/Destructor
	GenderGenerator() : FeatureGenerator() { generator_type = FTR_GEN_GENDER; names.push_back("Gender"); genderId = -1; req_signals.push_back("GENDER");}
	GenderGenerator(int _genderId) : FeatureGenerator() {generator_type = FTR_GEN_GENDER; names.push_back("Gender"); genderId = _genderId; req_signals.push_back("GENDER");}

	~GenderGenerator() {};

	// Name
	void set_names() { if (names.empty()) names.push_back("Gender"); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Serialization
	size_t get_size() { return 0; };
	size_t serialize(unsigned char *blob) { return 0; };
	size_t deserialize(unsigned char *blob) { set_names();  return 0; };
};

//.......................................................................................
//.......................................................................................
// BinnedLinearModels : Apply a set of liner models to generate features
//.......................................................................................
//.......................................................................................

struct BinnedLmEstimatesParams {
	vector<int> bin_bounds;
	int min_period ;
	int max_period ;
	float rfactor;

	vector<int> estimation_points;

};

class BinnedLmEstimates : public FeatureGenerator {
public:
	// Feature Descrption
	string signalName;
	int signalId, byearId, genderId;

	BinnedLmEstimatesParams params;
	vector<MedLM> models;
	vector<float> xmeans, xsdvs, ymeans, ysdvs;
	vector<float> means[2];
	
	// Naming 
	void set_names();

	// Constructor/Destructor
	BinnedLmEstimates() : FeatureGenerator() { signalName = ""; init_defaults(); };
	BinnedLmEstimates(string _signalName) : FeatureGenerator() { signalName = _signalName; init_defaults(); req_signals.push_back(signalName); set_names(); };
	BinnedLmEstimates(string _signalName, string init_string) : FeatureGenerator() {signalName = _signalName;init_defaults(); req_signals.push_back(signalName); init_from_string(init_string); set_names();};

	~BinnedLmEstimates() {};

	void set(string& _signalName);
	void set(string& _signalName, BinnedLmEstimatesParams* _params);

	// Init
	void init_defaults();
	int init(map<string, string>& mapper);

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors);

	// generate new feature(s)
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Number of features generated
	virtual int nfeatures() { return (int) params.estimation_points.size(); }

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

// Utilities
float get_last_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_first_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_last2_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_avg_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_max_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_min_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_std_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_last_delta_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_last_days_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);
float get_last2_days_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val);


#endif