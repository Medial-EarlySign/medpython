// FeatureProcess : Apply processing on features
// E.g. - Cleaning, Normalizing

#ifndef _FTR_GENERATOR_H_
#define _FTR_GENERATOR_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "Logger/Logger/Logger.h"
#include "MedProcessTools/MedProcessTools/RepProcess.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include <MedTime/MedTime/MedTime.h>
#include "InfraMed/InfraMed/MedRepositoryType.h"

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

	// SerialId

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
	FeatureGenerator() { learn_nthreads = DEFAULT_FEAT_GNRTR_NTHREADS; pred_nthreads = DEFAULT_FEAT_GNRTR_NTHREADS;  missing_val = MED_MAT_MISSING_VALUE; serial_id = ++global_serial_id_cnt; };
	~FeatureGenerator() {};

	// Required Signals
	vector<string> req_signals;
	vector<int> req_signal_ids; 
	void get_required_signal_ids(unordered_set<int>& signalIds, MedDictionarySections& dict);
	void get_required_signal_names(unordered_set<string>& signalNames);
	virtual void set_required_signal_ids(MedDictionarySections& dict);

	// Signal Ids
	virtual void set_signal_ids(MedDictionarySections& dict) { return; }

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

	virtual int init(void *generator_params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	// Number of features generated
	virtual int nfeatures() { return 1; }

	// Serialization
	size_t get_generator_size();
	size_t generator_serialize(unsigned char *blob);

	virtual void print() { fprintf(stderr, "Print Not Implemented for feature\n"); }

	static int global_serial_id_cnt;
	int serial_id;		// serial id of feature to 

};


/*
see http://stackoverflow.com/a/582456/574187
this factory creates feature generators by their class names
*/
template<typename T> FeatureGenerator * createFeatureGenerator() { return new T; }
struct FeatureGeneratorFactory {
	// create a generator
	static FeatureGenerator *make_generator(string name);
	static FeatureGenerator *make_generator(string name, string params);
	static FeatureGenerator *make_generator(FeatureGeneratorTypes type);
	static FeatureGenerator *make_generator(FeatureGeneratorTypes type, string params);
	static FeatureGenerator *create_generator(string &params); // must include fg_type	
	typedef std::map<std::string, FeatureGenerator*(*)()> map_type;
protected:
	static map_type * getMap() {
		// never delete'ed. (exist until program termination)
		// because we can't guarantee correct destruction order 
		if (!my_map) { my_map = new map_type; }
		return my_map;
	}
private:
	static map_type * my_map;
};

template<typename T>
struct DerivedRegister : FeatureGeneratorFactory {
	DerivedRegister(std::string const& s) {
		getMap()->insert(std::make_pair(s, &createFeatureGenerator<T>));
	}
};

FeatureGeneratorTypes ftr_generator_name_to_type(const string& generator_name);

#define DEC_FEATURE_GENERATOR(NAME) \
    static DerivedRegister<NAME> reg

#define DEF_FEATURE_GENERATOR(NAME) \
    DerivedRegister<NAME> NAME::reg(#NAME)

//..............................................................................................
// FeatureSingleChannel -
// This class is a mediator between FeatureGenerator and classes that generate
// Features on a single variable (not including age and gender) and in it in a single channel.
//..............................................................................................



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
	FTR_SLOPE_VALUE = 10,
	FTR_WIN_DELTA_VALUE = 11,
	FTR_CATEGORY_SET = 12,
	FTR_CATEGORY_SET_COUNT = 13,
	FTR_CATEGORY_SET_SUM = 14,
	FTR_NSAMPLES = 15,

	FTR_LAST
} BasicFeatureTypes;

class BasicFeatGenerator : public FeatureGenerator {
public:
	// Feature Descrption
	string signalName;
	int signalId;

	// parameters (should be serialized)
	BasicFeatureTypes type = FTR_LAST;
	int win_from = 0, win_to = 360000;			// time window for feature: date-win_to <= t < date-win_from
	int d_win_from = 360, d_win_to = 360000;	// delta time window for the FTR_WIN_DELTA_VALUE feature
	int time_unit_win = MedTime::Undefined;			// the time unit in which the windows are given. Default: Undefined
	int time_channel = 0;						// n >= 0 : use time channel n , default: 0.
	int val_channel = 0;						// n >= 0 : use val channel n , default : 0.
	int sum_channel = 1;						// for FTR_CETEGORY_SET_SUM
	vector<string> sets;						// for FTR_CATEGORY_SET_* , the list of sets 
	int time_unit_sig = MedTime::Undefined;		// the time init in which the signal is given. (set correctly from Repository in learn and Generate)

	// helpers
	vector<char> lut;							// to be used when generating FTR_CATEGORY_SET_*

	// Naming 
	void set_names();

	// Constructor/Destructor
	BasicFeatGenerator() : FeatureGenerator() { init_defaults(); };
	~BasicFeatGenerator() {};
	void set(string& _signalName, BasicFeatureTypes _type) { set(_signalName, _type, 0, 360000); req_signals.assign(1,signalName);}
	void set(string& _signalName, BasicFeatureTypes _type, int _time_win_from, int _time_win_to) { 
		signalName = _signalName;type = _type; win_from = _time_win_from; win_to = _time_win_to;
		set_names(); req_signals.assign(1, signalName);
	}

	BasicFeatureTypes name_to_type(const string &name);

	// Init
	int init(map<string, string>& mapper);
	void init_defaults();

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(signalName)].time_unit; return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);
	float get_value(PidDynamicRec& rec, int index, int date);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }

	// actual generators
	float uget_last(UniversalSigVec &usv, int time_point, int _win_from, int _win_to); // Added the win as needed to be called on different ones in uget_win_delta
	float uget_first(UniversalSigVec &usv, int time_point);
	float uget_last2(UniversalSigVec &usv, int time_point);
	float uget_avg(UniversalSigVec &usv, int time_point);
	float uget_max(UniversalSigVec &usv, int time_point);
	float uget_min(UniversalSigVec &usv, int time_point);
	float uget_std(UniversalSigVec &usv, int time_point);
	float uget_last_delta(UniversalSigVec &usv, int time_point);
	float uget_last_time(UniversalSigVec &usv, int time_point);
	float uget_last2_time(UniversalSigVec &usv, int time_point);
	float uget_slope(UniversalSigVec &usv, int time_point);
	float uget_win_delta(UniversalSigVec &usv, int time_point);
	float uget_category_set(PidDynamicRec &rec, UniversalSigVec &usv, int time_point);
	float uget_category_set_count(PidDynamicRec &rec, UniversalSigVec &usv, int time_point);
	float uget_category_set_sum(PidDynamicRec &rec, UniversalSigVec &usv, int time_point);
	float uget_nsamples(UniversalSigVec &usv, int time, int _win_from, int _win_to);
	// helpers

	// gets a [-_win_to, -_win_from] window in win time unit, and returns [_min_time, _max_time] window in signal time units relative to _win_time
	void get_window_in_sig_time(int _win_from, int _win_to, int _time_unit_win, int _time_unit_sig, int _win_time, int &_min_time, int &_max_time);


	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
	DEC_FEATURE_GENERATOR(BasicFeatGenerator); 
};

//.......................................................................................
//.......................................................................................
// Age
//.......................................................................................
//.......................................................................................
class AgeGenerator : public FeatureGenerator {
public:

	// Is Age Directly given ?
	bool directlyGiven;

	// Signak Id
	int signalId;

	// Constructor/Destructor
	AgeGenerator() : FeatureGenerator() { generator_type = FTR_GEN_AGE; names.push_back("Age"); signalId = -1; directlyGiven = med_rep_type.ageDirectlyGiven; }
	AgeGenerator(int _signalId) : FeatureGenerator() { generator_type = FTR_GEN_AGE; names.push_back("Age"); signalId = _signalId; directlyGiven = med_rep_type.ageDirectlyGiven;}
	~AgeGenerator() {};

	// Name
	void set_names() { if (names.empty()) names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".Age"); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { if (directlyGiven) signalId = dict.id("Age");  else signalId = dict.id("BYEAR"); }
	void set_required_signal_ids(MedDictionarySections& dict) {if (directlyGiven) req_signal_ids.assign(1, dict.id("Age"));  else req_signal_ids.assign(1, dict.id("BYEAR")); }

	// Serialization
	size_t get_size() { return MedSerialize::get_size(generator_type, names); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, names); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, names); }
	DEC_FEATURE_GENERATOR(AgeGenerator);
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
	GenderGenerator() : FeatureGenerator() { generator_type = FTR_GEN_GENDER; names.push_back("Gender"); genderId = -1; req_signals.assign(1, med_rep_type.genderSignalName);}
	GenderGenerator(int _genderId) : FeatureGenerator() {generator_type = FTR_GEN_GENDER; names.push_back("Gender"); genderId = _genderId; req_signals.assign(1, med_rep_type.genderSignalName);}

	~GenderGenerator() {};

	// Name
	void set_names() { if (names.empty()) names.push_back("FTR_" + int_to_string_digits(serial_id,6) + ".Gender"); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { genderId = dict.id(med_rep_type.genderSignalName); }
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.assign(1, dict.id(med_rep_type.genderSignalName)); }

	// Serialization
	size_t get_size() { return MedSerialize::get_size(generator_type, names); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, names); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, names); }
	DEC_FEATURE_GENERATOR(GenderGenerator);
};

//.......................................................................................
//.......................................................................................
// BinnedLinearModels : Apply a set of liner models to generate features
//.......................................................................................
//.......................................................................................

struct BinnedLmEstimatesParams {
	vector<int> bin_bounds ;
	int min_period ;
	int max_period ;
	float rfactor;

	vector<int> estimation_points;

};

class BinnedLmEstimates : public FeatureGenerator {
public:
	// Feature Descrption
	string signalName;
	int signalId, byearId, genderId, ageId;

	// Is age directly given (Hospital data) or through birth-year (Primary-care data)
	bool ageDirectlyGiven ;

	BinnedLmEstimatesParams params;
	vector<MedLM> models;
	vector<float> xmeans, xsdvs, ymeans, ysdvs;
	vector<float> means[2];

	int time_unit_periods = MedTime::Undefined;		// the time unit in which the periods are given. Default: Undefined
	int time_unit_sig = MedTime::Undefined;			// the time init in which the signal is given. Default: Undefined
	int time_channel = 0;						// n >= 0 : use time channel n , default: 0.
	int val_channel = 0;						// n >= 0 : use val channel n , default : 0.

	// Naming 
	void set_names();

	// Constructor/Destructor
	BinnedLmEstimates() : FeatureGenerator() { signalName = ""; init_defaults(); };
	BinnedLmEstimates(string _signalName) : FeatureGenerator() { signalName = _signalName; init_defaults(); req_signals.push_back(signalName); names.clear();  set_names(); };
	BinnedLmEstimates(string _signalName, string init_string) : FeatureGenerator() { signalName = _signalName; init_defaults(); req_signals.push_back(signalName); init_from_string(init_string); };

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

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict);

	// Number of features generated
	virtual int nfeatures() { return (int) params.estimation_points.size(); }

	// Age Related functions
	void prepare_for_age(PidDynamicRec& rec, UniversalSigVec& ageUsv, int &age, int &byear);
	void prepare_for_age(MedPidRepository& rep, int id, UniversalSigVec& ageUsv, int &age, int &byear);
	inline void get_age(int time, int time_unit_from, int& age, int byear);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
	DEC_FEATURE_GENERATOR(BinnedLmEstimates);
	// print 
	void print();
};

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(BasicFeatGenerator)
MEDSERIALIZE_SUPPORT(AgeGenerator)
MEDSERIALIZE_SUPPORT(GenderGenerator)
MEDSERIALIZE_SUPPORT(BinnedLmEstimates)

#endif