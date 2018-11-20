#ifndef _FTR_GENERATOR_H_
#define _FTR_GENERATOR_H_

#include <InfraMed/InfraMed/InfraMed.h>
#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/RepProcess.h>
#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <MedProcessTools/MedProcessTools/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedModelExceptions.h>
#include <MedTime/MedTime/MedTime.h>
#include <MedProcessTools/MedProcessTools/MedGlobals.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>

#define DEFAULT_FEAT_GNRTR_NTHREADS 8

// For ModelFeatureGenerator
class MedModel;

// TBD : Add wrapper for management of features list (read/write to file, etc.)

/** @enum
* Types of feature generators
*/
typedef enum {
	FTR_GEN_NOT_SET,
	FTR_GEN_BASIC, ///< "basic" - creates basic statistic on time windows - BasicFeatGenerator 
	FTR_GEN_AGE, ///< "age" - creating age feature - AgeGenerator 
	FTR_GEN_SINGLETON, ///< "singleton" - take the value of a time-less signale - SingletonGenerator
	FTR_GEN_GENDER, ///< "gender" - creating gender feature - GenderGenerator (special case of signleton)
	FTR_GEN_BINNED_LM, ///< "binnedLm" or "binnedLM" - creating linear model for esitmating feature in time points - BinnedLmEstimates
	FTR_GEN_SMOKING, ///< "smoking" - creating smoking feature - SmokingGenerator
	FTR_GEN_KP_SMOKING, ///< "kp_smoking" - creating smoking feature - KpSmokingGenerator
	FTR_GEN_RANGE, ///< "range" - creating RangeFeatGenerator
	FTR_GEN_DRG_INTAKE, ///< "drugIntake" - creating drugs feature coverage of prescription time - DrugIntakeGenerator
	FTR_GEN_ALCOHOL, ///< "alcohol" - creating alcohol feature - AlcoholGenerator
	FTR_GEN_MODEL, ///< "model" - creating ModelFeatGenerator
	FTR_GEN_TIME, ///< "time" - creating sample-time features (e.g. differentiate between times of day, season of year, days of the week, etc.)
	FTR_GEN_LAST
} FeatureGeneratorTypes;

/** @file
* FeatureGenerator : creating features from raw signals
*/
class FeatureGenerator : public SerializableObject {
public:

	/// Type
	FeatureGeneratorTypes generator_type;

	/// Feature name
	vector<string> names;

	// Threading
	int learn_nthreads, pred_nthreads;

	/// Missing value
	float missing_val;

	/// Tags - for defining labels or groups. may be used later for filtering for example
	vector<string> tags;

	/// Feature/Weights generator
	int iGenerateWeights = 0;

	// Naming
	virtual void set_names() { names.clear(); }

	// Helper - pointers to data vectors in MedFeatures (to save time in generatin)
	vector <float *> p_data;

	// Prepare for feature generation
	virtual void prepare(MedFeatures &features, MedPidRepository& rep, MedSamples& samples);
	virtual void get_p_data(MedFeatures& features);

	// Constructor/Destructor
	FeatureGenerator() { learn_nthreads = DEFAULT_FEAT_GNRTR_NTHREADS; pred_nthreads = DEFAULT_FEAT_GNRTR_NTHREADS;  missing_val = MED_MAT_MISSING_VALUE; serial_id = ++MedFeatures::global_serial_id_cnt; };
	virtual ~FeatureGenerator() { clear(); };
	virtual void clear() { };

	// Required Signals
	vector<string> req_signals;
	vector<int> req_signal_ids;

	void get_required_signal_names(unordered_set<string>& signalNames);
	virtual void set_required_signal_ids(MedDictionarySections& dict);
	void get_required_signal_ids(unordered_set<int>& signalIds);

	// Signal Ids
	virtual void set_signal_ids(MedDictionarySections& dict) { return; }

	// Init required tables
	virtual void init_tables(MedDictionarySections& dict) { return; }

	// Learn a generator
	virtual int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }
	int learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { set_names(); return _learn(rep, ids, processors); }
	int learn(MedPidRepository& rep, vector<RepProcessor *> processors) { set_names(); return _learn(rep, rep.pids, processors); }
	int learn(MedPidRepository& rep, vector<int>& ids) { set_names(); return _learn(rep, ids, vector<RepProcessor *>()); }
	int learn(MedPidRepository& rep) { set_names(); return _learn(rep, rep.pids, vector<RepProcessor *>()); }

	// generate feature data from repository
	// We assume the corresponding MedSamples have been inserted to MedFeatures : either at the end or at position index
	virtual int _generate(PidDynamicRec& in_rep, MedFeatures& features, int index, int num) { return 0; }
	int generate(PidDynamicRec& in_rep, MedFeatures& features, int index, int num) { return _generate(in_rep, features, index, num); }
	int generate(PidDynamicRec& in_rep, MedFeatures& features);
	int generate(MedPidRepository& rep, int id, MedFeatures& features);
	int generate(MedPidRepository& rep, int id, MedFeatures& features, int index, int num);

	// generate feature data from other features
	virtual int _generate(MedFeatures& features) { return 0; }
	int generate(MedFeatures& features) { return _generate(features); }

	// Init
	// create a generator
	static FeatureGenerator *make_generator(string name);
	static FeatureGenerator *make_generator(string name, string params);
	static FeatureGenerator *make_generator(FeatureGeneratorTypes type);
	static FeatureGenerator *make_generator(FeatureGeneratorTypes type, string params);

	static FeatureGenerator *create_generator(string &params); // must include fg_type

	virtual int init(void *generator_params) { return 0; };
	virtual int init(map<string, string>& mapper);
	virtual void init_defaults() {};

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *generator; }

	// Number of features generated
	virtual int nfeatures() { return (int)names.size(); }

	// Filter generated features according to a set. return number of valid features
	virtual int filter_features(unordered_set<string>& validFeatures);

	// Serialization
	size_t get_generator_size();
	size_t generator_serialize(unsigned char *blob);

	virtual void print() { fprintf(stderr, "Print Not Implemented for feature\n"); }


	// debug print for a feature generator. fg_flag can 
	void dprint(const string &pref, int fg_flag);


	int serial_id;		// serial id of feature
};

FeatureGeneratorTypes ftr_generator_name_to_type(const string& generator_name);

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

/** @enum
* BasicFeatGenerator types for calculating stats
*/
typedef enum {
	FTR_LAST_VALUE = 0, ///<"last" - Last Value in Window
	FTR_FIRST_VALUE = 1, ///<"first" - First Value in Window
	FTR_LAST2_VALUE = 2, ///<"last2" - One before last value in Window
	FTR_AVG_VALUE = 3, ///<"avg" - Mean value in Window
	FTR_MAX_VALUE = 4, ///<"max" - Max value in Window
	FTR_MIN_VALUE = 5, ///<"min" - Min value in Window
	FTR_STD_VALUE = 6, ///<"std" - Standart Dev. value in Window
	FTR_LAST_DELTA_VALUE = 7, ///<"last_delta" - Last delta. last-previous_last value
	FTR_LAST_DAYS = 8, ///<"last_time" - time diffrence from prediction time to last time has signal
	FTR_LAST2_DAYS = 9,///<"last_time2" - time diffrence from prediction time to one previous last time has signal
	FTR_SLOPE_VALUE = 10, ///<"slope" - calculating the slope over the points in the window
	FTR_WIN_DELTA_VALUE = 11, ///<"win_delta" - diffrence in value in two time windows (only if both exists, otherwise missing_value). value in [win_from,win_to] minus value in [d_win_from, d_win_to]
	FTR_CATEGORY_SET = 12, ///<"category_set" - boolean 0/1 if the signal has the value in the given lut (which initialized by the "sets" that can be specific single definition or name of set definition. the lookup is hierarchical)
	FTR_CATEGORY_SET_COUNT = 13,///<"category_set_count" - counts the number of appearnces of sets in the time window
	FTR_CATEGORY_SET_SUM = 14, ///<"category_set_sum" - sums the values of appearnces of sets in the time window
	FTR_NSAMPLES = 15, ///<"nsamples" - counts the number of times the signal apear in the time window
	FTR_EXISTS = 16, ///<"exists" - boolean 0/1 if the signal apears in the time window
	FTR_CATEGORY_SET_FIRST = 17, ///<"category_set_first" - boolean 0/1 if the signal apears in the time window and did not appear ever before the window
	FTR_MAX_DIFF = 18, ///<maximum diff in window
	FTR_FIRST_DAYS = 19, ///< time diffrence from prediction time to first time with signal
	FTR_LAST
} BasicFeatureTypes;

/** @enum
* TimeRangeTypes determines how the time window depends on a time-range determined by another signal
*/
typedef enum {
	TIME_RANGE_CURRENT = 0, ///<"current" - consider only the current time-range
	TIME_RANGE_BEFORE = 1, ///< "before" - consider anything before the current time-range
	TIME_RANGE_LAST
} TimeRangeTypes ;

/**
* A Basic Stats Generator for calcing simple statics on time window
*/
class BasicFeatGenerator : public FeatureGenerator {
private:
	// actual generators
	float uget_last(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime); // Added the win as needed to be called on different ones in uget_win_delta
	float uget_first(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_last2(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_avg(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_max(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_min(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_std(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_last_delta(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_last_time(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_last2_time(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_slope(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_win_delta(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int _d_win_from, int _d_win_to, int outcomeTime);
	float uget_category_set(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_category_set_count(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_category_set_sum(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_nsamples(UniversalSigVec &usv, int time, int _win_from, int _win_to, int outcomeTime);
	float uget_exists(UniversalSigVec &usv, int time, int _win_from, int _win_to, int outcomeTime);
	float uget_max_diff(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_first_time(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_category_set_first(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);

	// update time window according to time-range signal
	void get_updated_time_window(UniversalSigVec& time_range_usv, TimeRangeTypes type, int time, int& updated_win_from, int& updated_win_to, bool delta_win,
		int& updated_d_win_from, int& updated_d_win_to);
	void get_updated_time_window(TimeRangeTypes type, int range_from, int range_to, int time, int _win_from, int _win_to, int& updated_win_from, int& updated_win_to);

public:
	// Feature Descrption
	string signalName;
	int signalId;

	int version() { return 2; } ///< added "bound_outcomeTime" in version 1
								///< added time_range_signal in version 2

	// Signal to determine allowed time-range (e.g. current stay/admission for inpatients)
	string timeRangeSignalName = "";
	int timeRangeSignalId;
	TimeRangeTypes timeRangeType = TIME_RANGE_CURRENT;
	int time_unit_range_sig = MedTime::Undefined;		///< the time init in which the range signal is given. (set correctly from Repository in learn and _generate)

	// parameters (should be serialized)
	BasicFeatureTypes type = FTR_LAST;
	int win_from = 0;///< time window for feature: win_from is the minimal time before from the prediction time
	int win_to = 360000;///< time window for feature: win_to is the maximal time before the prediction time			 
	int d_win_from = 360; ///< delta time window for the FTR_WIN_DELTA_VALUE feature. the second time window
	int d_win_to = 360000;	///< delta time window for the FTR_WIN_DELTA_VALUE feature. the second time window
	int time_unit_win = MedTime::Undefined;			///< the time unit in which the windows are given. Default: Undefined
	int time_channel = 0;						///< n >= 0 : use time channel n , default: 0.
	int val_channel = 0;						///< n >= 0 : use val channel n , default : 0.
	int sum_channel = 1;						///< for FTR_CETEGORY_SET_SUM
	vector<string> sets;						///< for FTR_CATEGORY_SET_* , the list of sets 
	int time_unit_sig = MedTime::Undefined;		///< the time init in which the signal is given. (set correctly from Repository in learn and _generate)
	string in_set_name = "";					///< set name (if not given - take list of members)
	bool bound_outcomeTime; ///< If true will truncate time window till outcomeTime

	// helpers
	vector<char> lut;							///< to be used when generating FTR_CATEGORY_SET_*

	// Naming 
	void set_names();

	// Constructor/Destructor
	BasicFeatGenerator() : FeatureGenerator() { init_defaults(); };
	//~BasicFeatGenerator() {};
	void set(string& _signalName, BasicFeatureTypes _type) { 
		set(_signalName, _type, 0, 360000); 
		req_signals.assign(1, signalName); 
		if (timeRangeSignalName != "")
			req_signals.push_back(timeRangeSignalName);
	}
	
	void set(string& _signalName, BasicFeatureTypes _type, int _time_win_from, int _time_win_to) {
		signalName = _signalName; type = _type; win_from = _time_win_from; win_to = _time_win_to;
		set_names();
		req_signals.assign(1, signalName);
		if (timeRangeSignalName != "")
			req_signals.push_back(timeRangeSignalName);
	}

	/// Converts a name to type - please reffer to BasicFeatureTypes
	BasicFeatureTypes name_to_type(const string &name);

	/// Conversion between time-range type and name
	TimeRangeTypes time_range_name_to_type(const string& name);
	string time_range_type_to_name(TimeRangeTypes type);

	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp BasicFeatGenerator::init
	int init(map<string, string>& mapper);
	void init_defaults();

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<BasicFeatGenerator *>(generator)); }

	/// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) {
		time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(signalName)].time_unit;
		if (timeRangeSignalName != "")
			time_unit_range_sig = rep.sigs.Sid2Info[rep.sigs.sid(timeRangeSignalName)].time_unit;
		return 0;
	}

	/// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);
	float get_value(PidDynamicRec& rec, int index, int time, int outcomeTime);

	/// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); timeRangeSignalId = dict.id(timeRangeSignalName); }

	/// Init required tables
	void init_tables(MedDictionarySections& dict);

	// Serialization
	ADD_SERIALIZATION_FUNCS(generator_type, type, tags, serial_id, win_from, win_to, d_win_from, d_win_to,
		time_unit_win, time_channel, val_channel, sum_channel, signalName, sets,
		names, req_signals, in_set_name ,bound_outcomeTime, timeRangeSignalName, timeRangeType)

};

/**
* Age Generator
*/
class AgeGenerator : public FeatureGenerator {
public:

	string signalName;
	/// Signal Id
	int signalId;

	~AgeGenerator() { clear(); }
	void clear() { }

	// Constructor/Destructor
	AgeGenerator() {
		generator_type = FTR_GEN_AGE; names.push_back("Age"); signalId = -1; signalName = "BYEAR"; req_signals.assign(1, signalName);
	}
	//~AgeGenerator() {};

	// Name
	void set_names() { if (names.empty()) names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".Age"); tags.push_back("Age"); }

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<AgeGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }

	// Serialization
	int version() { return 1; }
	size_t get_size() { return MedSerialize::get_size(generator_type, names, tags, iGenerateWeights, signalName, req_signals); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, names, tags, iGenerateWeights, signalName, req_signals); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, names, tags, iGenerateWeights, signalName, req_signals); }
	virtual int init(map<string, string>& mapper);
};

/**
* Singleton
*/
class SingletonGenerator : public FeatureGenerator {
private:
	vector<char> lut;			///< to be used when generating sets*
public:

	/// Signal Id
	string signalName;
	int signalId;

	vector<string> sets = {};		/// list of sets 
	string in_set_name = "";

	// Constructor/Destructor
	SingletonGenerator() : FeatureGenerator() { generator_type = FTR_GEN_SINGLETON; names.push_back(signalName); signalId = -1; req_signals.assign(1, signalName); }
	SingletonGenerator(int _signalId) : FeatureGenerator() { generator_type = FTR_GEN_SINGLETON; names.push_back(signalName); signalId = _signalId; req_signals.assign(1, signalName); }

	// Name
	void set_names();

	// Init LUT for categorial variable
	void init_tables(MedDictionarySections& dict);
	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp SingletonGenerator::init
	int init(map<string, string>& mapper);

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<SingletonGenerator *>(generator)); }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.assign(1, dict.id(signalName)); }

	// Serialization
	ADD_SERIALIZATION_FUNCS(generator_type, req_signals, signalName, names, tags, iGenerateWeights, sets, lut)
};


/**
* Gender
*/
class GenderGenerator : public FeatureGenerator {
public:

	/// Gender Id
	int genderId;

	// Constructor/Destructor
	GenderGenerator() : FeatureGenerator() { generator_type = FTR_GEN_GENDER; names.push_back("Gender"); genderId = -1; req_signals.assign(1, "GENDER"); }
	GenderGenerator(int _genderId) : FeatureGenerator() { generator_type = FTR_GEN_GENDER; names.push_back("Gender"); genderId = _genderId; req_signals.assign(1, "GENDER"); }

	//~GenderGenerator() {};

	// Name
	void set_names() { if (names.empty()) names.push_back("Gender"); tags.push_back("Gender"); }

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<GenderGenerator *>(generator)); }

	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp SingletonGenerator::init
	int init(map<string, string>& mapper);

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { genderId = dict.id("GENDER"); }
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.assign(1, dict.id("GENDER")); }

	// Serialization
	size_t get_size() { return MedSerialize::get_size(generator_type, names, tags, iGenerateWeights); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, generator_type, names, tags, iGenerateWeights); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, generator_type, names, tags, iGenerateWeights); }
};

/**
* BinnedLinearModels : parameters
*/
struct BinnedLmEstimatesParams {
	vector<int> bin_bounds;
	int min_period;
	int max_period;
	float rfactor;

	vector<int> estimation_points;

};

/**
* BinnedLinearModels : Apply a set of liner models to generate features
*/
class BinnedLmEstimates : public FeatureGenerator {
public:
	// Feature Descrption
	string signalName;
	int signalId, byearId, genderId;

	BinnedLmEstimatesParams params;
	vector<MedLM> models;
	vector<float> xmeans, xsdvs, ymeans, ysdvs;
	vector<float> means[2];

	int time_unit_periods = MedTime::Undefined;		///< the time unit in which the periods are given. Default: Undefined
	int time_unit_sig = MedTime::Undefined;			///< the time init in which the signal is given. Default: Undefined
	int time_channel = 0;						///< n >= 0 : use time channel n , default: 0.
	int val_channel = 0;						///< n >= 0 : use val channel n , default : 0.

	/// Naming 
	void set_names();

	// Constructor/Destructor
	BinnedLmEstimates() : FeatureGenerator() { signalName = ""; init_defaults(); };
	BinnedLmEstimates(string _signalName) : FeatureGenerator() { signalName = _signalName; init_defaults(); req_signals.push_back(signalName); names.clear();  set_names(); };
	BinnedLmEstimates(string _signalName, string init_string) : FeatureGenerator() { signalName = _signalName; init_defaults(); req_signals.push_back(signalName); init_from_string(init_string); };

	//~BinnedLmEstimates() {};

	void set(string& _signalName);
	void set(string& _signalName, BinnedLmEstimatesParams* _params);

	void init_defaults();
	/// The parsed fields from init command.
	/// @snippet BinnedLmEstimates.cpp BinnedLmEstimates::init
	int init(map<string, string>& mapper);

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<BinnedLmEstimates *>(generator)); }

	/// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors);

	/// generate new feature(s)
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	/// Filter generated features according to a set. return number of valid features (does not affect single-feature genertors, just returns 1/0 if feature name in set)
	int filter_features(unordered_set<string>& validFeatures);

	// get pointers to data
	void get_p_data(MedFeatures& features);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict);

	// Age Related functions
	void prepare_for_age(PidDynamicRec& rec, UniversalSigVec& ageUsv, int &age, int &byear);
	void prepare_for_age(MedPidRepository& rep, int id, UniversalSigVec& ageUsv, int &age, int &byear);
	inline void get_age(int time, int time_unit_from, int& age, int byear);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// print 
	void print();
};


/** @enum
* RangeFeatGenerator types
*/
typedef enum {
	FTR_RANGE_CURRENT = 0, ///<"current" - finds the value of the time range signal that intersect with win_from. signal start_time is before this time and signal end_time is after this time point
	FTR_RANGE_LATEST = 1, ///<"latest" - finds the last value of the time range signal, that there is intersection of time signal range with the defined time window
	FTR_RANGE_MAX = 2, ///<"max" - finds the maximal value of the time range signal, that there is intersection of time signal range with the defined time window
	FTR_RANGE_MIN = 3, ///<"min" - finds the minimal value of the time range signal, that there is intersection of time signal range with the defined time window
	FTR_RANGE_EVER = 4,///<"ever" - boolean 0/1 - finds if there is intersection between signal time window and the defined time window with specific lut value. uses set.
	///"time_diff" - returns time diffrences between first intersection(if check_first is True) between signal time window and the defined time window with specific lut value. uses set.
	///if check_first is false returns the time diffrences between last intersection between signal time window and the defined time window. prediction time minus the last intersecting signal end time window.
	///if the last intersction if time ranges has no match to sets value and check_first is false will return -win_to value, otherwise missing value
	FTR_RANGE_TIME_DIFF = 5,
	FTR_RANGE_LAST
} RangeFeatureTypes;

/**
* RangeFeatGenerator : Generate features for a time range with value signal (for example drug)
*/
class RangeFeatGenerator : public FeatureGenerator {
private:
	// actual generators
	float uget_range_current(UniversalSigVec &usv, int time_point);
	float uget_range_latest(UniversalSigVec &usv, int time_point);
	float uget_range_min(UniversalSigVec &usv, int time_point);
	float uget_range_max(UniversalSigVec &usv, int time_point);
	float uget_range_ever(UniversalSigVec &usv, int time_point);
	float uget_range_time_diff(UniversalSigVec &usv, int time_point);
public:

	string signalName; ///< Signal to consider
	int signalId;
	vector<string> sets;						///< FTR_RANGE_EVER checks if the signal ever was in one of these sets/defs from the respective dict
	RangeFeatureTypes type; ///< Type of comorbidity index to generate
	int win_from = 0; ///< time window for feature: from is the minimal time before prediciton time
	int win_to = 360000;			///< time window for feature: to is the maximal time before prediciton time
	int time_unit_win = MedTime::Undefined;			///< the time unit in which the windows are given. Default: Undefined
	int time_unit_sig = MedTime::Undefined;		///< the time init in which the signal is given. (set correctly from Repository in learn and Generate)
	int val_channel = 0;						///< n >= 0 : use val channel n , default : 0.
	int check_first = 1;						///< if 1 choose first occurance of check_val otherwise choose last


	vector<char> lut;							///< to be used when generating FTR_RANGE_EVER



	// Constructor/Destructor
	RangeFeatGenerator() : FeatureGenerator() { init_defaults(); };
	//~RangeFeatGenerator() {};
	void set(string& _signalName, RangeFeatureTypes _type) { set(_signalName, _type, 0, 360000); req_signals.assign(1, signalName); }
	void set(string& _signalName, RangeFeatureTypes _type, int _time_win_from, int _time_win_to) {
		signalName = _signalName; type = _type; win_from = _time_win_from; win_to = _time_win_to;
		set_names(); req_signals.assign(1, signalName);
	}

	// Naming 
	void set_names();

	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp RangeFeatGenerator::init
	int init(map<string, string>& mapper);
	void init_defaults();
	RangeFeatureTypes name_to_type(const string &name); ///< please reffer to RangeFeatureTypes to understand the options
	void init_tables(MedDictionarySections& dict);
	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<RangeFeatGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(signalName)].time_unit; return 0; }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);
	float get_value(PidDynamicRec& rec, int index, int date);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }


	// Serialization
	// Serialization
	virtual int version() { return  1; };	// ihadanny 20171206 - added sets
	ADD_SERIALIZATION_FUNCS(generator_type, signalName, type, win_from, win_to, val_channel, names, tags, req_signals, sets, check_first)
};

/**
* Use a model to generate predictions to be used as features
*/
class ModelFeatGenerator : public FeatureGenerator {
public:

	string modelFile = ""; ///<  File for serialized model
	MedModel *model = NULL; ///< model
	string modelName = ""; ///< name of final feature
	int n_preds = 1;  ///< how many features to create
	int impute_existing_feature = 0; ///< If true will impute using model feature, otherwise using preds
	int use_overriden_predictions = 0;
	string medSamples_path = ""; ///< If provided will override predictions using samples
	/// A container for the predictions

	/// Naming 
	void set_names();

	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp ModelFeatGenerator::init
	int init(map<string, string>& mapper);
	int init_from_model(MedModel *_model);

	/// Hack - Instead of actually predicting with the model, load predictions from a MedSamples object. Compare to the models MedSamples (unless empty)
	void override_predictions(MedSamples& inSamples, MedSamples& modelSamples);

	/// Do the actual prediction prior to feature generation ...
	void prepare(MedFeatures & features, MedPidRepository& rep, MedSamples& samples);

	/// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// (De)Serialize
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	//dctor:
	~ModelFeatGenerator();
private:
	vector<float> preds;
	MedSamples _preloaded;

	void find_predictions(MedSamples& requestSamples, MedFeatures &features);
};


/**
* Time Feature Generator: creating sample-time features (e.g. differentiate between times of day, season of year, days of the week, etc.)
*/

typedef enum {
	FTR_TIME_YEAR = 0, ///< Year (as is)
	FTR_TIME_MONTH = 1, ///< Month of year (0-11)
	FTR_TIME_DAY_IN_MONTH = 2, ///< Day of the month (0-30)
	FTR_TIME_DAY_IN_WEEK = 3, ///< Day of the week (0-6)
	FTR_TIME_HOUR = 4, ///< Hour of the day (0-23)
	FTR_TIME_MINUTE = 5, ///< Minute of the hout (0-59)
	FTR_TIME_LAST,
} TimeFeatTypes;


class TimeFeatGenerator : public FeatureGenerator {
public:

	// Time Unit
	TimeFeatTypes time_unit = FTR_TIME_LAST;

	// Binning of time units
	vector<int> time_bins;
	vector<string> time_bin_names;
	
	// Constructor/Destructor
	TimeFeatGenerator() { generator_type = FTR_GEN_TIME; }
	~TimeFeatGenerator() {}

	// Naming 
	void set_names();

	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp TimeFeatGenerator::init
	int init(map<string, string>& mapper);
	int get_time_unit(string name);
	int get_time_bins(string& binsInfo);
	int get_nBins();
	void set_default_bins();
	string time_unit_to_string(TimeFeatTypes time_unit);

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<TimeFeatGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { return 0; }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);

	// Serialization
	ADD_SERIALIZATION_FUNCS(names, time_unit, time_bins, time_bin_names)
};


//=======================================
// Helpers
//=======================================

/// gets a [-_win_to, -_win_from] window in win time unit, and returns [_min_time, _max_time] window in signal time units relative to _win_time
void get_window_in_sig_time(int _win_from, int _win_to, int _time_unit_win, int _time_unit_sig, int _win_time, int &_min_time, int &_max_time);

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(BasicFeatGenerator)
MEDSERIALIZE_SUPPORT(AgeGenerator)
MEDSERIALIZE_SUPPORT(GenderGenerator)
MEDSERIALIZE_SUPPORT(BinnedLmEstimates)
MEDSERIALIZE_SUPPORT(RangeFeatGenerator)

#endif