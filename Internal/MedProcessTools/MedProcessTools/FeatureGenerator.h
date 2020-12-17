#ifndef _FTR_GENERATOR_H_
#define _FTR_GENERATOR_H_

#include <InfraMed/InfraMed/InfraMed.h>
#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/RepProcess.h>
#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedModelExceptions.h>
#include <MedTime/MedTime/MedTime.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedAlgo/MedAlgo/MedLM.h>
#include <cfloat>
#include <boost/regex.hpp>

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
	FTR_GEN_UNIFIED_SMOKING, ///< "unified_smoking" - creating smoking feature - UnifiedSmokingGenerator
	FTR_GEN_RANGE, ///< "range" - creating RangeFeatGenerator
	FTR_GEN_DRG_INTAKE, ///< "drugIntake" - creating drugs feature coverage of prescription time - DrugIntakeGenerator
	FTR_GEN_ALCOHOL, ///< "alcohol" - creating alcohol feature - AlcoholGenerator
	FTR_GEN_MODEL, ///< "model" - creating ModelFeatGenerator
	FTR_GEN_TIME, ///< "time" - creating sample-time features (e.g. differentiate between times of day, season of year, days of the week, etc.). Creates TimeFeatGenerator
	FTR_GEN_ATTR, ///< "attr" - creating features from samples attributes. Creates AttrFeatGenerator
	FTR_GEN_CATEGORY_DEPEND, ///< "category_depend" - creates features from categorical signal that have statistical strength in samples - CategoryDependencyGenerator
	FTR_GEN_EMBEDDING, ///< "embedding" - allows applying a pre trained embedding model to incorporate features into matrix. Creates EmbeddingGenerator
	FTR_GEN_EXTRACT_TBL, ///< "extract_tbl" - extract values from table with keys and rules to join with each patient. Creates FeatureGenExtractTable
	FTR_GEN_ELIXHAUSER, ///< Calculate Current Elixhauser given latest DRG and Diagnosis information. Creates ElixhauserGenerator
	FTR_GEN_DIABETES_FINDER, ///< "diabetes_finder" - Diabetes Finder feature. Creates DiabetesFinderGenerator
	FTR_GEN_LAST
} FeatureGeneratorTypes;

/** @file
* FeatureGenerator : creating features from raw signals
*/
class FeatureGenerator : public SerializableObject {
public:

	/// Type
	FeatureGeneratorTypes generator_type = FTR_GEN_LAST;

	/// Feature name
	vector<string> names;

	// Threading
	int learn_nthreads = 16, pred_nthreads = 16;

	/// Missing value
	float missing_val = (float)MED_MAT_MISSING_VALUE;

	/// Tags - for defining labels or groups. may be used later for filtering for example
	vector<string> tags;

	/// Feature/Weights generator
	int iGenerateWeights = 0;

	// Naming
	virtual void set_names() { names.clear(); }

	// Helper - pointers to data vectors in MedFeatures (to save time in generation)
	vector <float *> p_data;

	// Prepare for feature generation
	virtual void prepare(MedFeatures &features, MedPidRepository& rep, MedSamples& samples);
	virtual void get_p_data(MedFeatures& features, vector<float *> &_p_data);
	void get_p_data(MedFeatures& features) { get_p_data(features, p_data); }

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

	// generated features
	virtual void get_generated_features(unordered_set<string>& names_list) { for (auto &s : names) names_list.insert(s); }

	// Signal Ids
	virtual void set_signal_ids(MedSignals& sigs) { return; }

	// Init required tables
	virtual void init_tables(MedDictionarySections& dict) { return; }

	// Learn a generator
	virtual int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) { return 0; }
	int learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) { set_names(); return _learn(rep, samples, processors); }
	int learn(MedPidRepository& rep, const MedSamples& samples) { set_names(); return _learn(rep, samples, vector<RepProcessor *>()); }

	// generate feature data from repository
	// We assume the corresponding MedSamples have been inserted to MedFeatures : either at the end or at position index
	int _generate(PidDynamicRec& in_rep, MedFeatures& features, int index, int num) {return _generate(in_rep, features, index, num, p_data); }

	// the following is the MAIN generation routine to implement.
	// note that it is given a p_data of its own. This is in order to allow different records to write results to different places.
	// the default run will use it with the generator p_data.
	virtual int _generate(PidDynamicRec& in_rep, MedFeatures& features, int index, int num, vector<float *> &_p_data) { return 0; }

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

	/// returns for each used signal it's used categories
	virtual void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const {};

	// Filter generated features according to a set. return number of valid features
	virtual int filter_features(unordered_set<string>& validFeatures);

	///<summary>
	/// prints summary of generator job. optional, called after generate.
	/// for example - prints how many values were missing value
	///</summary>
	virtual void make_summary() {};

	// Serialization
	ADD_CLASS_NAME(FeatureGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, names, learn_nthreads, pred_nthreads, missing_val, tags, iGenerateWeights)
	void *new_polymorphic(string derived_class_name);

	size_t get_generator_size();
	size_t generator_serialize(unsigned char *blob);

	virtual void print() { fprintf(stderr, "Print Not Implemented for feature\n"); }


	// debug print for a feature generator. fg_flag can 
	virtual void dprint(const string &pref, int fg_flag);


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
	FTR_LAST_DAYS = 8, ///<"last_time" - time diffrence from prediction time to last time has signal in range of values
	FTR_LAST2_DAYS = 9,///<"last_time2" - time diffrence from prediction time to one previous last time has signal in range of values
	FTR_SLOPE_VALUE = 10, ///<"slope" - calculating the slope over the points in the window
	FTR_WIN_DELTA_VALUE = 11, ///<"win_delta" - diffrence in value in two time windows (only if both exists, otherwise missing_value). value in [win_from,win_to] minus value in [d_win_from, d_win_to]
	FTR_CATEGORY_SET = 12, ///<"category_set" - boolean 0/1 if the signal has the value in the given lut (which initialized by the "sets" that can be specific single definition or name of set definition. the lookup is hierarchical)
	FTR_CATEGORY_SET_COUNT = 13,///<"category_set_count" - counts the number of appearnces of sets in the time window
	FTR_CATEGORY_SET_SUM = 14, ///<"category_set_sum" - sums the values of appearnces of sets in the time window
	FTR_NSAMPLES = 15, ///<"nsamples" - counts the number of times the signal apear in the time window
	FTR_EXISTS = 16, ///<"exists" - boolean 0/1 if the signal apears in the time window
	FTR_CATEGORY_SET_FIRST = 17, ///<"category_set_first" - boolean 0/1 if the signal apears in the time window and did not appear ever before the window
	FTR_MAX_DIFF = 18, ///<"max_diff" maximum diff in window
	FTR_FIRST_DAYS = 19, ///<"first_time" time diffrence from prediction time to first time with signal
	FTR_RANGE_WIDTH = 20, ///<"range_width" maximal value - minimal value in a given window time frame
	FTR_CATEGORY_SET_FIRST_TIME = 21, ///<"category_set_first_time" - first time of category set found in the time window
	FTR_SUM_VALUE=22, ///<"sum" - sum of values in window
	FTR_LAST_NTH_VALUE = 23, ///<"last_nth" : (set also N_th parameter to use), get the last N_th in window, 0 is last, 1 is last2, etc.
	FTR_CATEGORY_SET_LAST_NTH = 24, ///<"category_set_last_nth" : (set also N_th parameter to use), check is the last N_th in window is in the given set
	FTR_TIME_SINCE_LAST_CHANGE = 25, ///<"time_since_last_change" : go over states signal, take last time since the value changed
	FTR_LAST
} BasicFeatureTypes;

/** @enum
* TimeRangeTypes determines how the time window depends on a time-range determined by another signal
*/
typedef enum {
	TIME_RANGE_CURRENT = 0, ///<"current" - consider only the current time-range
	TIME_RANGE_BEFORE = 1, ///< "before" - consider anything before the current time-range
	TIME_RANGE_LAST
} TimeRangeTypes;

/**
* A Basic Stats Generator for calcing simple statics on time window
*/
class BasicFeatGenerator : public FeatureGenerator {
private:
	// actual generators
	float uget_last(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime); // Added the win as needed to be called on different ones in uget_win_delta
	float uget_last_nth(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_first(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_last2(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_avg(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_max(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_min(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_sum(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_std(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_last_delta(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_last_time(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_last2_time(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_slope(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_win_delta(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int _d_win_from, int _d_win_to, int outcomeTime);
	float uget_category_set(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_category_set_last_nth(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_category_set_count(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_category_set_sum(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_nsamples(UniversalSigVec &usv, int time, int _win_from, int _win_to, int outcomeTime);
	float uget_exists(UniversalSigVec &usv, int time, int _win_from, int _win_to, int outcomeTime);
	float uget_range_width(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_max_diff(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_first_time(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_category_set_first(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_category_set_first_time(PidDynamicRec &rec, UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);
	float uget_time_since_last_change(UniversalSigVec &usv, int time_point, int _win_from, int _win_to, int outcomeTime);

	// Applying non FTR_CATEGORY_SET_* types to categorical data
	unordered_set<int> categ_require_dict = { FTR_LAST_VALUE, FTR_FIRST_VALUE, FTR_LAST2_VALUE, FTR_LAST_NTH_VALUE }; // Types that requriew dictionary if applied on categorical data
	unordered_set<int> categ_forbidden = { FTR_AVG_VALUE, FTR_MAX_VALUE, FTR_MIN_VALUE,  FTR_STD_VALUE, FTR_LAST_DELTA_VALUE,
		FTR_SLOPE_VALUE,FTR_WIN_DELTA_VALUE, FTR_MAX_DIFF,FTR_RANGE_WIDTH, FTR_SUM_VALUE }; // Types that are not allowed for categorical data
	map<string, int> categ_value2id;

public:
	// Feature Descrption
	string signalName;
	int signalId;

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
	bool bound_outcomeTime = false; ///< If true will truncate time window till outcomeTime
	float min_value = -FLT_MAX, max_value = FLT_MAX; ///< values range for FTR_LAST(2)_DAYS
	int N_th = 0; ///< used in last_nth and category_set_last_nth 
	int zero_missing = 0; ///< in some cases of category_set (or others) we may want to get 0 instead of missing_value, turn this on for that

	// helpers
	vector<char> lut;							///< to be used when generating FTR_CATEGORY_SET_*
	vector<float> categ_map;					///< to be used when applying non FTR_CATEGORY_SET_* types to categorical data
	bool apply_categ_map;

	// Naming 
	void set_names();

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

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
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) {
		time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(signalName)].time_unit;
		if (timeRangeSignalName != "")
			time_unit_range_sig = rep.sigs.Sid2Info[rep.sigs.sid(timeRangeSignalName)].time_unit;
		return 0;
	}

	/// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);
	float get_value(PidDynamicRec& rec, int index, int time, int outcomeTime);

	/// Signal Ids
	void set_signal_ids(MedSignals& sigs) { signalId = sigs.sid(signalName); timeRangeSignalId = sigs.sid(timeRangeSignalName); }

	/// Init required tables
	void init_tables(MedDictionarySections& dict);

	void prepare(MedFeatures &features, MedPidRepository& rep, MedSamples& samples);

	// Serialization
	ADD_CLASS_NAME(BasicFeatGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, type, tags, serial_id, win_from, win_to, d_win_from, d_win_to, time_unit_win, time_channel, val_channel, sum_channel, min_value, max_value, signalName, sets,
			names, req_signals, in_set_name, bound_outcomeTime, timeRangeSignalName, timeRangeType, time_unit_sig, N_th, zero_missing, missing_val, categ_value2id)

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
		generator_type = FTR_GEN_AGE; names.push_back("Age"); signalId = -1; signalName = "BDATE"; req_signals.assign(1, signalName);
	}
	//~AgeGenerator() {};

	// Name
	void set_names() { if (names.empty()) names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".Age"); tags.push_back("Age"); }

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<AgeGenerator *>(generator)); }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// Signal Ids
	void set_signal_ids(MedSignals& sigs) { signalId = sigs.sid(signalName); }

	// Serialization
	ADD_CLASS_NAME(AgeGenerator)
		ADD_SERIALIZATION_FUNCS(generator_type, names, tags, iGenerateWeights, signalName, req_signals)

		virtual int init(map<string, string>& mapper);
};

/**
* Singleton
*/
class SingletonGenerator : public FeatureGenerator {
private:
	vector<char> lut;			///< to be used when generating sets*
	unordered_map<string, float> name2Value; ///< Used for mapping dictionary strings to values (we don't rely on dictionary not to change)
	vector<float> id2Value; ///< mapping of dictionary id to value (rebuilt according to dictionary + name2Value)

	void get_id2Value(MedDictionarySections& dict);
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

	// learn generator (learning name2Value)
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors);

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// Preparation - just fill the value2Name attribute
	void prepare(MedFeatures &features, MedPidRepository& rep, MedSamples& samples);

	// Signal Ids
	void set_signal_ids(MedSignals& sigs) { signalId = sigs.sid(signalName); }
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.assign(1, dict.id(signalName)); }

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

	// Serialization
	ADD_CLASS_NAME(SingletonGenerator)
		ADD_SERIALIZATION_FUNCS(generator_type, req_signals, signalName, names, tags, iGenerateWeights, sets, lut, name2Value)
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
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// Signal Ids
	void set_signal_ids(MedSignals& sigs) { genderId = sigs.sid("GENDER"); }
	void set_required_signal_ids(MedDictionarySections& dict) { req_signal_ids.assign(1, dict.id("GENDER")); }

	// Serialization
	ADD_CLASS_NAME(GenderGenerator)
		ADD_SERIALIZATION_FUNCS(generator_type, names, tags, iGenerateWeights)
};

/**
* BinnedLinearModels : parameters
*/
struct BinnedLmEstimatesParams : public SerializableObject {
	vector<int> bin_bounds;
	int min_period;
	int max_period;
	float rfactor;

	vector<int> estimation_points;
	ADD_CLASS_NAME(BinnedLmEstimatesParams)
		ADD_SERIALIZATION_FUNCS(bin_bounds, min_period, max_period, rfactor, estimation_points)

};

/**
* BinnedLinearModels : which time-points to take
*/

typedef enum {
	BINNED_LM_TAKE_ALL = 0,
	BINNED_LM_STOP_AT_FIRST = 1,
	BINNED_LM_STOP_AT_LAST = 2,
	BINNED_LM_LAST
} BinnedLMSamplingStrategy;

/**
* BinnedLinearModels : Apply a set of liner models to generate features
*/
class BinnedLmEstimates : public FeatureGenerator {
public:
	// Feature Descrption
	string signalName;
	int signalId, bdateId, genderId;

	BinnedLmEstimatesParams params;
	BinnedLMSamplingStrategy sampling_strategy = BINNED_LM_TAKE_ALL;
	vector<MedLM> models;
	vector<float> xmeans, xsdvs, ymeans, ysdvs;
	vector<vector<float>> means = { {}, {} };

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
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors);

	/// generate new feature(s)
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	/// Filter generated features according to a set. return number of valid features (does not affect single-feature genertors, just returns 1/0 if feature name in set)
	int filter_features(unordered_set<string>& validFeatures);

	// get pointers to data
	void get_p_data(MedFeatures& features, vector<float *> &_p_data);

	// Signal Ids
	void set_signal_ids(MedSignals& sigs);

	// Sampling strategy
	void set_sampling_strategy(string& strategy);

	// Age Related functions
	void prepare_for_age(PidDynamicRec& rec, UniversalSigVec& ageUsv, int &age, int &byear);
	void prepare_for_age(MedPidRepository& rep, int id, UniversalSigVec& ageUsv, int &age, int &byear);
	inline void get_age(int time, int time_unit_from, int& age, int byear);

	void dprint(const string &pref, int fg_flag);

	// Serialization
	ADD_CLASS_NAME(BinnedLmEstimates)
		ADD_SERIALIZATION_FUNCS(generator_type, signalName, names, tags, req_signals, time_unit_periods, iGenerateWeights, params, xmeans, xsdvs, ymeans, means, models, time_unit_sig, sampling_strategy)

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
	///<"recurrence_count" - count the number of time the event occur shortly after a previous event, there is an intersection of the time signal range with the defined time window
	///previous event does not need to intersect the time window. 
	FTR_RANGE_RECURRENCE_COUNT = 6,
	/// "time_covered"  : give a time window, sum up all the times in ranges that intersect the time window
	FTR_RANGE_TIME_COVERED = 7,
	/// "last_nth_time_len"  : gives the length (in win_time_unit) of the last_n range in the window. If in middle of range, cuts to current time
	FTR_RANGE_LAST_NTH_TIME_LENGTH = 8,
	FTR_RANGE_TIME_DIFF_START = 9,
	FTR_RANGE_TIME_INSIDE = 10, ///<< "time_inside" : checks if the prediction time point is currently INSIDE a range, if not returns 0, if it is , then how long since the start.
	FTR_RANGE_LAST
} RangeFeatureTypes;

/**
* RangeFeatGenerator : Generate features for a time range with value signal (for example drug)
*/
class RangeFeatGenerator : public FeatureGenerator {
private:
	// actual generators
	float uget_range_current(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_latest(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_min(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_max(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_ever(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_time_diff(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_recurrence_count(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_time_covered(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_last_nth_time_len(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_time_diff_start(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);
	float uget_range_time_inside(UniversalSigVec &usv, int updated_win_from, int updated_win_to, int time);

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
	float div_factor = 1.0f;					/// dividing by this number in time_covered option

	vector<char> lut;							///< to be used when generating FTR_RANGE_EVER
	int recurrence_delta = 30 * 24 * 60;		///< maximum time for a subsequent range signal to be considered a recurrence in in window time units
	int min_range_time = -1;					///< if different from -1, the minimum length for a range to be considered valid in window time units (else not checked)
	int N_th = 0;								///< the index of the N-th range in order to consider in the last_nth_time_len option
	int zero_missing = 0;	///< in some cases we may want to get 0 instead of missing values
	int strict_times = 0;  ///< if on , will ignore cases in which the second time channel is after the prediction time
	int conditional_channel = -1; ///< in some cases (currently last_nth_len, and time_covered) we allow doing the calculation only on ranges passing the condition of being included in sets in this channel
	bool regex_on_sets= false;        ///< if on , regex is applied on .*sets[i].* and aggregated. 

	int first_evidence_time_channel = 1;	 ///< sometimes we have a different time channel stating WHEN the range started. We are strict and use the default of last time in range, but sometimes it is not like that and this can allow calculating if we are now IN the range, and how LONG since the start.

	// Signal to determine allowed time-range (e.g. current stay/admission for inpatients)
	string timeRangeSignalName = "";
	int timeRangeSignalId;
	TimeRangeTypes timeRangeType = TIME_RANGE_CURRENT;
	int time_unit_range_sig = MedTime::Undefined;		///< the time unit in which the range signal is given. (set correctly from Repository in learn and _generate)



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

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp RangeFeatGenerator::init
	int init(map<string, string>& mapper);
	void init_defaults();
	RangeFeatureTypes name_to_type(const string &name); ///< please reffer to RangeFeatureTypes to understand the options
	void init_tables(MedDictionarySections& dict);
	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<RangeFeatGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) { time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(signalName)].time_unit; return 0; }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);
	float get_value(PidDynamicRec& rec, int index, int date);

	// Signal Ids
	void set_signal_ids(MedSignals& sigs) { signalId = sigs.sid(signalName); }


	// Serialization
	ADD_CLASS_NAME(RangeFeatGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, signalName, type, win_from, win_to, val_channel, names, tags, req_signals, sets, check_first, timeRangeSignalName, timeRangeType, recurrence_delta, min_range_time,
			time_unit_sig, time_unit_win, div_factor, N_th, zero_missing, strict_times,conditional_channel, regex_on_sets)
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
	int impute_existing_feature = 0; ///< If true will use model to impute an existing feature (determined by model name. Otherwise - generate new feature(s)
	int use_overriden_predictions = 0; ///< Use a given vector of predictions instead of applying model
	int time_unit_win = global_default_windows_time_unit; ///< the time unit in which the times are given. Default: global_default_windows_time_unit
	int time_unit_sig = global_default_windows_time_unit; ///< the time init in which the signal is given. (set correctly from Repository in learn and Generate)
	vector<int> times;

	/// Naming 
	void set_names();

	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp ModelFeatGenerator::init
	int init(map<string, string>& mapper);
	int init_from_model(MedModel *_model);

	/// Use a given vector of predictions instead of applying model
	void override_predictions(MedSamples& inSamples, MedSamples& modelSamples);

	/// Do the actual prediction prior to feature generation ...
	void prepare(MedFeatures & features, MedPidRepository& rep, MedSamples& samples);

	/// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	void modifySampleTime(MedSamples& samples, int time);

	// (De)Serialize
	ADD_CLASS_NAME(ModelFeatGenerator)
		ADD_SERIALIZATION_HEADERS()

		//dctor:
		~ModelFeatGenerator();
private:
	vector<vector<vector<float>>> preds;
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
	FTR_TIME_DATE = 6, ///< Completete date (as is)
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

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// Serialization
	ADD_CLASS_NAME(TimeFeatGenerator)
		ADD_SERIALIZATION_FUNCS(generator_type, names, time_unit, time_bins, time_bin_names, tags)
};

/**
* Attribute Feature Generator: creating features from samples attributes
*/


class AttrFeatGenerator : public FeatureGenerator {
public:

	// Attribute to use
	string attribute;

	// Feature name (if empty - use attribute)
	string ftr_name;

	// Constructor/Destructor
	AttrFeatGenerator() { generator_type = FTR_GEN_ATTR; }
	~AttrFeatGenerator() {}

	// Naming 
	void set_names();

	/// The parsed fields from init command.
	/// @snippet FeatureGenerator.cpp TimeFeatGenerator::init
	int init(map<string, string>& mapper);

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<AttrFeatGenerator *>(generator)); }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// Serialization
	ADD_CLASS_NAME(AttrFeatGenerator);
	ADD_SERIALIZATION_FUNCS(generator_type, ftr_name, attribute, names);
};

enum class category_stat_test {
	chi_square = 1,
	mcnemar = 2
};

/**
* Creates multipal features based on categorical values and statistical dependency strength by Age,Gender groups
*/
class CategoryDependencyGenerator : public FeatureGenerator {
private:
	int bdate_sid;
	int gender_sid;
	map<int, vector<string>> categoryId_to_name; //for regex filter
	map<int, vector<int>> _member2Sets; //for hierarchy
	map<int, vector<int>> _set2Members; //for hierarchy
	unordered_map<int, vector<int>> _member2Sets_flat_cache; //for hierarchy cache in get_parents

	vector<string> top_codes;
	vector<vector<char>> luts;

	void get_parents(int codeGroup, vector<int> &parents, const boost::regex &reg_pat, const boost::regex & remove_reg_pat);

	void get_stats(const unordered_map<int, vector<vector<vector<int>>>> &categoryVal_to_stats,
		vector<int> &all_signal_values, vector<int> &signal_indexes, vector<double> &valCnts, vector<double> &posCnts,
		vector<double> &lift, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio, vector<int> &dof, const vector<vector<double>> &prior_per_bin) const;
public:
	string signalName; ///< the signal name
	int signalId;
	int time_channel;						///< n >= 0 : use time channel n , default: 0.
	int val_channel;						///< n >= 0 : use val channel n , default : 0.
	int win_from;///< time window for feature: win_from is the minimal time before from the prediction time
	int win_to;///< time window for feature: win_to is the maximal time before the prediction time			
	int time_unit_win;			///< the time unit in which the windows are given. Default: Undefined
	int min_age; ///< minimal age for testing statistical dependency
	int max_age; ///< maximal age for testing statistical dependency
	int age_bin; ///< age bin for testing statistical dependency
	string regex_filter; ///< regex filter for filtering categories in learn
	string remove_regex_filter; ///< remove regex filter for filtering categories in learn
	int min_code_cnt; ///< minimal number of occourences to consider signal
	float fdr; ///< the FDR value
	int take_top; ///< maximal number of features to create
	float lift_below; ///< filter lift to keep below it
	float lift_above; ///< filter lift to keep above it
	float filter_child_pval_diff; ///< below this threshold of pvalue diff change to remove child category (with AND condition on average lift change)
	float filter_child_lift_ratio; ///< below this threshold of lift change to remove child category
	float filter_child_count_ratio; ///< If child ratio count is too similar, small change from parent code - keep only paretn code
	float filter_child_removed_ratio; ///< If child removed ratio is beyond this and has other child taken - remove parent
	category_stat_test stat_metric; ///< statistical test
	float chi_square_at_least; ///< chi_square arg to test for at least that change in lift to measure bigger diffrence
	int minimal_chi_cnt; ///< chi_square arg to keep at least count to use row in calc
	int sort_by_chi = 0; ///< sort results by chi-square
	int max_depth; ///< maximal depth to go in heirarchy
	int max_parents; ///< controls maximum parents count
	bool use_fixed_lift; ///< If true will also sort be lifts below 1
	bool filter_hierarchy; /// Apply hierarchy filtering
	bool verbose; ///< in Learn will print selected features
	bool verbose_full; ///< If true will print a lot - table of all stats for each code
	string verbose_full_file; ///< output file for verbose_full debug in learn
	string feature_prefix; ///< additional prefix to add to name to describe the feature
	bool generate_with_counts; ///< If true will generate feature with counts not just as set

	void set_signal_ids(MedSignals& sigs);

	void init_tables(MedDictionarySections& dict);

	// Constructor/Destructor
	CategoryDependencyGenerator() : FeatureGenerator() { init_defaults(); };
	void init_defaults();

	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<CategoryDependencyGenerator *>(generator)); }

	/// The parsed fields from init command.
	/// @snippet CategoryDependencyGenerator.cpp CategoryDependencyGenerator::init
	int init(map<string, string>& mapper);

	void set_names();
	int filter_features(unordered_set<string>& validFeatures);

	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors);
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	int nfeatures();

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

	ADD_CLASS_NAME(CategoryDependencyGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, req_signals, top_codes, names, signalName, time_channel, val_channel, win_from, win_to, time_unit_win, feature_prefix, generate_with_counts, tags)
};

//=======================================
// Helpers
//=======================================

/// gets a [-_win_to, -_win_from] window in win time unit, and returns [_min_time, _max_time] window in signal time units relative to _win_time
/// boundOutcomeTime is used to future time windows when looking to the future to limit the time window till the outcomeTime
void get_window_in_sig_time(int _win_from, int _win_to, int _time_unit_win, int _time_unit_sig, int _win_time, int &_min_time, int &_max_time,
	bool boundOutcomeTime = false, int outcome_time = -1);

/// Conversion between time-range type and name
TimeRangeTypes time_range_name_to_type(const string& name);
string time_range_type_to_name(TimeRangeTypes type);

// update time window according to time-range signal
void get_updated_time_window(UniversalSigVec& time_range_usv, TimeRangeTypes type, int time_unit_range_sig, int time_unit_win, int time_unit_sig, int time,
	int win_from, int& updated_win_from, int win_to, int& updated_win_to, bool delta_win, int d_win_from, int& updated_d_win_from, int d_win_to, int& updated_d_win_to);
void get_updated_time_window(TimeRangeTypes type, int range_from, int range_to, int time, int _win_from, int _win_to, int& updated_win_from, int& updated_win_to);


//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(FeatureGenerator)
MEDSERIALIZE_SUPPORT(BasicFeatGenerator)
MEDSERIALIZE_SUPPORT(AgeGenerator)
MEDSERIALIZE_SUPPORT(GenderGenerator)
MEDSERIALIZE_SUPPORT(SingletonGenerator)
MEDSERIALIZE_SUPPORT(BinnedLmEstimatesParams)
MEDSERIALIZE_SUPPORT(BinnedLmEstimates)
MEDSERIALIZE_SUPPORT(RangeFeatGenerator)
MEDSERIALIZE_SUPPORT(ModelFeatGenerator)

#endif
