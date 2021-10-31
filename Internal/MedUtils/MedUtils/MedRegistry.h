/// @file
/// registry methods over MedRegistry Object
#ifndef __MED_REGISTRY_H__
#define __MED_REGISTRY_H__
#include <vector>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/RepProcess.h>
#include "MedSamplingStrategy.h"
#include "MedLabels.h"
#include "MedEnums.h"
#include "MedRegistryRecord.h"

using namespace std;

/**
* A class that holds all registry records on all patients.\n
* It has several ways to be initialized:\n
* 1. by reading from disk - binary format or text format\n
* 2. by creating registry using create_registry method. need to implement
* get_registry_records to handle single patient records.\n
* \n
* the class have also the ability to create contingency table with other signal:\n
* for each Gender,Age_bin - the 4 stats number of the registry with the appearances or not appearances of
* the signal value
*/
class MedRegistry : public SerializableObject
{
public:
	vector<MedRegistryRecord> registry_records; ///< the registry records vector
	int time_unit; ///< The time unit

	/// <summary>
	/// Writes the file to text file in tab delimeted format: PID, Start_Date, End_Date, min_allowed_date, max_allowed_date, Age, RegistryValue
	/// </summary>
	void write_text_file(const string &file_path) const;
	/// <summary>
	/// Reads the file in text format in tab delimeted
	/// </summary>
	void read_text_file(const string &file_path);

	/// <summary>
	/// Creates vector of registry using already initialized MedPidRepository with signals
	/// in parallel manner for each patient
	/// </summary>
	void create_registry(MedPidRepository &dataManager, medial::repository::fix_method method = medial::repository::fix_method::none, vector<RepProcessor *> *rep_processors = NULL);

	/// <summary>
	/// returns the signal codes used to create the registry - to load from rep
	/// </summary>
	void get_registry_creation_codes(vector<string> &signal_codes) const;
	/// <summary>
	/// returns the signal codes used to create the registry - to use in create_registry command as input
	/// </summary>
	void get_registry_use_codes(vector<string> &signal_codes) const;


	/// <summary>
	/// returns all patients ids from registry - unique patient ids
	/// @param pids the unique patient ids result vector
	/// </summary>
	void get_pids(vector<int> &pids) const;

	/// <summary>
	/// Merges registry record with same registry values (and continues in time)
	/// </summary>
	void merge_records();

	void *new_polymorphic(string dname);

	/// creates registry type and initialize it if init_str is not empty
	/// Use "binary" for MedRegistryCodesList and "categories" for MedRegistryCategories.
	/// @snippet MedRegistry.cpp MedRegistry::make_registry
	static MedRegistry *make_registry(const string &registry_type, const string &init_str = "");

	/// creates registry type and initialize it if init_str is not empty
	/// Use "binary" for MedRegistryCodesList and "categories" for MedRegistryCategories.
	static MedRegistry *make_registry(const string &registry_type, MedRepository &rep, const string &init_str = "");

	/// <summary>
	/// Creates vector of registry records - handles everything for you
	/// in parallel manner for each patient - uses create_registry
	/// </summary>
	static MedRegistry *create_registry_full(const string &registry_type, const string &init_str,
		const string &repository_path, MedModel &model_with_rep_processor, medial::repository::fix_method method = medial::repository::fix_method::none);

	/// Default Ctor
	MedRegistry() {
		need_bdate = false;
		time_unit = global_default_time_unit;
	}

	virtual ~MedRegistry() {};

	/// A function to clear creation variables that are on memory if needed
	virtual void clear_create_variables() {};

	virtual bool get_pid_records(PidRec &rec, int bDateCode, const vector<int> &used_sigs, vector<MedRegistryRecord> &results);

	/// Sets Repository object to initialize all registry object, if not given will try to use repository path
	/// to read and initialize repository
	void set_rep_for_init(MedRepository &rep) { rep_for_init = &rep; }

	ADD_CLASS_NAME(MedRegistry)
		ADD_SERIALIZATION_FUNCS(time_unit, registry_records)
protected:
	vector<string> signalCodes_names; ///< the signals codes by name
	bool need_bdate; ///< If true Bdate is also used in registry creation
	medial::repository::fix_method resolve_conlicts = medial::repository::fix_method::none; ///< resolve conflicts in registry method
	MedRepository *rep_for_init = NULL; ///< repository pointer to init dicts
private:
	virtual void get_registry_records(int pid, int bdate, vector<UniversalSigVec_mem> &usv, vector<MedRegistryRecord> &results) { throw logic_error("Not Implemented"); };
};

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	* \brief contingency_tables namespace
	*/
	namespace contingency_tables {
		/// \brief calc chi square distance for all groups with 4 value vector
		double calc_chi_square_dist(const map<float, vector<int>> &gender_sorted, int smooth_balls = 0,
			float allowed_error = 0, int minimal_balls = 0);

		/// \brief calc mcnemar square distance for all groups with 4 value vector
		double calc_mcnemar_square_dist(const map<float, vector<int>> &gender_sorted);

		/// \brief calc cmh  https://en.wikipedia.org/wiki/Cochran%E2%80%93Mantel%E2%80%93Haenszel_statistics
		double calc_cmh_square_dist(const map<float, vector<int>> *gender_sorted, const map<float, vector<int>> *gender_sorted2, bool &ok);

		/// \brief calcs chi square for full male, female and stores all the results stats and values in the vectors
		void calc_chi_scores(const map<float, map<float, vector<int>>> &male_stats,
			const map<float, map<float, vector<int>>> &female_stats,
			vector<int> &all_signal_values, vector<int> &signal_indexes,
			vector<double> &valCnts, vector<double> &posCnts,
			vector<double> &lift, vector<double> &scores,
			vector<double> &p_values, vector<double> &pos_ratio, int smooth_balls = 0, float allowed_error = 0,
			int minimal_balls = 0);

		/// \brief calcs mcnemar test square for full male, female and stores all the results stats and values in the vectors
		void calc_mcnemar_scores(const map<float, map<float, vector<int>>> &male_stats,
			const map<float, map<float, vector<int>>> &female_stats,
			vector<int> &all_signal_values, vector<int> &signal_indexes,
			vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
			, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio);

		/// \brief calc cmh  https://en.wikipedia.org/wiki/Cochran%E2%80%93Mantel%E2%80%93Haenszel_statistics
		void calc_cmh_scores(const map<float, map<float, vector<int>>> &male_stats,
			const map<float, map<float, vector<int>>> &female_stats,
			vector<int> &all_signal_values, vector<int> &signal_indexes,
			vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
			, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio);

		/// \brief filter by range
		void FilterRange(vector<int> &indexes, const vector<double> &vecCnts
			, double min_val, double max_val);

		/// \brief filter category code by hirarchy of categories. similar child code to parent will be removed.
		/// Similar is or:
		/// 1. small  change in count - below count_similarity cahnge in ratio
		/// 2. small change in P-Value And Average Lift - the change in P-Value is controlled with pValue_diff, and in lift (change factor child lift to paretn lift) is by lift_th
		void filterHirarchy(const map<int, vector<int>> &member2Sets, const map<int, vector<int>> &set2Members,
			vector<int> &indexes, const vector<int> &signal_values, const vector<double> &pVals,
			const vector<double> &valCnts, const vector<double> &lifts, const unordered_map<int, double> &code_unfiltered_cnts,
			float pValue_diff, float lift_th, float count_similarity, float child_fitlered_ratio, const map<int, vector<string>> *categoryId_to_name = NULL);

		/// \brief calc chi square probability from distance, DOF
		double chisqr(int Dof, double Cv);
		/// \brief serialize male,female stats
		void write_stats(const string &file_path,
			const map<float, map<float, vector<int>>> &males_stats, const map<float, map<float, vector<int>>> &females_stats);
		/// \brief deserialize male,female stats
		void read_stats(const string &file_path,
			map<float, map<float, vector<int>>> &males_stats, map<float, map<float, vector<int>>> &females_stats);
		/// \brief filter by FDR
		void FilterFDR(vector<int> &indexes,
			const vector<double> &scores, const vector<double> &p_vals, const vector<double> &lift,
			double filter_pval);
	}
	/*!
	*  \brief print namespace
	*/
	namespace print {
		/// \brief printing registry stats for labels inside of it.
		void print_reg_stats(const vector<MedRegistryRecord> &regRecords, const string &log_file = "");
	}
	/*!
	*  \brief registry namespace
	*/
	namespace registry {
		/// \brief completes control period for registry giving active period for patient. active_periods_registry - is time ranges for each patient (not looking on registry_value - like in censor registry)
		void complete_active_period_as_controls(vector<MedRegistryRecord> &registry,
			const vector<MedRegistryRecord> &active_periods_registry, bool unite_full_controls = true);
	}
}

/**
* A abstract class that represents a signal used to create registry and it's
* filter conditions to change outcome values based on current time point
* usefull if the registry depends seperatly by each signal / only one signal
*/
class RegistrySignal : public SerializableObject {
public:
	string signalName; ///< the signal name
	int duration_flag; ///< the duration for each positive to merge time ranges
	int buffer_duration; ///< a buffer duration between positive to negative
	bool take_only_first; ///< if True will take only first occournce
	int channel; ///< the channel number the rule operates on
	float outcome_value; ///< the outcome value when condition holds

	/// Default init ctor for object, that won't contain garbage when not initialized specifically
	RegistrySignal() {
		signalName = "";
		duration_flag = 0;
		buffer_duration = 0;
		take_only_first = false;
		channel = 0;
		outcome_value = 1;
	}

	/// a function that retrive current outcome based on new time point
	virtual bool get_outcome(const UniversalSigVec &s, int current_i, float &result) = 0;

	/// creates Registry rule. can have "set" for RegistrySignalSet and "range" for RegistrySignalRange.
	/// /// @snippet MedRegistry.cpp RegistrySignal::make_registry_signal
	static RegistrySignal *make_registry_signal(const string &type, MedRepository &rep, const string &path_to_cfg_file);
	/// creates Registry rule and uses init_string to initialize the type
	static RegistrySignal *make_registry_signal(const string &type, MedRepository &rep, const string &init_string, const string &path_to_cfg_file);

	/// <summary>
	/// parsing of registry signal rules - each line is new signal rule in this format:\n
	/// Each line is TAB seperated by RegistrySignal type and RegistrySignal init string calling 
	/// RegistrySignal::make_registry_signal 
	/// </summary>
	static void parse_registry_rules(const string &reg_cfg, MedRepository &rep,
		vector<RegistrySignal *> &result);

	/// The parsed fields from init command.\n
	/// @snippet MedRegistry.cpp RegistrySignal::init
	int init(map<string, string>& mapper);

	/// Each specific init function for pther arguments - called from init
	virtual void _init(const map<string, string>& mapper) {};

	virtual ~RegistrySignal() {};

	ADD_CLASS_NAME(RegistrySignal)
		ADD_SERIALIZATION_FUNCS(signalName, duration_flag, buffer_duration, take_only_first, channel, outcome_value)

		void *new_polymorphic(string dname);
};

/**
* A Class that condition a set of codes in dictionary.
* use "set" keyword to refernce this class
*/
class RegistrySignalSet : public RegistrySignal {
public:
	vector<string> sets;

	RegistrySignalSet(const string &sigName, int durr_time, int buffer_time, bool take_first,
		MedRepository &rep, const vector<string> &sets, const string &path_to_cfg_file, float outcome_val = 1, int chan = 0);
	RegistrySignalSet(const string &init_string, MedRepository &rep, const vector<string> &sets, const string &path_to_cfg_file, float outcome_val = 1);
	bool get_outcome(const UniversalSigVec &s, int current_i, float &result);

	/// The parsed fields from init command.\n
	/// @snippet MedRegistry.cpp RegistrySignalSet::_init
	void _init(const map<string, string>& mapper);

	/// Checks if has flags inside or it's empty one
	bool is_empty() { return Flags.empty(); }
private:
	vector<char> Flags;
	MedRepository *repo;
	string base_cfg_path;
};

/**
* A Class that condition a value range.
* use "range" keyword to refernce this class
*/
class RegistrySignalRange : public RegistrySignal {
public:
	float min_value; ///< the minimal value to turn control into case. greater than or equal
	float max_value; ///< the maximal value to turn control into case. smaller than or equal

	RegistrySignalRange(const string &sigName, int durr_time, int buffer_time, bool take_first,
		float min_range, float max_range, float outcome_val = 1, int chan = 0);
	bool get_outcome(const UniversalSigVec &s, int current_i, float &result);

	/// The parsed fields from init command.\n
	/// @snippet MedRegistry.cpp RegistrySignalRange::_init
	void _init(const map<string, string>& mapper);
private:

};

/**
* A Class that conditions nothing, just exising of the signal. usefull for DEATH signal
* Can have only time channel.
* use "aby" keyword to refernce this class
*/
class RegistrySignalAny : public RegistrySignal {

	bool get_outcome(const UniversalSigVec &s, int current_i, float &result);

	ADD_CLASS_NAME(RegistrySignalAny)
		ADD_SERIALIZATION_FUNCS(signalName, duration_flag, buffer_duration, take_only_first, channel, outcome_value)
};

/**
* A Registry operator to handle drugs with condition on drug type and dosage range
*/
class RegistrySignalDrug : public RegistrySignal {
public:
	vector<string> sets;
	RegistrySignalDrug(MedRepository &rep, const string &path_to_cfg_file);
	/// The parsed fields from init command.\n
	/// @snippet MedRegistry.cpp RegistrySignalDrug::_init
	void _init(const map<string, string>& mapper);

	/// Checks if has flags inside or it's empty one
	bool is_empty() { return Flags.empty(); }

	bool get_outcome(const UniversalSigVec &s, int current_i, float &result);
private:
	vector<char> Flags; ///< first if exists
	vector<pair<float, float>> Flags_range; ///< range for dosage
	MedRepository *repo;
	string base_path;
};

/**
* A Registry Signal class wrapper for AND condition on multiple Registry signal channels.
* it works only on same signal on the same time point
*/
class RegistrySignalAnd : public RegistrySignal {
public:
	vector<RegistrySignal *> conditions; ///< the list of conditions to calc AND on them

	RegistrySignalAnd(MedRepository &rep);
	/// The parsed fields from init command.\n
	/// @snippet MedRegistry.cpp RegistrySignalAnd::_init
	void _init(const map<string, string>& mapper);

	bool get_outcome(const UniversalSigVec &s, int current_i, float &result);

	~RegistrySignalAnd();
private:
	MedRepository * repo;
};

/**
* A Class which creates registry based on readcode lists.
*  Important: must be initialized by init_lists first
*/
class MedRegistryCodesList : public MedRegistry {
public:
	int start_buffer_duration; ///< the duration buffer form start
	int end_buffer_duration; ///< the duration buffer from last date
	int max_repo_date; ///< the maximal date for the repository
	bool allow_prediciton_in_case; ///< If True will allow to give prediciton after\in case time range
	bool seperate_cases; ///< If true will seperate each "case" time zone

	vector<RegistrySignal *> signal_filters; ///< the signal filters

	MedRegistryCodesList() {
		init_called = false;
		start_buffer_duration = 0;
		end_buffer_duration = 0;
		max_repo_date = 0;
		need_bdate = false;
		allow_prediciton_in_case = false;
		seperate_cases = false;
	}

	~MedRegistryCodesList() {
		clear_create_variables();
	}

	/// <summary>
	/// The init function in code API
	/// @param rep initialized repository with MedDictionry for initialization
	/// @param start_dur a minimal time for patient to enter registry from first signal after birth
	/// @param end_durr a minimal time for patient to leave registry from last signal
	/// @param max_repo the last date in the repositry - censor after this date
	/// @param signal_conditions vector of rules to calc when we turn patient into case
	/// @param skip_pid_file a file with blacklist of patient ids to skip
	/// @param pid_to_censor_dates an object to map between each patient and censor date for him
	/// </summary>
	void init(MedRepository &rep, int start_dur, int end_durr, int max_repo,
		const vector<RegistrySignal *> signal_conditions, const string &skip_pid_file = "",
		const unordered_map<int, int> *pid_to_censor_dates = NULL);

	/// <summary>
	/// the initializtion params. it has also "config_signals_rules", "pid_to_censor_dates", "rep" file paths.
	/// @param rep the repository path
	/// @param pid_to_censor_dates file path to pid censors. each line is pid TAB censor_date
	/// @param config_signals_rules file path to RegistrySignal rules. parsing is done with 
	/// MedRegistryCodesList::parse_registry_rules \n
	/// The parsed fields from init command.
	/// @snippet MedRegistry.cpp MedRegistryCodesList::init
	/// </summary>
	int init(map<string, string>& map);

	///clears the signal_filters
	void clear_create_variables();

	ADD_CLASS_NAME(MedRegistryCodesList)
		ADD_SERIALIZATION_FUNCS(time_unit, registry_records,
			start_buffer_duration, end_buffer_duration, max_repo_date, allow_prediciton_in_case,
			seperate_cases)
private:
	vector<bool> SkipPids; ///< black list of patients mask
	unordered_map<int, int> pid_to_max_allowed; ///< max date allowed to each pid constrain

	void get_registry_records(int pid, int bdate, vector<UniversalSigVec_mem> &usv, vector<MedRegistryRecord> &results);
	bool init_called; ///< a flag to mark that init was called
};

/**
* A Regsitry creator to create categoriezed outcome by signal rules.
* Esch signal is condition independence in the rest of the signals.
*/
class MedRegistryCategories : public MedRegistry {
public:
	int start_buffer_duration; ///< the duration buffer form start
	int end_buffer_duration; ///< the duration buffer from last date
	int max_repo_date; ///< the maximal date for the repository

	vector<vector<RegistrySignal *>> signals_rules; ///< the signal rules vectors, first index is signal id, second is list of rules

	/// Initialize class parameters - it also needs repository_path parameter which called "rep".
	/// @snippet MedRegistry.cpp MedRegistryCategories::init
	int init(map<string, string>& map);

	MedRegistryCategories() {
		start_buffer_duration = 0;
		end_buffer_duration = 0;
		max_repo_date = 0;
		need_bdate = false;
	}

	///clears the signals_rules
	void clear_create_variables();

	~MedRegistryCategories() {
		clear_create_variables();
	}

	ADD_CLASS_NAME(MedRegistryCategories)
		ADD_SERIALIZATION_FUNCS(time_unit, registry_records, start_buffer_duration, end_buffer_duration,
			max_repo_date)
private:
	unordered_map<int, int> pid_to_max_allowed; ///< max date allowed to each pid constrain

	void get_registry_records(int pid, int bdate, vector<UniversalSigVec_mem> &usv, vector<MedRegistryRecord> &results);
};

/**
* Keep Alive registry - use for censoring "dead" times
*/
class MedRegistryKeepAlive : public MedRegistry {
public:
	int duration; ///< the duration buffer form start
	int max_repo_date; ///< the maximal date for the repository
	int start_buffer_duration; ///< the buffer duration from first signal
	int secondry_start_buffer_duration; ///< the buffer duration for new region (after not active) - can be negative to look backward
	int end_buffer_duration; ///< the buffer duration from last signal
	vector<string> signal_list; ///< list of signals to fetch for keep alive time ranges

	MedRegistryKeepAlive() {
		duration = 0;
		max_repo_date = 0;
		need_bdate = false;
		start_buffer_duration = 0;
		end_buffer_duration = 0;
		secondry_start_buffer_duration = 0;
	}

	/// <summary>
	/// @snippet MedRegistry.cpp MedRegistryKeepAlive::init
	/// </summary>
	int init(map<string, string>& map);

	ADD_CLASS_NAME(MedRegistryKeepAlive)
		ADD_SERIALIZATION_FUNCS(time_unit, registry_records, duration, max_repo_date, start_buffer_duration,
			secondry_start_buffer_duration, end_buffer_duration, signal_list)
private:
	unordered_map<int, int> pid_to_max_allowed; ///< max date allowed to each pid constrain

	void get_registry_records(int pid, int bdate, vector<UniversalSigVec_mem> &usv, vector<MedRegistryRecord> &results);
};

MEDSERIALIZE_SUPPORT(MedRegistry)
MEDSERIALIZE_SUPPORT(MedRegistryCodesList)
MEDSERIALIZE_SUPPORT(MedRegistryCategories)
MEDSERIALIZE_SUPPORT(MedRegistryKeepAlive)

#endif
