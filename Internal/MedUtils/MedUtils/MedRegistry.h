#ifndef __MED_REGISTRY_H__
#define __MED_REGISTRY_H__
#include <vector>
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include <MedProcessTools/MedProcessTools/RepProcess.h>
#include "MedSamplingStrategy.h"

using namespace std;

class MedSamplingStrategy;

/**
* A class which represnt a registry record of patient in time range from start_date to end_date
* It has min_allowed and max_allowed dates for sampling and it has optional mark for
* current sample age
*/
class MedRegistryRecord : public SerializableObject
{
public:
	int pid; ///< patient ID
	//defines the registry value apply date range
	int start_date; ///< the start_date range for the record
	int end_date; ///< the end_date range for the record
	//defines the allowed sampling date range for the record
	int min_allowed_date; ///< min time for sampling
	int max_allowed_date; ///< max time for sampling

	int age; ///< optional mark for record age. the age in start_date
	float registry_value; ///< the registry value/state

	ADD_SERIALIZATION_FUNCS(pid, start_date, end_date, min_allowed_date, max_allowed_date, age, registry_value)
};

static unordered_set<float> default_empty_set;

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
	/// returns the signal codes used to create the registry
	/// </summary>
	void get_registry_creation_codes(vector<int> &signal_codes) const;

	/// <summary>
	/// calculates table statitics for interrsecting with registry of signal
	/// @param repository_path the repsitory path
	/// @param signalCode the signal code to calculate the stats of the registry with
	/// @param signalHirerchyType the Hirerchy type: [None, RC, ATC, BNF] for signals with hirerchy
	/// @param ageBinValue the age bin size
	/// @param time_window_from the minimal time before the event (registry start_time). may be negative to start after event
	/// @param time_window_to the maximal time before the event (registry start_time). may be negative to end after event start
	/// @param sampler the sampler for how to calc the non appearance of the signal with the registry. 
	/// You may use MedSamplingAge for example to sample each patient once in each age in the registry dates
	/// @param debug_file If provided the output path to write detailed results of the intersection of registry and signal
	/// @param debug_vals If not empty and has debug_file. will write the intersection(by the time window) of the registry with those values
	/// </summary>
	/// <returns>
	/// @param maleSignalToStats The stats for males. the first key in the dictionary is the signal_value.
	/// the second key is age_bin and the vector is always of size 4: [signal_not_appear&registry_is_false, signal_not_appear&registry_is_true, signal_appears&registry_is_false, signal_appears&registry_is_true]
	/// @param femaleSignalToStats The stats for the females. same format as males
	/// </returns>
	void calc_signal_stats(const string &repository_path, int signalCode,
		const string &signalHirerchyType, int ageBinValue, int time_window_from, int time_window_to,
		MedSamplingStrategy &sampler,
		map<float, map<float, vector<int>>> &maleSignalToStats,
		map<float, map<float, vector<int>>> &femaleSignalToStats,
		const string &debug_file = "", const unordered_set<float> &debug_vals = default_empty_set) const;

	/// <summary>
	/// returns all patients ids from registry - unique patient ids
	/// @param pids the unique patient ids result vector
	/// </summary>
	void get_pids(vector<int> &pids) const;

	/// <summary>
	/// calculate incidence and writes the result into file with old and new format
	/// @param file_path the output file path to write the results
	/// @param rep_path the repository path to calculate the incidence
	/// @param age_bin the age_bin for binning age groups for the incidence
	/// @param min_age the minimal age fro the incidence
	/// @param max_age the maximal age fro the incidence
	/// @param time_period the time period in days for the incidence
	/// @param use_kaplan_meir if True will calc using kaplan meier survivol rates
	/// @param sampler_name the sampler name for calculating incidence
	/// @param sampler_args the sampler args for calculating incidence - may control trail years for example
	/// </summary>
	void create_incidence_file(const string &file_path, const string &rep_path, int age_bin, int min_age,
		int max_age, int time_period = 365, bool use_kaplan_meir = false, const string &sampler_name = "yearly",
		const string &sampler_args = "conflict_method=max;use_allowed=1;day_jump=365;allowed_time_from=0;"
		"allowed_time_to=365;start_year=2007;end_year=2012") const;

	/// creates registry type and initialize it if init_str is not empty
	/// Use "binary" for MedRegistryCodesList and "categories" for MedRegistryCategories.
	/// @snippet MedRegistry.cpp MedRegistry::make_registry
	static MedRegistry *make_registry(const string &registry_type, const string &init_str = "");

	/// Default Ctor
	MedRegistry() {
		need_bdate = false;
		time_unit = global_default_time_unit;
	}

	/// A function to clear creation variables that are on memory if needed
	virtual void clear_create_variables() {};

	ADD_SERIALIZATION_FUNCS(registry_records)
protected:
	vector<int> signalCodes; ///< The signals codes to fetch in create_registry. will be used in get_registry_records
	bool need_bdate; ///< If true Bdate is also used in registry creation
private:
	virtual void get_registry_records(int pid, int bdate, vector<UniversalSigVec_mem> &usv, vector<MedRegistryRecord> &results) { throw logic_error("Not Implemented"); };
};

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	* \brief signal_hierarchy namespace
	*/
	namespace signal_hierarchy {
		/// \brief Getting signal with hirarchy options for each siganl
		void getRecords_Hir(int pid, vector<UniversalSigVec> &signals, MedDictionarySections &dict,
			const string &signalHirerchyType,
			vector<MedRegistryRecord> &res);
	}
	/*!
	* \brief contingency_tables namespace
	*/
	namespace contingency_tables {
		/// \brief calc chi square distance for all groups with 4 value vector
		double calc_chi_square_dist(const map<float, vector<int>> &gender_sorted, int smooth_balls = 0,
			float allowed_error = 0, int minimal_balls = 0);

		/// \brief calc mcnemar square distance for all groups with 4 value vector
		double calc_mcnemar_square_dist(const map<float, vector<int>> &gender_sorted);

		/// \brief calcs chi square for full male, female and stores all the results stats and values in the vectors
		void calc_chi_scores(const map<float, map<float, vector<int>>> &male_stats,
			const map<float, map<float, vector<int>>> &female_stats,
			vector<float> &all_signal_values, vector<int> &signal_indexes,
			vector<double> &valCnts, vector<double> &posCnts,
			vector<double> &lift, vector<double> &scores,
			vector<double> &p_values, vector<double> &pos_ratio, int smooth_balls = 0, float allowed_error = 0,
			int minimal_balls = 0);

		/// \brief calcs mcnemar test square for full male, female and stores all the results stats and values in the vectors
		void calc_mcnemar_scores(const map<float, map<float, vector<int>>> &male_stats,
			const map<float, map<float, vector<int>>> &female_stats,
			vector<float> &all_signal_values, vector<int> &signal_indexes,
			vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
			, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio);

		/// \brief filter by range
		void FilterRange(vector<int> &indexes, const vector<double> &vecCnts
			, double min_val, double max_val);

		/// \brief calc chi square probabilty from distance, DOF
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
	namespace io {
		void read_codes_file(const string &file_path, vector<string> &tokens);
	}
	namespace repository {
		/// \brief fetches the next date from all signals in patientFile by date order.
		/// the signalPointers is array of indexes of each signal. it also advances the right index
		/// returns the signal with the minimal date - "the next date"
		template<class T> int fetch_next_date(vector<T> &patientFile, vector<int> &signalPointers);
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
	float outcome_value; ///< the outcome value when condition holds

	/// a function that retrive current outcome based on new time point
	virtual bool get_outcome(UniversalSigVec &s, int current_i, float &result) = 0;

	/// creates Registry rule. can have "set" for RegistrySignalSet and "range" for RegistrySignalRange.
	/// /// @snippet MedRegistry.cpp RegistrySignal::make_registry_signal
	static RegistrySignal *make_registry_signal(const string &type, MedRepository &rep);
	/// creates Registry rule and uses init_string to initialize the type
	static RegistrySignal *make_registry_signal(const string &type, MedRepository &rep, const string &init_string);

	/// <summary>
	/// parsing of registry signal rules - each line is new signal rule in this format:\n
	/// Each line is TAB seperated by RegistrySignal type and RegistrySignal init string calling 
	/// RegistrySignal::make_registry_signal 
	/// </summary>
	static void parse_registry_rules(const string &reg_cfg, MedRepository &rep,
		vector<RegistrySignal *> &result);
};

/**
* A Class that condition a set of codes in dictionary.
* use "set" keyword to refernce this class
*/
class RegistrySignalSet : public RegistrySignal {
public:
	RegistrySignalSet(const string &sigName, int durr_time, int buffer_time, bool take_first,
		MedRepository &rep, const vector<string> &sets, float outcome_val = 1);
	RegistrySignalSet(const string &init_string, MedRepository &rep, const vector<string> &sets, float outcome_val = 1);
	bool get_outcome(UniversalSigVec &s, int current_i, float &result);

	/// The parsed fields from init command.\n
	/// @snippet MedRegistry.cpp RegistrySignalSet::init
	int init(map<string, string>& map);

	/// Checks if has flags inside or it's empty one
	bool is_empty() { return Flags.empty(); }
private:
	vector<char> Flags;
	MedRepository *repo;
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
		float min_range, float max_range, float outcome_val = 1);
	bool get_outcome(UniversalSigVec &s, int current_i, float &result);

	/// The parsed fields from init command.\n
	/// @snippet MedRegistry.cpp RegistrySignalRange::init
	int init(map<string, string>& map);
private:

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

	vector<RegistrySignal *> signal_filters; ///< the signal filters

	MedRegistryCodesList() {
		init_called = false;
		start_buffer_duration = 0;
		end_buffer_duration = 0;
		max_repo_date = 0;
		need_bdate = false;
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
private:
	void get_registry_records(int pid, int bdate, vector<UniversalSigVec_mem> &usv, vector<MedRegistryRecord> &results);
};

MEDSERIALIZE_SUPPORT(MedRegistryRecord)
MEDSERIALIZE_SUPPORT(MedRegistry)

#endif
