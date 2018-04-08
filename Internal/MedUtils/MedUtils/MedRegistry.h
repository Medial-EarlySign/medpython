#ifndef __MED_REGISTRY_H__
#define __MED_REGISTRY_H__
#include <vector>
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
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

	/// <summary>
	/// Writes the file to text file in tab delimeted format: PID, Start_Date, End_Date, min_allowed_date, max_allowed_date, Age, RegistryValue
	/// </summary>
	void write_text_file(const string &file_path);
	/// <summary>
	/// Reads the file in text format in tab delimeted
	/// </summary>
	void read_text_file(const string &file_path);

	/// <summary>
	/// Creates vector of registry using already initialized MedPidRepository with signals
	/// in parallel manner for each patient
	/// </summary>
	void create_registry(MedPidRepository &dataManager);

	/// <summary>
	/// returns the signal codes used to create the registry
	/// </summary>
	void get_registry_creation_codes(vector<int> &signal_codes);

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
		const string &debug_file = "", const unordered_set<float> &debug_vals = default_empty_set);

	ADD_SERIALIZATION_FUNCS(registry_records)
protected:
	vector<int> signalCodes; ///< The signals codes to fetch in create_registry. will be used in get_registry_records
private:
	virtual void get_registry_records(int pid, int bdate, vector<UniversalSigVec> &usv, vector<MedRegistryRecord> &results) { throw logic_error("Not Implemented"); };
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
		double calc_chi_square_dist(const map<float, vector<int>> &gender_sorted, int smooth_balls);

		/// \brief calcs chi square for full male, female and stores all the results stats and values in the vectors
		void calc_chi_scores(const map<float, map<float, vector<int>>> &male_stats,
			const map<float, map<float, vector<int>>> &female_stats,
			vector<float> &all_signal_values, vector<int> &signal_indexes,
			vector<double> &valCnts, vector<double> &posCnts,
			vector<double> &lift, vector<double> &scores,
			vector<double> &p_values, vector<double> &pos_ratio, int smooth_balls);

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
		int fetch_next_date(vector<UniversalSigVec> &patientFile, vector<int> &signalPointers);
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

	/// a function that retrive current outcome based on new time point
	virtual float get_outcome(UniversalSigVec &s, int current_i) = 0;

	static RegistrySignal *make_registry_signal(const string &type);
	static RegistrySignal *make_registry_signal(const string &type, const string &init_string);
};

/**
* A Class that condition a set of codes in dictionary
*/
class RegistrySignalSet : public RegistrySignal {
public:
	RegistrySignalSet(const string &sigName, int durr_time, int buffer_time, bool take_first,
		MedRepository &rep, const vector<string> &sets);
	RegistrySignalSet(const string &init_string, MedRepository &rep, const vector<string> &sets);
	float get_outcome(UniversalSigVec &s, int current_i);

	int init(map<string, string>& map);
private:
	vector<char> Flags;
};

/**
* A Class that condition a value range
*/
class RegistrySignalRange : public RegistrySignal {
public:
	RegistrySignalRange(const string &sigName, int durr_time, int buffer_time, bool take_first,
		float min_rnage, float max_range);
	float get_outcome(UniversalSigVec &s, int current_i);

	int init(map<string, string>& map);
private:
	float min_value;
	float max_value;
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

	/// <summary>
	/// The init function
	/// </summary>
	void init(MedRepository &rep, int start_dur, int end_durr, int max_repo,
		const vector<RegistrySignal *> signal_conditions, const string &skip_pid_file = "");
private:
	vector<bool> SkipPids; ///< black list of patients mask

	void get_registry_records(int pid, int bdate, vector<UniversalSigVec> &usv, vector<MedRegistryRecord> &results);
	bool init_called; ///< a flag to mark that init was called
};

MEDSERIALIZE_SUPPORT(MedRegistryRecord)
MEDSERIALIZE_SUPPORT(MedRegistry)

#endif
