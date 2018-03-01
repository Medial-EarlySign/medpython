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

namespace medial {
	namespace signal_hierarchy {
		///Getting signal with hirarchy options for each siganl
		void getRecords_Hir(int pid, vector<UniversalSigVec> &signals, MedDictionarySections &dict,
			const string &signalHirerchyType,
			vector<MedRegistryRecord> &res);

		///filtering hierarchy codes
		string filter_code_hierarchy(const vector<string> &vec, const string &signalHirerchyType);
		///getting parents in hierarchy codes
		vector<int> parents_code_hierarchy(MedDictionarySections &dict, const string &group, const string &signalHirerchyType);
		///getting sons in hierarchy codes
		vector<int> sons_code_hierarchy(MedDictionarySections &dict, const string &group, const string &signalHirerchyType);
		/// gets codes
		string get_readcode_code(MedDictionarySections &dict, int id, const string &signalHirerchyType);
	}
	namespace repository {
		///Helper function to calc diff between dates in years
		float DateDiff(int refDate, int dateSample);
		///Helper function to add days to date
		int DateAdd(int refDate, int daysAdd);
	}
	namespace contingency_tables {
		void calc_chi_scores(const map<float, map<float, vector<int>>> &male_stats,
			const map<float, map<float, vector<int>>> &female_stats,
			vector<float> &all_signal_values, vector<int> &signal_indexes,
			vector<double> &valCnts, vector<double> &posCnts,
			vector<double> &lift, vector<double> &scores,
			vector<double> &p_values, vector<double> &pos_ratio, int smooth_balls);

		void FilterRange(vector<int> &indexes, const vector<double> &vecCnts
			, double min_val, double max_val);

		double chisqr(int Dof, double Cv);
		void write_stats(const string &file_path, 
			const map<float, map<float, vector<int>>> &males_stats, const map<float, map<float, vector<int>>> &females_stats);
		void read_stats(const string &file_path, 
			map<float, map<float, vector<int>>> &males_stats, map<float, map<float, vector<int>>> &females_stats);

		void FilterFDR(vector<int> &indexes,
			const vector<double> &scores, const vector<double> &p_vals, const vector<double> &lift,
			double filter_pval);
	}
}

/**
* A Class which creates registry based on readcode lists.
*  Important: must be initialized by init_lists first
*/
class MedRegistryCodesList : public MedRegistry {
public:
	vector<char> RCFlags; ///< Readcode flags
	vector<bool> SkipPids; ///< black list of patients

	int duration_flag; ///< the duration for each positive to merge time ranges
	int buffer_duration; ///< a buffer duration between positive to negative
	bool take_only_first; ///< if True will take only first occournce
	int max_repo_date; ///< the maximal date for the repository

	/// <summary>
	/// The init function
	/// </summary>
	void init_lists(MedRepository &rep, int dur_flag, int buffer_dur, bool takeOnlyFirst,
		int max_repo, const vector<string> *rc_sets = NULL, const string &skip_pid_file = "");
private:
	void get_registry_records(int pid, int bdate, vector<UniversalSigVec> &usv, vector<MedRegistryRecord> &results);
	bool init_lists_called;
};

MEDSERIALIZE_SUPPORT(MedRegistryRecord)
MEDSERIALIZE_SUPPORT(MedRegistry)

#endif
