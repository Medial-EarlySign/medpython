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

/**
* A class that holds all registry records on all patients
*/
class MedRegistry : public SerializableObject
{
public:
	vector<MedRegistryRecord> registry_records; ///< the registry record vector
	vector<int> signalCodes; ///< The signalsCodes to fetch. make private

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
	/// </summary>
	void calc_dependency_mem(const string &repository_path, int signalCode,
		const string &signalHirerchyType, int ageBinValue, int time_window_from, int time_window_to,
		MedSamplingStrategy &sampler,
		map<float, map<float, vector<int>>> &maleSignalToStats,
		map<float, map<float, vector<int>>> &femaleSignalToStats);

	ADD_SERIALIZATION_FUNCS(registry_records)
private:
	virtual void get_registry_records(int pid, int bdate, vector<UniversalSigVec> &usv, vector<MedRegistryRecord> &results) { throw logic_error("Not Implemented"); };
};

///Helper function to calc diff between dates in years
float DateDiff(int refDate, int dateSample);
///Helper function to add days to date
int DateAdd(int refDate, int daysAdd);

/**
* A Class which creates registry based on readcode lists
*/
class MedRegistryCodesList : public MedRegistry {
public:
	vector<char> RCFlags; ///< Readcode flags
	vector<bool> SkipPids; ///< black list of patients

	int duration_flag; ///< the duration for each positive to merge time ranges
	int buffer_duration; ///< a buffer duration between positive to negative
	bool take_only_first; ///< if True will take only first occournce

	/// <summary>
	/// The init function
	/// </summary>
	void init_lists(MedRepository &rep, int dur_flag, int buffer_dur, bool takeOnlyFirst,
		const vector<string> *rc_sets = NULL, const string &skip_pid_file = "");
private:
	void get_registry_records(int pid, int bdate, vector<UniversalSigVec> &usv, vector<MedRegistryRecord> &results);
};

MEDSERIALIZE_SUPPORT(MedRegistryRecord)
MEDSERIALIZE_SUPPORT(MedRegistry)

#endif
