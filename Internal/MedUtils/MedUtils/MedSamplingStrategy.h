/// @file
/// Sampling methods over MedRegistry Object
#ifndef __MED_SAMPLING_STRATEGY_H__
#define __MED_SAMPLING_STRATEGY_H__

#include <vector>
#include <random>
#include <MedUtils/MedUtils/MedRegistry.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>

class MedRegistryRecord;

using namespace std;

/// @enum
/// Sampling options
enum SamplingMode {
	Case_Control = 0,///< "case_control" - cases need to pass start_time of registry. controls need to be within start_time to before end_time
	Pass = 1, ///< "pass" - cases need to pass start_time of registry. controls need to pass start_time. controls are like cases in case_control
	Within = 2, ///< "within" - cases need to be within start_time and end_time. controls need to be within start_time and end_time. cases are like control in case_control
};
extern vector<string> SamplingMode_to_name;
SamplingMode SamplingMode_name_to_type(const string& SamplingMode_name);

/// @enum
/// Conflicting options
enum ConflictMode {
	All = 0,///< "all" - take all
	Drop = 1, ///< "drop" - drop when conflcit
	Max = 2, ///< "max" - take max on conflict 
};
extern vector<string> ConflictMode_to_name;
ConflictMode ConflictMode_name_to_type(const string& ConflictMode_name);

/**
* An abstract class with sampling methods over registry records to convert to MedSamples
*/
class MedSamplingStrategy : public SerializableObject {
public:
	/// The sampling method to be implemented
	virtual void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples) = 0;

	/// @snippet MedSamplingStrategy.cpp MedSamplingStrategy::make_sampler
	static MedSamplingStrategy *make_sampler(const string &sampler_name);
	/// @snippet MedSamplingStrategy.cpp MedSamplingStrategy::make_sampler
	static MedSamplingStrategy *make_sampler(const string &sampler_name, const string &init_params);
};

/**
* A Class which samples records on registry for certain time window.
* For Each registry record it samples randomly time point which falls withing the time window
* from min_allowed to max_allowed or till max_allowed backward using the time window params
*/
class MedSamplingTimeWindow : public MedSamplingStrategy {
public:
	bool take_max; ///< If true will random sample between all time range of min_allowed to max_allowed
	int maximal_time_case; ///< If take_max is false. will sample till maximal_time_case duration from max_allowed for cases. In days
	int minimal_time_case; ///< Will treat at least minimal_time_case duration from max_allowed for cases. In days
	int maximal_time_control; ///< If take_max is false. Will sample till maximla_time_control duration from max_allowed for controls. In days
	int minimal_time_control; ///< Will treat at least minimal_time_control duration from max_allowed for controls. In days
	int sample_count; ///< how many samples to take in each time window

	///sample random using Environment variable. params: [Random_Duration, Back_Time_Window_Years, Jump_Time_Period_Years]
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples);

	int init(map<string, string>& map);

	MedSamplingTimeWindow() {
		gen = mt19937(rd()); sample_count = 1; take_max = false;
		maximal_time_case = 0; minimal_time_case = 0; minimal_time_control = 0; maximal_time_control = 0;
	}
private:
	random_device rd;
	mt19937 gen;
};

/**
* A Class which samples by year from year to year by jump and find match in registry.
* suitble for incidence calculation
*/
class MedSamplingYearly : public MedSamplingStrategy {
public:
	int start_year; ///< The start year to sample from
	int end_year; ///< The end year to sample from
	int prediction_month_day; ///< the prediciton month_day in each year
	int back_random_duration; ///< Random duration backward from prediciton month_day. to cancel use 0
	int day_jump; ///< the years bin, how many years to jump backward from each prediciton date
	SamplingMode mode; ///< sampling mode
	int time_from; ///< time window settings from
	int time_to; ///< time window settings to
	ConflictMode conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records

	///sample by year from year to year by jump and find match in registry
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples);

	int init(map<string, string>& map);

	MedSamplingYearly() {
		gen = mt19937(rd());
		conflict_method = Drop; //default
		prediction_month_day = 101; //deafult
		back_random_duration = 0; //default
		day_jump = 0;
		mode = Case_Control;
		start_year = 0;
		end_year = 0;
		time_to = 0;
		time_from = 0;
	}
private:
	random_device rd;
	mt19937 gen;
};

/**
* A Class which samples by age from age to age by jump and find match in registry.
* suitble for incidence calculation
*/
class MedSamplingAge : public MedSamplingStrategy {
public:
	int start_age; ///< The start age to sample from
	int end_age; ///< The end age to sample from
	int age_bin; ///< the age bin in years for jumping
	SamplingMode mode; ///< sampling mode
	ConflictMode conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records

	///sample by year from age to age by jump and find match in registry
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples);

	MedSamplingAge() {
		start_age = 0;
		end_age = 120;
		age_bin = 1;
		conflict_method = All;
		mode = Case_Control;
	}

	int init(map<string, string>& map);

private:

};

/**
* Samples between given dates for ech patient
*/
class MedSamplingDates : public MedSamplingStrategy {
public:
	int take_count; ///< How many samples to take in each date
	SamplingMode mode; ///< sampling mode
	int time_from; ///< time window settings from
	int time_to; ///< time window settings to
	ConflictMode conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records
	vector<vector<pair<int, int>>> samples_list_pid_dates; ///< All sample options for pid,date to sample from. row is sample with all options to sample from 

	///sample Take_Count samples for each record in samples_list_pid_dates.
	///each record is vector<pair<int, int>> which is list of all options to choose from
	/// each record in the options is (pid, prediction_time)				
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples);

	int init(map<string, string>& map);

	MedSamplingDates() {
		gen = mt19937(rd()); take_count = 1; mode = Case_Control;
		conflict_method = All; time_from = 0; time_to = 0;
	}

private:
	random_device rd;
	mt19937 gen;
};


/**
* A Class which samples from start_time to end_time by jump and find match in registry.
* also suitble for incidence calculation
*/
class MedSamplingFixedTime : public MedSamplingStrategy {
public:
	int start_time; ///< The start time to sample from. If 0 will use min time of pid
	int end_time; ///< The end time to sample from. If 0 will use max time of pid
	int back_random_duration; ///< Random duration backward from prediciton month_day. to cancel use 0
	int time_jump; ///< the time jump, how much jump from each prediciton date
	SamplingMode mode; ///< sampling mode
	int time_from; ///< time window settings from
	int time_to; ///< time window settings to
	ConflictMode conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records

	///sample by year from year to year by jump and find match in registry
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples);

	int init(map<string, string>& map);

	MedSamplingFixedTime() {
		gen = mt19937(rd());
		conflict_method = Drop; //default
		back_random_duration = 0; //default
		time_jump = 0;
		mode = Case_Control;
		start_time = 0;
		end_time = 0;
		time_to = 0;
		time_from = 0;
	}
private:
	random_device rd;
	mt19937 gen;
};

#endif