/// @file
/// Sampling methods over MedRegistry Object
#ifndef __MED_SAMPLING_STRATEGY_H__
#define __MED_SAMPLING_STRATEGY_H__

#include <vector>
#include <random>
#include "LabelParams.h"
#include <InfraMed/InfraMed/MedPidRepository.h>
#include "MedEnums.h"
#include "FilterParams.h"

using namespace std;

/**
* An abstract class with sampling methods over registry records to convert to MedSamples
*/
class MedSamplingStrategy : public SerializableObject {
public:
	FilterParams filtering_params; ///< the filtering constraints prior to sampling

	/// The sampler need repository for Age  filters if exist and some samplers also uses bdate.
	virtual void init_sampler(MedRepository &rep);

	/// The sampling options - calls _get_sampling_options, before that applies filters
	void get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
		unordered_map<int, vector<int>> &pid_options) const;

	/// The specific sampler get_options to implement
	virtual void _get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
		unordered_map<int, vector<int>> &pid_options) const = 0;

	/// @snippet MedSamplingStrategy.cpp MedSamplingStrategy::make_sampler
	static MedSamplingStrategy *make_sampler(const string &sampler_name);
	/// @snippet MedSamplingStrategy.cpp MedSamplingStrategy::make_sampler
	static MedSamplingStrategy *make_sampler(const string &sampler_name, const string &init_params);

	/// <summary>
	/// stores filters prior to sampling
	/// </summary>
	void set_filters(const string &filtering_str);

	virtual ~MedSamplingStrategy() {};
protected:
	unordered_map<int, int> pids_bdates;
private:
	bool apply_filter_params(unordered_map<int, vector<pair<int, int>>> &pid_time_ranges) const;
};

/**
* A Class which samples records on registry for certain time window.
* For Each registry record it samples randomly time point which falls withing the time window
* from min_allowed to max_allowed or till max_allowed backward using the time window params
*/
class MedSamplingTimeWindow : public MedSamplingStrategy {
public:
	bool take_max; ///< If true will random sample between all time range of min_allowed to max_allowed
	vector<int> minimal_times; ///< minimal times for window options
	vector<int> maximal_times; ///< maximal times for window options
	int sample_count; ///< how many samples to take in each time window

	/// sample random using Environment variable. params: [Random_Duration, Back_Time_Window_Years, Jump_Time_Period_Years]
	void _get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	int init(map<string, string>& map);

	MedSamplingTimeWindow() {
		sample_count = 1; take_max = false;
	}
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

	///sample by year from age to age by jump and find match in registry
	void _get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	MedSamplingAge() {
		start_age = 0;
		end_age = 120;
		age_bin = 1;
	}

	int init(map<string, string>& map);
};

/**
* Samples between given dates for ech patient
*/
class MedSamplingDates : public MedSamplingStrategy {
public:
	int take_count; ///< How many samples to take in each date
	vector<vector<pair<int, int>>> samples_list_pid_dates; ///< All sample options for pid,date to sample from. row is sample with all options to sample from 
	bool sample_with_filters; ///< If True will do sampling after time range filtering of years,age,censoring. otherwise total randomally choose times

	void get_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, 
		const vector<vector<pair<int, int>>> &samples_list_pid_opts,
		unordered_map<int, vector<int>> &pid_options) const;

	///sample Take_Count samples for each record in samples_list_pid_dates.
	///each record is vector<pair<int, int>> which is list of all options to choose from
	/// each record in the options is (pid, prediction_time)				
	void _get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	int init(map<string, string>& map);

	MedSamplingDates() {
		take_count = 1;
		sample_with_filters = true;
	}
};

/**
* A Class which samples from start_time to end_time by jump and find match in registry.
* also suitble for incidence calculation
*/
class MedSamplingFixedTime : public MedSamplingStrategy {
private:
	int add_time(int time, int add) const;
public:
	int start_time; ///< The start time to sample from. If 0 will use min time of pid
	int end_time; ///< The end time to sample from. If 0 will use max time of pid
	int back_random_duration; ///< Random duration backward from prediciton month_day. to cancel use 0
	int time_jump; ///< the time jump, how much jump from each prediciton date
	int time_range_unit; ///< the start_time,end_time unit
	int time_jump_unit; ///< the time jump unit

	///sample by year from year to year by jump and find match in registry
	void _get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	int init(map<string, string>& map);

	MedSamplingFixedTime() {
		back_random_duration = 0; //default
		time_jump = 0;
		start_time = 0;
		end_time = 0;
		time_jump_unit = global_default_windows_time_unit;
		time_range_unit = global_default_time_unit;
	}
};

/**
* DEPRECATED - A Class which samples by year from year to year by jump and find match in registry.
* suitble for incidence calculation
* uses MedSamplingFixedTime - just for backward compitablilty - DEPRECATED
*/
class MedSamplingYearly : public MedSamplingFixedTime {
public:
	int init(map<string, string>& map);
};

/**
* A Sampler to sample on one of signals test.
* You may also look at this example to create more complicated rules.
* All you need is to fetch and prepare the right data from repository with init_sampler and populate the values in samples_list_pid_dates
* Uses first time channel of each signal
*/
class MedSamplingStick : public MedSamplingDates {
public:
	vector<string> signal_list; ///< list of signals to take times for sampling on each patient
	int delta_time = 0; ///< delta time before the sticked signals. a date before should be negative
	int delta_time_unit = MedTime::Days;
	int minimal_time_between_samples = 0; ///< minimal time restriction between samples - starts from most recent and takes most recent in each time window

	/// Initialize samples_list_pid_dates by reading signals from repository
	void init_sampler(MedRepository &rep);

	int init(map<string, string>& map);

	void _get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	MedSamplingStick();
};

#endif