/// @file
/// Sampling methods over MedRegistry Object
#ifndef __MED_SAMPLING_STRATEGY_H__
#define __MED_SAMPLING_STRATEGY_H__

#include <vector>
#include <random>
#include "LabelParams.h"
#include <InfraMed/InfraMed/MedPidRepository.h>
#include "MedEnums.h"

using namespace std;

/**
* An abstract class with sampling methods over registry records to convert to MedSamples
*/
class MedSamplingStrategy : public SerializableObject {
public:
	virtual void init_sampler(MedRepository &rep) {};
	
	/// The sampling options to be implemented
	virtual void get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, 
		unordered_map<int, vector<int>> &pid_options) const = 0;

	/// @snippet MedSamplingStrategy.cpp MedSamplingStrategy::make_sampler
	static MedSamplingStrategy *make_sampler(const string &sampler_name);
	/// @snippet MedSamplingStrategy.cpp MedSamplingStrategy::make_sampler
	static MedSamplingStrategy *make_sampler(const string &sampler_name, const string &init_params);

	virtual ~MedSamplingStrategy() {};
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

	void init_sampler(MedRepository &rep);

	/// sample random using Environment variable. params: [Random_Duration, Back_Time_Window_Years, Jump_Time_Period_Years]
	void get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	int init(map<string, string>& map);

	MedSamplingTimeWindow() {
		sample_count = 1; take_max = false;
	}
private:
	unordered_map<int, int> pids_bdates;
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
	
	///sample by year from year to year by jump and find match in registry
	void get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	int init(map<string, string>& map);

	MedSamplingYearly() {
		prediction_month_day = 101; //deafult
		back_random_duration = 0; //default
		day_jump = 0;
		start_year = 0;
		end_year = 0;
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
	void init_sampler(MedRepository &rep);

	///sample by year from age to age by jump and find match in registry
	void get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	MedSamplingAge() {
		start_age = 0;
		end_age = 120;
		age_bin = 1;
	}

	int init(map<string, string>& map);

private:
	unordered_map<int, int> pids_bdates;
};

/**
* Samples between given dates for ech patient
*/
class MedSamplingDates : public MedSamplingStrategy {
public:
	int take_count; ///< How many samples to take in each date
	vector<vector<pair<int, int>>> samples_list_pid_dates; ///< All sample options for pid,date to sample from. row is sample with all options to sample from 

	///sample Take_Count samples for each record in samples_list_pid_dates.
	///each record is vector<pair<int, int>> which is list of all options to choose from
	/// each record in the options is (pid, prediction_time)				
	void get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	int init(map<string, string>& map);

	MedSamplingDates() {
		take_count = 1;
	}
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
	
	///sample by year from year to year by jump and find match in registry
	void get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, unordered_map<int, vector<int>> &pid_options) const;

	int init(map<string, string>& map);

	MedSamplingFixedTime() {
		back_random_duration = 0; //default
		time_jump = 0;
		start_time = 0;
		end_time = 0;
	}
};

#endif