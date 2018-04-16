#ifndef __MED_SAMPLING_STRATEGY_H__
#define __MED_SAMPLING_STRATEGY_H__

#include <vector>
#include <random>
#include "MedRegistry.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"

class MedRegistryRecord;

using namespace std;

/**
* An abstract class with sampling methods over registry records to convert to MedSamples
*/
class MedSamplingStrategy : public SerializableObject {
public:
	/// The sampling method to be implemented
	virtual void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples) = 0;

	static MedSamplingStrategy *make_sampler(const string &sampler_name);
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
	bool use_allowed; ///< If True will check for registry time window intersection with min_allowed=>max_allowed. instead start=>end
	int allowed_time_from; ///< time window settings whem use_allowed is on
	int allowed_time_to; ///< time window settings whem use_allowed is on
	string conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records

	///sample by year from year to year by jump and find match in registry
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples);

	int init(map<string, string>& map);

	MedSamplingYearly() {
		gen = mt19937(rd());
		conflict_method = "drop"; //default
		prediction_month_day = 101; //deafult
		back_random_duration = 0; //default
		int day_jump = 0;
		use_allowed = false;
		conflict_method = "all";
		start_year = 0;
		end_year = 0;
		allowed_time_to = 0;
		allowed_time_from = 0;
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
	bool use_allowed; ///< If True will check for registry time window intersection with min_allowed=>max_allowed. instead start=>end
	string conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records

	///sample by year from age to age by jump and find match in registry
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples);

	MedSamplingAge() {
		start_age = 0;
		end_age = 120;
		age_bin = 1;
		conflict_method = "all";
		use_allowed = false;
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
	bool use_allowed; ///< If True will check for registry time window intersection with min_allowed=>max_allowed. instead start=>end
	int allowed_time_from; ///< time window settings whem use_allowed is on
	int allowed_time_to; ///< time window settings whem use_allowed is on
	string conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records
	vector<vector<pair<int, int>>> samples_list_pid_dates; ///< All sample options for pid,date to sample from. row is sample with all options to sample from 

	///sample Take_Count samples for each record in samples_list_pid_dates.
	///each record is vector<pair<int, int>> which is list of all options to choose from
	/// each record in the options is (pid, prediction_time)				
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples);

	int init(map<string, string>& map);

	MedSamplingDates() {
		gen = mt19937(rd()); take_count = 1; use_allowed = false;
		conflict_method = "all"; allowed_time_from = 0; allowed_time_to = 0;
	}

private:
	random_device rd;
	mt19937 gen;
};


#endif