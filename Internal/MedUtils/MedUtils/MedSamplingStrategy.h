/// @file
/// Sampling methods over MedRegistry Object
#ifndef __MED_SAMPLING_STRATEGY_H__
#define __MED_SAMPLING_STRATEGY_H__

#include <vector>
#include <random>
#include <MedUtils/MedUtils/MedRegistry.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include "MedEnums.h"

class MedRegistryRecord;

using namespace std;

/**
* A warpper class for initializing rules for time window interaction
*/
class TimeWindowInteraction {
private:
	bool has_default_mode;
	bool has_default_range;
	TimeWindowMode default_modes[2];
	pair<float, float> default_intersection_range;

public:
	unordered_map<float, TimeWindowMode[2]> interaction_mode; ///< map from label value to interaction rules
	unordered_map<float, pair<float, float>> intersection_range_condition; ///< intersection rate range with registry condition

	TimeWindowInteraction() {
		has_default_mode = false;
		has_default_range = false;
	}

	TimeWindowMode *operator[] (float x) {
		if (interaction_mode.find(x) != interaction_mode.end() || !has_default_mode)
			return interaction_mode[x];
		//has_default_mode is true and not exist
		return default_modes;
	}

	bool find(const float x) const {
		return has_default_mode || interaction_mode.find(x) != interaction_mode.end();
	}

	const TimeWindowMode *at(float x) const {
		if (interaction_mode.find(x) != interaction_mode.end() || !has_default_mode)
			return interaction_mode.at(x);
		return default_modes;
	}

	void set_default(TimeWindowMode defaults_m[2]) {
		if (has_default_mode)
			HMTHROW_AND_ERR("Error - TimeWindowInteraction has already default\n");
		has_default_mode = true;
		default_modes[0] = defaults_m[0];
		default_modes[1] = defaults_m[1];
	}

	void set_default_range(float min_range, float max_range) {
		default_intersection_range.first = min_range;
		default_intersection_range.second = max_range;
		has_default_range = true;
	}

	bool get_inresection_range_cond(float x, float &min_range, float &max_range) const {
		if (intersection_range_condition.find(x) != intersection_range_condition.end()) {
			min_range = intersection_range_condition.at(x).first;
			max_range = intersection_range_condition.at(x).second;
			return true;
		}
		else if (has_default_range) {
			min_range = default_intersection_range.first;
			max_range = default_intersection_range.second;
			return true;
		}
		return false;
	}

	void reset_for_init() {
		has_default_mode = false;
		has_default_range = false;
		interaction_mode.clear();
		intersection_range_condition.clear();
	}
};

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	*  \brief process namespace
	*/
	namespace process {
		/// <summary>
		/// checks for time range intersection
		/// @param pred_date prediction time
		/// @param start_time start time of window
		/// @param end_time end time of window
		/// @param reverse If we are in reverse window - looking backward in time
		/// @param mode the intersection method test
		/// </summary>
		/// <returns>
		/// If has intersection with time window
		/// </returns>
		bool in_time_window_simple(int pred_date, int start_time, int end_time, bool reverse, TimeWindowMode mode);

		/// <summary>
		/// checks for time range intersection
		/// @param pred_date prediction time
		/// @param r_outcome the registry record for label
		/// @param r_censor all the patient registry records for censoring. if empty - no censoring
		/// @param time_from the time window from - to check with censoring date
		/// @param time_to the time window to - to check with outcome registry
		/// @param mode the intersection method test for outcome
		/// @param mode_prediction the intersection method test for censoring
		/// </summary>
		/// <returns>
		/// If has intersection with time window
		/// </returns>
		bool in_time_window(int pred_date, const MedRegistryRecord *r_outcome, const vector<const MedRegistryRecord *> &r_censor,
			int time_from, int time_to, const TimeWindowMode mode[2], const TimeWindowMode mode_prediction[2]);

		/// <summary>
		/// checks for time range intersection
		/// @param pred_date prediction time
		/// @param r_outcome the registry record for label
		/// @param r_censor all the patient registry records for censoring. if empty - no censoring
		/// @param time_from the time window from - to check with censoring date
		/// @param time_to the time window to - to check with outcome registry
		/// @param mode_outcome the intersection method test for outcome
		/// @param mode_censoring the intersection method test for censoring
		/// @param filter_no_censor what to do when no censoring record options are given
		/// </summary>
		/// <returns>
		/// If has intersection with time window
		/// </returns>
		bool in_time_window(int pred_date, const MedRegistryRecord *r_outcome, const vector<const MedRegistryRecord *> &r_censor,
			int time_from, int time_to, const TimeWindowInteraction &mode_outcome, const TimeWindowInteraction &mode_censoring,
			bool filter_no_censor = true);
	}

	/*!
	*  \brief sampling namespace
	*/
	namespace sampling {
		/// <summary>
		/// Supports reading of complex map object in format: registry_val:mode_for_start,mode_for_end. use "|" to seperate labels
		/// has ability to give 3rd tokens in format number-number to specify range for intersection rate condition
		/// Example: "0:all,before_end|1:before_start,after_start"
		/// for complex labels has keyword "all" to activate rule on all labels till specific override
		/// </summary>
		void init_time_window_mode(const string &init, TimeWindowInteraction &mode);

		/// <summary>
		/// checks for time range intersection
		/// @param pred_time prediction time
		/// @param pid_records the registry records of patient which are candidated for labeling
		/// @param r_censor all the patient registry records for censoring. if empty - no censoring
		/// @param time_from the time window from - to check with censoring date
		/// @param time_to the time window to - to check with outcome registry
		/// @param mode_outcome the intersection method test for outcome
		/// @param mode_censoring the intersection method test for censoring
		/// @param filter_no_censor what to do when no censoring record options are given
		/// </summary>
		/// <returns>
		/// If has intersection with time window
		/// </returns>
		void get_label_for_sample(int pred_time, const vector<const MedRegistryRecord *> &pid_records
			, const vector<const MedRegistryRecord *> &r_censor, int time_from, int time_to,
			const TimeWindowInteraction &mode_outcome, const TimeWindowInteraction &mode_censoring,
			ConflictMode conflict_mode, vector<MedSample> &idSamples,
			int &no_rule_found, int &conflict_count, int &done_count, bool filter_no_censor = true);
	}
}

/**
* An abstract class with sampling methods over registry records to convert to MedSamples
*/
class MedSamplingStrategy : public SerializableObject {
public:
	virtual void init_sampler(MedRepository &rep) {};
	/// The sampling method to be implemented
	virtual void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry = NULL) = 0;

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

	void init_sampler(MedRepository &rep);

	///sample random using Environment variable. params: [Random_Duration, Back_Time_Window_Years, Jump_Time_Period_Years]
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry = NULL);

	int init(map<string, string>& map);

	MedSamplingTimeWindow() {
		gen = mt19937(rd()); sample_count = 1; take_max = false;
		maximal_time_case = 0; minimal_time_case = 0; minimal_time_control = 0; maximal_time_control = 0;
	}
private:
	random_device rd;
	mt19937 gen;
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
	TimeWindowInteraction outcome_interaction_mode; ///< sampling mode per label - for start and end
	TimeWindowInteraction censor_interaction_mode; ///< sampling mode per label - for start and end for censor
	int time_from; ///< time window settings from
	int time_to; ///< time window settings to
	ConflictMode conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records

	///sample by year from year to year by jump and find match in registry
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry = NULL);

	int init(map<string, string>& map);

	MedSamplingYearly() {
		gen = mt19937(rd());
		conflict_method = ConflictMode::Drop; //default
		prediction_month_day = 101; //deafult
		back_random_duration = 0; //default
		day_jump = 0;
		medial::sampling::init_time_window_mode("0:all,before_end|1:before_start,after_start", outcome_interaction_mode);
		medial::sampling::init_time_window_mode("all:within,within", censor_interaction_mode);
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
	TimeWindowInteraction outcome_interaction_mode; ///< sampling mode per label - for start and end
	TimeWindowInteraction censor_interaction_mode; ///< sampling mode per label - for start and end for censor
	ConflictMode conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records
	void init_sampler(MedRepository &rep);

	///sample by year from age to age by jump and find match in registry
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry = NULL);

	MedSamplingAge() {
		start_age = 0;
		end_age = 120;
		age_bin = 1;
		conflict_method = ConflictMode::All;
		medial::sampling::init_time_window_mode("0:all,before_end|1:before_start,after_start", outcome_interaction_mode);
		medial::sampling::init_time_window_mode("all:within,within", censor_interaction_mode);
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
	TimeWindowInteraction outcome_interaction_mode; ///< sampling mode per label - for start and end
	TimeWindowInteraction censor_interaction_mode; ///< sampling mode per label - for start and end for censor
	int time_from; ///< time window settings from
	int time_to; ///< time window settings to
	ConflictMode conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records
	vector<vector<pair<int, int>>> samples_list_pid_dates; ///< All sample options for pid,date to sample from. row is sample with all options to sample from 

	///sample Take_Count samples for each record in samples_list_pid_dates.
	///each record is vector<pair<int, int>> which is list of all options to choose from
	/// each record in the options is (pid, prediction_time)				
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry = NULL);

	int init(map<string, string>& map);

	MedSamplingDates() {
		gen = mt19937(rd()); take_count = 1;
		medial::sampling::init_time_window_mode("0:all,before_end|1:before_start,after_start", outcome_interaction_mode);
		medial::sampling::init_time_window_mode("all:within,within", censor_interaction_mode);
		conflict_method = ConflictMode::All; time_from = 0; time_to = 0;
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
	TimeWindowInteraction outcome_interaction_mode; ///< sampling mode per label - for start and end
	TimeWindowInteraction censor_interaction_mode; ///< sampling mode per label - for start and end for censor
	int time_from; ///< time window settings from
	int time_to; ///< time window settings to
	ConflictMode conflict_method; ///< options: all,max,drop how to treat intesections with multiple registry records

	///sample by year from year to year by jump and find match in registry
	void do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry = NULL);

	int init(map<string, string>& map);

	MedSamplingFixedTime() {
		gen = mt19937(rd());
		conflict_method = ConflictMode::Drop; //default
		back_random_duration = 0; //default
		time_jump = 0;
		medial::sampling::init_time_window_mode("0:all,before_end|1:before_start,after_start", outcome_interaction_mode);
		medial::sampling::init_time_window_mode("all:within,within", censor_interaction_mode);
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