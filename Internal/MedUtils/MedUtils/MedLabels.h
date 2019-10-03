#ifndef __MED_LABELS_H__
#define __MED_LABELS_H__
/// @file
/// Labeling methods over MedRegistry Object

#include <vector>
#include <unordered_map>
#include "MedRegistryRecord.h"
#include "MedEnums.h"
#include "LabelParams.h"
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include "MedSamplingStrategy.h"

using namespace std;

class SamplingRes {
public:
	int done_cnt = 0;
	int conflict_cnt = 0;
	int no_rule_cnt = 0;
	int miss_pid_in_reg_cnt = 0;
};

class MedSamplingStrategy;
static unordered_set<int> default_empty_set;
/**
* A Class which represent time ranges of label values based on registry and labeling method.
* The main procedure is to quety label value using get_label method
*/
class MedLabels {
private:
	vector<MedRegistryRecord> all_reg_records; ///< a copy of all registry records
	vector<MedRegistryRecord> all_censor_records; ///< a copy of all censor records
	unordered_map<int, vector<const MedRegistryRecord *>> pid_reg_records; ///< registry records aggregated by pid
	unordered_map<int, vector<const MedRegistryRecord *>> pid_censor_records; ///< registry records aggregated by pid
public:
	LabelParams labeling_params; ///< the labeling parameters - problem definition

	/// <summary>.
	/// prepare MedLabels from registry (and censor if provided) based on labeling rules
	/// </summary>
	void prepare_from_registry(const vector<MedRegistryRecord> &reg_records, const vector<MedRegistryRecord> *censor_records = NULL);

	/// <summary>
	/// return true if censor registry was provided
	/// </summary>
	bool has_censor_reg() const;

	/// <summary>
	/// return true if found censor records for patient
	/// </summary>
	bool has_censor_reg(int pid) const;

	/// <summary>
	/// return all availbel pids from registry
	/// </summary>
	void get_pids(vector<int> &pids) const;

	/// <summary>
	/// get all registry and censor records for patient
	/// </summary>
	void get_records(int pid, vector<const MedRegistryRecord *> &reg_records, vector<const MedRegistryRecord *> &censor_records) const;

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
	void calc_signal_stats(const string &repository_path, const string &signal_name,
		const string &signalHirerchyType, int ageBinValue, MedSamplingStrategy &sampler,
		const LabelParams &inc_labeling_params, map<float, map<float, vector<int>>> &maleSignalToStats,
		map<float, map<float, vector<int>>> &femaleSignalToStats,
		const string &debug_file = "", const unordered_set<int> &debug_vals = default_empty_set) const;

	/// <summary>
	/// calculate incidence and writes the result into file with old and new format
	/// @param file_path the output file path to write the results
	/// @param rep_path the repository path to calculate the incidence
	/// @param age_bin the age_bin for binning age groups for the incidence
	/// @param min_age the minimal age fro the incidence
	/// @param max_age the maximal age fro the incidence
	/// @param use_kaplan_meir if True will calc using kaplan meier survivol rates
	/// @param sampler_name the sampler name for calculating incidence
	/// @param sampler_args the sampler args for calculating incidence - may control trail years for example
	/// </summary>
	void create_incidence_file(const string &file_path, const string &rep_path, int age_bin, int min_age,
		int max_age, bool use_kaplan_meir = false, const string &sampler_name = "yearly",
		const string &sampler_args = "day_jump=365;start_year=2007;end_year=2012;prediction_month_day=101",
		const string &debug_file = "") const;

	/// <summary>
	/// returns label value for time point
	/// </summary>
	SamplingRes get_samples(int pid, int time, vector<MedSample> &samples, bool show_conflicts = false) const;

	/// <summary>
	/// update outcome and outcomeTime keeping everything else
	/// </summary>
	SamplingRes get_samples(int pid, const vector<MedSample> &samples, vector<MedSample> &new_samples, bool show_conflicts = false) const;

	/// <summary>
	/// returns label value for time points
	/// </summary>
	SamplingRes get_samples(int pid, const vector<int> &times, vector<MedSample> &samples, bool show_conflicts = false) const;

	/// <summary>
	///Creates MedSamples using MedSampling
	/// </summary>
	void create_samples(const MedSamplingStrategy *sampler, MedSamples &samples, bool show_conflicts = true) const;

	/// <summary>
	/// relabels the samples besed on the labeler
	/// </summary>
	void relabel_samples(MedSamples &samples, bool show_conflicts = true) const;

	/// <summary>
	/// Ctor with labling params
	/// </summary>
	MedLabels(const LabelParams &params);

	/// <summary>
	/// Ctor with labling params as string
	/// </summary>
	MedLabels(const string &labling_params);
};

#endif
