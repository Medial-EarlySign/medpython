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

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	*  \brief sampling namespace
	*/
	namespace sampling {
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

	/*!
	* \brief signal_hierarchy namespace
	*/
	namespace signal_hierarchy {
		/// \brief Getting signal with hirarchy options for each siganl
		void getRecords_Hir(int pid, vector<UniversalSigVec> &signals, MedDictionarySections &dict,
			const string &signalHirerchyType,
			vector<MedRegistryRecord> &res);
	}
}

class SamplingRes {
public:
	int done_cnt = 0;
	int conflict_cnt = 0;
	int no_rule_cnt = 0;
	int miss_pid_in_reg_cnt = 0;
};

class MedSamplingStrategy;
static unordered_set<float> default_empty_set;
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
	void calc_signal_stats(const string &repository_path, int signalCode,
		const string &signalHirerchyType, int ageBinValue, MedSamplingStrategy &sampler,
		const LabelParams &inc_labeling_params, map<float, map<float, vector<int>>> &maleSignalToStats,
		map<float, map<float, vector<int>>> &femaleSignalToStats,
		const string &debug_file = "", const unordered_set<float> &debug_vals = default_empty_set) const;

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
	SamplingRes get_samples(int pid, int time, vector<MedSample> &samples) const;

	/// <summary>
	/// returns label value for time points
	/// </summary>
	SamplingRes get_samples(int pid, const vector<int> &times, vector<MedSample> &samples) const;

	void create_samples(const MedSamplingStrategy *sampler, MedSamples &samples) const;

	MedLabels(const LabelParams &params);

};

#endif
