#ifndef __MED_SAMPLING_HELPER__
#define __MED_SAMPLING_HELPER__

#include "FilterParams.h"
#include "MedRegistryRecord.h"
#include <MedProcessTools/MedProcessTools/MedSamples.h>

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
		/// @param time_from the time window from - to check with outcome registry
		/// @param time_to the time window to - to check with outcome registry
		/// @param censor_time_from the time window from - to check with censoring registry
		/// @param censor_time_to the time window to - to check with censoring registry
		/// @param mode the intersection method test for outcome
		/// @param mode_prediction the intersection method test for censoring
		/// </summary>
		/// <returns>
		/// If has intersection with time window
		/// </returns>
		bool in_time_window(int pred_date, const MedRegistryRecord *r_outcome, const vector<const MedRegistryRecord *> &r_censor,
			int time_from, int time_to, int censor_time_from, int  censor_time_to, 
			const TimeWindowMode mode[2], const TimeWindowMode mode_prediction[2]);

		/// <summary>
		/// checks for time range intersection
		/// @param pred_date prediction time
		/// @param r_outcome the registry record for label
		/// @param r_censor all the patient registry records for censoring. if empty - no censoring
		/// @param time_from the time window from - to check with outcome registry
		/// @param time_to the time window to - to check with outcome registry
		/// @param censor_time_from the time window from - to check with censoring registry
		/// @param censor_time_to the time window to - to check with censoring registry
		/// @param mode_outcome the intersection method test for outcome
		/// @param mode_censoring the intersection method test for censoring
		/// @param filter_no_censor what to do when no censoring record options are given
		/// </summary>
		/// <returns>
		/// If has intersection with time window
		/// </returns>
		bool in_time_window(int pred_date, const MedRegistryRecord *r_outcome, const vector<const MedRegistryRecord *> &r_censor,
			int time_from, int time_to, int censor_time_from, int  censor_time_to, const TimeWindowInteraction &mode_outcome, const TimeWindowInteraction &mode_censoring,
			bool filter_no_censor = true);

		/// <summary>
		/// checks for time range intersection
		/// @param pred_time prediction time
		/// @param pid_records the registry records of patient which are candidated for labeling
		/// @param r_censor all the patient registry records for censoring. if empty - no censoring
		/// @param time_from the time window from - to check with outcome registry
		/// @param time_to the time window to - to check with outcome registry
		/// @param censor_time_from the time window from - to check with censoring registry
		/// @param censor_time_to the time window to - to check with censoring registry
		/// @param mode_outcome the intersection method test for outcome
		/// @param mode_censoring the intersection method test for censoring
		/// @param filter_no_censor what to do when no censoring record options are given
		/// </summary>
		/// <returns>
		/// If has intersection with time window
		/// </returns>
		void get_label_for_sample(int pred_time, const vector<const MedRegistryRecord *> &pid_records
			, const vector<const MedRegistryRecord *> &r_censor, int time_from, int time_to, int censor_time_from, int  censor_time_to,
			const TimeWindowInteraction &mode_outcome, const TimeWindowInteraction &mode_censoring,
			ConflictMode conflict_mode, vector<MedSample> &idSamples,
			int &no_rule_found, int &conflict_count, int &done_count, bool treat_0_class_as_other_classes, bool filter_no_censor = true, bool show_conflicts = false);
	}
}

#endif // !__MED_SAMPLING_HELPER__

