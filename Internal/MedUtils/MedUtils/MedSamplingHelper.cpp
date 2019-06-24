#include "MedSamplingHelper.h"
#include <InfraMed/InfraMed/InfraMed.h>
#include <Logger/Logger/Logger.h>

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

bool medial::sampling::in_time_window_simple(int pred_date, int start_time, int end_time, bool reverse, TimeWindowMode mode) {
	switch (mode)
	{
	case TimeWindowMode::All_:
		return true;
	case TimeWindowMode::Before_End:
		if (reverse)
			return pred_date >= start_time;
		else
			return pred_date <= end_time;
	case TimeWindowMode::Before_Start:
		if (reverse)
			return pred_date >= end_time;
		else
			return pred_date <= start_time;
	case TimeWindowMode::After_Start:
		if (reverse)
			return (pred_date <= end_time);
		else
			return (pred_date >= start_time);
	case TimeWindowMode::Within:
		return  (pred_date >= start_time) && (pred_date <= end_time);
	default:
		MTHROW_AND_ERR("Error in in_time_window - unsupported mode - %d\n", mode);
	}
}

// testing for time_window - for specific registry_value. has rule for pred_date - which is from_time_window 
// time and rules for outcome
bool medial::sampling::in_time_window(int pred_date, const MedRegistryRecord *r_outcome, const vector<const MedRegistryRecord *> &r_censor,
	int time_from, int time_to, const TimeWindowMode mode[2], const TimeWindowMode mode_prediction[2]) {
	int sig_start_date = medial::repository::DateAdd(pred_date, time_from);
	int sig_end_date = medial::repository::DateAdd(pred_date, time_to);
	int reffer_date = sig_start_date, op_reffer = sig_end_date;
	if (time_from < 0) {//if looking backward force end_date to be in allowed
		reffer_date = sig_end_date;
		op_reffer = sig_start_date;
	}
	bool reverse = time_from < 0;
	//if (reffer_date > r->max_allowed_date || reffer_date < r->min_allowed_date)

	int idx_time = 0;
	bool can_have_pred = r_censor.empty();
	while (idx_time < r_censor.size() && !can_have_pred) {
		can_have_pred = in_time_window_simple(reffer_date, r_censor[idx_time]->start_date,
			r_censor[idx_time]->end_date, reverse, mode_prediction[0]);
		can_have_pred &= in_time_window_simple(op_reffer, r_censor[idx_time]->start_date,
			r_censor[idx_time]->end_date, reverse, mode_prediction[1]);
		++idx_time;
	}
	if (!can_have_pred)
		return false; //can't give prediction

	bool has_interact = in_time_window_simple(reffer_date, r_outcome->start_date, r_outcome->end_date, reverse, mode[0]);
	has_interact &= in_time_window_simple(op_reffer, r_outcome->start_date, r_outcome->end_date, reverse, mode[1]);
	return has_interact;
}

float interect_time_window(int pred_date, int time_from, int time_to,
	const MedRegistryRecord *r_outcome) {
	int sig_start_date = medial::repository::DateAdd(pred_date, time_from);
	int sig_end_date = medial::repository::DateAdd(pred_date, time_to);
	int start_window = min(sig_start_date, sig_end_date);
	int end_window = max(sig_start_date, sig_end_date);
	int window_size = abs(time_to - time_from);

	int max_start = max(start_window, r_outcome->start_date);
	int min_end = max(end_window, r_outcome->end_date);
	int interact_size = min_end - max_start;
	if (interact_size < 0)
		interact_size = 0;

	return float(interact_size) / window_size;
}

bool medial::sampling::in_time_window(int pred_date, const MedRegistryRecord *r_outcome, const vector<const MedRegistryRecord *> &r_censor,
	int time_from, int time_to, const TimeWindowInteraction &mode_outcome, const TimeWindowInteraction &mode_censoring,
	bool filter_no_censor) {
	const TimeWindowMode *mode = NULL;
	const TimeWindowMode  *mode_censor = NULL;
	if (mode_outcome.find(r_outcome->registry_value))
		mode = mode_outcome.at(r_outcome->registry_value);
	if (mode_censoring.find(r_outcome->registry_value))
		mode_censor = mode_censoring.at(r_outcome->registry_value);

	float min_range, max_range;
	bool has_interact = in_time_window(pred_date, r_outcome, r_censor, time_from, time_to, mode, mode_censor);
	if (mode_outcome.get_inresection_range_cond(r_outcome->registry_value, min_range, max_range)) {
		float intersect_rate = interect_time_window(pred_date, time_from, time_to, r_outcome);
		has_interact &= intersect_rate >= min_range && intersect_rate <= max_range;
	}
	if (mode_censoring.get_inresection_range_cond(r_outcome->registry_value, min_range, max_range)) {
		bool any = r_censor.empty() && !filter_no_censor;
		for (size_t i = 0; i < r_censor.size() && !any; ++i)
		{
			float intersect_rate = interect_time_window(pred_date, time_from, time_to, r_censor[i]);
			any = intersect_rate >= min_range && intersect_rate <= max_range;

		}
		has_interact &= any;
	}

	return has_interact;
}

void medial::sampling::get_label_for_sample(int pred_time, const vector<const MedRegistryRecord *> &pid_records
	, const vector<const MedRegistryRecord *> &r_censor, int time_from, int time_to,
	const TimeWindowInteraction &mode_outcome, const TimeWindowInteraction &mode_censoring,
	ConflictMode conflict_mode, vector<MedSample> &idSamples,
	int &no_rule_found, int &conflict_count, int &done_count, bool filter_no_censor, bool show_conflicts) {
	int curr_index = 0, final_selected = -1;
	float reg_val = -1;
	int reg_time = -1;
	if (pid_records.empty())
		return;
	MedSample smp;
	smp.time = pred_time;
	smp.id = pid_records.front()->pid;
	vector<const MedRegistryRecord *> matched_regs;

	//run on all matches:
	while (curr_index < pid_records.size()) {
		if (curr_index < pid_records.size()) {
			if (!mode_outcome.find(pid_records[curr_index]->registry_value)) {
#pragma omp atomic
				++no_rule_found;
				if (no_rule_found < 5)
					MWARN("Warning: missing rule for %f - skipping!!\n", pid_records[curr_index]->registry_value);
				++curr_index;
				continue;
			}

			if (!mode_censoring.find(pid_records[curr_index]->registry_value)) {
#pragma omp atomic
				++no_rule_found;
				if (no_rule_found < 5)
					MWARN("Warning: missing censor rule for %f - skipping!!\n", pid_records[curr_index]->registry_value);
				++curr_index;
				continue;
			}
		}
		if (curr_index < pid_records.size() &&
			!medial::sampling::in_time_window(pred_time, pid_records[curr_index], r_censor,
				time_from, time_to, mode_outcome, mode_censoring, filter_no_censor)) {
			++curr_index;
			continue;
		}
		if (curr_index >= pid_records.size())
			break; //skip if no match
				   //found match:
		if (reg_time == -1) { //first match
			reg_val = pid_records[curr_index]->registry_value;
			reg_time = pid_records[curr_index]->end_date;
			final_selected = curr_index;
			if (show_conflicts)
				matched_regs.push_back(pid_records[curr_index]);
		}
		else if (reg_val != pid_records[curr_index]->registry_value) {
			if (show_conflicts)
				matched_regs.push_back(pid_records[curr_index]);
			//if already found and conflicting:
			if (conflict_mode == ConflictMode::Drop) {
				reg_val = -1;
				reg_time = -1;
				final_selected = -1;
				break;
			}
			else if (conflict_mode == ConflictMode::Max) {
				if (reg_val < pid_records[curr_index]->registry_value) {
					reg_val = pid_records[curr_index]->registry_value;
					reg_time = pid_records[curr_index]->end_date;
					final_selected = curr_index;
				}
			}
			else if (conflict_mode == ConflictMode::All) {
				//insert current and update next:
				smp.outcomeTime = reg_val > 0 ? pid_records[curr_index]->start_date : reg_time;
				smp.outcome = reg_val;
				idSamples.push_back(smp);
#pragma omp atomic
				++done_count;
				reg_val = pid_records[curr_index]->registry_value;
				reg_time = pid_records[curr_index]->end_date;
				final_selected = curr_index;
			}
			else
				MTHROW_AND_ERR("Error in medial::sampling::get_label_for_sample - Unsupported conflict method %d\n", conflict_mode);
#pragma omp atomic
			++conflict_count;
			//break;
		}

		++curr_index;
	}

	if (reg_time != -1) {
		smp.outcomeTime = reg_val > 0 ? pid_records[final_selected]->start_date : reg_time;
		smp.outcome = reg_val;
		idSamples.push_back(smp);
		++done_count;
	}

	if (show_conflicts && matched_regs.size() >= 2) {
		string buffer_str = "";
		if (!matched_regs.empty())
			buffer_str += "Time:[" + to_string(matched_regs[0]->start_date) +
			"," + to_string(matched_regs[0]->end_date) + "],Label:" + to_string(matched_regs[0]->registry_value);
		for (size_t i = 1; i < matched_regs.size(); ++i)
			buffer_str += "|Time:[" + to_string(matched_regs[i]->start_date) +
			"," + to_string(matched_regs[i]->end_date) + "],Label:" + to_string(matched_regs[i]->registry_value);

		MWARN("Conflict example for pid %d on time %d - Registry:[%s] \n", smp.id, pred_time, buffer_str.c_str());
	}
}