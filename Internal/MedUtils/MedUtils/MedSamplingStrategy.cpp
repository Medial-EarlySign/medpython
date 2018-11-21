#include <MedUtils/MedUtils/MedSamplingStrategy.h>
#include <fstream>
#include <string>
#include <iostream>
#include <Logger/Logger/Logger.h>

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

vector<string> TimeWindow_to_name = { "before", "before_start" ,"pass", "within", "all" };
vector<string> ConflictMode_to_name = { "all", "drop", "max" };

TimeWindowMode TimeWindow_name_to_type(const string& TimeWindow_name) {
	for (int i = 0; i < TimeWindow_to_name.size(); ++i)
		if (TimeWindow_to_name[i] == TimeWindow_name) {
			return TimeWindowMode(i);
		}
	MTHROW_AND_ERR("Error in SamplingMode_name_to_type - Unsupported \"%s\". options are: %s\n",
		TimeWindow_name.c_str(), medial::io::get_list(TimeWindow_to_name).c_str());
}
ConflictMode ConflictMode_name_to_type(const string& ConflictMode_name) {
	for (int i = 0; i < ConflictMode_to_name.size(); ++i)
		if (ConflictMode_to_name[i] == ConflictMode_name) {
			return ConflictMode(i);
		}
	MTHROW_AND_ERR("Error in SamplingMode_name_to_type - Unsupported \"%s\". options are: %s\n",
		ConflictMode_name.c_str(), medial::io::get_list(ConflictMode_to_name).c_str());
}

int MedSamplingTimeWindow::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "take_max")
			take_max = stoi(it->second) > 0;
		else if (it->first == "minimal_time_case")
			minimal_time_case = stoi(it->second);
		else if (it->first == "maximal_time_case")
			maximal_time_case = stoi(it->second);
		else if (it->first == "minimal_time_control")
			minimal_time_control = stoi(it->second);
		else if (it->first == "maximal_time_control")
			maximal_time_control = stoi(it->second);
		else if (it->first == "sample_count")
			sample_count = stoi(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (minimal_time_case < 0 || maximal_time_case < 0
		|| minimal_time_control < 0 || maximal_time_control < 0)
		MTHROW_AND_ERR("negative time_back values aren't allowed\n");
	if (minimal_time_control <= maximal_time_case)
		MWARN("Warning: minimal_time_control <= maximal_time_case.\n"
			"can sample on same date control and case\n");
	return 0;
}

void get_bdates(MedRepository &rep, unordered_map<int, int> &bdates) {
	int bDateCode = rep.sigs.sid("BDATE");
	int bYearCode = rep.sigs.sid("BYEAR");
	int use_code = bDateCode;
	if (rep.pids.empty() || bDateCode <= 0)
		MTHROW_AND_ERR("Error MedSamplingStrategy::get_bdates - repository wasn't initialized and contains BDATE\n");
	if (!rep.index.index_table[bDateCode].is_loaded) {
		use_code = bYearCode;
		if (!rep.index.index_table[bYearCode].is_loaded)
			MTHROW_AND_ERR("Error MedSamplingStrategy::get_bdates - repository wasn't loaded with BDATE or BYEAR\n");
	}
	for (size_t i = 0; i < rep.pids.size(); ++i)
	{
		int pid = rep.pids[i];
		int bdate_val = medial::repository::get_value(rep, pid, use_code);
		if (use_code == bYearCode)
			bdate_val = med_time_converter.convert_years(global_default_time_unit, bdate_val);
		bdates[pid] = bdate_val;
	}
	MLOG_D("MedSamplingStrategy::get_bdates - loaded %zu patients\n", bdates.size());
}

void MedSamplingTimeWindow::init_sampler(MedRepository &rep) {
	get_bdates(rep, pids_bdates);
}

void MedSamplingTimeWindow::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry) {
	int random_back_dur = 1;
	int diff_window_cases = maximal_time_case - minimal_time_case;
	int diff_window_controls = maximal_time_control - minimal_time_control;
	bool use_random = !take_max && (diff_window_cases > 1 || diff_window_controls > 1);
	unordered_map<int, vector<const MedRegistryRecord *>> pid_censor_dates;
	if (censor_registry != NULL)
		for (const MedRegistryRecord &rec : *censor_registry)
			pid_censor_dates[rec.pid].push_back(&rec);
	else
		MWARN("Warning MedSamplingTimeWindow::do_sample - no censor registry\n");

	//create samples file:
	unordered_map<int, int> pid_to_ind;
	int skip_end_smaller_start = 0, skip_no_bdate = 0, example_pid = -1, no_censor = 0;
	for (const MedRegistryRecord &rec : registry)
	{
		vector<const MedRegistryRecord *> *pid_dates = NULL;
		if (pid_censor_dates.find(rec.pid) != pid_censor_dates.end())
			pid_dates = &pid_censor_dates[rec.pid];
		else
			++no_censor;
		int max_allowed_date = rec.end_date;
		int min_allowed_date = rec.start_date;
		if (pid_dates != NULL) {
			for (size_t i = 0; i < pid_dates->size(); ++i)
				if ((*pid_dates)[i]->end_date > max_allowed_date)
					max_allowed_date = (*pid_dates)[i]->end_date;
			for (size_t i = 0; i < pid_dates->size(); ++i)
				if ((*pid_dates)[i]->start_date < min_allowed_date)
					min_allowed_date = (*pid_dates)[i]->start_date;
		}

		bool addNew = false;
		int currDate = rec.end_date;
		if (currDate > max_allowed_date)
			currDate = max_allowed_date;
		int diff_window = diff_window_cases;
		if (rec.registry_value > 0)
			currDate = medial::repository::DateAdd(currDate, -minimal_time_case);
		else {
			currDate = medial::repository::DateAdd(currDate, -minimal_time_control);
			diff_window = diff_window_controls;
		}

		int pid_bdate = -1;
		if (pids_bdates.find(rec.pid) != pids_bdates.end())
			pid_bdate = pids_bdates.at(rec.pid);
		else {
			++skip_no_bdate;
			example_pid = rec.pid;
			continue;
		}
		float year_diff_to_first_pred;
		if (min_allowed_date <= 0) //has no limit - if "max" go back until date of birth
			year_diff_to_first_pred = medial::repository::DateDiff(pid_bdate, currDate);
		else
			year_diff_to_first_pred = medial::repository::DateDiff(min_allowed_date, currDate);
		if (year_diff_to_first_pred < 0 || rec.end_date <= rec.start_date || rec.end_date <= min_allowed_date) {
			++skip_end_smaller_start;
			if (skip_end_smaller_start < 5) {
				MLOG("Exampled Row Skipped: pid=%d, reg_dates=[%d => %d],  outcome=%f, age=%d\n",
					rec.pid, rec.start_date, rec.end_date, rec.registry_value,
					(int)medial::repository::DateDiff(pid_bdate, currDate));
			}
			continue;
		}
		int min_pred_date; //how many years to go back
		if (diff_window > 365 * year_diff_to_first_pred) //validate we wont go back too far
			diff_window = int(365 * year_diff_to_first_pred); //window passed max allowed - so cut in max
		min_pred_date = medial::repository::DateAdd(currDate, -diff_window); //how many years to go back
		MedIdSamples patient_samples(rec.pid);
		if (pid_to_ind.find(rec.pid) == pid_to_ind.end()) {
			pid_to_ind[rec.pid] = (int)samples.idSamples.size();
			addNew = true;
		}

		int rnd_days_diff = 0;
		if (take_max || use_random) {
			if (diff_window < 1) { //not enought time to sample - skip
				if (addNew && patient_samples.samples.empty()) //was new and haven't been added yet
					pid_to_ind.erase(rec.pid);
				continue;
			}
			uniform_int_distribution<> rand_int(0, diff_window);
			rnd_days_diff = (int)rand_int(gen);
		}

		for (size_t i = 0; i < sample_count; ++i)
		{
			int sample_pred_date = medial::repository::DateAdd(currDate, -rnd_days_diff);
			MedSample smp;
			smp.id = rec.pid;
			smp.outcome = rec.registry_value;
			smp.outcomeTime = rec.registry_value > 0 ? rec.start_date : rec.end_date;
			smp.split = 0;
			smp.time = sample_pred_date;
			patient_samples.samples.push_back(smp);
		}

		if (addNew) {
			if (!patient_samples.samples.empty())
				samples.idSamples.push_back(patient_samples);
		}
		else
			samples.idSamples[pid_to_ind[rec.pid]].samples.insert(samples.idSamples[pid_to_ind[rec.pid]].samples.end(),
				patient_samples.samples.begin(), patient_samples.samples.end());
	}
	samples.sort_by_id_date();

	if (no_censor > 0)
		MLOG("WARNING MedSamplingTimeWindow:do_sample - has %d samples with no censor dates\n", no_censor);
	if (skip_no_bdate > 0)
		MLOG("WARNING :: Skipped %d registry records because no bdate: example pid=%d\n", skip_no_bdate, example_pid);
	if (skip_end_smaller_start > 0)
		MLOG("WARNING :: Skipped %d registry records because end_date<start_date\n", skip_end_smaller_start);
}

int MedSamplingYearly::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "start_year") {
			start_year = stoi(it->second);
			if (start_year <= 1900 || start_year >= 2100)
				MTHROW_AND_ERR("start_year must be initialize between 1900 to 2100\n");
		}
		else if (it->first == "end_year") {
			end_year = stoi(it->second);
			if (end_year <= 1900 || end_year >= 2100)
				MTHROW_AND_ERR("end_year must be initialize between 1900 to 2100\n");
		}
		else if (it->first == "day_jump") {
			day_jump = stoi(it->second);
			if (day_jump <= 0)
				MTHROW_AND_ERR("day_jump must be positive > 0\n");
		}
		else if (it->first == "prediction_month_day")
			prediction_month_day = stoi(it->second);
		else if (it->first == "back_random_duration")
			back_random_duration = stoi(it->second);
		else if (it->first == "time_from")
			time_from = stoi(it->second);
		else if (it->first == "time_to")
			time_to = stoi(it->second);
		else if (it->first == "outcome_interaction_mode")
			medial::sampling::init_time_window_mode(it->second, outcome_interaction_mode);
		else if (it->first == "censor_interaction_mode")
			medial::sampling::init_time_window_mode(it->second, censor_interaction_mode);
		else if (it->first == "conflict_method")
			conflict_method = ConflictMode_name_to_type(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}

	if (prediction_month_day < 100 || prediction_month_day % 100 > 31)
		MTHROW_AND_ERR("prediction_month_day must be positive >= 100 <=1231\n");
	if (back_random_duration < 0)
		MTHROW_AND_ERR("back_random_duration must be positive\n");
	return 0;
}

bool medial::process::in_time_window_simple(int pred_date, int start_time, int end_time, bool reverse, TimeWindowMode mode) {
	switch (mode)
	{
	case TimeWindowMode::All_:
		return true;
	case TimeWindowMode::Before:
		if (reverse)
			return pred_date >= start_time;
		else
			return pred_date <= end_time;
	case TimeWindowMode::Before_Start:
		if (reverse)
			return pred_date >= end_time;
		else
			return pred_date <= start_time;
	case TimeWindowMode::Pass:
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
bool medial::process::in_time_window(int pred_date, const MedRegistryRecord *r_outcome, const vector<const MedRegistryRecord *> &r_censor,
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

bool medial::process::in_time_window(int pred_date, const MedRegistryRecord *r_outcome, const vector<const MedRegistryRecord *> &r_censor,
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

void medial::sampling::init_time_window_mode(const string &init, TimeWindowInteraction &mode) {
	mode.reset_for_init();
	vector<string> tokens;
	boost::split(tokens, init, boost::is_any_of("|"));
	for (size_t i = 0; i < tokens.size(); ++i)
	{
		vector<string> tokens_inner, tokens_rules, intersection_tokens;
		//Format of tokens[i] is: "label:start,end"
		boost::split(tokens_inner, tokens[i], boost::is_any_of(":"));
		if (tokens_inner.size() != 2)
			MTHROW_AND_ERR("Error in medial::sampling::init_time_window_mode - reading token \"%s\" and missing"
				" \":\". format should be label:start,end(,num-num as optional)\n", tokens[i].c_str());
		const string &label = tokens_inner[0];
		boost::split(tokens_rules, tokens_inner[1], boost::is_any_of(","));
		if (tokens_rules.size() != 2 && tokens_rules.size() != 3)
			MTHROW_AND_ERR("Error in medial::sampling::init_time_window_mode - reading token \"%s\" and missing"
				" \",\". format should be start,end. full_token = \"%s\"\n", tokens_inner[1].c_str(), tokens[i].c_str());
		if (label == "all") {
			//mode
			TimeWindowMode temp_mode[2];
			temp_mode[0] = TimeWindow_name_to_type(tokens_rules[0]);
			temp_mode[1] = TimeWindow_name_to_type(tokens_rules[1]);

			mode.set_default(temp_mode);
		}
		else {
			mode[med_stof(label)][0] = TimeWindow_name_to_type(tokens_rules[0]);
			mode[med_stof(label)][1] = TimeWindow_name_to_type(tokens_rules[1]);
		}
		if (tokens_rules.size() == 3) {
			//aditional args for intersection:
			boost::split(intersection_tokens, tokens_rules[2], boost::is_any_of("-"));
			if (intersection_tokens.size() != 2)
				MTHROW_AND_ERR("Error in medial::sampling::init_time_window_mode - reading token \"%s\" and missing"
					" \",\". format should be number-number. full_token = \"%s\"\n",
					tokens_rules[2].c_str(), tokens[i].c_str());
			if (label != "all") {
				mode.intersection_range_condition[med_stof(label)].first = med_stof(intersection_tokens[0]);
				mode.intersection_range_condition[med_stof(label)].second = med_stof(intersection_tokens[1]);
			}
			else
				mode.set_default_range(med_stof(intersection_tokens[0]), med_stof(intersection_tokens[1]));
		}
	}


}

void MedSamplingYearly::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples,
	const vector<MedRegistryRecord> *censor_registry) {
	if (day_jump <= 0)
		MTHROW_AND_ERR("day_jump must be positive > 0\n");
	if (end_year <= 1900 || end_year >= 2100 || start_year <= 1900 || start_year >= 2100)
		MTHROW_AND_ERR("start_year,end_year must be initialize between 1900 to 2100\n");
	int random_back_dur = 1;
	bool use_random = back_random_duration > 0;
	if (use_random)
		random_back_dur = back_random_duration;

	uniform_int_distribution<> rand_int(0, random_back_dur);
	unordered_map<int, int> pid_to_ind;
	vector<MedIdSamples> idSamples;
	vector<const MedRegistryRecord *> empty_censor;

	unordered_map<int, vector<const MedRegistryRecord *>> pid_to_regs, pid_to_censor;
	for (size_t i = 0; i < registry.size(); ++i)
		pid_to_regs[registry[i].pid].push_back(&registry[i]);
	if (censor_registry != NULL)
		for (size_t i = 0; i < censor_registry->size(); ++i)
			pid_to_censor[(*censor_registry)[i].pid].push_back(&(*censor_registry)[i]);

	int conflict_count = 0, done_count = 0, no_censor = 0, no_rule = 0;
	for (auto it = pid_to_regs.begin(); it != pid_to_regs.end(); ++it)
	{
		vector<const MedRegistryRecord *> *all_pid_records = &it->second;
		int min_date = start_year;
		int max_date = end_year;
		long start_date = min_date * 10000 + prediction_month_day;
		long end_date = max_date * 10000 + prediction_month_day;
		if (pid_to_ind.find(it->first) == pid_to_ind.end()) {
			pid_to_ind[it->first] = (int)idSamples.size();
			MedIdSamples pid_sample(it->first);
			idSamples.push_back(pid_sample);
		}
		vector<const MedRegistryRecord *> *r_censor = &empty_censor;
		if (pid_to_censor.find(it->first) != pid_to_censor.end())
			r_censor = &pid_to_censor[it->first];
		else
			++no_censor;
		if (censor_registry != NULL && pid_to_censor.find(it->first) == pid_to_censor.end())
			continue; //filter sample
		for (long date = start_date; date <= end_date; date = medial::repository::DateAdd(date, day_jump)) {
			//search for match in all regs:
			int pred_date = date;
			if (use_random)
				pred_date = medial::repository::DateAdd(pred_date, -rand_int(gen));

			MedSample smp;
			smp.id = it->first;
			smp.time = pred_date;
			int curr_index = 0, final_selected = -1;
			float reg_val = -1;
			int reg_time = -1;
			//run on all matches:
			while (curr_index < all_pid_records->size()) {
				if (curr_index < all_pid_records->size()) {
					if (!outcome_interaction_mode.find((*all_pid_records)[curr_index]->registry_value)) {
						++no_rule;
						if (no_rule < 5)
							MWARN("Warning: missing rule for %f - skipping!!\n", (*all_pid_records)[curr_index]->registry_value);
						++curr_index;
						continue;
					}

					if (!censor_interaction_mode.find((*all_pid_records)[curr_index]->registry_value)) {
						++no_rule;
						if (no_rule < 5)
							MWARN("Warning: missing censor rule for %f - skipping!!\n", (*all_pid_records)[curr_index]->registry_value);
						++curr_index;
						continue;
					}
				}
				if (curr_index < all_pid_records->size() &&
					!medial::process::in_time_window(pred_date, (*all_pid_records)[curr_index], *r_censor,
						time_from, time_to, outcome_interaction_mode, censor_interaction_mode)) {
					++curr_index;
					continue;
				}
				if (curr_index >= all_pid_records->size())
					break; //skip if no match
				//found match:
				if (reg_time == -1) { //first match
					reg_val = (*all_pid_records)[curr_index]->registry_value;
					reg_time = (*all_pid_records)[curr_index]->end_date;
					final_selected = curr_index;
				}
				else if (reg_val != (*all_pid_records)[curr_index]->registry_value) {
					//if already found and conflicting:
					if (conflict_method == Drop) {
						reg_val = -1;
						reg_time = -1;
						final_selected = -1;
						break;
					}
					else if (conflict_method == Max) {
						if (reg_val < (*all_pid_records)[curr_index]->registry_value) {
							reg_val = (*all_pid_records)[curr_index]->registry_value;
							reg_time = (*all_pid_records)[curr_index]->end_date;
							final_selected = curr_index;
						}
					}
					else {
						//insert current and update next:
						smp.outcomeTime = reg_val > 0 ? (*all_pid_records)[curr_index]->start_date : reg_time;
						smp.outcome = reg_val;
						idSamples[pid_to_ind.at(it->first)].samples.push_back(smp);
						++done_count;
						reg_val = (*all_pid_records)[curr_index]->registry_value;
						reg_time = (*all_pid_records)[curr_index]->end_date;
						final_selected = curr_index;
					}
					++conflict_count;
					break;
				}

				++curr_index;
			}

			if (reg_time != -1) {
				smp.outcomeTime = reg_val > 0 ? (*all_pid_records)[final_selected]->start_date : reg_time;
				smp.outcome = reg_val;
				idSamples[pid_to_ind.at(it->first)].samples.push_back(smp);
				++done_count;
			}
		}
	}

	if (no_rule > 0)
		MLOG("WARNING MedSamplingYearly:do_sample - has %d samples with no rules for time window\n", no_rule);
	if (no_censor > 0)
		if (censor_registry != NULL)
			MLOG("WARNING MedSamplingYearly:do_sample - has %d patients with no censor dates\n", no_censor);
		else
			MLOG("WARNING MedSamplingYearly:do_sample - no censoring time region was given\n");
	if (conflict_count > 0)
		MLOG("Sampled registry with %d conflicts. has %d registry records\n", conflict_count, done_count);
	//keep non empty pids:
	for (int i = (int)idSamples.size() - 1; i >= 0; --i)
		if (!idSamples[i].samples.empty())
			samples.idSamples.push_back(idSamples[i]);
	samples.sort_by_id_date();
}

int MedSamplingAge::init(map<string, string>& map) {
	conflict_method = Drop; //default
	age_bin = 1; //deafult
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "start_age")
			start_age = stoi(it->second);
		else if (it->first == "end_age")
			end_age = stoi(it->second);
		else if (it->first == "age_bin")
			age_bin = stoi(it->second);
		else if (it->first == "conflict_method")
			conflict_method = ConflictMode_name_to_type(it->second);
		else if (it->first == "outcome_interaction_mode")
			medial::sampling::init_time_window_mode(it->second, outcome_interaction_mode);
		else if (it->first == "censor_interaction_mode")
			medial::sampling::init_time_window_mode(it->second, censor_interaction_mode);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	return 0;
}

void MedSamplingAge::init_sampler(MedRepository &rep) {
	get_bdates(rep, pids_bdates);
}

void MedSamplingAge::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry) {
	if (start_age < 0 || end_age < 0 || end_age > 120 || start_age > 120)
		MTHROW_AND_ERR("start_age,end_age must be initialize between 0 to 120\n");
	if (age_bin <= 0)
		MTHROW_AND_ERR("age_bin must be positive > 0\n");
	unordered_map<int, vector<const MedRegistryRecord *>> pid_to_regs, pid_to_censor;
	for (size_t i = 0; i < registry.size(); ++i)
		pid_to_regs[registry[i].pid].push_back(&registry[i]);
	if (censor_registry != NULL)
		for (size_t i = 0; i < censor_registry->size(); ++i)
			pid_to_censor[(*censor_registry)[i].pid].push_back(&(*censor_registry)[i]);

	unordered_map<int, int> pid_to_ind;
	vector<MedIdSamples> idSamples;
	vector<const MedRegistryRecord *> empty_censor;

	int conflict_count = 0, done_count = 0, skip_no_bdate = 0, no_censor = 0, example_pid = -1, no_rule = 0;
	for (auto it = pid_to_regs.begin(); it != pid_to_regs.end(); ++it) {
		vector<const MedRegistryRecord *> *all_pid_records = &it->second;
		if (pid_to_ind.find(it->first) == pid_to_ind.end()) {
			pid_to_ind[it->first] = (int)idSamples.size();
			MedIdSamples pid_sample(it->first);
			idSamples.push_back(pid_sample);
		}
		int pid_bdate = -1;
		if (pids_bdates.find(it->first) != pids_bdates.end())
			pid_bdate = pids_bdates.at(it->first);
		else {
			++skip_no_bdate;
			example_pid = it->first;
			continue;
		}
		vector<const MedRegistryRecord *> *r_censor = &empty_censor;
		if (pid_to_censor.find(it->first) != pid_to_censor.end())
			r_censor = &pid_to_censor[it->first];
		else
			++no_censor;
		if (censor_registry != NULL && pid_to_censor.find(it->first) == pid_to_censor.end())
			continue; //filter sample
		for (int age = start_age; age <= end_age; age += age_bin) {
			//search for match in all regs:
			int pred_start_date = medial::repository::DateAdd(pid_bdate, 365 * age); //mark start date in age_bin to age
			int pred_end_date = medial::repository::DateAdd(pred_start_date, 365 * age_bin); //end date in age_bin

			MedSample smp;
			smp.id = it->first;
			smp.time = medial::repository::DateAdd(pred_start_date, age_bin * 365 / 2); //choose middle
			int curr_index = 0, final_selected = -1;
			float reg_val = -1;
			int reg_time = -1;
			//run on all matches:
			while (curr_index < all_pid_records->size()) {
				if (curr_index < all_pid_records->size()) {
					if (!outcome_interaction_mode.find((*all_pid_records)[curr_index]->registry_value)) {
						++no_rule;
						if (no_rule < 5)
							MWARN("Warning: missing rule for %f - skipping!!\n", (*all_pid_records)[curr_index]->registry_value);
						++curr_index;
						continue;
					}

					if (!censor_interaction_mode.find((*all_pid_records)[curr_index]->registry_value)) {
						++no_rule;
						if (no_rule < 5)
							MWARN("Warning: missing censor rule for %f - skipping!!\n", (*all_pid_records)[curr_index]->registry_value);
						++curr_index;
						continue;
					}
				}
				if (curr_index < all_pid_records->size() &&
					!medial::process::in_time_window(pred_start_date, (*all_pid_records)[curr_index], *r_censor,
						0, 365 * age_bin, outcome_interaction_mode, censor_interaction_mode)) {
					++curr_index;
					continue;
				}
				if (curr_index >= all_pid_records->size())
					break; //skip if no match
				//found match:
				if (reg_time == -1) { //first match
					reg_val = (*all_pid_records)[curr_index]->registry_value;
					reg_time = (*all_pid_records)[curr_index]->end_date;
					final_selected = curr_index;
				}
				else if (reg_val != (*all_pid_records)[curr_index]->registry_value) {
					//if already found and conflicting:
					if (conflict_method == Drop) {
						reg_val = -1;
						reg_time = -1;
						final_selected = -1;
						break;
					}
					else if (conflict_method == Max) {
						if (reg_val < (*all_pid_records)[curr_index]->registry_value) {
							reg_val = (*all_pid_records)[curr_index]->registry_value;
							reg_time = (*all_pid_records)[curr_index]->end_date;
							final_selected = curr_index;
						}
					}
					else {
						//insert current and update next:
						smp.outcomeTime = reg_val > 0 ? (*all_pid_records)[curr_index]->start_date : reg_time;
						smp.outcome = reg_val;
						idSamples[pid_to_ind.at(it->first)].samples.push_back(smp);
						++done_count;
						reg_val = (*all_pid_records)[curr_index]->registry_value;
						reg_time = (*all_pid_records)[curr_index]->end_date;
						final_selected = curr_index;
					}
					++conflict_count;
					break;
				}

				++curr_index;
			}

			if (reg_time != -1) {
				smp.outcomeTime = reg_val > 0 ? (*all_pid_records)[final_selected]->start_date : reg_time;
				smp.outcome = reg_val;
				idSamples[pid_to_ind.at(it->first)].samples.push_back(smp);
				++done_count;
			}
		}
	}

	if (no_rule > 0)
		MLOG("WARNING MedSamplingYearly:do_sample - has %d samples with no rules for time window\n", no_rule);
	if (no_censor > 0)
		if (censor_registry != NULL)
			MLOG("WARNING MedSamplingYearly:do_sample - has %d patients with no censor dates\n", no_censor);
		else
			MLOG("WARNING MedSamplingYearly:do_sample - no censoring time region was given\n");
	if (skip_no_bdate > 0)
		MLOG("WARNING :: Skipped %d registry records because no bdate: example pid=%d\n", skip_no_bdate, example_pid);
	if (conflict_count > 0)
		MLOG("Sampled registry with %d conflicts. has %d registry records\n", conflict_count, done_count);
	//keep non empty pids:
	for (int i = (int)idSamples.size() - 1; i >= 0; --i)
		if (!idSamples[i].samples.empty())
			samples.idSamples.push_back(idSamples[i]);
	samples.sort_by_id_date();
}

int MedSamplingDates::init(map<string, string>& map) {
	conflict_method = Drop; //default
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "take_count")
			take_count = stoi(it->second);
		else if (it->first == "outcome_interaction_mode")
			medial::sampling::init_time_window_mode(it->second, outcome_interaction_mode);
		else if (it->first == "censor_interaction_mode")
			medial::sampling::init_time_window_mode(it->second, censor_interaction_mode);
		else if (it->first == "time_from")
			time_from = stoi(it->second);
		else if (it->first == "time_to")
			time_to = stoi(it->second);
		else if (it->first == "conflict_method")
			conflict_method = ConflictMode_name_to_type(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (take_count <= 0)
		MTHROW_AND_ERR("take_count must be positive > 0\n");
	return 0;
}

void MedSamplingDates::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry) {
	unordered_map<int, vector<const MedRegistryRecord *>> pid_to_regs, pid_to_censor;
	for (size_t i = 0; i < registry.size(); ++i)
		pid_to_regs[registry[i].pid].push_back(&registry[i]);
	if (censor_registry != NULL)
		for (size_t i = 0; i < censor_registry->size(); ++i)
			pid_to_censor[(*censor_registry)[i].pid].push_back(&(*censor_registry)[i]);
	vector<const MedRegistryRecord *> empty_censor;

	unordered_map<int, MedIdSamples> map_pid_samples;
	int no_censor = 0, no_rule = 0;
	for (size_t i = 0; i < samples_list_pid_dates.size(); ++i)
	{
		const vector<pair<int, int>> &all_sample_options = samples_list_pid_dates[i];
		if (all_sample_options.empty())
			continue;
		uniform_int_distribution<> current_rand(0, (int)all_sample_options.size() - 1);
		for (size_t k = 0; k < take_count; ++k)
		{
			int choosed_index = current_rand(gen);
			const pair<int, int> &choosed_option = all_sample_options[choosed_index];
			int choosed_pid = choosed_option.first, choosed_time = choosed_option.second;
			if (pid_to_regs.find(choosed_pid) == pid_to_regs.end())
				continue;
			vector<const MedRegistryRecord *> &all_pid_records = pid_to_regs.at(choosed_pid);
			MedIdSamples sample_id(choosed_pid);
			if (map_pid_samples.find(choosed_pid) == map_pid_samples.end())
				map_pid_samples[choosed_pid] = sample_id;
			//find registry match for the selected option:
			vector<const MedRegistryRecord *> *r_censor = &empty_censor;
			if (pid_to_censor.find(choosed_pid) != pid_to_censor.end())
				r_censor = &pid_to_censor[choosed_pid];
			else
				++no_censor;
			if (censor_registry != NULL && pid_to_censor.find(choosed_pid) == pid_to_censor.end())
				continue; //filter sample

			int curr_index = 0;
			float reg_val = -1;
			int reg_time = -1;
			//run on all matches:
			while (curr_index < all_pid_records.size()) {
				if (curr_index < all_pid_records.size()) {
					if (!outcome_interaction_mode.find(all_pid_records[curr_index]->registry_value)) {
						++no_rule;
						if (no_rule < 5)
							MWARN("Warning: missing rule for %f - skipping!!\n", all_pid_records[curr_index]->registry_value);
						++curr_index;
						continue;
					}

					if (!censor_interaction_mode.find(all_pid_records[curr_index]->registry_value)) {
						++no_rule;
						if (no_rule < 5)
							MWARN("Warning: missing censor rule for %f - skipping!!\n", all_pid_records[curr_index]->registry_value);
						++curr_index;
						continue;
					}
				}
				if (curr_index < all_pid_records.size() && !medial::process::in_time_window(choosed_time, all_pid_records[curr_index],
					*r_censor, time_from, time_to, outcome_interaction_mode, censor_interaction_mode)) {
					++curr_index;
					continue;
				}
				if (curr_index >= all_pid_records.size())
					break; //skip if no match
						   //found match:
				if (reg_time == -1) { //first match
					reg_val = all_pid_records[curr_index]->registry_value;
					if (reg_val <= 0)
						reg_time = all_pid_records[curr_index]->end_date; //control take end_time
					else
						reg_time = all_pid_records[curr_index]->start_date; // case take start_time
				}
				else if (reg_val != all_pid_records[curr_index]->registry_value) {
					//if already found and conflicting:
					if (conflict_method == Drop) {
						reg_val = -1;
						reg_time = -1;
						break;
					}
					else if (conflict_method == Max) {
						if (reg_val < all_pid_records[curr_index]->registry_value) {
							reg_val = all_pid_records[curr_index]->registry_value;
							reg_time = all_pid_records[curr_index]->end_date;
						}
					}
					else {
						MedSample smp;
						smp.id = choosed_pid;
						smp.time = choosed_time;

						smp.outcomeTime = reg_time;
						smp.outcome = reg_val;
						map_pid_samples[choosed_pid].samples.push_back(smp);

						reg_val = all_pid_records[curr_index]->registry_value;
						reg_time = all_pid_records[curr_index]->end_date;
					}
				}

				++curr_index;
			}
			//if found has value in reg_val, reg_time else -1
			if (reg_time != -1) {

				MedSample smp;
				smp.id = choosed_pid;
				smp.time = choosed_time;

				smp.outcomeTime = reg_time;
				smp.outcome = reg_val;
				map_pid_samples[choosed_pid].samples.push_back(smp);
			}
		}
	}

	if (no_rule > 0)
		MLOG("WARNING MedSamplingYearly:do_sample - has %d samples with no rules for time window\n", no_rule);
	if (no_censor > 0)
		if (censor_registry != NULL)
			MLOG("WARNING MedSamplingYearly:do_sample - has %d patients with no censor dates\n", no_censor);
		else
			MLOG("WARNING MedSamplingYearly:do_sample - no censoring time region was given\n");

	for (auto it = map_pid_samples.begin(); it != map_pid_samples.end(); ++it)
		if (!it->second.samples.empty())
			samples.idSamples.push_back(it->second);
	samples.sort_by_id_date();
}

MedSamplingStrategy *MedSamplingStrategy::make_sampler(const string &sampler_name) {
	MedSamplingStrategy *sampler;

	//! [MedSamplingStrategy::make_sampler]
	if (sampler_name == "time_window")
		sampler = new MedSamplingTimeWindow;
	else if (sampler_name == "yearly")
		sampler = new MedSamplingYearly;
	else if (sampler_name == "age")
		sampler = new MedSamplingAge;
	else if (sampler_name == "dates")
		sampler = new MedSamplingDates;
	else if (sampler_name == "fixed_time")
		sampler = new MedSamplingFixedTime;
	else
		MTHROW_AND_ERR("Unsupported Sampling method %s\n", sampler_name.c_str());
	//! [MedSamplingStrategy::make_sampler]

	return sampler;
}

MedSamplingStrategy *MedSamplingStrategy::make_sampler(const string &sampler_name, const string &init_params) {
	MedSamplingStrategy *sampler = make_sampler(sampler_name);
	sampler->init_from_string(init_params);
	return sampler;
}

int MedSamplingFixedTime::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "start_time")
			start_time = stoi(it->second);
		else if (it->first == "end_time")
			end_time = stoi(it->second);
		else if (it->first == "time_jump") {
			time_jump = stoi(it->second);
			if (time_jump <= 0)
				MTHROW_AND_ERR("time_jump must be positive > 0\n");
		}
		else if (it->first == "back_random_duration")
			back_random_duration = stoi(it->second);
		else if (it->first == "time_from")
			time_from = stoi(it->second);
		else if (it->first == "time_to")
			time_to = stoi(it->second);
		else if (it->first == "outcome_interaction_mode")
			medial::sampling::init_time_window_mode(it->second, outcome_interaction_mode);
		else if (it->first == "censor_interaction_mode")
			medial::sampling::init_time_window_mode(it->second, censor_interaction_mode);
		else if (it->first == "conflict_method")
			conflict_method = ConflictMode_name_to_type(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (back_random_duration < 0)
		MTHROW_AND_ERR("back_random_duration must be positive\n");

	return 0;
}

void MedSamplingFixedTime::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples, const vector<MedRegistryRecord> *censor_registry) {
	if (time_jump <= 0)
		MTHROW_AND_ERR("time_jump must be positive > 0\n");

	int random_back_dur = 1;
	bool use_random = back_random_duration > 0;
	if (use_random)
		random_back_dur = back_random_duration;

	uniform_int_distribution<> rand_int(0, random_back_dur);
	unordered_map<int, int> pid_to_ind;
	vector<MedIdSamples> idSamples;
	vector<const MedRegistryRecord *> empty_censor;

	unordered_map<int, vector<const MedRegistryRecord *>> pid_to_regs, pid_to_censor;
	for (size_t i = 0; i < registry.size(); ++i)
		pid_to_regs[registry[i].pid].push_back(&registry[i]);
	if (censor_registry != NULL)
		for (size_t i = 0; i < censor_registry->size(); ++i)
			pid_to_censor[(*censor_registry)[i].pid].push_back(&(*censor_registry)[i]);

	int conflict_count = 0, done_count = 0, no_censor = 0, no_rule = 0;
	for (auto it = pid_to_regs.begin(); it != pid_to_regs.end(); ++it)
	{
		vector<const MedRegistryRecord *> *all_pid_records = &it->second;
		vector<const MedRegistryRecord *> *r_censor = &empty_censor;
		if (pid_to_censor.find(it->first) != pid_to_censor.end())
			r_censor = &pid_to_censor[it->first];
		else
			++no_censor;
		if (censor_registry != NULL && pid_to_censor.find(it->first) == pid_to_censor.end())
			continue; //filter sample

		long start_date = start_time;
		long end_date = end_time;

		//vector<const MedRegistryRecord *> *selected_p = all_pid_records;
		vector<const MedRegistryRecord *> *selected_p = r_censor;
		if (r_censor->empty()) {
			selected_p = all_pid_records;
			++no_censor;
		}

		if (start_date == 0 && !selected_p->empty()) {
			//select min{min_allowed on all_pid_records}
			start_date = selected_p->front()->start_date;
			for (size_t i = 1; i < selected_p->size(); ++i)
				if (start_date > selected_p->at(i)->start_date)
					start_date = selected_p->at(i)->start_date;
		}
		if (end_date == 0 && !selected_p->empty()) {
			//select max{max_allowed on all_pid_records}
			end_date = selected_p->front()->end_date;
			for (size_t i = 1; i < selected_p->size(); ++i)
				if (end_date < selected_p->at(i)->end_date)
					end_date = selected_p->at(i)->end_date;
		}



		if (pid_to_ind.find(it->first) == pid_to_ind.end()) {
			pid_to_ind[it->first] = (int)idSamples.size();
			MedIdSamples pid_sample(it->first);
			idSamples.push_back(pid_sample);
		}

		for (long date = start_date; date <= end_date; date = medial::repository::DateAdd(date, time_jump)) {
			//search for match in all regs:
			int pred_date = date;
			if (use_random)
				pred_date = medial::repository::DateAdd(pred_date, -rand_int(gen));

			MedSample smp;
			smp.id = it->first;
			smp.time = pred_date;
			int curr_index = 0, final_selected = -1;
			float reg_val = -1;
			int reg_time = -1;
			//run on all matches:
			while (curr_index < all_pid_records->size()) {
				if (curr_index < all_pid_records->size()) {
					if (!outcome_interaction_mode.find((*all_pid_records)[curr_index]->registry_value)) {
						++no_rule;
						if (no_rule < 5)
							MWARN("Warning: missing rule for %f - skipping!!\n", (*all_pid_records)[curr_index]->registry_value);
						++curr_index;
						continue;
					}

					if (!censor_interaction_mode.find((*all_pid_records)[curr_index]->registry_value)) {
						++no_rule;
						if (no_rule < 5)
							MWARN("Warning: missing censor rule for %f - skipping!!\n", (*all_pid_records)[curr_index]->registry_value);
						++curr_index;
						continue;
					}
				}
				if (curr_index < all_pid_records->size() &&
					!medial::process::in_time_window(pred_date, (*all_pid_records)[curr_index], *r_censor,
						time_from, time_to, outcome_interaction_mode, censor_interaction_mode)) {
					++curr_index;
					continue;
				}
				if (curr_index >= all_pid_records->size())
					break; //skip if no match
						   //found match:
				if (reg_time == -1) { //first match
					reg_val = (*all_pid_records)[curr_index]->registry_value;
					reg_time = (*all_pid_records)[curr_index]->end_date;
					final_selected = curr_index;
				}
				else if (reg_val != (*all_pid_records)[curr_index]->registry_value) {
					//if already found and conflicting:
					if (conflict_method == Drop) {
						reg_val = -1;
						reg_time = -1;
						final_selected = -1;
						break;
					}
					else if (conflict_method == Max) {
						if (reg_val < (*all_pid_records)[curr_index]->registry_value) {
							reg_val = (*all_pid_records)[curr_index]->registry_value;
							reg_time = (*all_pid_records)[curr_index]->end_date;
							final_selected = curr_index;
						}
					}
					else {
						//insert current and update next:
						smp.outcomeTime = reg_val > 0 ? (*all_pid_records)[curr_index]->start_date : reg_time;
						smp.outcome = reg_val;
						idSamples[pid_to_ind.at(it->first)].samples.push_back(smp);
						++done_count;
						reg_val = (*all_pid_records)[curr_index]->registry_value;
						reg_time = (*all_pid_records)[curr_index]->end_date;
						final_selected = curr_index;
					}
					++conflict_count;
					break;
				}

				++curr_index;
			}

			if (reg_time != -1) {
				smp.outcomeTime = reg_val > 0 ? (*all_pid_records)[final_selected]->start_date : reg_time;
				smp.outcome = reg_val;
				idSamples[pid_to_ind.at(it->first)].samples.push_back(smp);
				++done_count;
			}
		}
	}

	if (no_rule > 0)
		MLOG("WARNING MedSamplingYearly:do_sample - has %d samples with no rules for time window\n", no_rule);
	if (no_censor > 0)
		if (censor_registry != NULL)
			MLOG("WARNING MedSamplingYearly:do_sample - has %d patients with no censor dates\n", no_censor);
		else
			MLOG("WARNING MedSamplingYearly:do_sample - no censoring time region was given\n");
	if (conflict_count > 0)
		MLOG("Sampled registry with %d conflicts. has %d registry records\n", conflict_count, done_count);
	//keep non empty pids:
	for (int i = (int)idSamples.size() - 1; i >= 0; --i)
		if (!idSamples[i].samples.empty())
			samples.idSamples.push_back(idSamples[i]);
	samples.sort_by_id_date();
}