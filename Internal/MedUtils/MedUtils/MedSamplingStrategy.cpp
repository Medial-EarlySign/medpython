#include <MedUtils/MedUtils/MedSamplingStrategy.h>
#include <fstream>
#include <string>
#include <iostream>
#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int MedSamplingTimeWindow::init(map<string, string>& map) {
	vector<string> tokens;
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "take_max")
			take_max = stoi(it->second) > 0;
		else if (it->first == "minimal_times") {
			boost::split(tokens, it->second, boost::is_any_of(",;"));
			minimal_times.resize(tokens.size());
			for (size_t i = 0; i < tokens.size(); ++i)
				minimal_times[i] = med_stoi(tokens[i]);
		}
		else if (it->first == "maximal_times") {
			boost::split(tokens, it->second, boost::is_any_of(",;"));
			maximal_times.resize(tokens.size());
			for (size_t i = 0; i < tokens.size(); ++i)
				maximal_times[i] = med_stoi(tokens[i]);
		}
		else if (it->first == "sample_count")
			sample_count = stoi(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (minimal_times.empty() || maximal_times.empty())
		MTHROW_AND_ERR("Error in MedSamplingTimeWindow::init - empty time windows - please provide both minimal_times,maximal_times\n");
	if (minimal_times.size() <= maximal_times.size())
		MTHROW_AND_ERR("Error in MedSamplingTimeWindow::init - minimal_times.size()!=maximal_times.size()\n");
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

void MedSamplingTimeWindow::get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
	unordered_map<int, vector<int>> &pid_options) const {
	random_device rd;
	mt19937 gen(rd());

	int skip_end_smaller_start = 0, skip_no_bdate = 0, example_pid = -1;
	for (auto it = pid_time_ranges.begin(); it != pid_time_ranges.end(); ++it)
	{
		int pid = it->first;
		int pid_bdate = -1;
		if (pids_bdates.find(pid) != pids_bdates.end())
			pid_bdate = pids_bdates.at(pid);
		else {
			++skip_no_bdate;
			example_pid = pid;
			continue;
		}

		for (size_t i = 0; i < minimal_times.size(); ++i)
		{
			for (size_t i = 0; i < it->second.size(); ++i)
			{
				int max_allowed_date = it->second[i].first;
				int min_allowed_date = it->second[i].second;

				int currDate = max_allowed_date;
				int diff_window = maximal_times[i] - minimal_times[i];
				bool use_random = !take_max && (diff_window > 1);
				currDate = medial::repository::DateAdd(currDate, -minimal_times[i]);

				float year_diff_to_first_pred;
				if (min_allowed_date <= 0) //has no limit - if "max" go back until date of birth
					year_diff_to_first_pred = medial::repository::DateDiff(pid_bdate, currDate);
				else
					year_diff_to_first_pred = medial::repository::DateDiff(min_allowed_date, currDate);
				if (year_diff_to_first_pred < 0) {
					++skip_end_smaller_start;
					if (skip_end_smaller_start < 5) {
						MLOG("Exampled Row Skipped: pid=%d, pid_bdate=%d, min_allowed_date=%d, currDate=%d, age=%d\n",
							pid, pid_bdate, min_allowed_date, currDate, (int)medial::repository::DateDiff(pid_bdate, currDate));
					}
					continue;
				}
				year_diff_to_first_pred *= 365; //now time in days - convert to time unit
				year_diff_to_first_pred = med_time_converter.convert_days(global_default_windows_time_unit, (int)year_diff_to_first_pred);

				int min_pred_date; //how many years to go back
				if (diff_window > year_diff_to_first_pred) //validate we wont go back too far
					diff_window = int(year_diff_to_first_pred); //window passed max allowed - so cut in max
				min_pred_date = medial::repository::DateAdd(currDate, -diff_window); //how many years to go back

				int rnd_days_diff = 0;
				if (take_max || use_random) {
					if (diff_window < 1)  //not enought time to sample - skip
						continue;
					uniform_int_distribution<> rand_int(0, diff_window);
					rnd_days_diff = (int)rand_int(gen);
				}

				for (size_t i = 0; i < sample_count; ++i)
				{
					int sample_pred_date = medial::repository::DateAdd(currDate, -rnd_days_diff);
					pid_options[pid].push_back(sample_pred_date);
				}
			}
		}
	}
	if (skip_no_bdate > 0)
		MWARN("Warninig MedSamplingTimeWindow::get_sampling_options - had %d pid with no bdate. example pid=%d\n", skip_no_bdate, example_pid);
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
			outcome_interaction_mode.init_from_string(it->second);
		else if (it->first == "censor_interaction_mode")
			censor_interaction_mode.init_from_string(it->second);
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

void MedSamplingYearly::get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
	unordered_map<int, vector<int>> &pid_options) const {
	if (day_jump <= 0)
		MTHROW_AND_ERR("day_jump must be positive > 0\n");
	if (end_year <= 1900 || end_year >= 2100 || start_year <= 1900 || start_year >= 2100)
		MTHROW_AND_ERR("start_year,end_year must be initialize between 1900 to 2100\n");
	random_device rd;
	mt19937 gen(rd());
	int random_back_dur = 1;
	bool use_random = back_random_duration > 0;
	if (use_random)
		random_back_dur = back_random_duration;
	uniform_int_distribution<> rand_int(0, random_back_dur);

	for (auto it = pid_time_ranges.begin(); it != pid_time_ranges.end(); ++it)
	{
		int min_date = start_year;
		int max_date = end_year;
		long start_date = min_date * 10000 + prediction_month_day;
		long end_date = max_date * 10000 + prediction_month_day;

		for (long date = start_date; date <= end_date; date = medial::repository::DateAdd(date, day_jump)) {
			//search for match in all regs:
			int pred_date = date;
			if (use_random)
				pred_date = medial::repository::DateAdd(pred_date, -rand_int(gen));
			pid_options[it->first].push_back(pred_date);
		}
	}
}

int MedSamplingAge::init(map<string, string>& map) {
	conflict_method = ConflictMode::Drop; //default
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
			outcome_interaction_mode.init_from_string(it->second);
		else if (it->first == "censor_interaction_mode")
			censor_interaction_mode.init_from_string(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	return 0;
}

void MedSamplingAge::init_sampler(MedRepository &rep) {
	get_bdates(rep, pids_bdates);
}

void MedSamplingAge::get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
	unordered_map<int, vector<int>> &pid_options) const {
	if (start_age < 0 || end_age < 0 || end_age > 120 || start_age > 120)
		MTHROW_AND_ERR("start_age,end_age must be initialize between 0 to 120\n");
	if (age_bin <= 0)
		MTHROW_AND_ERR("age_bin must be positive > 0\n");

	int skip_no_bdate = 0, example_pid = -1;
	for (auto it = pid_time_ranges.begin(); it != pid_time_ranges.end(); ++it) {
		int pid = it->first;
		int pid_bdate = -1;
		if (pids_bdates.find(pid) != pids_bdates.end())
			pid_bdate = pids_bdates.at(pid);
		else {
			++skip_no_bdate;
			example_pid = pid;
			continue;
		}
		for (int age = start_age; age <= end_age; age += age_bin) {
			//search for match in all regs:
			int pred_start_date = medial::repository::DateAdd(pid_bdate, med_time_converter.convert_days(global_default_windows_time_unit, 365 * age)); //mark start date in age_bin to age

			pid_options[pid].push_back(pred_start_date);
		}
	}

	if (skip_no_bdate > 0)
		MLOG("WARNING :: Skipped %d registry records because no bdate: example pid=%d\n", skip_no_bdate, example_pid);
}

int MedSamplingDates::init(map<string, string>& map) {
	conflict_method = ConflictMode::Drop; //default
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "take_count")
			take_count = stoi(it->second);
		else if (it->first == "outcome_interaction_mode")
			outcome_interaction_mode.init_from_string(it->second);
		else if (it->first == "censor_interaction_mode")
			censor_interaction_mode.init_from_string(it->second);
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

void MedSamplingDates::get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
	unordered_map<int, vector<int>> &pid_options) const {
	random_device rd;
	mt19937 gen(rd());

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
			if (pid_time_ranges.find(choosed_pid) == pid_time_ranges.end())
				continue;
			//TODO: add option to check that it's legal - to sample after filtering:

			pid_options[choosed_pid].push_back(choosed_time);
		}
	}
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
			outcome_interaction_mode.init_from_string(it->second);
		else if (it->first == "censor_interaction_mode")
			censor_interaction_mode.init_from_string(it->second);
		else if (it->first == "conflict_method")
			conflict_method = ConflictMode_name_to_type(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (back_random_duration < 0)
		MTHROW_AND_ERR("back_random_duration must be positive\n");

	return 0;
}

void MedSamplingFixedTime::get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
	unordered_map<int, vector<int>> &pid_options) const {
	if (time_jump <= 0)
		MTHROW_AND_ERR("time_jump must be positive > 0\n");
	random_device rd;
	mt19937 gen(rd());
	int random_back_dur = 1;
	bool use_random = back_random_duration > 0;
	if (use_random)
		random_back_dur = back_random_duration;

	uniform_int_distribution<> rand_int(0, random_back_dur);

	for (auto it = pid_time_ranges.begin(); it != pid_time_ranges.end(); ++it)
	{
		int pid = it->first;
		long start_date = start_time;
		long end_date = end_time;
		const vector<pair<int, int>> &pid_dates = it->second;

		if (start_date == 0 && !pid_dates.empty()) {
			//select min{min_allowed on all_pid_records}
			start_date = pid_dates.front().first;
			for (size_t i = 1; i < pid_dates.size(); ++i)
				if (start_date > pid_dates[i].first)
					start_date = pid_dates[i].first;
		}
		if (end_date == 0 && !pid_dates.empty()) {
			//select max{max_allowed on all_pid_records}
			end_date = pid_dates.front().second;
			for (size_t i = 1; i < pid_dates.size(); ++i)
				if (end_date < pid_dates[i].second)
					end_date = pid_dates[i].second;
		}

		for (long date = start_date; date <= end_date; date = medial::repository::DateAdd(date, time_jump)) {
			//search for match in all regs:
			int pred_date = date;
			if (use_random)
				pred_date = medial::repository::DateAdd(pred_date, -rand_int(gen));
			pid_options[pid].push_back(pred_date);
		}
	}
}