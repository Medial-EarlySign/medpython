#include <MedUtils/MedUtils/MedSamplingStrategy.h>
#include <fstream>
#include <string>
#include <iostream>
#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string.hpp>
#include "MedSamplingHelper.h"

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
	if (minimal_times.size() != maximal_times.size())
		MTHROW_AND_ERR("Error in MedSamplingTimeWindow::init - minimal_times.size()[%zu]!=maximal_times.size()[%zu]\n",
			minimal_times.size(), maximal_times.size());
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

void MedSamplingTimeWindow::_get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
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
			for (size_t j = 0; j < it->second.size(); ++j)
			{
				int min_allowed_date = it->second[j].first;
				int max_allowed_date = it->second[j].second;

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
						MLOG("Exampled Row Skipped: pid=%d, pid_bdate=%d, min_allowed_date=%d, max_allowed_date=%d, currDate=%d, age=%d\n",
							pid, pid_bdate, min_allowed_date, max_allowed_date, currDate, (int)medial::repository::DateDiff(pid_bdate, currDate));
					}
					continue;
				}
				year_diff_to_first_pred *= 365; //now time in days - convert to time unit
				year_diff_to_first_pred = med_time_converter.convert_days(global_default_windows_time_unit, (int)year_diff_to_first_pred);

				if (diff_window > year_diff_to_first_pred) //validate we wont go back too far
					diff_window = int(year_diff_to_first_pred); //window passed max allowed - so cut in max

				int rnd_days_diff = 0;
				if (take_max || use_random) {
					if (diff_window < 1)  //not enought time to sample - skip
						continue;
					uniform_int_distribution<> rand_int(0, diff_window);
					rnd_days_diff = (int)rand_int(gen);
				}

				for (size_t k = 0; k < sample_count; ++k)
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
	int prediction_month_day = 101;
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "start_year" || it->first == "from_year") {
			start_time = stoi(it->second);
			if (start_time <= 1900 || start_time >= 2100)
				MTHROW_AND_ERR("start_year must be initialize between 1900 to 2100\n");
			start_time = start_time * 10000; //convert to DATE format
		}
		else if (it->first == "start_time")
			start_time = stoi(it->second);
		else if (it->first == "end_year" || it->first == "to_year") {
			end_time = stoi(it->second);
			if (end_time <= 1900 || end_time >= 2100)
				MTHROW_AND_ERR("end_year must be initialize between 1900 to 2100\n");
			end_time = end_time * 10000; //convert to DATE format
		}
		else if (it->first == "end_time")
			end_time = stoi(it->second);
		else if (it->first == "day_jump" || it->first == "time_jump") {
			time_jump = stoi(it->second);
			if (time_jump <= 0)
				MTHROW_AND_ERR("day_jump must be positive > 0\n");
		}
		else if (it->first == "time_jump_unit")
			time_jump_unit = med_time_converter.string_to_type(it->second);
		else if (it->first == "time_range_unit")
			time_range_unit = med_time_converter.string_to_type(it->second);
		else if (it->first == "prediction_month_day")
			prediction_month_day = stoi(it->second);
		else if (it->first == "back_random_duration")
			back_random_duration = stoi(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}

	if (prediction_month_day < 100 || prediction_month_day % 100 > 31)
		MTHROW_AND_ERR("prediction_month_day must be positive >= 100 <=1231\n");
	if (back_random_duration < 0)
		MTHROW_AND_ERR("back_random_duration must be positive\n");

	start_time += prediction_month_day;
	return 0;
}

int MedSamplingAge::init(map<string, string>& map) {
	age_bin = 1; //deafult
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "start_age")
			start_age = stoi(it->second);
		else if (it->first == "end_age")
			end_age = stoi(it->second);
		else if (it->first == "age_bin")
			age_bin = stoi(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	return 0;
}

void MedSamplingAge::_get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
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
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "take_count")
			take_count = med_stoi(it->second);
		else if (it->first == "sample_with_filters")
			sample_with_filters = med_stoi(it->second) > 0;
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (take_count <= 0)
		MTHROW_AND_ERR("take_count must be positive > 0\n");
	return 0;
}

void MedSamplingDates::_get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
	unordered_map<int, vector<int>> &pid_options) const {
	random_device rd;
	mt19937 gen(rd());
	//keep only "legal" dates by pid_time_ranges:
	vector<vector<pair<int, int>>> filtered_opts;
	filtered_opts.reserve(samples_list_pid_dates.size());
	if (sample_with_filters)
		for (size_t i = 0; i < samples_list_pid_dates.size(); ++i)
		{
			vector<pair<int, int>> filtered_grp_opts;
			for (size_t k = 0; k < samples_list_pid_dates[i].size(); ++k)
			{
				int pid = samples_list_pid_dates[i][k].first;
				int time = samples_list_pid_dates[i][k].second;

				bool is_legal = false;

				if (pid_time_ranges.find(pid) != pid_time_ranges.end()) {
					const vector<pair<int, int>> &pid_range = pid_time_ranges.at(pid);
					int pid_range_iter = 0;
					//check time is legal in pid_time_ranges:
					while (!is_legal && pid_range_iter < pid_range.size()) {
						int start_time = pid_range[pid_range_iter].first;
						int end_time = pid_range[pid_range_iter].second;
						is_legal = (time >= start_time && time <= end_time);
						++pid_range_iter;
					}
					if (is_legal)
						filtered_grp_opts.push_back(samples_list_pid_dates[i][k]);
				}
			}

			filtered_opts.push_back(filtered_grp_opts);
		}
	else
		filtered_opts = samples_list_pid_dates;

	for (size_t i = 0; i < filtered_opts.size(); ++i)
	{
		const vector<pair<int, int>> &all_sample_options = filtered_opts[i];
		if (all_sample_options.empty())
			continue;
		uniform_int_distribution<> current_rand(0, (int)all_sample_options.size() - 1);
		for (size_t k = 0; k < take_count; ++k)
		{
			int choosed_index = current_rand(gen);
			const pair<int, int> &choosed_option = all_sample_options[choosed_index];
			int choosed_pid = choosed_option.first, choosed_time = choosed_option.second;

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
	else if (sampler_name == "stick")
		sampler = new MedSamplingStick;
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
		else if (it->first == "time_jump_unit")
			time_jump_unit = med_time_converter.string_to_type(it->second);
		else if (it->first == "time_range_unit")
			time_range_unit = med_time_converter.string_to_type(it->second);
		else if (it->first == "back_random_duration")
			back_random_duration = stoi(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (back_random_duration < 0)
		MTHROW_AND_ERR("back_random_duration must be positive\n");

	return 0;
}

int MedSamplingFixedTime::add_time(int time, int add) const {
	int conv = med_time_converter.convert_times(time_range_unit, time_jump_unit, time) + add;
	return med_time_converter.convert_times(time_jump_unit, time_range_unit, conv);
}

bool validate_in_dates(const vector<pair<int, int>> &pid_dates, int pred_date, const TimeWindowMode &interaction) {
	bool inside = false;
	
	for (size_t i = 0; i < pid_dates.size() && !inside; ++i)
		inside = medial::sampling::in_time_window_simple(pred_date, 
			pid_dates[i].first, pid_dates[i].second, false, interaction);
		//inside = pred_date >= pid_dates[i].first && pred_date <= pid_dates[i].second;
	return inside;
}

void filter_opts(unordered_map<int, vector<int>> &pid_options, 
	const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges, const TimeWindowMode &interaction) {
	unordered_map<int, vector<int>> pid_filtered;
	for (auto it = pid_options.begin(); it != pid_options.end(); ++it)
	{
		int pid = it->first;
		if (pid_time_ranges.find(pid) == pid_time_ranges.end())
			continue;
		const vector<pair<int, int>> &pid_date = pid_time_ranges.at(pid);
		for (int pred_date : it->second)
			if (validate_in_dates(pid_date, pred_date, interaction))
				pid_filtered[pid].push_back(pred_date);
	}
	pid_options = move(pid_filtered);
}

void MedSamplingFixedTime::_get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
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

		for (long date = start_date; date <= end_date; date = add_time(date, time_jump)) {
			//search for match in all regs:
			int pred_date = date;
			if (use_random)
				pred_date = add_time(pred_date, -rand_int(gen));
			pid_options[pid].push_back(pred_date);
		}
	}
}

void MedSamplingStrategy::set_filters(const string &filtering_str) {
	filtering_params.init_from_string(filtering_str);
}

void MedSamplingStrategy::init_sampler(MedRepository &rep) {
	get_bdates(rep, pids_bdates);
}

template<class T> void commit_selection_vec(vector<T> &vec, const vector<int> &idx) {
	vector<T> filt(idx.size());
	for (size_t i = 0; i < idx.size(); ++i)
		filt[i] = vec[idx[i]];
	vec.swap(filt);
}

bool MedSamplingStrategy::apply_filter_params(unordered_map<int, vector<pair<int, int>>> &pid_time_ranges) const {
	if (filtering_params.max_age < 0 && filtering_params.max_time < 0 && filtering_params.min_age == 0 && filtering_params.min_time == 0)
		return false; //no filters

	for (auto it = pid_time_ranges.begin(); it != pid_time_ranges.end(); ++it)
	{
		int pid_bdate = 0; //Age need BDATE, will only be set if age filters are active
		if (filtering_params.min_age > 0 || filtering_params.max_age > 0) { //age filtering is needed, update bdate
			if (pids_bdates.find(it->first) == pids_bdates.end())
				MTHROW_AND_ERR("Error in MedSamplingStrategy::apply_filter_params - sampler isn'n initialized with init_sampler or pid %d"
					" is ilegel without bdate.\n", it->first);
			pid_bdate = pids_bdates.at(it->first);
		}
		vector<int> selected_idx;
		for (size_t i = 0; i < it->second.size(); ++i)
		{
			int start = pid_time_ranges[it->first][i].first;
			int end = pid_time_ranges[it->first][i].second;
			//apply Year filters:
			if (start < filtering_params.min_time)
				start = filtering_params.min_time;
			if (filtering_params.max_time > 0 && end > filtering_params.max_time)
				end = filtering_params.max_time;
			//apply age filters - If needed:
			if (filtering_params.min_age > 0 || filtering_params.max_age > 0) {
				float start_age = medial::repository::DateDiff(pid_bdate, start);
				float end_age = medial::repository::DateDiff(pid_bdate, end);
				if (start_age < filtering_params.min_age)
					start = medial::repository::DateAdd(pid_bdate, filtering_params.min_age * 365);
				if (filtering_params.max_age > 0 && end_age > filtering_params.max_age)
					end = medial::repository::DateAdd(pid_bdate, filtering_params.max_age * 365);
			}
			//set new range
			if (end > start) {
				pid_time_ranges[it->first][i].first = start;
				pid_time_ranges[it->first][i].second = end;
				selected_idx.push_back((int)i);
			}
		}

		commit_selection_vec(pid_time_ranges[it->first], selected_idx);
	}
	return true;
}

void MedSamplingStrategy::get_sampling_options(const unordered_map<int, vector<pair<int, int>>> &pid_time_ranges,
	unordered_map<int, vector<int>> &pid_options) const {
	//process filters
	unordered_map<int, vector<pair<int, int>>> pid_time_ranges_filtered = pid_time_ranges;
	bool has_filters = apply_filter_params(pid_time_ranges_filtered);

	_get_sampling_options(pid_time_ranges_filtered, pid_options);

	//force apply filters:
	if (has_filters)
		filter_opts(pid_options, pid_time_ranges_filtered, filtering_params.interaction_mode);
}

int MedSamplingStick::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "signal_list")
			boost::split(signal_list, it->second, boost::is_any_of(","));
		else if (it->first == "take_count")
			take_count = med_stoi(it->second);
		else if (it->first == "sample_with_filters")
			sample_with_filters = med_stoi(it->second) > 0;
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (signal_list.empty())
		MTHROW_AND_ERR("Error in MedSamplingStick::init - please provide \"signal_list\" init argument\n");
	return 0;
}

void MedSamplingStick::init_sampler(MedRepository &rep) {
	MedSamplingDates::init_sampler(rep);

	string rep_path = rep.config_fname;
	//check we have all signals otherwise load new rep (assume at least rep was initialized, will raise error in base inti_sampler otherwise):
	vector<int> sig_ids(signal_list.size());
	bool need_read = false;
	for (size_t i = 0; i < signal_list.size(); ++i) {
		sig_ids[i] = rep.sigs.sid(signal_list[i]);
		if (sig_ids[i] < 0)
			MTHROW_AND_ERR("Error in MedSamplingStick::init_sampler - Unknown signal %s in repository %s\n",
				signal_list[i].c_str(), rep_path.c_str());
		if (!rep.index.index_table[sig_ids[i]].is_loaded)
			need_read = true;
	}
	MedRepository *p_rep = &rep;
	MedRepository rep2;
	if (need_read) {
		if (rep2.read_all(rep_path, rep.pids, sig_ids) < 0)
			MTHROW_AND_ERR("Error in MedSamplingStick::init_sampler - can't read repository %s\n", rep_path.c_str());
		p_rep = &rep2;
	}

	//use p_rep to fetch signals as candidate dates for patient:
	for (size_t i = 0; i < p_rep->pids.size(); ++i)
	{
		int pid = p_rep->pids[i];
		set<int> possible_dates;
		for (size_t k = 0; k < sig_ids.size(); ++k)
		{
			UniversalSigVec usv;
			p_rep->uget(pid, sig_ids[k], usv);
			for (int l = 0; l < usv.len; ++l)
				possible_dates.insert(usv.Time(l));
		}

		if (!possible_dates.empty())
			samples_list_pid_dates.push_back(vector<pair<int, int>>());
		for (int candidate_date : possible_dates)
			samples_list_pid_dates.back().push_back(pair<int, int>(pid, candidate_date));
	}

}