#include "MedSamplingStrategy.h"
#include <fstream>
#include <string>
#include <iostream>
#include "Logger/Logger/Logger.h"

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int MedSamplingTimeWindow::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "take_max")
			take_max = stoi(it->second) > 0;
		else if (it->first == "minimal_time_back")
			minimal_time_back = stoi(it->second);
		else if (it->first == "maximal_time_back")
			maximal_time_back = stoi(it->second);
		else if (it->first == "sample_count")
			sample_count = stoi(it->second);
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (minimal_time_back < 0 || maximal_time_back < 0)
		MTHROW_AND_ERR("negative time_back values aren't allowed\n");
	return 0;
}

void MedSamplingTimeWindow::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples) {
	int random_back_dur = 1;
	if (!take_max)
		random_back_dur = maximal_time_back;
	bool use_random = maximal_time_back > 1;

	uniform_int_distribution<> rand_int(0, random_back_dur);
	//create samples file:
	unordered_map<int, int> pid_to_ind;
	int skip_end_smaller_start = 0;
	for (MedRegistryRecord rec : registry)
	{
		bool addNew = false;
		int currDate = rec.end_date;
		if (currDate > rec.max_allowed_date)
			currDate = rec.max_allowed_date;
		float year_diff_to_first_pred;
		if (rec.min_allowed_date <= 0) //has no limit - if "max" go back until date of birth
			year_diff_to_first_pred = DateDiff(rec.start_date, currDate) + rec.age;
		else
			year_diff_to_first_pred = DateDiff(rec.min_allowed_date, currDate);
		if (year_diff_to_first_pred < 0 || rec.end_date <= rec.start_date || rec.end_date <= rec.min_allowed_date) {
			++skip_end_smaller_start;
			if (skip_end_smaller_start < 5) {
				MLOG("Exampled Row Skipped: pid=%d, reg_dates=[%d => %d], pred_dates=[%d => %d], outcome=%f, age=%2.1f\n",
					rec.pid, rec.start_date, rec.end_date, rec.min_allowed_date, rec.max_allowed_date,
					rec.registry_value, rec.age);
			}
			continue;
		}
		int min_pred_date; //how many years to go back
		if (minimal_time_back + maximal_time_back > 365 * year_diff_to_first_pred) //validate we wont go back too far
			min_pred_date = DateAdd(currDate, -int(365 * year_diff_to_first_pred));
		else
			min_pred_date = DateAdd(currDate, -minimal_time_back - maximal_time_back); //how many years to go back
		MedIdSamples patient_samples(rec.pid);
		if (pid_to_ind.find(rec.pid) == pid_to_ind.end()) {
			pid_to_ind[rec.pid] = (int)samples.idSamples.size();
			addNew = true;
		}

		int rnd_days_diff = 0;
		if (take_max) {
			float curr_year_diff;
			if (rec.min_allowed_date <= 0) //has no limit - if "max" go back until date of birth
				curr_year_diff = DateDiff(rec.start_date, currDate) + rec.age;
			else
				curr_year_diff = DateDiff(rec.min_allowed_date, currDate);
			int max_diff = int(365 * curr_year_diff) - 1;
			if (max_diff < 2) {
				if (addNew && patient_samples.samples.empty()) //was new and haven't been added yet
					pid_to_ind.erase(rec.pid);
				break;
			}
			uniform_int_distribution<> rand_int_p = uniform_int_distribution<>(0, max_diff);
			rnd_days_diff = (int)rand_int_p(gen);
		}
		else
			if (use_random)
				rnd_days_diff = (int)rand_int(gen);

		for (size_t i = 0; i < sample_count; ++i)
		{
			int sample_pred_date = DateAdd(currDate, -rnd_days_diff);
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

	if (skip_end_smaller_start > 0)
		MLOG("WARNING :: Skipped %d registry records because end_date<start_date\n", skip_end_smaller_start);
}

int MedSamplingYearly::init(map<string, string>& map) {
	conflict_method = "drop"; //default
	prediction_month_day = 101; //deafult
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
		else if (it->first == "use_allowed")
			use_allowed = stoi(it->second) > 0;
		else if (it->first == "conflict_method")
			conflict_method = it->second;
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}

	if (!(conflict_method == "drop" || conflict_method == "max" ||
		conflict_method == "all"))
		MTHROW_AND_ERR("Unsuported conflcit method - please choose: drop,all,max\n");
	if (prediction_month_day < 100 || prediction_month_day % 100 > 31)
		MTHROW_AND_ERR("prediction_month_day must be positive >= 100 <=1231\n");
	if (back_random_duration < 0)
		MTHROW_AND_ERR("back_random_duration must be positive\n");
	return 0;
}

void MedSamplingYearly::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples) {
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

	unordered_map<int, vector<const MedRegistryRecord *>> pid_to_regs;
	for (size_t i = 0; i < registry.size(); ++i)
		pid_to_regs[registry[i].pid].push_back(&registry[i]);

	int conflict_count = 0, done_count = 0;
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

		for (long date = start_date; date <= end_date; date = DateAdd(date, day_jump)) {
			//search for match in all regs:
			int pred_date = date;
			if (use_random)
				pred_date = DateAdd(pred_date, -rand_int(gen));

			MedSample smp;
			smp.id = it->first;
			smp.time = pred_date;
			int curr_index = 0, final_selected = -1;
			float reg_val = -1;
			int reg_time = -1;
			//run on all matches:
			while (curr_index < all_pid_records->size()) {
				if (!use_allowed)
					while (curr_index < all_pid_records->size() &&
						(pred_date < (*all_pid_records)[curr_index]->start_date ||
							pred_date >(*all_pid_records)[curr_index]->end_date))
						++curr_index;
				else
					while (curr_index < all_pid_records->size() &&
						(pred_date < (*all_pid_records)[curr_index]->min_allowed_date ||
							pred_date >(*all_pid_records)[curr_index]->max_allowed_date))
						++curr_index;
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
					if (conflict_method == "drop") {
						reg_val = -1;
						reg_time = -1;
						final_selected = -1;
						break;
					}
					else if (conflict_method == "max") {
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

	if (conflict_count > 0)
		MLOG("Sampled registry with %d conflicts. has %d registry records\n", conflict_count, done_count);
	//keep non empty pids:
	for (int i = (int)idSamples.size() - 1; i >= 0; --i)
		if (!idSamples[i].samples.empty())
			samples.idSamples.push_back(idSamples[i]);
	samples.sort_by_id_date();
}

int MedSamplingAge::init(map<string, string>& map) {
	conflict_method = "drop"; //default
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
			conflict_method = it->second;
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (!(conflict_method == "drop" || conflict_method == "max" ||
		conflict_method == "all"))
		MTHROW_AND_ERR("Unsuported conflcit method - please choose: drop,all,max\n");
	return 0;
}

void MedSamplingAge::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples) {
	if (start_age < 0 || end_age < 0 || end_age > 120 || start_age > 120)
		MTHROW_AND_ERR("start_age,end_age must be initialize between 0 to 120\n");
	if (age_bin <= 0)
		MTHROW_AND_ERR("age_bin must be positive > 0\n");
	unordered_map<int, vector<const MedRegistryRecord *>> pid_to_regs;
	for (size_t i = 0; i < registry.size(); ++i)
		pid_to_regs[registry[i].pid].push_back(&registry[i]);

	unordered_map<int, int> pid_to_ind;
	vector<MedIdSamples> idSamples;

	int conflict_count = 0, done_count = 0;
	for (auto it = pid_to_regs.begin(); it != pid_to_regs.end(); ++it) {
		vector<const MedRegistryRecord *> *all_pid_records = &it->second;
		if (pid_to_ind.find(it->first) == pid_to_ind.end()) {
			pid_to_ind[it->first] = (int)idSamples.size();
			MedIdSamples pid_sample(it->first);
			idSamples.push_back(pid_sample);
		}
		for (int age = start_age; age <= end_age; age += age_bin) {
			//search for match in all regs:
			int pred_start_date = DateAdd(all_pid_records->front()->start_date, -365 * (all_pid_records->front()->age - age)); //mark start date in age_bin to age
			int pred_end_date = DateAdd(pred_start_date, 365 * age_bin); //end date in age_bin

			MedSample smp;
			smp.id = it->first;
			smp.time = DateAdd(pred_start_date, age_bin * 365 / 2); //choose middle
			int curr_index = 0;
			float reg_val = -1;
			int reg_time = -1;
			//run on all matches:
			while (curr_index < all_pid_records->size()) {
				while (curr_index < all_pid_records->size() &&
					(pred_end_date < (*all_pid_records)[curr_index]->start_date ||
						pred_start_date >(*all_pid_records)[curr_index]->end_date))
					++curr_index;
				if (curr_index >= all_pid_records->size())
					break; //skip if no match
						   //found match:
				if (reg_time == -1) { //first match
					reg_val = (*all_pid_records)[curr_index]->registry_value;
					reg_time = (*all_pid_records)[curr_index]->end_date;
				}
				else if (reg_val != (*all_pid_records)[curr_index]->registry_value) {
					//if already found and conflicting:
					if (conflict_method == "drop") {
						reg_val = -1;
						reg_time = -1;
						break;
					}
					else if (conflict_method == "max") {
						if (reg_val < (*all_pid_records)[curr_index]->registry_value) {
							reg_val = (*all_pid_records)[curr_index]->registry_value;
							reg_time = (*all_pid_records)[curr_index]->end_date;
						}
					}
					else {
						//insert current and update next:
						smp.outcomeTime = reg_time;
						smp.outcome = reg_val;
						idSamples[pid_to_ind.at(it->first)].samples.push_back(smp);
						++done_count;
						reg_val = (*all_pid_records)[curr_index]->registry_value;
						reg_time = (*all_pid_records)[curr_index]->end_date;
					}
					++conflict_count;
					break;
				}

				++curr_index;
			}

			if (reg_time != -1) {
				smp.outcomeTime = reg_time;
				smp.outcome = reg_val;
				idSamples[pid_to_ind.at(it->first)].samples.push_back(smp);
				++done_count;
			}
		}
	}

	if (conflict_count > 0)
		MLOG("Sampled registry with %d conflicts. has %d registry records\n", conflict_count, done_count);
	//keep non empty pids:
	for (int i = (int)idSamples.size() - 1; i >= 0; --i)
		if (!idSamples[i].samples.empty())
			samples.idSamples.push_back(idSamples[i]);
	samples.sort_by_id_date();
}

int MedSamplingDates::init(map<string, string>& map) {
	conflict_method = "drop"; //default
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "take_count")
			take_count = stoi(it->second);
		else if (it->first == "use_allowed")
			use_allowed = stoi(it->second) > 0;
		else if (it->first == "conflict_method")
			conflict_method = it->second;
		else
			MTHROW_AND_ERR("Unsupported parameter %s for Sampler\n", it->first.c_str());
	}
	if (take_count <= 0)
		MTHROW_AND_ERR("take_count must be positive > 0\n");
	if (!(conflict_method == "drop" || conflict_method == "max" ||
		conflict_method == "all"))
		MTHROW_AND_ERR("Unsuported conflcit method - please choose: drop,all,max\n");
	return 0;
}

void MedSamplingDates::do_sample(const vector<MedRegistryRecord> &registry, MedSamples &samples) {
	unordered_map<int, vector<const MedRegistryRecord *>> pid_to_regs;
	for (size_t i = 0; i < registry.size(); ++i)
		pid_to_regs[registry[i].pid].push_back(&registry[i]);

	unordered_map<int, MedIdSamples> map_pid_samples;
	for (size_t i = 0; i < samples_list_pid_dates.size(); ++i)
	{
		const vector<pair<int, int>> &all_sample_options = samples_list_pid_dates[i];
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

			int curr_index = 0;
			float reg_val = -1;
			int reg_time = -1;
			//run on all matches:
			while (curr_index < all_pid_records.size()) {
				if (!use_allowed)
					while (curr_index < all_pid_records.size() &&
						(choosed_time < all_pid_records[curr_index]->start_date ||
							choosed_time >all_pid_records[curr_index]->end_date))
						++curr_index;
				else
					while (curr_index < all_pid_records.size() &&
						(choosed_time < all_pid_records[curr_index]->min_allowed_date ||
							choosed_time >all_pid_records[curr_index]->max_allowed_date))
						++curr_index;
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
					if (conflict_method == "drop") {
						reg_val = -1;
						reg_time = -1;
						break;
					}
					else if (conflict_method == "max") {
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

	for (auto it = map_pid_samples.begin(); it != map_pid_samples.end(); ++it)
		samples.idSamples.push_back(it->second);
	samples.sort_by_id_date();
}

MedSamplingStrategy *MedSamplingStrategy::make_sampler(const string &sampler_name) {
	MedSamplingStrategy *sampler;

	if (sampler_name == "time_window")
		sampler = new MedSamplingTimeWindow;
	else if (sampler_name == "yearly")
		sampler = new MedSamplingYearly;
	else if (sampler_name == "age")
		sampler = new MedSamplingAge;
	else if (sampler_name == "dates")
		sampler = new MedSamplingDates;
	else
		MTHROW_AND_ERR("Unsupported Sampling method %s\n", sampler_name.c_str());

	return sampler;
}

MedSamplingStrategy *MedSamplingStrategy::make_sampler(const string &sampler_name, const string &init_params) {
	MedSamplingStrategy *sampler = make_sampler(sampler_name);
	sampler->init_from_string(init_params);
	return sampler;
}