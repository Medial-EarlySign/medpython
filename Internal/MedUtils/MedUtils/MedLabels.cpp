#include "MedLabels.h"
#include <algorithm>
#include <Logger/Logger/Logger.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedUtils/MedUtils/MedUtils.h>

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

MedLabels::MedLabels(const LabelParams &params) {
	labeling_params = params;
}

MedLabels::MedLabels(const string &labling_params) {
	labeling_params.init_from_string(labling_params);
}

void MedLabels::prepare_from_registry(const vector<MedRegistryRecord> &reg_records, const vector<MedRegistryRecord> *censor_records) {
	all_reg_records = reg_records;
	for (size_t i = 0; i < all_reg_records.size(); ++i)
		pid_reg_records[all_reg_records[i].pid].push_back(&all_reg_records[i]);

	if (censor_records != NULL)
		all_censor_records = *censor_records;
	for (size_t i = 0; i < all_censor_records.size(); ++i)
		pid_censor_records[all_censor_records[i].pid].push_back(&all_censor_records[i]);
}

bool MedLabels::has_censor_reg() const {
	return !all_censor_records.empty();
}

bool MedLabels::has_censor_reg(int pid) const {
	return pid_censor_records.find(pid) != pid_censor_records.end();
}

void MedLabels::get_pids(vector<int> &pids) const {
	pids.reserve(pid_reg_records.size());
	for (auto it = pid_reg_records.begin(); it != pid_reg_records.end(); ++it)
		pids.push_back(it->first);
}

SamplingRes MedLabels::get_samples(int pid, int time, vector<MedSample> &samples, bool show_conflicts) const {
	if (pid_reg_records.empty())
		MTHROW_AND_ERR("Error in MedLabels::get_samples - please init MedLabels by calling prepare_from_registry\n");
	//search where time falls inside records - assume no conflicts (not checking for this)
	vector<const MedRegistryRecord *> empty_censor;
	const vector<const MedRegistryRecord *> *censor_p = &empty_censor;
	SamplingRes r;
	if (pid_censor_records.find(pid) != pid_censor_records.end())
		censor_p = &pid_censor_records.at(pid);
	if (pid_reg_records.find(pid) != pid_reg_records.end()) {
		const vector<const MedRegistryRecord *> &pid_recs = pid_reg_records.at(pid);

		medial::sampling::get_label_for_sample(time, pid_recs, *censor_p, labeling_params.time_from, labeling_params.time_to,
			labeling_params.censor_time_from, labeling_params.time_to, labeling_params.label_interaction_mode,
			labeling_params.censor_interaction_mode, labeling_params.conflict_method, samples, r.no_rule_cnt, r.conflict_cnt,
			r.done_cnt, false, show_conflicts);
	}
	else
		++r.miss_pid_in_reg_cnt;

	return r;
}

SamplingRes MedLabels::get_samples(int pid, const vector<int> &times, vector<MedSample> &samples, bool show_conflicts) const {
	if (pid_reg_records.empty())
		MTHROW_AND_ERR("Error in MedLabels::get_samples - please init MedLabels by calling prepare_from_registry\n");
	vector<const MedRegistryRecord *> empty_censor;
	const vector<const MedRegistryRecord *> *censor_p = &empty_censor;
	if (pid_censor_records.find(pid) != pid_censor_records.end())
		censor_p = &pid_censor_records.at(pid);
	SamplingRes r;
	if (pid_reg_records.find(pid) != pid_reg_records.end()) {
		const vector<const MedRegistryRecord *> &pid_recs = pid_reg_records.at(pid);
		for (size_t i = 0; i < times.size(); ++i) {
			medial::sampling::get_label_for_sample(times[i], pid_recs, *censor_p, labeling_params.time_from, labeling_params.time_to,
				labeling_params.censor_time_from, labeling_params.censor_time_to, labeling_params.label_interaction_mode,
				labeling_params.censor_interaction_mode, labeling_params.conflict_method, samples, r.no_rule_cnt, r.conflict_cnt,
				r.done_cnt, false, show_conflicts);
		}

	}
	else
		++r.miss_pid_in_reg_cnt;
	return r;
}

void MedLabels::get_records(int pid, vector<const MedRegistryRecord *> &reg_records, vector<const MedRegistryRecord *> &censor_records) const {
	if (pid_reg_records.find(pid) != pid_reg_records.end())
		reg_records = pid_censor_records.at(pid);
	if (pid_censor_records.find(pid) != pid_censor_records.end())
		censor_records = pid_censor_records.at(pid);
}

void update_loop(int pos, int ageBin_index, float ageBin, int signal_val,
	map<float, map<float, vector<int>>> &signalToStats, vector<unordered_map<float, int>> &val_seen_pid_pos) {
	if ((pos == 3 && val_seen_pid_pos[ageBin_index][signal_val] / 2 > 0) ||
		(pos == 2 && val_seen_pid_pos[ageBin_index][signal_val] % 2 > 0)) {
		return; //continue;
	}
	val_seen_pid_pos[ageBin_index][signal_val] += pos - 1;
	//update cnts:
#pragma omp critical 
	{
		vector<int> *cnts = &(signalToStats[signal_val][ageBin]);
		if (cnts->empty())
			cnts->resize(4); // first time
		++(*cnts)[pos];
	}
}

static void _get_parents(int codeGroup, vector<int> &parents, bool has_regex, const regex &reg_pat,
	int max_depth, int max_parents, const map<int, vector<int>> &_member2Sets, const map<int, vector<string>> &categoryId_to_name) {
	vector<int> last_parents = { codeGroup };
	if (last_parents.front() < 0)
		return; //no parents
	parents = {};

	for (size_t k = 0; k < max_depth; ++k) {
		vector<int> new_layer;
		for (int par : last_parents)
			if (_member2Sets.find(par) != _member2Sets.end()) {
				new_layer.insert(new_layer.end(), _member2Sets.at(par).begin(), _member2Sets.at(par).end());
				parents.insert(parents.end(), _member2Sets.at(par).begin(), _member2Sets.at(par).end()); //aggregate all parents
			}
		if (parents.size() >= max_parents)
			break;
		new_layer.swap(last_parents);
		if (last_parents.empty())
			break; //no more parents to loop up
	}

	if (has_regex) {
		vector<int> filtered_p;
		filtered_p.reserve(parents.size());
		for (int code : parents)
		{
			if (categoryId_to_name.find(code) == categoryId_to_name.end())
				MTHROW_AND_ERR("CategoryDependencyGenerator::post_learn_from_samples - code %d wasn't found in dict\n", code);
			const vector<string> &names = categoryId_to_name.at(code);
			int nm_idx = 0;
			bool pass_regex_filter = false;
			while (!pass_regex_filter && nm_idx < names.size())
			{
				pass_regex_filter = regex_match(names[nm_idx], reg_pat);
				++nm_idx;
			}
			if (pass_regex_filter)
				filtered_p.push_back(code);
		}
		parents.swap(filtered_p);
	}
}

static void propogate_hir(map<float, map<float, vector<int>>> &categoryVal_to_stats, bool has_regex, const regex &reg_pat
	, const map<int, vector<int>> &_member2Sets, const map<int, vector<string>> &categoryId_to_name) {
	for (auto it = categoryVal_to_stats.begin(); it != categoryVal_to_stats.end(); ++it) {
		int base_code = it->first;
		vector<int> all_parents;
		_get_parents(base_code, all_parents, has_regex, reg_pat, 50, 500, _member2Sets, categoryId_to_name);

		const map<float, vector<int>> &base_code_stats = categoryVal_to_stats.at(base_code);
		for (int code : all_parents)
		{
			//process request for code - aggregate stats from [2+0] [2+1] to code:
			map<float, vector<int>> &code_stats = categoryVal_to_stats[code]; //age, 4 counts per state
			for (auto jt = base_code_stats.begin(); jt != base_code_stats.end(); ++jt) {
				if (code_stats[jt->first].empty())
					code_stats[jt->first].resize(4);
				code_stats[jt->first][2 + 0] += jt->second[2 + 0];
				code_stats[jt->first][2 + 1] += jt->second[2 + 1];
			}
		}

	}
}

static void filter_regex_hir(map<float, map<float, vector<int>>> &categoryVal_to_stats, const regex &reg_pat
	, const map<int, vector<string>> &categoryId_to_name) {
	for (auto it = categoryVal_to_stats.begin(); it != categoryVal_to_stats.end();) {
		int base_code = (int)it->first;
		bool found_match = false;
		const vector<string> &names = categoryId_to_name.at(base_code);
		int pos_i = 0;
		while (pos_i < names.size() && !found_match) {
			found_match = regex_match(names[pos_i], reg_pat);
			++pos_i;
		}
		if (found_match)
			++it;
		else
			it = categoryVal_to_stats.erase(it);
	}
}

void MedLabels::calc_signal_stats(const string &repository_path, const string &signal_name,
	const string &signalHirerchyType, int ageBinValue, MedSamplingStrategy &sampler, const LabelParams &inc_labeling_params,
	map<float, map<float, vector<int>>> &maleSignalToStats,
	map<float, map<float, vector<int>>> &femaleSignalToStats,
	const string &debug_file, const unordered_set<int> &debug_vals) const {
	int sig_val_channel = 0;
	int sig_time_channel = 0;
	regex reg_pat;
	if (!signalHirerchyType.empty())
		reg_pat = regex(signalHirerchyType);
	bool using_regex_filter = !signalHirerchyType.empty() && signalHirerchyType != "None";
	MedRepository dataManager;
	time_t start = time(NULL);
	int duration;

	int time_window_to = labeling_params.time_to;
	int time_window_from = labeling_params.time_from;
	if (time_window_from > time_window_to)
		MTHROW_AND_ERR("Error in MedLabels::calc_signal_stats - you gave time window params in wrong order [%d, %d]\n"
			, time_window_from, time_window_to);

	vector<int> pids;
	get_pids(pids);

	vector<string> readSignals = { "GENDER" , "BDATE" };
	readSignals.push_back(signal_name);

	MLOG("Fetching signal %s using repository %s\n", signal_name.c_str(), repository_path.c_str());
	if (dataManager.read_all(repository_path, pids, readSignals) < 0)
		MTHROW_AND_ERR("error reading from repository %s\n", repository_path.c_str());
	int genderCode = dataManager.sigs.sid("GENDER");
	int bdateCode = dataManager.sigs.sid("BDATE");
	int signalCode = dataManager.sigs.sid(signal_name);
	if (genderCode < 0 || bdateCode < 0 || signalCode < 0)
		MTHROW_AND_ERR("Error in MedLabels::calc_signal_stats - can't find on of signals: GENDER,BDATE,%s in repository %s\n",
			signal_name.c_str(), repository_path.c_str());
	int sectionId = 0;
	if (dataManager.dict.SectionName2Id.find(signal_name) == dataManager.dict.SectionName2Id.end())
		MTHROW_AND_ERR("Error in MedLabels::calc_signal_stats - signal %s has no section, not categorical?\n",
			signal_name.c_str());
	sectionId = dataManager.dict.SectionName2Id.at(signal_name);
	vector<int> &all_pids = dataManager.pids;

	MLOG("Sampling for incidence stats...\n");
	MedSamples incidence_samples;
	sampler.init_sampler(dataManager);
	MedLabels inc_labeler(inc_labeling_params);
	inc_labeler.prepare_from_registry(all_reg_records, &all_censor_records);
	inc_labeler.create_samples(&sampler, incidence_samples);
	duration = (int)difftime(time(NULL), start);
	MLOG("Done in %d seconds with %zu patient ids!\n", duration, incidence_samples.idSamples.size());

	start = time(NULL);

	unordered_map<float, vector<int>> male_total_prevalence; //key=age
	unordered_map<float, vector<int>> female_total_prevalence; //key=age
	vector<unordered_map<float, unordered_set<int>>> male_pid_seen(2);
	vector<unordered_map<float, unordered_set<int>>> female_pid_seen(2);
	int unknown_gender = 0, min_age = 200, max_age = 0;
	for (MedIdSamples idSample : incidence_samples.idSamples)
		for (MedSample rec : idSample.samples)
		{
			int ind = rec.outcome > 0;
			int gend = medial::repository::get_value(dataManager, rec.id, genderCode);
			int bdate = medial::repository::get_value(dataManager, rec.id, bdateCode);
			if (gend == -1) {
				++unknown_gender;
				continue;
			}
			double curr_age = medial::repository::DateDiff(bdate, rec.time);

			float ageBin = float(ageBinValue * floor(curr_age / ageBinValue));
			if (gend == GENDER_MALE) {
				if (male_pid_seen[ind][ageBin].find(rec.id) == male_pid_seen[ind][ageBin].end()) {
					male_pid_seen[ind][ageBin].insert(rec.id);
					if (male_total_prevalence[ageBin].size() == 0)
						male_total_prevalence[ageBin].resize(2);
					++male_total_prevalence[ageBin][ind];
				}
			}
			else {
				if (female_pid_seen[ind][ageBin].find(rec.id) == female_pid_seen[ind][ageBin].end()) {
					female_pid_seen[ind][ageBin].insert(rec.id);
					if (female_total_prevalence[ageBin].size() == 0)
						female_total_prevalence[ageBin].resize(2);
					++female_total_prevalence[ageBin][ind];
				}
			}
			if (ageBin < min_age)
				min_age = (int)ageBin;
			if (ageBin > max_age)
				max_age = (int)ageBin;
		}

	if (unknown_gender > 0)
		MWARN("Has %d Unknown genders.\n", unknown_gender);
	if (!debug_file.empty()) {
		ofstream dbg_file_totals;
		dbg_file_totals.open(debug_file + ".totals");
		if (!dbg_file_totals.good())
			MTHROW_AND_ERR("IOError: Can't open debug file %s to write",
			(debug_file + ".totals").c_str());
		dbg_file_totals << "PID" << "\t" << "gender" << "\t" << "age_bin" <<
			"\t" << "registry_value" << endl;
		for (size_t i = 0; i < male_pid_seen.size(); ++i)
			for (auto it = male_pid_seen[i].begin(); it != male_pid_seen[i].end(); ++it)
				for (int pid_in : it->second)
					dbg_file_totals << pid_in << "\t" << GENDER_MALE
					<< "\t" << int(it->first) << "\t" << i << "\n";
		for (size_t i = 0; i < female_pid_seen.size(); ++i)
			for (auto it = female_pid_seen[i].begin(); it != female_pid_seen[i].end(); ++it)
				for (int pid_in : it->second)
					dbg_file_totals << pid_in << "\t" << GENDER_FEMALE
					<< "\t" << int(it->first) << "\t" << i << "\n";
		dbg_file_totals.close();
	}

	ofstream dbg_file;
	if (!debug_file.empty()) {
		dbg_file.open(debug_file);
		if (!dbg_file.good())
			MTHROW_AND_ERR("IOError: Cann't open debug file %s to write", debug_file.c_str());
		dbg_file << "PID" << "\t" << "signal_date" << "\t" << "signal_value" <<
			"\t" << "registry_start_date" << "\t" << "registry_end_date"
			<< "\t" << "gender" << "\t" << "age_bin" << "\t" << "registry_value" << endl;
	}

	duration = (int)difftime(time(NULL), start);
	MLOG("Done prep registry in %d seconds. min_age=%d, max_age=%d\n", duration, min_age, max_age);
	start = time(NULL);

	int age_bin_count = (max_age - min_age) / ageBinValue + 1;
	time_t last_time_print = start;
	int prog_pid = 0, no_rule = 0, conflict_count = 0;
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < all_pids.size(); ++i) {
		int pid = all_pids[i];
		if (has_censor_reg() && !has_censor_reg(pid))
			continue;
		//calcs on the fly pid records:
		int gender = medial::repository::get_value(dataManager, pid, genderCode);
		int BDate = medial::repository::get_value(dataManager, pid, bdateCode);
		UniversalSigVec patientFile;
		dataManager.uget(pid, signalCode, patientFile);

		vector<unordered_map<float, int>> val_seen_pid_pos(age_bin_count); //for age bin index and value (it's for same pid so gender doesnt change) - if i saw the value already
		for (int j = 0; j < patientFile.len; ++j)
		{
			int signal_val = (int)patientFile.Val(j, sig_val_channel);
			int signal_time = patientFile.Time(j, sig_time_channel);
			if (signal_val <= 0)
				continue;

			vector<int> all_vals;
			_get_parents(signal_val, all_vals, using_regex_filter, reg_pat, 50, 500,
				dataManager.dict.dict(sectionId)->Member2Sets, dataManager.dict.dict(sectionId)->Id2Names);
			all_vals.push_back(signal_val);

			int pos;
			vector<int> cnts;
			float ageBin;
			int ageBin_index;

			ageBin = float(ageBinValue * floor(double(medial::repository::DateDiff(BDate, signal_time)) / ageBinValue));
			ageBin_index = int((ageBin - min_age) / ageBinValue);
			if (ageBin < min_age || ageBin > max_age)
				continue; //skip out of range...

			vector<MedSample> found_samples;
			get_samples(pid, signal_time, found_samples);
			for (const MedSample &smp : found_samples)
			{
				pos = 2;
				//pos += 1; //registry_value > 0 - otherwise skip this
				pos += int(smp.outcome > 0);
				for (int v : all_vals)
				{
					if (gender == GENDER_MALE)
						update_loop(pos, ageBin_index, ageBin, v, maleSignalToStats, val_seen_pid_pos);
					else
						update_loop(pos, ageBin_index, ageBin, v, femaleSignalToStats, val_seen_pid_pos);
				}

				if (!debug_file.empty() && debug_vals.find(signal_val) != debug_vals.end()) {
#pragma omp critical
					dbg_file << pid << "\t" << signal_time << "\t" << signal_val
						<< "\t" << smp.time << "\t" << smp.outcomeTime
						<< "\t" << gender << "\t" << ageBin << "\t" << smp.outcome
						<< "\n";
				}
			}
		}

#pragma omp atomic
		++prog_pid;

		if (prog_pid % 10000 == 0 && (int)difftime(time(NULL), last_time_print) >= 60) {
			last_time_print = time(NULL);
			float time_elapsed = (float)difftime(time(NULL), start);
			float estimate_time = float(all_pids.size() - prog_pid) / prog_pid * time_elapsed;
			cout << "Processed " << prog_pid << " out of " << all_pids.size() << "(" << round(10000.0*(prog_pid / float(all_pids.size()))) / 100.0
				<< "%) time elapsed: " << round(time_elapsed / 6) / 10 << " Minutes, estimate time to finish " << round(10 * estimate_time / 60.0) / 10 << " Minutes"
				<< endl;
		}
	}

	if (no_rule > 0)
		MWARN("Warning has %d records with no rules for labels\n", no_rule);
	if (conflict_count > 0)
		MWARN("has %d records with conflicts\n", conflict_count);
	if (!debug_file.empty())
		dbg_file.close();
	unordered_set<float> vals;
	for (auto it = maleSignalToStats.begin(); it != maleSignalToStats.end(); ++it)
		vals.insert(it->first);
	for (auto it = femaleSignalToStats.begin(); it != femaleSignalToStats.end(); ++it)
		vals.insert(it->first);

	//filter regex:
	if (using_regex_filter) {
		filter_regex_hir(maleSignalToStats, reg_pat, dataManager.dict.dict(sectionId)->Id2Names);
		filter_regex_hir(femaleSignalToStats, reg_pat, dataManager.dict.dict(sectionId)->Id2Names);
	}

	//update values prevalence
	int warn_cnt = 0; int max_warns = 5, tot_problems = 0;
	for (auto it = vals.begin(); it != vals.end(); ++it)
	{
		if (maleSignalToStats.find(*it) != maleSignalToStats.end())
			for (auto jt = maleSignalToStats[*it].begin(); jt != maleSignalToStats[*it].end(); ) {
				if (male_total_prevalence.find(jt->first) == male_total_prevalence.end()) {
					++tot_problems;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedLabels::calc_signal_stats - Sample is too small, no incidences for age_bin=%d in males (value was=%f, cnts=[%d, %d])\n"
							, int(jt->first), *it, maleSignalToStats[*it][jt->first][2], maleSignalToStats[*it][jt->first][3]);
					}
					jt = maleSignalToStats[*it].erase(jt);
					continue;
				}
				maleSignalToStats[*it][jt->first][0] = male_total_prevalence[jt->first][0] - maleSignalToStats[*it][jt->first][2];
				maleSignalToStats[*it][jt->first][1] = male_total_prevalence[jt->first][1] - maleSignalToStats[*it][jt->first][3];
				if (maleSignalToStats[*it][jt->first][0] < 0) {
					maleSignalToStats[*it][jt->first][0] = 0;
					++tot_problems;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedLabels::calc_signal_stats - Control Male age_bin=%d, signal_value=%f, total=%d, signal=%d\n",
							int(jt->first), *it, male_total_prevalence[jt->first][0], maleSignalToStats[*it][jt->first][2]);
					}
				}
				if (maleSignalToStats[*it][jt->first][1] < 0) {
					maleSignalToStats[*it][jt->first][1] = 0;
					++tot_problems;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedLabels::calc_signal_stats - Cases Male age_bin=%d, signal_value=%f, total=%d, signal=%d\n",
							int(jt->first), *it, male_total_prevalence[jt->first][1], maleSignalToStats[*it][jt->first][3]);
					}
				}
				++jt;
			}
		if (femaleSignalToStats.find(*it) != femaleSignalToStats.end())
			for (auto jt = femaleSignalToStats[*it].begin(); jt != femaleSignalToStats[*it].end();) {
				if (female_total_prevalence.find(jt->first) == female_total_prevalence.end()) {
					++tot_problems;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedLabels::calc_signal_stats - Sample is too small, no incidences for age_bin=%d in females (value was=%f, cnts=[%d, %d])\n"
							, int(jt->first), *it, femaleSignalToStats[*it][jt->first][2], femaleSignalToStats[*it][jt->first][3]);
					}
					jt = femaleSignalToStats[*it].erase(jt);
					continue;
				}
				femaleSignalToStats[*it][jt->first][0] = female_total_prevalence[jt->first][0] - femaleSignalToStats[*it][jt->first][2];
				femaleSignalToStats[*it][jt->first][1] = female_total_prevalence[jt->first][1] - femaleSignalToStats[*it][jt->first][3];
				if (femaleSignalToStats[*it][jt->first][0] < 0) {
					femaleSignalToStats[*it][jt->first][0] = 0;
					++tot_problems;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedLabels::calc_signal_stats - Control Female age_bin=%d, signal_value=%f, total=%d, signal=%d\n",
							int(jt->first), *it, female_total_prevalence[jt->first][0], femaleSignalToStats[*it][jt->first][2]);
					}
				}
				if (femaleSignalToStats[*it][jt->first][1] < 0) {
					femaleSignalToStats[*it][jt->first][1] = 0;
					++tot_problems;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedLabels::calc_signal_stats - Cases Female age_bin=%d, signal_value=%f, total=%d, signal=%d\n",
							int(jt->first), *it, female_total_prevalence[jt->first][1], femaleSignalToStats[*it][jt->first][3]);
					}
				}
				++jt;
			}
	}
	if (tot_problems > 0)
		MWARN("Warning: MedLabels::calc_signal_stats - total miss matches: %d\n", tot_problems);

	duration = (int)difftime(time(NULL), start);
	MLOG("Finished in %d seconds with %d records in males and %d records in females\n",
		duration, (int)maleSignalToStats.size(), (int)femaleSignalToStats.size());
}

void MedLabels::create_incidence_file(const string &file_path, const string &rep_path, int age_bin, int min_age,
	int max_age, bool use_kaplan_meir, const string &sampler_name, const string &sampler_args, const string &debug_file) const {
	MedSamplingStrategy *sampler = MedSamplingStrategy::make_sampler(sampler_name, sampler_args);

	MedRepository rep;
	vector<int> pids;
	get_pids(pids);
	vector<string> signal_to_read = { "BYEAR", "GENDER" };
	if (rep.read_all(rep_path, pids, signal_to_read) < 0)
		MTHROW_AND_ERR("FAILED reading repository %s\n", rep_path.c_str());
	min_age = int(min_age / age_bin) * age_bin;
	MedSamples incidence_samples;
	sampler->init_sampler(rep);
	MLOG("Sampling for incidence stats...\n");
	create_samples(sampler, incidence_samples);
	MLOG("Done...\n");
	delete sampler;
	ofstream fw_debug;
	if (!debug_file.empty())
		fw_debug.open(debug_file);
	int time_period = labeling_params.time_to - labeling_params.time_from;

	vector<int> all_cnts = { 0,0 };
	int bin_counts = (max_age - min_age) / age_bin + 1;
	vector<pair<int, int>> counts(bin_counts), male_counts(bin_counts), female_counts(bin_counts);
	vector<vector<vector<pair<int, int>>>> filters_idxs(2); //gender, age_bin, vector of idx
	if (use_kaplan_meir) {
		for (size_t i = 0; i < filters_idxs.size(); ++i)
			filters_idxs[i].resize(bin_counts);
	}
	for (int i = min_age; i < max_age; i += age_bin)
		counts[(i - min_age) / age_bin] = pair<int, int>(0, 0);
	int byear_sid = rep.sigs.sid("BYEAR");
	int gender_sid = rep.sigs.sid("GENDER");
	int len;
	for (int i = 0; i < incidence_samples.idSamples.size(); ++i)
		for (int j = 0; j < incidence_samples.idSamples[i].samples.size(); ++j) {
			int pid = incidence_samples.idSamples[i].samples[j].id;
			int byear = (int)((((SVal *)rep.get(pid, byear_sid, len))[0]).val);
			int age = int(incidence_samples.idSamples[i].samples[j].time / 10000) - byear;
			int gender = (int)((((SVal *)rep.get(pid, gender_sid, len))[0]).val);
			//int bin = age_bin*(age / age_bin);
			int age_index = (age - min_age) / age_bin;
			if (age < min_age || age > max_age || age_index < 0 || age_index >= counts.size())
				continue;

			++counts[age_index].first;
			if (gender == GENDER_MALE)
				++male_counts[age_index].first;
			else if (gender == GENDER_FEMALE)
				++female_counts[age_index].first;
			all_cnts[0]++;
			if (incidence_samples.idSamples[i].samples[j].outcome > 0) {
				++counts[age_index].second;
				all_cnts[1]++;
				if (gender == GENDER_MALE)
					++male_counts[age_index].second;
				else if (gender == GENDER_FEMALE)
					++female_counts[age_index].second;
				/*if (age_index*age_bin + min_age == 95)
				MLOG("DEBUG:: pid=%d, time=%d, outcomeTime=%d\n", pid,
				incidence_samples.idSamples[i].samples[j].time,
				incidence_samples.idSamples[i].samples[j].outcomeTime);*/
			}
			if (!debug_file.empty()) {
				//Debug: pid, year, outcome, age, gender
				fw_debug << incidence_samples.idSamples[i].samples[j].id << "\t"
					<< incidence_samples.idSamples[i].samples[j].time << "\t"
					<< incidence_samples.idSamples[i].samples[j].outcome << "\t"
					<< age << "\t" << gender << "\n";
			}

			if (use_kaplan_meir)
				filters_idxs[gender - 1][age_index].push_back(pair<int, int>(i, j));
		}
	if (!debug_file.empty())
		fw_debug.close();

	if (use_kaplan_meir) {
		int kaplan_meier_controls_count = 100000;
		//for each group - Age, Age+Gender... whatever
		ofstream of_new;
		if (file_path != "/dev/null") {
			of_new.open(file_path + ".new_format");
			if (!of_new.good())
				MTHROW_AND_ERR("IO Error: can't write \"%s\"\n", (file_path + ".new_format").c_str());
			of_new << "AGE_BIN" << "\t" << age_bin << "\n";
			of_new << "AGE_MIN" << "\t" << min_age << "\n";
			of_new << "AGE_MAX" << "\t" << max_age << "\n";
			of_new << "OUTCOME_VALUE" << "\t" << "0.0" << "\n";
			of_new << "OUTCOME_VALUE" << "\t" << "1.0" << "\n";

			for (size_t gender = 0; gender < 2; ++gender) {
				string gender_str = gender + 1 == GENDER_MALE ? "MALE" : "FEMALE";
				for (int c = 0; c < bin_counts; ++c) {
					vector<pair<int, int>> &filters = filters_idxs[gender][c];
					double prob = medial::stats::kaplan_meir_on_samples(incidence_samples, time_period, &filters);
					if (prob > 0 && prob < 1) {
						int age = c * age_bin + min_age;
						//print to file:
						MLOG("%s:Ages[%d - %d]:%d :: %2.2f%% (size=%zu) (kaplan meier)\n", gender_str.c_str(), age, age + age_bin,
							age + age_bin / 2, 100 * prob, filters.size());

						kaplan_meier_controls_count = (int)filters.size();
						if (age >= min_age && age <= max_age) {
							of_new << "STATS_ROW" << "\t" << gender_str << "\t" <<
								age + age_bin / 2 << "\t" << "0.0" << "\t" << int(kaplan_meier_controls_count * (1 - prob)) << "\n";
							of_new << "STATS_ROW" << "\t" << gender_str << "\t" <<
								age + age_bin / 2 << "\t" << "1.0" << "\t" << int(kaplan_meier_controls_count * prob) << "\n";
						}
					}
				}
			}
			of_new.close();
		}
	}
	else {
		//regular inc calc
		MLOG("Total counts: 0: %d 1: %d : inc %f\n", all_cnts[0], all_cnts[1],
			(float)all_cnts[1] / all_cnts[0]);
		int nlines = 0;
		for (int c = 0; c < counts.size(); ++c) {
			int age = c * age_bin + min_age;
			int n0 = counts[c].first;
			int n1 = counts[c].second;

			if (age >= min_age && age < max_age) nlines++;

			if (n0 > 0)
				MLOG("Ages: %d - %d : %d : 0: %d 1: %d : %f\n", age, age + age_bin,
					age + age_bin / 2, n0, n1, (n0 > 0) ? (float)n1 / n0 : 0);
		}

		ofstream of(file_path);
		if (!of.good())
			MTHROW_AND_ERR("IO Error: can't write \"%s\"\n", file_path.c_str());

		of << "KeySize 1\n";
		of << "Nkeys " << nlines << "\n";
		of << "1.0\n";
		for (int c = 0; c < counts.size(); ++c) {
			int age = c * age_bin + min_age;
			int n0 = counts[c].first;
			int n1 = counts[c].second;

			if (age >= min_age && age < max_age)
				of << age + age_bin / 2 << " " << n1 << " " << n0 - n1 << "\n";
		}
		of.close();

		//New Format:
		ofstream of_new;
		if (file_path != "/dev/null") {
			of_new.open(file_path + ".new_format");
			if (!of_new.good())
				MTHROW_AND_ERR("IO Error: can't write \"%s\"\n", (file_path + ".new_format").c_str());

			of_new << "AGE_BIN" << "\t" << age_bin << "\n";
			of_new << "AGE_MIN" << "\t" << min_age << "\n";
			of_new << "AGE_MAX" << "\t" << max_age << "\n";
			of_new << "OUTCOME_VALUE" << "\t" << "0.0" << "\n";
			of_new << "OUTCOME_VALUE" << "\t" << "1.0" << "\n";

			for (int c = 0; c < counts.size(); ++c) {

				int age = c * age_bin + min_age;
				int male_n0 = male_counts[c].first;
				int male_n1 = male_counts[c].second;

				if (age >= min_age && age <= max_age && male_n0 > 0) {
					of_new << "STATS_ROW" << "\t" << "MALE" << "\t" <<
						age + age_bin / 2 << "\t" << "0.0" << "\t" << male_n0 - male_n1 << "\n";
					of_new << "STATS_ROW" << "\t" << "MALE" << "\t" <<
						age + age_bin / 2 << "\t" << "1.0" << "\t" << male_n1 << "\n";
				}

				int female_n0 = female_counts[c].first;
				int female_n1 = female_counts[c].second;

				if (age >= min_age && age <= max_age && female_n0 > 0) {
					of_new << "STATS_ROW" << "\t" << "FEMALE" << "\t" <<
						age + age_bin / 2 << "\t" << "0.0" << "\t" << female_n0 - female_n1 << "\n";
					of_new << "STATS_ROW" << "\t" << "FEMALE" << "\t" <<
						age + age_bin / 2 << "\t" << "1.0" << "\t" << female_n1 << "\n";
				}
			}
			of_new.close();
		}
	}
}

void MedLabels::create_samples(const MedSamplingStrategy *sampler, MedSamples &samples, bool show_conflicts) const {
	if (pid_reg_records.empty())
		MTHROW_AND_ERR("Error in MedLabels::get_samples - please init MedLabels by calling prepare_from_registry\n");
	unordered_map<int, vector<pair<int, int>>> pid_time_ranges;
	//create availible times for each pid:
	if (has_censor_reg()) {
		for (auto it = pid_censor_records.begin(); it != pid_censor_records.end(); ++it)
		{
			if (pid_reg_records.find(it->first) == pid_reg_records.end())
				continue; //not in registry skip
			for (size_t i = 0; i < it->second.size(); ++i)
			{
				pair<int, int> tm(it->second[i]->start_date, it->second[i]->end_date);
				pid_time_ranges[it->first].push_back(tm);
			}

		}
	}
	else {
		MWARN("Warning MedLabels::create_samples - no censor registry\n");
		//will take minimal start_time and maximal end_time as availble time range
		for (auto it = pid_reg_records.begin(); it != pid_reg_records.end(); ++it)
		{
			if (!it->second.empty()) {
				pair<int, int> tm(it->second.front()->start_date, it->second.front()->end_date);
				for (size_t i = 1; i < it->second.size(); ++i)
				{
					if (it->second[i]->start_date < tm.first)
						tm.first = it->second[i]->start_date;
					if (it->second[i]->end_date > tm.second)
						tm.second = it->second[i]->end_date;

				}
				pid_time_ranges[it->first].push_back(tm);
			}
		}
	}
	unordered_map<int, vector<int>> pid_times;
	sampler->get_sampling_options(pid_time_ranges, pid_times);
	int conflict_count = 0, done_count = 0, no_censor = 0, no_rule = 0;
	int max_to_shown = 5;

	for (auto it = pid_times.begin(); it != pid_times.end(); ++it)
	{
		if (has_censor_reg() && !has_censor_reg(it->first)) {
			++no_censor;
			continue; //filter sample
		}
		vector<int> &times = it->second;
		MedIdSamples smp_id(it->first);
		SamplingRes r = get_samples(it->first, times, smp_id.samples, show_conflicts);
		done_count += r.done_cnt;  no_rule += r.no_rule_cnt; conflict_count += r.conflict_cnt;
		if (r.conflict_cnt > 0 && max_to_shown > 0) {
			--max_to_shown;
			if (max_to_shown <= 0)
				show_conflicts = false;
		}

		if (!smp_id.samples.empty())
			samples.idSamples.push_back(smp_id);
	}
	if (no_rule > 0)
		MLOG("WARNING MedLabels::create_samples - has %d samples with no rules for time window\n", no_rule);
	if (no_censor > 0) {
		if (has_censor_reg())
			MLOG("WARNING MedLabels::create_samples - has %d patients with no censor dates\n", no_censor);
		else
			MLOG("WARNING MedLabels::create_samples - no censoring time region was given\n");
	}
	if (conflict_count > 0)
		MLOG("Sampled registry with %d conflicts. has %d registry records\n", conflict_count, done_count);
	//do sort:
	samples.sort_by_id_date();
}

void MedLabels::relabel_samples(MedSamples &samples, bool show_conflicts) const {
	if (pid_reg_records.empty())
		MTHROW_AND_ERR("Error in MedLabels::relabel_samples - please init MedLabels by calling prepare_from_registry\n");
	int conflict_count = 0, done_count = 0, no_censor = 0, no_rule = 0;
	int max_to_shown = 5;

	vector<int> keep_ids;
	keep_ids.reserve(samples.idSamples.size());
	for (size_t i = 0; i < samples.idSamples.size(); ++i) {
		int pid = samples.idSamples[i].id;
		vector<int> times(samples.idSamples[i].samples.size());
		samples.idSamples[i].get_times(times);
		samples.idSamples[i].samples.clear();
		SamplingRes r = get_samples(pid, times, samples.idSamples[i].samples, show_conflicts);
		done_count += r.done_cnt;  no_rule += r.no_rule_cnt; conflict_count += r.conflict_cnt;
		if (r.conflict_cnt > 0 && max_to_shown > 0) {
			--max_to_shown;
			if (max_to_shown <= 0)
				show_conflicts = false;
		}
		if (!samples.idSamples[i].samples.empty())
			keep_ids.push_back((int)i);
	}
	//filter only keep:
	MedSamples filtered;
	filtered.time_unit = samples.time_unit;
	filtered.raw_format = samples.raw_format;
	filtered.idSamples.resize(keep_ids.size());
	for (size_t i = 0; i < keep_ids.size(); ++i)
	{
		filtered.idSamples[i].id = samples.idSamples[keep_ids[i]].id;
		filtered.idSamples[i].split = samples.idSamples[keep_ids[i]].split;
		filtered.idSamples[i].samples = move(samples.idSamples[keep_ids[i]].samples);
	}
	samples = move(filtered);

	if (no_rule > 0)
		MLOG("WARNING MedLabels::relabel_samples - has %d samples with no rules for time window\n", no_rule);
	if (no_censor > 0) {
		if (has_censor_reg())
			MLOG("WARNING MedLabels::relabel_samples - has %d patients with no censor dates\n", no_censor);
		else
			MLOG("WARNING MedLabels::relabel_samples - no censoring time region was given\n");
	}
	if (conflict_count > 0)
		MLOG("Sampled registry with %d conflicts. has %d registry records\n", conflict_count, done_count);
	//do sort:
	samples.sort_by_id_date();
}