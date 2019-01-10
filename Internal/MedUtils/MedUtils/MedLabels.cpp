#include "MedLabels.h"
#include <algorithm>
#include <Logger/Logger/Logger.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedUtils/MedUtils/MedUtils.h>

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
	int &no_rule_found, int &conflict_count, int &done_count, bool filter_no_censor) {
	int curr_index = 0, final_selected = -1;
	float reg_val = -1;
	int reg_time = -1;
	if (pid_records.empty())
		return;
	MedSample smp;
	smp.time = pred_time;
	smp.id = pid_records.front()->pid;

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
		}
		else if (reg_val != pid_records[curr_index]->registry_value) {
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
			else {
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
#pragma omp atomic
			++conflict_count;
			break;
		}

		++curr_index;
	}

	if (reg_time != -1) {
		smp.outcomeTime = reg_val > 0 ? pid_records[final_selected]->start_date : reg_time;
		smp.outcome = reg_val;
		idSamples.push_back(smp);
		++done_count;
	}
}

void medial::signal_hierarchy::getRecords_Hir(int pid, vector<UniversalSigVec> &signals, MedDictionarySections &dict,
	const string &signalHirerchyType,
	vector<MedRegistryRecord> &res) {
	UniversalSigVec &signalVal = signals[0];
	int max_search_depth = 3;

	for (int i = 0; i < signalVal.len; ++i)
	{
		MedRegistryRecord rec;
		rec.pid = pid;
		rec.start_date = signalVal.Time(i);
		rec.end_date = medial::repository::DateAdd(rec.start_date, 1);
		rec.registry_value = signalVal.Val(i);
		res.push_back(rec);
		if (signalVal.Val(i) <= 0)
			continue; //has no hirerachy
					  //take care of hirerachy:

		if (signalHirerchyType.empty() || signalHirerchyType == "None")
			continue;
		string s = get_readcode_code(dict, (int)signalVal.Val(i), signalHirerchyType);
		if (s.empty())
			continue;

		vector<int> nums = parents_code_hierarchy(dict, s, signalHirerchyType, max_search_depth);
		for (size_t k = 0; k < nums.size(); ++k)
		{
			if (nums[k] <= 0)
				continue;
			MedRegistryRecord rec2;

			rec2.pid = pid;
			rec2.start_date = signalVal.Time(i);
			rec2.end_date = medial::repository::DateAdd(rec2.start_date, 1);
			rec2.registry_value = (float)nums[k];
			res.push_back(rec2);
		}
	}
}

MedLabels::MedLabels(const LabelParams &params) {
	labeling_params = params;
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

SamplingRes MedLabels::get_samples(int pid, int time, vector<MedSample> &samples) const {
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
			labeling_params.label_interaction_mode, labeling_params.censor_interaction_mode, labeling_params.conflict_method,
			samples, r.no_rule_cnt, r.conflict_cnt, r.done_cnt, false);
	}
	else
		++r.miss_pid_in_reg_cnt;

	return r;
}

SamplingRes MedLabels::get_samples(int pid, const vector<int> &times, vector<MedSample> &samples) const {
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
				labeling_params.label_interaction_mode, labeling_params.censor_interaction_mode, labeling_params.conflict_method,
				samples, r.no_rule_cnt, r.conflict_cnt, r.done_cnt, false);
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

void update_loop(int pos, int ageBin_index, float ageBin, const MedRegistryRecord &sigRec,
	map<float, map<float, vector<int>>> &signalToStats, vector<unordered_map<float, int>> &val_seen_pid_pos) {
	if ((pos == 3 && val_seen_pid_pos[ageBin_index][sigRec.registry_value] / 2 > 0) ||
		(pos == 2 && val_seen_pid_pos[ageBin_index][sigRec.registry_value] % 2 > 0)) {
		return; //continue;
	}
	val_seen_pid_pos[ageBin_index][sigRec.registry_value] += pos - 1;
	//update cnts:
#pragma omp critical 
	{
		vector<int> *cnts = &(signalToStats[sigRec.registry_value][ageBin]);
		if (cnts->empty())
			cnts->resize(4); // first time
		++(*cnts)[pos];
	}
}

void MedLabels::calc_signal_stats(const string &repository_path, const string &signal_name,
	const string &signalHirerchyType, int ageBinValue, MedSamplingStrategy &sampler, const LabelParams &inc_labeling_params,
	map<float, map<float, vector<int>>> &maleSignalToStats,
	map<float, map<float, vector<int>>> &femaleSignalToStats,
	const string &debug_file, const unordered_set<float> &debug_vals) const {
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
		vector<UniversalSigVec> patientFile(1);
		dataManager.uget(pid, signalCode, patientFile[0]);

		vector<MedRegistryRecord> signal_vals;
		medial::signal_hierarchy::getRecords_Hir(pid, patientFile, dataManager.dict, signalHirerchyType, signal_vals);

		vector<unordered_map<float, int>> val_seen_pid_pos(age_bin_count); //for age bin index and value (it's for same pid so gender doesnt change) - if i saw the value already
		for (MedRegistryRecord sigRec : signal_vals)
		{
			if (sigRec.registry_value <= 0) {
				continue;
			}
			int pos;
			vector<int> cnts;
			float ageBin;
			int ageBin_index;

			ageBin = float(ageBinValue * floor(double(medial::repository::DateDiff(BDate, sigRec.start_date)) / ageBinValue));
			ageBin_index = int((ageBin - min_age) / ageBinValue);
			if (ageBin < min_age || ageBin > max_age)
				continue; //skip out of range...

			vector<MedSample> found_samples;
			get_samples(pid, sigRec.start_date, found_samples);
			for (const MedSample &smp : found_samples)
			{
				pos = 2;
				//pos += 1; //registry_value > 0 - otherwise skip this
				pos += int(smp.outcome > 0);
				if (gender == GENDER_MALE)
					update_loop(pos, ageBin_index, ageBin, sigRec, maleSignalToStats, val_seen_pid_pos);
				else
					update_loop(pos, ageBin_index, ageBin, sigRec, femaleSignalToStats, val_seen_pid_pos);
				if (!debug_file.empty() && debug_vals.find(sigRec.registry_value) != debug_vals.end()) {
#pragma omp critical
					dbg_file << pid << "\t" << sigRec.start_date << "\t" << sigRec.registry_value
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

	//update values prevalence
	int warn_cnt = 0; int max_warns = 5;
	for (auto it = vals.begin(); it != vals.end(); ++it)
	{
		if (maleSignalToStats.find(*it) != maleSignalToStats.end())
			for (auto jt = maleSignalToStats[*it].begin(); jt != maleSignalToStats[*it].end(); ) {
				if (male_total_prevalence.find(jt->first) == male_total_prevalence.end()) {
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedRegistry::calc_signal_stats - Sample is too small, no incidences for age_bin=%d in males (value was=%f, cnts=[%d, %d])\n"
							, int(jt->first), *it, maleSignalToStats[*it][jt->first][2], maleSignalToStats[*it][jt->first][3]);
					}
					jt = maleSignalToStats[*it].erase(jt);
					continue;
				}
				maleSignalToStats[*it][jt->first][0] = male_total_prevalence[jt->first][0] - maleSignalToStats[*it][jt->first][2];
				maleSignalToStats[*it][jt->first][1] = male_total_prevalence[jt->first][1] - maleSignalToStats[*it][jt->first][3];
				if (maleSignalToStats[*it][jt->first][0] < 0) {
					maleSignalToStats[*it][jt->first][0] = 0;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedRegistry::calc_signal_stats - Control Male age_bin=%d, signal_value=%f, total=%d, signal=%d\n",
							int(jt->first), *it, male_total_prevalence[jt->first][0], maleSignalToStats[*it][jt->first][2]);
					}
				}
				if (maleSignalToStats[*it][jt->first][1] < 0) {
					maleSignalToStats[*it][jt->first][1] = 0;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedRegistry::calc_signal_stats - Cases Male age_bin=%d, signal_value=%f, total=%d, signal=%d\n",
							int(jt->first), *it, male_total_prevalence[jt->first][1], maleSignalToStats[*it][jt->first][3]);
					}
				}
				++jt;
			}
		if (femaleSignalToStats.find(*it) != femaleSignalToStats.end())
			for (auto jt = femaleSignalToStats[*it].begin(); jt != femaleSignalToStats[*it].end();) {
				if (female_total_prevalence.find(jt->first) == female_total_prevalence.end()) {

					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedRegistry::calc_signal_stats - Sample is too small, no incidences for age_bin=%d in females (value was=%f, cnts=[%d, %d])\n"
							, int(jt->first), *it, femaleSignalToStats[*it][jt->first][2], femaleSignalToStats[*it][jt->first][3]);
					}
					jt = femaleSignalToStats[*it].erase(jt);
					continue;
				}
				femaleSignalToStats[*it][jt->first][0] = female_total_prevalence[jt->first][0] - femaleSignalToStats[*it][jt->first][2];
				femaleSignalToStats[*it][jt->first][1] = female_total_prevalence[jt->first][1] - femaleSignalToStats[*it][jt->first][3];
				if (femaleSignalToStats[*it][jt->first][0] < 0) {
					femaleSignalToStats[*it][jt->first][0] = 0;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedRegistry::calc_signal_stats - Control Female age_bin=%d, signal_value=%f, total=%d, signal=%d\n",
							int(jt->first), *it, female_total_prevalence[jt->first][0], femaleSignalToStats[*it][jt->first][2]);
					}
				}
				if (femaleSignalToStats[*it][jt->first][1] < 0) {
					femaleSignalToStats[*it][jt->first][1] = 0;
					if (warn_cnt < max_warns) {
						++warn_cnt;
						MWARN("Warning: MedRegistry::calc_signal_stats - Cases Female age_bin=%d, signal_value=%f, total=%d, signal=%d\n",
							int(jt->first), *it, female_total_prevalence[jt->first][1], femaleSignalToStats[*it][jt->first][3]);
					}
				}
				++jt;
			}
	}

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
	vector<vector<int>> sorted_times(bin_counts);
	vector<vector<vector<pair<int, int>>>> times_indexes(bin_counts);
	vector<vector<bool>> all_times;
	if (use_kaplan_meir) {
		all_times.resize(bin_counts);
		times_indexes.resize(bin_counts);
		sorted_times.resize(bin_counts);
		for (size_t i = 0; i < bin_counts; ++i)
		{
			sorted_times[i].reserve(time_period + 1);
			all_times[i].resize(time_period + 1, false);
		}
	}
	for (int i = min_age; i < max_age; i += age_bin)
		counts[(i - min_age) / age_bin] = pair<int, int>(0, 0);
	int byear_sid = rep.sigs.sid("BYEAR");
	int gender_sid = rep.sigs.sid("GENDER");
	int len;
	for (size_t i = 0; i < incidence_samples.idSamples.size(); ++i)
		for (size_t j = 0; j < incidence_samples.idSamples[i].samples.size(); ++j) {
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

			if (use_kaplan_meir) {
				int time_diff = int(365 * medial::repository::DateDiff(incidence_samples.idSamples[i].samples[j].time,
					incidence_samples.idSamples[i].samples[j].outcomeTime));
				if (time_diff > time_period)
					time_diff = time_period;
				if (!all_times[age_index][time_diff]) {
					sorted_times[age_index].push_back(time_diff);
					all_times[age_index][time_diff] = true;
				}
			}
		}
	if (use_kaplan_meir) {
		for (int c = 0; c < sorted_times.size(); ++c)
		{
			sort(sorted_times[c].begin(), sorted_times[c].end());
			times_indexes[c].resize(sorted_times[c].size());
		}
		//prepare times_indexes:
		for (size_t i = 0; i < incidence_samples.idSamples.size(); ++i)
			for (size_t j = 0; j < incidence_samples.idSamples[i].samples.size(); ++j) {
				int pid = incidence_samples.idSamples[i].samples[j].id;
				int byear = (int)((((SVal *)rep.get(pid, byear_sid, len))[0]).val);
				int age = int(incidence_samples.idSamples[i].samples[j].time / 10000) - byear;
				int age_index = (age - min_age) / age_bin;
				if (age < min_age || age > max_age || age_index < 0 || age_index >= counts.size())
					continue;

				int time_diff = int(365 * medial::repository::DateDiff(incidence_samples.idSamples[i].samples[j].time,
					incidence_samples.idSamples[i].samples[j].outcomeTime));
				int original_time = time_diff;
				if (time_diff > time_period)
					time_diff = time_period;
				int ind = medial::process::binary_search_index(sorted_times[age_index].data(),
					sorted_times[age_index].data() + sorted_times[age_index].size() - 1, time_diff);
				if (incidence_samples.idSamples[i].samples[j].outcome <= 0 ||
					original_time <= time_period)
					times_indexes[age_index][ind].push_back(pair<int, int>((int)i, (int)j));
			}
	}

	if (!debug_file.empty())
		fw_debug.close();

	if (use_kaplan_meir) {
		bool warn_shown = false;
		int kaplan_meier_controls_count = 100000;
		//for each group - Age, Age+Gender... whatever
		ofstream of_new(file_path + ".new_format");
		if (!of_new.good())
			MTHROW_AND_ERR("IO Error: can't write \"%s\"\n", (file_path + ".new_format").c_str());
		of_new << "AGE_BIN" << "\t" << age_bin << "\n";
		of_new << "AGE_MIN" << "\t" << min_age << "\n";
		of_new << "AGE_MAX" << "\t" << max_age << "\n";
		of_new << "OUTCOME_VALUE" << "\t" << "0.0" << "\n";
		of_new << "OUTCOME_VALUE" << "\t" << "1.0" << "\n";

		for (int c = 0; c < sorted_times.size(); ++c)
		{
			double total_controls_all = 0;
			for (size_t sort_ind = 0; sort_ind < sorted_times[c].size(); ++sort_ind) {
				const vector<pair<int, int>> &index_order = times_indexes[c][sort_ind];
				//update only controls count for group - do for all groups
				//to get total count from all time windows kaplan meier
				for (const pair<int, int> &p_i_j : index_order)
					if (incidence_samples.idSamples[p_i_j.first].samples[p_i_j.second].outcome <= 0)
						++total_controls_all;
			}

			double controls = 0, cases = 0, prob = 1;
			for (size_t sort_ind = 0; sort_ind < sorted_times[c].size(); ++sort_ind) {
				const vector<pair<int, int>> &index_order = times_indexes[c][sort_ind];
				for (const pair<int, int> &p_i_j : index_order) {
					//keep update kaplan meir in time point
					if (incidence_samples.idSamples[p_i_j.first].samples[p_i_j.second].outcome > 0)
						++cases;
					else
						++controls;
				}
				//reset kaplan meir - flash last time prob
				if (!warn_shown && total_controls_all < 10) {
					MWARN("the kaplan_meir left with small amount of controls - "
						" try increasing the sampling / use smaller time window because the"
						" registry has not so long period of tracking patients\n");
					warn_shown = true;
				}
				if (total_controls_all > 0 || cases > 0)
					prob *= total_controls_all / (cases + total_controls_all);
				total_controls_all -= controls; //remove controls from current time-window - they are now censored
				controls = 0; cases = 0;
			}
			prob = 1 - prob;
			if (prob > 0 && prob < 1) {
				int age = c * age_bin + min_age;
				//print to file:
				MLOG("Ages[%d - %d]:%d :: %2.2f%% (kaplan meier)\n", age, age + age_bin,
					age + age_bin / 2, 100 * prob);

				if (age >= min_age && age <= max_age) {
					of_new << "STATS_ROW" << "\t" << "MALE" << "\t" <<
						age + age_bin / 2 << "\t" << "0.0" << "\t" << int(kaplan_meier_controls_count * (1 - prob)) << "\n";
					of_new << "STATS_ROW" << "\t" << "MALE" << "\t" <<
						age + age_bin / 2 << "\t" << "1.0" << "\t" << int(kaplan_meier_controls_count * prob) << "\n";

					of_new << "STATS_ROW" << "\t" << "FEMALE" << "\t" <<
						age + age_bin / 2 << "\t" << "0.0" << "\t" << int(kaplan_meier_controls_count * (1 - prob)) << "\n";
					of_new << "STATS_ROW" << "\t" << "FEMALE" << "\t" <<
						age + age_bin / 2 << "\t" << "1.0" << "\t" << int(kaplan_meier_controls_count * prob) << "\n";
				}
			}
		}
		of_new.close();
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
		ofstream of_new(file_path + ".new_format");
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

void MedLabels::create_samples(const MedSamplingStrategy *sampler, MedSamples &samples) const {
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

	for (auto it = pid_times.begin(); it != pid_times.end(); ++it)
	{
		if (has_censor_reg() && !has_censor_reg(it->first)) {
			++no_censor;
			continue; //filter sample
		}
		vector<int> &times = it->second;
		MedIdSamples smp_id(it->first);
		SamplingRes r = get_samples(it->first, times, smp_id.samples);
		done_count += r.done_cnt;  no_rule += r.no_rule_cnt; conflict_count = r.conflict_cnt;
		if (!smp_id.samples.empty())
			samples.idSamples.push_back(smp_id);
	}
	if (no_rule > 0)
		MLOG("WARNING MedSamplingYearly:do_sample - has %d samples with no rules for time window\n", no_rule);
	if (no_censor > 0) {
		if (has_censor_reg())
			MLOG("WARNING MedSamplingYearly:do_sample - has %d patients with no censor dates\n", no_censor);
		else
			MLOG("WARNING MedSamplingYearly:do_sample - no censoring time region was given\n");
	}
	if (conflict_count > 0)
		MLOG("Sampled registry with %d conflicts. has %d registry records\n", conflict_count, done_count);
	//do sort:
	samples.sort_by_id_date();
}