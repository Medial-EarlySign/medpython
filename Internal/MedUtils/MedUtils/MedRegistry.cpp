#include "MedRegistry.h"
#include <fstream>
#include <string>
#include <iostream>
#include <random>
#include <map>
#include "Logger/Logger/Logger.h"
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void MedRegistry::read_text_file(const string &file_path) {
	ifstream fp(file_path);
	if (!fp.good())
		MTHROW_AND_ERR("IOError: can't read file %s\n", file_path.c_str());
	registry_records.clear();
	string line;
	int lineNum = 0;
	vector<string> tokens;
	tokens.reserve(7);
	//Format: [ PID, Start_Date, End_Date, min_allowed_date, max_allowed_date, Age, RegistryValue ]
	while (getline(fp, line))
	{
		++lineNum;
		if (line.size() < 3 || line.at(0) == '#') {
			continue; //empty or comment line
		}
		tokens.clear();
		boost::split(tokens, line, boost::is_any_of("\t"));
		if (tokens.size() != 7) {
			cerr << "Bad File Format in line " << lineNum << " got \"" << line << "\"" << " parsed " << tokens.size() << " tokens" << endl;
			throw out_of_range("File has bad format");
		}


		MedRegistryRecord pr;
		pr.pid = stoi(tokens[0]);
		pr.start_date = stoi(tokens[1]);
		pr.end_date = stoi(tokens[2]);
		pr.min_allowed_date = stoi(tokens[3]);
		pr.max_allowed_date = stoi(tokens[4]);
		pr.age = stoi(tokens[5]);
		pr.registry_value = stof(tokens[6]);
		registry_records.push_back(pr);

	}
	fp.close();
	MLOG("[Read %d registry records from %s]\n", (int)registry_records.size(), file_path.c_str());
}

void MedRegistry::write_text_file(const string &file_path) {
	char delim = '\t';
	ofstream fw(file_path);
	if (!fw.good())
		MTHROW_AND_ERR("IOError: can't write file %s\n", file_path.c_str());
	fw << "# Format: PID, Start_Date, End_Date, min_allowed_date, max_allowed_date, Age, RegistryValue" << endl;
	fw << "# Created By Script - Insert Comments following #..." << endl;
	fw << endl;
	for (size_t i = 0; i < registry_records.size(); ++i)
		fw << registry_records[i].pid << delim << registry_records[i].start_date << delim << registry_records[i].end_date << delim
		<< registry_records[i].min_allowed_date << delim << registry_records[i].max_allowed_date << delim << registry_records[i].age
		<< delim << registry_records[i].registry_value << "\n";

	fw.flush();
	fw.close();
	MLOG("[Wrote %d registry records to %s]\n", (int)registry_records.size(), file_path.c_str());
}

void MedRegistry::create_registry(MedPidRepository &dataManager) {

	MLOG_D("Creating registry...\n");

	time_t start = time(NULL);
	time_t last_time_print = start;
	double duration;
	int prog_pid = 0;
	int bDateCode = dataManager.sigs.sid("BDATE");
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < dataManager.pids.size(); ++i)
	{
		vector<UniversalSigVec> sig_vec((int)signalCodes.size());
		for (size_t k = 0; k < sig_vec.size(); ++k)
			dataManager.uget(dataManager.pids[i], signalCodes[k], sig_vec[k]);
		int birth = medial::repository::get_value(dataManager, dataManager.pids[i], bDateCode);
		vector<MedRegistryRecord> vals;
		get_registry_records(dataManager.pids[i], birth, sig_vec, vals);

		if (registry_records.size() >= 1000000000 && registry_records.size() + vals.size() > registry_records.capacity())
			MWARN("BIG DICTIONARY SIZE BEFORE %d may crash\n", (int)registry_records.size());
#pragma omp critical 
		{
			registry_records.insert(registry_records.end(), vals.begin(), vals.end());
			++prog_pid;
		}

		if (prog_pid % 10000 == 0 && (int)difftime(time(NULL), last_time_print) >= 15) {
			last_time_print = time(NULL);
			float time_elapsed = (float)difftime(time(NULL), start);
			float estimate_time = float(dataManager.pids.size() - prog_pid) / prog_pid * time_elapsed;
			MLOG("Processed %d out of %d(%2.2f%) time elapsed: %2.1f Minutes, estimate time to finish %2.1f Minutes. has %d records\n",
				prog_pid, dataManager.pids.size(), 100.0*(prog_pid / float(dataManager.pids.size())), time_elapsed / 60, estimate_time / 60.0, registry_records.size());
		}
	}

	duration = difftime(time(NULL), start);
	MLOG("Finished creating registy in %d seconds\n", (int)duration);
	dataManager.clear();
}

void MedRegistry::get_registry_creation_codes(vector<int> &signal_codes)
{
	signal_codes = signalCodes;
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
		rec.age = -1;
		rec.start_date = signalVal.Date(i);
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
			rec2.age = -1;
			rec2.start_date = signalVal.Date(i);
			rec2.end_date = medial::repository::DateAdd(rec2.start_date, 1);
			rec2.registry_value = (float)nums[k];
			res.push_back(rec2);
		}
	}
}

/// test that sig_start is in allowed range (means signal date is in allowed time).
/// If looking backward - checks that end_date is in allowed range (backward search).
/// use_whole changes the check to "contain" - which is relevents for controls to not passed reg_end_date.
/// In cases we don't need it
bool date_intersection(int min_allowed_date, int max_allowed_date, int reg_start, int reg_end, int sig_date,
	int time_window_from, int time_window_to, bool use_whole) {
	//Registry, Signal
	int sig_start_date = medial::repository::DateAdd(sig_date, time_window_from);
	int sig_end_date = medial::repository::DateAdd(sig_date, time_window_to);
	int reffer_date = sig_start_date;
	if (time_window_from < 0) //if looking backward force end_date to be in allowed
		reffer_date = sig_end_date;
	if (reffer_date > max_allowed_date || reffer_date < min_allowed_date)
		return false;

	if (time_window_from < 0)
		return (sig_start_date <= reg_end) && (!use_whole || sig_start_date >= reg_start);
	else
		return (sig_end_date >= reg_start) && (!use_whole || sig_end_date <= reg_end);
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

void MedRegistry::calc_signal_stats(const string &repository_path, int signalCode,
	const string &signalHirerchyType, int ageBinValue, int time_window_from, int time_window_to,
	MedSamplingStrategy &sampler,
	map<float, map<float, vector<int>>> &maleSignalToStats,
	map<float, map<float, vector<int>>> &femaleSignalToStats,
	const string &debug_file, const unordered_set<float> &debug_vals) {
	MedRepository dataManager;
	time_t start = time(NULL);
	int duration;

	if (time_window_from > time_window_to) {
		MWARN("Warning: you gave time window params in wrong order [%d, %d]."
			" switching them for you\n", time_window_from, time_window_to);
		int tmp = time_window_to;
		time_window_to = time_window_from;
		time_window_from = tmp;
	}

	vector<int> pids;
	MedSamples incidence_samples;

	MLOG("Sampling for incidence stats...\n");
	sampler.do_sample(registry_records, incidence_samples);
	incidence_samples.get_ids(pids);
	duration = (int)difftime(time(NULL), start);
	MLOG("Done in %d seconds!\n", duration);

	vector<int> readSignals;
	readSignals.push_back(signalCode);
	int curr_level = global_logger.levels[LOG_REP];
	global_logger.levels[LOG_REP] = MAX_LOG_LEVEL;
	dataManager.init(repository_path);
	global_logger.levels[LOG_REP] = curr_level;
	int genderCode = dataManager.sigs.sid("GENDER");
	int bdateCode = dataManager.sigs.sid("BDATE");
	readSignals.push_back(genderCode);
	readSignals.push_back(bdateCode);
	MLOG("Fetching signal %s using repository %s\n", dataManager.sigs.name(signalCode).c_str(),
		repository_path.c_str());
	if (dataManager.read_all(repository_path, pids, readSignals) < 0)
		MTHROW_AND_ERR("error reading from repository %s\n", repository_path.c_str());
	vector<int> &all_pids = dataManager.pids;

	start = time(NULL);

	unordered_map<int, vector<int>> registryPidToInds;
	for (int i = 0; i < registry_records.size(); ++i)
		registryPidToInds[registry_records[i].pid].push_back(i);

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
			if (gend == 1) {
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
	int prog_pid = 0;
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < all_pids.size(); ++i) {
		int pid = all_pids[i];
		vector<int> *registry_inds = NULL;
		if (registryPidToInds.find(pid) != registryPidToInds.end()) {
			registry_inds = &registryPidToInds.at(pid);
		}
		else {
#pragma omp atomic
			++prog_pid;
			continue; // show only stats for registry=1,0
		}
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
			bool has_intr = false;
			for (int indReg : *registry_inds)
			{
				MedRegistryRecord &regRec = registry_records.at(indReg);
				if (regRec.registry_value < 0) {
					has_intr = true;//censored out - mark as done, no valid registry records for pid
					break;
				}
				if (regRec.age != -1 && regRec.registry_value > 0)
					ageBin = float(ageBinValue * floor(double(regRec.age) / ageBinValue));
				else
					ageBin = float(ageBinValue * floor(double(medial::repository::DateDiff(BDate, sigRec.start_date)) / ageBinValue));
				ageBin_index = int((ageBin - min_age) / ageBinValue);
				if (ageBin < min_age || ageBin > max_age)
					continue; //skip out of range...
				/*if (ageBin < min_age)
					ageBin_index = 0;
				if (ageBin > max_age)
					ageBin_index = age_bin_count - 1;*/
				pos = 2;

				bool intersect = date_intersection(regRec.min_allowed_date, regRec.max_allowed_date,
					regRec.start_date, regRec.end_date, sigRec.start_date,
					time_window_from, time_window_to, regRec.registry_value <= 0);
				if (intersect) {
					has_intr = true;
					//pos += 1; //registry_value > 0 - otherwise skip this
					pos += int(regRec.registry_value > 0);
					if (gender == 1)
						update_loop(pos, ageBin_index, ageBin, sigRec, maleSignalToStats, val_seen_pid_pos);
					else
						update_loop(pos, ageBin_index, ageBin, sigRec, femaleSignalToStats, val_seen_pid_pos);
					if (!debug_file.empty() && debug_vals.find(sigRec.registry_value) != debug_vals.end()) {
#pragma omp critical
						dbg_file << pid << "\t" << sigRec.start_date << "\t" << sigRec.registry_value
							<< "\t" << regRec.start_date << "\t" << regRec.end_date
							<< "\t" << gender << "\t" << ageBin << "\t" << regRec.registry_value
							<< "\n";
					}
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
	if (!debug_file.empty())
		dbg_file.close();
	unordered_set<float> vals;
	for (auto it = maleSignalToStats.begin(); it != maleSignalToStats.end(); ++it)
		vals.insert(it->first);
	for (auto it = femaleSignalToStats.begin(); it != femaleSignalToStats.end(); ++it)
		vals.insert(it->first);

	//update values prevalence
	for (auto it = vals.begin(); it != vals.end(); ++it)
	{
		if (maleSignalToStats.find(*it) != maleSignalToStats.end())
			for (auto jt = maleSignalToStats[*it].begin(); jt != maleSignalToStats[*it].end(); ++jt) {
				if (male_total_prevalence.find(jt->first) == male_total_prevalence.end())
					MTHROW_AND_ERR("Sample is too small - no incidences for age_bin=%f in males\n", jt->first);
				/*if (male_total_prevalence[jt->first][1] < maleSignalToStats[*it][jt->first][3])
					MTHROW_AND_ERR("Bug in sampling by age or registry - not found[1,3]. Males, age_bin=%f"
						, jt->first);
				if (male_total_prevalence[jt->first][0] < maleSignalToStats[*it][jt->first][2])
					MTHROW_AND_ERR("Bug in sampling by age or registry - not found[0,2]. Males, age_bin=%f"
						, jt->first);*/
				maleSignalToStats[*it][jt->first][0] = male_total_prevalence[jt->first][0] - maleSignalToStats[*it][jt->first][2];
				maleSignalToStats[*it][jt->first][1] = male_total_prevalence[jt->first][1] - maleSignalToStats[*it][jt->first][3];
				if (maleSignalToStats[*it][jt->first][0] < 0)
					maleSignalToStats[*it][jt->first][0] = 0;
				if (maleSignalToStats[*it][jt->first][1] < 0)
					maleSignalToStats[*it][jt->first][1] = 0;

			}
		if (femaleSignalToStats.find(*it) != femaleSignalToStats.end())
			for (auto jt = femaleSignalToStats[*it].begin(); jt != femaleSignalToStats[*it].end(); ++jt) {
				if (female_total_prevalence.find(jt->first) == female_total_prevalence.end())
					MTHROW_AND_ERR("Sample is too small - no incidences for age_bin=%f in females\n", jt->first);
				/*if (female_total_prevalence[jt->first][1] < femaleSignalToStats[*it][jt->first][3])
					MTHROW_AND_ERR("Bug in sampling by age or registry - not found[1,3]. Females, age_bin=%f"
						, jt->first);
				if (female_total_prevalence[jt->first][0] < femaleSignalToStats[*it][jt->first][2])
					MTHROW_AND_ERR("Bug in sampling by age or registry - not found[0,2]. Females, age_bin=%f"
						, jt->first);*/
				femaleSignalToStats[*it][jt->first][0] = female_total_prevalence[jt->first][0] - femaleSignalToStats[*it][jt->first][2];
				femaleSignalToStats[*it][jt->first][1] = female_total_prevalence[jt->first][1] - femaleSignalToStats[*it][jt->first][3];
				if (femaleSignalToStats[*it][jt->first][0] < 0)
					femaleSignalToStats[*it][jt->first][0] = 0;
				if (femaleSignalToStats[*it][jt->first][1] < 0)
					femaleSignalToStats[*it][jt->first][1] = 0;
			}
	}

	duration = (int)difftime(time(NULL), start);
	MLOG("Finished in %d seconds with %d records in males and %d records in females\n",
		duration, (int)maleSignalToStats.size(), (int)femaleSignalToStats.size());

}

void MedRegistry::get_pids(vector<int> &pids) {
	pids.clear();
	unordered_set<int> seen_pid;
	for (size_t i = 0; i < registry_records.size(); ++i)
		seen_pid.insert(registry_records[i].pid);
	pids.insert(pids.end(), seen_pid.begin(), seen_pid.end());
}

void MedRegistry::create_incidence_file(const string &file_path, const string &rep_path,
	int age_bin, int min_age, int max_age, int time_period, bool use_kaplan_meir,
	const string &sampler_name, const string &sampler_args) {
	MedSamplingStrategy *sampler = MedSamplingStrategy::make_sampler(sampler_name, sampler_args);
	MedSamples incidence_samples;
	MLOG("Sampling for incidence stats...\n");
	sampler->do_sample(registry_records, incidence_samples);
	MLOG("Done...\n");
	delete sampler;
	MedRepository rep;
	vector<int> pids;
	get_pids(pids);
	vector<string> signal_to_read = { "BYEAR", "GENDER" };
	if (rep.read_all(rep_path, pids, signal_to_read) < 0)
		MTHROW_AND_ERR("FAILED reading repository %s\n", rep_path.c_str());
	min_age = int(min_age / age_bin) * age_bin;

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

inline void init_list(const string &reg_path, vector<bool> &list) {
	list.resize(16000000);
	ifstream file;
	file.open(reg_path);
	if (!file.is_open())
		MTHROW_AND_ERR("Unable to open test indexes file\n%s", reg_path.c_str());
	string line;
	//getline(file, line); //ignore first line
	while (getline(file, line)) {
		boost::trim(line);
		if (line.empty())
			continue;
		if (line.at(0) == '#')
			continue;
		list[stoi(line) - 5000000] = true;
	}
	file.close();
}

RegistrySignalSet::RegistrySignalSet(const string &sigName, int durr_time, int buffer_time, bool take_first,
	MedRepository &rep, const vector<string> &sets) {
	signalName = sigName;
	buffer_duration = buffer_time;
	duration_flag = durr_time;
	take_only_first = take_first;
	if (!sets.empty()) {
		int section_id = rep.dict.section_id(sigName);
		rep.dict.curr_section = section_id;
		rep.dict.default_section = section_id;
		rep.dict.prep_sets_lookup_table(section_id, sets, Flags);
	}
}

float RegistrySignalSet::get_outcome(UniversalSigVec &s, int current_i) {
	if (current_i < 0 || current_i >= s.len
		|| s.Val(current_i) < 0 || s.Val(current_i) >= Flags.size())
		return 0;
	return Flags[(int)s.Val(current_i)];
}

int RegistrySignalSet::init(map<string, string>& map) {

	string sets_arg = "";
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		//! [RegistrySignalSet::init]
		if (it->first == "signalName")
			signalName = it->second;
		else if (it->first == "duration_flag")
			duration_flag = stoi(it->second);
		else if (it->first == "buffer_duration")
			buffer_duration = stoi(it->second);
		else if (it->first == "take_only_first")
			take_only_first = stoi(it->second) > 0;
		else if (it->first == "sets")
			sets_arg = it->second;
		else
			MTHROW_AND_ERR("unsupported element \"%s\"\n",
				it->first.c_str());
		//! [RegistrySignalSet::init]
	}
	//save to end
	if (!sets_arg.empty()) {
		std::map<string, string> set_init;
		//rep=;sets=;
		init_map_from_string(sets_arg, set_init);
		if (set_init.find("rep") == set_init.end())
			MTHROW_AND_ERR("sets should contain rep\n");
		if (set_init.find("sets") == set_init.end())
			MTHROW_AND_ERR("sets should contain sets for file path\n");
		vector<string> sets;
		medial::io::read_codes_file(set_init["sets"], sets);
		MedRepository rep;
		if (rep.init(set_init["rep"]) < 0)
			MTHROW_AND_ERR("unabale to init rep by %s\n", set_init["rep"].c_str());
		if (!sets.empty()) {
			int section_id = rep.dict.section_id(signalName);
			rep.dict.curr_section = section_id;
			rep.dict.default_section = section_id;
			rep.dict.prep_sets_lookup_table(section_id, sets, Flags);
		}
	}
	return 0;
}

RegistrySignalSet::RegistrySignalSet(const string &init_string, MedRepository &rep, const vector<string> &sets) {
	init_from_string(init_string);
	if (!sets.empty()) {
		int section_id = rep.dict.section_id(signalName);
		rep.dict.curr_section = section_id;
		rep.dict.default_section = section_id;
		rep.dict.prep_sets_lookup_table(section_id, sets, Flags);
	}
}

void MedRegistryCodesList::init(MedRepository &rep, int start_dur, int end_durr, int max_repo,
	const vector<RegistrySignal *> signal_conditions, const string &skip_pid_file,
	const unordered_map<int, int> *pid_to_censor_dates) {
	if (signal_conditions.empty())
		MTHROW_AND_ERR("must be initialize with something\n");
	init_called = true;
	start_buffer_duration = start_dur;
	end_buffer_duration = end_durr;
	max_repo_date = max_repo;
	if (!skip_pid_file.empty())
		init_list(skip_pid_file, SkipPids);
	if (pid_to_censor_dates != NULL)
		pid_to_max_allowed = *pid_to_censor_dates;
	signalCodes.clear();
	for (size_t i = 0; i < signal_conditions.size(); ++i)
		signalCodes.push_back(rep.sigs.sid(signal_conditions[i]->signalName));
	//the user called init for this signal_conditions
	signal_filters = signal_conditions;
}

RegistrySignalRange::RegistrySignalRange(const string &sigName, int durr_time, int buffer_time,
	bool take_first, float min_range, float max_range) {
	signalName = sigName;
	duration_flag = durr_time;
	buffer_duration = buffer_time;
	take_only_first = take_first;

	min_value = min_range;
	max_value = max_range;
}

float RegistrySignalRange::get_outcome(UniversalSigVec &s, int current_i) {
	return current_i < s.len && s.Val(current_i) >= min_value && s.Val(current_i) <= max_value;
}

int RegistrySignalRange::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		//! [RegistrySignalRange::init]
		if (it->first == "signalName")
			signalName = it->second;
		else if (it->first == "duration_flag")
			duration_flag = stoi(it->second);
		else if (it->first == "buffer_duration")
			buffer_duration = stoi(it->second);
		else if (it->first == "take_only_first")
			take_only_first = stoi(it->second) > 0;
		else if (it->first == "min_value")
			min_value = stof(it->second);
		else if (it->first == "max_value")
			max_value = stof(it->second);
		else
			MTHROW_AND_ERR("unsupported element \"%s\"\n",
				it->first.c_str());
		//! [RegistrySignalRange::init]
	}
	return 0;
}

inline int Date_wrapper(UniversalSigVec &signal, int i) {
	if (signal.get_type() != T_Value)
		return signal.Date(i);
	else
		return (int)signal.Val(i);
}

int medial::repository::fetch_next_date(vector<UniversalSigVec> &patientFile, vector<int> &signalPointers) {
	int minDate = -1, minDate_index = -1;
	for (size_t i = 0; i < patientFile.size(); ++i)
	{
		UniversalSigVec &data = patientFile[i];
		if (signalPointers[i] >= data.len)
			continue; //already reached the end for this signal
		if (data.get_type() == T_Value) {
			if (minDate_index == -1 || data.Val(signalPointers[i]) < minDate) {
				minDate = (int)data.Val(signalPointers[i]);
				minDate_index = (int)i;
			}
		}
		else {
			if (minDate_index == -1 || data.Date(signalPointers[i]) < minDate) {
				minDate = data.Date(signalPointers[i]);
				minDate_index = (int)i;
			}
		}
	}
	if (minDate_index >= 0)
		++signalPointers[minDate_index];
	return minDate_index;
}

void MedRegistryCodesList::get_registry_records(int pid,
	int bdate, vector<UniversalSigVec> &usv, vector<MedRegistryRecord> &results) {
	if (!init_called)
		MTHROW_AND_ERR("Must be initialized by init before use\n");
	vector<int> signals_indexes_pointers(signal_filters.size()); //all in 0

	int max_allowed_date = max_repo_date;
	if (pid_to_max_allowed.find(pid) != pid_to_max_allowed.end())
		max_allowed_date = pid_to_max_allowed[pid];

	MedRegistryRecord r;
	r.pid = pid;
	int start_date = -1, last_date = -1;
	int signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
	while (signal_index >= 0)
	{
		UniversalSigVec &signal = usv[signal_index];
		RegistrySignal *signal_prop = signal_filters[signal_index];
		int i = signals_indexes_pointers[signal_index] - 1; //the current signal time
		//find first date if not marked already
		if (start_date == -1) {
			if (Date_wrapper(signal, i) >= bdate) {
				start_date = Date_wrapper(signal, i);
				r.start_date = medial::repository::DateAdd(start_date, start_buffer_duration);
			}
			else {
				signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
				continue;
			}
		}
		int min_date = medial::repository::DateAdd(start_date, start_buffer_duration);
		r.min_allowed_date = min_date; //at least 1 year data
		r.age = int(medial::repository::DateDiff(bdate, r.start_date));
		r.registry_value = 0;
		//I have start_date
		if (Date_wrapper(signal, i) > max_allowed_date)
			break;
		last_date = Date_wrapper(signal, i);

		if (signal_prop->get_outcome(signal, i) > 0) {
			//flush buffer
			int last_date_c = medial::repository::DateAdd(Date_wrapper(signal, i), -signal_prop->buffer_duration);
			r.end_date = last_date_c;
			r.max_allowed_date = last_date_c;
			if (r.end_date > r.start_date)
				results.push_back(r);

			//start new record
			//r.pid = pid;
			r.min_allowed_date = min_date;
			r.max_allowed_date = Date_wrapper(signal, i);
			r.start_date = Date_wrapper(signal, i);
			r.age = (int)medial::repository::DateDiff(bdate, Date_wrapper(signal, i));
			r.registry_value = 1;
			if (signal_prop->take_only_first) {
				r.end_date = 30000000;
				results.push_back(r);
				return;
			}
			else
				r.end_date = medial::repository::DateAdd(Date_wrapper(signal, i), signal_prop->duration_flag);
			int max_search = medial::repository::DateAdd(r.end_date,
				signal_prop->buffer_duration - 1);
			//advanced till passed end_date + buffer with no reapeating RC:
			while (signal_index >= 0 && Date_wrapper(signal, i) < max_search) {
				if (signal_prop->get_outcome(signal, i) > 0) {
					r.end_date = medial::repository::DateAdd(Date_wrapper(signal, i), signal_prop->duration_flag);
					max_search = medial::repository::DateAdd(r.end_date, signal_prop->buffer_duration - 1);
				}

				signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
				if (signal_index < 0)
					break;
				i = signals_indexes_pointers[signal_index] - 1; //the current signal time
				signal_prop = signal_filters[signal_index];
				signal = usv[signal_index];
			}
			results.push_back(r);
			if (signal_index < 0) {
				r.start_date = 30000000; //no more control times, reached the end
				break;
			}
			//prepare for next:
			r.min_allowed_date = min_date;
			r.registry_value = 0;
			r.start_date = Date_wrapper(signal, i); //already after duration and buffer. can start new control
			r.age = int(medial::repository::DateDiff(bdate, r.start_date));
			continue; //dont call fetch_next again
		}

		signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
	}

	r.end_date = last_date;
	last_date = medial::repository::DateAdd(last_date, -end_buffer_duration);
	r.max_allowed_date = last_date;
	if (r.end_date > r.start_date)
		results.push_back(r);
}

double medial::contingency_tables::calc_chi_square_dist(const map<float, vector<int>> &gender_sorted,
	int smooth_balls, float allowed_error, int minimal_balls) {
	//calc over all ages
	double regScore = 0;
	for (auto i = gender_sorted.begin(); i != gender_sorted.end(); ++i) { //iterate over age bins
		const vector<int> &probs_tmp = i->second; //the forth numbers
		if (!(probs_tmp[0] >= minimal_balls && probs_tmp[1] >= minimal_balls
			&& probs_tmp[2] >= minimal_balls && probs_tmp[3] >= minimal_balls))
			continue; //skip row with minimal balls
		vector<double> probs((int)probs_tmp.size()); //the forth numbers - float with fix
		double totCnt = 0;
		vector<double> R(2);
		vector<double> C(2);

		C[0] = probs_tmp[0] + probs_tmp[2]; //how much controls
		C[1] = probs_tmp[1] + probs_tmp[1 + 2]; //how much cases
		totCnt = C[0] + C[1];
		for (size_t j = 0; j < probs_tmp.size(); ++j)
			probs[j] = probs_tmp[j] + (smooth_balls * C[j % 2] / totCnt);  /* add smooth addition */

		totCnt = 0;
		R[0] = probs[0] + probs[1];
		R[1] = probs[2 + 0] + probs[2 + 1];
		C[0] = probs[0] + probs[2]; //how much controls
		C[1] = probs[1] + probs[1 + 2]; //how much cases
		for (size_t j = 0; j < probs.size(); ++j)
			totCnt += probs[j];

		for (size_t j = 0; j < probs.size(); ++j)
		{
			double	Qij = probs[j];
			double Eij = (R[j / 2] * C[j % 2]) / totCnt;
			double Dij = abs(Qij - Eij) - (allowed_error / 100) * Eij;
			if (Dij < 0)
				Dij = 0;

			if (Eij > 0)
				regScore += (Dij * Dij) / (Eij); //Chi-square
		}

	}
	return regScore;
}

double medial::contingency_tables::calc_mcnemar_square_dist(const map<float, vector<int>> &gender_sorted) {
	//calc over all ages
	double regScore = 0;
	for (auto i = gender_sorted.begin(); i != gender_sorted.end(); ++i) { //iterate over age bins
		const vector<int> &counts = i->second; //the forth numbers
		double totCnt = 0;
		vector<double> R(2);
		R[0] = counts[0] + counts[1];
		R[1] = counts[2 + 0] + counts[2 + 1];
		totCnt = R[0] + R[1];

		double p_min;
		double b, c;
		int min_ind = 0;
		//Mathcing to lower:
		p_min = R[0];
		if (R[1] < p_min) {
			p_min = R[1];
			++min_ind;
		}
		if (p_min = 0)
			continue; //no matching possible 0's in both cells
		if (min_ind == 0) {
			//matching to first whos lower(no sig appear before is lower):
			b = counts[0 + 1];
			c = counts[2 + 1] * R[0] / R[1];
		}
		else {
			//matching to second whos lower:
			b = counts[2 + 1];
			c = counts[0 + 1] * R[1] / R[0];
		}

		if (b + c > 0)
			regScore += (b - c) * (b - c) / (b + c); //Mcnemar
	}
	return regScore;
}

double medial::contingency_tables::chisqr(int Dof, double Cv)
{
	if (Dof < 1 || Cv <= 0) {
		return 1;
	}
	boost::math::chi_squared dist(Dof);
	return (1.0 - boost::math::cdf(dist, Cv));
}

int _count_legal_rows(const  map<float, vector<int>> &m, int minimal_balls) {
	int res = 0;
	for (auto it = m.begin(); it != m.end(); ++it)
	{
		int ind = 0;
		bool all_good = true;
		while (all_good && ind < it->second.size()) {
			all_good = it->second[ind] >= minimal_balls;
			++ind;
		}
		res += int(all_good);
	}
	return res;
}

void medial::contingency_tables::calc_chi_scores(const map<float, map<float, vector<int>>> &male_stats,
	const map<float, map<float, vector<int>>> &female_stats,
	vector<float> &all_signal_values, vector<int> &signal_indexes,
	vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
	, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio, int smooth_balls
	, float allowed_error, int minimal_balls) {

	unordered_set<float> all_vals;
	for (auto i = male_stats.begin(); i != male_stats.end(); ++i)
		all_vals.insert(i->first);
	for (auto i = female_stats.begin(); i != female_stats.end(); ++i)
		all_vals.insert(i->first);
	all_signal_values.insert(all_signal_values.end(), all_vals.begin(), all_vals.end());
	signal_indexes.resize(all_signal_values.size());
	for (size_t i = 0; i < signal_indexes.size(); ++i)
		signal_indexes[i] = (int)i;
	posCnts.resize(all_signal_values.size());
	valCnts.resize(all_signal_values.size());
	lift.resize(all_signal_values.size());
	scores.resize(all_signal_values.size());
	p_values.resize(all_signal_values.size());
	pos_ratio.resize(all_signal_values.size());

	for (int index : signal_indexes)
	{
		float signalVal = all_signal_values[index];
		//check chi-square for this value:
		double totCnt = 0;
		double sig_sum = 0;
		double sum_noSig_reg = 0;
		double sum_noSig_tot = 0;

		if (male_stats.find(signalVal) != male_stats.end())
			for (auto jt = male_stats.at(signalVal).begin(); jt != male_stats.at(signalVal).end(); ++jt) {
				totCnt += jt->second[2] + jt->second[3];
				posCnts[index] += jt->second[1 + 2];
				sig_sum += jt->second[0 + 2];
				sum_noSig_reg += jt->second[1 + 0];
				sum_noSig_tot += jt->second[1 + 0] + jt->second[0 + 0];
			}
		if (female_stats.find(signalVal) != female_stats.end())
			for (auto jt = female_stats.at(signalVal).begin(); jt != female_stats.at(signalVal).end(); ++jt) {
				totCnt += jt->second[2] + jt->second[3];
				posCnts[index] += jt->second[1 + 2];
				sig_sum += jt->second[0 + 2];
				sum_noSig_reg += jt->second[1 + 0];
				sum_noSig_tot += jt->second[1 + 0] + jt->second[0 + 0];
			}
		if (totCnt == 0)
			continue;
		valCnts[index] = totCnt; //for signal apeareance
		sig_sum += posCnts[index];
		if (sig_sum > 0 && sum_noSig_reg > 0)
			lift[index] = (posCnts[index] / sig_sum) / (sum_noSig_reg / sum_noSig_tot);
		if (sig_sum > 0 && sum_noSig_reg <= 0)
			lift[index] = 2 * posCnts[index]; //maximum lift
		pos_ratio[index] = posCnts[index] / totCnt;

		double regScore = 0;
		if (male_stats.find(signalVal) != male_stats.end())
			regScore += calc_chi_square_dist(male_stats.at(signalVal), smooth_balls, allowed_error, minimal_balls); //Males
		if (female_stats.find(signalVal) != female_stats.end())
			regScore += calc_chi_square_dist(female_stats.at(signalVal), smooth_balls, allowed_error, minimal_balls); //Females

		scores[index] = (float)regScore;
		int dof = -1;
		if (male_stats.find(signalVal) != male_stats.end())
			dof += _count_legal_rows(male_stats.at(signalVal), minimal_balls);
		if (female_stats.find(signalVal) != female_stats.end())
			dof += _count_legal_rows(female_stats.at(signalVal), minimal_balls);
		double pv = chisqr(dof, regScore);
		p_values[index] = pv;
	}
}

void medial::contingency_tables::calc_mcnemar_scores(const map<float, map<float, vector<int>>> &male_stats,
	const map<float, map<float, vector<int>>> &female_stats,
	vector<float> &all_signal_values, vector<int> &signal_indexes,
	vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
	, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio) {

	unordered_set<float> all_vals;
	for (auto i = male_stats.begin(); i != male_stats.end(); ++i)
		all_vals.insert(i->first);
	for (auto i = female_stats.begin(); i != female_stats.end(); ++i)
		all_vals.insert(i->first);
	all_signal_values.insert(all_signal_values.end(), all_vals.begin(), all_vals.end());
	signal_indexes.resize(all_signal_values.size());
	for (size_t i = 0; i < signal_indexes.size(); ++i)
		signal_indexes[i] = (int)i;
	posCnts.resize(all_signal_values.size());
	valCnts.resize(all_signal_values.size());
	lift.resize(all_signal_values.size());
	scores.resize(all_signal_values.size());
	p_values.resize(all_signal_values.size());
	pos_ratio.resize(all_signal_values.size());

	for (int index : signal_indexes)
	{
		float signalVal = all_signal_values[index];
		//check chi-square for this value:
		double totCnt = 0;
		double sig_sum = 0;
		double sum_noSig_reg = 0;
		double sum_noSig_tot = 0;

		if (male_stats.find(signalVal) != male_stats.end())
			for (auto jt = male_stats.at(signalVal).begin(); jt != male_stats.at(signalVal).end(); ++jt) {
				totCnt += jt->second[2] + jt->second[3];
				posCnts[index] += jt->second[1 + 2];
				sig_sum += jt->second[0 + 2];
				sum_noSig_reg += jt->second[1 + 0];
				sum_noSig_tot += jt->second[1 + 0] + jt->second[0 + 0];
			}
		if (female_stats.find(signalVal) != female_stats.end())
			for (auto jt = female_stats.at(signalVal).begin(); jt != female_stats.at(signalVal).end(); ++jt) {
				totCnt += jt->second[2] + jt->second[3];
				posCnts[index] += jt->second[1 + 2];
				sig_sum += jt->second[0 + 2];
				sum_noSig_reg += jt->second[1 + 0];
				sum_noSig_tot += jt->second[1 + 0] + jt->second[0 + 0];
			}
		if (totCnt == 0)
			continue;
		valCnts[index] = totCnt; //for signal apeareance
		sig_sum += posCnts[index];
		if (sig_sum > 0 && sum_noSig_reg > 0)
			lift[index] = (posCnts[index] / sig_sum) / (sum_noSig_reg / sum_noSig_tot);
		if (sig_sum > 0 && sum_noSig_reg <= 0)
			lift[index] = 2 * posCnts[index]; //maximum lift
		pos_ratio[index] = posCnts[index] / totCnt;

		double regScore = 0;
		if (male_stats.find(signalVal) != male_stats.end())
			regScore += calc_mcnemar_square_dist(male_stats.at(signalVal)); //Males
		if (female_stats.find(signalVal) != female_stats.end())
			regScore += calc_mcnemar_square_dist(female_stats.at(signalVal)); //Females

		scores[index] = (float)regScore;
		int dof = -1;
		if (male_stats.find(signalVal) != male_stats.end())
			dof += _count_legal_rows(male_stats.at(signalVal), 0);
		if (female_stats.find(signalVal) != female_stats.end())
			dof += _count_legal_rows(female_stats.at(signalVal), 0);
		double pv = chisqr(dof, regScore);
		p_values[index] = pv;
	}
}

void medial::contingency_tables::FilterRange(vector<int> &indexes, const vector<double> &vecCnts
	, double min_val, double max_val) {
	vector<int> filtered_indexes;
	filtered_indexes.reserve(indexes.size());
	for (size_t i = 0; i < indexes.size(); ++i)
		if (vecCnts[indexes[i]] >= min_val && vecCnts[indexes[i]] <= max_val)
			filtered_indexes.push_back(indexes[i]);
	indexes.swap(filtered_indexes);
}

void read_gender_val(ifstream &fr, map<float, vector<int>> &vec) {
	int dict_size;
	fr.read((char *)&dict_size, sizeof(int));
	for (size_t i = 0; i < dict_size; ++i)
	{
		float ageBin;
		fr.read((char *)&ageBin, sizeof(float));
		vec[ageBin] = vector<int>(4);
		for (size_t j = 0; j < vec[ageBin].size(); ++j)
		{
			fr.read((char *)&vec[ageBin][j], sizeof(int));
		}
	}
}

void write_gender_val(ofstream &fw, const map<float, vector<int>> &gender_stats) {
	int dict_size = (int)gender_stats.size();
	fw.write((char *)&dict_size, sizeof(int));

	for (auto it = gender_stats.begin(); it != gender_stats.end(); ++it)
	{
		float ageBin = it->first;
		fw.write((char*)&ageBin, sizeof(float));
		if (it->second.size() != 4) {
			throw logic_error("validation failed for stats vector of size 4");
		}
		for (size_t i = 0; i < it->second.size(); ++i) //fixed size - 4
		{
			int num = it->second[i];
			fw.write((char *)&num, sizeof(int));
		}
	}
}

void write_gender(const map<float, map<float, vector<int>>> &dict, const string &file_path) {
	ofstream fw(file_path, ios::binary);
	//Format is dictionary Size then each rowL float, male_vector, feamle_vector. gender_vector = map_size, for each row: float, 4 int vector numbers for stats
	int dict_size = (int)dict.size();
	fw.write((char *)&dict_size, sizeof(int));
	for (auto it = dict.begin(); it != dict.end(); ++it)
	{
		float signalValue = it->first;
		const map<float, vector<int>> &stats = it->second;

		fw.write((char *)&signalValue, sizeof(float));
		write_gender_val(fw, stats);
		//lets serialize male and then female:
	}

	fw.close();
}

void read_gender(const string &file_path, map<float, map<float, vector<int>>> &dict) {
	ifstream fr(file_path, ios::binary);
	int dict_size;
	fr.read((char *)&dict_size, sizeof(int));

	for (size_t i = 0; i < dict_size; ++i)
	{
		float signalValue;
		fr.read((char *)&signalValue, sizeof(float));
		read_gender_val(fr, dict[signalValue]);
	}

	fr.close();
}

void medial::contingency_tables::write_stats(const string &file_path,
	const map<float, map<float, vector<int>>> &males_stats, const map<float, map<float, vector<int>>> &females_stats) {

	ofstream fw(file_path, ios::binary);
	if (!fw.good())
		MTHROW_AND_ERR("IOError: can't open %s for writing.\n", file_path.c_str());
	//Format is dictionary Size then each rowL float, male_vector, feamle_vector. gender_vector = map_size, for each row: float, 4 int vector numbers for stats
	int dict_size = (int)males_stats.size();
	fw.write((char *)&dict_size, sizeof(int));
	for (auto it = males_stats.begin(); it != males_stats.end(); ++it)
	{
		float signalValue = it->first;
		const map<float, vector<int>> &stats = it->second;

		fw.write((char *)&signalValue, sizeof(float));
		write_gender_val(fw, stats);
	}

	dict_size = (int)females_stats.size();
	fw.write((char *)&dict_size, sizeof(int));
	for (auto it = females_stats.begin(); it != females_stats.end(); ++it)
	{
		float signalValue = it->first;
		const map<float, vector<int>> &stats = it->second;

		fw.write((char *)&signalValue, sizeof(float));
		write_gender_val(fw, stats);
	}

	fw.close();
	MLOG("wrote [%d] keys on both male and female.\n", int(males_stats.size() + females_stats.size()));
}

void medial::contingency_tables::read_stats(const string &file_path,
	map<float, map<float, vector<int>>> &males_stats, map<float, map<float, vector<int>>> &females_stats) {
	ifstream fr(file_path, ios::binary);
	int dict_size;
	fr.read((char *)&dict_size, sizeof(int));
	for (size_t i = 0; i < dict_size; ++i)
	{
		float signalValue;
		fr.read((char *)&signalValue, sizeof(float));
		read_gender_val(fr, males_stats[signalValue]);
	}
	fr.read((char *)&dict_size, sizeof(int));
	for (size_t i = 0; i < dict_size; ++i)
	{
		float signalValue;
		fr.read((char *)&signalValue, sizeof(float));
		read_gender_val(fr, females_stats[signalValue]);
	}

	fr.close();
	MLOG("read [%d] records on both male and female stats.\n",
		int(males_stats.size() + females_stats.size()));
}

void medial::contingency_tables::FilterFDR(vector<int> &indexes,
	const vector<double> &scores, const vector<double> &p_vals, const vector<double> &lift,
	double filter_pval) {
	//sort by  pVal (if equal than -score (Floating point round and they are almost all has same dof)) also use posCnts/ valCnts:
	int num_of_init_test = (int)indexes.size();
	double cm = 0;
	for (size_t i = 0; i < num_of_init_test; ++i)
		cm += 1 / (i + double(1));

	double num_test_factor = num_of_init_test * cm;//dependence correction
	vector<pair<int, vector<double>>> keysSorted(indexes.size());

	for (int i = 0; i < indexes.size(); ++i) {
		vector<double> vec = { p_vals[indexes[i]],
			-lift[indexes[i]] , -scores[indexes[i]] };
		keysSorted[i] = pair<int, vector<double>>(indexes[i], vec);
	}

	sort(keysSorted.begin(), keysSorted.end(), [](pair<int, vector<double>> a, pair<int, vector<double>> b) {
		int pos = 0;
		while (pos < a.second.size() &&
			a.second[pos] == b.second[pos])
			++pos;
		if (pos >= a.second.size())
			return false;
		return b.second[pos] > a.second[pos];
	});

	double normAlpha = filter_pval / num_test_factor;
	int totSum = 0;
	int stop_index = 0;
	while (stop_index < keysSorted.size() && normAlpha * (stop_index + 1) >= keysSorted[stop_index].second[0])
		++stop_index;

	//Keep only filtered indexes
	indexes.resize(stop_index);
	for (size_t i = 0; i < stop_index; ++i)
		indexes[i] = keysSorted[i].first;
}

void medial::print::print_reg_stats(const vector<MedRegistryRecord> &regRecords, const string &log_file) {
	map<float, int> histCounts, histCounts_All;
	vector<unordered_set<int>> pid_index(2);
	for (size_t k = 0; k < regRecords.size(); ++k)
	{
		if (pid_index[regRecords[k].registry_value > 0].find(regRecords[k].pid) == pid_index[regRecords[k].registry_value > 0].end()) {
			if (histCounts.find(regRecords[k].registry_value) == histCounts.end()) {
				histCounts[regRecords[k].registry_value] = 0;
			}
			++histCounts[regRecords[k].registry_value];
		}
		pid_index[regRecords[k].registry_value > 0].insert(regRecords[k].pid);
		++histCounts_All[regRecords[k].registry_value];
	}
	cout << "Registry has " << regRecords.size() << " records. [";
	string delim = " ";
	if (histCounts.size() > 2)
		delim = "\n";
	for (auto it = histCounts.begin(); it != histCounts.end(); ++it)
		cout << delim << (int)it->first << "=" << it->second;
	cout << " ]" << " All = [ ";
	for (auto it = histCounts_All.begin(); it != histCounts_All.end(); ++it)
		cout << delim << (int)it->first << "=" << it->second;
	cout << " ]" << endl;
	if (!log_file.empty()) {
		ofstream fo(log_file, ios::app);
		fo << "Registry has " << regRecords.size() << " records. [";
		if (histCounts.size() > 2)
			delim = "\n";
		for (auto it = histCounts.begin(); it != histCounts.end(); ++it)
			fo << delim << (int)it->first << "=" << it->second;
		fo << " ]" << " All = [ ";
		for (auto it = histCounts_All.begin(); it != histCounts_All.end(); ++it)
			fo << delim << (int)it->first << "=" << it->second;
		fo << " ]" << endl;
		fo.close();
	}
}

void medial::io::read_codes_file(const string &file_path, vector<string> &tokens) {
	tokens.clear();
	ifstream file;
	file.open(file_path);
	if (!file.is_open())
		throw logic_error("Unable to open test indexes file:\n" + file_path);
	string line;
	//getline(file, line); //ignore first line
	while (getline(file, line)) {
		boost::trim(line);
		if (line.empty())
			continue;
		if (line.at(0) == '#')
			continue;
		if (line.find('\t') != string::npos)
			line = line.substr(0, line.find('\t'));
		tokens.push_back(line);
	}
	file.close();
}

RegistrySignal *RegistrySignal::make_registry_signal(const string &type) {
	MedRepository rep;
	vector<string> empty_arr;
	string empty_str = "";
	if (type == "set")
		return new RegistrySignalSet(empty_str, 0, 0, false, rep, empty_arr);
	else if (type == "range")
		return new RegistrySignalRange(empty_str, 0, 0, false, 0, 0);
	else
		MTHROW_AND_ERR("Unsupported type %s\n", type.c_str());
}

RegistrySignal *RegistrySignal::make_registry_signal(const string &type, const string &init_string) {
	RegistrySignal *reg = make_registry_signal(type);
	reg->init_from_string(init_string);
	return reg;
}

void MedRegistryCodesList::parse_registry_rules(const string &reg_cfg, const string &rep_path,
	vector<RegistrySignal *> &result) {
	ifstream fr(reg_cfg);
	if (!fr.good())
		MTHROW_AND_ERR("IOError: can't read registry rules config from %s\n", reg_cfg.c_str());
	string line;
	result.clear();
	while (getline(fr, line))
	{
		boost::trim(line);
		if (line.empty() || line[0] == '#')
			continue; //skip line
		vector<string> tokens;
		boost::split(tokens, line, boost::is_any_of("\t"));
		if (tokens.size() != 2)
			MTHROW_AND_ERR("IOERROR: bad file format excepting one [TAB]. got:\n%s\n", line.c_str());
		string type = tokens[0];
		string ini = tokens[1];
		boost::replace_all(ini, "$REP", rep_path);
		RegistrySignal *reg = RegistrySignal::make_registry_signal(type, ini);
		result.push_back(reg);
	}
	fr.close();
}

int MedRegistryCodesList::init(map<string, string>& map) {
	MedRepository repo;
	string rep_path = "";
	string registry_file_path = "";
	for (auto it = map.begin(); it != map.end(); ++it)
		//! [MedRegistryCodesList::init]
		if (it->first == "rep")
			rep_path = it->second;
		else if (it->first == "start_buffer_duration")
			start_buffer_duration = stoi(it->second);
		else if (it->first == "end_buffer_duration")
			end_buffer_duration = stoi(it->second);
		else if (it->first == "max_repo_date")
			max_repo_date = stoi(it->second);
		else if (it->first == "pid_to_censor_dates") {
			ifstream fr(it->second);
			if (!fr.good())
				MTHROW_AND_ERR("Error in MedRegistryCodesList::init - unable to open %s for reading.",
					it->second.c_str());
			string line;
			while (getline(fr, line)) {
				vector<string> tokens;
				boost::split(tokens, line, boost::is_any_of("\t"));
				if (tokens.size() != 2)
					MTHROW_AND_ERR("Error in MedRegistryCodesList::init - in parsing pid_to_censor_dates file"
						" - %s excpeting TAB for line:\n%s\n", it->second.c_str(), line.c_str());
				pid_to_max_allowed[stoi(tokens[0])] = stoi(tokens[1]);
			}
			fr.close();
		}
		else if (it->first == "config_signals_rules")
			registry_file_path = it->second;
		else
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - Unsupported init param \"%s\"\n", it->first.c_str());
		//! [MedRegistryCodesList::init]
		if (rep_path.empty())
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - please provide rep param to init function\n");
		if (max_repo_date == 0)
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - please provide max_repo_date param to init function\n");
		if (registry_file_path.empty())
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - please provide config_signals_rules param to init function\n");

		if (repo.init(rep_path) < 0)
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - Unable to init repositrory from path %s\n", rep_path.c_str());

		parse_registry_rules(registry_file_path, rep_path, signal_filters);
		init_called = true;

		signalCodes.clear();
		for (size_t i = 0; i < signal_filters.size(); ++i)
			signalCodes.push_back(repo.sigs.sid(signal_filters[i]->signalName));

		return 0;
}