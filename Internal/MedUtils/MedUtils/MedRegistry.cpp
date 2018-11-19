#include "MedRegistry.h"
#include <fstream>
#include <string>
#include <iostream>
#include <random>
#include <map>
#include "Logger/Logger/Logger.h"
#include <MedUtils/MedUtils/MedUtils.h>
#include <MedProcessTools/MedProcessTools/MedProcessUtils.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <omp.h>

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int read_time(int time_unit, const string &p) {
	int max_date_mark = 30000000;
	if (time_unit != MedTime::Date)
		max_date_mark = 2000000000;
	if (p == to_string(max_date_mark))
		return max_date_mark;

	return med_time_converter.convert_datetime_safe(time_unit, p, 2);
}
string write_time(int time_unit, int time) {
	int max_date_mark = 30000000;
	if (time_unit != MedTime::Date)
		max_date_mark = 2000000000;
	if (time == max_date_mark)
		return to_string(max_date_mark);

	return med_time_converter.convert_times_S(time_unit, MedTime::DateTimeString, time);
}

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
		if (tokens.size() == 2 && tokens[0] == "TIME_UNIT") {
			MLOG("MedRegistry TIME_UNIT=%s\n", tokens[1].c_str());
			time_unit = med_time_converter.string_to_type(tokens[1]);
			continue;
		}

		if (tokens.size() != 7) {
			cerr << "Bad File Format in line " << lineNum << " got \"" << line << "\"" << " parsed " << tokens.size() << " tokens" << endl;
			throw out_of_range("File has bad format");
		}


		MedRegistryRecord pr;
		pr.pid = stoi(tokens[0]);
		pr.start_date = read_time(time_unit, tokens[1]);
		pr.end_date = read_time(time_unit, tokens[2]);
		pr.min_allowed_date = read_time(time_unit, tokens[3]);
		pr.max_allowed_date = read_time(time_unit, tokens[4]);
		pr.registry_value = stof(tokens[5]);
		registry_records.push_back(pr);

	}
	fp.close();
	MLOG("[Read %d registry records from %s]\n", (int)registry_records.size(), file_path.c_str());
}

void MedRegistry::write_text_file(const string &file_path) const {
	char delim = '\t';
	ofstream fw(file_path);
	if (!fw.good())
		MTHROW_AND_ERR("IOError: can't write file %s\n", file_path.c_str());
	fw << "# Format: PID, Start_Date, End_Date, min_allowed_date, max_allowed_date, Age, RegistryValue" << endl;
	fw << "# Created By Script - Insert Comments following #..." << endl;
	fw << "TIME_UNIT" << delim << med_time_converter.type_to_string(time_unit) << endl;
	fw << endl;
	for (size_t i = 0; i < registry_records.size(); ++i)
		fw << registry_records[i].pid << delim <<
		write_time(time_unit, registry_records[i].start_date) << delim <<
		write_time(time_unit, registry_records[i].end_date) << delim
		<< write_time(time_unit, registry_records[i].min_allowed_date) << delim
		<< write_time(time_unit, registry_records[i].max_allowed_date) << delim
		<< delim << registry_records[i].registry_value << "\n";

	fw.flush();
	fw.close();
	MLOG("[Wrote %d registry records to %s]\n", (int)registry_records.size(), file_path.c_str());
}

void MedRegistry::create_registry(MedPidRepository &dataManager, medial::repository::fix_method method, vector<RepProcessor *> *rep_processors) {
	MLOG_D("Creating registry...\n");
	vector<int> used_sigs;

	time_t start = time(NULL);
	time_t last_time_print = start;
	double duration;
	int prog_pid = 0;
	int bDateCode = dataManager.sigs.sid("BDATE");
	//update using rep_processors:
	vector<string> physical_signals;
	vector<string> sig_names_use = medial::repository::prepare_repository(dataManager, signalCodes_names,
		physical_signals, rep_processors);
	vector<int> final_sigs_to_read(sig_names_use.size()), physical_ids(physical_signals.size());
	for (size_t i = 0; i < sig_names_use.size(); ++i) {
		int sid = dataManager.sigs.sid(sig_names_use[i]);
		if (sid < 0)
			MTHROW_AND_ERR("Error in MedRegistry::create_registry - Couldn't find signal %s in repository or virtual\n",
				sig_names_use[i].c_str());
		final_sigs_to_read[i] = sid;
	}
	for (size_t i = 0; i < physical_ids.size(); ++i) {
		int sid = dataManager.sigs.sid(physical_signals[i]);
		if (sid < 0)
			MTHROW_AND_ERR("Error in MedRegistry::create_registry - Couldn't find signal %s in repository or virtual\n",
				physical_signals[i].c_str());
		physical_ids[i] = sid;
	}
	vector<int> signalCodes(signalCodes_names.size());
	for (size_t i = 0; i < signalCodes_names.size(); ++i)
		signalCodes[i] = dataManager.sigs.sid(signalCodes_names[i]);

	if (!dataManager.index.index_table[bDateCode].is_loaded)
		MTHROW_AND_ERR("Error in MedRegistry::create_registry - you haven't loaded BDATE for repository which is needed\n");
	for (size_t i = 0; i < signalCodes.size(); ++i)
		if (!dataManager.index.index_table[signalCodes[i]].is_loaded && !dataManager.sigs.Sid2Info[signalCodes[i]].virtual_sig)
			MTHROW_AND_ERR("Error in MedRegistry::create_registry - you haven't loaded %s for repository which is needed\n",
				dataManager.sigs.name(signalCodes[i]).c_str());
	for (size_t i = 0; i < physical_signals.size(); ++i)
		if (!dataManager.index.index_table[physical_ids[i]].is_loaded)
			MTHROW_AND_ERR("Error in MedRegistry::create_registry - you haven't loaded %s for repository which is needed by rep_processor!\n",
				physical_signals[i].c_str());

	used_sigs.reserve(signalCodes.size());
	if (need_bdate)
		used_sigs = signalCodes;
	else
		for (size_t i = 0; i < signalCodes.size(); ++i)
			if (dataManager.sigs.name(signalCodes[i]) != "BDATE")
				used_sigs.push_back(signalCodes[i]);

	int N_tot_threads = omp_get_max_threads();
	vector<PidDynamicRec> idRec(N_tot_threads);

	resolve_conlicts = method;
	int fixed_cnt = 0; int example_pid = -1;
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < dataManager.pids.size(); ++i)
	{
		int n_th = omp_get_thread_num();
		if (idRec[n_th].init_from_rep(std::addressof(dataManager), dataManager.pids[i], final_sigs_to_read, 1) < 0)
			MTHROW_AND_ERR("Unable to read repository\n");

		if (rep_processors != NULL && !rep_processors->empty()) {
			MedIdSamples pid_samples(dataManager.pids[i]);
			MedSample smp;
			smp.id = pid_samples.id; smp.time = INT_MAX;
			pid_samples.samples.push_back(smp);
			for (unsigned int i = 0; i < rep_processors->size(); ++i) {
				(*rep_processors)[i]->apply(idRec[n_th], pid_samples);
			}
		}

		vector<UniversalSigVec_mem> sig_vec((int)used_sigs.size());
		for (size_t k = 0; k < sig_vec.size(); ++k) {
			UniversalSigVec vv;
			idRec[n_th].uget(used_sigs[k], 0, vv);
			bool did_something = medial::repository::fix_contradictions(vv, medial::repository::fix_method::none, sig_vec[k]);
			if (did_something) {
#pragma omp atomic
				++fixed_cnt;
#pragma omp critical
				example_pid = dataManager.pids[i];
			}
		}
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
	string fixed_count_str = "";
	if (fixed_cnt > 0)
		fixed_count_str = "(fixed " + to_string(fixed_cnt) + " patient signals. example patient id=" +
		to_string(example_pid) + ")";
	MLOG("Finished creating registy in %d seconds%s\n", (int)duration, fixed_count_str.c_str());
	dataManager.clear();
}

void MedRegistry::get_registry_creation_codes(vector<string> &signal_codes) const
{
	signal_codes = signalCodes_names;
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
	const string &debug_file, const unordered_set<float> &debug_vals) const {
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
				const MedRegistryRecord &regRec = registry_records.at(indReg);
				if (regRec.registry_value < 0) {
					has_intr = true;//censored out - mark as done, no valid registry records for pid
					break;
				}

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

void MedRegistry::get_pids(vector<int> &pids) const {
	pids.clear();
	unordered_set<int> seen_pid;
	for (size_t i = 0; i < registry_records.size(); ++i)
		seen_pid.insert(registry_records[i].pid);
	pids.insert(pids.end(), seen_pid.begin(), seen_pid.end());
}

void MedRegistry::create_incidence_file(const string &file_path, const string &rep_path,
	int age_bin, int min_age, int max_age, int time_period, bool use_kaplan_meir,
	const string &sampler_name, const string &sampler_args) const {
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
	MedRepository &rep, const vector<string> &sets, float outcome_val, int chan) {
	signalName = sigName;
	buffer_duration = buffer_time;
	duration_flag = durr_time;
	take_only_first = take_first;
	outcome_value = outcome_val;
	channel = chan;
	repo = &rep;
	if (!sigName.empty()) {
		int sid = rep.sigs.sid(sigName);
		if (sid < 0)
			MTHROW_AND_ERR("Error in RegistrySignalSet::RegistrySignalSet - couldn't find signal \"%s\" in repo. maybe repo not initialized?\n",
				sigName.c_str());
		int max_chan = rep.sigs.Sid2Info[sid].n_val_channels;
		if (channel >= max_chan)
			MTHROW_AND_ERR("Error in RegistrySignalSet::RegistrySignalSet - channel %d not exists in signal \"%s\"\n",
				channel, signalName.c_str());
	}
	if (!sets.empty()) {
		int section_id = rep.dict.section_id(sigName);
		rep.dict.curr_section = section_id;
		rep.dict.default_section = section_id;
		rep.dict.prep_sets_lookup_table(section_id, sets, Flags);
	}
}

bool RegistrySignalSet::get_outcome(const UniversalSigVec &s, int current_i, float &result) {
	bool is_active = false;
	result = 0;
	is_active = !(current_i < 0 || current_i >= s.len
		|| s.Val(current_i, channel) < 0 || s.Val(current_i, channel) >= Flags.size());
	is_active = is_active && Flags[(int)s.Val(current_i, channel)];
	if (is_active)
		result = outcome_value;
	return is_active;
}

int RegistrySignalSet::init(map<string, string>& map) {

	string sets_path = "";
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
		else if (it->first == "outcome_value")
			outcome_value = stof(it->second);
		else if (it->first == "channel")
			channel = stoi(it->second);
		else if (it->first == "sets_path") //should contain "sets_path=" which points to file with list of codes
			sets_path = it->second;
		else
			MTHROW_AND_ERR("Error in RegistrySignalSet::init - unsupported element \"%s\"\n",
				it->first.c_str());
		//! [RegistrySignalSet::init]
	}
	int sid = repo->sigs.sid(signalName);
	if (sid < 0)
		MTHROW_AND_ERR("Error in RegistrySignalSet::init - couldn't find signal \"%s\" in repo. maybe repo not initialized?\n",
			signalName.c_str());
	int max_chan = repo->sigs.Sid2Info[sid].n_val_channels;
	if (channel >= max_chan)
		MTHROW_AND_ERR("Error in RegistrySignalSet::init - channel %d not exists in signal \"%s\"\n",
			channel, signalName.c_str());
	//save to end
	if (!sets_path.empty()) {
		vector<string> sets;
		medial::io::read_codes_file(sets_path, sets);
		if (!sets.empty()) {
			int section_id = repo->dict.section_id(signalName);
			repo->dict.curr_section = section_id;
			repo->dict.default_section = section_id;
			try {
				repo->dict.prep_sets_lookup_table(section_id, sets, Flags);
			}
			catch (const exception &) {
				MERR("ERROR in RegistrySignalSet::init - for signal %s(%d) sets_path=%s\n",
					signalName.c_str(), sid, sets_path.c_str());
				throw;
			}
		}
	}
	return 0;
}

RegistrySignalSet::RegistrySignalSet(const string &init_string, MedRepository &rep, const vector<string> &sets, float outcome_val) {
	repo = &rep;
	init_from_string(init_string);
	outcome_value = outcome_val;
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
	allow_prediciton_in_case = false;
	if (!skip_pid_file.empty())
		init_list(skip_pid_file, SkipPids);
	if (pid_to_censor_dates != NULL)
		pid_to_max_allowed = *pid_to_censor_dates;
	signalCodes_names.clear();
	for (size_t i = 0; i < signal_conditions.size(); ++i)
		signalCodes_names.push_back(signal_conditions[i]->signalName);
	signalCodes_names.push_back("BDATE");
	for (size_t i = 0; i < signal_conditions.size(); ++i)
		if (signal_conditions[i]->signalName == "BDATE") {
			need_bdate = true;
			break;
		}
	//the user called init for this signal_conditions
	signal_filters = signal_conditions;
}

RegistrySignalRange::RegistrySignalRange(const string &sigName, int durr_time, int buffer_time,
	bool take_first, float min_range, float max_range, float outcome_val, int chan) {
	signalName = sigName;
	duration_flag = durr_time;
	buffer_duration = buffer_time;
	take_only_first = take_first;

	min_value = min_range;
	max_value = max_range;
	outcome_value = outcome_val;
	channel = chan;
}

bool RegistrySignalRange::get_outcome(const UniversalSigVec &s, int current_i, float &result) {
	bool is_active = current_i < s.len && s.Val(current_i, channel) >= min_value && s.Val(current_i, channel) <= max_value;
	if (is_active)
		result = outcome_value;
	return is_active;
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
		else if (it->first == "outcome_value")
			outcome_value = stof(it->second);
		else if (it->first == "channel")
			channel = stoi(it->second);
		else
			MTHROW_AND_ERR("Error in RegistrySignalRange::init - unsupported element \"%s\"\n",
				it->first.c_str());
		//! [RegistrySignalRange::init]
	}
	return 0;
}

RegistrySignalDrug::RegistrySignalDrug(MedRepository &rep) {
	repo = &rep;
	signalName = "";
	duration_flag = 0;
	buffer_duration = 0;
	take_only_first = false;
	outcome_value = 1;
}

int RegistrySignalDrug::init(map<string, string>& map) {

	string sets_path = "";
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		//! [RegistrySignalDrug::init]
		if (it->first == "signalName")
			signalName = it->second;
		else if (it->first == "duration_flag")
			duration_flag = stoi(it->second);
		else if (it->first == "buffer_duration")
			buffer_duration = stoi(it->second);
		else if (it->first == "take_only_first")
			take_only_first = stoi(it->second) > 0;
		else if (it->first == "outcome_value")
			outcome_value = stof(it->second);
		else if (it->first == "sets_path") //should contain "sets=" which points to file with list of codes with TAB min_dosage_range TAB max_dosage_range
			sets_path = it->second;
		else
			MTHROW_AND_ERR("Error in RegistrySignalDrug::init - unsupported element \"%s\"\n",
				it->first.c_str());
		//! [RegistrySignalDrug::init]
	}
	//save to end
	if (!sets_path.empty()) {
		vector<string> sets;
		medial::io::read_codes_file(sets_path, sets);
		vector<int> matched_ids(sets.size());
		int max_id = 0;
		if (!sets.empty()) {
			int sid = repo->sigs.sid(signalName);
			if (sid < 0)
				MTHROW_AND_ERR("ERROR in RegistrySignalDrug::init - can't find signal %s in repository\n",
					signalName.c_str());
			int section_id = repo->dict.section_id(signalName);
			repo->dict.curr_section = section_id;
			repo->dict.default_section = section_id;

			try {
				repo->dict.prep_sets_lookup_table(section_id, sets, Flags);
			}
			catch (const exception &) {
				MERR("ERROR in RegistrySignalDrug::init - for signal %s(%d) sets_path=%s\n",
					signalName.c_str(), sid, sets_path.c_str());
				throw;
			}

			for (size_t i = 0; i < matched_ids.size(); ++i) {
				matched_ids[i] = repo->dict.id(section_id, sets[i]);
				if (matched_ids[i] > max_id)
					max_id = matched_ids[i];
			}
		}
		Flags_range.resize(max_id + 1);
		//now parse range part:
		ifstream file;
		file.open(sets_path);
		if (!file.good())
			MTHROW_AND_ERR("Error in RegistrySignalDrug::init - Unable to open test indexes file:\n%s\n", sets_path.c_str());
		string line;
		//getline(file, line); //ignore first line
		int set_id = 0;
		while (getline(file, line)) {
			boost::trim(line);
			if (line.empty())
				continue;
			if (line.at(0) == '#')
				continue;
			vector<string> tokens;
			boost::split(tokens, line, boost::is_any_of("\t"));
			if (tokens.size() != 3)
				MTHROW_AND_ERR("Error in RegistrySignalDrug::init - parsing %s file where each line should contain 3 tokens seprated by TAB. got line:\n%s\n",
					sets_path.c_str(), line.c_str());
			Flags_range[matched_ids[set_id]].first = stof(tokens[1]);
			Flags_range[matched_ids[set_id]].second = stof(tokens[2]);
			++set_id;
		}
		file.close();
	}
	return 0;
}

bool RegistrySignalDrug::get_outcome(const UniversalSigVec &s, int current_i, float &result) {
	bool is_active = false;
	result = 0;
	is_active = !(current_i < 0 || current_i >= s.len
		|| s.Val(current_i) < 0 || s.Val(current_i) >= Flags.size());
	is_active = is_active && Flags[(int)s.Val(current_i)]; //has the drug in set
	is_active = is_active && s.Val(current_i, 1) >= Flags_range[(int)s.Val(current_i)].first; //dosage check minimal
	is_active = is_active && s.Val(current_i, 1) <= Flags_range[(int)s.Val(current_i)].second; //dosage check maximal
	if (is_active)
		result = outcome_value;
	return is_active;
}

bool RegistrySignalAnd::get_outcome(const UniversalSigVec &s, int current_i, float &result) {
	bool is_active = true;
	result = 0;
	float temp;
	for (size_t i = 0; i < conditions.size() && is_active; ++i)
		is_active = conditions[i]->get_outcome(s, current_i, temp);
	if (is_active)
		result = outcome_value;
	return is_active;
}

RegistrySignalAnd::RegistrySignalAnd(MedRepository &rep) {
	repo = &rep;
	signalName = "";
	duration_flag = 0;
	buffer_duration = 0;
	take_only_first = false;
	outcome_value = 1;
	channel = 0;
}

int RegistrySignalAnd::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		//! [RegistrySignalAnd::init]
		if (it->first == "signalName")
			MWARN("Warning in RegistrySignalAnd::init - ignoring signalName argument. this is wrapper operation\n");
		else if (it->first == "conditions")  //not checking for infinite loop
			RegistrySignal::parse_registry_rules(it->second, *repo, conditions);
		else
			MTHROW_AND_ERR("ERROR in RegistrySignalAnd::init - Unsupported Argument %s\n", it->first.c_str());
		//! [RegistrySignalAnd::init]
	}
	if (conditions.empty())
		MTHROW_AND_ERR("ERROR in RegistrySignalAnd::init - conditions is empty. please use conditions to reffer to file with and conditions of signals\n");
	return 0;
}

RegistrySignalAnd::~RegistrySignalAnd() {
	for (size_t i = 0; i < conditions.size(); ++i)
		delete conditions[i];
	conditions.clear();
}

inline int Date_wrapper(const UniversalSigVec &signal, int i) {
	if (signal.get_type() != T_Value)
		return signal.Time(i);
	else
		return (int)signal.Val(i);
}

void MedRegistryCodesList::get_registry_records(int pid,
	int bdate, vector<UniversalSigVec_mem> &usv, vector<MedRegistryRecord> &results) {
	if (!init_called)
		MTHROW_AND_ERR("Must be initialized by init before use\n");
	vector<int> signals_indexes_pointers(signal_filters.size()); //all in 0
	int max_date_mark = 30000000;
	if (time_unit != MedTime::Date)
		max_date_mark = 2000000000;

	int max_allowed_date = max_repo_date;
	if (pid_to_max_allowed.find(pid) != pid_to_max_allowed.end())
		max_allowed_date = pid_to_max_allowed[pid];

	MedRegistryRecord r;
	r.pid = pid;
	int start_date = -1, last_date = -1;
	int signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
	while (signal_index >= 0)
	{
		const UniversalSigVec *signal = &usv[signal_index];
		RegistrySignal *signal_prop = signal_filters[signal_index];
		int i = signals_indexes_pointers[signal_index] - 1; //the current signal time
		//find first date if not marked already
		if (start_date == -1) {
			if (Date_wrapper(*signal, i) >= bdate) {
				start_date = Date_wrapper(*signal, i);
				r.start_date = medial::repository::DateAdd(start_date, start_buffer_duration);
			}
			else {
				signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
				continue;
			}
		}
		int min_date = medial::repository::DateAdd(start_date, start_buffer_duration);
		r.min_allowed_date = min_date; //at least 1 year data
		r.registry_value = 0;
		//I have start_date
		if (Date_wrapper(*signal, i) > max_allowed_date)
			break;
		last_date = Date_wrapper(*signal, i);
		float registry_outcome_result;
		if (signal_prop->get_outcome(*signal, i, registry_outcome_result)) {
			//flush buffer
			int last_date_c = medial::repository::DateAdd(Date_wrapper(*signal, i), -signal_prop->buffer_duration);
			r.end_date = last_date_c;
			r.max_allowed_date = last_date_c;
			if (r.end_date > r.start_date)
				results.push_back(r);

			//start new record
			//r.pid = pid;
			r.min_allowed_date = min_date;
			r.max_allowed_date = Date_wrapper(*signal, i);
			r.start_date = Date_wrapper(*signal, i);
			r.registry_value = registry_outcome_result;
			if (signal_prop->take_only_first) {
				r.end_date = max_date_mark;
				if (allow_prediciton_in_case)
					r.max_allowed_date = r.end_date;
				results.push_back(r);
				return;
			}
			else
				r.end_date = medial::repository::DateAdd(Date_wrapper(*signal, i), signal_prop->duration_flag);
			int max_search = medial::repository::DateAdd(r.end_date,
				signal_prop->buffer_duration - 1);
			//advanced till passed end_date + buffer with no reapeating RC:
			while (signal_index >= 0 && Date_wrapper(*signal, i) < max_search) {
				if (signal_prop->get_outcome(*signal, i, registry_outcome_result)) {
					r.end_date = medial::repository::DateAdd(Date_wrapper(*signal, i), signal_prop->duration_flag);
					max_search = medial::repository::DateAdd(r.end_date, signal_prop->buffer_duration - 1);
				}

				signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
				if (signal_index < 0)
					break;
				i = signals_indexes_pointers[signal_index] - 1; //the current signal time
				signal_prop = signal_filters[signal_index];
				signal = &usv[signal_index];
			}
			if (allow_prediciton_in_case)
				r.max_allowed_date = r.end_date;
			results.push_back(r);
			if (signal_index < 0) {
				r.start_date = max_date_mark; //no more control times, reached the end
				break;
			}
			//prepare for next:
			r.min_allowed_date = min_date;
			r.registry_value = 0;
			r.start_date = Date_wrapper(*signal, i); //already after duration and buffer. can start new control
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
	ofstream fo;
	if (!log_file.empty()) {
		fo.open(log_file);
		if (!fo.good())
			MWARN("Warning: can log into file %s\n", log_file.c_str());
	}
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
	string delim = ", ";
	if (histCounts.size() > 2)
		delim = "\n";
	int total = 0, total_all = 0;
	for (auto it = histCounts.begin(); it != histCounts.end(); ++it)
		total += it->second;
	for (auto it = histCounts_All.begin(); it != histCounts_All.end(); ++it)
		total_all += it->second;

	if (histCounts.size() > 2)
		log_with_file(fo, "Registry has %zu records:\n", regRecords.size());
	else if (!regRecords.empty())
		log_with_file(fo, "Registry has %zu records. [", regRecords.size());
	else {
		log_with_file(fo, "Registry is empty.\n");
		if (fo.good())
			fo.close();
		return;
	}

	auto iter = histCounts.begin();
	if (!histCounts.empty())
		log_with_file(fo, "%d=%d(%2.2f%%)", (int)iter->first, iter->second,
			100.0 * iter->second / float(total));
	++iter;
	for (; iter != histCounts.end(); ++iter)
		log_with_file(fo, "%s%d=%d(%2.2f%%)", delim.c_str(), (int)iter->first, iter->second,
			100.0 * iter->second / float(total));
	if (histCounts.size() > 2)
		log_with_file(fo, "\nAll Records:\n");
	else
		log_with_file(fo, "] All = [");

	iter = histCounts_All.begin();
	if (!histCounts_All.empty())
		log_with_file(fo, "%d=%d(%2.2f%%)", (int)iter->first, iter->second,
			100.0 * iter->second / float(total_all));
	++iter;
	for (; iter != histCounts_All.end(); ++iter)
		log_with_file(fo, "%s%d=%d(%2.2f%%)", delim.c_str(), (int)iter->first, iter->second,
			100.0 * iter->second / float(total_all));
	if (histCounts.size() > 2)
		log_with_file(fo, "\n");
	else
		log_with_file(fo, "]\n");
	if (fo.good())
		fo.close();
}

RegistrySignal *RegistrySignal::make_registry_signal(const string &type, MedRepository &rep) {
	vector<string> empty_arr;
	string empty_str = "";
	//! [RegistrySignal::make_registry_signal]
	if (type == "set")
		return new RegistrySignalSet(empty_str, 0, 0, false, rep, empty_arr);
	else if (type == "range")
		return new RegistrySignalRange(empty_str, 0, 0, false, 0, 0);
	else if (type == "drug")
		return new RegistrySignalDrug(rep);
	else if (type == "and")
		return new RegistrySignalAnd(rep);
	else
		MTHROW_AND_ERR("Error: Unsupported type \"%s\" for RegistrySignal::make_registry_signal\n", type.c_str());
	//! [RegistrySignal::make_registry_signal]
}

RegistrySignal *RegistrySignal::make_registry_signal(const string &type, MedRepository &rep, const string &init_string) {
	RegistrySignal *reg = make_registry_signal(type, rep);
	reg->init_from_string(init_string);
	return reg;
}

void RegistrySignal::parse_registry_rules(const string &reg_cfg, MedRepository &rep,
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
		boost::replace_all(ini, "$REP", rep.config_fname);
		RegistrySignal *reg = RegistrySignal::make_registry_signal(type, rep, ini);
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
		else if (it->first == "allow_prediciton_in_case")
			allow_prediciton_in_case = stoi(it->second) > 0;
		else if (it->first == "pid_to_censor_dates") {
			ifstream fr(it->second);
			if (!fr.good())
				MTHROW_AND_ERR("Error in MedRegistryCodesList::init - unable to open %s for reading.",
					it->second.c_str());
			string line;
			while (getline(fr, line)) {
				boost::trim(line);
				if (line.empty() || line.at(0) == '#')
					continue;
				vector<string> tokens;
				boost::split(tokens, line, boost::is_any_of("\t"));
				if (tokens.size() != 2)
					MTHROW_AND_ERR("Error in MedRegistryCodesList::init - in parsing pid_to_censor_dates file"
						" - %s excpeting TAB for line:\n%s\n", it->second.c_str(), line.c_str());
				if (pid_to_max_allowed.find(stoi(tokens[0])) == pid_to_max_allowed.end() ||
					pid_to_max_allowed[stoi(tokens[0])] > stoi(tokens[1]))
					pid_to_max_allowed[stoi(tokens[0])] = stoi(tokens[1]);
			}
			fr.close();
		}
		else if (it->first == "config_signals_rules")
			registry_file_path = it->second;
		else
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - Unsupported init param \"%s\"\n", it->first.c_str());
		//! [MedRegistryCodesList::init]
		if (rep_path.empty() && rep_for_init == NULL)
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - please provide rep param to init function\n");
		if (max_repo_date == 0)
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - please provide max_repo_date param to init function\n");
		if (registry_file_path.empty())
			MTHROW_AND_ERR("Error in MedRegistryCodesList::init - please provide config_signals_rules param to init function\n");

		if (rep_for_init == NULL) {
			rep_for_init = &repo;
			if (rep_for_init->init(rep_path) < 0)
				MTHROW_AND_ERR("Error in MedRegistryCodesList::init - Unable to init repositrory from path %s\n", rep_path.c_str());
		}

		RegistrySignal::parse_registry_rules(registry_file_path, *rep_for_init, signal_filters);
		init_called = true;

		signalCodes_names.clear();
		for (size_t i = 0; i < signal_filters.size(); ++i)
			signalCodes_names.push_back(signal_filters[i]->signalName);
		signalCodes_names.push_back("BDATE");
		for (size_t i = 0; i < signal_filters.size(); ++i)
			if (signal_filters[i]->signalName == "BDATE") {
				need_bdate = true;
				break;
			}

		return 0;
}

MedRegistry *MedRegistry::make_registry(const string &registry_type, const string &init_str) {
	MedRegistry *registry;
	//! [MedRegistry::make_registry]
	if (registry_type == "binary")
		registry = new MedRegistryCodesList;
	else if (registry_type == "categories")
		registry = new MedRegistryCategories;
	else
		MTHROW_AND_ERR("Error: Unsupported type \"%s\" for MedRegistry::make_registry\n",
			registry_type.c_str());
	//! [MedRegistry::make_registry]
	if (!init_str.empty())
		registry->init_from_string(init_str);

	return registry;
}

MedRegistry *MedRegistry::make_registry(const string &registry_type, MedRepository &rep, const string &init_str) {
	MedRegistry *registry = make_registry(registry_type, "");
	registry->rep_for_init = &rep;

	if (!init_str.empty())
		registry->init_from_string(init_str);

	return registry;
}

int MedRegistryCategories::init(map<string, string>& map) {
	string repository_path, registry_cfg_path;
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		//! [MedRegistryCategories::init]
		if (it->first == "rep")
			repository_path = it->second;
		else if (it->first == "start_buffer_duration")
			start_buffer_duration = stoi(it->second);
		else if (it->first == "end_buffer_duration")
			end_buffer_duration = stoi(it->second);
		else if (it->first == "max_repo_date")
			max_repo_date = stoi(it->second);
		else if (it->first == "config_signals_rules")
			registry_cfg_path = it->second;
		else if (it->first == "pid_to_censor_dates") {
			ifstream fr(it->second);
			if (!fr.good())
				MTHROW_AND_ERR("Error in MedRegistryCategories::init - unable to open %s for reading.",
					it->second.c_str());
			string line;
			while (getline(fr, line)) {
				boost::trim(line);
				if (line.empty() || line.at(0) == '#')
					continue;
				vector<string> tokens;
				boost::split(tokens, line, boost::is_any_of("\t"));
				if (tokens.size() != 2)
					MTHROW_AND_ERR("Error in MedRegistryCategories::init - in parsing pid_to_censor_dates file"
						" - %s excpeting TAB for line:\n%s\n", it->second.c_str(), line.c_str());
				if (pid_to_max_allowed.find(stoi(tokens[0])) == pid_to_max_allowed.end() ||
					pid_to_max_allowed[stoi(tokens[0])] > stoi(tokens[1]))
					pid_to_max_allowed[stoi(tokens[0])] = stoi(tokens[1]);
			}
			fr.close();
		}
		else
			MTHROW_AND_ERR("Error in MedRegistryCategories::init - Unsupported init param \"%s\"\n",
				it->first.c_str());
		//! [MedRegistryCategories::init]
	}

	if (repository_path.empty() && rep_for_init == NULL)
		MTHROW_AND_ERR("Error in MedRegistryCategories::init - please provide repository param to init function\n");
	if (registry_cfg_path.empty())
		MTHROW_AND_ERR("Error in MedRegistryCategories::init - please provide config_signals_rules param to init function\n");

	MedPidRepository repo;
	if (rep_for_init == NULL) {
		rep_for_init = &repo;
		if (rep_for_init->init(repository_path) < 0)
			MTHROW_AND_ERR("Error in MedRegistryCategories::init - Unable to init repositrory from path %s\n", repository_path.c_str());
	}

	vector<RegistrySignal *> all_rules;
	RegistrySignal::parse_registry_rules(registry_cfg_path, *rep_for_init, all_rules);
	//transpose to all_rules -> signals_rules:
	unordered_map<string, int> signal_name_to_idx;
	for (size_t i = 0; i < all_rules.size(); ++i)
	{
		int current_signal_idx = (int)signals_rules.size();
		if (signal_name_to_idx.find(all_rules[i]->signalName) == signal_name_to_idx.end()) {
			signal_name_to_idx[all_rules[i]->signalName] = current_signal_idx;
			signals_rules.resize(current_signal_idx + 1); //open new empty signal rules list
		}
		else
			current_signal_idx = signal_name_to_idx[all_rules[i]->signalName];
		signals_rules[current_signal_idx].push_back(all_rules[i]);
	}

	signalCodes_names.clear();
	for (size_t i = 0; i < signals_rules.size(); ++i)
		signalCodes_names.push_back(signals_rules[i][0]->signalName);
	signalCodes_names.push_back("BDATE");
	for (size_t i = 0; i < signals_rules.size(); ++i)
		if (signals_rules[i][0]->signalName == "BDATE") {
			need_bdate = true;
			break;
		}

	return 0;
}

void MedRegistryCategories::get_registry_records(int pid, int bdate, vector<UniversalSigVec_mem> &usv,
	vector<MedRegistryRecord> &results) {
	if (signals_rules.empty())
		MTHROW_AND_ERR("Must be initialized by init before use\n");
	vector<int> signals_indexes_pointers(signals_rules.size()); //all in 0
	int max_allowed_date = max_repo_date;
	if (pid_to_max_allowed.find(pid) != pid_to_max_allowed.end())
		if (pid_to_max_allowed[pid] > 0)
			max_allowed_date = pid_to_max_allowed[pid];
		else
			return;

	unordered_set<float> outcomes_may_not_use;
	int last_buffer_duration = -1;
	int max_date_mark = 30000000;
	if (time_unit != MedTime::Date)
		max_date_mark = 2000000000;

	MedRegistryRecord r;
	r.pid = pid;
	r.registry_value = -1;
	int start_date = -1, last_date = -1;
	int signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
	//fetch till passing bdate
	while (signal_index >= 0 && start_date == -1) {
		UniversalSigVec &signal = usv[signal_index];
		int i = signals_indexes_pointers[signal_index] - 1; //the current signal time
		//find first date if not marked already
		if (Date_wrapper(signal, i) >= bdate)
			start_date = Date_wrapper(signal, i);
		else
			signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
	}
	int min_date = medial::repository::DateAdd(start_date, start_buffer_duration);
	r.min_allowed_date = min_date;
	r.start_date = min_date;

	bool same_date = false;
	float rule_activated = 0;
	bool is_rule_active = false;
	bool mark_no_match = true;
	while (signal_index >= 0)
	{
		const UniversalSigVec &signal = usv[signal_index];
		vector<RegistrySignal *> *all_signal_prop = &signals_rules[signal_index];
		int i = signals_indexes_pointers[signal_index] - 1; //the current signal time

		//I have start_date
		int curr_date = Date_wrapper(signal, i);
		if (max_allowed_date > 0 && curr_date > max_allowed_date)
			break;


		if (!same_date) {
			rule_activated = 0;
			is_rule_active = false;
		}
		for (size_t rule_idx = 0; rule_idx < all_signal_prop->size(); ++rule_idx)
		{
			RegistrySignal *signal_prop = (*all_signal_prop)[rule_idx];

			float signal_prop_outcome;
			if (signal_prop->get_outcome(signal, i, signal_prop_outcome)) {
				if (is_rule_active && rule_activated != signal_prop_outcome) {//validates no contradicted rule passes this condition
					if (resolve_conlicts == medial::repository::fix_method::none) {
						MTHROW_AND_ERR("Error in MedRegistryCategories - specific signal \"%s\" has contradicted"
							" rules in same time point with diffrent outcomes(pid=%d, time=%d, value=%2.3f)\n",
							signal_prop->signalName.c_str(), pid, curr_date, signal.Val(i));
					}
					else if (resolve_conlicts == medial::repository::fix_method::take_first) {
						break;
					}
					else if (resolve_conlicts == medial::repository::fix_method::take_max) {
						if (rule_activated > signal_prop_outcome)
							continue;
					}
					else if (resolve_conlicts == medial::repository::fix_method::take_min) {
						if (rule_activated < signal_prop_outcome)
							continue;
					}
					else if (resolve_conlicts == medial::repository::fix_method::drop) {
						mark_no_match = true;
						break;
					}
					else if (resolve_conlicts == medial::repository::fix_method::take_last) {
						//do nothing - continue
					}
					else
						MTHROW_AND_ERR("Resolve Conflict mode %d isn't supported\n", resolve_conlicts);
				}
				rule_activated = signal_prop_outcome;
				is_rule_active = true;

				//check if we need to merge this outcome with previous one current state or open new one:
				if (r.registry_value == signal_prop_outcome && !mark_no_match) {
					//same outcome - update end_time:
					if (curr_date < medial::repository::DateAdd(r.end_date,
						last_buffer_duration - 1)) {

						int new_end_date = medial::repository::DateAdd(curr_date, signal_prop->duration_flag);
						if (new_end_date > r.end_date) {
							int prev_search = medial::repository::DateAdd(r.end_date, last_buffer_duration - 1);
							r.end_date = new_end_date;
							int new_search_date = medial::repository::DateAdd(r.end_date, signal_prop->buffer_duration - 1);
							last_buffer_duration = signal_prop->buffer_duration;
							if (new_search_date < prev_search)
								last_buffer_duration += med_time_converter.convert_times(global_default_time_unit, global_default_windows_time_unit, prev_search) -
								med_time_converter.convert_times(global_default_time_unit, global_default_windows_time_unit, new_search_date);
						}
					}
					else {
						if (signal_prop->take_only_first && outcomes_may_not_use.find(signal_prop_outcome) != outcomes_may_not_use.end())
							continue;
						//finished time - flush and open new registry with 0 outcome:
						if (signal_prop->take_only_first) {
							r.end_date = max_date_mark;
							//results.push_back(r);
							outcomes_may_not_use.insert(signal_prop_outcome); //if happens again will ignore and skip
							last_buffer_duration = -1; //no buffer duration
							results.push_back(r);
							mark_no_match = true;
							continue;
						}

						if (r.end_date > r.start_date && r.max_allowed_date > r.min_allowed_date)
							results.push_back(r);
						//start new record with 0 outcome:
						//r.registry_value = signal_prop_outcome; //left the same no need
						r.start_date = curr_date; //continue from where left
						//r.end_date = medial::repository::DateAdd(r.start_date, 1); //let's start from 1 day
						r.max_allowed_date = curr_date;
						r.end_date = medial::repository::DateAdd(curr_date, signal_prop->duration_flag);
						last_buffer_duration = signal_prop->buffer_duration;
					}
				}
				else { //diffrent outcome - no contradiction in same time point:
					//flush last 
					int last_date_c = medial::repository::DateAdd(curr_date, -signal_prop->buffer_duration);
					if (!mark_no_match && last_date_c < r.end_date)
						r.end_date = last_date_c;
					//if (r.registry_value == 0)
					//	r.max_allowed_date = last_date_c;

					if (!mark_no_match && r.end_date > r.start_date && r.max_allowed_date > r.min_allowed_date)
						results.push_back(r);

					//skip if may not use
					if (signal_prop->take_only_first && outcomes_may_not_use.find(signal_prop_outcome) != outcomes_may_not_use.end())
						continue;
					//start new record
					//r.pid = pid;
					//r.min_allowed_date = min_date;
					r.max_allowed_date = curr_date;
					r.start_date = curr_date;
					r.registry_value = signal_prop_outcome;

					if (signal_prop->take_only_first) {
						r.end_date = max_date_mark;
						//results.push_back(r);
						outcomes_may_not_use.insert(signal_prop_outcome); //if happens again will ignore and skip
						last_buffer_duration = -1; //no buffer duration
						results.push_back(r);
						mark_no_match = true;
						continue; //finished handling!
					}

					r.end_date = medial::repository::DateAdd(curr_date, signal_prop->duration_flag);
					last_buffer_duration = signal_prop->buffer_duration;

					//results.push_back(r);
				}
				mark_no_match = false;
			}
		}

		if (!same_date && !is_rule_active) {
			//check if need to close buffer - no rule happend in this time and has outcome in buffer
			if (curr_date > r.end_date) { //if need to skip or has passed buffer
				if (!mark_no_match && r.end_date > r.start_date && r.max_allowed_date > r.min_allowed_date)
					results.push_back(r);
				//start new record with 0 outcome:
				r.registry_value = -1;
				//r.start_date = r.end_date; //continue from where left
				//r.age = (int)medial::repository::DateDiff(bdate, r.start_date);
				mark_no_match = true;
				//r.end_date = medial::repository::DateAdd(r.start_date, 1); //let's start from 1 day
				//r.max_allowed_date = r.end_date;
			}
		}

		last_date = curr_date;
		signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
		if (signal_index >= 0)
			same_date = last_date == Date_wrapper(usv[signal_index], signals_indexes_pointers[signal_index] - 1);
	}

	if (mark_no_match)
		r.end_date = last_date;
	last_date = medial::repository::DateAdd(last_date, -end_buffer_duration);
	r.max_allowed_date = last_date;
	if (r.end_date > r.start_date && r.max_allowed_date > r.min_allowed_date && !mark_no_match)
		results.push_back(r);
}

void MedRegistryCategories::clear_create_variables() {
	for (size_t i = 0; i < signals_rules.size(); ++i) {
		for (size_t j = 0; j < signals_rules[i].size(); ++j)
			delete signals_rules[i][j];
		signals_rules[i].clear();
	}
	signals_rules.clear();
}

void MedRegistryCodesList::clear_create_variables() {
	for (size_t i = 0; i < signal_filters.size(); ++i)
		delete signal_filters[i];
	signal_filters.clear();
}