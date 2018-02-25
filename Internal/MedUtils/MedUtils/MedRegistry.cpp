#include "MedRegistry.h"
#include <fstream>
#include <string>
#include <iostream>
#include <random>
#include "Logger/Logger/Logger.h"
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

float DateDiff(int refDate, int dateSample) {
	return float((med_time_converter.convert_date(MedTime::Days, dateSample) -
		med_time_converter.convert_date(MedTime::Days, refDate)) / 365.0);
}

int DateAdd(int refDate, int daysAdd) {
	return med_time_converter.convert_days(MedTime::Date,
		med_time_converter.convert_date(MedTime::Days, refDate) + daysAdd);
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
}

int get_value(MedRepository &rep, int pid, int sigCode) {
	int len, gend = -1;
	void *data = rep.get(pid, sigCode, len);
	if (len > 0)
		gend = (int)(*(SVal *)data).val;
	return gend;
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
		int birth = get_value(dataManager, dataManager.pids[i], bDateCode);
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

#pragma region Readcodes hirerachy

bool is_code_chr(char ch) {
	return (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') || (ch >= '0' && ch <= '9') || ch == '.';
}

bool all_code_chars(const string &str) {
	int i = 0;
	while (i < str.size() && is_code_chr(str.at(i)))
		++i;
	return i >= str.size();
}

string filter_g_code(const vector<string> &vec) {
	for (string s : vec)
	{
		if (s.size() == 7 && s.at(0) == 'G' && s.at(1) == '_') {
			return s;
		}
		if (s.size() == 7 && all_code_chars(s)) {
			//return "G_" + s.substr(0, 5);
			return s;
		}
	}
	return "";
}

vector<int> get_parents_rc(MedDictionarySections &dict, const string &group) {
	int sectionId = 1;
	vector<int> res(5);
	static vector<string> constStrings = { "",".", "..", "...", "...." };

	int ii = 0;
	string lookupVal = group;
	if (all_code_chars(group)) {
		int code = dict.dicts[sectionId].id(group);
		res.insert(res.begin(), code);
		++ii;
		lookupVal = "G_" + group.substr(0, 5);
	}

	for (size_t i = 0; i < 5; ++i)
	{
		if (lookupVal.at(lookupVal.size() - 1 - i) == '.') {
			res.pop_back();
			continue; //skip - will do next
		}
		string manip = lookupVal.substr(0, lookupVal.size() - i) + constStrings[i];
		int code = dict.id(sectionId, manip);

		res[ii] = code;
		++ii;
	}

	return res;
}

#pragma endregion

#pragma region ATC hirerachy

bool is_numeric(const string &str) {
	int i = 2;
	while (i < str.size() && (str[i] >= '0' && str[i] <= '9')) {
		++i;
	}
	return i >= str.size();
}

string filter_atc_code(const vector<string> &vec) {
	for (string s : vec)
	{
		if (s.size() >= 12 && s.at(0) == 'A' && s.at(1) == 'T' && s.at(2) == 'C' && s.at(3) == '_') {
			if (s.size() == 12)
				return s;
			else
				return s.substr(0, 12);
		}
		if (s.size() == 10 && s.at(0) == 'd' && s.at(1) == 'c' && is_numeric(s)) {
			return s;
		}
	}
	return "";
}

vector<int> get_parents_atc(MedDictionarySections &dict, const string &group) {
	int sectionId = dict.section_id("Drug");
	int maxH = 5;
	vector<int> res(maxH);
	static vector<string> constStrings = { "", "__", "___", /* skip _ in last*/ "_____", "_______" };
	static vector<int> indexLookup = { 1, 3, 5, 7 };

	string lookupCode = group;
	int ii = 0;
	//If dc - find it's ATC (ATc si size 12, dc is 10)
	if (group.size() == 10) {
		int code = dict.id(sectionId, group);
		res.insert(res.begin(), code);
		++ii;
		++maxH;

		//search all options for ATC:
		vector<int> sets;
		dict.get_member_sets(sectionId, code, sets);
		bool found_atc = false;
		for (int set_i : sets)
		{
			if (dict.dicts[sectionId].Id2Names.find(set_i) == dict.dicts[sectionId].Id2Names.end()) {
				continue;
			}
			vector<string> set_names = dict.dicts[sectionId].Id2Names.at(set_i);
			for (string set_name : set_names) {
				if (set_name.size() == 12 && set_name.at(0) == 'A' && set_name.at(1) == 'T' && set_name.at(2) == 'C' && set_name.at(3) == '_') {
					found_atc = true;
					lookupCode = set_name;
					break;
				}
			}
			if (found_atc)
			{
				break;
			}
		}
	}

	for (size_t i = 0; i < indexLookup.size(); ++i)
	{
		if (lookupCode.at(lookupCode.size() - indexLookup[i]) == '_') {
			res.pop_back();
			continue; //skip - will do next anyway
		}
		string manip = lookupCode.substr(0, lookupCode.size() - constStrings[i].size()) + constStrings[i];
		int code = dict.id(sectionId, manip);

		res[ii] = code;
		++ii;
	}
	string manip = lookupCode.substr(0, lookupCode.size() - constStrings[constStrings.size() - 1].size()) + constStrings[constStrings.size() - 1];
	int code = dict.id(sectionId, manip);

	res[ii] = code;

	return res;
}

#pragma endregion

#pragma region BNF hirerachy

string filter_bnf_code(const vector<string> &vec) {
	for (string s : vec)
	{
		if (s.size() == 15 && s.at(0) == 'B' && s.at(1) == 'N' && s.at(2) == 'F' && s.at(4) == '_') {
			return s;
		}
		if (s.size() == 10 && s.at(0) == 'd' && s.at(1) == 'c' && is_numeric(s)) {
			return s;
		}
	}
	return "";
}

vector<int> get_parents_bnf(MedDictionarySections &dict, const string &group) {
	int sectionId = dict.section_id("Drug");
	int maxH = 4;
	vector<int> res(maxH);
	static vector<int> indexLookup = { 1, 4, 7 };
	static vector<string> constStrings = { "", "00", "00.00", "00.00.00" };

	string lookupCode = group;
	int ii = 0;
	//If dc - find it's BNF (BNF size 15, dc is 10)
	if (group.size() == 10) {
		int code = dict.id(sectionId, group);
		res.insert(res.begin(), code);
		++ii;
		++maxH;

		//search all options for BNF:
		vector<int> sets;
		dict.get_member_sets(sectionId, code, sets);
		bool found_bnf = false;
		for (int set_i : sets)
		{
			if (dict.dicts[sectionId].Id2Names.find(set_i) == dict.dicts[sectionId].Id2Names.end()) {
				continue;
			}
			vector<string> set_names = dict.dicts[sectionId].Id2Names.at(set_i);
			for (string set_name : set_names) {
				if (set_name.size() == 15 && set_name.at(0) == 'B' && set_name.at(1) == 'N' && set_name.at(2) == 'F' && set_name.at(3) == '_') {
					found_bnf = true;
					lookupCode = set_name;
					break;
				}
			}
			if (found_bnf)
			{
				break;
			}
		}
	}

	for (size_t i = 0; i < indexLookup.size(); ++i)
	{
		if (lookupCode.at(lookupCode.size() - indexLookup[i]) == '0' && lookupCode.at(lookupCode.size() - indexLookup[i] - 1) == '0') {
			res.pop_back();
			continue; //skip - will do next
		}
		string manip = lookupCode.substr(0, lookupCode.size() - constStrings[i].size()) + constStrings[i];
		int code = dict.id(sectionId, manip);

		res[ii] = code;
		++ii;
	}

	string manip = lookupCode.substr(0, lookupCode.size() - constStrings[constStrings.size() - 1].size()) + constStrings[constStrings.size() - 1];
	int code = dict.id(sectionId, manip);
	res[ii] = code;

	return res;
}

#pragma endregion

string get_readcode_code(MedDictionarySections &dict, int id, string(*filterCode)(const vector<string> &), int sectionId) {
	string res = "";
	vector<int> sets;
	dict.get_member_sets(sectionId, id, sets);
	if (dict.dicts[sectionId].Id2Names.find(id) != dict.dicts[sectionId].Id2Names.end()) {
		vector<string> nms = dict.dicts[sectionId].Id2Names.at(id);
		string filterd = filterCode(nms);
		if (!filterd.empty()) {
			res = filterd;
			return res;
		}
	}
	for (int s : sets)
	{
		vector<string> names;
		if (dict.dicts[sectionId].Id2Names.find(s) != dict.dicts[sectionId].Id2Names.end()) {
			names = dict.dicts[sectionId].Id2Names.at(s);
		}
		string name = filterCode(names);
		if (!name.empty() && !res.empty()) {
			//more than one option - for now not supported:
			cerr << "not supported" << endl;
			throw logic_error("not supported");
		}
		if (!name.empty()) {
			res = name;
		}
	}

	return res;
}

void MedRegistry::getRecords_Hir(int pid, vector<UniversalSigVec> &signals, MedDictionarySections &dict,
	const string &signalHirerchyType,
	vector<MedRegistryRecord> &res) {
	UniversalSigVec &signalVal = signals[0];

	for (int i = 0; i < signalVal.len; ++i)
	{
		MedRegistryRecord rec;
		rec.pid = pid;
		rec.age = -1;
		rec.start_date = signalVal.Date(i);
		rec.end_date = DateAdd(rec.start_date, 1);
		rec.registry_value = signalVal.Val(i);
		res.push_back(rec);
		if (signalVal.Val(i) <= 0)
			continue; //has no hirerachy
		//take care of hirerachy:
		auto func_filter = filter_g_code;
		auto func_parents = get_parents_rc;
		int sectionId = 1;
		if (signalHirerchyType == "None") {
			//do nothing - no hirarchy
		}
		else if (signalHirerchyType == "RC") {
			func_filter = filter_g_code;
			func_parents = get_parents_rc;
			sectionId = 1;
		}
		else if (signalHirerchyType == "ATC") {
			func_filter = filter_atc_code;
			func_parents = get_parents_atc;
			sectionId = 2;
		}
		else if (signalHirerchyType == "BNF")
		{
			func_filter = filter_bnf_code;
			func_parents = get_parents_bnf;
			sectionId = 2;
		}
		else
			MTHROW_AND_ERR("Signal_Hirerchy_Type not suppoted %s. options are: ATC,BNF,RC,None\n",
				signalHirerchyType.c_str());


		string s = get_readcode_code(dict, (int)signalVal.Val(i), func_filter, sectionId);
		if (s.empty())
			continue;

		vector<int> nums = func_parents(dict, s);
		for (size_t k = 1; k < nums.size() && k <= 3; ++k) //take till 3
		{
			if (nums[k] <= 0)
				continue;
			MedRegistryRecord rec2;

			rec2.pid = pid;
			rec2.age = -1;
			rec2.start_date = signalVal.Date(i);
			rec2.end_date = DateAdd(rec2.start_date, 1);
			rec2.registry_value = (float)nums[k];
			res.push_back(rec2);
		}
	}
}

/// Check if signal date intersect with registry record in time window: min_dur, duration
bool date_intersection(int min_allowed_date, int max_allowed_date, int reg_start, int reg_end, int sig_date,
	int time_window_from, int time_window_to) {
	//Registry, Signal
	int sig_start_date = DateAdd(sig_date, time_window_from);
	int sig_end_date = DateAdd(sig_date, time_window_to);
	if (sig_start_date > max_allowed_date || sig_start_date < min_allowed_date)
		return false;

	return (sig_end_date >= reg_start) && (sig_start_date <= reg_end);
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
			int gend = get_value(dataManager, rec.id, genderCode);
			int bdate = get_value(dataManager, rec.id, bdateCode);
			if (gend == -1) {
				++unknown_gender;
				continue;
			}
			double curr_age = DateDiff(bdate, rec.time);

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
		int gender = get_value(dataManager, pid, genderCode);
		int BDate = get_value(dataManager, pid, bdateCode);
		vector<UniversalSigVec> patientFile(1);
		dataManager.uget(pid, signalCode, patientFile[0]);

		vector<MedRegistryRecord> signal_vals;
		getRecords_Hir(pid, patientFile, dataManager.dict, signalHirerchyType, signal_vals);

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
				if (regRec.age != -1)
					ageBin = float(ageBinValue * floor(double(regRec.age) / ageBinValue));
				else
					ageBin = float(ageBinValue * floor(double(DateDiff(BDate, sigRec.start_date)) / ageBinValue));
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
					time_window_from, time_window_to);
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

	//update values prevalence
	for (auto it = maleSignalToStats.begin(); it != maleSignalToStats.end(); ++it)
	{
		for (auto jt = maleSignalToStats[it->first].begin(); jt != maleSignalToStats[it->first].end(); ++jt) {
			maleSignalToStats[it->first][jt->first][0] = male_total_prevalence[jt->first][0] - maleSignalToStats[it->first][jt->first][2];
			maleSignalToStats[it->first][jt->first][1] = male_total_prevalence[jt->first][1] - maleSignalToStats[it->first][jt->first][3];
		}
		for (auto jt = femaleSignalToStats[it->first].begin(); jt != femaleSignalToStats[it->first].end(); ++jt) {
			femaleSignalToStats[it->first][jt->first][0] = female_total_prevalence[jt->first][0] - femaleSignalToStats[it->first][jt->first][2];
			femaleSignalToStats[it->first][jt->first][1] = female_total_prevalence[jt->first][1] - femaleSignalToStats[it->first][jt->first][3];
		}
	}

	duration = (int)difftime(time(NULL), start);
	MLOG("Finished in %d seconds with %d records in males and %d records in females\n",
		duration, (int)maleSignalToStats.size(), (int)femaleSignalToStats.size());

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

void MedRegistryCodesList::init_lists(MedRepository &rep, int dur_flag, int buffer_dur, bool takeOnlyFirst,
	int max_repo, const vector<string> *rc_sets, const string &skip_pid_file) {
	if (rc_sets != NULL) {
		int section_id = rep.dict.section_id("RC");
		rep.dict.curr_section = section_id;
		rep.dict.default_section = section_id;
		rep.dict.prep_sets_lookup_table(section_id, *rc_sets, RCFlags);
	}

	duration_flag = dur_flag;
	buffer_duration = buffer_dur;
	take_only_first = takeOnlyFirst;
	max_repo_date = max_repo;

	if (!skip_pid_file.empty())
		init_list(skip_pid_file, SkipPids);
	signalCodes.clear();
	signalCodes.push_back(rep.sigs.sid("RC"));
}

void MedRegistryCodesList::get_registry_records(int pid,
	int bdate, vector<UniversalSigVec> &usv, vector<MedRegistryRecord> &results) {

	UniversalSigVec &signal = usv[0]; //RC signal
	if (signal.len <= 0)
		return;
	int start_date = signal.Date(0);
	for (int i = 1; i < signal.len && start_date < bdate; ++i)
		start_date = signal.Date(i); //finds first legal date

	int min_date = DateAdd(start_date, 365);

	MedRegistryRecord r;
	r.pid = pid;
	r.min_allowed_date = min_date; //at least 1 year data
	r.start_date = min_date;
	r.age = int(DateDiff(bdate, r.start_date));
	r.registry_value = 0;

	int last_date = start_date;
	for (int i = 0; i < signal.len; ++i)
	{
		if (signal.Date(i) > max_repo_date)
			break;
		last_date = signal.Date(i);
		if (signal.Val(i) > 0 && RCFlags[(int)signal.Val(i)]) {
			//flush buffer
			int last_date = DateAdd(signal.Date(i), -buffer_duration);
			r.end_date = last_date;
			r.max_allowed_date = last_date;
			if (r.end_date > r.start_date)
				results.push_back(r);

			//start new record
			//r.pid = pid;
			//r.male = gender == 1;
			r.min_allowed_date = last_date;
			r.max_allowed_date = signal.Date(i);
			r.start_date = signal.Date(i);
			r.age = (int)round(DateDiff(bdate, signal.Date(i)));
			r.registry_value = 1;
			if (take_only_first) {
				r.end_date = 30000000;
				results.push_back(r);
				return;
			}
			else
				r.end_date = DateAdd(signal.Date(i), duration_flag);
			int max_search = DateAdd(r.end_date, buffer_duration - 1);
			//advanced till passed end_date + buffer with no reapeating RC:
			while (i < signal.len && signal.Date(i) < max_search) {
				if (signal.Val(i) > 0 && RCFlags[(int)signal.Val(i)])
					r.end_date = DateAdd(signal.Date(i), duration_flag);
				++i;
			}
			results.push_back(r);
			if (i >= signal.len) {
				r.start_date = 30000000; //no more control times, reached the end
				break;
			}
			//prepare for next:
			start_date = DateAdd(signal.Date(i), buffer_duration); //next time don't predict before last record
			if (start_date < min_date)
				start_date = min_date;
			r.min_allowed_date = start_date;
			r.registry_value = 0;
			r.start_date = start_date;
			r.age = int(DateDiff(bdate, r.start_date));
			--i; //return to previous one...
		}
	}

	r.end_date = last_date;
	last_date = DateAdd(last_date, -buffer_duration);
	r.max_allowed_date = last_date;
	if (r.end_date > r.start_date)
		results.push_back(r);
}