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
#include <boost/filesystem.hpp>
#include <omp.h>
#include <MedMat/MedMat/MedMat.h>

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
	tokens.reserve(4);
	//Format: [ PID, Start_Date, End_Date, RegistryValue ]
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

		if (tokens.size() != 4) {
			cerr << "Bad File Format in line " << lineNum << " got \"" << line << "\"" << " parsed " << tokens.size() << " tokens" << endl;
			throw out_of_range("File has bad format");
		}

		MedRegistryRecord pr;
		pr.pid = stoi(tokens[0]);
		pr.start_date = read_time(time_unit, tokens[1]);
		pr.end_date = read_time(time_unit, tokens[2]);
		pr.registry_value = stof(tokens[3]);
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
	fw << "# Format: PID, Start_Date, End_Date, RegistryValue" << endl;
	fw << "# Created By Script - Insert Comments following #..." << endl;
	fw << "TIME_UNIT" << delim << med_time_converter.type_to_string(time_unit) << endl;
	fw << endl;
	for (size_t i = 0; i < registry_records.size(); ++i)
		fw << registry_records[i].pid << delim <<
		write_time(time_unit, registry_records[i].start_date) << delim <<
		write_time(time_unit, registry_records[i].end_date) << delim
		<< registry_records[i].registry_value << "\n";

	fw.flush();
	fw.close();
	MLOG("[Wrote %d registry records to %s]\n", (int)registry_records.size(), file_path.c_str());
}

void MedRegistry::create_registry(MedPidRepository &dataManager, medial::repository::fix_method method, vector<RepProcessor *> *rep_processors) {
	MLOG_D("Creating registry...\n");
	vector<int> used_sigs;

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
	MedProgress prog_pid("MedRegistry::create_registry", (int)dataManager.pids.size(), 15);
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
			idRec[n_th].uget(used_sigs[k], idRec[n_th].get_n_versions() > 0 ? idRec[n_th].get_n_versions() - 1 : 0, vv); //get last version
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
		registry_records.insert(registry_records.end(), vals.begin(), vals.end());
		prog_pid.update();
	}

	string fixed_count_str = "";
	if (fixed_cnt > 0)
		fixed_count_str = " (fixed " + to_string(fixed_cnt) + " patient signals. example patient id=" +
		to_string(example_pid) + ")";
	MLOG("Finished creating registy%s\n", fixed_count_str.c_str());
	dataManager.clear();
}

void MedRegistry::get_registry_creation_codes(vector<string> &signal_codes) const
{
	signal_codes = signalCodes_names;
}

void MedRegistry::get_registry_use_codes(vector<string> &signal_codes) const
{
	if (need_bdate)
		signal_codes = signalCodes_names;
	else {
		signal_codes.clear();
		for (size_t i = 0; i < signalCodes_names.size(); ++i)
			if (signalCodes_names[i] != "BDATE")
				signal_codes.push_back(signalCodes_names[i]);
	}
}


void MedRegistry::get_pids(vector<int> &pids) const {
	pids.clear();
	unordered_set<int> seen_pid;
	for (size_t i = 0; i < registry_records.size(); ++i)
		seen_pid.insert(registry_records[i].pid);
	pids.insert(pids.end(), seen_pid.begin(), seen_pid.end());
}

void *MedRegistry::new_polymorphic(string dname) {
	//SERIALIZATION NOT SUPPORTED YET:
	//CONDITIONAL_NEW_CLASS(dname, MedRegistryCodesList);
	//CONDITIONAL_NEW_CLASS(dname, MedRegistryCategories);

	CONDITIONAL_NEW_CLASS(dname, MedRegistryKeepAlive);

	MTHROW_AND_ERR("Error in MedRegistry::new_polymorphic - Unsupported class %s\n", dname.c_str());
}

void *RegistrySignal::new_polymorphic(string dname) {
	CONDITIONAL_NEW_CLASS(dname, RegistrySignalAny);

	MTHROW_AND_ERR("Error in RegistrySignal::new_polymorphic - Unsupported class %s\n", dname.c_str());
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
	MedRepository &rep, const vector<string> &sets, const string &path_to_cfg_file,
	float outcome_val, int chan) {
	signalName = sigName;
	buffer_duration = buffer_time;
	duration_flag = durr_time;
	take_only_first = take_first;
	outcome_value = outcome_val;
	channel = chan;
	repo = &rep;
	boost::filesystem::path p(path_to_cfg_file);
	base_cfg_path = p.parent_path().string();
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

int RegistrySignal::init(map<string, string>& mapper) {
	map<string, string> rest_arguments;
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RegistrySignal::init]
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
		else
			rest_arguments[it->first] = it->second;
		//! [RegistrySignal::init]
	}

	_init(rest_arguments);
	return 0;
}

void RegistrySignalSet::_init(const map<string, string>& mapper) {

	string sets_path = "";
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RegistrySignalSet::_init]
		if (it->first == "sets_path") //should contain "sets_path=" which points to file with list of codes
			sets_path = it->second;
		else
			MTHROW_AND_ERR("Error in RegistrySignalSet::init - unsupported element \"%s\"\n",
				it->first.c_str());
		//! [RegistrySignalSet::_init]
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
	if (boost::starts_with(sets_path, "path_rel:") && !base_cfg_path.empty()) {
		//use relative path to the current config file
		sets_path = sets_path.substr(9);
		//use base_paht as fix to relative path:
		if (!((sets_path.size() > 2 && (sets_path[0] == '/' || sets_path[1] == ':')))) 
			sets_path = base_cfg_path + path_sep() + sets_path;
	}
	if (!sets_path.empty()) {
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
}

RegistrySignalSet::RegistrySignalSet(const string &init_string, MedRepository &rep,
	const vector<string> &sets, const string &path_to_cfg_file, float outcome_val) {
	repo = &rep;
	init_from_string(init_string);
	outcome_value = outcome_val;
	boost::filesystem::path p(path_to_cfg_file);
	base_cfg_path = p.parent_path().string();
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

void RegistrySignalRange::_init(const map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RegistrySignalRange::_init]
		if (it->first == "min_value")
			min_value = stof(it->second);
		else if (it->first == "max_value")
			max_value = stof(it->second);
		else
			MTHROW_AND_ERR("Error in RegistrySignalRange::init - unsupported element \"%s\"\n",
				it->first.c_str());
		//! [RegistrySignalRange::_init]
	}
}

RegistrySignalDrug::RegistrySignalDrug(MedRepository &rep, const string &path_to_cfg_file) {
	repo = &rep;
	signalName = "";
	duration_flag = 0;
	buffer_duration = 0;
	take_only_first = false;
	outcome_value = 1;
	boost::filesystem::path p(path_to_cfg_file);
	base_path = p.parent_path().string();
}

void RegistrySignalDrug::_init(const map<string, string>& mapper) {

	string sets_path = "";
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RegistrySignalDrug::_init]
		if (it->first == "sets_path") //should contain "sets=" which points to file with list of codes with TAB min_dosage_range TAB max_dosage_range
			sets_path = it->second;
		else
			MTHROW_AND_ERR("Error in RegistrySignalDrug::init - unsupported element \"%s\"\n",
				it->first.c_str());
		//! [RegistrySignalDrug::_init]
	}
	//save to end
	if (boost::starts_with(sets_path, "path_rel:") && !base_path.empty()) {
		//use relative path to the current config file
		sets_path = sets_path.substr(9);
		//use base_paht as fix to relative path:
		if (!((sets_path.size() > 2 && (sets_path[0] == '/' || sets_path[1] == ':'))))
			sets_path = base_path + path_sep() + sets_path;
	}
	if (!sets_path.empty()) {
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

void RegistrySignalAnd::_init(const map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RegistrySignalAnd::_init]
		if (it->first == "conditions")  //not checking for infinite loop
			RegistrySignal::parse_registry_rules(it->second, *repo, conditions);
		else
			MTHROW_AND_ERR("ERROR in RegistrySignalAnd::init - Unsupported Argument %s\n", it->first.c_str());
		//! [RegistrySignalAnd::_init]
	}
	if (conditions.empty())
		MTHROW_AND_ERR("ERROR in RegistrySignalAnd::init - conditions is empty. please use conditions to reffer to file with and conditions of signals\n");
	//given signalName argument
	if (!signalName.empty())
		MWARN("Warning in RegistrySignalAnd::init - ignoring signalName argument. this is wrapper operation\n");
}

RegistrySignalAnd::~RegistrySignalAnd() {
	for (size_t i = 0; i < conditions.size(); ++i)
		delete conditions[i];
	conditions.clear();
}

bool RegistrySignalAny::get_outcome(const UniversalSigVec &s, int current_i, float &result) {
	bool is_active = current_i < s.len;
	if (is_active)
		result = outcome_value;
	return is_active;
}

inline int Date_wrapper(const UniversalSigVec &signal, int i) {
	if (signal.n_time_channels() != 0)
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
			if (r.end_date > r.start_date)
				results.push_back(r);

			//start new record
			//r.pid = pid;
			r.start_date = Date_wrapper(*signal, i);
			r.registry_value = registry_outcome_result;
			if (signal_prop->take_only_first) {
				r.end_date = max_date_mark;
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
					if (!seperate_cases) {
						r.end_date = medial::repository::DateAdd(Date_wrapper(*signal, i), signal_prop->duration_flag);
						max_search = medial::repository::DateAdd(r.end_date, signal_prop->buffer_duration - 1);
					}
					else {
						results.push_back(r);
						r.end_date = medial::repository::DateAdd(Date_wrapper(*signal, i), signal_prop->duration_flag);
						r.start_date = Date_wrapper(*signal, i);
						max_search = medial::repository::DateAdd(r.end_date, signal_prop->buffer_duration - 1);
					}
				}

				signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
				if (signal_index < 0)
					break;
				i = signals_indexes_pointers[signal_index] - 1; //the current signal time
				signal_prop = signal_filters[signal_index];
				signal = &usv[signal_index];
			}
			results.push_back(r);
			if (signal_index < 0) {
				r.start_date = max_date_mark; //no more control times, reached the end
				break;
			}
			//prepare for next:
			r.registry_value = 0;
			r.start_date = Date_wrapper(*signal, i); //already after duration and buffer. can start new control
			continue; //dont call fetch_next again
		}

		signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
	}

	r.end_date = last_date;
	if (last_date > 0) {
		last_date = medial::repository::DateAdd(last_date, -end_buffer_duration);
		if (r.end_date > r.start_date)
			results.push_back(r);
	}
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

/// it's not mcnemar it's more like: https://en.wikipedia.org/wiki/Cochran%E2%80%93Mantel%E2%80%93Haenszel_statistics
double medial::contingency_tables::calc_cmh_square_dist(const map<float, vector<int>> *gender_sorted, \
	const map<float, vector<int>> *gender_sorted2, bool &ok) {
	//calc over all ages
	double regScore = 0;
	double up = 0, down = 0;
	if (gender_sorted != NULL)
		for (auto i = gender_sorted->begin(); i != gender_sorted->end(); ++i) { //iterate over age bins
			const vector<int> &counts = i->second; //the forth numbers
			double totCnt = 0;
			vector<double> R(2), C(2);
			R[0] = counts[0] + counts[1]; //total no sig
			R[1] = counts[2 + 0] + counts[2 + 1]; //total with sig
			C[0] = counts[0] + counts[2]; //total controls
			C[1] = counts[0 + 1] + counts[2 + 1]; //total cases
			totCnt = R[0] + R[1];

			double b, c;
			b = counts[2 + 1];
			c = R[1] * C[1] / totCnt;
			if (totCnt < 1 || R[0] == 0 || R[1] == 0 || C[0] == 0 || C[1] == 0)
				continue;

			up += (b - c) * (b - c); //Cochran–Mantel–Haenszel statistics
			down += C[0] * C[1] * R[0] * R[1] / totCnt / totCnt / (totCnt - 1); //Cochran–Mantel–Haenszel statistics

		}
	if (gender_sorted2 != NULL)
		for (auto i = gender_sorted2->begin(); i != gender_sorted2->end(); ++i) { //iterate over age bins
			const vector<int> &counts = i->second; //the forth numbers
			double totCnt = 0;
			vector<double> R(2), C(2);
			R[0] = counts[0] + counts[1]; //total no sig
			R[1] = counts[2 + 0] + counts[2 + 1]; //total with sig
			C[0] = counts[0] + counts[2]; //total controls
			C[1] = counts[0 + 1] + counts[2 + 1]; //total cases
			totCnt = R[0] + R[1];

			double b, c;
			b = counts[2 + 1];
			c = R[1] * C[1] / totCnt;
			if (totCnt < 1 || R[0] == 0 || R[1] == 0 || C[0] == 0 || C[1] == 0)
				continue;

			up += (b - c) * (b - c); //Cochran–Mantel–Haenszel statistics
			down += C[0] * C[1] * R[0] * R[1] / totCnt / totCnt / (totCnt - 1); //Cochran–Mantel–Haenszel statistics
		}
	ok = false;
	if (down > 0) {
		ok = true;
		regScore += up / down; //Cochran–Mantel–Haenszel statistics
	}
	return regScore;
}

double medial::contingency_tables::calc_mcnemar_square_dist(const map<float, vector<int>> &gender_sorted) {
	//calc over all ages
	double regScore = 0;
	for (auto i = gender_sorted.begin(); i != gender_sorted.end(); ++i) { //iterate over age bins
		const vector<int> &counts = i->second; //the forth numbers
		vector<double> R(2);
		R[0] = counts[0] + counts[1]; //total no sig
		R[1] = counts[2 + 0] + counts[2 + 1]; //total with sig

		double p_min;
		double b, c;
		int min_ind = 0;
		//Mathcing to lower:
		p_min = R[0];
		if (R[1] < p_min) {
			p_min = R[1]; //matching to upper - almost always the case
			++min_ind;
		}
		if (p_min == 0)
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
			regScore += ((b - c) * (b - c)) / (b + c); //Mcnemar
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

void collect_stats(const map<float, map<float, vector<int>>> &male_stats,
	const map<float, map<float, vector<int>>> &female_stats, vector<int> &all_signal_values, vector<int> &signal_indexes,
	vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
	, vector<double> &pos_ratio) {
	unordered_set<int> all_vals;
	for (auto i = male_stats.begin(); i != male_stats.end(); ++i)
		all_vals.insert((int)i->first);
	for (auto i = female_stats.begin(); i != female_stats.end(); ++i)
		all_vals.insert((int)i->first);
	all_signal_values.insert(all_signal_values.end(), all_vals.begin(), all_vals.end());
	signal_indexes.resize(all_signal_values.size());
	for (size_t i = 0; i < signal_indexes.size(); ++i)
		signal_indexes[i] = (int)i;
	posCnts.resize(all_signal_values.size());
	valCnts.resize(all_signal_values.size());
	lift.resize(all_signal_values.size());
	pos_ratio.resize(all_signal_values.size());

	for (int index : signal_indexes)
	{
		float signalVal = all_signal_values[index];
		//check chi-square for this value:
		double totCnt = 0, weighted_lift = 0;
		/*unordered_map<float, double> prior_lift_males, prior_lift_females; //for each age_bin - prior outcome

		if (male_stats.find(signalVal) != male_stats.end())
			for (auto jt = male_stats.at(signalVal).begin(); jt != male_stats.at(signalVal).end(); ++jt)
				if (jt->second[1] + jt->second[1 + 2] > 0)
					prior_lift_males[jt->first] = (jt->second[1] + jt->second[1 + 2])
					/ (jt->second[1] + jt->second[1 + 2] + jt->second[0] + jt->second[0 + 2]);
		if (female_stats.find(signalVal) != female_stats.end())
			for (auto jt = female_stats.at(signalVal).begin(); jt != female_stats.at(signalVal).end(); ++jt)
				if (jt->second[1] + jt->second[1 + 2] > 0)
					prior_lift_females[jt->first] = (jt->second[1] + jt->second[1 + 2])
					/ (jt->second[1] + jt->second[1 + 2] + jt->second[0] + jt->second[0 + 2]);*/


		if (male_stats.find(signalVal) != male_stats.end())
			for (auto jt = male_stats.at(signalVal).begin(); jt != male_stats.at(signalVal).end(); ++jt) {
				totCnt += jt->second[2] + jt->second[3];
				posCnts[index] += jt->second[1 + 2];
				if (jt->second[1 + 0] > 0)
					weighted_lift += double(jt->second[1 + 2]) / double(jt->second[1 + 0]) * double(jt->second[1 + 0] + jt->second[0 + 0]);
				//weighted_lift += jt->second[1 + 2] / prior_lift_males[jt->first];
			}
		if (female_stats.find(signalVal) != female_stats.end())
			for (auto jt = female_stats.at(signalVal).begin(); jt != female_stats.at(signalVal).end(); ++jt) {
				totCnt += jt->second[2] + jt->second[3];
				posCnts[index] += jt->second[1 + 2];
				if (jt->second[1 + 0] > 0)
					weighted_lift += double(jt->second[1 + 2]) / double(jt->second[1 + 0]) * double(jt->second[1 + 0] + jt->second[0 + 0]);
			}
		if (totCnt == 0)
			continue;
		valCnts[index] = totCnt; //for signal apeareance
		lift[index] = weighted_lift / totCnt;

		pos_ratio[index] = posCnts[index] / totCnt;
	}
}

void medial::contingency_tables::calc_chi_scores(const map<float, map<float, vector<int>>> &male_stats,
	const map<float, map<float, vector<int>>> &female_stats,
	vector<int> &all_signal_values, vector<int> &signal_indexes,
	vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
	, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio, int smooth_balls
	, float allowed_error, int minimal_balls) {

	collect_stats(male_stats, female_stats, all_signal_values, signal_indexes, valCnts, posCnts, lift, pos_ratio);
	scores.resize(all_signal_values.size());
	p_values.resize(all_signal_values.size());

	for (int index : signal_indexes)
	{
		float signalVal = all_signal_values[index];
		//check chi-square for this value:
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

/// it's not mcnemar it's more like: https://en.wikipedia.org/wiki/Cochran%E2%80%93Mantel%E2%80%93Haenszel_statistics
void medial::contingency_tables::calc_cmh_scores(const map<float, map<float, vector<int>>> &male_stats,
	const map<float, map<float, vector<int>>> &female_stats,
	vector<int> &all_signal_values, vector<int> &signal_indexes,
	vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
	, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio) {

	collect_stats(male_stats, female_stats, all_signal_values, signal_indexes, valCnts, posCnts, lift, pos_ratio);
	scores.resize(all_signal_values.size());
	p_values.resize(all_signal_values.size());

	for (int index : signal_indexes)
	{
		float signalVal = all_signal_values[index];
		//check chi-square for this value:

		double regScore = 0;
		const map<float, vector<int>> *p1 = NULL, *p2 = NULL;
		if (male_stats.find(signalVal) != male_stats.end())
			p1 = &male_stats.at(signalVal);
		if (female_stats.find(signalVal) != female_stats.end())
			p2 = &female_stats.at(signalVal);
		bool is_ok;
		regScore = calc_cmh_square_dist(p1, p2, is_ok);
		if (is_ok) {
			scores[index] = (float)regScore;
			int dof = 1;
			double pv = chisqr(dof, regScore);
			p_values[index] = pv;
		}
	}
}

void medial::contingency_tables::calc_mcnemar_scores(const map<float, map<float, vector<int>>> &male_stats,
	const map<float, map<float, vector<int>>> &female_stats,
	vector<int> &all_signal_values, vector<int> &signal_indexes,
	vector<double> &valCnts, vector<double> &posCnts, vector<double> &lift
	, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio) {

	collect_stats(male_stats, female_stats, all_signal_values, signal_indexes, valCnts, posCnts, lift, pos_ratio);
	scores.resize(all_signal_values.size());
	p_values.resize(all_signal_values.size());

	for (int index : signal_indexes)
	{
		float signalVal = all_signal_values[index];
		//check chi-square for this value:
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

void medial::contingency_tables::filterHirarchy(const map<int, vector<int>> &member2Sets, const map<int, vector<int>> &set2Members,
	vector<int> &indexes, const vector<int> &signal_values, const vector<double> &pVals,
	const vector<double> &valCnts, const vector<double> &lifts, const unordered_map<int, double> &code_unfiltered_cnts,
	float pValue_diff, float lift_th, float count_similarity, float child_fitlered_ratio,
	const map<int, vector<string>> *categoryId_to_name) {
	double minPerc = lift_th; //max allowed diff in lift

	unordered_map<int, double> pVals_d, lifts_d, valCnts_d;
	for (int index : indexes)
	{
		pVals_d[signal_values[index]] = pVals[index];
		lifts_d[signal_values[index]] = lifts[index];
		valCnts_d[signal_values[index]] = valCnts[index];
	}
	unordered_set<int> signal_values_set(signal_values.begin(), signal_values.end());

	unordered_set<int> toRemove;
	unordered_set<int> removingParents;
	for (int index : indexes)
	{
		int keyVal = signal_values[index];
		double currCnt = valCnts[index];
		double pV = pVals[index];
		double currLift = lifts[index];
		if (member2Sets.find(keyVal) == member2Sets.end())
			continue; // no parents

		vector<int> parents = member2Sets.at(keyVal);
		vector<int> filtered;
		for (size_t i = 0; i < parents.size(); ++i)
			if (signal_values_set.find(parents[i]) != signal_values_set.end())
				filtered.push_back(parents[i]);
		filtered.swap(parents);

		for (int parentId : parents) { //has parents with similarity, so remove child
			double parentLift = 1, parentCnt = 0, parentPv = 1;
			if (lifts_d.find(parentId) != lifts_d.end()) {
				parentLift = lifts_d.at(parentId);
				parentCnt = valCnts_d.at(parentId);
				parentPv = pVals_d.at(parentId);
			}
			if (pVals_d.find(parentId) != pVals_d.end()) { //has parent to compare to
				if (abs(parentCnt - currCnt) / parentCnt <= count_similarity) { //less than diff, remove child:
					if (categoryId_to_name != NULL)
						MLOG("DEBUG: remove key %s, parent has similar count:%d and current:%d\n",
							categoryId_to_name->at(keyVal).back().c_str(), (int)parentCnt, (int)currCnt);
					toRemove.insert(keyVal);
					removingParents.insert(parentId);
					break;
				}

				if (abs(parentPv - pV) <= pValue_diff) {
					double cmp = -1;
					if (currLift > 0 && parentLift > currLift)
						cmp = abs(parentLift / currLift - 1);
					else if (parentLift > 0 && currLift > parentLift)
						cmp = abs(currLift / parentLift - 1);
					if ((parentLift == 0 && currLift == 0) || (cmp > 0 && cmp <= minPerc)) { //less than 5% diff, remove child:
						if (categoryId_to_name != NULL)
							MLOG("DEBUG: remove key %s, parent has similar lift:%2.3f and current:%2.3f\n",
								categoryId_to_name->at(keyVal).back().c_str(), parentLift, currLift);
						toRemove.insert(keyVal);
						removingParents.insert(parentId);
						break;
					}
				}
			}
		}
	}

	vector<int> keep_indexes;
	keep_indexes.reserve(indexes.size());
	for (int index : indexes)
		if (toRemove.find(signal_values[index]) == toRemove.end())
			keep_indexes.push_back(index);
	indexes.swap(keep_indexes);

	//remove parents that has childs removed( beyond percantage) and child that are kept (means the parent is not 
	// right resolution, so remove
	toRemove.clear();
	unordered_set<int> filter_idx;
	for (int index : indexes)
		filter_idx.insert(signal_values[index]);

	for (int index : indexes)
	{
		int keyVal = signal_values[index];

		// If this parent has caused the removal of one of it's children, we can't remove it !
		if (removingParents.find(keyVal) != removingParents.end())
			continue;

		double currCnt = valCnts[index];
		//test if that's parent that need to be removed - has child that has been removed, 
		//and at least onr child that haven't moved:
		if (set2Members.find(keyVal) == set2Members.end())
			continue; // no childs - it's leaf
		const vector<int> &childs = set2Members.at(keyVal);
		double removed_child_counts = 0;
		bool has_keep_child = false;
		stringstream log_desc;
		for (int child_id : childs) {
			if (filter_idx.find(child_id) != filter_idx.end()) {
				has_keep_child = true;
				if (categoryId_to_name != NULL)
					log_desc << "\n\t###" << "Kept_Child:" << child_id << "(" << categoryId_to_name->at(child_id).back() << ")"
					<< " Count:" << valCnts_d.at(child_id);
			}
			else if (code_unfiltered_cnts.find(child_id) != code_unfiltered_cnts.end()) {
				removed_child_counts += code_unfiltered_cnts.at(child_id);
				if (categoryId_to_name != NULL)
					log_desc << "\n\t###" << "Removed_Child:" << child_id << "(" << categoryId_to_name->at(child_id).back() << ")"
					<< " Count:" << code_unfiltered_cnts.at(child_id);
			}
		}
		if (has_keep_child && removed_child_counts > 0 && removed_child_counts / currCnt >= child_fitlered_ratio) {
			if (categoryId_to_name != NULL)
				MLOG("DEBUG: remove parent code %d(%s) with %d, has big removed childs:%s\n",
					keyVal, categoryId_to_name->at(keyVal).back().c_str(), (int)currCnt, log_desc.str().c_str());
			toRemove.insert(keyVal);
		}
	}

	keep_indexes.clear();
	for (int index : indexes)
		if (toRemove.find(signal_values[index]) == toRemove.end())
			keep_indexes.push_back(index);
	indexes.swap(keep_indexes);
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
	if (num_of_init_test == 0)
		return;
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
	int stop_index = 0;
	for (unsigned int i = 0; i < keysSorted.size(); i++) {
		if (keysSorted[i].second[0] <= normAlpha * (i + 1))
			stop_index = i + 1;
	}

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

RegistrySignal *RegistrySignal::make_registry_signal(const string &type, MedRepository &rep,
	const string &path_to_cfg_file) {
	vector<string> empty_arr;
	string empty_str = "";
	//! [RegistrySignal::make_registry_signal]
	if (type == "set")
		return new RegistrySignalSet(empty_str, 0, 0, false, rep, empty_arr, path_to_cfg_file);
	else if (type == "range")
		return new RegistrySignalRange(empty_str, 0, 0, false, 0, 0);
	else if (type == "drug")
		return new RegistrySignalDrug(rep, path_to_cfg_file);
	else if (type == "and")
		return new RegistrySignalAnd(rep);
	else if (type == "any")
		return new RegistrySignalAny;
	else
		MTHROW_AND_ERR("Error: Unsupported type \"%s\" for RegistrySignal::make_registry_signal\n", type.c_str());
	//! [RegistrySignal::make_registry_signal]
}

RegistrySignal *RegistrySignal::make_registry_signal(const string &type, MedRepository &rep,
	const string &init_string, const string &path_to_cfg_file) {
	RegistrySignal *reg = make_registry_signal(type, rep, path_to_cfg_file);
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
		RegistrySignal *reg = RegistrySignal::make_registry_signal(type, rep, ini, reg_cfg);
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
		else if (it->first == "seperate_cases")
			seperate_cases = stoi(it->second) > 0;
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
	else if (registry_type == "keep_alive")
		registry = new MedRegistryKeepAlive;
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

void MedRegistry::merge_records()
{
	vector<MedRegistryRecord> merged_records;
	merged_records.push_back(registry_records[0]);

	for (int i = 1; i < registry_records.size(); i++)
	{
		if ((registry_records[i].pid == merged_records.back().pid) && (registry_records[i].start_date <= merged_records.back().end_date + 1) && (registry_records[i].registry_value == merged_records.back().registry_value))
		{
			merged_records.back().end_date = registry_records[i].end_date;
		}
		else
			merged_records.push_back(registry_records[i]);
	}
	registry_records = move(merged_records);
}

MedRegistry *MedRegistry::create_registry_full(const string &registry_type, const string &init_str,
	const string &repository_path, MedModel &model_with_rep_processor, medial::repository::fix_method method) {
	MedPidRepository rep;
	if (rep.init(repository_path) < 0)
		MTHROW_AND_ERR("Can't read repository %s", repository_path.c_str());
	model_with_rep_processor.collect_and_add_virtual_signals(rep);

	MedRegistry *registry = MedRegistry::make_registry(registry_type, rep, init_str);

	vector<int> pids;
	vector<string> sig_codes;
	registry->get_pids(pids);
	registry->get_registry_creation_codes(sig_codes);
	vector<string> physical_signal;
	medial::repository::prepare_repository(rep, sig_codes, physical_signal, &model_with_rep_processor.rep_processors);

	if (rep.read_all(repository_path, pids, physical_signal) < 0)
		MTHROW_AND_ERR("Can't read repository %s", repository_path.c_str());
	registry->create_registry(rep, method, &model_with_rep_processor.rep_processors);
	registry->clear_create_variables();

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
	if (pid_to_max_allowed.find(pid) != pid_to_max_allowed.end()) {
		if (pid_to_max_allowed[pid] > 0)
			max_allowed_date = pid_to_max_allowed[pid];
		else
			return;
	}

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

						if (r.end_date > r.start_date)
							results.push_back(r);
						//start new record with 0 outcome:
						//r.registry_value = signal_prop_outcome; //left the same no need
						r.start_date = curr_date; //continue from where left
						//r.end_date = medial::repository::DateAdd(r.start_date, 1); //let's start from 1 day
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

					if (!mark_no_match && r.end_date > r.start_date)
						results.push_back(r);

					//skip if may not use
					if (signal_prop->take_only_first && outcomes_may_not_use.find(signal_prop_outcome) != outcomes_may_not_use.end())
						continue;
					//start new record
					//r.pid = pid;
					//r.min_allowed_date = min_date;
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
				if (!mark_no_match && r.end_date > r.start_date)
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
	if (last_date > 0) {
		last_date = medial::repository::DateAdd(last_date, -end_buffer_duration);
		if (r.end_date > r.start_date && !mark_no_match)
			results.push_back(r);
	}
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

int MedRegistryKeepAlive::init(map<string, string>& map) {
	string repository_path;
	for (auto it = map.begin(); it != map.end(); ++it) {
		//! [MedRegistryKeepAlive::init]
		if (it->first == "duration")
			duration = med_stoi(it->second);
		else if (it->first == "rep")
			repository_path = it->second;
		else if (it->first == "max_repo_date")
			max_repo_date = med_stoi(it->second);
		else if (it->first == "secondry_start_buffer_duration")
			secondry_start_buffer_duration = med_stoi(it->second);
		else if (it->first == "start_buffer_duration")
			start_buffer_duration = med_stoi(it->second);
		else if (it->first == "end_buffer_duration")
			end_buffer_duration = med_stoi(it->second);
		else if (it->first == "signal_list")
			boost::split(signal_list, it->second, boost::is_any_of(",;"));
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
		else
			MTHROW_AND_ERR("Error in MedRegistryKeepAlive::init - unknown parameter %s\n", it->first.c_str());
		//! [MedRegistryKeepAlive::init]
	}

	signalCodes_names = signal_list;
	signalCodes_names.push_back("BDATE");
	for (size_t i = 0; i < signal_list.size(); ++i)
		if (signal_list[i] == "BDATE") {
			need_bdate = true;
			break;
		}

	return 0;
}

void MedRegistryKeepAlive::get_registry_records(int pid, int bdate, vector<UniversalSigVec_mem> &usv, vector<MedRegistryRecord> &results) {
	vector<int> signals_indexes_pointers(usv.size()); //all in 0

	int max_allowed_date = max_repo_date;
	if (pid_to_max_allowed.find(pid) != pid_to_max_allowed.end())
		max_allowed_date = pid_to_max_allowed[pid];

	MedRegistryRecord r;
	r.pid = pid;
	r.registry_value = 1; //we mark only the legal time range for sampling
	r.start_date = 0;
	int start_date = -1, last_date = -1;
	int signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
	while (signal_index >= 0)
	{
		const UniversalSigVec *signal = &usv[signal_index];
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
		if (r.start_date <= 0) {
			r.start_date = min_date;
			r.end_date = medial::repository::DateAdd(start_date, duration);
		}
		//I have start_date
		int curr_date = Date_wrapper(*signal, i);
		if (max_allowed_date != 0 && curr_date > max_allowed_date)
			break;

		int new_start = medial::repository::DateAdd(curr_date, secondry_start_buffer_duration);
		if (curr_date < r.end_date || (secondry_start_buffer_duration < 0 && new_start < r.end_date))
			r.end_date = medial::repository::DateAdd(curr_date, duration);
		else {
			//has dead region - close buffer and open new one:
			if (r.end_date > r.start_date)
				results.push_back(r);
			r.start_date = new_start;
			r.end_date = medial::repository::DateAdd(curr_date, duration);
		}

		last_date = curr_date;
		signal_index = medial::repository::fetch_next_date(usv, signals_indexes_pointers);
	}

	if (last_date > 0) {
		r.end_date = medial::repository::DateAdd(last_date, duration - end_buffer_duration);
		if (start_date > 0 && r.end_date > r.start_date)
			results.push_back(r);
	}
}

void medial::registry::complete_active_period_as_controls(vector<MedRegistryRecord> &registry,
	const vector<MedRegistryRecord> &active_periods_registry, bool unite_full_controls) {
	unordered_map<int, vector<const MedRegistryRecord *>> pid_to_periods;
	for (size_t i = 0; i < active_periods_registry.size(); ++i)
		pid_to_periods[active_periods_registry[i].pid].push_back(&active_periods_registry[i]);
	unordered_map<int, vector<MedRegistryRecord>> pid_to_regs;
	for (size_t i = 0; i < registry.size(); ++i)
		pid_to_regs[registry[i].pid].push_back(registry[i]);
	vector<MedRegistryRecord> new_reg;

	if (!unite_full_controls) {
		for (auto &rec : pid_to_regs)
		{
			int pid = rec.first;
			if (pid_to_periods.find(pid) == pid_to_periods.end())
				continue; //not found - no completion
			vector<MedRegistryRecord> &reg_recs = rec.second;
			vector<const MedRegistryRecord *> &active_pr = pid_to_periods.at(pid);
			sort(reg_recs.begin(), reg_recs.end(), [](const MedRegistryRecord &a, const MedRegistryRecord &b)
			{ return a.start_date < b.start_date; });
			sort(active_pr.begin(), active_pr.end(), [](const MedRegistryRecord *a, const MedRegistryRecord *b)
			{ return a->start_date < b->start_date; });

			//both sorted - now "join" them
			int active_i = 0, reg_i = 0;
			int curr_time = 0;
			vector<MedRegistryRecord> added_recs;
			while (active_i < active_pr.size())
			{
				const MedRegistryRecord *curr_active = active_pr[active_i];
				const MedRegistryRecord *curr_reg = NULL;
				if (reg_i < reg_recs.size())
					curr_reg = &reg_recs[reg_i];
				if (curr_reg == NULL) {
					if (curr_time == 0)
						curr_time = curr_active->start_date;
					MedRegistryRecord reg_rec;
					reg_rec.pid = pid;
					reg_rec.registry_value = 0;
					reg_rec.start_date = curr_time;
					reg_rec.end_date = curr_active->end_date;
					if (reg_rec.end_date > reg_rec.start_date) // add if not equal:
						added_recs.push_back(reg_rec);

					curr_time = 0;
					++active_i; //finished active period
					continue;
				}
				if (curr_time > 0 && curr_time > curr_active->end_date) {
					++active_i;
					continue;
				}

				if (curr_active->start_date < curr_reg->start_date) {
					if (curr_time == 0 || curr_time < curr_active->start_date)
						curr_time = curr_active->start_date;
					//complete till active_end or start_date of reg:
					MedRegistryRecord reg_rec;
					reg_rec.pid = pid;
					reg_rec.registry_value = 0;
					reg_rec.start_date = curr_time;
					reg_rec.end_date = curr_active->end_date;
					if (curr_reg->start_date < reg_rec.end_date) {
						reg_rec.end_date = curr_reg->start_date;
						curr_time = curr_reg->end_date; //update to current reg time - there is intersection, skip till curr_reg.end_date
						++reg_i; //read and skip reg record can move on...
						if (curr_reg->end_date > curr_active->end_date)
							++active_i; //finished active period
					}
					else {
						++active_i; //finished active period
						curr_time = 0; //need to test again for curr_time
					}

					if (reg_rec.end_date > reg_rec.start_date) // add if not equal:
						added_recs.push_back(reg_rec);
				}
				else {
					curr_time = curr_reg->end_date;
					++reg_i;
					//for efficancy - will work anyway
					if (curr_active->end_date < curr_reg->end_date)
						++active_i;
				}

			}

			//add and sort new control records:
			reg_recs.insert(reg_recs.end(), added_recs.begin(), added_recs.end());
			sort(reg_recs.begin(), reg_recs.end(), [](const MedRegistryRecord &a, const MedRegistryRecord &b)
			{ return a.start_date < b.start_date; });
		}

		//add as controls missings from reg:
		for (auto &rec : active_periods_registry)
			if (pid_to_regs.find(rec.pid) == pid_to_regs.end()) {
				MedRegistryRecord new_rec;
				new_rec.pid = rec.pid;
				new_rec.registry_value = 0;
				new_rec.start_date = rec.start_date;
				new_rec.end_date = rec.end_date;
				pid_to_regs[rec.pid].push_back(new_rec);
			}

		//commit to reg with unite:

		MedRegistryRecord rec_temp;
		for (const auto rec : pid_to_regs) {
			//sort(rec.second.begin(), rec.second.end(), [](const MedRegistryRecord &a, const MedRegistryRecord &b)
			//{ return a.start_date < b.start_date; });
			//unite control times:
			if (!rec.second.empty()) {
				rec_temp.pid = -1;
				for (size_t i = 0; i < rec.second.size(); ++i)
				{
					if (rec.second[i].registry_value > 0) {
						//close buffer:
						if (rec_temp.pid > 0) {
							new_reg.push_back(rec_temp);
							rec_temp.pid = -1;
						}
						new_reg.push_back(rec.second[i]);
					}
					else {
						if (rec_temp.pid == -1) {
							rec_temp = rec.second[i];
						}
						else if (rec.second[i].start_date <= rec_temp.end_date) {
							//unite buffers
							rec_temp.end_date = rec.second[i].end_date;
						}
						else {
							//close buffer & start new one:
							if (rec_temp.pid > 0) {
								new_reg.push_back(rec_temp);
								rec_temp = rec.second[i];
							}
						}
					}
				}
				//close buffer:
				if (rec_temp.pid > 0)
					new_reg.push_back(rec_temp);

				//new_reg.insert(new_reg.end(), rec.second.begin(), rec.second.end());
			}
		}
	}
	else { // just add all as controls
		for (auto &rec : active_periods_registry) {
			MedRegistryRecord new_rec;
			new_rec.pid = rec.pid;
			new_rec.registry_value = 0;
			new_rec.start_date = rec.start_date;
			new_rec.end_date = rec.end_date;
			pid_to_regs[rec.pid].push_back(new_rec);
		}

		for (auto &rec : pid_to_regs) {
			sort(rec.second.begin(), rec.second.end(), [](const MedRegistryRecord &a, const MedRegistryRecord &b)
			{ return a.start_date < b.start_date; });
			new_reg.insert(new_reg.end(), rec.second.begin(), rec.second.end());
		}
	}




	registry = move(new_reg);
}

bool MedRegistry::get_pid_records(PidRec &rec, int bDateCode, const vector<int> &used_sigs,
	vector<MedRegistryRecord> &results) {

	vector<UniversalSigVec_mem> sig_vec((int)used_sigs.size());
	bool pid_fix = false;
	for (size_t k = 0; k < sig_vec.size(); ++k) {
		UniversalSigVec vv;
		rec.uget(used_sigs[k], vv);
		bool did_something = medial::repository::fix_contradictions(vv, medial::repository::fix_method::none, sig_vec[k]);
		if (did_something)
			pid_fix = did_something;
	}
	int birth = medial::repository::get_value(rec, bDateCode);

	get_registry_records(rec.pid, birth, sig_vec, results);
	return pid_fix;
}