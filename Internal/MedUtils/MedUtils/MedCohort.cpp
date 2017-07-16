#include "MedCohort.h"

#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string.hpp>
#include <InfraMed/InfraMed/InfraMed.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <fstream>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=====================================================================================
// CohortRec
//=====================================================================================
//-------------------------------------------------------------------------------------
int CohortRec::init(map<string, string>& map)
{
	for (auto &m : map) {
		if (m.first == "pid") pid = stoi(m.second);
		else if (m.first == "from") from = stoi(m.second);
		else if (m.first == "to") to = stoi(m.second);
		else if (m.first == "outcome_date") outcome_date = stoi(m.second);
		else if (m.first == "outcome") outcome = stof(m.second);
		else {
			MERR("Unknown variable %s in CohortRec\n", m.first.c_str());
		}

	}
	return 0;
}

//-------------------------------------------------------------------------------------
int CohortRec::get_string(string &to_str)
{
	to_str = to_string(pid) + "\t" + to_string(from) + "\t" + to_string(to) + "\t" + to_string(outcome_date) + "\t" + to_string(outcome);
	return 0;
}

//-------------------------------------------------------------------------------------
int CohortRec::from_string(string &from_str)
{
	vector<string> fields;
	boost::split(fields, from_str, boost::is_any_of("\t"));

	if (fields.size() >= 5) {
		pid = stoi(fields[0]);
		from = stoi(fields[1]);
		to = stoi(fields[2]);
		outcome_date = stoi(fields[3]);
		outcome = stof(fields[4]);
	}
	else
		return -1;

	return 0;
}

//=====================================================================================
// SamplingParams
//=====================================================================================
//-------------------------------------------------------------------------------------
int SamplingParams::init(map<string, string>& map)
{
	for (auto &m : map) {
		if (m.first == "min_control" || m.first == "min_control_years") min_control_years = stof(m.second);
		else if (m.first == "max_control" || m.first == "max_control_years") max_control_years = stof(m.second);
		else if (m.first == "min_case" || m.first == "min_case_years") min_case_years = stof(m.second);
		else if (m.first == "max_case" || m.first == "max_case_years") max_case_years = stof(m.second);
		else if (m.first == "is_continous") is_continous = stoi(m.second);
		else if (m.first == "min_days_from_outcome" || m.first == "min_days") min_days_from_outcome = stoi(m.second);
		else if (m.first == "jump_days") jump_days = stoi(m.second);
		else if (m.first == "min_year") min_year = stoi(m.second);
		else if (m.first == "max_year") max_year = stoi(m.second);
		else if (m.first == "gender_mask") gender_mask = stoi(m.second);
		else if (m.first == "train_mask") train_mask = stoi(m.second);
		else if (m.first == "min_age") min_age = stoi(m.second);
		else if (m.first == "max_age") max_age = stoi(m.second);
		else if (m.first == "rep") rep_fname = m.second;
		else if (m.first == "stick_to" || m.first == "stick_to_sigs") {
			boost::split(stick_to_sigs, m.second, boost::is_any_of(","));
		}
		else {
			MERR("Unknown variable %s in SamplingParams\n", m.first.c_str());
		}

	}
	return 0;
}

//=====================================================================================
// IncidenceParams
//=====================================================================================
//-------------------------------------------------------------------------------------
int IncidenceParams::init(map<string, string>& map)
{
	for (auto &m : map) {
		if (m.first == "age_bin") age_bin = stoi(m.second);
		else if (m.first == "min_samples_in_bin") min_samples_in_bin = stoi(m.second);
		else if (m.first == "from_year") from_year = stoi(m.second);
		else if (m.first == "to_year") to_year = stoi(m.second);
		else if (m.first == "gender_mask") gender_mask = stoi(m.second);
		else if (m.first == "train_mask") train_mask = stoi(m.second);
		else if (m.first == "from_age") from_age = stoi(m.second);
		else if (m.first == "to_age") to_age = stoi(m.second);
		else if (m.first == "rep") rep_fname = m.second;
		else {
			MERR("Unknown variable %s in IncidenceParams\n", m.first.c_str());
		}

	}
	return 0;
}

//=====================================================================================
// MedCohort
//=====================================================================================
//-------------------------------------------------------------------------------------
int MedCohort::read_from_file(string fname)
{
	ifstream inf(fname);

	MLOG("MedCohort: reading %s\n", fname.c_str());
	if (!inf) {
		MERR("MedCohort: can't open file %s for read\n", fname.c_str());
		return -1;
	}

	string curr_line;

	recs.clear();
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			CohortRec cr;
			if (cr.from_string(curr_line) < 0) {
				MERR("Bad Cohort Line: %s\n", curr_line.c_str());
				inf.close();
				return -1;
			}

			recs.push_back(cr);
		}
	}
	inf.close();
	MLOG("Read %d cohort lines from %s\n", recs.size(), fname.c_str());
	return 0;
}

//-------------------------------------------------------------------------------------
int MedCohort::write_to_file(string fname)
{
	ofstream of(fname);

	MLOG("MedCohort: writing to %s\n", fname.c_str());
	if (!of) {
		MERR("MedCohort: can't open file %s for writing\n", fname.c_str());
		return -1;
	}

	for (auto &rc: recs) {
		string sout;
		rc.get_string(sout);
		of << sout << endl;
	}

	MLOG("wrote [%d] records in cohort file %s\n", recs.size(), fname.c_str());
	of.close();
	return 0;
}

//-------------------------------------------------------------------------------------
int MedCohort::get_pids(vector<int> &pids)
{
	pids.clear();
	for (auto &cr : recs) pids.push_back(cr.pid);
	return 0;
}

//-------------------------------------------------------------------------------------
int MedCohort::create_incidence_file(IncidenceParams &i_params, string out_file)
{
	//string inc_params; // from_to_fname,pids_to_use_fname,from_year,to_year,min_age,max_age,bin_size,inc_file

	vector<int> train_to_take ={ 0,0,0,0 };
	if (i_params.train_mask &0x1) train_to_take[1] = 1;
	if (i_params.train_mask &0x2) train_to_take[2] = 1;
	if (i_params.train_mask &0x4) train_to_take[3] = 1;

	vector<int> pids;
	get_pids(pids);

	// read byears, gender and TRAIN
	MedRepository rep;
	if (rep.read_all(i_params.rep_fname, pids, { "BYEAR", "GENDER", "TRAIN" }) < 0) {
		MERR("FAILED reading repository %s\n", i_params.rep_fname.c_str());
		return -1;
	}

	// actual sampling
	unordered_set<int> pids_set;
	for (auto pid : pids) pids_set.insert(pid);

	map<int, pair<int, int>> counts;

	// age bins init
	for (int i = 0; i < 200; i += i_params.age_bin) counts[i] = pair<int, int>(0, 0);


	int byear_sid = rep.sigs.sid("BYEAR");
	int gender_sid = rep.sigs.sid("GENDER");
	int train_sid = rep.sigs.sid("TRAIN");
	int len;

	vector<int> all_cnts ={ 0,0 };

	for (auto &crec : recs) {
		int fyear = crec.from / 10000;
		int to_date = crec.to;
		if (crec.outcome != 0) to_date = crec.outcome_date;
		int tyear = to_date / 10000;
		int byear = (int)((((SVal *)rep.get(crec.pid, byear_sid, len))[0]).val);
		int gender = (int)((((SVal *)rep.get(crec.pid, gender_sid, len))[0]).val);
		int train = (int)((((SVal *)rep.get(crec.pid, train_sid, len))[0]).val);

		if ((gender & i_params.gender_mask) && (train_to_take[train]))
			for (int year = fyear; year <= tyear; year++) {
				if (year >= i_params.from_year && year <= i_params.to_year) {
					int age = year - byear;
					int bin = i_params.age_bin*(age / i_params.age_bin);
					counts[bin].first++; all_cnts[0]++;
					if (year == tyear && (crec.outcome!=0)) { counts[bin].second++; all_cnts[1]++; }
				}
			}
	}

	MLOG("Total counts: 0: %d 1: %d : inc %f\n", all_cnts[0], all_cnts[1], (float)all_cnts[1] / all_cnts[0]);

	int nlines = 0;
	for (auto &c : counts) {

		int age = c.first;
		int n0 = c.second.first;
		int n1 = c.second.second;

		if (age >= i_params.from_age && age <= i_params.to_age) nlines++;

		if (n0 > 0)
			MLOG("Ages: %d - %d : %d : 0: %d 1: %d : %f\n", age, age + i_params.age_bin, age + i_params.age_bin / 2, n0, n1, (n0 > 0) ? (float)n1 / n0 : 0);
	}

	ofstream of(out_file);

	of << "KeySize 1\n";
	of << "Nkeys " << nlines << "\n";
	of << "1.0\n";

	for (auto &c : counts) {

		int age = c.first;
		int n0 = c.second.first;
		int n1 = c.second.second;

		if (age >= i_params.from_age && age <= i_params.to_age) {
			of << age + i_params.age_bin / 2 << " " << n1 << " " << n0 - n1 << "\n";
		}

	}

	of.close();

	return 0;
}


//-------------------------------------------------------------------------------------
int MedCohort::create_sampling_file(SamplingParams &s_params, string out_sample_file)
{
	if (s_params.is_continous == 0)
		return create_sampling_file_sticked(s_params, out_sample_file);

	vector<int> train_to_take ={ 0,0,0,0 };
	if (s_params.train_mask &0x1) train_to_take[1] = 1;
	if (s_params.train_mask &0x2) train_to_take[2] = 1;
	if (s_params.train_mask &0x4) train_to_take[3] = 1;
	vector<int> pids;
	get_pids(pids);
	MedRepository rep;
	if (rep.read_all(s_params.rep_fname, pids, { "BYEAR", "GENDER", "TRAIN" }) < 0) {
		MERR("FAILED reading repository %s\n", s_params.rep_fname.c_str());
		return -1;
	}

	MedSamples samples;
	int byear_sid = rep.sigs.sid("BYEAR");
	int gender_sid = rep.sigs.sid("GENDER");
	int train_sid = rep.sigs.sid("TRAIN");
	int len;

	int nsamp = 0;

	for (auto &rc : recs) {
		int byear = (int)((((SVal *)rep.get(rc.pid, byear_sid, len))[0]).val);
		int gender = (int)((((SVal *)rep.get(rc.pid, gender_sid, len))[0]).val);
		int train = (int)((((SVal *)rep.get(rc.pid, train_sid, len))[0]).val);

		//MLOG("s: %d outcome %d %d from-to %d %d byear %d gender %d (mask %d) train %d (mask %d)\n", rc.pid, (int)rc.outcome, rc.outcome_date, rc.from, rc.to, byear, gender, s_params.gender_mask, train, s_params.train_mask);
		if ((gender & s_params.gender_mask) && (train_to_take[train])) {

			//MLOG("pid %d passed masks\n", rc.pid);
			if (rc.from <= rc.outcome_date && rc.outcome_date <= rc.to && rc.from <= rc.to) {

				//MLOG("pid %d passed from-to\n", rc.pid);
				// first moving to work with days
				int to_days = med_time_converter.convert_date(MedTime::Days, rc.to);
				int from_days = med_time_converter.convert_date(MedTime::Days, rc.from);
				int outcome_days = med_time_converter.convert_date(MedTime::Days, rc.outcome_date);

				// then - adjusting from - to to be within the frame for outcomes 0 and 1
				if (rc.outcome != 0) {
					// case
					to_days = outcome_days - int(365.0f * s_params.min_case_years);
					from_days = max(from_days, from_days - int(365.0f * s_params.max_case_years));
				}
				else {
					// control
					to_days = outcome_days - int(365.0f * s_params.min_control_years);
					from_days = max(from_days, from_days - int(365.0f * s_params.max_control_years));
				}

				to_days = min(to_days, outcome_days - s_params.min_days_from_outcome);

				int delta = to_days - from_days;
				//MLOG("pid %d to_days %d outcome_days %d from_days %d delta %d\n", rc.pid, to_days, outcome_days, from_days, delta);

				MedIdSamples mis;
				mis.id = rc.pid;

				while (delta >= 0) {

					int range = min(delta, s_params.jump_days);
					int r = rand_N(range);

					int rand_date = to_days - r;
					rand_date = med_time_converter.convert_days(MedTime::Date, rand_date);
					MedSample ms;
					ms.id = rc.pid;
					ms.outcome = rc.outcome;
					ms.outcomeTime = rc.outcome_date;
					ms.time = rand_date;
					nsamp++;

					int age = (rand_date/10000) - byear;

					//MLOG("pid %d age %d delta %d \n", rc.pid, age, delta);
					if (age >= s_params.min_age && age <= s_params.max_age)
						mis.samples.push_back(ms);

					to_days -= s_params.jump_days;
					delta -= s_params.jump_days;

				}

				if (mis.samples.size() > 0)
					samples.idSamples.push_back(mis);

			}


		}


	}

	samples.sort_by_id_date();
	if (samples.write_to_file(out_sample_file) < 0) {
		MERR("FAILED writing samples file %s\n", out_sample_file.c_str());
		return -1;
	}
	MLOG("Created samples file %s : %d samples for %d ids\n", out_sample_file.c_str(), nsamp, samples.idSamples.size());

	return 0;
}
