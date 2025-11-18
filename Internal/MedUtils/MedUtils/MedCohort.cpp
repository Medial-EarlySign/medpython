#include <MedUtils/MedUtils/MedCohort.h>

#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string.hpp>
#include <MedUtils/MedUtils/MedGenUtils.h>
#include <MedUtils/MedUtils/MedGlobalRNG.h>
#include <fstream>
#include <algorithm>
#include <fstream>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

// small helper for sampling
struct SamplesBucket {
	int from_date;
	int to_date;
	int from_days;
	int to_days;
	vector<int> sample_dates;
	vector<int> sample_days;
	void clear() { sample_dates.clear(); sample_days.clear(); }

};


//=====================================================================================
// CohortRec
//=====================================================================================
// Init from map
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

// Represent a cohort as a tab-delimited string
//-------------------------------------------------------------------------------------
void CohortRec::get_string(string &to_str)
{
	to_str = to_string(pid) + "\t" + to_string(from) + "\t" + to_string(to) + "\t" + to_string(outcome_date) + "\t" + to_string(outcome);
	if (comments != "") to_str += "\t" + comments;
}

// Get a cohort from a tab-delimited string. Return -1 if wrong number of fields, 0 upon success
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
	else if (fields.size() >= 6) {
		comments = fields[5];
	}
	else
		return -1;

	return 0;
}

//=====================================================================================
// SamplingParams
//=====================================================================================
// Init from map
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
		else if (m.first == "max_samples_per_id") max_samples_per_id = stoi(m.second);
		else if (m.first == "max_samples_per_id_method") max_samples_per_id_method = m.second;
		else if (m.first == "take_closest") take_closest = stoi(m.second);
		else if (m.first == "take_all") take_all = stoi(m.second);
		else if (m.first == "stick_to" || m.first == "stick_to_sigs") {
			boost::split(stick_to_sigs, m.second, boost::is_any_of(","));
			is_continous = 0;
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
// Init from map
//-------------------------------------------------------------------------------------
int IncidenceParams::init(map<string, string>& map)
{
	for (auto &m : map) {
		if (m.first == "age_bin") age_bin = stoi(m.second);
		else if (m.first == "min_samples_in_bin") min_samples_in_bin = stoi(m.second);
		else if (m.first == "from_year") from_year = stoi(m.second);
		else if (m.first == "to_year") to_year = stoi(m.second);
		else if (m.first == "start_date") start_date = stoi(m.second);
		else if (m.first == "gender_mask") gender_mask = stoi(m.second);
		else if (m.first == "train_mask") train_mask = stoi(m.second);
		else if (m.first == "from_age") from_age = stoi(m.second);
		else if (m.first == "to_age") to_age = stoi(m.second);
		else if (m.first == "rep") rep_fname = m.second;
		else if (m.first == "incidence_years_window") incidence_years_window = stoi(m.second);
		else if (m.first == "incidence_days_win") incidence_days_win = stoi(m.second);
		else {
			MERR("Unknown variable %s in IncidenceParams\n", m.first.c_str());
		}

	}
	return 0;
}

//=====================================================================================
// MedCohort
//=====================================================================================
// Read from tab-delimeted file. Return -1 if fail to open file
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

// Write to tab-delimeted file. Return -1 if fail to open file
//-------------------------------------------------------------------------------------
int MedCohort::write_to_file(string fname)
{
	ofstream of(fname);

	MLOG("MedCohort: writing to %s\n", fname.c_str());
	if (!of) {
		MERR("MedCohort: can't open file %s for writing\n", fname.c_str());
		return -1;
	}

	for (auto &rc : recs) {
		string sout;
		rc.get_string(sout);
		of << sout << endl;
	}

	MLOG("wrote [%d] records in cohort file %s\n", recs.size(), fname.c_str());
	of.close();
	return 0;
}

// Get all pids
//-------------------------------------------------------------------------------------
void MedCohort::get_pids(vector<int> &pids)
{
	pids.clear();
	for (auto &cr : recs) pids.push_back(cr.pid);
}

// Generate an incidence file from cohort + incidence-params
// Check all patient-years within cohort that fit to IncidenceParams and count positive outcomes within i_params.incidence_years_window
// Outcome - incidence per age-bin - is written to file
// Return 0 upon success. -1 upon failre to read repository
//-------------------------------------------------------------------------------------
int MedCohort::create_incidence_file(IncidenceParams &i_params, string out_file, const string &debug_file)
{
	//string inc_params; // from_to_fname,pids_to_use_fname,from_year,to_year,min_age,max_age,bin_size,inc_file

	vector<int> train_to_take = { 0,0,0,0 };
	if (i_params.train_mask & 0x1) train_to_take[1] = 1;
	if (i_params.train_mask & 0x2) train_to_take[2] = 1;
	if (i_params.train_mask & 0x4) train_to_take[3] = 1;

	vector<int> pids;
	get_pids(pids);

	// read byears, gender and TRAIN
	MedRepository rep;
	if (rep.read_all(i_params.rep_fname, pids, { "BDATE", "GENDER", "TRAIN" }) < 0) {
		MERR("FAILED reading repository %s\n", i_params.rep_fname.c_str());
		return -1;
	}

	// actual sampling
	unordered_set<int> pids_set;
	for (auto pid : pids) pids_set.insert(pid);

	map<int, pair<int, int>> counts, male_counts, female_counts;

	// age bins init
	for (int i = 0; i < 200; i += i_params.age_bin) counts[i] = pair<int, int>(0, 0);
	for (int i = 0; i < 200; i += i_params.age_bin) male_counts[i] = pair<int, int>(0, 0);
	for (int i = 0; i < 200; i += i_params.age_bin) female_counts[i] = pair<int, int>(0, 0);


	int bdate_sid = rep.sigs.sid("BDATE");
	int gender_sid = rep.sigs.sid("GENDER");
	int train_sid = rep.sigs.sid("TRAIN");

	vector<int> all_cnts = { 0,0 };

	//
	// To Estimate the annual statistics for the given time window we do the following:
	// We try to estimate the incidence at each 1.1.YYYY and then to weight average over all different years.
	// To do that we look at a cohort record that has a from_date and to_date , where in outcome==0 (controls) the to_date
	// is the last KNOWN date to be 0, and in cases (outcome!=0) it is the first date of conversion from case to control.
	//
	// Therefore:
	// (1) If the 1.1.YYYY is contained in the 0 period we count it as a 0 sample.
	// (2) If the 1.1.YYYY is contained in the 0 period AND the outcomedate for 1 is IN YYYY (or after if measuring longer periods) we count it as 1.
	//

	ofstream fw_debug;
	if (!debug_file.empty())
		fw_debug.open(debug_file);

	if (i_params.start_date == 101) {

		// old code for compatability

		for (auto &crec : recs) {
			int fyear = crec.from / 10000;
			if ((crec.from % 10000) > 101) fyear++;
			int to_date = crec.to;
			if (crec.outcome != 0) to_date = crec.outcome_date;
			else {
				if (i_params.incidence_days_win < 0)
					to_date -= i_params.incidence_years_window * 10000;
				else {
					int days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, to_date);
					days -= i_params.incidence_days_win;
					to_date = med_time_converter.convert_times(MedTime::Days, MedTime::Date, days);
				}
			}
			int tyear = to_date / 10000;
			int bdate = medial::repository::get_value(rep, crec.pid, bdate_sid);
			int byear = int(bdate / 10000);
			int gender = medial::repository::get_value(rep, crec.pid, gender_sid);
			int train = medial::repository::get_value(rep, crec.pid, train_sid);

			if ((gender & i_params.gender_mask) && (train_to_take[train]))
				for (int year = fyear; year <= tyear; year++) {
					if (year >= i_params.from_year && year <= i_params.to_year) {

						int age = year - byear;
						int bin = i_params.age_bin*(age / i_params.age_bin);
						counts[bin].first++; all_cnts[0]++;
						if (gender == GENDER_MALE)
							++male_counts[bin].first;
						else if (gender == GENDER_FEMALE)
							++female_counts[bin].first;

						bool count_this_year = crec.outcome == 0;
						if (crec.outcome != 0) {
							// handlind cases
							// first we will calculate if 1.1.year is indeed at most incidence_years or incindence_days BEFORE the outcome date


							if (i_params.incidence_days_win < 0) {
								// case1 : we use years:
								count_this_year = (year > tyear - i_params.incidence_years_window);
							}
							else {
								// case2 : we use days (remember to_date now is the date of the case event)
								int to_days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, to_date);
								int curr_days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, year * 10000 + 101);
								count_this_year = (to_days - curr_days <= i_params.incidence_days_win);
							}

							if (count_this_year) {
								counts[bin].second++; all_cnts[1]++;
								if (gender == GENDER_MALE)
									++male_counts[bin].second;
								else if (gender == GENDER_FEMALE)
									++female_counts[bin].second;
							}


						}

						if (!debug_file.empty() && count_this_year && age >= i_params.from_age && age <= i_params.to_age) {
							//Debug: pid, year, outcome, age, gender
							fw_debug << crec.pid << "\t" << year << "\t" << crec.outcome << "\t" << age << "\t" << gender
								<< "\n";
						}
					}
				}
		}
	}
	else {

		// general purpose code : gets a start_date, and works with incidence_days
		int incidence_days = i_params.incidence_years_window * 365;
		if (i_params.incidence_days_win > 0)
			incidence_days = i_params.incidence_days_win;

		MLOG("INCIDENCE_DAYS IS %d (in params: days %d years %d)\n", incidence_days, i_params.incidence_days_win, i_params.incidence_years_window);

		for (auto &crec : recs) {

			// the eligible dates are the dates of the form
			// YYYYMMDD , where MMDD is our date and YYYYMMDD is in the range [crec.from, crec.to]
			// We calculate the eligible dates and the relevant Y (outcome) for each
			// Then we can simply go over them , filter some more conditions and sum

			vector<int> edates, Y;
			int to_date = crec.to;
			if (crec.outcome != 0) to_date = crec.outcome_date;
			else {
				int days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, to_date);
				days -= incidence_days;
				to_date = med_time_converter.convert_times(MedTime::Days, MedTime::Date, days);
			}

			int year = crec.from / 10000;
			int edate = 0;
			while (edate < to_date) {
				edate = year * 10000 + i_params.start_date;
				if (edate >= crec.from && edate < to_date) {
					edates.push_back(edate);
					if (crec.outcome == 0)
						Y.push_back(0);
					else {
						// have to check the window for incidence_days
						int to_days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, to_date);
						int curr_days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, edate);
						//MLOG("edate %d to %d to_days %d curr_days %d diff %d incidence_days %d\n", edate, crec.to, to_days, curr_days, to_days-curr_days, incidence_days);
						if (to_days - curr_days <= incidence_days)
							Y.push_back(1);
						else
							Y.push_back(0);
					}
				}
				year++;
			}

			//MLOG("Cohort: pid %d : %d - %d : %d , %f : to_date %d : edates :", crec.pid, crec.from, crec.to, crec.outcome_date, crec.outcome, to_date);
			//for (int i=0; i<edates.size(); i++) MLOG(" %d,%d", edates[i], Y[i]);
			//MLOG("\n");

			// go over dates, filter and count
			int bdate = medial::repository::get_value(rep, crec.pid, bdate_sid);
			int byear = int(bdate / 10000);
			int gender = medial::repository::get_value(rep, crec.pid, gender_sid);
			int train = medial::repository::get_value(rep, crec.pid, train_sid);
			for (int i = 0; i < edates.size(); i++) {
				int year = edates[i] / 10000;
				if ((gender & i_params.gender_mask) && (train_to_take[train]))
					if (year >= i_params.from_year && year <= i_params.to_year) {

						int age = year - byear;
						int bin = i_params.age_bin*(age / i_params.age_bin);

						// outcome 0 cases - always counted
						counts[bin].first++; all_cnts[0]++;
						if (gender == GENDER_MALE)
							++male_counts[bin].first;
						else if (gender == GENDER_FEMALE)
							++female_counts[bin].first;


						// outcome 1 cases
						if (Y[i]) {
							counts[bin].second++; all_cnts[1]++;
							if (gender == GENDER_MALE)
								++male_counts[bin].second;
							else if (gender == GENDER_FEMALE)
								++female_counts[bin].second;
						}
					}
			}

		}

	}

	if (!debug_file.empty())
		fw_debug.close();

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
	if (!of.good())
		MTHROW_AND_ERR("IO Error: can't write \"%s\"\n", out_file.c_str());

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

	//New Format:
	ofstream of_new(out_file + ".new_format");
	if (!of_new.good()) {
		MERR("IO Error: can't write \"%s\"\n", (out_file + ".new_format").c_str());
		return -1;
	}

	of_new << "AGE_BIN" << "\t" << i_params.age_bin << "\n";
	of_new << "AGE_MIN" << "\t" << i_params.from_age << "\n";
	of_new << "AGE_MAX" << "\t" << i_params.to_age << "\n";
	of_new << "OUTCOME_VALUE" << "\t" << "0.0" << "\n";
	of_new << "OUTCOME_VALUE" << "\t" << "1.0" << "\n";

	for (auto it = counts.begin(); it != counts.end(); ++it) {

		int age = it->first;
		int male_n0 = male_counts[it->first].first;
		int male_n1 = male_counts[it->first].second;

		if (age >= i_params.from_age && age <= i_params.to_age && male_n0 > 0) {
			of_new << "STATS_ROW" << "\t" << "MALE" << "\t" <<
				age + i_params.age_bin / 2 << "\t" << "0.0" << "\t" << male_n0 - male_n1 << "\n";
			of_new << "STATS_ROW" << "\t" << "MALE" << "\t" <<
				age + i_params.age_bin / 2 << "\t" << "1.0" << "\t" << male_n1 << "\n";
		}

		int female_n0 = female_counts[it->first].first;
		int female_n1 = female_counts[it->first].second;

		if (age >= i_params.from_age && age <= i_params.to_age && female_n0 > 0) {
			of_new << "STATS_ROW" << "\t" << "FEMALE" << "\t" <<
				age + i_params.age_bin / 2 << "\t" << "0.0" << "\t" << female_n0 - female_n1 << "\n";
			of_new << "STATS_ROW" << "\t" << "FEMALE" << "\t" <<
				age + i_params.age_bin / 2 << "\t" << "1.0" << "\t" << female_n1 << "\n";
		}
	}

	of_new.close();

	return 0;
}

// Generate a samples file from cohort + sampling-params
// Generate samples within cohort times that fit SampleingParams criteria and windows.
// Sample dates are selected randomly for each window of s_params.jump_days in the legal period, and written to file
// Return 0 upon success. -1 upon failre to read repository
//-------------------------------------------------------------------------------------
int MedCohort::create_sampling_file(SamplingParams &s_params, string out_sample_file)
{

	if (s_params.is_continous == 0)
		return create_sampling_file_sticked(s_params, out_sample_file);

	vector<int> pids;
	get_pids(pids);
	MedRepository rep;
	if (rep.read_all(s_params.rep_fname, pids, { "BDATE", "GENDER", "TRAIN" }) < 0) {
		MERR("FAILED reading repository %s\n", s_params.rep_fname.c_str());
		return -1;
	}

	MedSamples samples;
	create_samples(rep, s_params, samples);

	if (samples.write_to_file(out_sample_file) < 0) {
		MERR("FAILED writing samples file %s\n", out_sample_file.c_str());
		return -1;
	}
	MLOG("Created samples file %s : %d samples for %d ids\n", out_sample_file.c_str(), samples.nSamples(), samples.idSamples.size());
	return 0;
}

int MedCohort::create_samples(MedRepository& rep, SamplingParams &s_params, MedSamples& samples)
{

	vector<int> train_to_take = { 0,0,0,0 };
	if (s_params.train_mask & 0x1) train_to_take[1] = 1;
	if (s_params.train_mask & 0x2) train_to_take[2] = 1;
	if (s_params.train_mask & 0x4) train_to_take[3] = 1;

	int bdate_sid = rep.sigs.sid("BDATE");
	int gender_sid = rep.sigs.sid("GENDER");
	int train_sid = rep.sigs.sid("TRAIN");

	int nsamp = 0;

	int min_date = s_params.min_year * 10000 + 0101;
	int max_date = s_params.max_year * 10000 + 1230;

	for (auto &rc : recs) {

		if (rc.from < min_date) rc.from = min_date;
		if (rc.to > max_date) rc.from = max_date;

		int bdate = medial::repository::get_value(rep, rc.pid, bdate_sid);
		int byear = int(bdate / 10000);
		int gender = medial::repository::get_value(rep, rc.pid, gender_sid);
		int train = medial::repository::get_value(rep, rc.pid, train_sid);

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
					from_days = max(from_days, outcome_days - int(365.0f * s_params.max_case_years));
				}
				else {
					// control
					to_days = outcome_days - int(365.0f * s_params.min_control_years);
					from_days = max(from_days, outcome_days - int(365.0f * s_params.max_control_years));
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
					mis.samples.push_back(ms);

					to_days -= s_params.jump_days;
					delta -= s_params.jump_days;

				}
				vector<MedSample>::iterator last_it = mis.samples.end();
				if (mis.samples.size() > s_params.max_samples_per_id) {
					// randomizing
					if (s_params.max_samples_per_id_method == "rand") {
						shuffle(mis.samples.begin(), mis.samples.end(), globalRNG::get_engine());
					}
					// It is assumed that samples are ordered from last to first, so other s_params.max_samples_per_id_method = 'last' should do nothing before setting end iterator.
					last_it = mis.samples.begin() + s_params.max_samples_per_id;
				}

				MedIdSamples mis_new;
				for (auto ms = mis.samples.begin(); ms != last_it; ++ms) {
					int age = (ms->time / 10000) - byear;

					//MLOG("pid %d age %d delta %d \n", rc.pid, age, delta);
					if (age >= s_params.min_age && age <= s_params.max_age)
						mis_new.samples.push_back(*ms);
				}

				mis = mis_new;
				mis.id = rc.pid;

				if (mis.samples.size() > 0)
					samples.idSamples.push_back(mis);
			}
		}
	}

	samples.normalize();


	return 0;
}

// Generate a samples file from cohort + sampling-params
// Generate samples within cohort times that fit SampleingParams criteria and windows.
// Sample dates are those with the required signals for each window of s_params.jump_days in the legal period (if existing), and written to file
//-------------------------------------------------------------------------------------
int MedCohort::create_sampling_file_sticked(SamplingParams &s_params, string out_sample_file)
{

	vector<int> pids;
	get_pids(pids);
	MedRepository rep;
	vector<string> sigs = { "BDATE", "GENDER", "TRAIN" };
	sigs.insert(sigs.end(), s_params.stick_to_sigs.begin(), s_params.stick_to_sigs.end());
	if (rep.read_all(s_params.rep_fname, pids, sigs) < 0) {
		MERR("FAILED reading repository %s\n", s_params.rep_fname.c_str());
		return -1;
	}

	MedSamples samples;
	create_samples_sticked(rep, s_params, samples);

	if (samples.write_to_file(out_sample_file) < 0) {
		MERR("FAILED writing samples file %s\n", out_sample_file.c_str());
		return -1;
	}
	MLOG("Created samples file %s : %d samples for %d ids\n", out_sample_file.c_str(), samples.nSamples(), samples.idSamples.size());
	return 0;
}
int MedCohort::create_samples_sticked(MedRepository& rep, SamplingParams &s_params, MedSamples& samples)
{
	vector<int> train_to_take = { 0,0,0,0 };
	if (s_params.train_mask & 0x1) train_to_take[1] = 1;
	if (s_params.train_mask & 0x2) train_to_take[2] = 1;
	if (s_params.train_mask & 0x4) train_to_take[3] = 1;


	int bdate_sid = rep.sigs.sid("BDATE");
	int gender_sid = rep.sigs.sid("GENDER");
	int train_sid = rep.sigs.sid("TRAIN");
	//vector<int> sids_to_stick;
	//for (auto &sig : sigs) sids_to_stick.push_back(rep.sigs.sid(sig));

	vector<int> dates_to_take;

	int from0_days = (int)(s_params.max_control_years * 365.0f);
	int to0_days = (int)(s_params.min_control_years * 365.0f);
	int from1_days = (int)(s_params.max_case_years * 365.0f);
	int to1_days = (int)(s_params.min_case_years * 365.0f);
	int min_date = s_params.min_year * 10000 + 0101;
	int max_date = s_params.max_year * 10000 + 1230;

	int nsamp = 0;

	for (auto &rc : recs) {
		bool print = false;
		int bdate = medial::repository::get_value(rep, rc.pid, bdate_sid);
		int byear = int(bdate / 10000);
		int gender = medial::repository::get_value(rep, rc.pid, gender_sid);
		int train = medial::repository::get_value(rep, rc.pid, train_sid);

		if (print) MLOG("pid %d from %d to %d outcome %f outcome %d : byear %d gender %d train %d\n", rc.pid, rc.from, rc.to, rc.outcome, rc.outcome_date, byear, gender, train);
		if (print) MLOG("pid %d from0 %d to0 %d from1 %d to1 %d\n", rc.pid, from0_days, to0_days, from1_days, to1_days);

		if ((gender & s_params.gender_mask) && (train_to_take[train])) {

			vector<int> dates_with_sigs;
			rep.get_dates_with_signal(rc.pid, s_params.stick_to_sigs, dates_with_sigs);

			int outcome_days = med_time_converter.convert_date(MedTime::Days, rc.outcome_date);

			if (print) MLOG("pid %d outcome_days %d dates_with_sig size: %d \n", rc.pid, outcome_days, dates_with_sigs.size());
			// split dates_with_sigs into buckets
			map<int, vector<int>> buckets;
			for (auto date : dates_with_sigs) {
				if (date < min_date || date > max_date)
					continue;
				int days = med_time_converter.convert_date(MedTime::Days, date);
				int relative_days = outcome_days - days;
				int age = date / 10000 - byear;
				bool date_is_in_cohort = (date >= rc.from && date <= rc.to);
				bool date_is_in_sample_window = (((rc.outcome == 0) && (relative_days <= from0_days) && (relative_days > to0_days)) ||
					((rc.outcome != 0) && (relative_days <= from1_days) && (relative_days > to1_days)));
				bool date_is_in_age_range = ((age >= s_params.min_age) && (age <= s_params.max_age));
				int bucket_num = 0;
				if (relative_days > 0)
					bucket_num = relative_days / s_params.jump_days;
				else
					bucket_num = -((-relative_days) / s_params.jump_days) - 1;
				if (print) MLOG("pid %d date %d days %d relative_days %d age %d in_cohort %d in_win %d in age %d bucket %d\n",
					rc.pid, date, days, relative_days, age, date_is_in_cohort, date_is_in_sample_window, date_is_in_age_range, bucket_num);

				if (date_is_in_cohort && date_is_in_sample_window && date_is_in_age_range) {
					// Important: only samples that can be considered are put into a bucket.
					if (buckets.find(bucket_num) == buckets.end()) {
						buckets[bucket_num] = vector<int>();
					}
					buckets[bucket_num].push_back(date);
				}
			}

			// Now going over buckets and actually sampling
			dates_to_take = {};
			for (auto &bucket : buckets) {
				if (s_params.take_all)
					dates_to_take.insert(dates_to_take.end(), bucket.second.begin(), bucket.second.end());
				else if (s_params.take_closest)
					dates_to_take.push_back(bucket.second.back());
				else {
					int r = rand_N((int)bucket.second.size());
					dates_to_take.push_back(bucket.second[r]);
				}
			}

			if (dates_to_take.size() > s_params.max_samples_per_id) {
				if (s_params.max_samples_per_id_method == "last") {
					sort(dates_to_take.begin(), dates_to_take.end(), std::greater<int>());
				}
				if (s_params.max_samples_per_id_method == "rand") {
					shuffle(dates_to_take.begin(), dates_to_take.end(), globalRNG::get_engine());
				}
			}

			// Now simply pushing into our MedSamples
			MedIdSamples mis;
			mis.id = rc.pid;

			int samples_cnt = 0;
			for (auto date : dates_to_take) {
				if (samples_cnt >= s_params.max_samples_per_id) { break; }
				MedSample ms;
				ms.id = rc.pid;
				ms.outcome = rc.outcome;
				ms.outcomeTime = rc.outcome_date;
				ms.time = date;
				nsamp++;
				samples_cnt++;
				mis.samples.push_back(ms);
			}
			samples.idSamples.push_back(mis);

		}
	}

	samples.sort_by_id_date();

	return 0;
}
