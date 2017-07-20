#pragma once
#ifndef __MED_COHORT_H__
#define __MED_COHORT_H__

//===================================================================================
// Data Structures and helpers to deal with a cohort.
// A cohort is simply a list of:
// - pid
// - follow up time : from , to
// - outcome date
// - outcome
//
// Major functionalities needed are:
// (1) read/write from/to file
// (2) Sample and convert to sampling file
// (3) Create incidence file
// (4) Calculate life time risk
//

#include <vector>
#include <MedProcessTools/MedProcessTools/SerializableObject.h>
using namespace std;


//===================================================================================
struct CohortRec : SerializableObject {
	int pid = -1;
	int from = 0;
	int to = 0;
	int outcome_date = 0;
	float outcome = -1;

	CohortRec() {};
	CohortRec(int _pid, int _from, int _to, int _outcome_date, float _outcome) {
		pid = _pid; from = _from; to = _to; outcome_date = _outcome_date; outcome = _outcome;
	}

	int init(map<string, string>& map);
	int get_string(string &to_str);
	int from_string(string &from_str);
};

//===================================================================================
struct SamplingParams : SerializableObject {
	float min_control_years = 0;
	float max_control_years = 10;
	float min_case_years = 0;
	float max_case_years = 1;
	int is_continous = 1;			// continous mode of sampling vs. stick to
	int min_days_from_outcome = 30;
	int jump_days = 180;
	int min_year = 1900;
	int max_year = 2100;
	int gender_mask = 0x3;
	int train_mask = 0x7;
	int min_age = 0;
	int max_age = 200;
	string rep_fname;

	// sticking related
	vector<string> stick_to_sigs;		// only use time points with these signals tested
	int take_closest = 0;
	int take_all = 0;

	int init(map<string, string>& map);
};

//===================================================================================
struct IncidenceParams : SerializableObject {
	int from_year = 2007;
	int to_year = 2013;
	int from_age = 30;
	int to_age = 90;
	int age_bin = 5;
	int min_samples_in_bin = 20;
	int gender_mask = 0x3;
	int train_mask = 0x7;
	string rep_fname;

	int init(map<string, string>& map);
};

//===================================================================================
class MedCohort : SerializableObject {

 public:

	vector<CohortRec> recs;

	void insert(int pid, int from, int to, int outcome_date, float outcome) { recs.push_back(CohortRec(pid, from, to, outcome_date, outcome)); }
	int read_from_file(string fname);
	int write_to_file(string fname);
	int read_from_bin_file(string fname) { return SerializableObject::read_from_file(fname); }
	int write_to_bin_file(string fname) { return SerializableObject::write_to_file(fname); }

	int get_pids(vector<int> &pids);

	//int print_general_stats();
	int create_incidence_file(IncidenceParams &i_params, string out_file);
	int create_sampling_file(SamplingParams &s_params, string out_sample_file);
	int create_sampling_file_sticked(SamplingParams &s_params, string out_sample_file);

};


//===================================================================================
// A few more MedSamples Helpers

// Scanner ::
// Given a MedSamples file , allows defining a sub-sample of it (can be all), 
// And define a list of tests and a list of base_tests
//
// The Scanner then allows the following:
// (1) Count for every test / base_test how many had at least N tests in a window W.
//     This helps in considering only variables that HAVE data.
// (2) Train a model M for each of:
//     - base_tests only
//     - base_tests + single test 
//	     -> for all train/test group
//       -> only for the subgroup of points that HAS no missing values (and compare to the base just on those)
//


#endif
