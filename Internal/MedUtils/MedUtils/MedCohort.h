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
// CohortRec : a single entry within a cohort; includes:
//	- pid
// - follow up time : from , to
// - outcome date
// - outcome
//
//===================================================================================
struct CohortRec : SerializableObject {
	int pid = -1;			// Patient Id
	int from = 0;			// Followup start
	int to = 0;				// Followup end
	int outcome_date = 0;	// Date(Time) at which outcome is given
	float outcome = -1;		// Outcome

	// Constructors Initialization
	CohortRec() {};
	CohortRec(int _pid, int _from, int _to, int _outcome_date, float _outcome) {
		pid = _pid; from = _from; to = _to; outcome_date = _outcome_date; outcome = _outcome;
	}
	int init(map<string, string>& map);

	// Transfer to/from string
	int get_string(string &to_str);
	int from_string(string &from_str);
};

//===================================================================================
// SamplingParams : Parameters for sampling from repostory + cohort
//===================================================================================
struct SamplingParams : SerializableObject {
	float min_control_years = 0;
	float max_control_years = 10;
	float min_case_years = 0;
	float max_case_years = 1;
	int is_continous = 1;			// continous mode of sampling vs. stick to (0 = stick)
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
// IncidenceParams: Parameters for calculating incidence from repostory + cohort
//===================================================================================
struct IncidenceParams : SerializableObject {
	int from_year = 2007;
	int to_year = 2013;
	int from_age = 30;
	int to_age = 90;
	int age_bin = 5;
	int incidence_years_window = 1; // how many years ahead do we consider an outcome?
	int min_samples_in_bin = 20;
	int gender_mask = 0x3;
	int train_mask = 0x7;
	string rep_fname;

	int init(map<string, string>& map);
};

//===================================================================================
// Cohort : a vector of CohortRecs - eac
//===================================================================================
class MedCohort : SerializableObject {

 public:

	vector<CohortRec> recs;

	// Add a record
	void insert(int pid, int from, int to, int outcome_date, float outcome) { recs.push_back(CohortRec(pid, from, to, outcome_date, outcome)); }
	// Read/Write to/from files
	int read_from_file(string fname);
	int write_to_file(string fname);
	int read_from_bin_file(string fname) { return SerializableObject::read_from_file(fname); }
	int write_to_bin_file(string fname) { return SerializableObject::write_to_file(fname); }

	// Get all pids
	int get_pids(vector<int> &pids);

	//int print_general_stats();
	// Generate an incidence file from cohort + incidence-params
	// Check all patient-years within cohort that fit to IncidenceParams and count positive outcomes within i_params.incidence_years_window
	// Outcome - incidence per age-bin - is written to file
	int create_incidence_file(IncidenceParams &i_params, string out_file);

	// Generate a samples file from cohort + sampling-params
	// Generate samples within cohort times that fit SampleingParams criteria and windows.
	// Sample dates are selected randomly for each window of s_params.jump_days in the legal period, and written to file
	int create_sampling_file(SamplingParams &s_params, string out_sample_file);

	// Generate a samples file from cohort + sampling-params
	// Generate samples within cohort times that fit SampleingParams criteria and windows.
	// Sample dates are those with the required signals for each window of s_params.jump_days in the legal period (if existing), and written to file
	//-------------------------------------------------------------------------------------
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
