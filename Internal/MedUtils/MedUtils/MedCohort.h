#pragma once
#ifndef __MED_COHORT_H__
#define __MED_COHORT_H__

//===================================================================================
/** Data Structures and helpers to deal with a cohort. <br>
* A cohort is simply a list of: <br>
* - pid <br>
* - follow up time : from , to <br>
* - outcome date <br>
* - outcome <br>
*
* Major functionalities needed are: <br>
* (1) read/write from/to file <br>
* (2) Sample and convert to sampling file <br>
* (3) Create incidence file <br>
* (4) Calculate life time risk <br>
*/
//===================================================================================


#include <vector>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <InfraMed/InfraMed/InfraMed.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>

using namespace std;


//===================================================================================
/** CohortRec : a single entry within a cohort; includes: <br>
*	- pid <br>
* - follow up time : from , to <br>
* - outcome date <br>
* - outcome <br>
*/
//===================================================================================
struct CohortRec : SerializableObject {
	int pid = -1;			///< Patient Id
	int from = 0;			///< Followup start
	int to = 0;				///< Followup end
	int outcome_date = 0;	///< Date(Time) at which outcome is given
	float outcome = -1;		///< Outcome
	string comments = "";	///< additional option for comments

	/// <summary> empty constructor </summary>
	CohortRec() {};
	/// <summary> Constructor with initialization </summary>
	CohortRec(int _pid, int _from, int _to, int _outcome_date, float _outcome) {
		pid = _pid; from = _from; to = _to; outcome_date = _outcome_date; outcome = _outcome;
	}

	/// <summary> Initialize from a map </summary>
	int init(map<string, string>& map);

	/// <summary> Represent a cohort as a tab-delimited string </summary>
	void get_string(string &to_str);
	/// <summary> Get a cohort rec from a tab-delimited string </summary>
	/// <returns> -1 if wrong number of fields, 0 upon success </returns>
	int from_string(string &from_str);
};

//===================================================================================
/** SamplingParams : Parameters for sampling from repostory + cohort
*/
//===================================================================================
struct SamplingParams : SerializableObject {
	float min_control_years = 0;	///< minimal number of years before outcome for controls
	float max_control_years = 10;	///< maximal number of years before outcome for controls
	float min_case_years = 0;	///< minimal number of years before outcome for cases
	float max_case_years = 1;	///< maximal number of years before outcome for cases
	int is_continous = 1;			///< continous mode of sampling vs. stick to (0 = stick)
	int min_days_from_outcome = 30;		///< minimal number of days before outcome
	int jump_days = 180;	///< days to jump between sampling periods
	int min_year = 1900;	///< first year for sampling
	int max_year = 2100;	///< last year for sampling
	int gender_mask = 0x3;	///< mask for gender specification (rightmost bit on for male, second for female)
	int train_mask = 0x7;	///< mask for TRAIN-value specification (three rightmost bits for TRAIN = 1,2,3)
	int min_age = 0;	///< minimum age for sampling
	int max_age = 200;	///< maximum age for sampling
	string rep_fname;	///< Repository configration file

	/// sticking related. if none of take_closest/take_all is on, a random sample with requrired-signal within each sampling period is selected
	vector<string> stick_to_sigs;		///< only use time points with these signals tested
	int take_closest = 0;	///< flag: take the sample with requrired-signals that is closest to each target sampling-date
	int take_all = 0; ///< flag: take all samples with requrired-signal within each sampling period is selected
	int max_samples_per_id = (int)1e9; ///< maximum samples per ID
	string max_samples_per_id_method = "last"; ///< determine how to pick samples - 'last' or 'rand'

	/// <summary> Initialize from a map </summary>
	int init(map<string, string>& map);
};

//===================================================================================
/** IncidenceParams: Parameters for calculating incidence from repostory + cohort
*/
//===================================================================================
struct IncidenceParams : SerializableObject {
	int from_year = 2007;	///< first year to consider in calculating incidence
	int to_year = 2013;	///< last year to consider in calculating incidence
	int start_date = 101; /// the date in each year to start looking from (default is 0101), format is MMDD
	int from_age = 30;	///< minimal age to consider
	int to_age = 90;	///< maximal age to consider
	int age_bin = 5;	///< binning of ages
	int incidence_years_window = 1; ///< how many years ahead do we consider an outcome?
	int incidence_days_win = -1;  /// if -1: using incidence_years_window
	int min_samples_in_bin = 20; ///< minimal required samples to estimate incidence per bin
	int gender_mask = 0x3;	///< mask for gender specification (rightmost bit on for male, second for female)
	int train_mask = 0x7;	///< mask for TRAIN-value specification (three rightmost bits for TRAIN = 1,2,3)
	string rep_fname;	///< Repository configration file


	/// <summary> Initialize from a map </summary>
	int init(map<string, string>& map);
};

//===================================================================================
/** MedCohort - a vector of CohortRec's
*/
//===================================================================================
class MedCohort : SerializableObject {

 public:
	vector<CohortRec> recs; ///< Cohort information

	/// <summary>
	/// Add a record
	/// </summary>
	void insert(int pid, int from, int to, int outcome_date, float outcome) { recs.push_back(CohortRec(pid, from, to, outcome_date, outcome)); }

	/// <summary> Read to tab-delimited file </summary>
	/// <returns> 1- if fail to open, 0 upon success </returns>
	int read_from_file(string fname);
	/// <summary>  Write from tab-delimited file </summary>
	/// <returns> 1- if fail to open, 0 upon success </returns>
	int write_to_file(string fname);

	/// <summary> Read to binary file </summary>
	/// <returns> 1- if fail to open, 0 upon success </returns>
	int read_from_bin_file(string fname) { return SerializableObject::read_from_file(fname); }
	/// <summary> Write from binary file </summary>
	/// <returns>1- if fail to open, 0 upon success </returns>
	int write_to_bin_file(string fname) { return SerializableObject::write_to_file(fname); }

	/// <summary> Get all pids </summary>
	void get_pids(vector<int> &pids);

	//int print_general_stats();

	/// <summary>
	/// Generate an incidence file from cohort + incidence-params <br>
	/// Check all patient-years within cohort that fit to IncidenceParams and count positive outcomes within i_params.incidence_years_window <br>
	/// Outcome - incidence per age-bin - is written to file 
	/// </summary>
	/// <returns> -1 if writing to file failed, 0 upon success </returns>
	int create_incidence_file(IncidenceParams &i_params, string out_file, const string &debug_file = "");

	/// <summary>
	/// Generate a samples file from cohort + sampling-params <br>
	/// Generate samples within cohort times that fit SampleingParams criteria and windows. <br>
	/// Sample dates are selected randomly for each window of s_params.jump_days in the legal period, and written to file <br>
	/// </summary>
	/// <returns> 0 upon success. -1 upon failre to read repository </returns>
	int create_sampling_file(SamplingParams &s_params, string out_sample_file);
	int create_samples(MedRepository& rep, SamplingParams &s_params, MedSamples& samples);

	/// <summary>
	/// Generate a samples file from cohort + sampling-params <br>
	/// Generate samples within cohort times that fit SampleingParams criteria and windows. <br>
	/// Sample dates are those with the required signals for each window of s_params.jump_days in the legal period (if existing), and written to file <br>
	/// </summary>
	/// <returns>  0 upon success. -1 upon failre to read repository </returns>
	int create_sampling_file_sticked(SamplingParams &s_params, string out_sample_file);
	int create_samples_sticked(MedRepository& rep, SamplingParams &s_params, MedSamples& samples);
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
