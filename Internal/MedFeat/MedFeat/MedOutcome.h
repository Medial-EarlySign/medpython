//
// MedOutcome - classes and tools to help managing medical outcomes ("Y").
//

#ifndef __MEDOUTCOME_H__
#define __MEDOUTCOME_H__

//#include "MedSplit.h"
#include "MedFeat.h"
#include <InfraMed/InfraMed/InfraMed.h>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#define MAX_OUTCOME_LENGTH		1000000
using namespace std;

// OLD CODE NOT NEEDED

#if 0
class OutcomeRecord {
	public:
		int pid;
		int date; // 0 date for an outcome means the category for the relevant id could be at any time point
		unsigned long long time;
		unsigned long long raw_time;
		int raw_date;
		float categ;
		int length;
		OutcomeRecord() { pid=0; date=0; time=0; raw_time=0; raw_date=0; categ=0; length=0; }
		//inline bool operator< (const OutcomeRecord &or2) const { if (pid == or2.pid) return (time < or2.time); return (pid < or2.pid); };
		inline bool operator< (const OutcomeRecord &or2) const { if (pid == or2.pid) return (date < or2.date); return (pid < or2.pid); };
		inline bool operator== (const OutcomeRecord &or2) const { if (pid == or2.pid && date == or2.date && time == or2.time && raw_time == or2.raw_time && raw_date == or2.raw_date && categ == or2.categ &&  length == or2.length) return true; return false; };
};


class MedOutcome {
	public:

		// core vals
		string name;
		string desc;
		string type; // single, binary, multi, regression
		vector<OutcomeRecord> out;
		vector<int> ignore_pids;
		int n_categ; // 0 for regression

		// calcuated
		map<int,int> pid2use; // a map from pid->int containing -1 for ignored pids, and 1 for event pids (control pids are considered to be all the rest, but can be inside too (as 0))

		void clear() {name=""; desc=""; type=""; out.clear(); ignore_pids.clear(); pid2use.clear();}

		// file io
		int read_from_file(const string &fname);
		int write_to_file(const string &fname);

		void shuffle() {random_shuffle(out.begin(),out.end());} // shuffles out vector randomly
		void sort() { std::sort(out.begin(), out.end()); }
		void sort_unique() { MedOutcome::sort(); auto it = std::unique(out.begin(), out.end()); out.resize(std::distance(out.begin(), it)); }

		void subsample(float prob);

};


// given an outcome we often need to filter it and throw some ids due to some basic rules (dates, age, gender, etc)
class OutcomeFilter {
	public:
		// filtering related
		int filter_flag;			// if 0 toggles off filtering (helps in some APIs)
		int min_date,max_date;		// min_date < 0 means no filtering on dates
		unsigned int train_to_take; // TRAIN values are 1,2,3 (train,internal,external) , and this is a bit_mask to tell which to take (0x7 is all, 0x4 is just external, etc...)
		int gender;					// 0x3 both, 0x1 men, 0x2 women
		float min_age, max_age;
		int filter_just_demographics; // filtering only date, train, gender, age

		int history_time_before;	// how much days of history has the sample got (at least)
		vector<string> clock_sigs;	// signals by which to measure the history_time_before attribute

		vector<string> sigs;
		int sigs_win_min, sigs_win_max; // negative values are AFTER the event....
		int min_sigs_time_points_in_win;

		// moving to closest point
		int stick_to_sigs;			// will move date to closest sigs time point in window (in case window restrictions are on)
		int n_points_to_take_in_win;

		// matching related
		int match_flag;
		int match_age_flag;
		float age_bin_width; // in years
		int match_do_shuffle_flag;
		int match_const_flag;
		int match_gender_flag;
		int match_dates_flag;
		float dates_bin_width; // in days
		int match_location_flag;
		int match_signals_flag; // matching to the signals values of the closest before test
		vector<string> sigs_to_match;
		vector<float> sig_bin_width;
		int match_min_size; // min number of events in a group to consider for ratio calculations
		int match_max_ratio;
		float match_event_case_price_ratio;	// when calculating the optimal match ratio, how many controls are we willing to lose to save one case?
		int match_max_samples;

		OutcomeFilter() {filter_flag = 0; min_date = -1; train_to_take=0x7; gender=0x3; min_age=0; max_age=185; filter_just_demographics=0; sigs.clear(); 
					     sigs_win_min=0; sigs_win_max=360; min_sigs_time_points_in_win=1; history_time_before = 0; clock_sigs.clear();
						 stick_to_sigs = 0; n_points_to_take_in_win = 1;
						 match_flag = 0; match_age_flag = 0; age_bin_width = 5; match_gender_flag = 0; match_const_flag = 0; match_do_shuffle_flag = 1; match_dates_flag = 0; dates_bin_width = 365;
						 match_location_flag = 0; match_signals_flag = 0; sigs_to_match.clear();
						 match_min_size = 10; match_max_ratio = 4; match_event_case_price_ratio = 100.0; match_max_samples = 1000000;
		}

		// the following is very useful when in need of init from command_line or file (sits on top of defaults)
		void init_from_string(const string &s);

};

// filtering a single case
int filter_outcome_record(MedRepository &rep, OutcomeFilter &filter, OutcomeRecord &orec);

// filtering a whole outcome
void get_filtered_outcome(MedRepository &rep, MedOutcome &mout, OutcomeFilter &filter, MedOutcome &mout_f);

// adding randomized control points to a single outcome
void simple_single_to_binary_outcome(MedRepository &rep, MedOutcome &mout_s, MedOutcome &mout_b);

// adding matched control points to a single outcome
void matched_single_to_binary_outcome(MedRepository &rep, OutcomeFilter filter, MedOutcome &mout_s, MedOutcome &mout_b);

void match_binary_outcome(MedRepository &rep, const string &match_params, MedOutcome &mout_in, MedOutcome &mout_out);
void match_binary_outcome(MedRepository &rep, OutcomeFilter &filter, MedOutcome &mout_in, MedOutcome &mout_out);

// general way to create features given a repository, an outcome, a list of raw signals, and feature types
// (a work in progress)
//int collect_features_for_outcome(MedRepository &rep, MedOutcome &mout, vector<string> &raw_sigs, vector<string> &feat_types, float missing, MedFeaturesData &mf);
int collect_features_for_outcome(MedRepository &rep, MedOutcome &mout, vector<string> &sigs, float missing, MedFeaturesData &mf);

// for each entry in mout_in, adds ALL the samples in the range of 0-days (with the same categ)
// this is done to enrich samples to contain samples in all periods before the event
void expand_binary_outcome_prefix(MedRepository &rep,int days, const string &sig, MedOutcome &mout_in, MedOutcome &mout_out);



// a few more classes needed for easier management

// prepare a balanced split given a repository, and an outcome
// the result is a split to the pids in the outcome, such that in each split there's the same distribution of categories
int get_balanced_split(MedSplit &spl, int nsplits, MedRepository &rep, MedOutcome &mout);

// mode - "train" : get all pids not in split , "test" : get all pids in split
int get_split_from_outcome(MedOutcome &mout, MedSplit &spl, int isplit, const string &mode, MedOutcome &mout_s);

// next class holds the parameters needed in order to create a feature
enum FeatureType {F_NONE, F_Age_Calculated, F_Age_Param, F_Const, F_Last, F_TLast, F_Avg, F_Max, F_Min, F_Delta, F_Slope , F_NumV };

class FeatParams {
	public:
		string full_name;
		string sig_name;
		string feat_name;
		int sid;
		int sig_type;
		FeatureType f_type;
		int ind;
		int from_day;
		int to_day;
		unsigned long long from_time;
		unsigned long long to_time;

		FeatParams() { 
			full_name=""; sig_name=""; feat_name="";  sid = -1; sig_type = -1; f_type = F_NONE; ind = 1; from_day = 0; to_day = 365000; 
			from_time = 0; to_time = ((unsigned long long)1)<<45; 
		}

		int init(const string &in_name);
		int init_sid(MedRepository &rep);
};

class MedFeatList {
	private:
		void push_feature(string feat);
	public:
		string desc;
		vector<string> features;
		vector<FeatParams> feats;
		vector<string> signals; // keeping also a list of raw signals, to enable a faster read of a repository - only on the asked for signals.

		MedFeatList() {clear();}
		void clear() { features.clear(); signals.clear(); feats.clear(); }

		int add_sid_to_feats(MedRepository &rep) { for (auto &fp : feats) fp.init_sid(rep); }


		int read_from_file(const string &fname);
		int write_to_file(const string &fname);

};


//-----------------------------------------------------------------------
// following: very simple ways to keep prediction results

struct PredictedRes {
	int pid;
	long long date;				// can hold time too
	float pred;
	long long outcome_date;		// the date of the outcome allows us to calculate the time between the prediction and the event
	float outcome;			// actual outcome we wanted to predict
	int	split;
	int nsplits;			// is split == nsplits , we know we have a true validation set (out of train and cv)

	void format(string &s);	  // formating a line ready to be printed
	int parse(string &line);  // parsing a line containing a record

	inline bool operator< (const PredictedRes &pr2) const {
		if (split == pr2.split) {
			if (pid == pr2.pid) return (date < pr2.date);
			return (pid < pr2.pid);
		}
		return (split < pr2.split);
	};
};
inline bool comp_pred(const PredictedRes &pr1, const PredictedRes &pr2) {
	return pr1.pred < pr2.pred;
}
class PredictionList {
	public:
		vector<PredictedRes> preds;

		void sort() { std::sort(preds.begin(), preds.end()); }
		void sort_by_pred() { std::sort(preds.begin(), preds.end(), comp_pred); }
		void clear() { preds.clear(); }
		int write_to_file(const string &fname);
		int read_from_file(const string &fname);
		int get_splits(vector<int> splits, PredictionList &pl);					// create a list containing only the given splits
		void get_preds_and_y(vector<float> &predictions, vector<float> &y);		// get preds and y ready to be used in performance
		void get_split_preds_and_y(int split, vector<float> &predictions, vector<float> &y); // same for a specific split

};

#endif

#endif