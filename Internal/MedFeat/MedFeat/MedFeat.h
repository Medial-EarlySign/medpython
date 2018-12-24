//
// MedFeat.h - Tools to help build features, and for feature selection
//

#ifndef __MEDFEAT_H__
#define __MEDFEAT_H__

#if 0

#include "Logger/Logger/Logger.h"
#include "InfraMed/InfraMed/MedSignals.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "MedStat/MedStat/MedStat.h"

//===========================================================
// SDateVal tools and features
//===========================================================

// SDateVal time query
int sdv_get_first_before(SDateVal *sdv, int len, int date, int before); // returns -1 or the index in the sdv array that is the maximal that is <= date-before_days in its date

// SDateVal get sub list
SDateVal *sdv_get_sub_list(SDateVal *sdv, int len, int date, int days_before_start, int days_before_end, int &lout);

// SDateVal general Features
float sdv_get_min(SDateVal *sdv, int len); // returns -1 if len is 0
float sdv_get_max(SDateVal *sdv, int len); // returns -1 if len is 0
float sdv_get_avg(SDateVal *sdv, int len);
float sdv_get_std(SDateVal *sdv, int len);
float sdv_get_fraction_below(SDateVal *sdv, int len, float bound);
float sdv_get_fraction_above(SDateVal *sdv, int len, float bound);
float sdv_get_slope(SDateVal *sdv, int len); // slope is 0 (= no change) if there are less than 2 samples
float sdv_get_linear_delta(SDateVal *sdv, int len, int days);
float sdv_get_linear_val(SDateVal *sdv, int len, int days);

// Next Class is a basic class to hold a list of "named" features.

class Sample_Id {

	public:
		string name;
		int id;
		int date;
};

class MedFeaturesData {
	public:

		int nsplits ;
		map<string,vector<float>> data ;			// from features name -> feature vector
		vector<float> label ;						// y 
		map<string,vector<MedCleaner>> cleaners;
		vector<MedCleaner> label_cleaners;			// y cleaner
		string label_name;							// y name

		vector<string> signals;						// feature names to take
		vector<int> splits ;

		vector<Sample_Id> row_ids;					// optional place to hold id of each sample

		// Predictions 
		int n_preds_per_sample;
		bool predict_on_train;
		vector<vector<float> > split_preds ;
		vector<vector<float> > split_preds_on_train;
		vector<float> preds ;

		//MedFeaturesData(const int& _nsplits) { nsplits = _nsplits; n_preds_per_sample = 1; predict_on_train = false; }
		MedFeaturesData(const int _nsplits) { nsplits = _nsplits; n_preds_per_sample = 1; predict_on_train = false; }
		MedFeaturesData() { MedFeaturesData(1); }

		void clear() { data.clear(); label.clear(); cleaners.clear(); label_cleaners.clear(); label_name.clear(); signals.clear(); splits.clear(); split_preds.clear(); split_preds_on_train.clear(); preds.clear(); }
		int loadFromRepository(MedRepository &rep, string label_name, vector<string> sigs_name);

		// Functions 
		int get_learning_nrows(int isplit) ;
		int get_testing_nrows(int isplit) ;

		void setup_cleaners_for_all_folds(const string& signal, bool trim_flag, bool remove_flag, bool replace_missing_to_mean_flag, bool normalize_flag, float missing_v = MED_DEFAULT_MISSING_VALUE);
		void get_normalization(const string& signal) ;
		void get_label_normalization() ;
		void get_normalization_and_cleaning(const string& signal, float missing_v) ;
		void get_normalization_and_cleaning(const string& signal) {return get_normalization_and_cleaning(signal,-1);}
		void get_all_normalization_and_cleaning(float missing_v);

		int apply_clean_and_normalize(int i_split, float missing_v, string &sig, int normalize_only_flag); // for a single split - if == nsplits then on all
		int apply_clean_and_normalize(float missing_v, string &sig, int i_cleaner, int normalize_only_flag); // on each split

		const static int Transpose = 0x0001;
		const static int Normalize = 0x0002;
		const static int Clean = 0x0004;
		const static int Split_Equal = 0x0008;
		// translating Features to a MedMat matrix with the following options:
		// clean_flag: will remove(that replace with mean), trim, and complete (put mean).
		// transpose_flag: load matrix transposed
		// normalize_flag: normalize with mean,std in cleaner
		// split_eq: take equal or different from split (predict or learn)
		// _split: to take
		void get_features_as_mat(MedMat<float> &x, int _isplit, int flags) {get_features_as_mat(signals,x,_isplit,flags);}
		void get_features_as_mat(vector<string> &f_names, MedMat<float> &x, int _isplit, int flags);
		void get_features_as_mat(const string &f_name, MedMat<float> &x, int _isplit, int flags) {vector<string> vs; vs.push_back(f_name); get_features_as_mat(vs,x,_isplit,flags);}
		void get_label_as_mat(MedMat<float> &y, int _isplit, int flags);

		// note that when write/read cleaners one has to make sure to call the matching read with exactly the same nsplits and exactly the same order of signals in the vector sigs.
		int write_cleaners_to_file(const string &fname, vector<string> &sigs);
		int read_cleaners_from_file(const string &fname, vector<string> &sigs);
		int write_label_cleaners_to_file(const string &fname);
		int read_label_cleaners_to_file(const string &fname);
} ;

#endif

#endif

