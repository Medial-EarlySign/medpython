#ifndef _MED_VALUE_CLEANER_H_
#define _MED_VALUE_CLEANER_H_

#define NUMERICAL_CORRECTION_EPS 1e-8
#undef max
#undef min

/* @enum
* Basic Cleaner types
*/
typedef enum {
	VAL_CLNR_ITERATIVE, ///<"iterative"
	VAL_CLNR_QUANTILE, ///<"quantile"
	VAL_CLNR_LAST,
} ValueCleanerType;

class ValueCleanerParams {
public:
	ValueCleanerType type;

	// General
	int take_log;
	float missing_value;
	float range_min = (float)-1e20;
	float range_max = (float)1e20;
	float trim_range_min = (float)-1e21;
	float trim_range_max = (float)1e21;

	// Iterative
	float trimming_sd_num, removing_sd_num, nbrs_sd_num ;

	// Quantile
	float quantile, trimming_quantile_factor, removing_quantile_factor, nbrs_quantile_factor;

	// Application
	bool doTrim;
	bool doRemove;

	/// Utility : maximum number of samples to take for moments calculations
	int max_samples = 10000;

	ValueCleanerParams() {
		//defautls
		quantile = 0;
		removing_quantile_factor = 1;
		nbrs_quantile_factor = 0;
		trimming_quantile_factor = 1;
		doTrim = doRemove = true;
		take_log = 0;
		trimming_sd_num = 7;
		removing_sd_num = 14;
		nbrs_sd_num = 0;
		missing_value = -65336;
	}

};

/** @file
*  A parent class for single-value cleaners
*/
class MedValueCleaner {
public:

	/// Learning parameters
	ValueCleanerParams params;

	/// Thresholds for trimming
	float trimMax, trimMin;

	/// Thresholds for removing
	float removeMax, removeMin;

	/// Thresholds for neighbors
	float nbrsMax, nbrsMin;

	int num_samples_after_cleaning;

	// Functions
	/// Learning 
	int get_quantile_min_max(vector<float>& values);
	int get_iterative_min_max(vector<float>& values);

	// Init
	virtual void init_defaults() { return; }
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedValueCleaner.cpp MedValueCleaner::init
	int init(map<string, string>& mapper);
	
	/// Get Type
	ValueCleanerType get_cleaner_type(string name);

	///default ctor:
	MedValueCleaner() {
		trimMin = numeric_limits<float>().min();
		trimMax = numeric_limits<float>().max();
		removeMin = numeric_limits<float>().min();
		removeMax = numeric_limits<float>().max();
		nbrsMin = numeric_limits<float>().min();
		nbrsMax = numeric_limits<float>().max();
	}
};

#endif

