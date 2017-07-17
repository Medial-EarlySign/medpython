// A parent class for single-value cleaners

#ifndef _MED_VALUE_CLEANER_H_
#define _MED_VALUE_CLEANER_H_

#define NUMERICAL_CORRECTION_EPS 1e-8

typedef enum {
	VAL_CLNR_ITERATIVE,
	VAL_CLNR_QUANTILE,
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

	// Iterative
	float trimming_sd_num, removing_sd_num, nbrs_sd_num ;

	// Quantile
	float quantile, trimming_quantile_factor, removing_quantile_factor, nbrs_quantile_factor;

	// Application
	bool doTrim;
	bool doRemove;

	// Utility : maximum number of samples to take for moments calculations
	int max_samples = 10000;

};

class MedValueCleaner {
public:

	// Learning parameters
	ValueCleanerParams params;

	// Thresholds for trimming
	float trimMax, trimMin;

	// Thresholds for removing
	float removeMax, removeMin;

	// Thresholds for neighbors
	float nbrsMax, nbrsMin;

	// Functions
	// Learning 
	int get_quantile_min_max(vector<float>& values);
	int get_iterative_min_max(vector<float>& values);

	// Init
	virtual void init_defaults() { return; }
	int init(void *params);
	int init(map<string, string>& mapper);
	
	// Get Type
	ValueCleanerType get_cleaner_type(string name);
};

#endif

