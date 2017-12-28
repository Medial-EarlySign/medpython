//
// MedStat - Statistics utilities :
//			 1. MedCleaner used for normalization and handling outliers
//			 2. MedPerformance used for analyzing performance
//

#ifndef _MED_STAT_H_
#define _MED_STAT_H_

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include "assert.h"
#include "math.h"

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cstring>

#include "string.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"


#define MED_CLEANER_MAX_Z 15
#define MED_CLEANER_EPSILON 0.0001

using namespace std ;

//=============================================
// general helper routines
//=============================================

// get a vector of values, a vector of probabilities, and returning a matching vector of values such that Prob(x<=pvals[i])=p[i]
void get_percentiles(vector<float> &vals, vector<float> &p, vector<float> &pvals);
void get_percentiles(vector<float> &vals, vector<float> &p, vector<float> &pvals, int only_positive_flag);

// gets a vector of values, and checks the best way to round it at the 10^-3 to 10^3 (at 10x steps) range.
float get_best_rounding(vector<float> &vals, vector<float>& res, vector<int> &counts, float missing_value = -1);
//float get_best_rounding(vector<float> &vals) { vector<float> res; vector<int> counts; return(get_best_rounding(vals, res, counts)); }

// Moments
void get_mean(vector<float> &vals, float &mean);

template <typename T> void get_median(const vector<T>& vals, T &median) {
	vector<T> tempValues = vals;
	sort_and_get_median(tempValues, median);
}

template <typename T> void sort_and_get_median(vector<T>& vals, T &median) {
	if (vals.size() == 0)
		median = T();
	else {
		nth_element(vals.begin(), vals.begin() + vals.size() / 2, vals.end());
		median = vals[vals.size() / 2];
	}
}

void get_common(vector<float> &vals, float &common);

void get_mean_and_std(vector<float> &vals, float &mean, float &std);
void get_mean_and_std(vector<float> &vals, float &mean, float &std, float missing_val);


// rounding
#define ROUNDING_EPSILON 0.00001
inline float roundf(float val, float rounder) {
	return ((float)((int)((val+(float)ROUNDING_EPSILON)/rounder)) * rounder);
}


//======================================================================================
// MedCleaner - class to handle cleaning/normalization of data
//======================================================================================

#define MED_DEFAULT_MISSING_VALUE		-1
#define MED_DEFAULT_MIN_TRIM			-1e9
#define MED_DEFAULT_MAX_TRIM			 1e9
// Cleaner class : Normalizing and cleaning of outliers
class MedCleaner : SerializableObject {
public:

	float missing_value ;
	float min_trim;

	int n,nvals,most_common_count,nzeros ;
	float median,q1,q3,iqr,mean,sdv,skew,min,max ;
	float most_common_value ;

	bool trim_flag,remove_flag,normalize_flag,replace_missing_to_mean_flag;
	float trim_min,trim_max ;
	float remove_min,remove_max ;

	float sk ;
	int skew_sign ;

	MedCleaner() {missing_value = MED_DEFAULT_MISSING_VALUE; remove_flag = trim_flag = normalize_flag = replace_missing_to_mean_flag = true; n=0; nvals=0;
				  most_common_count=0; nzeros=0; median=q1=q3=iqr=MED_DEFAULT_MISSING_VALUE; median = MED_DEFAULT_MISSING_VALUE;mean=MED_DEFAULT_MISSING_VALUE; sdv=0; skew=0;
				  most_common_value = MED_DEFAULT_MISSING_VALUE; trim_min=MED_DEFAULT_MISSING_VALUE;  min_trim = MED_DEFAULT_MIN_TRIM;
				  trim_max=MED_DEFAULT_MISSING_VALUE; remove_min=MED_DEFAULT_MISSING_VALUE; remove_max=-1; sk=0; skew_sign=0;}

	void print(const string& prefix) ;
	void print_short(const string& prefix) ;
	void calculate(vector<float> &values) ; 
	void get_mean_and_sdv(vector<float> &values, bool take_missing_into_account = false) ;
	void get_cleaning_range (vector<float>& values, float& min_val, float& max_val, float std_mult = MED_CLEANER_MAX_Z) ;
	void get_limits_iteratively(vector<float> values, float std_mult = MED_CLEANER_MAX_Z) ;
	void get_cleaning_params(vector<float> values) ;
	int clear(vector<float>& values) ;
	int clean(vector<float>& values) {return clear(values);};
	void remove_trim_replace(vector<float> &values);
	
	void normalize(vector<float>& values) ;

	bool is_valid(float value) ;
	float get_trimmed(float value) ;
	float get_value(float value) ;
	int trim(float& value) ;
	void single_remove_trim_replace(float &val);
	void single_normalize(float &val);

	size_t get_size() ;
	size_t serialize(unsigned char *buffer) ;
	size_t deserialize(unsigned char *buffer) ;
} ;

MEDSERIALIZE_SUPPORT(MedCleaner)


//======================================================================================
// MedHist - simple class to get histograms
//======================================================================================
class MedHist {

	public:

	float missing_value;
	int positive_only_flag;
	float rounder; // to round values with, -1: no rounding
	float from_val, to_val; // histogram boundaries (first cell will be <from_val, last cell will be >= to_val), if not set will be the 0.001 - 0.999 parts of the distribution
	int ncells;	// if too big due to rounding values, will be shrinked automatically
	float min_percentile, max_percentile;
	vector<int> hist;
	vector<float> from;
	vector<float> to;

	void clear() {hist.clear(); from.clear(); to.clear(); missing_value = (float)MED_DEFAULT_MISSING_VALUE; from_val = 0; to_val = -1; ncells = 200; rounder = -1; positive_only_flag = 1;
				  min_percentile = (float)0.001; max_percentile = (float)0.999;}
	MedHist() {clear();}

	int get_hist(vector<float> &vals);

	int get_cdf(vector<float> &cdf);
	int sum();
	void print(string &prefix);
};


#endif
