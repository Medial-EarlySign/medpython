#ifndef _MED_HIST_H_
#define _MED_HIST_H_

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include <assert.h>
#include <math.h>

#include <vector>
#include <map>
#include <string>
#include <algorithm>

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


	MedHist() { clear(); }
	
	void clear();
	int get_hist(vector<float> &vals);

	int get_cdf(vector<float> &cdf);
	int sum();
	void print(string &prefix);
};

#endif