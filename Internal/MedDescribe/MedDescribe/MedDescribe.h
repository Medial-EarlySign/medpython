#pragma once

//
// MedDescribe
//
// Contains tools to study the followings:
// (1) Distributions of signals.
// (2) Biases of signals.
// (3) Joint distribution of signals.
// (4) Joint distributions of outcomes, and signals vs. outcomes.
//

#include <string>
#include <map>
#include <vector>
#include <InfraMed/InfraMed/InfraMed.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>

using namespace std;

class VecMoments {

public:
	vector<float> quantiles = { (float)0.01, (float)0.1, (float)0.25, (float)0.5, (float)0.75, (float)0.9, (float)0.99 };
	float min_val = (float)-1e6;
	float max_val = (float)1e6;
	vector<float> quantiles_vals;
	float mean;
	float std;
	float median;
	int N;
	int N_pos, N_neg;
	int N0;
	int N_diff_vals;
	int N_out_of_range;
	float vmin, vmax;

	int get_for_vec(vector<float> &v);
	int safe_get_for_vec(vector<float> v) { return get_for_vec(v); }
	void print(const string prefix);
};


class MedMutualDist : public SerializableObject {

public:
	// signal names
	string sig1;
	string sig2;

	vector<int> pids_to_check;
	MedMutualDist() {};
	~MedMutualDist() {};

	// universal sigs
	int sig1_time_ch = 0;
	int sig1_val_ch = 0;
	int sig2_time_ch = 0;
	int sig2_val_ch = 0;

	// binning
	int sig1_is_categorial = 0;
	int auto_range = 1;
	float min_val1 = (float)-1e10, max_val1 = (float)1e10;
	float min_val2 = (float)-1e10, max_val2 = (float)1e10;
	float bin_size = (float)0.1;

	// general restrictions
	int min_time = 0;
	int max_time = 1999999999;
	int gender_mask = 0x3;

	// sampling params
	int jump1_time = 180;			// minimal jump between sig1 samples, one random sample will be taken in each bucket.
	int win_from = 60, win_to = -60;	// window around sampling point to search for sig2
	int win_time_unit = MedTime::Days;

	// collected values
	map<float, vector<float>> values;

	// functions

	// init_from_string :: implemented in SerializableObject
	int init(map<string, string>& mapper);
	int collect_values(MedRepository &rep);
	int collect_values(const string &rep_fname);

};