#ifndef __MEDBOOTSTRAP_ANALYSIS_H__
#define __MEDBOOTSTRAP_ANALYSIS_H__
#include <unordered_map>
#include "bootstrap.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"

map<string, map<string, float>> booststrap_base(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
	const map<string, vector<float>> &additional_info,
	const map<string, vector<Filter_Param>> &filter_cohort,
	const ROC_Params &function_params, float sample_ratio, int sample_per_pid, int loopCnt);

map<string, map<string, float>> booststrap(MedSamples &samples,
	map<string, vector<float>> &additional_info,
	const map<string, vector<Filter_Param>> &filter_cohort,
	const ROC_Params &roc_params, float sample_ratio = (float)0.7, int sample_per_pid = 1,
	int loopCnt = 500);

map<string, map<string, float>> booststrap(MedFeatures &features,
	const map<string, vector<Filter_Param>> &filter_cohort,
	const ROC_Params &roc_params, float sample_ratio = (float)0.7, int sample_per_pid = 1,
	int loopCnt = 500);

map<string, map<string, float>> booststrap(MedSamples &samples, const string &rep_path,
	const map<string, vector<Filter_Param>> &filter_cohort,
	const ROC_Params &roc_params, float sample_ratio = (float)0.7, int sample_per_pid = 1,
	int loopCnt = 500);

void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedSamples &samples);
void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedFeatures &features);
void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedFeatures &features);
void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedSamples &samples);

void change_sample_autosim(MedSamples &samples, int min_time, int max_time, MedSamples &new_samples);
void change_sample_autosim(MedFeatures &features, int min_time, int max_time, MedFeatures &new_features);

#endif