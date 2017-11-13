#ifndef __MEDBOOTSTRAP_ANALYSIS_H__
#define __MEDBOOTSTRAP_ANALYSIS_H__
#include <unordered_map>
#include "bootstrap.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"

class MedBootstrap : public SerializableObject {
public:
	ROC_Params roc_Params;
	map<string, vector<Filter_Param>> filter_cohort;
	float sample_ratio;
	int sample_per_pid;
	int loopCnt;

	MedBootstrap();

	map<string, map<string, float>> booststrap(MedFeatures &features);
	map<string, map<string, float>> booststrap(MedSamples &samples, map<string, vector<float>> &additional_info);
	map<string, map<string, float>> booststrap(MedSamples &samples, const string &rep_path);

	void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedSamples &samples);
	void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedFeatures &features);
	void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedFeatures &features);
	void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedSamples &samples);

	void change_sample_autosim(MedSamples &samples, int min_time, int max_time, MedSamples &new_samples);
	void change_sample_autosim(MedFeatures &features, int min_time, int max_time, MedFeatures &new_features);

	ADD_SERIALIZATION_FUNCS(sample_ratio, sample_per_pid, loopCnt, roc_Params, filter_cohort);

private:
	map<string, map<string, float>> booststrap_base(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
		const map<string, vector<float>> &additional_info);
	bool use_time_window();
};

class MedBootstrapResult : public SerializableObject {
public:
	MedBootstrap bootstrap_params;
	map<string, map<string, float>> bootstrap_results;

	void booststrap(MedFeatures &features);
	void booststrap(MedSamples &samples, map<string, vector<float>> &additional_info);
	void booststrap(MedSamples &samples, const string &rep_path);

	void find_working_points(const map<string, float> &bootstrap_cohort,
		vector<float> &sens_points, vector<float> &pr_points);

	void explore_score(float score, map<string, float> &score_measurements,
		const string &string_cohort = "All", float max_search_range = 0.1);

	ADD_SERIALIZATION_FUNCS(bootstrap_params, bootstrap_results)
private:
	bool find_in_range(const vector<float> &vec, float search, float th);
	void explore_measure(const string &measure_name, float value, map<string, float> &score_measurements,
		const string &string_cohort = "All", float max_search_range = 0.1);
};

MEDSERIALIZE_SUPPORT(MedBootstrap)
MEDSERIALIZE_SUPPORT(MedBootstrapResult)

#endif