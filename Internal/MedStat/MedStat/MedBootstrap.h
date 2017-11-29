#ifndef __MEDBOOTSTRAP_ANALYSIS_H__
#define __MEDBOOTSTRAP_ANALYSIS_H__
#include <unordered_map>
#include "bootstrap.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "InfraMed/InfraMed/MedPidRepository.h"

class MedBootstrap : public SerializableObject {
public:
	ROC_Params roc_Params;
	map<string, vector<Filter_Param>> filter_cohort;
	float sample_ratio; //the sample ratio of the patients out of all patients in each bootstrap
	int sample_per_pid; //how many samples to take for each patients. 0 - means no sampling take all sample for patient
	bool sample_patient_label; //if true will treat patient+label as the "id" for the sampling
	int loopCnt; //the bootstrap count

	void add_filter_cohorts(const map<string, vector<pair<float, float>>> &parameters_ranges);
	void add_filter_cohorts(const vector<vector<Filter_Param>> &parameters_ranges);
	
	void parse_cohort_file(const string &cohorts_path);

	MedBootstrap();
	//init_string format: paramter_name=value;... roc_Params is the init_string of roc_Params
	//filter_cohort - it's a file where each line is new cohort. "cohort_name TAB paramter_name:min_range,max_range;..."
	MedBootstrap(const string &init_string);

	void clean_feature_name_prefix(map<string, vector<float>> &features);

	//if results_per_split is NULL (not provided) will not calculate results by each split. the return value is on all values without splits
	void prepare_bootstrap(MedFeatures &features, vector<float> &preds, vector<float> &y, vector<int> &pids,
		map<string, vector<float>> &final_additional_info, unordered_map<int, vector<int>> *splits_inds = NULL);
	void prepare_bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info, vector<float> &preds, vector<float> &y, vector<int> &pids,
		unordered_map<int, vector<int>> *splits_inds = NULL);
	map<string, map<string, float>> bootstrap(MedFeatures &features, map<int, map<string, map<string, float>>> *results_per_split = NULL); 
	map<string, map<string, float>> bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info, map<int, map<string, map<string, float>>> *results_per_split = NULL);
	map<string, map<string, float>> bootstrap(MedSamples &samples, const string &rep_path, map<int, map<string, map<string, float>>> *results_per_split = NULL);
	map<string, map<string, float>> bootstrap(MedSamples &samples, MedPidRepository &rep, map<int, map<string, map<string, float>>> *results_per_split = NULL);

	void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedSamples &samples);
	void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedFeatures &features);
	void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedFeatures &features);
	void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedSamples &samples);

	void change_sample_autosim(MedSamples &samples, int min_time, int max_time, MedSamples &new_samples);
	void change_sample_autosim(MedFeatures &features, int min_time, int max_time, MedFeatures &new_features);

	ADD_SERIALIZATION_FUNCS(sample_ratio, sample_per_pid, loopCnt, roc_Params, filter_cohort);

private:
	map<string, map<string, float>> bootstrap_base(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
		const map<string, vector<float>> &additional_info);
	void add_splits_results(const vector<float> &preds, const vector<float> &y,
		const vector<int> &pids, const map<string, vector<float>> &data,
		const unordered_map<int, vector<int>> &splits_inds,
		map<int, map<string, map<string, float>>> &results_per_split);
	bool use_time_window();
};

class MedBootstrapResult : public SerializableObject {
public:
	MedBootstrap bootstrap_params;
	map<string, map<string, float>> bootstrap_results;

	void bootstrap(MedFeatures &features);
	void bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info);
	void bootstrap(MedSamples &samples, const string &rep_path);

	void find_working_points(const map<string, float> &bootstrap_cohort,
		vector<float> &sens_points, vector<float> &pr_points);

	void explore_score(float score, map<string, float> &score_measurements,
		const string &string_cohort = "All", float max_search_range = 0.1);

	void write_results_to_text_file(const string &path, bool pivot_format = true);
	void read_results_to_text_file(const string &path, bool pivot_format = true);

	ADD_SERIALIZATION_FUNCS(bootstrap_params, bootstrap_results)
private:
	bool find_in_range(const vector<float> &vec, float search, float th);
	void explore_measure(const string &measure_name, float value, map<string, float> &score_measurements,
		const string &string_cohort = "All", float max_search_range = 0.1);
};

MEDSERIALIZE_SUPPORT(MedBootstrap)
MEDSERIALIZE_SUPPORT(MedBootstrapResult)

#endif