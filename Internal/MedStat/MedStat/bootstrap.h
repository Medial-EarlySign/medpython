#ifndef __BOOTSTRAP_ANALYSIS_H__
#define __BOOTSTRAP_ANALYSIS_H__

#include "MedUtils/MedUtils/MedUtils.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include <vector>
#include <string>
#include <map>
#include "Logger/Logger/Logger.h"

using namespace std;

static MedTime med_time;

#pragma region Measurements Fucntions
map<string, float> calc_npos_nneg(const vector<float> &preds, const vector<float> &y, void *function_params);
map<string, float> calc_only_auc(const vector<float> &preds, const vector<float> &y, void *function_params);
map<string, float> calc_roc_measures(const vector<float> &preds, const vector<float> &y, void *function_params); //SENS, SPEC, SCORE. no ppv, no PR
map<string, float> calc_roc_measures_full(const vector<float> &preds, const vector<float> &y, void *function_params); //with PPV and PR
map<string, float> calc_roc_measures_with_inc(const vector<float> &preds, const vector<float> &y, void *function_params); //with PPV and PR
//For example we can put here statistical measures for regression problem or more measurements for classification..
#pragma endregion

//The Incident Object:
class Incident_Stats : public SerializableObject {
public:
	//age bin config:
	int age_bin_years;
	float min_age, max_age;
	//outcome_labels - sorted:
	vector<float> sorted_outcome_labels;
	//male:
	vector<vector<double>> male_labels_count_per_age; //for each age_bin, histogram of outcome labels
	//female:
	vector<vector<double>> female_labels_count_per_age; //for each age_bin, histogram of outcome labels

	ADD_SERIALIZATION_FUNCS(age_bin_years, min_age, max_age, sorted_outcome_labels, 
		male_labels_count_per_age, female_labels_count_per_age)
};

//Parameter object for calc_roc_measures fucntions:
class ROC_Params : public SerializableObject {
public:
	vector<float> working_point_FPR;
	vector<float> working_point_SENS;
	vector<float> working_point_PR;
	bool use_score_working_points; //if true will calculate all roc measurements based on scores working points
	float max_diff_working_point; //the maximal diff in calculated working point to requested working point to drop
	int score_bins;
	Incident_Stats inc_stats; //the incedince data if provided for general population
	ROC_Params() {
		max_diff_working_point = (float)0.05;
		use_score_working_points = false;
		working_point_FPR = { (float)0.1, 1, 5, 10,20,30,40,50,55,60,65,70,75,80,85,90,95 };
		score_bins = 0;
		incidence_fix = 0;
	}
	double incidence_fix;
	
	ADD_SERIALIZATION_FUNCS(working_point_FPR, working_point_SENS, working_point_PR, use_score_working_points,
		max_diff_working_point, score_bins, inc_stats)
};

#pragma region Cohort Fucntions
bool filter_range_param(const map<string, vector<float>> &record_info, int index, void *cohort_params); //on single param
bool filter_range_params(const map<string, vector<float>> &record_info, int index, void *cohort_params); //on vector of params
#pragma endregion

//Parameter object for filter_params fucntion:
class Filter_Param : public SerializableObject { //for example Age and range for filter
public:
	string param_name;
	float min_range, max_range;

	ADD_SERIALIZATION_FUNCS(param_name, min_range, max_range)
};

//Infra
typedef map<string, float>(*MeasurementFunctions)(const vector<float> &preds, const vector<float> &y, void *function_params);
typedef bool(*FilterCohortFunc)(const map<string, vector<float>> &record_info, int index, void *cohort_params);
typedef void(*ProcessMeasurementParamFunc)(const map<string, vector<float>> &additional_info, const vector<float> &y, const vector<int> &pids, FilterCohortFunc cohort_def, void *cohort_params, void *function_params);

#pragma region Process Measurement Param Functions
void fix_cohort_sample_incidence(const map<string, vector<float>> &additional_info,
	const vector<float> &y, const vector<int> &pids, FilterCohortFunc cohort_def,
	void *cohort_params, void *function_params);
#pragma endregion

//The bootstrap function process
map<string, map<string, float>> booststrap_analyze(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
	const map<string, vector<float>> &additional_info, const map<string, FilterCohortFunc> &filter_cohort,
	const vector<MeasurementFunctions> &meas_functions = { calc_roc_measures_with_inc },
	const map<string, void *> *cohort_params = NULL, const vector<void *> *function_params = NULL,
	ProcessMeasurementParamFunc process_measurments_params = NULL,
	float sample_ratio = (float)1.0, int sample_per_pid = 1,
	int loopCnt = 500, bool binary_outcome = true);

//will output the bootstrap results  into file
void PrintMeasurement(const map<string, map<string, float>> &all_cohorts_measurments, const string &file_name);

MEDSERIALIZE_SUPPORT(Incident_Stats)
MEDSERIALIZE_SUPPORT(ROC_Params)
MEDSERIALIZE_SUPPORT(Filter_Param)

#endif // !__BOOTSTRAP_ANALYSIS_H__

