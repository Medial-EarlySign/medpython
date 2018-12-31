#ifndef __BOOTSTRAP_ANALYSIS_H__
#define __BOOTSTRAP_ANALYSIS_H__
#include <vector>
#include <string>
#include <map>
#include <random>
#include <MedTime/MedTime/MedTime.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <Logger/Logger/Logger.h>


using namespace std;

/**
* @file
* This is the infrastracture of bootstrap
*/

static MedTime med_time;

/**
* A class which fetches the samples in bootstrap manner in lazy way. \n
* The object doesn't creates vectors for labels,scores in each bootstrap but
* fetches the samples in randomization without creating memory allocation
*/
class Lazy_Iterator {
public:
	/// <summary>
	/// The Ctor
	/// @param p_pids a reference to the ids vector
	/// @param p_preds a reference to the predicitons vector
	/// @param p_y a reference to the labels vector
	/// @param p_sample_ratio a sample ratio parameter for the bootstrap (0-1]
	/// @param p_sample_per_pid a sample count per patient
	/// @param max_loops a maximal number of accessing the bootstrap 
	/// (num of threads or the bootstrap loop count).
	/// @param seed if 0 will use random device to select seed for randomization
	/// </summary>
	Lazy_Iterator(const vector<int> *p_pids, const vector<float> *p_preds,
		const vector<float> *p_y, const vector<float> *p_w, float p_sample_ratio, int p_sample_per_pid, int max_loops, int seed);

	/// <summary>
	/// Inline function to fetch next pred,label couple in the bootstrap process
	/// </summary>
	inline bool fetch_next(int thread, float &ret_y, float &ret_pred, float &weight);

	/// <summary>
	/// external function to fetch next pred,label couple in the bootstrap process for external implementitions
	/// </summary>
	bool fetch_next_external(int thread, float &ret_y, float &ret_pred, float &weight);

	/// <summary>
	/// to restart the iterator
	/// </summary>
	void restart_iterator(int thread);
	/// <summary>
	/// set the bootstrap to retrieve those vectors p_y,p_preds with no randomizations
	/// @param p_y a pointer to array labels
	/// @param p_preds a pointer to array predictions
	/// @param thread_num an access point to the bootstrap state - thread_numbeer or bootstrap loop count
	/// </summary>
	void set_static(const vector<float> *p_y, const vector<float> *p_preds, const vector<float> *p_w, int thread_num);

	~Lazy_Iterator();

	//sampling params:
	float sample_ratio; ///<the sample ratio of the patients out of all patients in each bootstrap
	int sample_per_pid; ///<how many samples to take for each patients. 0 - means no sampling take all sample for patient
	bool sample_all_no_sampling; ///<for calcing Obs if true
private:
	//internal structure - one time init
	static random_device rd;
	vector<mt19937> rd_gen;
	uniform_int_distribution<> rand_pids;
	vector<int> ind_to_pid;
	vector<vector<int>> pid_index_to_indexes; //for each pid_index retrieve the indexes in the original vectors
	vector<uniform_int_distribution<>> internal_random;
	int cohort_size;
	int min_pid_start;
	//init each time again
	//save for each Thread!
	vector<int> current_pos;
	vector<int> inner_pos; //only used when sample_per_pid==0
	vector<int> sel_pid_index; //only used when sample_per_pid==0
	vector<int> vec_size;
	vector<const float *> vec_y;
	vector<const float *> vec_preds;
	vector<const float *> vec_weights;

	//original vectors
	const float *preds;
	const float *y;
	const float *weights;
	const vector<int> *pids;


	//threading:
	int maxThreadCount;
};

#pragma region Measurements Fucntions
/// <summary>
/// A Function to calculate only NPOS,NNEG (already calculated in calc_roc_measures_with_inc). \n
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from each measurement name("NPOS" or "NNEG") to it's value
/// </returns>
map<string, float> calc_npos_nneg(Lazy_Iterator *iterator, int thread_num, void *function_params);
/// <summary>
/// A Function to calculate only AUC (already calculated in calc_roc_measures_with_inc). \n
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from measurement name "AUC" to it's value
/// </returns>
map<string, float> calc_only_auc(Lazy_Iterator *iterator, int thread_num, void *function_params);
/// <summary>
/// A Function to calculate all roc measurements- AUC, Sensitivity, speceficity 
/// positive rate, ppv...\n
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from each measurement name to it's value
/// </returns>
map<string, float> calc_roc_measures_with_inc(Lazy_Iterator *iterator, int thread_num, void *function_params); //with PPV and PR
/// <summary>
/// A Function to calculate calc_kandel_tau
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from measurement name "Kendell-Tau" to it's value
/// </returns>
map<string, float> calc_kandel_tau(Lazy_Iterator *iterator, int thread_num, void *function_params);
//For example we can put here statistical measures for regression problem or more measurements for classification..
#pragma endregion

/**
* The Incident Object which holds the gender, age incidence stats
*/
class Incident_Stats : public SerializableObject {
public:
	//age bin config:
	int age_bin_years; ///<age bin size in years
	float min_age; ///<the minimal age in the file
	float max_age; ///<the maximal age in the file
	///outcome_labels - sorted:
	vector<float> sorted_outcome_labels;
	//male:
	vector<vector<double>> male_labels_count_per_age; ///<for each age_bin, histogram of outcome labels
	//female:
	vector<vector<double>> female_labels_count_per_age; ///<for each age_bin, histogram of outcome labels

	/// Reading the file. the file format is: \n
	///   - a line with AGE_BIN[TAB]{NUMBER} \n
	///     AGE_BIN is keyword, and {NUMBER} is the age bin value numeric
	///   - a line with AGE_MIN[TAB]{NUMBER} \n
	///     AGE_MIN is keyword, and {NUMBER} is the age minimal value
	///   - a line with AGE_MAX[TAB]{NUMBER} \n
	///     AGE_MAX is keyword, and {NUMBER} is the age maximal value
	///   - a line with OUTCOME_VALUE[TAB]{NUMBER} \n
	///     OUTCOME_VALUE is keyword, and {NUMBER} is a possible outcome value.
	///     binary bootstrap will contain 2 lines with 0 and 1 values
	///   - a line with STATS_ROW[TAB]{MALE|FEMALE}[TAB]{AGE_NUMBER}[TAB]{OUTCOME_VALUE}[TAB]{NUMBER_COUNT} \n
	///     STATS_ROW is keyword, and than you provide either "MALE" or "FEMALE", TAB the age TAB the outcome
	///     value (in binary 0 for control, 1 for case) TAB the count. the incidence will calulate the 
	///     incidence rate in each group.
	void read_from_text_file(const string &text_file);
	/// Writing the file. please refer to read_from_text_file for the file format
	void write_to_text_file(const string &text_file);

	ADD_SERIALIZATION_FUNCS(age_bin_years, min_age, max_age, sorted_outcome_labels,
		male_labels_count_per_age, female_labels_count_per_age)
};

/**
* Parameter object for calc_roc_measures fucntions. this object
* stores the working point, and other parameters for the roc measurments
* bootstrap calculations.
*/
class ROC_Params : public SerializableObject {
public:
	vector<float> working_point_FPR; ///< The False Positive rate working point definition
	vector<float> working_point_SENS; ///< The True Positive rate working point definition
	vector<float> working_point_PR; ///< The Positive rate working point definition
	bool use_score_working_points; ///< If true will calculate all roc measurements based on scores working points
	float max_diff_working_point; ///< The maximal diff in calculated working point to requested working point to drop
	int score_bins; ///< score bin count for speed up calculation. 0 means no binning
	int score_min_samples; ///< score bin min sample count for speed up calculation. 0 means no limit
	float score_resolution; ///< score resultion to contorl bining for speed up calculation. 0 means no binning resulotion
	bool fix_label_to_binary; ///< If True will change label value to be binary 0,1 (default is True)
	Incident_Stats inc_stats; ///< the incedince data if provided for general population. look for Incident_Stats for more info
	/// <summary>
	/// Default Ctor
	/// </summary>
	ROC_Params() {
		max_diff_working_point = (float)0.05;
		use_score_working_points = false;
		working_point_FPR = { (float)0.1, 1, 5, 10,20,30,40,50,55,60,65,70,75,80,85,90,95 };
		score_bins = 0;
		score_resolution = 0;
		incidence_fix = 0;
		score_min_samples = 0;
		fix_label_to_binary = true;
	}
	/// <summary>
	/// Initializing each parameter from string in format: "parameter_name=value;...". \n
	/// for vectors values use "," between numbers
	/// </summary>
	ROC_Params(const string &init_string);
	/// <summary>
	/// Initializing each parameter from string in format: "parameter_name=value;...". \n
	/// for vectors values use "," between numbers
	/// @snippet bootstrap.cpp ROC_Params::init
	/// </summary>
	int init(map<string, string>& map);

	double incidence_fix; ///< The final incidence calculation on the cohort (will be calcuated)
	ADD_SERIALIZATION_FUNCS(working_point_FPR, working_point_SENS, working_point_PR, use_score_working_points,
		max_diff_working_point, score_bins, score_resolution, score_min_samples, fix_label_to_binary, inc_stats)
};

#pragma region Cohort Fucntions
/// <summary>
/// A function to filter samples based on single Filter_Param object. it's a FilterCohortFunc signature
/// </summary>
bool filter_range_param(const map<string, vector<float>> &record_info, int index, void *cohort_params); //on single param
/// <summary>
/// A function to filter samples based on multipal Filter_Param objects - in a vector with and condition
/// between each parameter range. it's a FilterCohortFunc signature
/// </summary>
bool filter_range_params(const map<string, vector<float>> &record_info, int index, void *cohort_params); //on vector of params
#pragma endregion

/**
* Parameter object for filter_params fucntions
*/
class Filter_Param : public SerializableObject { //for example Age and range for filter
public:
	string param_name; ///< The parameter name for the filtering
	float min_range; ///< the minimal range for the parameter
	float max_range; ///< the maximal range for the parameter

	/// <summary>
	/// initializing object in format: "PARAM_NAME:MIN_RANGE,MAX_RANGE". \n
	/// For example: \n
	/// Age:40,80 \n
	/// will create param_name="Age" in range 40 till 80.
	/// </summary>
	Filter_Param(const string &init_string);

	/// <summary>
	/// initializing object in format: "PARAM_NAME:MIN_RANGE,MAX_RANGE". \n
	/// For example: \n
	/// Age:40,80 \n
	/// will create param_name="Age" in range 40 till 80.
	/// </summary>
	int init_from_string(string init_string);

	/// default init function for each parameter. not the same as init_from_string!!!
	/// @snippet bootstrap.cpp Filter_Param::init
	int init(map<string, string>& map);

	/// <summary>
	/// default Ctor
	/// </summary>
	Filter_Param() {}

	ADD_SERIALIZATION_FUNCS(param_name, min_range, max_range)
};

struct ROC_And_Filter_Params {
	ROC_Params *roc_params;
	vector<Filter_Param> *filter;
};

//Infra
///Function which recieves Lazy_Iterator and the thread num for iterating the predictions and labels.
/// it also recieves function_params which are additional arguments for the function (can be working points
/// defintions for example)
typedef map<string, float>(*MeasurementFunctions)(Lazy_Iterator *iterator, int thread_num, void *function_params);
///Function which recieves map from feature name to vector of all samples value, sample index and cohort
/// definition params and return true\false if to include the sample in the cohort.
typedef bool(*FilterCohortFunc)(const map<string, vector<float>> &record_info, int index, void *cohort_params);
/// a function to process and maniplulate function params based on the given cohort - for example sotring
/// incedince information for the cohort
typedef void(*ProcessMeasurementParamFunc)(const map<string, vector<float>> &additional_info, const vector<float> &y, const vector<int> &pids, void *function_params,
	const vector<int> &filtered_indexes, const vector<float> &y_full, const vector<int> &pids_full);
/// a funtion to preprocess the prediction scores (binning for example to speed up bootstrap).
/// the function manipulate preds based on function_params
typedef void(*PreprocessScoresFunc)(vector<float> &preds, void *function_params);

#pragma region Process Measurement Param Functions
/// <summary>
/// a function to calculate the incidence in each cohort - preprocessing of function_params
/// and storing the incidence inside of it.
/// </summary>
void fix_cohort_sample_incidence(const map<string, vector<float>> &additional_info,
	const vector<float> &y, const vector<int> &pids, void *function_params,
	const vector<int> &filtered_indexes, const vector<float> &y_full, const vector<int> &pids_full);

/// <summary>
/// a function to calculate the incidence in each cohort - preprocessing of function_params
/// and storing the incidence inside of it. The old has same implementation as old bootstrap
/// only averaging incidence over the controls in the sample based on incidence in each group(age+gender)
/// </summary>
void fix_cohort_sample_incidence_old(const map<string, vector<float>> &additional_info,
	const vector<float> &y, const vector<int> &pids, void *function_params,
	const vector<int> &filtered_indexes, const vector<float> &y_full, const vector<int> &pids_full);
#pragma endregion

#pragma region Process Scores Functions
/// <summary>
/// Binning function of scores based on ROC_Params. look at score_bins,score_resolution
/// </summary>
void preprocess_bin_scores(vector<float> &preds, void *function_params);
#pragma endregion

/// <summary>
/// The main bootstrap function to run all bootstrap process with all the arguments
/// @param preds the predicitions vector
/// @param y the labels vector
/// @param pids the ids vector to sample different ids
/// @param additional_info The features to filter and create cohorts. map from feature name to 
/// it's values for each sample
/// @param filter_cohort The cohorts definition from the name to the filtering function
/// @param meas_functions The measurements to calculate for each bootstrap
/// @param cohort_params Additional parameters for the filtering cohort function. each key
/// is coresponding to the same key in filter_cohort
/// @param function_params Additional parameters for the measurements functions. in the same
/// order and size of meas_functions
/// @param process_measurments_params Function to process the function_params before running on each
/// cohort (helps to calc incidence for example)
/// @param preprocess_scores A function to preprocess all scores - for example binning the scores
/// to speedup bootstrap
/// @param preprocess_scores_params Additional parameters for the preprocess_scores function
/// @param sample_ratio A number in range (0,1] for subsampling the samples
/// @param sample_per_pid How many samples to sample on each id
/// @param loopCnt How many bootstrap to do
/// @param seed The random seed. If 0 will use random_device to create random seed
/// @param binary_outcome A flag to indicate the labels are binary (used to validate
/// The input labels)
/// </summary>
/// <returns>
/// Returns a map from each cohort name to the measurments results. each measurments results
/// is also a map from each measurement name to it's value
/// </returns>
map<string, map<string, float>> booststrap_analyze(const vector<float> &preds, const vector<float> &y, const vector<float> *weights
	, const vector<int> &pids, const map<string, vector<float>> &additional_info, const map<string, FilterCohortFunc> &filter_cohort,
	const vector<MeasurementFunctions> &meas_functions = { calc_roc_measures_with_inc },
	const map<string, void *> *cohort_params = NULL, const vector<void *> *function_params = NULL,
	ProcessMeasurementParamFunc process_measurments_params = NULL,
	PreprocessScoresFunc preprocess_scores = NULL, void *preprocess_scores_params = NULL,
	float sample_ratio = (float)1.0, int sample_per_pid = 1,
	int loopCnt = 500, int seed = 0, bool binary_outcome = true);

/// <summary>
/// Will output the bootstrap results into file in TAB delimeted format. each line is cohort and the
/// The columns are the measurements
/// </summary>
void write_bootstrap_results(const string &file_name, const map<string, map<string, float>> &all_cohorts_measurments, const string& run_id = "");
/// <summary>
/// Will read the bootstrap results from file in TAB delimeted format. each line is cohort and the
/// The columns are the measurements
/// </summary>
void read_bootstrap_results(const string &file_name, map<string, map<string, float>> &all_cohorts_measurments);

/// <summary>
/// Will output the bootstrap results into file with the new format with columns: 
/// "Cohort$Measure_Name", "Value"
/// </summary>
void write_pivot_bootstrap_results(const string &file_name, const map<string, map<string, float>> &all_cohorts_measurments, const string& run_id = "");
/// <summary>
/// Will read the bootstrap results into file with the new format with columns: 
/// "Cohort$Measure_Name", "Value"
/// <//summary>
void read_pivot_bootstrap_results(const string &file_name, map<string, map<string, float>> &all_cohorts_measurments);

MEDSERIALIZE_SUPPORT(Incident_Stats)
MEDSERIALIZE_SUPPORT(ROC_Params)
MEDSERIALIZE_SUPPORT(Filter_Param)

#endif // !__BOOTSTRAP_ANALYSIS_H__

