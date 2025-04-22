#ifndef __BOOTSTRAP_ANALYSIS_H__
#define __BOOTSTRAP_ANALYSIS_H__
#include <vector>
#include <string>
#include <map>
#include <random>
#include <MedTime/MedTime/MedTime.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <Logger/Logger/Logger.h>
#include <MedStat/MedStat/MedStat.h>

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
class Lazy_Iterator
{
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
				  const vector<float> *p_y, const vector<float> *p_w, float p_sample_ratio, int p_sample_per_pid, int max_loops, int seed, const vector<int> *p_preds_order = NULL);

	void init(const vector<int> *p_pids, const vector<float> *p_preds,
			  const vector<float> *p_y, const vector<float> *p_w, float p_sample_ratio, int p_sample_per_pid, int max_loops, int seed);

	/// <summary>
	/// Inline function to fetch next pred,label couple in the bootstrap process
	/// </summary>
	inline bool fetch_next(int thread, float &ret_y, float &ret_pred, float &weight);
	inline bool fetch_next(int thread, float &ret_y, const float *&ret_pred, float &weight, const int *&preds_order);
	/// <summary>
	/// external function to fetch next pred,label couple in the bootstrap process for external implementitions
	/// </summary>
	bool fetch_next_external(int thread, float &ret_y, float &ret_pred, float &weight);
	/// <summary>
	/// external function to fetch next pred,label couple in the bootstrap process for external implementitions
	/// </summary>
	bool fetch_next_external(int thread, float &ret_y, float &ret_pred, float &weight, const int *&preds_order);

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
	void set_static(const vector<float> *p_y, const vector<float> *p_preds, const vector<float> *p_w, const vector<int> *p_preds_order, int thread_num);

	void set_thread_sample_all(int thread);

	~Lazy_Iterator();

	// sampling params:
	float sample_ratio;			 ///< the sample ratio of the patients out of all patients in each bootstrap
	int sample_per_pid;			 ///< how many samples to take for each patients. 0 - means no sampling take all sample for patient
	bool sample_all_no_sampling; ///< for calcing Obs if true
	size_t num_categories;		 ///< number of categories (inferred)
private:
	// internal structure - one time init
	static random_device rd;
	vector<mt19937> rd_gen;
	uniform_int_distribution<> rand_pids;
	vector<int> ind_to_pid;
	vector<vector<int>> pid_index_to_indexes; // for each pid_index retrieve the indexes in the original vectors
	vector<uniform_int_distribution<>> internal_random;
	int cohort_size;
	int min_pid_start;
	// init each time again
	// save for each Thread!
	vector<int> current_pos;
	vector<int> inner_pos;	   // only used when sample_per_pid==0
	vector<int> sel_pid_index; // only used when sample_per_pid==0
	vector<int> vec_size;
	vector<const float *> vec_y;
	vector<const float *> vec_preds;
	vector<const float *> vec_weights;
	vector<const int *> vec_preds_order;
	// original vectors
	const float *preds;
	const float *y;
	const float *weights;
	const vector<int> *pids;
	const int *preds_order;

	// threading:
	int maxThreadCount;
};

/**
 * A class which fetches the samples in bootstrap manner in memort way. \n
 * The object selects indexes in each bootstrap
 * fetches the samples in randomization by index vector
 * in the future - choose dynamically the faster Iterator to use Lazy_Iterator or Mem_Iterator
 */
class Mem_Iterator
{
public:
	// no calc
	Mem_Iterator() {}

	/// <summary>
	/// The Ctor
	/// @param pids a reference to pids vector - without selection
	/// @param cohort_indexes a reference to selected indexes from cohort filter to bootstrap indexes from
	/// @param p_sample_ratio a sample ratio parameter for the bootstrap (0-1]
	/// @param p_sample_per_pid a sample count per patient
	/// (num of threads or the bootstrap loop count).
	/// @param seed if 0 will use random device to select seed for randomization
	/// </summary>
	Mem_Iterator(const vector<int> &pids, const vector<int> &cohort_indexes, float p_sample_ratio, int p_sample_per_pid, int seed);

	/// <summary>
	/// Inline function to fetch indexes
	/// </summary>
	inline void fetch_selection(vector<int> &indexes) const;

	/// <summary>
	/// external function to fetch indexes
	/// </summary>
	void fetch_selection_external(vector<int> &indexes) const;

	/// <summary>
	/// Inline function to fetch indexes - for multi thread - provide random generator
	/// </summary>
	inline void fetch_selection(mt19937 &rd_gen, vector<int> &indexes) const;

	/// <summary>
	/// external function to fetch indexes  - for multi thread - provide random generator
	/// </summary>
	void fetch_selection_external(mt19937 &rd_gen, vector<int> &indexes) const;

	// sampling params:
	float sample_ratio; ///< the sample ratio of the patients out of all patients in each bootstrap
	int sample_per_pid; ///< how many samples to take for each patients. 0 - means no sampling take all sample for patient
private:
	// internal structure - one time init
	mt19937 _rd_gen;

	int cohort_size;
	int tot_rec_cnt;
	vector<vector<int>> pid_to_inds;
	vector<int> ind_to_pid;
	vector<int> cohort_idx;
};

/**
 * A base class for measurements parameter
 */
class Measurement_Params : public SerializableObject
{
public:
	bool show_warns; ///< If True will show warnnings
	Measurement_Params();
	virtual ~Measurement_Params() {};
};

#pragma region Measurements Functions
/// <summary>
/// A Function to calculate only NPOS,NNEG (already calculated in calc_roc_measures_with_inc). \n
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from each measurement name("NPOS" or "NNEG") to it's value
/// </returns>
map<string, float> calc_npos_nneg(Lazy_Iterator *iterator, int thread_num, Measurement_Params *function_params);
/// <summary>
/// A Function to calculate only AUC (already calculated in calc_roc_measures_with_inc). \n
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from measurement name "AUC" to it's value
/// </returns>
map<string, float> calc_only_auc(Lazy_Iterator *iterator, int thread_num, Measurement_Params *function_params);
/// <summary>
/// A Function to calculate all roc measurements- AUC, Sensitivity, speceficity
/// positive rate, ppv...\n
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from each measurement name to it's value
/// </returns>
map<string, float> calc_roc_measures_with_inc(Lazy_Iterator *iterator, int thread_num, Measurement_Params *function_params); // with PPV and PR
/// <summary>
/// A Function to calculate calc_kandel_tau
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from measurement name "Kendell-Tau" to it's value
/// </returns>
map<string, float> calc_kandel_tau(Lazy_Iterator *iterator, int thread_num, Measurement_Params *function_params);
/// <summary>
/// A Function to calculate performance measurements for multicategory
/// Implements MeasurementFunctions signature function
/// </summary>
/// <returns>
/// A map from measurement name to it's value
/// </returns>
map<string, float> calc_multi_class(Lazy_Iterator *iterator, int thread_num, Measurement_Params *function_params);

/// <summary>
/// A Function to calculate performance measurements for harrell c statistic
/// Implements MeasurementFunctions signature function
/// Encoding:
/// Case/Control => effect outcome/y sign. positive is case, negative controls. Can't handle event in time zero.
/// Time to event => abs value of outcome/y
/// Score => the prediction
/// </summary>
/// <returns>
/// A map from measurement name to it's value
/// </returns>
map<string, float> calc_harrell_c_statistic(Lazy_Iterator *iterator, int thread_num, Measurement_Params *function_params);

/// <summary>
/// A Function to calculate performance measurements for regression problems
/// Implements MeasurementFunctions signature function
/// Accepted Regresion_Params
/// </summary>
/// <returns>
/// A map from measurement name to it's value
/// </returns>
map<string, float> calc_regression(Lazy_Iterator *iterator, int thread_num, Measurement_Params *function_params);

// For example we can put here statistical measures for regression problem or more measurements for classification..
#pragma endregion

/**
 * The Incident Object which holds the gender, age incidence stats
 */
class Incident_Stats : public SerializableObject
{
public:
	// age bin config:
	int age_bin_years; ///< age bin size in years
	float min_age;	   ///< the minimal age in the file
	float max_age;	   ///< the maximal age in the file
	/// outcome_labels - sorted:
	vector<float> sorted_outcome_labels;
	// male:
	vector<vector<double>> male_labels_count_per_age; ///< for each age_bin, histogram of outcome labels
	// female:
	vector<vector<double>> female_labels_count_per_age; ///< for each age_bin, histogram of outcome labels

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
 * Parameter object for calc_roc_measures functions. this object
 * stores the working point, and other parameters for the roc measurments
 * bootstrap calculations.
 */
class ROC_Params : public Measurement_Params
{
public:
	vector<float> working_point_FPR;		///< The False Positive rate working point definition
	vector<int> working_point_TOPN;			///< The Top N working points
	vector<float> working_point_SENS;		///< The True Positive rate working point definition
	vector<float> working_point_PR;			///< The Positive rate working point definition
	vector<float> working_point_auc;		///< The partial auc working points definition
	vector<float> working_point_Score;		///< The Scores workin point definition
	bool use_score_working_points;			///< If true will calculate all roc measurements based on scores working points
	float max_diff_working_point;			///< The maximal diff in calculated working point to requested working point to drop
	int score_bins;							///< score bin count for speed up calculation. 0 means no binning
	int score_min_samples;					///< score bin min sample count for speed up calculation. 0 means no limit
	float score_resolution;					///< score resultion to contorl bining for speed up calculation. 0 means no binning resulotion
	bool fix_label_to_binary;				///< If True will change label value to be binary 0,1 (default is True)
	int min_score_quants_to_force_score_wp; ///< The minimal count of unique score to force fetching scores only by score cutoffs
	Incident_Stats inc_stats;				///< the incedince data if provided for general population. look for Incident_Stats for more info

	/// <summary>
	/// Default Ctor
	/// </summary>
	ROC_Params()
	{
		max_diff_working_point = (float)0.05;
		use_score_working_points = false;
		working_point_FPR = {(float)0.1, 1, 5, 10, 20, 30, 40, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95};
		score_bins = 0;
		score_resolution = 0;
		incidence_fix = 0;
		score_min_samples = 0;
		fix_label_to_binary = true;
		show_warns = true;
		min_score_quants_to_force_score_wp = 10;
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
	int init(map<string, string> &map);

	double incidence_fix; ///< The final incidence calculation on the cohort (will be calcuated)
	ADD_SERIALIZATION_FUNCS(working_point_FPR, working_point_SENS, working_point_PR, working_point_TOPN, working_point_Score, working_point_auc, use_score_working_points,
							max_diff_working_point, score_bins, score_resolution, score_min_samples, fix_label_to_binary, show_warns, inc_stats, min_score_quants_to_force_score_wp)
};

/**
 * Parameter object for Multiclass measure functions
 */
class Multiclass_Params : public Measurement_Params
{
public:
	vector<int> top_n;				   ///< when looking on top predictions, this is the maximal index
	int n_categ;					   ///< Number of categories
	vector<float> dist_weights;		   ///< Vector of weights - for index i : dist_weights[i] = 1/sum(dist(i,k))
	vector<vector<float>> dist_matrix; /// dist(i,j)
	string dist_name = "JACCARD";
	string dist_file;
	bool do_class_auc = false;

	Multiclass_Params()
	{
		top_n = {1, 5};
		n_categ = 1;
	}
	Multiclass_Params(const string &init_string);

	int init(map<string, string> &map);
	void read_dist_matrix_from_file(const string &fileName);

	ADD_CLASS_NAME(Multiclass_Params)
	ADD_SERIALIZATION_FUNCS(top_n, n_categ, dist_weights, dist_file, dist_matrix, dist_name, do_class_auc)
};

/**
 * Parameter object for Regression measure functions
 */
class Regression_Params : public Measurement_Params
{
public:
	bool do_logloss = false;
	double epsilon = 1e-5;
	vector<float> coverage_quantile_percentages;

	/// <summary>
	/// Initializing each parameter from string in format: "parameter_name=value;...". \n
	/// for vectors values use "," between numbers
	/// @snippet bootstrap.cpp Regression_Params::init
	/// </summary>
	int init(map<string, string> &mapper);

	ADD_CLASS_NAME(Regression_Params)
	ADD_SERIALIZATION_FUNCS(do_logloss, coverage_quantile_percentages, epsilon)
};

#pragma region Cohort Functions
/// <summary>
/// A function to filter samples based on single Filter_Param object. it's a FilterCohortFunc signature
/// </summary>
bool filter_range_param(const map<string, vector<float>> &record_info, int index, void *cohort_params); // on single param
/// <summary>
/// A function to filter samples based on multipal Filter_Param objects - in a vector with and condition
/// between each parameter range. it's a FilterCohortFunc signature
/// </summary>
bool filter_range_params(const map<string, vector<float>> &record_info, int index, void *cohort_params); // on vector of params
#pragma endregion

/**
 * Parameter object for filter_params functions
 */
class Filter_Param : public SerializableObject
{ // for example Age and range for filter
public:
	string param_name; ///< The parameter name for the filtering
	float min_range;   ///< the minimal range for the parameter
	float max_range;   ///< the maximal range for the parameter

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
	int init(map<string, string> &map);

	/// <summary>
	/// default Ctor
	/// </summary>
	Filter_Param() { param_name = ""; }

	ADD_SERIALIZATION_FUNCS(param_name, min_range, max_range)
};

struct ROC_And_Filter_Params : public Measurement_Params
{
	ROC_Params *roc_params;
	vector<Filter_Param> *filter;
};

// Infra
/// Function which recieves Lazy_Iterator and the thread num for iterating the predictions and labels.
///  it also recieves function_params which are additional arguments for the function (can be working points
///  defintions for example)
typedef map<string, float> (*MeasurementFunctions)(Lazy_Iterator *iterator, int thread_num, Measurement_Params *function_params);
/// Function which recieves map from feature name to vector of all samples value, sample index and cohort
///  definition params and return true\false if to include the sample in the cohort.
typedef bool (*FilterCohortFunc)(const map<string, vector<float>> &record_info, int index, void *cohort_params);
/// a function to process and maniplulate function params based on the given cohort - for example sotring
/// incedince information for the cohort
typedef void (*ProcessMeasurementParamFunc)(const map<string, vector<float>> &additional_info, const vector<float> &y, const vector<int> &pids, Measurement_Params *function_params,
											const vector<int> &filtered_indexes, const vector<float> &y_full, const vector<int> &pids_full);
/// a funtion to preprocess the prediction scores (binning for example to speed up bootstrap).
/// the function manipulate preds based on function_params
typedef void (*PreprocessScoresFunc)(vector<float> &preds, Measurement_Params *function_params);

#pragma region Process Measurement Param Functions
/// <summary>
/// a function to calculate the incidence in each cohort - preprocessing of function_params
/// and storing the incidence inside of it.
/// </summary>
void fix_cohort_sample_incidence(const map<string, vector<float>> &additional_info,
								 const vector<float> &y, const vector<int> &pids, Measurement_Params *function_params,
								 const vector<int> &filtered_indexes, const vector<float> &y_full, const vector<int> &pids_full);

/// <summary>
/// a function to calculate the incidence in each cohort - preprocessing of function_params
/// and storing the incidence inside of it. The old has same implementation as old bootstrap
/// only averaging incidence over the controls in the sample based on incidence in each group(age+gender)
/// </summary>
void fix_cohort_sample_incidence_old(const map<string, vector<float>> &additional_info,
									 const vector<float> &y, const vector<int> &pids, Measurement_Params *function_params,
									 const vector<int> &filtered_indexes, const vector<float> &y_full, const vector<int> &pids_full);
#pragma endregion

#pragma region Process Scores Functions
/// <summary>
/// Binning function of scores based on ROC_Params. look at score_bins,score_resolution
/// </summary>
void preprocess_bin_scores(vector<float> &preds, Measurement_Params *function_params);
#pragma endregion

/// <summary>
/// Format out measurement
/// </summary>
inline string format_working_point(const string &init_str, float wp, bool perc = true)
{
	char res[100];
	if (perc)
		wp *= 100;
	snprintf(res, sizeof(res), "%s_%06.3f", init_str.c_str(), wp);
	return string(res);
}

inline string format_working_point_topn(const string &init_str, int wp, bool perc = true)
{
	char res[100];
	snprintf(res, sizeof(res), "%s_%d", init_str.c_str(), wp);
	return string(res);
}

/// <summary>
/// to run bootstrap on single cohort
/// </summary>
map<string, float> booststrap_analyze_cohort(const vector<float> &preds, const vector<int> &preds_order, const vector<float> &y,
											 const vector<int> &pids, float sample_ratio, int sample_per_pid, int loopCnt,
											 const vector<MeasurementFunctions> &meas_functions, const vector<Measurement_Params *> &function_params,
											 ProcessMeasurementParamFunc process_measurments_params,
											 const map<string, vector<float>> &additional_info, const vector<float> &y_full,
											 const vector<int> &pids_full, const vector<float> *weights, const vector<int> &filter_indexes, FilterCohortFunc cohort_def,
											 void *cohort_params, int &warn_cnt, const string &cohort_name, int seed = 0);

/// @brief The main bootstrap function to run all bootstrap process with all the arguments
/// @param preds vector of predictions
/// @param y labels
/// @param pids patient ids [used for random draws]
/// @param additional_info Dictionary used to keep values of features [for filtering]. Key is the feature name. Val contains vector of feature values.
/// @param filter_cohort Dictionary, where key is name of cohort, value is a function which performs the filtering [keeps only entries that belong to this cohort]

/// @param meas_functions Vector of metrics to calculate per one bootstrap (?) experiment
/// @param cohort_params Key: name of a cohort, Value: additional parameters which are passed to functions from filter_cohort.values
/// @param function_params Configuration parameters passed to "meas_functions" (like 2 for F2 metric).
/// @param process_measurments_params Function to process the function_params before running on each cohort (helps to calc incidence for example)
/// @param preprocess_scores A function to preprocess all scores - for example binning the scores (can sometimes speedup metrics calculation)
/// @param preprocess_scores_params Additional parameters for the preprocess_scores function
/// @param sample_ratio A number in range (0,1] for subsampling the samples [in order to speed-up the bootstrap]
/// @param sample_per_pid How many samples to sample on each id [sampling with replacement]
/// @param loopCnt How many bootstrap experiments(?) to do
/// @param seed The random seed. If 0 will use random_device to create random seed
/// @param binary_outcome A flag to indicate whether the labels are binary (used to validate the input labels)
/// <returns>
/// Returns a map from each cohort name to the measurments results. each measurments results
/// is also a map from each measurement name to it's value
/// </returns>
map<string, map<string, float>> booststrap_analyze(
	const vector<float> &preds,
	const vector<int> &preds_order,
	const vector<float> &y,
	const vector<float> *weights,
	const vector<int> &pids,
	const map<string, vector<float>> &additional_info,
	const map<string, FilterCohortFunc> &filter_cohort,
	const vector<MeasurementFunctions> &meas_functions = {calc_roc_measures_with_inc},
	const map<string, void *> *cohort_params = NULL,
	const vector<Measurement_Params *> *function_params = NULL,
	ProcessMeasurementParamFunc process_measurments_params = NULL,
	PreprocessScoresFunc preprocess_scores = NULL,
	Measurement_Params *preprocess_scores_params = NULL,
	float sample_ratio = (float)1.0,
	int sample_per_pid = 1,
	int loopCnt = 500,
	int seed = 0,
	bool binary_outcome = true);

/// <summary>
/// @param pids the pids vector
/// @param additional_info the data vector for filtering
/// @param filter_cohort The cohorts definition - the filtering function
/// @param cohort_params Additional parameters for the filtering cohort function
/// @param sample_ratio A number in range (0,1] for subsampling the samples
/// @param sample_per_pid How many samples to sample on each id
/// @param seed The random seed. If 0 will use random_device to create random seed
/// @param indexes the selected indexes results for the bootstrap
/// <returns>
/// Returns indexes vector
/// </returns>
/// </summary>
void prepare_for_bootstrap(const vector<int> &pids,
						   const map<string, vector<float>> &additional_info, FilterCohortFunc &filter_cohort, void *cohort_params, float sample_ratio, int sample_per_pid, int seed, vector<int> &indexes);

/// <summary>
/// Will output the bootstrap results into file in TAB delimeted format. each line is cohort and the
/// The columns are the measurements
/// </summary>
void write_bootstrap_results(const string &file_name, const map<string, map<string, float>> &all_cohorts_measurments, const string &run_id = "");
/// <summary>
/// Will read the bootstrap results from file in TAB delimeted format. each line is cohort and the
/// The columns are the measurements
/// </summary>
void read_bootstrap_results(const string &file_name, map<string, map<string, float>> &all_cohorts_measurments);

/// <summary>
/// Will output the bootstrap results into file with the new format with columns:
/// "Cohort$Measure_Name", "Value"
/// </summary>
void write_pivot_bootstrap_results(const string &file_name, const map<string, map<string, float>> &all_cohorts_measurments, const string &run_id = "");
/// <summary>
/// Will read the bootstrap results into file with the new format with columns:
/// "Cohort$Measure_Name", "Value"
/// <//summary>
void read_pivot_bootstrap_results(const string &file_name, map<string, map<string, float>> &all_cohorts_measurments);

MEDSERIALIZE_SUPPORT(Incident_Stats)
MEDSERIALIZE_SUPPORT(ROC_Params)
MEDSERIALIZE_SUPPORT(Filter_Param)
MEDSERIALIZE_SUPPORT(Regression_Params)
MEDSERIALIZE_SUPPORT(Multiclass_Params)

#endif // !__BOOTSTRAP_ANALYSIS_H__
