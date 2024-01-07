#ifndef __MEDBOOTSTRAP_ANALYSIS_H__
#define __MEDBOOTSTRAP_ANALYSIS_H__
#include <unordered_map>
#include <MedStat/MedStat/bootstrap.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedUtils/MedUtils/MedRegistry.h>

/**
* A class holder for optional arguments for bootstrap to use MedRegistry
* If provided will calculated more accurate the incidence rate on each filter cohort
* based on kaplan-meier(default option) and will use registry to censor samples when
* using with sim_time_window mode
*/
/// @enum
/// Measurment types
enum class MeasurmentFunctionType {
	calc_npos_nneg = 0,
	calc_only_auc = 1,
	calc_roc_measures_with_inc = 2, 
	calc_multi_class = 3,
	calc_kandel_tau = 4,
	calc_harrell_c_statistic = 5,
	calc_regression = 6
};
class with_registry_args {
public:
	MedRegistry *registry; ///< the registry of records
	MedRegistry *registry_censor; ///< the registry censor of records
	MedSamplingYearly *sampler; ///< the sampler for calculating incidence for example yearly from year to year
	bool do_kaplan_meir; ///< If true will do kaplan meier
	string rep_path; ///< repository path
	string json_model; ///< The json_model path to create matrix to calc incidence and filter cohort
	LabelParams labeling_params;
	with_registry_args() {
		do_kaplan_meir = true;
		registry_censor = NULL;
		registry = NULL;
		sampler = NULL;
	}
};

/**
* Bootstrap wrapper for Medila Infrastructure objects, simplify the parameters
* and the input, output process. \n
* for more control and lower level interface please refer to bootstrap.h
*/
class MedBootstrap : public SerializableObject {
public:
	ROC_Params roc_Params; ///< Controling the roc parameters: sensitivity, specificity...
	Regression_Params regression_params; ///< params for regerssion
	Multiclass_Params multiclass_params; ///< Controling the multi class parameters: top n...
	map<string, vector<Filter_Param>> filter_cohort; ///< the cohorts definitions. name to parameters range to intersect
	map<string, FilterCohortFunc> additional_cohorts; ///< not Serializable! additional cohorts given by function
	float sample_ratio; ///<the sample ratio of the patients out of all patients in each bootstrap
	int sample_per_pid; ///<how many samples to take for each patients. 0 - means no sampling take all sample for patient
	bool sample_patient_label; ///<if true will treat patient+label as the "id" for the sampling
	int sample_seed; ///<if 0 will use random_device
	int loopCnt; ///<the bootstrap count
	bool is_binary_outcome; ///< only used for validating bootstrap input
	bool use_time_control_as_case; ///< if True will use time window condition for controls same as cases.
	///Time window simulation (in cohorts with Time-Window filtering) - instead of censoring cases out of time range
	///, treat them as controls
	bool simTimeWindow;
	float censor_time_factor;
	bool sort_preds_in_multicategory;
	size_t num_categories; ///< number of categories
	vector<pair<MeasurementFunctions, Measurement_Params *>> measurements_with_params;  ///<not Serializable! the measurements with the params

	/// <summary>
	/// parsing specific line. please refer to parse_cohort_file for full spec
	/// </summary>
	void parse_cohort_line(const string &line);

	/// <summary>
	/// A function which reads a single cohort definition from the command line and parses it. \n
	/// Please refer to parse_cohort_file for full spec.
	/// </summary>
	void get_cohort_from_arg(const string &single_cohort);


	/// <summary>
	/// a function which reads cohorts file and stores it in filter_cohort.
	/// The file format may be in 2 options:
	///   -# COHORT_NAME[TAB]PARAMETERS_DEF - cohort name is string representing cohort \n
	///      name. PARAMETER_DEF is in format: "PARAMETER_NAME:MIN_RANGE,MAX_RANGE;..." \n
	///      the format can repeat itself with ";" between each parameter. the cohort \n
	///      will consist of intersection between all parameters ranges with "and" condition. \n
	///      there is single tab betwwen the name and the defenition. \n
	///      Example Line: \n
	///      1 year back & age 40-80	Time-Window:0,365;Age:40,80 \n
	///      will create cohort called "1 year back & age 40-80" and will filter out records \n
	///      with (Time-Window>=0 and Time-Window<=365) and (Age>=40 and Age<=80) \n
	///   -# MULTI[TAB]PARAMETERS_DEF[TAB]...PARAMETERS_DEF[TAB] - this definition with \n
	///      line starting with MULTI keyword will create all the cartesain options for each \n
	///      parameter definition with the each parameter definition in the next TABs. \n
	///      PARAMETERS_DEF - is same as option 1 format. \n
	///      Example Line: \n
	///      MULTI Time-Window:0,30;Time-Window:30,180	Age:40,60;Age:60,80;Age:40,80	Gender:1,1;Gender:2,2 \n
	///      will create 2*3*2=12 cohorts for each Time-Window, Age, and Gender option \n
	/// </summary>
	void parse_cohort_file(const string &cohorts_path);

	/// <summary>
	/// defualt Ctor. look for ROC_Params defaults. cohorts consists of 1 cohort called "All" with
	/// not filtering
	/// </summary>
	MedBootstrap();

	/// <summary>
	/// Initialization string with format "parameter_name=value;..."
	/// each paramter_name is same as the class name field. filter_cohort is path to file
	/// roc_Params is the init string for ROC_PARAMS
	/// @snippet MedBootstrap.cpp MedBootstrap::init
	/// </summary>
	int init(map<string, string>& map);

	/// <summary>
	/// cleans the initiale "FTR_" from the feature names in MedFeatures created by the infra pipeline
	/// </summary>
	void clean_feature_name_prefix(map<string, vector<float>> &features);

	/// <summary>
	///prepares the required vectors for bootstrap from MedFeatures &features
	/// </summary>
	/// <returns>
	/// updates - preds, y, pids, final_additional_info with the information from MedFeatures.
	/// if splits_inds is provided (and not NULL) it will fill a mapping from split_index to the 
	/// indexes in the samples vector correspond to each split value
	/// </returns>
	void prepare_bootstrap(const MedFeatures &features, vector<float> &preds, vector<float> &y, vector<int> &pids,
		map<string, vector<float>> &final_additional_info, vector<int> &preds_order, unordered_map<int, vector<int>> *splits_inds = NULL);
	/// <summary>
	///prepares the required vectors for bootstrap from samples, additional_info
	/// </summary>
	/// <returns>
	/// updates - preds, y, pids, final_additional_info with the information from samples, additional_info.
	/// if splits_inds is provided (and not NULL) it will fill a mapping from split_index to the 
	/// indexes in the samples vector correspond to each split value
	/// </returns>
	void prepare_bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info, vector<float> &preds, vector<float> &y, vector<int> &pids, vector<int> &preds_order,
		unordered_map<int, vector<int>> *splits_inds = NULL);
	/// <summary>
	/// Will run the bootstraping process on all cohorts and measurements.
	/// MedFeatures need to contains also the information for the cohorts defenitions.
	/// for example: if there is Age:40-80, MedFeatures should contain Age Feature
	/// </summary>
	/// <returns>
	/// the bootstrap results in map from cohort_name to all cohort measurements(a map).
	/// Each measurement is key,value in the map from measurement name to it's value
	/// if splits_inds is not NULL and mapping from each split value to it's coresponding
	/// indexes in the samples are provided - it will return also results for each split
	/// the higest level in the map is the split value
	/// </returns>
	map<string, map<string, float>> bootstrap(const MedFeatures &features, map<int, map<string, map<string, float>>> *results_per_split = NULL, with_registry_args *registry_args = NULL);
	/// <summary>
	/// Will run the bootstraping process on all cohorts and measurements.
	/// the input is samples, additional_info. additional_info is provided for filtering 
	/// and creating the cohorts. for example - Age:40-80 and Males
	/// </summary>
	/// <returns>
	/// the bootstrap results in map from cohort_name to all cohort measurements(a map).
	/// Each measurement is key,value in the map from measurement name to it's value
	/// if splits_inds is not NULL and mapping from each split value to it's coresponding
	/// indexes in the samples are provided - it will return also results for each split
	/// the higest level in the map is the split value
	/// </returns>
	map<string, map<string, float>> bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info, map<int, map<string, map<string, float>>> *results_per_split = NULL, with_registry_args *registry_args = NULL);
	/// <summary>
	/// Will run the bootstraping process on all cohorts and measurements.
	/// the input is samples, and rep_path. The rep_path is path to the repository which 
	/// adds Age,Gender signals for creating the cohorts definitions. it's simple overload
	/// for convention
	/// </summary>
	/// <returns>
	/// the bootstrap results in map from cohort_name to all cohort measurements(a map).
	/// Each measurement is key,value in the map from measurement name to it's value
	/// if splits_inds is not NULL and mapping from each split value to it's coresponding
	/// indexes in the samples are provided - it will return also results for each split
	/// the higest level in the map is the split value
	/// </returns>
	map<string, map<string, float>> bootstrap(MedSamples &samples, const string &rep_path, map<int, map<string, map<string, float>>> *results_per_split = NULL, with_registry_args *registry_args = NULL);
	/// <summary>
	/// Will run the bootstraping process on all cohorts and measurements.
	/// the input is samples, and rep. The rep is the repository which 
	/// adds Age,Gender signals for creating the cohorts definitions. it's simple overload
	/// for convention
	/// </summary>
	/// <returns>
	/// the bootstrap results in map from cohort_name to all cohort measurements(a map).
	/// Each measurement is key,value in the map from measurement name to it's value
	/// if splits_inds is not NULL and mapping from each split value to it's coresponding
	/// indexes in the samples are provided - it will return also results for each split
	/// the higest level in the map is the split value
	/// </returns>
	map<string, map<string, float>> bootstrap(MedSamples &samples, MedPidRepository &rep, map<int, map<string, map<string, float>>> *results_per_split = NULL, with_registry_args *registry_args = NULL);

	/// <summary>
	/// censors samples from samples based on time_range provided in pid_censor_dates.
	/// the format is map from pid to max_date the after that date the sample is filtered.
	/// </summary>
	/// <returns>
	/// update samples - changes outcomeDate for controls to censor date.
	/// </returns>
	void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedSamples &samples);
	/// <summary>
	/// censors samples from features based on time_range provided in pid_censor_dates.
	/// the format is map from pid to max_date the after that date the sample is filtered.
	/// </summary>
	/// <returns>
	/// update features - changes outcomeDate for controls to censor date.
	/// </returns>
	void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedFeatures &features);
	/// <summary>
	/// censors samples from features based on time_range provided in pids,censor_dates.
	/// pids and censor_dates are same sizes.
	/// for each pid and the coresponding date in censor_dates, filtering pid's samples after
	/// that date.
	/// </summary>
	/// <returns>
	/// update features - changes outcomeDate for controls to censor date.
	/// </returns>
	void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedFeatures &features);
	/// <summary>
	/// censors samples from samples based on time_range provided in pids,censor_dates.
	/// pids and censor_dates are same sizes.
	/// for each pid and the coresponding date in censor_dates, filtering pid's samples after
	/// that date.
	/// </summary>
	/// <returns>
	/// update samples - changes outcomeDate for controls to censor date.
	/// </returns>
	void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedSamples &samples);

	/// <summary>
	/// changing the samples to be auto-simulations - taking max score in the time window for each
	/// pid
	/// </summary>
	/// <returns>
	/// updates new_samples from samples
	/// </returns>
	void change_sample_autosim(MedSamples &samples, int min_time, int max_time, MedSamples &new_samples);
	/// <summary>
	/// changing the samples to be auto-simulations - taking max score in the time window for each
	/// pid
	/// </summary>
	/// <returns>
	/// updates new_features from features
	/// </returns>
	void change_sample_autosim(MedFeatures &features, int min_time, int max_time, MedFeatures &new_features);
	/// <summary>
	/// convert measurement function name to type 
	/// </summary>
	/// <returns>
	/// MeasurmentFunctionType 
	/// </returns>
	MeasurmentFunctionType measurement_function_name_to_type(const string& measurement_function_name);

	static unordered_map<string, MeasurmentFunctionType> measurement_function_name_map;

	/// <summary>
	/// commit bootstrap cohort filter on a given matrix
	/// @param features - matrix
	/// @param bt_cohort - a single line cohort (no support for MULTI) without the cohort name. 
	/// only the filter definition. no tabs in the string.
	/// </summary>
	/// <returns>
	/// filter rows from features by cohort definition
	/// </returns>
	static void filter_bootstrap_cohort(MedFeatures &features, const string &bt_cohort);
	/// <summary>
	/// commit bootstrap cohort filter on a given samples
	/// @param bt_repository - repository that was initialized for applying the bt_filters model
	/// to generate matrix
	/// @param bt_filters - the model to generate matrix for filtering the bootstrap cohort
	/// @param curr_samples - the samples to filter
	/// @param bt_cohort - a single line cohort (no support for MULTI) without the cohort name. 
	/// only the filter definition. no tabs in the string.
	/// </summary>
	/// <returns>
	/// filter samples from curr_samples by cohort definition
	/// </returns>
	static void filter_bootstrap_cohort(MedPidRepository &bt_repository, MedModel &bt_filters,
		MedSamples &curr_samples, const string &bt_cohort);
	/// <summary>
	/// commit bootstrap cohort filter on a given samples
	/// @param rep - repository path
	/// @param bt_json - the json model to generate matrix for filtering the bootstrap cohort
	/// Automatically Age,Gender are added
	/// @param curr_samples - the samples to filter
	/// @param bt_cohort - a single line cohort (no support for MULTI) without the cohort name. 
	/// only the filter definition. no tabs in the string.
	/// </summary>
	/// <returns>
	/// filter samples from curr_samples by cohort definition
	/// </returns>
	static void filter_bootstrap_cohort(const string &rep, const string &bt_json,
		MedSamples &curr_samples, const string &bt_cohort);

	ADD_CLASS_NAME(MedBootstrap)
	ADD_SERIALIZATION_FUNCS(sample_ratio, sample_per_pid, sample_patient_label, sample_seed, loopCnt, roc_Params, filter_cohort, simTimeWindow, multiclass_params, censor_time_factor)

private:
	map<string, map<string, float>> bootstrap_base(const vector<float> &preds, const vector<int> &preds_order,  const vector<float> &y, const vector<int> &pids,
		const vector<float> *weights, const map<string, vector<float>> &additional_info);
	map<string, map<string, float>> bootstrap_using_registry(const MedFeatures &features_mat,
		const with_registry_args& args, map<int, map<string, map<string, float>>> *results_per_split = NULL);
	void add_splits_results(const vector<float> &preds, vector<int> &preds_order, const vector<float> &y,
		const vector<int> &pids, const vector<float> *weights, const map<string, vector<float>> &data,
		const unordered_map<int, vector<int>> &splits_inds,
		map<int, map<string, map<string, float>>> &results_per_split);
	bool use_time_window();
	void add_filter_cohorts(const map<string, vector<pair<float, float>>> &parameters_ranges);
	void add_filter_cohorts(const vector<vector<Filter_Param>> &parameters_ranges);
	void sort_index_only(const vector<float> &vec, std::vector<int>::iterator ind_start, std::vector<int>::iterator ind_end);
};

/**
* A wrapper class which contains the MedBootstrap object and the results for later
* quering the scores for preformance details
*/
class MedBootstrapResult : public SerializableObject {
public:
	MedBootstrap bootstrap_params; ///<The boostrap parameters
	///The bootstrap results. map from cohort_name to all cohort measurements(a map). 
	///each measurement is key,value in the map from measurement name to it's value
	map<string, map<string, float>> bootstrap_results;

	/// <summary>
	/// run the bootstrap - look at MedBootstrap.bootstrap documentition and stores
	/// the results in bootstrap_results
	/// </summary>
	void bootstrap(MedFeatures &features, map<int, map<string, map<string, float>>> *results_per_split = NULL, with_registry_args *registry_args = NULL);
	/// <summary>
	/// run the bootstrap - look at MedBootstrap.bootstrap documentition and stores
	/// the results in bootstrap_results
	/// </summary>
	void bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info, map<int, map<string, map<string, float>>> *results_per_split = NULL, with_registry_args *registry_args = NULL);
	/// <summary>
	/// run the bootstrap - look at MedBootstrap.bootstrap documentition and stores
	/// the results in bootstrap_results
	/// </summary>
	void bootstrap(MedSamples &samples, const string &rep_path, map<int, map<string, map<string, float>>> *results_per_split = NULL, with_registry_args *registry_args = NULL);

	/// <summary>
	/// searches for the sensitivty(sens) and positive rate(pr) in the sepcific bootstrap_cohort results
	/// which is map from measure name to it's value. the bootstrap_cohort results is in binary/categorial
	/// (which working points are defined @SCORE)
	/// </summary>
	/// <returns>
	/// updates sens_points, pr_points vectors for SCORE defined working points bootstrap
	/// </returns>
	void find_working_points(const map<string, float> &bootstrap_cohort,
		vector<float> &sens_points, vector<float> &pr_points);

	/// <summary>
	/// searches for sepcific score value the corresonding measurments in that working
	/// point based on the bootstrap result.
	/// you need to run or load object with bootstrap_results which is not empty.
	/// @param score the score that defines the working point
	/// @param string_cohort the cohort name to search for the measurments with
	/// @param max_search_range max interval for searching working points measurements
	/// in the bootstrap. if the working point in the bootstrap results are too far away
	/// from the score working point it will not return any result
	/// </summary>
	/// <returns>
	/// updates score_measurements with the bootstrap corresponding measurements
	/// </returns>
	void explore_score(float score, map<string, float> &score_measurements,
		const string &string_cohort = "All", float max_search_range = 0.1);

	/// <summary>
	/// writes the results to file with TAB delimeted manner. you can also
	/// pivot the format
	/// </summary>
	void write_results_to_text_file(const string &path, bool pivot_format = true, const string& run_id = "");
	/// <summary>
	/// reads the results from file with TAB delimeted manner. you can also
	/// read pivot format
	/// </summary>
	void read_results_to_text_file(const string &path, bool pivot_format = true);

	ADD_CLASS_NAME(MedBootstrapResult)
	ADD_SERIALIZATION_FUNCS(bootstrap_params, bootstrap_results)
private:
	bool find_in_range(const vector<float> &vec, float search, float th);
	void explore_measure(const string &measure_name, float value, map<string, float> &score_measurements,
		const string &string_cohort = "All", float max_search_range = 0.1);
};

namespace medial {
	namespace process {
		/// \brief converts cohort to sim time window- changes cases that in long term will only
		/// turn in to cases into controls
		void make_sim_time_window(const string &cohort_name, const vector<Filter_Param> &filter_p,
			const vector<float> &y, const map<string, vector<float>> &additional_info,
			vector<float> &y_changed, map<string, vector<float>> &cp_info,
			map<string, FilterCohortFunc> &cohorts_t, map<string, void *> &cohort_params_t, float censor_time_factor = 2);
	}
}

MEDSERIALIZE_SUPPORT(MedBootstrap)
MEDSERIALIZE_SUPPORT(MedBootstrapResult)

#endif