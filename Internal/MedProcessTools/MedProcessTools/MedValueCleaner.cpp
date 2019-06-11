#include "InfraMed/InfraMed/InfraMed.h"
#include "Logger/Logger/Logger.h"
#include "MedValueCleaner.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"

#define LOCAL_SECTION LOG_VALCLNR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// MedValueCleaner
//=======================================================================================
// Quantile cleaning
int MedValueCleaner::get_quantile_min_max(vector<float>& values) {

	if (values.size() == 0)
		MTHROW_AND_ERR("Trying to get quantiles of an empty array\n");
	if (params.take_log) {
		for (unsigned int i = 0; i < values.size(); i++) {
			if (values[i] <= 0)
				values[i] = params.missing_value;
			else if (values[i] != params.missing_value)
				values[i] = log(values[i]);
		}
	}
	sort(values.begin(), values.end());

	float median = values[(int)(values.size() * 0.5)];
	int access_ind = (int)(values.size() * (1.0 - params.quantile));
	if (access_ind >= values.size()) //for example may happen if params.quantile == 0
		access_ind = (int)values.size() - 1;
	float upper = values[access_ind];
	int lower_access_ind = (int)(values.size() * params.quantile);
	float lower = values[lower_access_ind];
	num_samples_after_cleaning = access_ind - lower_access_ind;

	if (params.take_log) {
		if (median <= 0.0 || lower <= 0.0 || upper <= 0.0) {
			MERR("Cannot take log of non-positive quantile\n");
			return -1;
		}
		median = log(median);
		lower = log(lower);
		upper = log(upper);
	}

	trimMax = median + (upper - median)*params.trimming_quantile_factor; if (params.take_log) trimMax = exp(trimMax);
	trimMin = median + (lower - median)*params.trimming_quantile_factor; if (params.take_log) trimMin = exp(trimMin);
	removeMax = median + (upper - median)*params.removing_quantile_factor; if (params.take_log) removeMax = exp(removeMax);
	removeMin = median + (lower - median)*params.removing_quantile_factor; if (params.take_log) removeMin = exp(removeMin);
	nbrsMax = median + (upper - median)*params.nbrs_quantile_factor; if (params.take_log) nbrsMax = exp(nbrsMax);
	nbrsMin = median + (lower - median)*params.nbrs_quantile_factor; if (params.take_log) nbrsMin = exp(nbrsMin);

	return 0;
};

//.......................................................................................
// Iterative cleaning
int MedValueCleaner::get_iterative_min_max(vector<float>& values) {

	float max_range = params.range_max;
	float min_range = params.range_min;
	float trim_max_range = params.trim_range_max;
	float trim_min_range = params.trim_range_min;

	// Take Log if required
	if (params.take_log) {
		for (unsigned int i = 0; i < values.size(); i++) {
			if (values[i] <= 0)
				values[i] = params.missing_value;
			else if (values[i] != params.missing_value)
				values[i] = log(values[i]);
		}
		if (max_range <= 0) max_range = 1;
		if (min_range <= 0) min_range = (float)0.01;
		max_range = log(max_range);
		min_range = log(min_range);

		if (trim_max_range <= 0) trim_max_range = 1;
		if (trim_min_range <= 0) trim_min_range = (float)0.01;
		trim_max_range = log(trim_max_range);
		trim_min_range = log(trim_min_range);

	}

	bool need_to_clean = true;
	float mean, sd, vmin, vmax;

	while (need_to_clean) {
		need_to_clean = false;
		medial::stats::get_mean_and_std(values, params.missing_value, num_samples_after_cleaning, mean, sd);
		if (num_samples_after_cleaning == 0) {
			MWARN("EMPTY_VECTOR:: learning cleaning parameters from an empty vector\n");
			trimMax = 0;
			trimMin = 0;
			removeMax = 0;
			removeMin = 0;
			nbrsMax = 0;
			nbrsMin = 0;

			return 0;
		}

		vmax = min(min(mean + params.trimming_sd_num * sd, max_range), trim_max_range);
		vmin = max(max(mean - params.trimming_sd_num * sd, min_range), trim_min_range);

		// Clean
		need_to_clean = false;
		for (unsigned int i = 0; i < values.size(); i++) {
			if (values[i] != params.missing_value) {
				if (values[i] > vmax) {
					need_to_clean = true;
					values[i] = vmax;
				}
				else if (values[i] < vmin) {
					need_to_clean = true;
					values[i] = vmin;
				}
			}
		}
	}

	trimMax = vmax; if (params.take_log) trimMax = exp(trimMax);
	trimMin = vmin; if (params.take_log) trimMin = exp(trimMin);
	removeMax = mean + params.removing_sd_num * sd; if (params.take_log) removeMax = exp(removeMax);
	removeMin = mean - params.removing_sd_num * sd; if (params.take_log) removeMin = exp(removeMin);
	nbrsMax = mean + params.nbrs_sd_num * sd; if (params.take_log) nbrsMax = exp(nbrsMax);
	nbrsMin = mean - params.nbrs_sd_num * sd; if (params.take_log) nbrsMin = exp(nbrsMin);

	trimMax = min(trimMax, params.trim_range_max);
	trimMin = max(trimMin, params.trim_range_min);
	removeMax = min(removeMax, params.range_max);
	removeMin = max(removeMin, params.range_min);
	nbrsMax = min(nbrsMax, params.range_max);
	nbrsMin = max(nbrsMin, params.range_min);



	return 0;
}

// Init
//.......................................................................................
int MedValueCleaner::init(void *_in_params)
{

	ValueCleanerParams *in_params = (ValueCleanerParams *)_in_params;

	params.type = in_params->type;
	params.take_log = in_params->take_log;
	params.missing_value = in_params->missing_value;
	params.trimming_sd_num = in_params->trimming_sd_num;
	params.removing_sd_num = in_params->removing_sd_num;
	params.quantile = in_params->quantile;
	params.trimming_quantile_factor = in_params->trimming_quantile_factor;
	params.removing_quantile_factor = in_params->removing_quantile_factor;
	params.doTrim = in_params->doTrim;
	params.doRemove = in_params->doRemove;


	return 0;
}

//..............................................................................
int MedValueCleaner::init(map<string, string>& mapper) {

	unordered_set<string> remove_me = { "verbose_file" ,"fp_type", "rp_type", "unconditional", "signal", "time_channel", "val_channel", "nrem_attr", "ntrim_attr", "nrem_suff",
		"ntrim_suff", "time_unit", "nbr_time_unit", "nbr_time_width",  "tag", "conf_file", "clean_method","signals", "addRequiredSignals",
		"consideredRules", "print_summary", "print_summary_critical_cleaned" };


	for (auto entry : mapper) {
		string field = entry.first;
		//! [MedValueCleaner::init]
		if (field == "type") params.type = get_cleaner_type(entry.second);
		else if (field == "take_log") params.take_log = med_stoi(entry.second);
		else if (field == "missing_value") params.missing_value = med_stof(entry.second);
		else if (field == "trimming_sd_num") params.trimming_sd_num = med_stof(entry.second);
		else if (field == "removing_sd_num") params.removing_sd_num = med_stof(entry.second);
		else if (field == "nbrs_sd_num") params.nbrs_sd_num = med_stof(entry.second);
		else if (field == "quantile") params.quantile = med_stof(entry.second);
		else if (field == "trimming_quantile_factor") params.trimming_quantile_factor = med_stof(entry.second);
		else if (field == "removing_quantile_factor") params.removing_quantile_factor = med_stof(entry.second);
		else if (field == "nbrs_quantile_factor") params.nbrs_quantile_factor = med_stof(entry.second);
		else if (field == "doTrim") params.doTrim = (med_stoi(entry.second) != 0);
		else if (field == "doRemove") params.doRemove = (med_stoi(entry.second) != 0);
		else if (field == "range_min") params.range_min = med_stof(entry.second);
		else if (field == "range_max") params.range_max = med_stof(entry.second);
		else if (field == "trim_range_min") params.trim_range_min = med_stof(entry.second);
		else if (field == "trim_range_max") params.trim_range_max = med_stof(entry.second);
		else if (field == "max_samples") params.max_samples = med_stoi(entry.second);
		else if (remove_me.find(field) == remove_me.end()) MWARN("MedValueCleaner:: Warn Unknown param \"%s\"\n", field.c_str());
		//! [MedValueCleaner::init]

	}

	return 0;
}

//..............................................................................
// Get Type
ValueCleanerType MedValueCleaner::get_cleaner_type(string name) {

	boost::algorithm::to_lower(name);
	if (name == "iterative")
		return VAL_CLNR_ITERATIVE;
	else if (name == "quantile")
		return VAL_CLNR_QUANTILE;
	else
		return VAL_CLNR_LAST;
}


