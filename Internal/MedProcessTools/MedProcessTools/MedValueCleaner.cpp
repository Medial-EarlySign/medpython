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
	float upper = values[(int)(values.size() * (1.0 - params.quantile))];
	float lower = values[(int)(values.size() * params.quantile)];

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

	// Take Log if required
	if (params.take_log) {
		for (unsigned int i = 0; i < values.size(); i++) {
			if (values[i] <= 0)
				values[i] = params.missing_value;
			else if (values[i] != params.missing_value)
				values[i] = log(values[i]);
		}
	}

	bool need_to_clean = true;
	float mean, sd, min, max;

	vector<float> wgts(values.size(), 1.0);

	while (need_to_clean) {
		need_to_clean = false;

		int n = get_moments(values, wgts, params.missing_value, mean, sd);
		if (n == 0) {
			MERR("Cannot learn cleaning parameters from an empty vector\n");
			return -1;
		}

		max = mean + params.trimming_sd_num * sd;
		min = mean - params.trimming_sd_num * sd;

		// Clean
		need_to_clean = false;
		for (unsigned int i = 0; i < values.size(); i++) {
			if (values[i] != params.missing_value) {
				if (values[i] > max) {
					need_to_clean = true;
					values[i] = max;
				}
				else if (values[i] < min) {
					need_to_clean = true;
					values[i] = min;
				}
			}
		}
	}


	trimMax = max; if (params.take_log) trimMax = exp(trimMax);
	trimMin = min; if (params.take_log) trimMin = exp(trimMin);
	removeMax = mean + params.removing_sd_num * sd; if (params.take_log) removeMax = exp(removeMax);
	removeMin = mean - params.removing_sd_num * sd; if (params.take_log) removeMin = exp(removeMin);
	nbrsMax = mean + params.nbrs_sd_num * sd; if (params.take_log) nbrsMax = exp(nbrsMax);
	nbrsMin = mean - params.nbrs_sd_num * sd; if (params.take_log) nbrsMin = exp(nbrsMin);

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

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "type") params.type = get_cleaner_type(entry.second);
		else if (field == "take_log") params.take_log = stoi(entry.second);
		else if (field == "missing_value") params.missing_value = stof(entry.second);
		else if (field == "trimming_sd_num") params.trimming_sd_num = stof(entry.second);
		else if (field == "removing_sd_num") params.removing_sd_num = stof(entry.second);
		else if (field == "nbrs_sd_num") params.nbrs_sd_num = stof(entry.second);
		else if (field == "quantile") params.quantile = stof(entry.second);
		else if (field == "trimming_quantile_factor") params.trimming_quantile_factor = stof(entry.second);
		else if (field == "trimming_quantile_factor") params.trimming_quantile_factor = stof(entry.second);
		else if (field == "nbrs_quantile_factor") params.nbrs_quantile_factor = stof(entry.second);
		else if (field == "doTrim") params.doTrim = (stoi(entry.second) != 0);
		else if (field == "doRemove") params.doRemove = (stoi(entry.second) != 0);
		// next are in ignore ... used in level above
		else if (field != "signal" && field != "time_unit" && field != "time_channel" && field != "fp_type" &&
				 field != "val_channel" && field != "nbr_time_unit" && field != "nbr_time_width" && field != "rp_type")
			MLOG("Unknonw parameter \'%s\' for MedValueCleaner\n", field.c_str());

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


