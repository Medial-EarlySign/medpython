#include <unordered_set>
#include "MedBootstrap.h"
#include "Logger/Logger/Logger.h"

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL 

map<string, map<string, float>> booststrap_base(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
	const map<string, vector<float>> &additional_info,
	const map<string, vector<Filter_Param>> &filter_cohort,
	const ROC_Params &roc_params, float sample_ratio, int sample_per_pid, int loopCnt) {

	map<string, FilterCohortFunc> cohorts;
	vector<MeasurementFunctions> measures = { calc_npos_nneg, calc_roc_measures_full };
	map<string, void *> cohort_params;
	vector<void *> measurements_params = { NULL, (void *)&roc_params };

	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it) {
		cohorts[it->first] = filter_range_params;
		cohort_params[it->first] = (void *)&it->second;
	}

	return booststrap_analyze(preds, y, pids, additional_info, cohorts,
		measures, &cohort_params, &measurements_params, fix_cohort_sample_incidence,
		sample_ratio, sample_per_pid, loopCnt);
}

map<string, map<string, float>> booststrap(MedSamples &samples,
	map<string, vector<float>> &additional_info,
	const map<string, vector<Filter_Param>> &filter_cohort,
	const ROC_Params &roc_params, float sample_ratio, int sample_per_pid,
	int loopCnt) {
	vector<float> preds;
	vector<float> y;
	vector<int> pids(samples.nSamples());
	samples.get_y(y);
	samples.get_preds(preds);
	int c = 0;
	for (size_t i = 0; i < samples.idSamples.size(); ++i)
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
			pids[c++] = samples.idSamples[i].samples[j].id;

	bool uses_time_window = true;

	if (uses_time_window) {
		additional_info["Time-Window"].resize(c); //negative value to mark control
		MedTime tm;
		tm.init_time_tables();
		c = 0;
		for (size_t i = 0; i < samples.idSamples.size(); ++i)
			for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j) {
				int diff_days = (tm.convert_date(MedTime::Days,
					samples.idSamples[i].samples[j].outcomeTime)
					- tm.convert_date(MedTime::Days, samples.idSamples[i].samples[j].time));
				if (samples.idSamples[i].samples[j].outcome <= 0)
					diff_days = -diff_days;
				additional_info["Time-Window"][c] = (float)diff_days;
				++c;
			}
	}

	return booststrap_base(preds, y, pids, additional_info,
		filter_cohort, roc_params, sample_ratio, sample_per_pid, loopCnt);
}

map<string, map<string, float>> booststrap(MedFeatures &features,
	const map<string, vector<Filter_Param>> &filter_cohort,
	const ROC_Params &roc_params, float sample_ratio, int sample_per_pid,
	int loopCnt) {
	vector<float> preds((int)features.samples.size());
	vector<float> y((int)features.samples.size());
	vector<int> pids((int)features.samples.size());
	map<string, vector<float>> data = features.data;
	bool uses_time_window = true;

	if (uses_time_window)
		data["Time-Window"].resize((int)features.samples.size());
	MedTime tm;
	tm.init_time_tables();
	for (size_t i = 0; i < features.samples.size(); ++i)
	{
		pids[i] = features.samples[i].id;
		y[i] = features.samples[i].outcome;
		preds[i] = features.samples[i].prediction[0];
		if (uses_time_window) {
			int diff_days = (tm.convert_date(MedTime::Days, features.samples[i].outcomeTime)
				- tm.convert_date(MedTime::Days, features.samples[i].time));
			if (y[i] <= 0)
				diff_days = -diff_days;
			data["Time-Window"][i] = (float)diff_days;
		}
	}

	return booststrap_base(preds, y, pids, data, filter_cohort, roc_params, sample_ratio,
		sample_per_pid, loopCnt);
}

void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedSamples &samples) {
	vector<int> remove_indexes_pid, remove_indexes_sample;
	vector<int> remove_complete_pids;
	int tot_cnt = 0;
	for (size_t i = 0; i < samples.idSamples.size(); ++i) {
		bool keep_pid = false;
		tot_cnt += (int)samples.idSamples[i].samples.size();
		if (pid_censor_dates.find(samples.idSamples[i].id) == pid_censor_dates.end())
			continue; //no censor date
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
		{
			int censor_date = pid_censor_dates.at(samples.idSamples[i].samples[j].id);
			if (samples.idSamples[i].samples[j].time >= censor_date) {
				remove_indexes_pid.push_back((int)i);
				remove_indexes_sample.push_back((int)j);
			}
			else
				keep_pid = true;
		}
		if (!keep_pid)
			remove_complete_pids.push_back((int)i);
	}

	if (remove_indexes_pid.empty())
		return;

	//remove or keep- if remove less than 5% than remove directly, otherwise create new
	bool better_create_new = (double(remove_indexes_pid.size()) / double(tot_cnt)) > 0.05;
	if (better_create_new) {
		int remove_index = 0, remove_index_complete = 0;
		MedSamples new_samples;
		new_samples.time_unit = samples.time_unit;
		//allocate memory:
		new_samples.idSamples.reserve(samples.idSamples.size() - remove_complete_pids.size());
		//samples - filtered:
		for (size_t i = 0; i < samples.idSamples.size(); ++i) {
			if (remove_index_complete < remove_complete_pids.size() &&
				i == remove_complete_pids[remove_index_complete]) {
				++remove_index_complete;
				while (remove_index < remove_indexes_pid.size() &&
					remove_indexes_pid[remove_index] == i)
					++remove_index;
				continue;
			}
			for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
			{
				//remove selected row when matched:
				if (remove_index < remove_indexes_pid.size()
					&& i == remove_indexes_pid[remove_index]
					&& j == remove_indexes_sample[remove_index]) {
					++remove_index;
					continue;
				}
				//insert new sample for pid if needed
				if (new_samples.idSamples.empty() ||
					new_samples.idSamples.back().id != samples.idSamples[i].id)
					new_samples.idSamples.push_back(MedIdSamples(samples.idSamples[i].id));
				new_samples.idSamples.back().samples.push_back(samples.idSamples[i].samples[j]);
			}
		}

		//new_samples.sort_by_id_date(); //not needed
		samples = new_samples;
	}
	else {
		unordered_set<int> removed_pid_set(remove_complete_pids.begin(), remove_complete_pids.end());
		for (int i = (int)remove_indexes_pid.size() - 1; i >= 0; --i) {
			if (removed_pid_set.find(remove_indexes_pid[i]) != removed_pid_set.end())
				continue; //will removed all later
			samples.idSamples[remove_indexes_pid[i]].samples.erase(samples.idSamples[remove_indexes_pid[i]].samples.begin() +
				remove_indexes_sample[i]);
		}
		for (int i = (int)remove_complete_pids.size() - 1; i >= 0; --i)
			samples.idSamples.erase(samples.idSamples.begin() + remove_complete_pids[i]);

		//samples.sort_by_id_date(); //not needed
	}
}
void apply_censor(const unordered_map<int, int> &pid_censor_dates, MedFeatures &features) {
	vector<int> remove_indexes;
	for (size_t i = 0; i < features.samples.size(); ++i)
	{
		if (pid_censor_dates.find(features.samples[i].id) == pid_censor_dates.end())
			continue; //no censor date
		int censor_date = pid_censor_dates.at(features.samples[i].id);
		if (features.samples[i].time >= censor_date)
			remove_indexes.push_back((int)i);
	}
	if (remove_indexes.empty())
		return;

	//remove or keep- if remove less than 5% than remove directly, otherwise create new
	bool better_create_new = (double(remove_indexes.size()) / double(features.samples.size())) > 0.05;
	if (better_create_new) {
		int curr_ind = 0;
		MedFeatures new_feature;
		new_feature.time_unit = features.time_unit;
		new_feature.tags = features.tags;
		new_feature.attributes = features.attributes;
		//alocate memory:
		new_feature.samples.reserve(features.samples.size() - remove_indexes.size());
		for (auto iit = features.data.begin(); iit != features.data.end(); ++iit)
			new_feature.data[iit->first].reserve(features.samples.size() - remove_indexes.size());
		//data and samples - filtered:
		for (int i = 0; i < features.samples.size(); ++i)
		{
			//remove selected row when matched:
			if (curr_ind < remove_indexes.size() && i == remove_indexes[curr_ind]) {
				++curr_ind;
				continue;
			}
			new_feature.samples.push_back(features.samples[i]);
			for (auto iit = features.data.begin(); iit != features.data.end(); ++iit)
				new_feature.data[iit->first].push_back(iit->second[i]);
			if (!features.weights.empty())
				new_feature.weights.push_back(features.weights[i]);
		}

		new_feature.init_pid_pos_len();
		features = new_feature;
	}
	else {
		for (int i = (int)remove_indexes.size() - 1; i >= 0; --i)
		{
			//remove selected rows:
			features.samples.erase(features.samples.begin() + remove_indexes[i]);
			for (auto iit = features.data.begin(); iit != features.data.end(); ++iit)
				features.data[iit->first].erase(features.data[iit->first].begin() + remove_indexes[i]);
			if (!features.weights.empty())
				features.weights.erase(features.weights.begin() + remove_indexes[i]);
		}
		features.init_pid_pos_len();
	}
}
void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedSamples &samples) {
	unordered_map<int, int> pid_censor_dates;
	pid_censor_dates.reserve((int)pids.size());
	for (size_t i = 0; i < pids.size(); ++i)
		pid_censor_dates[pids[i]] = censor_dates[i];
	apply_censor(pid_censor_dates, samples);
}
void apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedFeatures &features) {
	unordered_map<int, int> pid_censor_dates;
	pid_censor_dates.reserve((int)pids.size());
	for (size_t i = 0; i < pids.size(); ++i)
		pid_censor_dates[pids[i]] = censor_dates[i];
	apply_censor(pid_censor_dates, features);
}

void change_sample_autosim(MedSamples &samples, int min_time, int max_time, MedSamples &new_samples) {
	MedTime tm;
	tm.init_time_tables();
	new_samples.time_unit = samples.time_unit;
	for (size_t i = 0; i < samples.idSamples.size(); ++i) {
		int keep_case_index = -1, keep_control_index = -1;
		float max_pred_case = 0, max_pred_control = 0;
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
		{
			int diff_days = (tm.convert_date(MedTime::Days,
				samples.idSamples[i].samples[j].outcomeTime)
				- tm.convert_date(MedTime::Days, samples.idSamples[i].samples[j].time));
			//check if valid for time_window:
			if (samples.idSamples[i].samples[j].outcome >0) {
				if ((diff_days >= min_time && diff_days <= max_time) &&
					(keep_case_index < 0 || max_pred_case < samples.idSamples[i].samples[j].prediction[0])) {
					max_pred_case = samples.idSamples[i].samples[j].prediction[0];
					keep_case_index = (int)j;
				}
			}
			else {
				if ((samples.idSamples[i].samples[j].outcome <= 0 && diff_days >= max_time) &&
					(keep_control_index < 0 || max_pred_control < samples.idSamples[i].samples[j].prediction[0])) {
					max_pred_control = samples.idSamples[i].samples[j].prediction[0];
					keep_control_index = (int)j;
				}
			}

		}
		if (keep_control_index >= 0) {
			new_samples.idSamples.push_back(MedIdSamples(samples.idSamples[i].id));
			new_samples.idSamples.back().samples.push_back(samples.idSamples[i].samples[keep_control_index]);
		}
		if (keep_case_index >= 0) {
			if (keep_control_index < 0)
				new_samples.idSamples.push_back(MedIdSamples(samples.idSamples[i].id));
			new_samples.idSamples.back().samples.push_back(samples.idSamples[i].samples[keep_case_index]);
		}
	}
	//new_samples.sort_by_id_date(); //not needed
}

void change_sample_autosim(MedFeatures &features, int min_time, int max_time, MedFeatures &new_features) {
	MedTime tm;
	tm.init_time_tables();
	new_features.attributes = features.attributes;
	new_features.tags = features.tags;
	new_features.time_unit = features.time_unit;
	unordered_map<int, vector<int>> pid_indexes;
	for (size_t i = 0; i < features.samples.size(); ++i)
		pid_indexes[features.samples[i].id].push_back((int)i);

	for (auto it = pid_indexes.begin(); it != pid_indexes.end(); ++it) {
		int keep_case_index = -1, keep_control_index = -1;
		float max_pred_case = 0, max_pred_control = 0;

		//scan pid indexes
		for (int i : it->second) {
			int diff_days = (tm.convert_date(MedTime::Days, features.samples[i].outcomeTime)
				- tm.convert_date(MedTime::Days, features.samples[i].time));
			//check if valid for time_window:
			if (features.samples[i].outcome > 0) {
				if ((diff_days >= min_time && diff_days <= max_time) &&
					(keep_case_index < 0 || max_pred_case < features.samples[i].prediction[0])) {
					max_pred_case = features.samples[i].prediction[0];
					keep_case_index = (int)i;
				}
			}
			else {
				if ((features.samples[i].outcome <= 0 && diff_days >= max_time) &&
					(keep_control_index < 0 || max_pred_control < features.samples[i].prediction[0])) {
					max_pred_control = features.samples[i].prediction[0];
					keep_control_index = (int)i;
				}
			}
		}

		if (keep_control_index >= 0) {
			new_features.samples.push_back(features.samples[keep_control_index]);
			for (auto iit = features.data.begin(); iit != features.data.end(); ++iit)
				new_features.data[iit->first].push_back(iit->second[keep_control_index]);
			if (!features.weights.empty())
				new_features.weights.push_back(features.weights[keep_control_index]);
		}
		if (keep_case_index >= 0) {
			new_features.samples.push_back(features.samples[keep_case_index]);
			for (auto iit = features.data.begin(); iit != features.data.end(); ++iit)
				new_features.data[iit->first].push_back(iit->second[keep_case_index]);
			if (!features.weights.empty())
				new_features.weights.push_back(features.weights[keep_case_index]);
		}
	}
	new_features.init_pid_pos_len();
}