#include <unordered_set>
#include "MedBootstrap.h"
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include "Logger/Logger/Logger.h"
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL 

MedBootstrap::MedBootstrap() {
	sample_ratio = (float)1.0;
	sample_per_pid = 1;
	loopCnt = 500;
	filter_cohort["All"] = {};
}

MedBootstrap::MedBootstrap(const string &init_string) {
	sample_ratio = (float)1.0;
	sample_per_pid = 1;
	loopCnt = 500;
	filter_cohort["All"] = {};

	//now read init_string to override default:
	vector<string> tokens;
	boost::split(tokens, init_string, boost::is_any_of(";"));
	for (string token : tokens)
	{
		if (token.find('=') == string::npos)
			MTHROW_AND_ERR("Wrong token. has no value \"%s\"\n", token.c_str());
		string param_name = token.substr(0, token.find('='));
		string param_value = token.substr(token.find('=') + 1);
		boost::to_lower(param_name);

		if (param_name == "sample_ratio") {
			sample_ratio = stof(param_value);
			if (sample_ratio > 1.0 || sample_ratio < 0)
				MTHROW_AND_ERR("sample_ratio should be between 0-1, got %2.3f\n", sample_ratio);
		}
		else if (param_name == "sample_per_pid")
			sample_per_pid = stoi(param_value);
		else if (param_name == "loopCnt")
			loopCnt = stoi(param_value);
		else if (param_name == "roc_params")
			roc_Params = ROC_Params(param_value);
		else if (param_name == "filter_cohort") {
			ifstream of(param_value);
			string line;
			while (getline(of, line)) {
				if (line.empty() || boost::starts_with(line, "#"))
					continue;
				if (line.find('\t') == string::npos)
					MTHROW_AND_ERR("filter_cohort file \"%s\" is in wing format. line=\"%s\"\n",
						param_value.c_str(), line.c_str());
				string cohort_name = line.substr(0, line.find('\t'));
				string cohort_definition = line.substr(line.find('\t') + 1);
				vector<string> params;
				boost::split(params, cohort_definition, boost::is_any_of(";"));
				vector<Filter_Param> convert_params((int)params.size());
				for (size_t i = 0; i < params.size(); ++i)
					convert_params[i] = Filter_Param(params[i]);
				filter_cohort[cohort_name] = convert_params;
			}
			of.close();
		}
		else
			MTHROW_AND_ERR("Unknown paramter \"%s\" for ROC_Params\n", param_name.c_str());
	}
}

map<string, map<string, float>> MedBootstrap::booststrap_base(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
	const map<string, vector<float>> &additional_info) {

	map<string, FilterCohortFunc> cohorts;
	vector<MeasurementFunctions> measures = { calc_roc_measures_with_inc };
	map<string, void *> cohort_params;
	vector<void *> measurements_params = { NULL, (void *)&roc_Params };

	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it) {
		cohorts[it->first] = filter_range_params;
		cohort_params[it->first] = (void *)&it->second;
	}

	return booststrap_analyze(preds, y, pids, additional_info, cohorts,
		measures, &cohort_params, &measurements_params, fix_cohort_sample_incidence,
		sample_ratio, sample_per_pid, loopCnt);
}

bool MedBootstrap::use_time_window() {
	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it)
		for (size_t i = 0; i < it->second.size(); ++i)
			if (boost::to_lower_copy(it->second[i].param_name) == "time-window")
				return true;
	return false;
}

void MedBootstrap::clean_feature_name_prefix(map<string, vector<float>> &features) {
	map<string, vector<float>> name_data;
	for (auto it = features.begin(); it != features.end(); ++it)
	{
		if (boost::starts_with(it->first, "FTR_") && it->first.find('.') != string::npos) {
			string new_name = it->first.substr(it->first.find('.') + 1);
			name_data[new_name].swap(it->second);
		}
	}
	name_data.swap(features); //feature_data will steal "name_data" data
}

map<string, map<string, float>> MedBootstrap::booststrap(MedFeatures &features) {
	vector<float> preds((int)features.samples.size());
	vector<float> y((int)features.samples.size());
	vector<int> pids((int)features.samples.size());
	map<string, vector<float>> data = features.data;
	clean_feature_name_prefix(data);
	bool uses_time_window = use_time_window();

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

	return booststrap_base(preds, y, pids, data);
}

map<string, map<string, float>> MedBootstrap::booststrap(MedSamples &samples,
	map<string, vector<float>> &additional_info) {
	vector<float> preds;
	vector<float> y;
	vector<int> pids(samples.nSamples());
	samples.get_y(y);
	samples.get_preds(preds);
	int c = 0;
	for (size_t i = 0; i < samples.idSamples.size(); ++i)
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
			pids[c++] = samples.idSamples[i].samples[j].id;

	bool uses_time_window = use_time_window();

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

	return booststrap_base(preds, y, pids, additional_info);
}

map<string, map<string, float>> MedBootstrap::booststrap(MedSamples &samples, const string &rep_path) {
	MedModel mdl;
	mdl.add_age();
	mdl.add_gender();
	vector<int> pids_to_take;
	samples.get_ids(pids_to_take);

	unordered_set<string> req_names;
	mdl.get_required_signal_names(req_names);
	vector<string> sigs = { "BYEAR", "GENDER", "TRAIN" };
	for (string s : req_names)
		sigs.push_back(s);
	sort(sigs.begin(), sigs.end());
	auto it = unique(sigs.begin(), sigs.end());
	sigs.resize(std::distance(sigs.begin(), it));

	MedPidRepository rep;
	if (rep.read_all(rep_path, pids_to_take, sigs) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", rep_path.c_str());

	if (mdl.apply(rep, samples, MedModelStage::MED_MDL_APPLY_FTR_GENERATORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS) < 0)
		MTHROW_AND_ERR("Error creating age,gender for samples\n");

	return booststrap(mdl.features);
}

void MedBootstrap::apply_censor(const unordered_map<int, int> &pid_censor_dates, MedSamples &samples) {
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
void MedBootstrap::apply_censor(const unordered_map<int, int> &pid_censor_dates, MedFeatures &features) {
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
void MedBootstrap::apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedSamples &samples) {
	unordered_map<int, int> pid_censor_dates;
	pid_censor_dates.reserve((int)pids.size());
	for (size_t i = 0; i < pids.size(); ++i)
		pid_censor_dates[pids[i]] = censor_dates[i];
	apply_censor(pid_censor_dates, samples);
}
void MedBootstrap::apply_censor(const vector<int> &pids, const vector<int> &censor_dates, MedFeatures &features) {
	unordered_map<int, int> pid_censor_dates;
	pid_censor_dates.reserve((int)pids.size());
	for (size_t i = 0; i < pids.size(); ++i)
		pid_censor_dates[pids[i]] = censor_dates[i];
	apply_censor(pid_censor_dates, features);
}

void MedBootstrap::change_sample_autosim(MedSamples &samples, int min_time, int max_time, MedSamples &new_samples) {
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

void MedBootstrap::change_sample_autosim(MedFeatures &features, int min_time, int max_time, MedFeatures &new_features) {
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

void MedBootstrapResult::booststrap(MedFeatures &features) {
	bootstrap_results = bootstrap_params.booststrap(features);
}
void MedBootstrapResult::booststrap(MedSamples &samples, map<string, vector<float>> &additional_info) {
	bootstrap_results = bootstrap_params.booststrap(samples, additional_info);
}
void MedBootstrapResult::booststrap(MedSamples &samples, const string &rep_path) {
	bootstrap_results = bootstrap_params.booststrap(samples, rep_path);
}

bool MedBootstrapResult::find_in_range(const vector<float> &vec, float search, float th) {
	for (size_t i = 0; i < vec.size(); ++i)
	{
		float diff = abs(vec[i] - search);
		if (diff < th)
			return true;
	}
	return false;
}

void MedBootstrapResult::find_working_points(const map<string, float> &bootstrap_cohort,
	vector<float> &sens_points, vector<float> &pr_points) {
	unordered_set<string> all_columns_uniq;
	float min_dif = (float)1.0; //remove too close working points

	for (auto it = bootstrap_cohort.begin(); it != bootstrap_cohort.end(); ++it)
		all_columns_uniq.insert(it->first);
	vector<string> all_columns(all_columns_uniq.begin(), all_columns_uniq.end());
	//search for @SCORE when SENS, and PR are not 0 or 1:
	for (string column_name : all_columns)
	{
		if (column_name.size() > 7 && column_name.substr(0, 7) == "SENS@SC" && column_name.substr(column_name.size() - 4) == "_Obs")  //sens@score
			if (bootstrap_cohort.at(column_name) != MED_MAT_MISSING_VALUE && bootstrap_cohort.at(column_name) != 0
				&& bootstrap_cohort.at(column_name) != 1) {
				float cand_val = round(10000 * bootstrap_cohort.at(column_name)) / 100;
				if (!find_in_range(sens_points, cand_val, min_dif))
					sens_points.push_back(cand_val);
			}

		if (column_name.size() > 5 && column_name.substr(0, 5) == "PR@SC"  && column_name.substr(column_name.size() - 4) == "_Obs")  //sens@score
			if (bootstrap_cohort.at(column_name) != MED_MAT_MISSING_VALUE && bootstrap_cohort.at(column_name) != 0
				&& bootstrap_cohort.at(column_name) != 1) {
				float cand_val = round(10000 * bootstrap_cohort.at(column_name)) / 100;
				if (!find_in_range(pr_points, cand_val, min_dif))
					pr_points.push_back(cand_val);
			}
	}

}

void MedBootstrapResult::explore_measure(const string &measure_name, float value, map<string, float> &score_measurements,
	const string &string_cohort, float max_search_range) {
	if (bootstrap_results.empty())
		MTHROW_AND_ERR("Please run bootstrap first\n");
	if (bootstrap_results.find(string_cohort) == bootstrap_results.end())
		MTHROW_AND_ERR("Couldn't find \"%s\" cohort in bootstrap\n", string_cohort.c_str());
	map<string, float> all_measures = bootstrap_results.at(string_cohort);
	string value_att = "@" + boost::to_upper_copy(measure_name) + "_";

	if (!bootstrap_params.roc_Params.use_score_working_points)
		MWARN("Warnning:: You have ran bootstrap without score preformance!\n");
	//search for closest score
	unordered_map<string, float> lower_bound, upper_bound;
	unordered_map<string, float> lower_val, upper_val;
	unordered_set<string> all_measures_opts;
	for (auto it = all_measures.begin(); it != all_measures.end(); ++it)
	{
		string column_name = it->first;
		if (column_name.find(value_att) != string::npos && boost::ends_with(column_name, "_Mean")) {
			string measure = column_name.substr(0, column_name.find_first_of('@'));
			all_measures_opts.insert(measure);
			if (it->second != MED_MAT_MISSING_VALUE) {
				float score_val = stof(column_name.substr(column_name.find(value_att) + value_att.length()));

				if (abs(score_val - value) <= max_search_range) {
					if (lower_bound.find(measure) == lower_bound.end() ||
						(lower_bound[measure] > score_val && score_val < value)) {
						lower_bound[measure] = score_val;
						lower_val[measure] = it->second;
					}
					if (upper_bound.find(measure) == upper_bound.end() ||
						(upper_bound[measure] < score_val && score_val > value)) {
						upper_bound[measure] = score_val;
						upper_val[measure] = it->second;
					}
				}
			}
		}
	}

	//calc interp for the score measurements:
	for (auto it = all_measures_opts.begin(); it != all_measures_opts.end(); ++it)
	{
		if (lower_bound.find(*it) == lower_bound.end() && upper_bound.find(*it) == upper_bound.end()) {
			MWARN("Couldn't find measure \"%s\" in search_range %2.3f\n",
				(*it).c_str(), max_search_range);
			continue;
		}

		//only upper exists:
		if (lower_bound.find(*it) == lower_bound.end()) {
			score_measurements[*it] = upper_val[*it];
			continue;
		}
		//only lower exists:
		if (upper_bound.find(*it) == upper_bound.end()) {
			score_measurements[*it] = lower_val[*it];
			continue;
		}
		//interp between both:
		float lower_diff = value - lower_bound[*it];
		float upper_diff = upper_bound[*it] - value;
		float tot_diff = lower_diff + upper_diff;
		float final_val = (upper_diff / tot_diff) * lower_val[*it] +
			(lower_diff / tot_diff) * upper_val[*it]; //in inverse order to distance
		score_measurements[*it] = final_val;
	}
}

void MedBootstrapResult::explore_score(float score, map<string, float> &score_measurements,
	const string &string_cohort, float max_search_range) {
	explore_measure("SCORE", score, score_measurements, string_cohort, max_search_range);
}

void MedBootstrapResult::write_results_to_text_file(const string &path) {
	write_bootstrap_results(path, bootstrap_results);
}

void MedBootstrapResult::read_results_to_text_file(const string &path) {
	read_bootstrap_results(path, bootstrap_results);
}