#include <unordered_set>
#include "MedBootstrap.h"
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include "Logger/Logger/Logger.h"
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL 

void MedBootstrap::parse_cohort_file(const string &cohorts_path) {
	ifstream of(cohorts_path);
	if (!of.good())
		MTHROW_AND_ERR("IO Error: can't read \"%s\"\n", cohorts_path.c_str());
	string line;
	while (getline(of, line)) {
		if (line.empty() || boost::starts_with(line, "#"))
			continue;
		if (line.find('\t') == string::npos)
			MTHROW_AND_ERR("filter_cohort file \"%s\" is in wing format. line=\"%s\"\n",
				cohorts_path.c_str(), line.c_str());
		string cohort_name = line.substr(0, line.find('\t'));
		string cohort_definition = line.substr(line.find('\t') + 1);
		if (cohort_name != "MULTI") {
			vector<string> params;
			boost::split(params, cohort_definition, boost::is_any_of(";"));
			vector<Filter_Param> convert_params((int)params.size());
			for (size_t i = 0; i < params.size(); ++i)
				convert_params[i] = Filter_Param(params[i]);
			filter_cohort[cohort_name] = convert_params;
		}
		else {
			vector<string> tokens_p;
			boost::split(tokens_p, cohort_definition, boost::is_any_of("\t"));
			vector<vector<Filter_Param>> multi_line;
			for (size_t i = 0; i < tokens_p.size(); ++i)
			{
				vector<string> params;
				boost::split(params, tokens_p[i], boost::is_any_of(";"));
				vector<Filter_Param> convert_params((int)params.size());
				for (size_t k = 0; k < params.size(); ++k)
					convert_params[k] = Filter_Param(params[k]);
				multi_line.push_back(convert_params);
			}

			add_filter_cohorts(multi_line);
		}
	}
	of.close();
}

MedBootstrap::MedBootstrap()
{
	sample_ratio = (float)1.0;
	sample_per_pid = 1;
	sample_patient_label = false;
	loopCnt = 500;
	filter_cohort["All"] = {};
}

MedBootstrap::MedBootstrap(const string &init_string) {
	sample_ratio = (float)1.0;
	sample_per_pid = 1;
	sample_patient_label = false;
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
		else if (param_name == "loopcnt")
			loopCnt = stoi(param_value);
		else if (param_name == "sample_patient_label")
			sample_patient_label = stoi(param_value) > 0;
		else if (param_name == "roc_params")
			roc_Params = ROC_Params(param_value);
		else if (param_name == "filter_cohort") {
			filter_cohort.clear();
			parse_cohort_file(param_value);
		}
		else
			MTHROW_AND_ERR("Unknown paramter \"%s\" for ROC_Params\n", param_name.c_str());
	}
}

map<string, map<string, float>> MedBootstrap::bootstrap_base(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
	const map<string, vector<float>> &additional_info) {

	map<string, FilterCohortFunc> cohorts;
	vector<MeasurementFunctions> measures = { calc_roc_measures_with_inc };
	map<string, void *> cohort_params;
	vector<void *> measurements_params = { (void *)&roc_Params };

	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it) {
		cohorts[it->first] = filter_range_params;
		cohort_params[it->first] = (void *)&it->second;
	}

	vector<int> new_ids((int)pids.size());
	if (sample_patient_label)
	{
		for (size_t i = 0; i < pids.size(); ++i)
			new_ids[i] = pids[i] * 10 + (int)y[i];
		return booststrap_analyze(preds, y, new_ids, additional_info, cohorts,
			measures, &cohort_params, &measurements_params, fix_cohort_sample_incidence,
			preprocess_bin_scores, &roc_Params, sample_ratio, sample_per_pid, loopCnt);
	}
	return booststrap_analyze(preds, y, pids, additional_info, cohorts,
		measures, &cohort_params, &measurements_params, fix_cohort_sample_incidence,
		preprocess_bin_scores, &roc_Params, sample_ratio, sample_per_pid, loopCnt);
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
		else
			name_data[it->first].swap(it->second);
	}
	name_data.swap(features); //feature_data will steal "name_data" data
}

void MedBootstrap::add_splits_results(const vector<float> &preds, const vector<float> &y,
	const vector<int> &pids, const map<string, vector<float>> &data,
	const unordered_map<int, vector<int>> &splits_inds,
	map<int, map<string, map<string, float>>> &results_per_split) {

	results_per_split.clear();
	for (auto it = splits_inds.begin(); it != splits_inds.end(); ++it)
	{
		int split_id = it->first;
		vector<float> split_preds((int)it->second.size());
		vector<float> split_y((int)it->second.size());
		vector<int> split_pids((int)it->second.size());
		map<string, vector<float>> split_data;
		for (auto jt = data.begin(); jt != data.end(); ++jt)
			split_data[jt->first].resize((int)it->second.size());
		for (size_t i = 0; i < split_preds.size(); ++i)
		{
			split_preds[i] = preds[it->second[i]];
			split_y[i] = y[it->second[i]];
			split_pids[i] = pids[it->second[i]];
			for (auto jt = data.begin(); jt != data.end(); ++jt)
				split_data[jt->first][i] = jt->second[i];
		}
		results_per_split[split_id] =
			bootstrap_base(split_preds, split_y, split_pids, split_data);
	}

}

void rec_filter_cohorts(const map<string, vector<pair<float, float>>> &parameters_ranges,
	map<string, vector<Filter_Param>> &current_filters) {
	if (parameters_ranges.empty())
		return;
	auto curr = parameters_ranges.begin();
	if (current_filters.empty()) {
		for (size_t i = 0; i < curr->second.size(); ++i)
		{
			char buff[1000];
			snprintf(buff, sizeof(buff), "%s:%1.3f-%1.3f", curr->first.c_str(),
				curr->second[i].first, curr->second[i].second);
			Filter_Param fp;
			fp.param_name = curr->first;
			fp.min_range = curr->second[i].first;
			fp.max_range = curr->second[i].second;
			current_filters[string(buff)].push_back(fp);
		}
		++curr;
		map<string, vector<pair<float, float>>> rest(curr, parameters_ranges.end());
		rec_filter_cohorts(rest, current_filters);
		return;
	}

	map<string, vector<Filter_Param>> new_step;
	for (auto it = current_filters.begin(); it != current_filters.end(); ++it) {
		for (size_t i = 0; i < curr->second.size(); ++i)
		{
			//add to it curr->second[i] 
			char buff[1000];
			snprintf(buff, sizeof(buff), "%s,%s:%1.3f-%1.3f", it->first.c_str(), curr->first.c_str(),
				curr->second[i].first, curr->second[i].second);
			Filter_Param fp;
			fp.param_name = curr->first;
			fp.min_range = curr->second[i].first;
			fp.max_range = curr->second[i].second;
			string str_name = string(buff);
			new_step[str_name] = it->second;
			new_step[str_name].push_back(fp);
		}
	}

	++curr;
	map<string, vector<pair<float, float>>> rest(curr, parameters_ranges.end());
	current_filters = new_step;
	rec_filter_cohorts(rest, current_filters);
}

void MedBootstrap::add_filter_cohorts(const map<string, vector<pair<float, float>>> &parameters_ranges) {
	map<string, vector<Filter_Param>> filters;
	rec_filter_cohorts(parameters_ranges, filters);
	//add filters to object:
	for (auto it = filters.begin(); it != filters.end(); ++it)
		filter_cohort[it->first] = it->second;
}
void rec_filter_cohorts(const vector<vector<Filter_Param>> &parameters_ranges,
	map<string, vector<Filter_Param>> &current_filters) {
	if (parameters_ranges.empty())
		return;
	vector<Filter_Param> curr = parameters_ranges.front();
	if (current_filters.empty()) {
		for (size_t i = 0; i < curr.size(); ++i)
		{
			char buff[1000];
			snprintf(buff, sizeof(buff), "%s:%1.3f-%1.3f", curr[i].param_name.c_str(),
				curr[i].min_range, curr[i].max_range);
			current_filters[string(buff)].push_back(curr[i]);
		}
	}
	else {
		map<string, vector<Filter_Param>> new_step;
		for (auto it = current_filters.begin(); it != current_filters.end(); ++it) {
			for (size_t i = 0; i < curr.size(); ++i)
			{
				//add to it curr->second[i] 
				char buff[1000];
				snprintf(buff, sizeof(buff), "%s,%s:%1.3f-%1.3f", it->first.c_str(), curr[i].param_name.c_str(),
					curr[i].min_range, curr[i].max_range);
				string str_name = string(buff);
				new_step[str_name] = it->second;
				new_step[str_name].push_back(curr[i]);
			}
		}
		current_filters = new_step;
	}

	vector<vector<Filter_Param>> rest(parameters_ranges.begin() + 1, parameters_ranges.end());
	rec_filter_cohorts(rest, current_filters);

}
void MedBootstrap::add_filter_cohorts(const vector<vector<Filter_Param>> &parameters_ranges) {
	map<string, vector<Filter_Param>> filters;
	rec_filter_cohorts(parameters_ranges, filters);
	//add filters to object:
	for (auto it = filters.begin(); it != filters.end(); ++it)
		filter_cohort[it->first] = it->second;
}

void MedBootstrap::prepare_bootstrap(MedFeatures &features, vector<float> &preds, vector<float> &y, vector<int> &pids,
	map<string, vector<float>> &final_additional_info, unordered_map<int, vector<int>> *splits_inds) {
	preds.resize((int)features.samples.size());
	y.resize((int)features.samples.size());
	pids.resize((int)features.samples.size());
	final_additional_info = features.data;
	clean_feature_name_prefix(final_additional_info);
	bool uses_time_window = use_time_window();

	if (uses_time_window) {
		final_additional_info["Time-Window"].resize((int)features.samples.size());
		final_additional_info["Label"].resize((int)features.samples.size());
	}
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
			final_additional_info["Time-Window"][i] = (float)diff_days;
			final_additional_info["Label"][i] = y[i];
		}
		if (splits_inds != NULL)
			(*splits_inds)[features.samples[i].split].push_back((int)i);
	}
}
map<string, map<string, float>> MedBootstrap::bootstrap(MedFeatures &features,
	map<int, map<string, map<string, float>>> *results_per_split) {

	vector<float> preds, y;
	vector<int> pids;
	map<string, vector<float>> data;
	unordered_map<int, vector<int>> splits_inds;

	if (results_per_split != NULL)
		prepare_bootstrap(features, preds, y, pids, data, &splits_inds);
	else
		prepare_bootstrap(features, preds, y, pids, data);

	if (results_per_split != NULL)
		add_splits_results(preds, y, pids, data, splits_inds, *results_per_split);

	return bootstrap_base(preds, y, pids, data);
}

void MedBootstrap::prepare_bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info, vector<float> &preds, vector<float> &y, vector<int> &pids,
	unordered_map<int, vector<int>> *splits_inds) {
	pids.resize(samples.nSamples());
	samples.get_y(y);
	samples.get_preds(preds);
	int c = 0;
	for (size_t i = 0; i < samples.idSamples.size(); ++i)
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j) {
			pids[c++] = samples.idSamples[i].samples[j].id;
			if (splits_inds != NULL)
				(*splits_inds)[samples.idSamples[i].samples[j].split].push_back((int)c);
		}

	bool uses_time_window = use_time_window();

	if (uses_time_window) {
		additional_info["Time-Window"].resize(c); //negative value to mark control
		additional_info["Label"].resize(c); //negative value to mark control
		MedTime tm;
		tm.init_time_tables();
		c = 0;
		for (size_t i = 0; i < samples.idSamples.size(); ++i)
			for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j) {
				int diff_days = (tm.convert_date(MedTime::Days,
					samples.idSamples[i].samples[j].outcomeTime)
					- tm.convert_date(MedTime::Days, samples.idSamples[i].samples[j].time));
				additional_info["Time-Window"][c] = (float)diff_days;
				additional_info["Label"][c] = samples.idSamples[i].samples[j].outcome;
				++c;
			}
	}
}
map<string, map<string, float>> MedBootstrap::bootstrap(MedSamples &samples,
	map<string, vector<float>> &additional_info, map<int, map<string, map<string, float>>> *results_per_split) {
	vector<float> preds, y;
	vector<int> pids;
	unordered_map<int, vector<int>> splits_inds;
	if (results_per_split == NULL)
		prepare_bootstrap(samples, additional_info, preds, y, pids);
	else
		prepare_bootstrap(samples, additional_info, preds, y, pids, &splits_inds);

	if (results_per_split != NULL)
		add_splits_results(preds, y, pids, additional_info, splits_inds, *results_per_split);
	return bootstrap_base(preds, y, pids, additional_info);
}

map<string, map<string, float>> MedBootstrap::bootstrap(MedSamples &samples, const string &rep_path, map<int, map<string, map<string, float>>> *results_per_split) {
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
	int curr_level = global_logger.levels.front();
	global_logger.init_all_levels(LOG_DEF_LEVEL);
	if (rep.read_all(rep_path, pids_to_take, sigs) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", rep_path.c_str());
	global_logger.init_all_levels(curr_level);

	if (mdl.apply(rep, samples, MedModelStage::MED_MDL_APPLY_FTR_GENERATORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS) < 0)
		MTHROW_AND_ERR("Error creating age,gender for samples\n");

	return bootstrap(mdl.features, results_per_split);
}

map<string, map<string, float>> MedBootstrap::bootstrap(MedSamples &samples, MedPidRepository &rep, map<int, map<string, map<string, float>>> *results_per_split) {
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
	if (rep.load(sigs, pids_to_take) < 0)
		MTHROW_AND_ERR("Couldn't load signals in given repository\n");

	if (mdl.apply(rep, samples, MedModelStage::MED_MDL_APPLY_FTR_GENERATORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS) < 0)
		MTHROW_AND_ERR("Error creating age,gender for samples\n");

	return bootstrap(mdl.features, results_per_split);
}

void MedBootstrap::apply_censor(const unordered_map<int, int> &pid_censor_dates, MedSamples &samples) {
	for (size_t i = 0; i < samples.idSamples.size(); ++i) {
		if (pid_censor_dates.find(samples.idSamples[i].id) == pid_censor_dates.end())
			continue; //no censor date
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
		{
			int censor_date = pid_censor_dates.at(samples.idSamples[i].samples[j].id);
			//act only on controls:
			if (samples.idSamples[i].samples[j].outcome <= 0 &&
				samples.idSamples[i].samples[j].outcomeTime > censor_date) {
				samples.idSamples[i].samples[j].outcomeTime = censor_date;
			}
		}
	}
}
void MedBootstrap::apply_censor(const unordered_map<int, int> &pid_censor_dates, MedFeatures &features) {
	for (size_t i = 0; i < features.samples.size(); ++i)
	{
		if (pid_censor_dates.find(features.samples[i].id) == pid_censor_dates.end())
			continue; //no censor date
		int censor_date = pid_censor_dates.at(features.samples[i].id);
		if (features.samples[i].outcome <= 0 && features.samples[i].outcomeTime > censor_date)
			features.samples[i].outcomeTime = censor_date;
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

void MedBootstrapResult::bootstrap(MedFeatures &features) {
	bootstrap_results = bootstrap_params.bootstrap(features);
}
void MedBootstrapResult::bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info) {
	bootstrap_results = bootstrap_params.bootstrap(samples, additional_info);
}
void MedBootstrapResult::bootstrap(MedSamples &samples, const string &rep_path) {
	bootstrap_results = bootstrap_params.bootstrap(samples, rep_path);
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
				float cand_val = round(100 * bootstrap_cohort.at(column_name)) / 100;
				if (!find_in_range(sens_points, cand_val, min_dif))
					sens_points.push_back(cand_val);
			}

		if (column_name.size() > 5 && column_name.substr(0, 5) == "PR@SC"  && column_name.substr(column_name.size() - 4) == "_Obs")  //sens@score
			if (bootstrap_cohort.at(column_name) != MED_MAT_MISSING_VALUE && bootstrap_cohort.at(column_name) != 0
				&& bootstrap_cohort.at(column_name) != 1) {
				float cand_val = round(100 * bootstrap_cohort.at(column_name)) / 100;
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

void MedBootstrapResult::write_results_to_text_file(const string &path, bool pivot_format) {
	if (pivot_format)
		write_pivot_bootstrap_results(path, bootstrap_results);
	else
		write_bootstrap_results(path, bootstrap_results);
}

void MedBootstrapResult::read_results_to_text_file(const string &path, bool pivot_format) {
	if (pivot_format)
		read_pivot_bootstrap_results(path, bootstrap_results);
	else
		read_bootstrap_results(path, bootstrap_results);
}