#include <unordered_set>
#include "MedBootstrap.h"
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include "Logger/Logger/Logger.h"
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL 

void MedBootstrap::get_cohort_from_arg(const string &single_cohort) {
	string cohort_line_with_name = single_cohort + "\t" + single_cohort;
	parse_cohort_line(cohort_line_with_name);
}

void MedBootstrap::parse_cohort_line(const string &line) {
	if (line.find('\t') == string::npos)
		MTHROW_AND_ERR("line is in wrong format:\n%s\n",
			line.c_str());

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

void MedBootstrap::parse_cohort_file(const string &cohorts_path) {
	ifstream of(cohorts_path);
	if (!of.good())
		MTHROW_AND_ERR("IO Error: can't read \"%s\"\n", cohorts_path.c_str());
	string line;
	while (getline(of, line)) {
		boost::trim(line);
		if (line.empty() || boost::starts_with(line, "#"))
			continue;
		parse_cohort_line(line);
	}
	of.close();
}

MedBootstrap::MedBootstrap()
{
	sample_ratio = (float)1.0;
	sample_per_pid = 1;
	sample_patient_label = false;
	sample_seed = 0;
	loopCnt = 500;
	filter_cohort["All"] = {};
	simTimeWindow = false;
}

MedBootstrap::MedBootstrap(const string &init_string) {
	sample_ratio = (float)1.0;
	sample_per_pid = 1;
	sample_patient_label = false;
	sample_seed = 0;
	loopCnt = 500;
	filter_cohort["All"] = {};
	simTimeWindow = false;

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
		else if (param_name == "sample_seed")
			sample_seed = stoi(param_value);
		else if (param_name == "sample_patient_label")
			sample_patient_label = stoi(param_value) > 0;
		else if (param_name == "roc_params")
			roc_Params = ROC_Params(param_value);
		else if (param_name == "filter_cohort") {
			filter_cohort.clear();
			parse_cohort_file(param_value);
		}
		else if (param_name == "simTimeWindow")
			simTimeWindow = stoi(param_value) > 0;
		else
			MTHROW_AND_ERR("Unknown paramter \"%s\" for ROC_Params\n", param_name.c_str());
	}
}

template<class T> bool has_element(const vector<T> &vec, const T &val) {
	for (size_t i = 0; i < vec.size(); ++i)
		if (vec[i] == val)
			return true;
	return false;
}
bool has_element(const vector<Filter_Param> &vec, const string &val) {
	for (size_t i = 0; i < vec.size(); ++i)
		if (vec[i].param_name == val)
			return true;
	return false;
}
void init_model(MedModel &mdl, MedRepository& rep, const string &json_model,
	const string &rep_path, const vector<int> &pids_to_take) {
	if (!json_model.empty()) {
		MLOG("Adding new features using %s\n", json_model.c_str());
		mdl.init_from_json_file(json_model);
	}

	bool need_age = true, need_gender = true;
	for (FeatureGenerator *generator : mdl.generators) {
		if (generator->generator_type == FTR_GEN_AGE)
			need_age = false;
		if (generator->generator_type == FTR_GEN_GENDER)
			need_gender = false;
	}
	if (need_age)
		mdl.add_age();
	if (need_gender)
		mdl.add_gender();

	unordered_set<string> req_names;
	for (FeatureGenerator *generator : mdl.generators)
		generator->get_required_signal_names(req_names);
	vector<string> sigs = { "BYEAR", "GENDER" };
	for (string s : req_names)
		sigs.push_back(s);
	sort(sigs.begin(), sigs.end());
	auto it = unique(sigs.begin(), sigs.end());
	sigs.resize(std::distance(sigs.begin(), it));
	bool clean_rep = medial::process::clean_uneeded_rep_processors(mdl, sigs);
	if (clean_rep)
		MWARN("WARNING TT: Cleaned unused Rep_Processors\n");

	int curr_level = global_logger.levels.front();
	global_logger.init_all_levels(LOG_DEF_LEVEL);
	if (rep.read_all(rep_path, pids_to_take, sigs) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", rep_path.c_str());
	global_logger.init_all_levels(curr_level);
}
void get_data_for_filter(const string &json_model, const string &rep_path,
	MedBootstrap &single_cohort, const vector<MedRegistryRecord> &registry_records,
	MedSamplingStrategy &sampler, map<string, vector<float>> &data_for_filtering,
	vector<MedSample> &inc_smps) {
	MedSamples inc_samples;
	sampler.do_sample(registry_records, inc_samples);
	MLOG("Done sampling for incidence by year. has %d patients\n",
		(int)inc_samples.idSamples.size());
	vector<int> pids_to_take;
	inc_samples.get_ids(pids_to_take);
	MedModel mdl;
	MedPidRepository rep;
	init_model(mdl, rep, json_model, rep_path, pids_to_take);

	if (mdl.learn(rep, &inc_samples, MedModelStage::MED_MDL_LEARN_REP_PROCESSORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS) < 0)
		MTHROW_AND_ERR("Error creating age,gender for samples\n");

	vector<float> preds, y;
	vector<int> pids;
	unordered_map<int, vector<int>> splits_inds;
	for (size_t i = 0; i < mdl.features.samples.size(); ++i)
		mdl.features.samples[i].prediction.resize(1, 0); //set to 0 that won't throw error in prepare
	single_cohort.prepare_bootstrap(mdl.features, preds, y, pids, data_for_filtering);
	inc_smps.swap(mdl.features.samples);
}

void medial::process::make_sim_time_window(const string &cohort_name, const vector<Filter_Param> &filter_p,
	const vector<float> &y, const map<string, vector<float>> &additional_info,
	vector<float> &y_changed, map<string, vector<float>> &cp_info,
	map<string, FilterCohortFunc> &cohorts_t, map<string, void *> &cohort_params_t) {
	string search_term = "Time-Window";
	cohorts_t[cohort_name] = filter_range_params;
	cohort_params_t[cohort_name] = (void *)&filter_p;
	//lets find the element:
	Filter_Param time_filter;
	for (size_t i = 0; i < filter_p.size(); ++i)
		if (filter_p[i].param_name == search_term)
		{
			time_filter = filter_p[i];
			break;
		}

	//update: y_changed based on y
	cp_info = additional_info; //a copy - i will change those each time
	y_changed.insert(y_changed.end(), y.begin(), y.end());
	//update: cp_info["Time-Window"] based on additional_info["Time-Window"] for cases not in time window
	int cases_censored = 0, cases_changed_to_controls = 0;
	for (size_t i = 0; i < y_changed.size(); ++i)
	{
		if (y[i] > 0 && !filter_range_param(additional_info, (int)i, &time_filter)) {
			// cases which are long before the outcome (>2*max_range) are considered as controls:

			// ----------------max_range*2--------max_range---------min_range---event---------------
			// -------------------^-------------------^----------------^----------^-----------------
			// ------control------|------censor-------|------case------|--------censor--------------

			if (additional_info.at(search_term)[i] > time_filter.max_range * 2) {
				y_changed[i] = 0;
				cp_info[search_term][i] = time_filter.min_range; // bogus time - wont filter because in time range
				cases_changed_to_controls++;
			}
			else
				cases_censored++;
		}
	}
	MLOG("make_sim_time_window for cohort [%s] censored %d cases and changed %d cases to controls as they were %f days before outcome date\n",
		cohort_name.c_str(), cases_censored, cases_changed_to_controls, time_filter.max_range * 2);
}

map<string, map<string, float>> MedBootstrap::bootstrap_base(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
	const map<string, vector<float>> &additional_info) {

	map<string, FilterCohortFunc> cohorts;
	map<string, void *> cohort_params;
	if (measurements_with_params.empty()) { //if not touched(than empty) - this is the default!
		measurements_with_params.resize(1);
		measurements_with_params[0] = pair<MeasurementFunctions, void*>(calc_roc_measures_with_inc, &roc_Params);
	}
	vector<MeasurementFunctions> measures((int)measurements_with_params.size());
	vector<void *> measurements_params((int)measurements_with_params.size());
	for (size_t i = 0; i < measurements_with_params.size(); ++i)
	{
		measures[i] = measurements_with_params[i].first;
		measurements_params[i] = measurements_with_params[i].second;
	}

	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it) {
		cohorts[it->first] = filter_range_params;
		cohort_params[it->first] = (void *)&it->second;
	}
	for (auto it = additional_cohorts.begin(); it != additional_cohorts.end(); ++it)
		if (cohorts.find(it->first) == cohorts.end())
			cohorts[it->first] = it->second;
		else
			MWARN("Cohort \"%s\" name already exists - skip additional\n");
	const vector<int> *rep_ids = &pids;
	vector<int> new_ids;
	if (sample_patient_label)
	{
		new_ids.resize((int)pids.size());
		for (size_t i = 0; i < pids.size(); ++i)
			new_ids[i] = pids[i] * 10 + (int)y[i];
		rep_ids = &new_ids;
	}

	if (!simTimeWindow) {
		return booststrap_analyze(preds, y, *rep_ids, additional_info, cohorts,
			measures, &cohort_params, &measurements_params, fix_cohort_sample_incidence,
			preprocess_bin_scores, &roc_Params, sample_ratio, sample_per_pid, loopCnt, sample_seed);
	}
	else {
		//split by each time window after editing each time window samples:
		string search_term = "Time-Window";
		cohorts.clear();
		cohort_params.clear();
		map<string, map<string, float>> all_results;
		for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it) {
			if (!has_element(it->second, search_term)) {
				cohorts[it->first] = filter_range_params;
				cohort_params[it->first] = (void *)&it->second;
			}
		}
		for (auto it = additional_cohorts.begin(); it != additional_cohorts.end(); ++it)
			if (cohorts.find(it->first) == cohorts.end())
				cohorts[it->first] = it->second;
			else
				MWARN("Cohort \"%s\" name already exists - skip additional\n");

		if (!cohorts.empty())
			all_results = booststrap_analyze(preds, y, *rep_ids, additional_info, cohorts,
				measures, &cohort_params, &measurements_params, fix_cohort_sample_incidence,
				preprocess_bin_scores, &roc_Params, sample_ratio, sample_per_pid, loopCnt, sample_seed);

		//now lets add each time window result:
		for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it) {
			if (has_element(it->second, search_term)) {
				vector<float> y_changed;
				map<string, vector<float>> cp_info;
				map<string, FilterCohortFunc> cohorts_t;
				map<string, void *> cohort_params_t;
				medial::process::make_sim_time_window(it->first, it->second,
					y, additional_info, y_changed, cp_info, cohorts_t, cohort_params_t);

				auto agg_res = booststrap_analyze(preds, y_changed, *rep_ids, cp_info, cohorts_t,
					measures, &cohort_params_t, &measurements_params, fix_cohort_sample_incidence,
					preprocess_bin_scores, &roc_Params, sample_ratio, sample_per_pid, loopCnt, sample_seed);
				if (!agg_res.empty()) // if the cohort is too small it does not return results
					all_results.insert(*agg_res.begin()); //has one key
			}
		}

		return all_results;
	}
}

map<string, map<string, float>> MedBootstrap::bootstrap_using_registry(MedFeatures &features_mat,
	const with_registry_args& args, map<int, map<string, map<string, float>>> *results_per_split) {
	MedBootstrap single_cohort = *this; //copy
	MedRegistry *registry = args.registry;
	string time_window_term = "Time-Window";
	map<string, map<string, float>> full_results;

	MLOG("Done reading %d record for registry\n", (int)registry->registry_records.size());
	MedSamplingYearly *sampler_year = args.sampler;
	//sample for each diffrent window_width in 3month res:
	unordered_map<int, map<string, vector<float>>> window_to_data;
	unordered_map<int, vector<MedSample>> window_to_smps;
	unordered_set<int> all_windows;
	unordered_map<string, int> cohort_to_time_res, cohort_to_time_filter_index;
	unordered_map<int, vector<int>> sorted_times;
	unordered_map<int, vector<vector<int>>> times_indexes;
	unordered_map<int, vector<MedRegistryRecord *>> pid_to_reg;
	MedFeatures *final_features = &features_mat;
	if (simTimeWindow) {
		for (size_t i = 0; i < registry->registry_records.size(); ++i)
			pid_to_reg[registry->registry_records[i].pid].push_back(&registry->registry_records[i]);
	}
	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it) {
		int from_w = 0, to_w = 0;
		bool found = false;
		for (size_t i = 0; i < it->second.size() && !found; ++i)
			if (it->second[i].param_name == time_window_term) {
				from_w = (int)it->second[i].min_range;
				to_w = (int)it->second[i].max_range;
				found = true;
				cohort_to_time_filter_index[it->first] = (int)i;
			}
		if (found) {
			//int time_res = ((to_w - from_w) / 90) * 90; //time window resultaion in 3 month's
			int time_res = to_w - from_w;
			all_windows.insert(time_res);
			cohort_to_time_res[it->first] = time_res;
		}
		else
			MWARN("Warning: No Time-Window in cohort \"%s\"\n", it->first.c_str());
	}

	MLOG("has %d time_window ranges\n", (int)all_windows.size());
	for (auto it = all_windows.begin(); it != all_windows.end(); ++it)
	{
		unordered_set<int> all_times;
		int time_res = *it;
		MLOG("Running Incidence calc for time Res %d\n", time_res);
		sampler_year->day_jump = time_res;
		get_data_for_filter(args.json_model, args.rep_path, single_cohort,
			registry->registry_records, *sampler_year, window_to_data[time_res], window_to_smps[time_res]);
		MLOG("Done preparing matrix of incidence for filtering with %d width...\n", time_res);
		if (args.do_kaplan_meir) {
			for (size_t i = 0; i < window_to_smps[time_res].size(); ++i) {
				int time_diff = (int)window_to_data[time_res][time_window_term][i];
				if (time_diff > time_res)
					time_diff = time_res;
				if (all_times.find(time_diff) == all_times.end()) {
					sorted_times[time_res].push_back(time_diff);
					all_times.insert(time_diff);
				}
			}
			sort(sorted_times[time_res].begin(), sorted_times[time_res].end());
			times_indexes[time_res].resize(sorted_times[time_res].size());
			for (size_t i = 0; i < window_to_smps[time_res].size(); ++i)
			{
				int time_diff = (int)window_to_data[time_res][time_window_term][i];
				if (time_diff > time_res)
					time_diff = time_res;
				int ind = medial::process::binary_search_index(sorted_times[time_res].data(),
					sorted_times[time_res].data() + sorted_times[time_res].size() - 1, time_diff);
				//if (ind < 0 || ind >= sorted_times.size())
				//	MTHROW_AND_ERR("BUG: bug in binary search\n");
				//skip cases after time window - will not occour, so skip
				if (window_to_smps[time_res][i].outcome <= 0 || (int)window_to_data[time_res][time_window_term][i] <= time_res)
					times_indexes[time_res][ind].push_back((int)i);
			}
		}
	}
	bool warn_shown = false;
	for (auto ii = filter_cohort.begin(); ii != filter_cohort.end(); ++ii) {
		final_features = &features_mat;
		single_cohort.filter_cohort.clear();
		single_cohort.filter_cohort[ii->first] = ii->second;
		Filter_Param time_filter;
		time_filter.param_name = "";
		for (size_t i = 0; i < ii->second.size(); ++i)
			if (ii->second[i].param_name == time_window_term)
			{
				time_filter = ii->second[i];
				break;
			}

		//save incidence_fix in ROC_Params
		if (cohort_to_time_res.find(ii->first) != cohort_to_time_res.end()) {
			int time_res = cohort_to_time_res[ii->first];
			int time_filter_index = cohort_to_time_filter_index[ii->first];
			double controls = 0, cases = 0, prob = 1;
			if (args.do_kaplan_meir) {
				double total_controls_all = 0;
				for (size_t sort_ind = 0; sort_ind < sorted_times[time_res].size(); ++sort_ind) {
					const vector<int> &index_order = times_indexes[time_res][sort_ind];
					ii->second[time_filter_index].max_range = (float)sorted_times[time_res][sort_ind];
					//update only controls count for time window - do for all time windows
					//to get total count from all time windows kaplna meir
					for (int i : index_order)
						if (window_to_smps[time_res][i].outcome <= 0 &&
							filter_range_params(window_to_data[time_res], (int)i, &ii->second))
							++total_controls_all;
				}

				for (size_t sort_ind = 0; sort_ind < sorted_times[time_res].size(); ++sort_ind) {
					const vector<int> &index_order = times_indexes[time_res][sort_ind];
					ii->second[time_filter_index].max_range = (float)sorted_times[time_res][sort_ind];
					for (int i : index_order) {
						if (filter_range_params(window_to_data[time_res], (int)i, &ii->second)) {
							//keep update kaplan meir in time point
							if (window_to_smps[time_res][i].outcome > 0)
								++cases;
							else
								++controls;
						}
					}
					//reset kaplan meir - flash last time prob
					if (!warn_shown && total_controls_all < 10) {
						MWARN("the kaplan_meir left with small amount of controls - "
							" try increasing the sampling / use smaller time window because the"
							" registry has not so long period of tracking patients\n");
						warn_shown = true;
					}
					if (total_controls_all > 0 || cases > 0)
						prob *= total_controls_all / (cases + total_controls_all);
					total_controls_all -= controls; //remove controls from current time-window - they are now censored
					controls = 0; cases = 0;
				}
				if (prob > 0 && prob < 1)
					single_cohort.roc_Params.incidence_fix = 1 - prob;
				MLOG("Incidence for %s cohort is %2.2f%% (kaplan meir)\n",
					ii->first.c_str(), single_cohort.roc_Params.incidence_fix * 100);
			}
			else {
				for (size_t i = 0; i < window_to_smps[time_res].size(); ++i)
					if (filter_range_params(window_to_data[time_res], (int)i, &ii->second))
						if (window_to_smps[time_res][i].outcome > 0)
							++cases;
						else
							++controls;
				if (cases > 0)
					single_cohort.roc_Params.incidence_fix = cases / (cases + controls);
				MLOG("Incidence for %s cohort is %2.2f%% (%d, %d)\n",
					ii->first.c_str(), single_cohort.roc_Params.incidence_fix * 100,
					(int)controls, (int)cases);
			}
		}
		else //no time window
			single_cohort.roc_Params.incidence_fix = 0;

		map<string, map<string, float>> part_res;
		MedFeatures sim_features;
		if (simTimeWindow) {
			MLOG("Censoring using MedRegistry for sim_time_window\n");
			sim_features = features_mat; //full copy
											//check for features_final - MedRegistry allowed if simTime - and filter
			vector<int> selected_rows;
			selected_rows.reserve(sim_features.samples.size());
			for (size_t i = 0; i < sim_features.samples.size(); ++i)
			{
				if (sim_features.samples[i].outcome <= 0 || time_filter.param_name.empty()) {
					selected_rows.push_back((int)i);
					continue;
				}
				int time_df = int(365 * (medial::repository::DateDiff(sim_features.samples[i].time,
					sim_features.samples[i].outcomeTime))); //time to turn into case
				if (time_df > time_filter.max_range) {
					//search for intersection:
					const vector<MedRegistryRecord *> &reg_records = pid_to_reg[sim_features.samples[i].id];
					bool is_legal = false;
					for (size_t k = 0; k < reg_records.size(); ++k)
					{
						if (reg_records[i]->registry_value > 0)
							continue;
						if (sim_features.samples[i].time >= reg_records[i]->min_allowed_date &&
							sim_features.samples[i].time <= reg_records[i]->max_allowed_date) {
							is_legal = true;
							break;
						}
					}
					if (is_legal)
						selected_rows.push_back((int)i); //add if legal
				}
				selected_rows.push_back((int)i);
			}
			medial::process::filter_row_indexes(sim_features, selected_rows);
			final_features = &sim_features;
		}

		map<int, map<string, map<string, float>>> part_results_per_split;
		map<int, map<string, map<string, float>>> *part_results_per_split_p = NULL;
		if (results_per_split != NULL)
			part_results_per_split_p = &part_results_per_split;
		part_res = single_cohort.bootstrap(*final_features, part_results_per_split_p);
		//Aggregate results:
		if (!part_res.empty()) {
			full_results[part_res.begin()->first] = part_res.begin()->second; //has one cohort
			if (results_per_split != NULL)
				for (auto ii = part_results_per_split.begin(); ii != part_results_per_split.end(); ++ii)
					(*results_per_split)[ii->first][ii->second.begin()->first] = ii->second.begin()->second;
		}
	}
	return full_results;
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
				split_data[jt->first][i] = jt->second[it->second[i]];
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
		if (features.samples[i].prediction.empty())
			MTHROW_AND_ERR("MedFeautres Prediciton is empty. need to run on MedFeatures with predictions\n");
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
	map<int, map<string, map<string, float>>> *results_per_split, with_registry_args *registry_args) {
	if (registry_args != NULL)
		return bootstrap_using_registry(features, *registry_args, results_per_split);

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
	map<string, vector<float>> &additional_info, map<int, map<string, map<string, float>>> *results_per_split, with_registry_args *registry_args) {
	MedFeatures features;
	features.data = additional_info;
	features.samples.reserve(samples.nSamples());
	for (size_t i = 0; i < samples.idSamples.size(); ++i)
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
			features.samples.push_back(samples.idSamples[i].samples[j]);
	features.init_pid_pos_len();
	return bootstrap(features, results_per_split, registry_args);
}

map<string, map<string, float>> MedBootstrap::bootstrap(MedSamples &samples, const string &rep_path, map<int, map<string, map<string, float>>> *results_per_split, with_registry_args *registry_args) {
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

	return bootstrap(mdl.features, results_per_split, registry_args);
}

map<string, map<string, float>> MedBootstrap::bootstrap(MedSamples &samples, MedPidRepository &rep, map<int, map<string, map<string, float>>> *results_per_split, with_registry_args *registry_args) {
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

	return bootstrap(mdl.features, results_per_split, registry_args);
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

void MedBootstrapResult::bootstrap(MedFeatures &features, map<int, map<string, map<string, float>>> *results_per_split, with_registry_args *registry_args) {
	bootstrap_results = bootstrap_params.bootstrap(features, results_per_split, registry_args);
}
void MedBootstrapResult::bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info, map<int, map<string, map<string, float>>> *results_per_split, with_registry_args *registry_args) {
	bootstrap_results = bootstrap_params.bootstrap(samples, additional_info, results_per_split, registry_args);
}
void MedBootstrapResult::bootstrap(MedSamples &samples, const string &rep_path, map<int, map<string, map<string, float>>> *results_per_split, with_registry_args *registry_args) {
	bootstrap_results = bootstrap_params.bootstrap(samples, rep_path, results_per_split, registry_args);
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

void MedBootstrapResult::write_results_to_text_file(const string &path, bool pivot_format, const string& run_id) {
	if (pivot_format)
		write_pivot_bootstrap_results(path, bootstrap_results, run_id);
	else
		write_bootstrap_results(path, bootstrap_results, run_id);
}

void MedBootstrapResult::read_results_to_text_file(const string &path, bool pivot_format) {
	if (pivot_format)
		read_pivot_bootstrap_results(path, bootstrap_results);
	else
		read_bootstrap_results(path, bootstrap_results);
}