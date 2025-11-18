#include <unordered_set>
#include "MedBootstrap.h"
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include "Logger/Logger/Logger.h"
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL 

void MedBootstrap::get_cohort_from_arg(const string &single_cohort) {
	string cohort_line_with_name = single_cohort + "\t" + boost::replace_all_copy(single_cohort, " ", "\t");
	if (boost::starts_with(single_cohort, "MULTI"))
		cohort_line_with_name = boost::replace_all_copy(single_cohort, " ", "\t");
	parse_cohort_line(cohort_line_with_name);
}

void MedBootstrap::parse_cohort_line(const string &line) {
	if (line.find('\t') == string::npos)
		MWARN("Warning MedBootstrap::parse_cohort_line - line has no filters only name: \"%s\"\n",
			line.c_str());

	string cohort_name = line.substr(0, line.find('\t'));
	string cohort_definition = "";
	if (line.find('\t') != string::npos)
		cohort_definition = line.substr(line.find('\t') + 1);
	if (cohort_name != "MULTI") {
		vector<string> params;
		if (!cohort_definition.empty())
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
	censor_time_factor = 2;
	filter_cohort["All"] = {};
	simTimeWindow = false;
	is_binary_outcome = true;
	use_time_control_as_case = false;
	num_categories = 1;
	sort_preds_in_multicategory = false;
}

int MedBootstrap::init(map<string, string>& map) {
	//now read init_string to override default:
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		const string &param_name = boost::to_lower_copy(it->first);
		const string &param_value = it->second;
		//! [MedBootstrap::init]
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
		else if (param_name == "simtimewindow")
			simTimeWindow = stoi(param_value) > 0;
		else if (param_name == "censor_time_factor")
			censor_time_factor = stof(param_value);
		else if (param_name == "is_binary_outcome")
			is_binary_outcome = stoi(param_value) > 0;
		else if (param_name == "use_time_control_as_case")
			use_time_control_as_case = stoi(param_value) > 0;
		else if (param_name == "sort_preds_in_multicategory")
			sort_preds_in_multicategory = stoi(param_value) > 0;
		else if (param_name == "measurement") {
			measurements_with_params.clear();
			//parse and create this object tuples need to pass in {}: format measure_name|param_init_string;...
			vector<string> tokens;
			boost::split(tokens, param_value, boost::is_any_of(";"));
			for (size_t i = 0; i < tokens.size(); ++i)
			{
				vector<string> parts;
				boost::split(parts, tokens[i], boost::is_any_of("|"));
				if (parts.size() != 2)
					MTHROW_AND_ERR("Error MedBootstrap::init - expecting 2 tokens by | delimeter\n");
				MeasurmentFunctionType measure_type = measurement_function_name_to_type(parts[0]);
				MeasurementFunctions func = NULL;
				Measurement_Params *pr_m = NULL;
				string init_m_params = parts[1];
				switch (measure_type)
				{
				case MeasurmentFunctionType::calc_npos_nneg:
					func = calc_npos_nneg;
					break;
				case MeasurmentFunctionType::calc_only_auc:
					func = calc_only_auc;
					break;
				case MeasurmentFunctionType::calc_roc_measures_with_inc:
					func = calc_roc_measures_with_inc;
					pr_m = &roc_Params;
					break;
				case MeasurmentFunctionType::calc_multi_class:
					func = calc_multi_class;
					pr_m = &multiclass_params;
					multiclass_params.init_from_string(init_m_params);
					is_binary_outcome = false;
					break;
				case MeasurmentFunctionType::calc_kandel_tau:
					func = calc_kandel_tau;
					break;
				case MeasurmentFunctionType::calc_harrell_c_statistic:
					func = calc_harrell_c_statistic;
					break;
				case MeasurmentFunctionType::calc_regression:
					func = calc_regression;
					pr_m = &regression_params;
					regression_params.init_from_string(init_m_params);
					is_binary_outcome = false;
					break;
				default:
					MTHROW_AND_ERR("Error MedBootstrap::init unsupported measure %d\n", (int)measure_type);
				}

				measurements_with_params.push_back(pair<MeasurementFunctions, Measurement_Params*>(func, pr_m));
			}

		}
		//! [MedBootstrap::init]
		else
			MTHROW_AND_ERR("Unknown paramter \"%s\" for MedBootstrap::init\n", param_name.c_str());
	}


	return 0;
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
void init_model(MedModel &mdl, MedPidRepository& rep, const string &json_model,
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
	//read only dicts of rep:
	if (rep.init(rep_path) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", rep_path.c_str());
	mdl.fit_for_repository(rep);
	mdl.filter_rep_processors();

	unordered_set<string> req_names;
	mdl.get_required_signal_names(req_names);
	vector<string> sigs = { "BDATE", "GENDER" };
	for (string s : req_names)
		sigs.push_back(s);
	sort(sigs.begin(), sigs.end());
	auto it = unique(sigs.begin(), sigs.end());
	sigs.resize(std::distance(sigs.begin(), it));


	int curr_level = global_logger.levels.front();
	global_logger.init_all_levels(LOG_DEF_LEVEL);
	if (rep.read_all(rep_path, pids_to_take, sigs) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", rep_path.c_str());
	global_logger.init_all_levels(curr_level);
}
void get_data_for_filter(const string &json_model, const string &rep_path,
	MedBootstrap &single_cohort, const vector<MedRegistryRecord> &registry_records,
	MedSamplingStrategy &sampler, map<string, vector<float>> &data_for_filtering,
	vector<MedSample> &inc_smps, const LabelParams &par, const MedRegistry *censor_registry = NULL) {
	MedSamples inc_samples;
	MedPidRepository rep;
	vector<int> pids_to_take;
	unordered_set<int> seen_id;
	for (size_t i = 0; i < registry_records.size(); ++i)
		seen_id.insert(registry_records[i].pid);
	pids_to_take.insert(pids_to_take.end(), seen_id.begin(), seen_id.end());
	MedModel mdl;
	init_model(mdl, rep, json_model, rep_path, pids_to_take);

	sampler.init_sampler(rep);

	MedLabels labeler(par);
	labeler.prepare_from_registry(registry_records, censor_registry == NULL ? NULL : &censor_registry->registry_records);
	labeler.create_samples(&sampler, inc_samples);
	MLOG("Done sampling for incidence by year. has %d patients\n",
		(int)inc_samples.idSamples.size());

	if (mdl.learn(rep, &inc_samples, MedModelStage::MED_MDL_LEARN_REP_PROCESSORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS) < 0)
		MTHROW_AND_ERR("Error creating age,gender for samples\n");

	vector<float> preds, y;
	vector<int> pids;
	unordered_map<int, vector<int>> splits_inds;
	for (size_t i = 0; i < mdl.features.samples.size(); ++i)
		mdl.features.samples[i].prediction.resize(1, 0); //set to 0 that won't throw error in prepare

	vector<int> preds_order;
	single_cohort.prepare_bootstrap(mdl.features, preds, y, pids, data_for_filtering, preds_order);
	inc_smps.swap(mdl.features.samples);
}

void medial::process::make_sim_time_window(
	const string &cohort_name, 
	const vector<Filter_Param> &filter_p,
	const vector<float> &y, 
	const map<string, vector<float>> &additional_info,
	vector<float> &y_changed, 
	map<string, vector<float>> &cp_info,
	map<string, FilterCohortFunc> &cohorts_t, 
	map<string, void *> &cohort_params_t, 
	float censor_time_factor) 
{
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
	
	if (time_filter.param_name == "")
		MTHROW_AND_ERR("Cohort has no Time-Window, can't use make_sim_time_window\n");
	//update: y_changed based on y
	cp_info = additional_info; //a copy - i will change those each time
	y_changed.insert(y_changed.end(), y.begin(), y.end());
	//update: cp_info["Time-Window"] based on additional_info["Time-Window"] for cases not in time window
	int cases_censored = 0, cases_changed_to_controls = 0, cases_left_as_cases = 0;
	int max_range = (int)time_filter.max_range;
	for (size_t i = 0; i < y_changed.size(); ++i)
	{
		if (y[i] > 0) {
			if (!filter_range_param(additional_info, (int)i, &time_filter)) {
				//when censor_time_factor=2
				// cases which are long before the outcome (>2*max_range) are considered as controls:

				// ----------------max_range*2--------max_range---------min_range---event---------------
				// -------------------^-------------------^----------------^----------^-----------------
				// ------control------|------censor-------|------case------|--------censor--------------

				if (additional_info.at(search_term)[i] > max_range * censor_time_factor) {
					y_changed[i] = 0;
					cp_info[search_term][i] = time_filter.min_range; // bogus time - wont filter because in time range
					cases_changed_to_controls++;
				}
				else
					cases_censored++;
			}
			else cases_left_as_cases++;
		}
	}
	MLOG("make_sim_time_window for cohort [%s] censored %d cases with %d<gap<%d and changed %d cases to controls with %d<gap, left with %d cases\n",
		cohort_name.c_str(), 
		cases_censored, 
		max_range, 
		max_range * (int)censor_time_factor, 
		cases_changed_to_controls, 
		max_range * (int)censor_time_factor, 
		cases_left_as_cases
	);
}

map<string, map<string, float>> MedBootstrap::bootstrap_base(const vector<float> &preds, const vector<int> &preds_order, const vector<float> &y, const vector<int> &pids,
	const vector<float> *weights, const map<string, vector<float>> &additional_info) {
	if (weights != NULL && !weights->empty() && weights->size() != y.size())
		MTHROW_AND_ERR("Bootstrap weights error - samples is of size %zu, weights.size()=%zu\n",
			y.size(), weights->size());

	map<string, FilterCohortFunc> cohorts;
	map<string, void *> cohort_params;
	if (measurements_with_params.empty()) { //if not touched(than empty) - this is the default!
		measurements_with_params.resize(1);
		measurements_with_params[0] = pair<MeasurementFunctions, Measurement_Params*>(calc_roc_measures_with_inc, &roc_Params);
	}
	vector<MeasurementFunctions> measures((int)measurements_with_params.size());
	vector<Measurement_Params *> measurements_params((int)measurements_with_params.size());
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
			MWARN("Cohort \"%s\" name already exists - skip additional\n", it->first.c_str());
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
		return booststrap_analyze(preds, preds_order, y, weights, *rep_ids, additional_info, cohorts,
			measures, &cohort_params, &measurements_params, fix_cohort_sample_incidence,
			preprocess_bin_scores, &roc_Params, sample_ratio, sample_per_pid, loopCnt, sample_seed, is_binary_outcome);
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
				MWARN("Cohort \"%s\" name already exists - skip additional\n", it->first.c_str());

		if (!cohorts.empty())
			all_results = booststrap_analyze(preds, preds_order, y, weights, *rep_ids, additional_info, cohorts,
				measures, &cohort_params, &measurements_params, fix_cohort_sample_incidence,
				preprocess_bin_scores, &roc_Params, sample_ratio, sample_per_pid, loopCnt, sample_seed, is_binary_outcome);

		//now lets add each time window result:
		for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it) {
			if (has_element(it->second, search_term)) {
				vector<float> y_changed;
				map<string, vector<float>> cp_info;
				map<string, FilterCohortFunc> cohorts_t;
				map<string, void *> cohort_params_t;
				medial::process::make_sim_time_window(it->first, it->second,
					y, additional_info, y_changed, cp_info, cohorts_t, cohort_params_t, censor_time_factor);

				auto agg_res = booststrap_analyze(preds, preds_order, y_changed, weights, *rep_ids, cp_info, cohorts_t,
					measures, &cohort_params_t, &measurements_params, fix_cohort_sample_incidence,
					preprocess_bin_scores, &roc_Params, sample_ratio, sample_per_pid, loopCnt, sample_seed, is_binary_outcome);
				if (!agg_res.empty()) // if the cohort is too small it does not return results
					all_results.insert(*agg_res.begin()); //has one key
			}
		}

		return all_results;
	}
}

struct pair_hash {
	template <class T1, class T2>
	std::size_t operator () (const std::pair<T1, T2> &p) const {
		auto h1 = std::hash<T1>{}(p.first);
		auto h2 = std::hash<T2>{}(p.second);

		// Mainly for demonstration purposes, i.e. works but is overly simple
		// In the real world, use sth. like boost.hash_combine
		return h1 ^ h2;
	}
};
map<string, map<string, float>> MedBootstrap::bootstrap_using_registry(const MedFeatures &features_mat,
	const with_registry_args& args, map<int, map<string, map<string, float>>> *results_per_split) {
	MedBootstrap single_cohort = *this; //copy
	MedRegistry *registry = args.registry;
	string time_window_term = "Time-Window";
	map<string, map<string, float>> full_results;

	MLOG("Done reading %d record for registry\n", (int)registry->registry_records.size());
	MedSamplingYearly *sampler_year = args.sampler;
	//sample for each diffrent window_width in 3month res:
	unordered_map<pair<int, int>, map<string, vector<float>>, pair_hash> window_to_data;
	unordered_map<pair<int, int>, vector<MedSample>, pair_hash> window_to_smps;
	unordered_set<pair<int, int>, pair_hash> all_windows;
	unordered_map<string, pair<int, int>> cohort_to_time_res;
	unordered_map<string, int> cohort_to_time_filter_index;
	unordered_map<int, vector<MedRegistryRecord *>> pid_to_reg, pid_to_censor;
	const MedFeatures *final_features = &features_mat;
	if (simTimeWindow) {
		for (size_t i = 0; i < registry->registry_records.size(); ++i)
			pid_to_reg[registry->registry_records[i].pid].push_back(&registry->registry_records[i]);
		if (args.registry_censor != NULL)
			for (size_t i = 0; i < args.registry_censor->registry_records.size(); ++i)
				pid_to_censor[args.registry_censor->registry_records[i].pid].push_back(&(args.registry_censor->registry_records)[i]);
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
			//int time_res = to_w - from_w;
			all_windows.insert(pair<int, int>(from_w, to_w));
			cohort_to_time_res[it->first] = pair<int, int>(from_w, to_w);
		}
		else
			MWARN("Warning: No Time-Window in cohort \"%s\"\n", it->first.c_str());
	}

	MLOG("has %d time_window ranges\n", (int)all_windows.size());
	for (auto it = all_windows.begin(); it != all_windows.end(); ++it)
	{
		unordered_set<int> all_times;
		pair<int, int> time_res = *it;
		MLOG("Running Incidence calc for time Res %d to %d\n", time_res.first, time_res.second);
		sampler_year->time_jump = time_res.second - time_res.first;
		sampler_year->time_jump_unit = global_default_windows_time_unit;
		LabelParams lbl_params = args.labeling_params;
		lbl_params.time_from = time_res.first;
		lbl_params.time_to = time_res.second;
		get_data_for_filter(args.json_model, args.rep_path, single_cohort,
			registry->registry_records, *sampler_year, window_to_data[time_res], window_to_smps[time_res],
			lbl_params, args.registry_censor);
		MLOG("Done preparing matrix of incidence for filtering with [%d, %d] window...\n", time_res.first, time_res.second);
	}
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
		vector<Filter_Param> without_time_filter = ii->second;
		without_time_filter.erase(without_time_filter.begin() + cohort_to_time_filter_index[ii->first]);

		//save incidence_fix in ROC_Params
		if (cohort_to_time_res.find(ii->first) != cohort_to_time_res.end()) {
			const pair<int, int> &time_res = cohort_to_time_res[ii->first];
			double controls = 0, cases = 0;
			if (args.do_kaplan_meir) {
				vector<int> sel_idx;
				//update only controls count for time window - do for all time windows
				//to get total count from all time windows kaplna meir
				for (int i = 0; i < window_to_smps[time_res].size(); ++i)
					if (filter_range_params(window_to_data[time_res], (int)i, &without_time_filter))
						sel_idx.push_back(i);
				double prob = medial::stats::kaplan_meir_on_samples(window_to_smps[time_res], features_mat.time_unit, time_res.second - time_res.first,
					&sel_idx);

				if (prob > 0 && prob < 1)
					single_cohort.roc_Params.incidence_fix = prob;
				else
					MWARN("In MedBootstrap::bootstrap_using_registry - Has no cases/controls in Incidence - please check registry or sampling args\n");

				MLOG("Incidence for %s cohort is %2.4f%% (kaplan meir) - total in cohort: (%zu)\n",
					ii->first.c_str(), single_cohort.roc_Params.incidence_fix * 100, sel_idx.size());
			}
			else {
				for (size_t i = 0; i < window_to_smps[time_res].size(); ++i)
					if (filter_range_params(window_to_data[time_res], (int)i, &ii->second)) {
						if (window_to_smps[time_res][i].outcome > 0)
							++cases;
						else
							++controls;
					}
				if (cases > 0)
					single_cohort.roc_Params.incidence_fix = cases / (cases + controls);
				else
					MWARN("In MedBootstrap::bootstrap_using_registry - Has no cases in Incidence - please check registry or sampling args\n");
				MLOG("Incidence for %s cohort is %2.2f%% (%d, %d)\n",
					ii->first.c_str(), single_cohort.roc_Params.incidence_fix * 100,
					(int)controls, (int)cases);
			}
		}
		else //no time window
			single_cohort.roc_Params.incidence_fix = 0;

		map<string, map<string, float>> part_res;
		MedFeatures sim_features;
		vector<MedRegistryRecord *> empty_arr;
		if (simTimeWindow) {
			MLOG("Censoring using MedRegistry for sim_time_window\n");
			sim_features = features_mat; //full copy
											//check for features_final - MedRegistry allowed if simTime - and filter
			vector<int> selected_rows;
			selected_rows.reserve(sim_features.samples.size());
			int no_censoring = 0;
			for (size_t i = 0; i < sim_features.samples.size(); ++i)
			{
				if (sim_features.samples[i].outcome <= 0 || time_filter.param_name.empty()) {
					selected_rows.push_back((int)i);
					continue;
				}
				int time_df = (med_time_converter.convert_times(sim_features.time_unit, MedTime::Days, sim_features.samples[i].outcomeTime)
					- med_time_converter.convert_times(sim_features.time_unit, MedTime::Days, sim_features.samples[i].time));
				if (time_df > time_filter.max_range) {
					//search for intersection:
					const vector<MedRegistryRecord *> *reg_censor = &empty_arr;
					if (pid_to_censor.find(sim_features.samples[i].id) != pid_to_censor.end())
						reg_censor = &pid_to_censor[sim_features.samples[i].id];
					else
						++no_censoring;
					bool is_legal = reg_censor->empty();
					for (size_t k = 0; k < reg_censor->size(); ++k)
					{
						int diff_to_allowed = med_time_converter.convert_days(global_default_windows_time_unit,
							int(365 * (medial::repository::DateDiff(sim_features.samples[i].time, reg_censor->at(k)->end_date))));

						if (diff_to_allowed >= time_filter.max_range && sim_features.samples[i].time >= reg_censor->at(k)->start_date) {
							is_legal = true;
							break;
						}
					}
					if (is_legal)
						selected_rows.push_back((int)i); //add if legal
				}
				else
					selected_rows.push_back((int)i);
			}
			medial::process::filter_row_indexes(sim_features, selected_rows);
			final_features = &sim_features;
			if (no_censoring > 0)
				MWARN("Warning in bootstarp has %d pids without censor registry\n", no_censoring);
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

void MedBootstrap::add_splits_results(const vector<float> &preds, vector<int> &preds_order, const vector<float> &y,
	const vector<int> &pids, const vector<float> *weights, const map<string, vector<float>> &data,
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
		vector<float> split_weights;

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
		if (weights != NULL && !weights->empty()) {
			split_weights.resize((int)it->second.size());
			for (size_t i = 0; i < split_preds.size(); ++i)
				split_weights[i] = weights->at(it->second[i]);
		}
		results_per_split[split_id] =
			bootstrap_base(split_preds, preds_order, split_y, split_pids, &split_weights, split_data);
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

void MedBootstrap::sort_index_only(const vector<float> &vec, vector<int>::iterator ind_start, vector<int>::iterator ind_end)
{
	iota(ind_start, ind_end, 0);
	stable_sort(ind_start, ind_end, [&vec](size_t i1, size_t i2) {return vec[i1] > vec[i2]; });
}

void MedBootstrap::prepare_bootstrap(const MedFeatures &features, vector<float> &preds, vector<float> &y, vector<int> &pids,
	map<string, vector<float>> &final_additional_info, vector<int> &preds_order, unordered_map<int, vector<int>> *splits_inds) {
	if (features.samples.empty() || features.samples[0].prediction.empty())
		MTHROW_AND_ERR("prepare_bootstrap - sample size %d or number of predictions %d are not valid \n", (int)features.samples.size(), (int)features.samples[0].prediction.size());
	num_categories = features.samples[0].prediction.size();

	preds.resize((int)features.samples.size()*num_categories);
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

	preds_order.assign((int)features.samples.size()*num_categories, 0);

	for (size_t i = 0; i < features.samples.size(); ++i)
	{
		pids[i] = features.samples[i].id;
		y[i] = features.samples[i].outcome;
		if (features.samples[i].prediction.empty())
			MTHROW_AND_ERR("MedFeautres Prediciton is empty. need to run on MedFeatures with predictions\n");

		for (size_t j = 0; j < num_categories; j++)
		{
			preds[i*num_categories + j] = features.samples[i].prediction[j];
			preds_order[i*num_categories + j] = (int)j;
		}
		if (sort_preds_in_multicategory)
		{
			sort_index_only(features.samples[i].prediction, preds_order.begin() + i * num_categories, preds_order.begin() + (i + 1) * num_categories);

		}

		if (uses_time_window) {
			int diff_days = (tm.convert_times(features.time_unit, MedTime::Days, features.samples[i].outcomeTime)
				- tm.convert_times(features.time_unit, MedTime::Days, features.samples[i].time));
			final_additional_info["Time-Window"][i] = (float)diff_days;
			final_additional_info["Label"][i] = !use_time_control_as_case ? y[i] : 1;
		}
		if (splits_inds != NULL)
			(*splits_inds)[features.samples[i].split].push_back((int)i);
	}
}
map<string, map<string, float>> MedBootstrap::bootstrap(const MedFeatures &features,
	map<int, map<string, map<string, float>>> *results_per_split, with_registry_args *registry_args) {
	if (registry_args != NULL)
		return bootstrap_using_registry(features, *registry_args, results_per_split);

	vector<float> preds, y;
	vector<int> pids;
	vector<int> preds_order;
	map<string, vector<float>> data;
	unordered_map<int, vector<int>> splits_inds;

	if (results_per_split != NULL)
		prepare_bootstrap(features, preds, y, pids, data, preds_order, &splits_inds);
	else
		prepare_bootstrap(features, preds, y, pids, data, preds_order);

	//check we have all signals ececpt Time,Label (will be completed in prepare):
	vector<string> all_names_ls;
	all_names_ls.reserve(data.size());
	for (auto &it : data)
		all_names_ls.push_back(it.first);
	unordered_set<string> print_mp;
	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it)
		for (Filter_Param &fp : it->second) {
			if (fp.param_name == "Time-Window" || fp.param_name == "Label")
				continue;
			if (data.find(fp.param_name) == data.end()) {
				//try fix:
				int fn_pos = find_in_feature_names(all_names_ls, fp.param_name, false);
				if (fn_pos < 0)
					MTHROW_AND_ERR("ERROR in MedBootstrap::bootstrap - missing "
						"filter_cohort parameter \"%s\" in cohort \"%s\" in input features.\n"
						"Please provide the feature in the input for filtering or remove cohort filter\n",
						fp.param_name.c_str(), it->first.c_str());
				//passed:
				string found_nm = all_names_ls[fn_pos];
				fp.param_name = found_nm; //change to name
				if (print_mp.find(fp.param_name) == print_mp.end()) {
					print_mp.insert(fp.param_name);
					MLOG("Info::bootstrap Mapped %s => %s\n", fp.param_name.c_str(), found_nm.c_str());
				}
			}
		}

	if (results_per_split != NULL)
		add_splits_results(preds, preds_order, y, pids, &features.weights, data, splits_inds, *results_per_split);

	return bootstrap_base(preds, preds_order, y, pids, &features.weights, data);
}

void MedBootstrap::prepare_bootstrap(MedSamples &samples, map<string, vector<float>> &additional_info, vector<float> &preds, vector<float> &y, vector<int> &pids, vector<int> &preds_order,
	unordered_map<int, vector<int>> *splits_inds) {
	pids.resize(samples.nSamples());
	samples.get_y(y);
	samples.get_preds(preds);
	if (y.empty() || preds.empty())
		MTHROW_AND_ERR("prepare_bootstrap - sample size %d or number of predictions %d are not valid \n", (int)y.size(), (int)preds.size());
	num_categories = preds.size() / y.size();

	preds_order.assign((int)preds.size()*num_categories, 0);

	int i = 0;
	if (num_categories > 1)
	{
		for (auto& idSamples : samples.idSamples)
		{
			for (auto& samples : idSamples.samples)
			{
				if (sort_preds_in_multicategory)
				{
					sort_index_only(samples.prediction, preds_order.begin() + i * num_categories, preds_order.begin() + (i + 1) * num_categories);
				}
				else
					iota(preds_order.begin() + i * num_categories, preds_order.begin() + (i + 1) * num_categories, 0);
				i++;
			}
		}
	}

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
				int diff_days = (tm.convert_times(samples.time_unit, MedTime::Days,
					samples.idSamples[i].samples[j].outcomeTime)
					- tm.convert_times(samples.time_unit, MedTime::Days, samples.idSamples[i].samples[j].time));
				additional_info["Time-Window"][c] = (float)diff_days;
				additional_info["Label"][c] = !use_time_control_as_case ? samples.idSamples[i].samples[j].outcome : 1;
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
	vector<string> sigs = { "BDATE", "GENDER", "TRAIN" };
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
	vector<string> sigs = { "BDATE", "GENDER", "TRAIN" };
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
			int diff_days = (tm.convert_times(samples.time_unit, MedTime::Days,
				samples.idSamples[i].samples[j].outcomeTime)
				- tm.convert_times(samples.time_unit, MedTime::Days, samples.idSamples[i].samples[j].time));
			//check if valid for time_window:
			if (samples.idSamples[i].samples[j].outcome > 0) {
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
			MedIdSamples smp_pid(samples.idSamples[i].id);
			smp_pid.split = samples.idSamples[i].split;
			new_samples.idSamples.push_back(smp_pid);
			new_samples.idSamples.back().samples.push_back(samples.idSamples[i].samples[keep_control_index]);
		}
		if (keep_case_index >= 0) {
			if (keep_control_index < 0) {
				MedIdSamples smp_pid(samples.idSamples[i].id);
				smp_pid.split = samples.idSamples[i].split;
				new_samples.idSamples.push_back(smp_pid);
			}
			new_samples.idSamples.back().samples.push_back(samples.idSamples[i].samples[keep_case_index]);
		}
	}
	//new_samples.sort_by_id_date(); //not needed
}

unordered_map<string, MeasurmentFunctionType> MedBootstrap::measurement_function_name_map = {
	{ "calc_npos_nneg",MeasurmentFunctionType::calc_npos_nneg },
{ "calc_only_auc",MeasurmentFunctionType::calc_only_auc },
{ "calc_roc_measures_with_inc",MeasurmentFunctionType::calc_roc_measures_with_inc },
{ "calc_multi_class",MeasurmentFunctionType::calc_multi_class },
{ "calc_kandel_tau", MeasurmentFunctionType::calc_kandel_tau },
{ "calc_harrell_c_statistic", MeasurmentFunctionType::calc_harrell_c_statistic },
{ "calc_regression", MeasurmentFunctionType::calc_regression }
};

MeasurmentFunctionType MedBootstrap::measurement_function_name_to_type(const string& measurement_function_name) {
	if (measurement_function_name_map.find(measurement_function_name) == measurement_function_name_map.end()) {
		string options = medial::io::get_list(measurement_function_name_map, "\n");
		MTHROW_AND_ERR("measurement_function_name_to_type: unknown name %s\nOptions:%s\n",
			measurement_function_name.c_str(), options.c_str());
	}
	return measurement_function_name_map[measurement_function_name];
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
			int diff_days = (tm.convert_times(features.time_unit, MedTime::Days, features.samples[i].outcomeTime)
				- tm.convert_times(features.time_unit, MedTime::Days, features.samples[i].time));
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

void map_bootstrap_names(vector<Filter_Param> &filters_test,
	const map<string, vector<float>> &data) {
	bool all_valid = true;
	vector<string> all_names_ls;
	for (auto &it : data)
		all_names_ls.push_back(it.first);

	for (Filter_Param param_name : filters_test)
		if (data.find(param_name.param_name) == data.end()) {
			//try and fix first:
			int fn_pos = find_in_feature_names(all_names_ls, param_name.param_name, false);
			if (fn_pos < 0) {
				all_valid = false;
				MERR("ERROR:: Wrong use in \"%s\" as filter params. the parameter is missing\n", param_name.param_name.c_str());
			}
			else {
				//fix name:
				string found_nm = all_names_ls[fn_pos];
				MLOG("Mapped %s => %s\n", param_name.param_name.c_str(), found_nm.c_str());
				for (Filter_Param &pp : filters_test)
					if (pp.param_name == param_name.param_name)
						pp.param_name = found_nm;
			}
		}

	if (!all_valid) {
		MLOG("Feature Names availible for cohort filtering:\n");
		for (auto it = data.begin(); it != data.end(); ++it)
			MLOG("%s\n", it->first.c_str());
		MLOG("Time-Window\n");
		MLOG("\n");

		MTHROW_AND_ERR("Cohort file has wrong paramter names. look above for all avaible params\n");
	}
}


void MedBootstrap::filter_bootstrap_cohort(MedFeatures &features, const string &bt_cohort) {
	//create bt object with cohort filter
	MedBootstrap mbr;
	mbr.filter_cohort.clear();
	string cohort_format = "bs\t" + bt_cohort;
	mbr.parse_cohort_line(cohort_format);

	vector<float> preds, labelsOrig;
	vector<int> pidsOrig;
	bool has_no_preds = false;
	// Create Dummy predictions (needed for prepare)
	for (MedSample &sample : features.samples)
		if (sample.prediction.empty()) {
			sample.prediction = { 0.0 };
			has_no_preds = true;
		}

	vector<int> preds_order;
	map<string, vector<float>> bt_data;
	mbr.prepare_bootstrap(features, preds, labelsOrig, pidsOrig, bt_data, preds_order);
	//check that all params are known:
	vector<Filter_Param> &filters_test = mbr.filter_cohort["bs"];
	map_bootstrap_names(filters_test, bt_data);
	
	//commit filter on: test_bt_features and store to curr_samples
	vector<int> sel;
	for (int i = 0; i < features.samples.size(); ++i)
		if (filter_range_params(bt_data, i, &mbr.filter_cohort["bs"]))
			sel.push_back(i);

	medial::process::filter_row_indexes(features, sel);
	//clear predictions:
	if (has_no_preds)
		for (size_t i = 0; i < features.samples.size(); ++i)
			features.samples[i].prediction.clear();
}

void MedBootstrap::filter_bootstrap_cohort(MedPidRepository &bt_repository, MedModel &bt_filters,
	MedSamples &curr_samples, const string &bt_cohort) {

	if (bt_filters.learn(bt_repository, &curr_samples,
		MedModelStage::MED_MDL_LEARN_REP_PROCESSORS, MedModelStage::MED_MDL_APPLY_FTR_PROCESSORS) < 0)
		MTHROW_AND_ERR("Error creating bootstrap features for test samples\n");
	MedFeatures &test_bt_features = bt_filters.features;
	//create bt object with cohort filter
	filter_bootstrap_cohort(test_bt_features, bt_cohort);
	//commit test_bt_features to curr_samples
	curr_samples.import_from_sample_vec(test_bt_features.samples);
}

void MedBootstrap::filter_bootstrap_cohort(const string &rep, const string &bt_json,
	MedSamples &curr_samples, const string &bt_cohort) {
	MedModel bt_filters;
	if (!bt_json.empty())
		bt_filters.init_from_json_file(bt_json);
	//check if has Age:
	bool fnd_age = false, fnd_gender = false;
	for (const FeatureGenerator *f : bt_filters.generators) {
		if (f->generator_type == FTR_GEN_AGE)
			fnd_age = true;
		if (f->generator_type == FTR_GEN_GENDER)
			fnd_gender = true;
		if (fnd_age && fnd_gender)
			break;
	}
	if (!fnd_age)
		bt_filters.add_age();
	if (!fnd_gender)
		bt_filters.add_gender();

	MedPidRepository bt_repository;
	bt_filters.load_repository(rep, bt_repository, true);

	filter_bootstrap_cohort(bt_repository, bt_filters, curr_samples, bt_cohort);
}