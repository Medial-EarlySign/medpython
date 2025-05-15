#pragma once

#include <string>
#include <InfraMed/InfraMed/InfraMed.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/ExplainWrapper.h>
#include <MedStat/MedStat/MedBootstrap.h>
#include "InputTesters.h"
#include "AlgoMarkerErr.h"
#include <cmath>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

class Explainer_record_config : public SerializableObject {
public:
	string contributer_group_name = ""; ///< name of explainer group
	string signal_name = ""; ///< name of signal to fetch
	int max_count = 1; ///< limit of maximal last values
	int max_time_window = 0; ///< limit on maximal time before prediction
	int time_channel = 0; ///< time channel to filter
	int time_unit = MedTime::Days; ///< time unit for max_time_window
	int val_channel = 0; ///< val channel to filter sets
	vector<string> sets; ///< sets to filter categorical signal

	// Not to be initialized:
	vector<char> lut;
	int init(map<string, string>& mapper) {
		for (auto &it : mapper)
		{
			if (it.first == "contributer_group_name")
				contributer_group_name = it.second;
			else if (it.first == "signal_name")
				signal_name = it.second;
			else if (it.first == "max_count")
				max_count = med_stoi(it.second);
			else if (it.first == "max_time_window")
				max_time_window = med_stoi(it.second);
			else if (it.first == "time_channel")
				time_channel = med_stoi(it.second);
			else if (it.first == "val_channel")
				val_channel = med_stoi(it.second);
			else if (it.first == "time_unit")
				time_unit = med_time_converter.string_to_type(it.second);
			else if (it.first == "sets")
				boost::split(sets, it.second, boost::is_any_of(","));
			else
				HMTHROW_AND_ERR("Error in Explainer_record_config::init - unknown parameter \"%s\"\n",
					it.first.c_str());
		}
		if (contributer_group_name.empty())
			HMTHROW_AND_ERR("Error in Explainer_record_config::init - contributer_group_name must be given\n");
		if (signal_name.empty())
			HMTHROW_AND_ERR("Error in Explainer_record_config::init - signal_name must be given\n");
		return  0;
	}

	ADD_CLASS_NAME(Explainer_record_config)
		ADD_SERIALIZATION_FUNCS(contributer_group_name, signal_name, max_count, max_time_window, time_channel, time_unit, val_channel, sets)
};

class Explainer_description_config : public SerializableObject {
public:
	map<string, Explainer_record_config> records;

	void read_cfg_file(const string &file);

	int init(map<string, string>& mapper) {
		for (auto &it : mapper)
		{
			if (it.first == "records")
				read_cfg_file(it.second);
			else
				HMTHROW_AND_ERR("Error in Explainer_description_config::init - unknown parameter \"%s\"\n",
					it.first.c_str());
		}
		return  0;
	}

	ADD_CLASS_NAME(Explainer_description_config)
		ADD_SERIALIZATION_FUNCS(records)
};

class Explainer_parameters : public SerializableObject {
public:
	float max_threshold = 0; ///< control max threshold
	int num_groups = 3; ///< control how much binning to present
	bool use_perc = false; ///< control if binning on absolute value or on percentage
	Explainer_description_config cfg; ///< file to configure fetching signal to present
	unordered_set<string> ignore_groups_list; ///< name list of groups to alwaya ignore
	int total_max_reasons = -1; ///< if bigger than zero max limit for all reasons
	int total_max_pos_reasons = -1; ///< if bigger than zero max limit for pos reasons
	int total_max_neg_reasons = -1; ///< if bigger than zero max limit for neg reasons
	float threshold_abs = -1; ///< absolute thershold if bigger than 0
	float threshold_percentage = -1; ///< percentage thershold if bigger than 0
	vector<string> static_features_info; ///< config of information to fetch for every patient.

	// to be init before:
	string base_dir = "";

	int init(map<string, string>& mapper) {
		for (auto &it : mapper)
		{
			if (it.first == "max_threshold")
				max_threshold = med_stof(it.second);
			else if (it.first == "num_groups")
				num_groups = med_stoi(it.second);
			else if (it.first == "total_max_reasons")
				total_max_reasons = med_stoi(it.second);
			else if (it.first == "total_max_pos_reasons")
				total_max_pos_reasons = med_stoi(it.second);
			else if (it.first == "total_max_neg_reasons")
				total_max_neg_reasons = med_stoi(it.second);
			else if (it.first == "threshold_abs")
				threshold_abs = med_stof(it.second);
			else if (it.first == "threshold_percentage")
				threshold_percentage = med_stof(it.second);
			else if (it.first == "use_perc")
				use_perc = med_stoi(it.second) > 0;
			else if (it.first == "static_features_info")
				boost::split(static_features_info, it.second, boost::is_any_of(","));
			else if (it.first == "ignore_groups_list") {
				vector<string> tokens;
				boost::split(tokens, it.second, boost::is_any_of(","));
				ignore_groups_list.insert(tokens.begin(), tokens.end());
			}
			else if (it.first == "cfg") {
				if (it.second != "" && it.second[0] != '/' && it.second[0] != '\\' && !base_dir.empty())
					cfg.read_cfg_file(base_dir + path_sep() + it.second);
				else
					cfg.read_cfg_file(it.second);
			}
			else
				HMTHROW_AND_ERR("Error in Explainer_parameters::init - unknown parameter \"%s\"\n",
					it.first.c_str());
		}
		if (max_threshold < 0)
			HMTHROW_AND_ERR("Error in Explainer_parameters::init - max_threshold should be positive\n");

		return 0;
	}

	ADD_CLASS_NAME(Explainer_parameters)
		ADD_SERIALIZATION_FUNCS(max_threshold, num_groups, cfg, ignore_groups_list, total_max_reasons, total_max_pos_reasons, total_max_neg_reasons, threshold_abs, threshold_percentage, static_features_info)
};

//===============================================================================
// MedAlgoMarkerInternal - a mid-way API class : hiding all details of 
// implementation that are specific to the base classes (MedRepository, MedSamples, MedModel)
// that we use today.
// All functions assume c style to allow for easy export to C#/.NET
//===============================================================================
class  MedAlgoMarkerInternal {
private:
	// we force working ONLY using the API

	MedPidRepository rep;
	MedModel model;
	MedSamples samples;
	unordered_map<int, unordered_map<string, unordered_set<string>>> unknown_codes;
	Explainer_parameters explainer_params;
	//InputSanityTester ist;
	map<string, map<string, float>> mbr; ///< read bootstrap cohort, then measure and then value
	string default_threshold = ""; ///< deafult trehsold

	string name;
	string model_fname;
	string rep_fname;
	vector<int> pids;
	int model_end_stage = MED_MDL_END;
	bool model_init_done = false;
	bool model_rep_done = false;
public:

	MedPidRepository & get_rep() { return rep; }
	//========================================================
	// Initializations
	//========================================================

	// init name
	void set_name(const char *_name) { name = string(_name); }
	void set_model_end_stage(int _model_end_stage) { model_end_stage = _model_end_stage; };

	// init repository config
	int init_rep_config(const char *config_fname) {
		rep.switch_to_in_mem_mode();
		if (rep.MedRepository::init(string(config_fname)) < 0) return -1;

		return 0;
	}

	// set time_unit env for repositories and models
	int set_time_unit_env(int time_unit) {
		global_default_time_unit = time_unit;
		return 0;
	}

	// init pids
	void set_pids(int *_pids, int npids) { pids.clear(); pids.assign(_pids, _pids + npids); }

	// init rep , model , samples
	int init_rep_with_file_data(const char *_rep_fname) {
		rep.clear();
		rep_fname = string(_rep_fname);
		vector<string> sigs = {};
		return (rep.read_all(rep_fname, pids, sigs));
	}

	// init model
	int init_model_from_file(const char *_model_fname) { model.clear();	model.verbosity = 0; return (model.read_from_file(string(_model_fname))); }
	int model_check_required_signals() {
		int ret = 0;
		vector<string> req_sigs;
		model.get_required_signal_names(req_sigs);
		for (const auto& s : req_sigs)
			if (0 == rep.sigs.Name2Sid.count(s)) {
				ret = -1;
				fprintf(stderr, "ERROR: AM model requires signal '%s' but signal does not exist in AM repository .signals file\n", s.c_str());
			}
		return ret;
	}

	// init model for apply
	int init_model_for_apply() {
		global_logger.log(LOG_APP, LOG_DEF_LEVEL, "Init MedModel for Apply\n");
		model_init_done = true;
		return model.init_model_for_apply(rep, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_END);
	}

	void fit_model_to_rep() {
		model.fit_for_repository(rep);
	}

	int init_model_for_rep() {
		//global_logger.log(LOG_APP, LOG_DEF_LEVEL, "Init MedModel for Rep\n");
		if (!model_rep_done) {
			model_rep_done = true;
			return model.init_model_for_apply(rep, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_APPLY_FTR_PROCESSORS);
		}
		return 0;
	}

	unordered_map<string, unordered_set<string>> *get_unknown_codes(int pid) {
		return &unknown_codes[pid];
	}
	// init samples
	int init_samples(int *pids, int *times, int n_samples) { clear_samples(); int rc = insert_samples(pids, times, n_samples); samples.normalize(); return rc; }
	int init_samples(int pid, int time) { return init_samples(&pid, &time, 1); }  // single prediction point initiation 

	// init input_tester
	//int init_input_tester(const char *_fname) { return ist.read_config(string(_fname)); }

	void add_json_dict(json &js) { rep.dict.add_json(js); }

	bool model_initiated() { return model_init_done; }

	//========================================================
	// Loading data to rep
	//========================================================

	 // init loading : actions that must be taken BEFORE any loading starts
	int data_load_init() { unknown_codes.clear(); rep.switch_to_in_mem_mode(); return 0; }

	// load n_elems for a pid,sig
	int data_load_pid_sig(int pid, const char *sig_name, int *times, float *vals, int n_elems) {
		int sid = rep.sigs.Name2Sid[string(sig_name)];
		if (sid < 0) return -1; // no such signal
		int n_times = n_elems * rep.sigs.Sid2Info[sid].n_time_channels, n_vals = n_elems * rep.sigs.Sid2Info[sid].n_val_channels;
		if (times == NULL) n_times = 0;
		if (vals == NULL) n_vals = 0;
		return rep.in_mem_rep.insertData(pid, sid, times, vals, n_times, n_vals);
	}

	// load pid,sig with vectors of times and vals
	int data_load_pid_sig(int pid, const char *sig_name, int *times, int n_times, float *vals, int n_vals,
	map<pair<int, int>, pair<int, vector<char>>> *data = NULL) {
		int sid = rep.sigs.Name2Sid[string(sig_name)];
		if (sid < 0) return -1; // no such signal
		if (data == NULL)
			data = &rep.in_mem_rep.data;
		return rep.in_mem_rep.insertData_to_buffer(pid, sid, times, vals, n_times, n_vals, rep.sigs, *data);
	}

	// load a single element for a pid,sig
	int data_load_pid_sig(int pid, const char *sig_name, int *times, float *vals) { return data_load_pid_sig(pid, sig_name, times, vals, 1); }

	// end loading : actions that must be taken AFTER all loading was done, and BEFORE we calculate the predictions
	int data_load_end() { return rep.in_mem_rep.sortData(); }

	void get_rep_signals(unordered_set<string> &sigs)
	{
		for (auto &sig : rep.sigs.signals_names)
		{
			sigs.insert(sig);
		}
	}
	// returns the available signals

	//========================================================
	// Samples
	//========================================================

	// clear prediction points BEFORE a new set of predictions is done using the same instance
	void clear_samples() { samples.clear(); }

	// insert prediction points
	int insert_samples(int *pids, int *times, int n_samples) {
		for (int i = 0; i < n_samples; i++)
			samples.insertRec(pids[i], times[i]);
		return 0;
	}

	int insert_sample(int pid, int time) { return insert_samples(&pid, &time, 1); }

	//MedSample *get_sample(int idx) { if (idx>=0 && idx<samples.get_size()) return & }

	// normalize samples must be called after finishing inserting all samples.
	int normalize_samples() { samples.normalize(); return 0; }

	MedSamples *get_samples_ptr() { return &samples; }

	//========================================================
	// Calculate predictions
	//========================================================
	// note that if (_pids,times) are not sorted, they will be changed and sorted.
	int get_preds(int *_pids, int *times, float *preds, int n_samples) {

		// init_samples
		init_samples(_pids, times, n_samples);

		return get_raw_preds(_pids, times, preds);
	}

	int get_preds(int *_pids, int *times, float *preds, int n_samples,
		const vector<Effected_Field> &requested_fields, MedPidRepository *_rep=NULL) {

		// init_samples
		init_samples(_pids, times, n_samples);
		if (_rep == NULL)
			_rep = &this->rep;
		return get_raw_preds(_pids, times, preds, requested_fields, _rep);
	}

	int get_raw_preds(int *_pids, int *times, float *preds,
		const vector<Effected_Field> &requested_fields, MedPidRepository *_rep) {

		try {

			try {
				// run model to calculate predictions
				if (!samples.idSamples.empty())
					model.no_init_apply_partial(*_rep, samples, requested_fields);
			}
			catch (...) {
				fprintf(stderr, "Caught an exception in no_init_apply_partial\n");
				return -1;
			}

			// export pids, times and preds to c arrays
			int j = 0;
			if (preds != NULL) {
				for (auto& idSample : samples.idSamples)
					for (auto& sample : idSample.samples) {
						_pids[j] = sample.id;
						times[j] = sample.time;
						preds[j] = sample.prediction.size() > 0 ? sample.prediction[0] : (float)AM_UNDEFINED_VALUE; // This is Naive - but works for simple predictors giving the Raw score.
						j++;
					}
			}

			return 0;
		}
		catch (int &exception_code) {
			fprintf(stderr, "Caught an exception code: %d\n", exception_code);
			return -1; // exception_code;
		}
		catch (...) {
			fprintf(stderr, "Caught Something...\n");
			return -1;
		}
	}


	int get_raw_preds(int *_pids, int *times, float *preds) {

		try {

			try {
				// run model to calculate predictions
				if (!samples.idSamples.empty())
					if (model.no_init_apply(rep, samples, (MedModelStage)0, (MedModelStage)model_end_stage) < 0) {
						fprintf(stderr, "ERROR: MedAlgoMarkerInternal::get_preds FAILED.");
						return -1;
					}
			}
			catch (...) {
				fprintf(stderr, "Caught an exception in no_init_apply\n");
				return -1;
			}

			// export pids, times and preds to c arrays
			int j = 0;
			if (preds != NULL) {
				for (auto& idSample : samples.idSamples)
					for (auto& sample : idSample.samples) {
						_pids[j] = sample.id;
						times[j] = sample.time;
						preds[j] = sample.prediction.size() > 0 ? sample.prediction[0] : (float)AM_UNDEFINED_VALUE; // This is Naive - but works for simple predictors giving the Raw score.
						j++;
					}
			}

			return 0;
		}
		catch (int &exception_code) {
			fprintf(stderr, "Caught an exception code: %d\n", exception_code);
			return -1; // exception_code;
		}
		catch (...) {
			fprintf(stderr, "Caught Something...\n");
			return -1;
		}
	}

	int get_preds(MedSamples &_samples, float *preds) {

		samples = _samples;

		// run model to calculate predictions
		if (model.no_init_apply(rep, samples, (MedModelStage)0, (MedModelStage)model_end_stage) < 0) {
			fprintf(stderr, "ERROR: MedAlgoMarkerInternal::get_preds FAILED.");
			return -1;
		}

		// export pids, times and preds to c arrays
		int j = 0;
		for (auto& idSample : samples.idSamples)
			for (auto& sample : idSample.samples) {
				preds[j++] = sample.prediction[0]; // This is Naive - but works for simple predictors giving the Raw score.
			}
		return 0;
	}

	int get_pred(int *pid, int *time, float *pred) { return get_preds(pid, time, pred, 1); }


	//========================================================
	// Clearing - freeing mem
	//========================================================
	void clear() { unknown_codes.clear(); pids.clear(); model.clear(); samples.clear(); rep.in_mem_rep.clear(); rep.clear(); }

	// clear_data() : leave model up, leave repository config up, but get rid of data and samples
	void clear_data() {
		samples.clear(); rep.in_mem_rep.clear(); unknown_codes.clear();
	}


	//========================================================
	// a few more needed APIs
	//========================================================
	const char *get_name() { return name.c_str(); }

	void write_features_mat(const string &feat_mat) { model.write_feature_matrix(feat_mat); }
	void add_features_mat(const string &feat_mat) { model.write_feature_matrix(feat_mat, false, true); }

	void get_signal_structure(string &sig, int &n_time_channels, int &n_val_channels, int* &is_categ)
	{
		int sid = this->rep.sigs.sid(sig);
		if (sid <= 0) {
			n_time_channels = 0;
			n_val_channels = 0;
		}
		else {
			n_time_channels = this->rep.sigs.Sid2Info[sid].n_time_channels;
			n_val_channels = this->rep.sigs.Sid2Info[sid].n_val_channels;
			is_categ = &(this->rep.sigs.Sid2Info[sid].is_categorical_per_val_channel[0]);
		}
	}

	void model_apply_verbose(bool flag) {
		if ((model.verbosity > 0) ^ flag) {
			model.verbosity = int(flag);

			string full_log_format = "$timestamp\t$level\t$section\t%s";
			global_logger.init_format(LOG_APP, full_log_format);
			global_logger.init_format(LOG_DEF, full_log_format);
			global_logger.init_format(LOG_MED_MODEL, full_log_format);
			global_logger.init_format(LOG_MEDALGO, full_log_format);
			MLOG("Activated logging, Version Info:\n%s\n", medial::get_git_version().c_str());
		}
	}

	string model_version_info() const {
		return model.version_info;
	}

	void get_model_signals_info(vector<string> &sigs,
		unordered_map<string, vector<string>> &res_categ) const {
		model.get_required_signal_names(sigs);
		model.get_required_signal_categories(res_categ);
	}

	void get_explainer_params(Explainer_parameters &out) const {
		out = explainer_params;
	}

	void get_explainer_output_options(vector<string> &opts) {
		vector<const PostProcessor *> flat;
		for (const PostProcessor *pp : model.post_processors) {
			if (pp->processor_type == PostProcessorTypes::FTR_POSTPROCESS_MULTI)
			{
				const MultiPostProcessor *multi = static_cast<const MultiPostProcessor *>(pp);
				for (const PostProcessor *m_pp : multi->post_processors)
					flat.push_back(m_pp);
			}
			else
				flat.push_back(pp);
		}

		for (const PostProcessor *pp : flat)
		{
			const ModelExplainer *explainer_m = dynamic_cast<const ModelExplainer *>(pp);
			if (explainer_m != NULL) {
				for (const string &grp : explainer_m->processing.groupNames)
					opts.push_back(grp);
				break;
			}
		}
	}

	void set_explainer_params(const string &params, const string &base_dir) {
		explainer_params.base_dir = base_dir;
		explainer_params.init_from_string(params);
	}

	void set_threshold_leaflet(const string &init_string, const string &base_dir) {
		map<string, string> params;
		if (MedSerialize::init_map_from_string(init_string, params) < 0)
			MTHROW_AND_ERR("Error Init from String %s\n", init_string.c_str());
		string bt_file_path = "";
		map<string, string> rename_cohorts;
		for (const auto &it : params)
		{
			if (it.first == "bootstrap_file_path")
				bt_file_path = it.second;
			else if (it.first == "rename_cohorts") {
				vector<string> tokens;
				boost::split(tokens, it.second, boost::is_any_of("#"));
				for (const string &tk : tokens)
				{
					vector<string> src_target;
					boost::split(src_target, tk, boost::is_any_of("|"));
					if (src_target.size() != 2)
						MTHROW_AND_ERR("Error expecting 2 tokens, recieved \"%s\"\n", tk.c_str());
					mes_trim(src_target[1]);
					mes_trim(src_target[1]);
					rename_cohorts[src_target[0]] = src_target[1];
				}
			}
			else if (it.first == "default_threshold") {
				default_threshold = it.second;
				mes_trim(default_threshold);
			}
			else
				MTHROW_AND_ERR("Error unknown param %s\n", it.first.c_str());
		}
		if (bt_file_path.empty())
			MTHROW_AND_ERR("Error must provide bootstrap_file_path in THRESHOLD_LEAFLET\n");

		if (bt_file_path != "" && bt_file_path[0] != '/' && bt_file_path[0] != '\\' && !base_dir.empty())
			bt_file_path = base_dir + path_sep() + bt_file_path;

		if (default_threshold.empty())
			MTHROW_AND_ERR("Error - must have default_threshold\n");

		map<string, map<string, float>> mbr_before;
		read_pivot_bootstrap_results(bt_file_path, mbr_before);

		//commit rename:
		for (auto &it : mbr_before)
		{
			string cohort = it.first;
			if (rename_cohorts.find(cohort) != rename_cohorts.end())
				cohort = rename_cohorts[cohort];
			//Filter to take only "SCORE@" prefix - and SKIP missing values
			map<string, float> &filt = mbr[cohort];
			for (const auto &jt : it.second)
				if (boost::starts_with(jt.first, "SCORE@") && boost::ends_with(jt.first, "_Mean") && jt.second != MED_MAT_MISSING_VALUE)
					filt[jt.first.substr(6, jt.first.length() - 11)] = jt.second;
		}

		//Test default is OK:
		string err_c;
		fetch_threshold(default_threshold, err_c);
		if (!err_c.empty()) {
			vector<string> opts;
			fetch_all_thresholds(opts);
			for (const string & s : opts)
				MLOG("Option: \"%s\"\n", s.c_str());
			MTHROW_AND_ERR("Error default_threshold is invalid - please select one in format as COHORT$MEASURE_NUMERIC\n");
		}
	}

	bool has_threshold_settings() const {
		return !mbr.empty();
	}

	string get_default_threshold() const { return default_threshold; }

	void fetch_all_thresholds(vector<string> &opts) const {
		for (const auto &it : mbr)
		{
			for (const auto &jt : it.second)
			{
				string res = it.first + "$" + jt.first;
				opts.push_back(res);
			}
		}
	}

	float fetch_threshold(const string &threshold, string &err_msg) const {
		vector<string> tokens;
		err_msg = "";
		boost::split(tokens, threshold, boost::is_any_of("$"));
		if (tokens.size() != 2) {
			err_msg = "(" + to_string(AM_THRESHOLD_ERROR_NON_FATAL) + ")Error flag_threshold should contain $";
			return MED_MAT_MISSING_VALUE;
		}
		mes_trim(tokens[0]);
		mes_trim(tokens[1]);
		if (mbr.find(tokens[0]) == mbr.end()) {
			err_msg = "(" + to_string(AM_THRESHOLD_ERROR_NON_FATAL) + ")Error flag_threshold doesn't contain threshold settings for " + tokens[0];
			return MED_MAT_MISSING_VALUE;
		}
		const map<string, float> &fnd = mbr.at(tokens[0]);
		//Search numericly:
		vector<string> meas_tokens;
		boost::split(meas_tokens, tokens[1], boost::is_any_of("_"));
		if (meas_tokens.size() != 2) {
			err_msg = "(" + to_string(AM_THRESHOLD_ERROR_NON_FATAL) + ")Error flag_threshold doesn't should contain _ in the cutoff setting part";
			return MED_MAT_MISSING_VALUE;
		}
		float num_val;
		try {
			num_val = stof(meas_tokens[1]);
		}
		catch (...) {
			err_msg = "(" + to_string(AM_THRESHOLD_ERROR_NON_FATAL) + ")Error flag_threshold search cutoff isn't numeric";
			return MED_MAT_MISSING_VALUE;
		}

		float res = MED_MAT_MISSING_VALUE;
		for (const auto &jt : fnd)
		{
			string cand = jt.first;
			vector<string> cand_tokens;
			boost::split(cand_tokens, cand, boost::is_any_of("_"));
			if (cand_tokens.size() != 2)
				continue;
			if (cand_tokens[0] != meas_tokens[0])
				continue;
			//Need to compare numericaly: cand_tokens[1] == meas_tokens[1]
			float num_val_cmp;
			try {
				num_val_cmp = stof(cand_tokens[1]);
			}
			catch (...) {
				continue;
			}
			if (abs(num_val_cmp - num_val) <= 1e-6) {
				res = jt.second;
				break; //found
			}
		}

		if (res == MED_MAT_MISSING_VALUE)
			err_msg = "(" + to_string(AM_THRESHOLD_ERROR_NON_FATAL) + ")Error flag_threshold doesn't contain threshold for " + tokens[1];
		return res;
	}
};

MEDSERIALIZE_SUPPORT(Explainer_record_config)
MEDSERIALIZE_SUPPORT(Explainer_description_config)
MEDSERIALIZE_SUPPORT(Explainer_parameters)
