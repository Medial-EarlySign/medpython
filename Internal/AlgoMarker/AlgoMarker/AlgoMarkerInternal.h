#pragma once

#include <string>
#include <InfraMed/InfraMed/InfraMed.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include "InputTesters.h"
#include "AlgoMarkerErr.h"

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
	//InputSanityTester ist;

	string name;
	string model_fname;
	string rep_fname;
	vector<int> pids;
	int model_end_stage = MED_MDL_END;
	bool model_init_done = false;

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
	int data_load_init() { rep.switch_to_in_mem_mode(); return 0; }

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
	int data_load_pid_sig(int pid, const char *sig_name, int *times, int n_times, float *vals, int n_vals) {
		int sid = rep.sigs.Name2Sid[string(sig_name)];
		if (sid < 0) return -1; // no such signal
		return rep.in_mem_rep.insertData(pid, sid, times, vals, n_times, n_vals);
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

	int get_raw_preds(int *_pids, int *times, float *preds) {

		try {

			try {
				// run model to calculate predictions
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
	}

	int get_pred(int *pid, int *time, float *pred) { return get_preds(pid, time, pred, 1); }


	//========================================================
	// Clearing - freeing mem
	//========================================================
	void clear() { pids.clear(); model.clear(); samples.clear(); rep.in_mem_rep.clear(); rep.clear(); }

	// clear_data() : leave model up, leave repository config up, but get rid of data and samples
	void clear_data() { samples.clear(); rep.in_mem_rep.clear(); }


	//========================================================
	// a few more needed APIs
	//========================================================
	const char *get_name() { return name.c_str(); }

	void write_features_mat(const string feat_mat) { model.write_feature_matrix(feat_mat); }

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

};

