#pragma once

#include <string>
#include <InfraMed/InfraMed/InfraMed.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>

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

	string name;
	string model_fname;
	string rep_fname;
	vector<int> pids;

public:
	//========================================================
	// Initializations
	//========================================================

	// init name
	void set_name(const char *_name) { name = string(_name); }

	// init repository config
	int init_rep_config(const char *config_fname) { if (rep.init(string(config_fname)) < 0) return -1; rep.switch_to_in_mem_mode(); return 0; }

	// init pids
	void set_pids(int *_pids, int npids) { pids.clear(); pids.assign(_pids, _pids + npids); }

	// init rep , model , samples
	int init_rep_with_file_data(const char *_rep_fname) {
		rep.clear();
		rep_fname = string(_rep_fname);
		vector<string> sigs ={};
		return (rep.read_all(rep_fname, pids, sigs));
	}

	// init model
	int init_model_from_file(const char *_model_fname) { model.clear();	model.verbosity = 0; return (model.read_from_file(string(_model_fname))); }


	// init samples
	int init_samples(int *pids, int *times, int n_samples) { clear_samples(); int rc = insert_samples(pids, times, n_samples); samples.normalize(); return rc; }
	int init_samples(int pid, int time) { return init_samples(&pid, &time, 1); }  // single prediction point initiation 


   //========================================================
   // Loading data to rep
   //========================================================
	
    // init loading : actions that must be taken BEFORE any loading starts
	int data_load_init() { rep.switch_to_in_mem_mode(); return 0; }

	// load n_elems for a pid,sig
	int data_load_pid_sig(int pid, const char *sig_name, int *times, float *vals, int n_elems) {
		int sid = rep.sigs.Name2Sid[string(sig_name)];
		if (sid < 0) return -1; // no such signal
		int n_times = n_elems * rep.sigs.Sid2Info[sid].n_time_channels, n_vals = n_elems  * rep.sigs.Sid2Info[sid].n_val_channels;
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


	//========================================================
	// Samples
	//========================================================

	// clear prediction points BEFORE a new set of predictions is done using the same instance
	void clear_samples() { samples.clear(); }

	// insert prediction points
	int insert_samples(int *pids, int *times, int n_samples) {
		for (int i=0; i<n_samples; i++)
			if (samples.insertRec(pids[i], times[i]) < 0)
				return -1;
		return 0;
	}

	int insert_sample(int pid, int time) { return insert_samples(&pid, &time, 1); }

	// normalize samples must be called after finishing inserting all samples.
	int normalize_samples() { samples.normalize(); return 0; }

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
		// run model to calculate predictions
		if (model.apply(rep, samples) < 0) {
			fprintf(stderr, "ERROR: MedAlgoMarkerInternal::get_preds FAILED.");
			return -1;
		}

		// export pids, times and preds to c arrays
		int j = 0;
		for (auto& idSample : samples.idSamples)
			for (auto& sample : idSample.samples) {
				_pids[j] = sample.id;
				times[j] = sample.time;
				preds[j++] = sample.prediction[0]; // This is Naive - but works for simple predictors giving the Raw score.
			}

		return 0;
	}

	int get_preds(MedSamples &_samples, float *preds) {

		samples = _samples;
		int n_samples = samples.nSamples();

		// run model to calculate predictions
		if (model.apply(rep, samples) < 0) {
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
	const char *get_name()	{ return name.c_str();	}
};
