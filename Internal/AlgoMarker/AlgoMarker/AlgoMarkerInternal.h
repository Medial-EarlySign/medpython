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
	void set_name(char *_name) { name = string(_name); }

	// init repository config
	int init_rep_config(char *config_fname) { if (rep.read_config(string(config_fname)) < 0) return -1; rep.switch_to_in_mem_mode(); }

	// init pids
	void set_pids(int *_pids, int npids) { pids.clear(); pids.assign(_pids, _pids + npids); }

	// init rep , model , samples
	int init_rep_with_file_data(char *_rep_fname) {
		rep.clear();
		rep_fname = string(_rep_fname);
		vector<string> sigs ={};
		return (rep.read_all(rep_fname, pids, sigs));
	}

	// init model
	int init_model_from_file(char *_model_fname) {	model.clear();	return (model.read_from_file(string(_model_fname))); }


	// init samples
	int init_samples(int *pids, int *times, int n_samples) { clear_samples(); return insert_samples(pids, times, n_samples); }
	int init_samples(int pid, int time) { return init_samples(&pid, &time, 1); }  // single prediction point initiation 


   //========================================================
   // Loading data to rep
   //========================================================
	
    // init loading : actions that must be taken BEFORE any loading starts
	int data_load_init() { rep.switch_to_in_mem_mode(); }

	// load n_elems for a pid,sig
	int data_load_pid_sig(int pid, char *sig_name, int *times, float *vals, int n_elems) {
		int n_times = n_elems, n_vals = n_elems;
		if (times == NULL) n_times = 0;
		if (vals == NULL) n_vals = 0;
		return rep.in_mem_rep.insertData(pid, sig_name, times, vals, n_times, n_vals);
	}

	// load a single element for a pid,sig
	int data_load_pid_sig(int pid, char *sig_name, int *times, float *vals) { return data_load_pid_sig(pid, sig_name, times, vals, 1); }

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

	//========================================================
	// Calculate predictions
	//========================================================
	int get_preds(int *pids, int *times, float *preds, int n_samples) {

		// init_samples
		init_samples(pids, times, n_samples);

		// run model to calculate predictions
		if (model.apply(rep, samples) < 0) {
			fprintf(stderr, "ERROR: MedAlgoMarkerInternal::get_preds FAILED.");
			return -1;
		}

		// export pids, times and preds to c arrays
		int j = 0;
		for (auto& idSample : samples.idSamples)
			for (auto& sample : idSample.samples) {
				pids[j] = sample.id;
				times[j] = sample.time;
				preds[j++] = sample.prediction[0]; // This Naive - but works for simple predictors giving the Raw score.
			}
	}

	int get_pred(int *pid, int *time, float *pred) { return get_preds(pid, time, pred, 1); }


	//========================================================
	// Clearing - freeing mem
	//========================================================
	void clear() { pids.clear(); model.clear(); samples.clear(); rep.in_mem_rep.clear(); rep.clear(); }


	//========================================================
	// a few more needed APIs
	//========================================================
	const char *get_name()	{ return name.c_str();	}
};
