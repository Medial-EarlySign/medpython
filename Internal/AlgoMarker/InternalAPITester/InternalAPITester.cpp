//
// Some basic tests for the in-mem API
// general plan:
// get a model file, a repository and a samples file
// get predictions using the regular API vs loading and predicting via the on mem API
//

#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedMat/MedMat/MedMat.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedProcessTools/MedProcessTools/SampleFilter.h>
#include <AlgoMarker/AlgoMarker/AlgoMarker.h>

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;
namespace po = boost::program_options;

//=========================================================================================================
int read_run_params(int argc, char *argv[], po::variables_map& vm) {
	po::options_description desc("Program options");

	try {
		desc.add_options()
			("help", "produce help message")
			("rep", po::value<string>()->default_value("/home/Repositories/THIN/thin_mar2017/thin.repository"), "repository file name")
			("samples", po::value<string>()->default_value(""), "medsamples file to use")
			("model", po::value<string>()->default_value(""), "model file to use")
			("test_filter", "flag to run a test on filters")
			;


		po::store(po::parse_command_line(argc, argv, desc), vm);
		if (vm.count("help")) {
			cerr << desc << "\n";
			exit(-1);

		}
		po::notify(vm);

		MLOG("=========================================================================================\n");
		MLOG("Command Line:");
		for (int i=0; i<argc; i++) MLOG(" %s", argv[i]);
		MLOG("\n");
		MLOG("..........................................................................................\n");
	}
	catch (exception& e) {
		cerr << "error: " << e.what() << "; run with --help for usage information\n";
		return -1;
	}
	catch (...) {
		cerr << "Exception of unknown type!\n";
		return -1;
	}

	return 0;
}

//=================================================================================================================
int get_preds_in_mem_loading(string rep_conf, MedPidRepository &rep, MedModel &model, MedSamples &samples, vector<int> &pids, vector<string> &sigs, vector<MedSample> &res)
{
	MedPidRepository rep_in_mem;

	if (rep_in_mem.init(rep_conf) < 0) return -1;

	MLOG("After initializing rep_in_mem\n");
	UniversalSigVec usv;

	int max_vals = 100000;
	vector<int> times(max_vals);
	vector<float> vals(max_vals);

	rep_in_mem.switch_to_in_mem_mode();
	for (auto pid : pids)
		for (auto &sig : sigs) {
			rep.uget(pid, sig, usv);
			int nelem = usv.len;
			if (nelem > 0) {
				int *p_times = &times[0];
				float *p_vals = &vals[0];
				int i_time = 0;
				int i_val = 0;

				if (usv.n_time_channels() > 0) {
					for (int i=0; i<nelem; i++)
						for (int j=0; j<usv.n_time_channels(); j++)
							p_times[i_time++] = usv.Time(i, j);
				}
				else
					p_times = NULL;

				if (usv.n_val_channels() > 0) {
					for (int i=0; i<nelem; i++)
						for (int j=0; j<usv.n_val_channels(); j++)
							p_vals[i_val++] = usv.Val(i, j);
				}
				else
					p_vals = NULL;

				//MLOG("Inserting pid %d sig %s nelem %d i_time %d i_val %d\n", pid, sig.c_str(), nelem, i_time, i_val);
				rep_in_mem.in_mem_rep.insertData(pid, sig.c_str(), p_times, p_vals, i_time, i_val);
			}
		}

	// finish rep loading 
	MLOG("Before Sorting in_mem_rep\n");
	rep_in_mem.in_mem_rep.sortData();
	MLOG("After Sorting in_mem_rep\n");

	// actual prediction
	model.apply(rep_in_mem, samples);

	// get res out
	res.clear();
	samples.export_to_sample_vec(res);

	return 0;
}

//=================================================================================================================
int get_preds_with_medalgomarker_internal(string rep_conf, MedPidRepository &rep, MedSamples &samples, string model_fname, vector<int> &pids, vector<string> &sigs, vector<float> &res)
{
	MedAlgoMarkerInternal marker;

	marker.set_name("Testing");

	// init rep config
	MLOG("init marker with rep_conf %s\n", rep_conf.c_str());
	if (marker.init_rep_config(rep_conf.c_str()) < 0) return -1;

	MLOG("init model from file %s\n", model_fname.c_str());
	// init model
	if (marker.init_model_from_file(model_fname.c_str()) < 0) return -2;

	// init samples
	vector<int> sample_pids, sample_times;
	for (auto &id : samples.idSamples)
		for (auto &s : id.samples) {
			sample_pids.push_back(s.id);
			sample_times.push_back(s.time);
		}
	//MLOG("Initializeing with %d samples\n", sample_pids.size());
	//if (marker.init_samples(&sample_pids[0], &sample_times[0], (int)sample_pids.size()) < 0) return -3;

	// init data
	marker.data_load_init();
	int max_vals = 100000;
	vector<int> times(max_vals);
	vector<float> vals(max_vals);

	UniversalSigVec usv;
	for (auto pid : pids)
		for (auto &sig : sigs) {
			rep.uget(pid, sig, usv);
			int nelem = usv.len;
			if (nelem > 0) {
				int *p_times = &times[0];
				float *p_vals = &vals[0];
				int i_time = 0;
				int i_val = 0;

				if (usv.n_time_channels() > 0) {
					for (int i=0; i<nelem; i++)
						for (int j=0; j<usv.n_time_channels(); j++)
							p_times[i_time++] = usv.Time(i, j);
				}
				else
					p_times = NULL;

				if (usv.n_val_channels() > 0) {
					for (int i=0; i<nelem; i++)
						for (int j=0; j<usv.n_val_channels(); j++)
							p_vals[i_val++] = usv.Val(i, j);
				}
				else
					p_vals = NULL;

				//MLOG("Inserting pid %d sig %s nelem %d i_time %d i_val %d\n", pid, sig.c_str(), nelem, i_time, i_val);
				if (marker.data_load_pid_sig(pid, sig.c_str(), p_times, p_vals, nelem) < 0) return -4;
			}
		}

	if (marker.data_load_end() < 0) return -5;


	// actual prediction
	res.resize(sample_pids.size(), 0);
	marker.get_preds(&sample_pids[0], &sample_times[0], &res[0], (int)sample_pids.size());
	//marker.get_preds(samples, &res[0]);

	return 0;
}


//=================================================================================================================
int get_preds_with_medalgomarker_internal_loop(string rep_conf, MedPidRepository &rep, MedSamples &samples, string model_fname, vector<int> &pids, vector<string> &sigs, vector<float> &res)
{
	MedAlgoMarkerInternal marker;

	marker.set_name("Testing");

	// init rep config
	MLOG("init marker with rep_conf %s\n", rep_conf.c_str());
	if (marker.init_rep_config(rep_conf.c_str()) < 0) return -1;

	MLOG("init model from file %s\n", model_fname.c_str());
	// init model
	if (marker.init_model_from_file(model_fname.c_str()) < 0) return -2;

	// init samples
	vector<int> sample_pids, sample_times;
	for (auto &id : samples.idSamples)
		for (auto &s : id.samples) {
			sample_pids.push_back(s.id);
			sample_times.push_back(s.time);
		}

	// init data
	marker.data_load_init();
	int max_vals = 100000;
	vector<int> times(max_vals);
	vector<float> vals(max_vals);

	UniversalSigVec usv;
	res.resize(sample_pids.size(), 0);
	MedTimer timer;
	timer.start();
	for (int k=0; k<sample_pids.size(); k++) {
		if (k%10000 == 0) { 
			timer.take_curr_time();
			double per_sec = (double)k/timer.diff_sec();
			MLOG(">>>>>>> k = %d     %f calls/sec <<<<<<<<<\n", k, per_sec); 
		}
		int pid = sample_pids[k];
		if (k==0 || pid != sample_pids[k-1]) {
			marker.clear_data();
			for (auto &sig : sigs) {
				rep.uget(pid, sig, usv);
				int nelem = usv.len;
				if (nelem > 0) {
					int *p_times = &times[0];
					float *p_vals = &vals[0];
					int i_time = 0;
					int i_val = 0;

					if (usv.n_time_channels() > 0) {
						for (int i=0; i<nelem; i++)
							for (int j=0; j<usv.n_time_channels(); j++)
								p_times[i_time++] = usv.Time(i, j);
					}
					else
						p_times = NULL;

					if (usv.n_val_channels() > 0) {
						for (int i=0; i<nelem; i++)
							for (int j=0; j<usv.n_val_channels(); j++)
								p_vals[i_val++] = usv.Val(i, j);
					}
					else
						p_vals = NULL;

					//MLOG("Inserting pid %d sig %s nelem %d i_time %d i_val %d\n", pid, sig.c_str(), nelem, i_time, i_val);
					if (marker.data_load_pid_sig(pid, sig.c_str(), p_times, p_vals, nelem) < 0) return -4;
				}
			}

			if (marker.data_load_end() < 0) return -5;
		}
		int time = sample_times[k];
		marker.get_pred(&pid, &time, &res[k]);

	}


	return 0;
}

//========================================================================================================
int read_files(po::variables_map &vm, MedModel &model, vector<string> &sigs, MedSamples &samples, vector<int> &pids, MedPidRepository &rep)
{
	// read model file
	if (model.read_from_file(vm["model"].as<string>()) < 0) {
		MERR("FAILED reading model file %s\n", vm["model"].as<string>().c_str());
		return -1;
	}

	unordered_set<string> sigs_set;
	//vector<string> sigs;
	model.get_required_signal_names(sigs_set);

	MLOG("Reuired signals:");
	for (auto &sig : sigs_set) {
		MLOG(" %s", sig.c_str());
		sigs.push_back(sig);
	}
	MLOG("\n");

	// read samples file
	if (samples.read_from_file(vm["samples"].as<string>())) {
		MERR("FAILES reading samples file %s\n", vm["samples"].as<string>().c_str());
		return -1;
	}

	samples.get_ids(pids);

	MLOG("Read samples file %s with %d samples from %d pids\n", vm["samples"].as<string>().c_str(), samples.nSamples(), pids.size());


	// read rep

	if (rep.read_all(vm["rep"].as<string>(), pids, sigs) < 0) return -1;

	return 0;
}

//=================================================================================================
int test_filter(po::variables_map &vm)
{
	
	MedModel model;
	vector<string> sigs;
	MedSamples samples;
	vector<int> pids;
	MedPidRepository rep;

	if (read_files(vm, model, sigs, samples, pids, rep) < 0) return -1;

	vector<SanitySimpleFilter> bfilters;

	SanitySimpleFilter bf;

	bf.init_from_string("sig=Glucose,win_from=0,win_to=730,min_NVals=1,max_Nvals=3");	bfilters.push_back(bf);
	bf.init_from_string("sig=Glucose,win_from=0,win_to=365,min_val=50,max_val=150,max_outliers=0,time_unit=Days");	bfilters.push_back(bf);

	MLOG("initialized %d filters\n", bfilters.size());

	vector<MedSample> svec;
	samples.export_to_sample_vec(svec);

	for (auto &s : svec) {
		for (int j=0; j<bfilters.size(); j++) {
			int nvals, noutliers, rc;
			if ((rc = bfilters[j].test_filter(s, rep, nvals, noutliers)) != SanitySimpleFilter::Passed) {
				MLOG("failed filter %d on sample %d,%d with code %d nvals %d noutliers %d\n", j, s.id, s.time, rc, nvals, noutliers);
			}
		}
	}


	return 0;

}


//========================================================================================
// MAIN
//========================================================================================

int main(int argc, char *argv[])
{
	int rc = 0;
	po::variables_map vm;

	// Running Parameters
	MLOG("Reading params\n");
	rc = read_run_params(argc, argv, vm);
	assert(rc >= 0);


	if (vm.count("test_filter")) return test_filter(vm);

	MedModel model;
	vector<string> sigs;
	MedSamples samples;
	vector<int> pids;
	MedPidRepository rep;

	if (read_files(vm, model, sigs, samples, pids, rep) < 0) return -1;

	MedSamples samples2 = samples;

/*
	// read model file
	MedModel model;
	if (model.read_from_file(vm["model"].as<string>()) < 0) {
		MERR("FAILED reading model file %s\n", vm["model"].as<string>().c_str());
		return -1;
	}

	unordered_set<string> sigs_set;
	vector<string> sigs;
	model.get_required_signal_names(sigs_set);

	MLOG("Reuired signals:");
	for (auto &sig : sigs_set) {
		MLOG(" %s", sig.c_str());
		sigs.push_back(sig);
	}
	MLOG("\n");

	// read samples file
	MedSamples samples, samples2;

	if (samples.read_from_file(vm["samples"].as<string>())) {
		MERR("FAILES reading samples file %s\n", vm["samples"].as<string>().c_str());
		return -1;
	}
	samples2 = samples;

	vector<int> pids;
	samples.get_ids(pids);

	MLOG("Read samples file %s with %d samples from %d pids\n", vm["samples"].as<string>().c_str(), samples.nSamples(), pids.size());


	// read rep
	MedPidRepository rep;

	if (rep.read_all(vm["rep"].as<string>(), pids, sigs) < 0) return -1;

*/
	// apply model (+ print top 50 scores)
	model.apply(rep, samples);

	// printing
	vector<MedSample> res1;
	samples.export_to_sample_vec(res1);
	for (int i=0; i<min(50, (int)res1.size()); i++) {
		MLOG("#Res1 :: pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0]);
	}


	//===============================================================================
	// TEST1: testing internal in_mem in a repository
	//===============================================================================
	// get preds via direct loading of signals
	vector<MedSample> res2;
	get_preds_in_mem_loading(vm["rep"].as<string>(), rep, model, samples2, pids, sigs, res2);
	for (int i=0; i<min(50, (int)res1.size()); i++) {
		MLOG("#Res1 :: pid %d time %d pred %f #Res2 pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], res2[i].id, res2[i].time, res2[i].prediction[0]);
	}

	// test results
	int nbad = 0;
	for (int i=0; i<res1.size(); i++) {
		if (res1[i].prediction[0] != res2[i].prediction[0]) {
			MLOG("ERROR !!!: #Res1 :: pid %d time %d pred %f #Res2 pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], res2[i].id, res2[i].time, res2[i].prediction[0]);
			nbad++;
		}
	}


	MLOG(">>>>>TEST1: test internal memory rep in repository: ");
	if (nbad == 0) MLOG("PASSED\n"); else MLOG("FAILED\n");



	//===============================================================================
	// Testing with internal MedAlgo API
	//===============================================================================

	MLOG("Direct test of Internal MedAlgoMarker API\n");
	vector<float> preds;
	int rc_test;
	if ((rc_test = get_preds_with_medalgomarker_internal(vm["rep"].as<string>(), rep, samples, vm["model"].as<string>(), pids, sigs, preds)) < 0) {
		MLOG("FAILED with rc %d\n", rc_test);
		return -1;
	}
	// test results
	nbad = 0;
	for (int i=0; i<res1.size(); i++) {
		if (res1[i].prediction[0] != preds[i]) {
			MLOG("ERROR !!!: #Res1 :: pid %d time %d pred %f !=  pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], preds[i]);
			nbad++;
		}
	}


	MLOG(">>>>>TEST2: test internal MedAlgoMarker class: ");
	if (nbad == 0) MLOG("PASSED\n"); else MLOG("FAILED\n");


	//===============================================================================
	// Testing with internal MedAlgo API loop
	//===============================================================================

	MLOG("Direct test of Internal MedAlgoMarker API Loop\n");
	preds.clear();
	if ((rc_test = get_preds_with_medalgomarker_internal_loop(vm["rep"].as<string>(), rep, samples, vm["model"].as<string>(), pids, sigs, preds)) < 0) {
		MLOG("FAILED with rc %d\n", rc_test);
		return -1;
	}
	// test results
	nbad = 0;
	for (int i=0; i<res1.size(); i++) {
		if (res1[i].prediction[0] != preds[i]) {
			MLOG("ERROR !!!: #Res1 :: pid %d time %d pred %f !=  pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], preds[i]);
			nbad++;
		}
	}


	MLOG(">>>>>TEST3: test internal MedAlgoMarker class loop: ");
	if (nbad == 0) MLOG("PASSED\n"); else MLOG("FAILED\n");
}



