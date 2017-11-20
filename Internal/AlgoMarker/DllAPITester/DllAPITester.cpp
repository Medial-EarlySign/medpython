//
// Test Program to the Dll API
//
// General Plan :
// 
// Compare same data/model/points prediction using the infrastructure directly AND using the DLL.
//

#define AM_DLL_IMPORT

#include "AlgoMarker/AlgoMarker/AlgoMarker.h"

#include <string>
#include <iostream>
#include <boost/program_options.hpp>


#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
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
			("amconfig" , po::value<string>()->default_value(""), "algo marker configuration file")
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
int get_preds_from_algomarker(AlgoMarker *am, string rep_conf, MedPidRepository &rep, MedModel &model, MedSamples &samples, vector<int> &pids, vector<string> &sigs, vector<MedSample> &res)
{
	UniversalSigVec usv;

	int max_vals = 100000;
	vector<long> times(max_vals);
	vector<float> vals(max_vals);

	AM_API_ClearData(am);

	//MLOG("Going over %d pids\n", pids.size());
	for (auto pid : pids)
		for (auto &sig : sigs) {
			rep.uget(pid, sig, usv);
			int nelem = usv.len;
			if (nelem > 0) {
				long *p_times = &times[0];
				float *p_vals = &vals[0];
				int i_time = 0;
				int i_val = 0;

				if (usv.n_time_channels() > 0) {
					for (int i=0; i<nelem; i++)
						for (int j=0; j<usv.n_time_channels(); j++)
							p_times[i_time++] = (long)usv.Time(i, j);
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

				//MLOG("Adding data: pid %d sig %s n_times %d n_vals %d\n", pid, sig.c_str(), i_time, i_val);
				AM_API_AddData(am, pid, sig.c_str(), i_time, p_times, i_val, p_vals);
			}
		}


	// finish rep loading 
	char *stypes[] ={ "Raw" };
	vector<int> _pids;
	vector<long> _timestamps;
	vector<MedSample> _vsamp;
	samples.export_to_sample_vec(_vsamp);
	for (auto &s : _vsamp) {
		_pids.push_back(s.id);
		_timestamps.push_back(s.time);
	}

	//MLOG("Before CreateRequest\n");
	// prep request
	AMRequest *req;
	int req_create_rc = AM_API_CreateRequest("test_request", stypes, 1, &_pids[0], &_timestamps[0], (int)_pids.size(), &req);
	if (req == NULL)
		MLOG("ERROR: Got a NULL request !!\n");
	AMResponses *resp;

	// calculate scores
	//MLOG("Before Calculate\n");
	int calc_rc = AM_API_Calculate(am, req, &resp);
	//MLOG("After Calculate: rc = %d\n", calc_rc);

	// go over reponses and pack them to a MesSample vector
	int n_resp = AM_API_GetResponsesNum(resp);
	MLOG("Got %d responses\n", n_resp);
	res.clear();
	int n_scr = 0;
	float *_scr = NULL;
	int pid;
	long ts;
	char **_scr_types;
	for (int i=0; i<n_resp; i++) {
		//MLOG("Getting response no. %d\n", i);
		int resp_rc = AM_API_GetResponse(resp, i, &pid, &ts, &n_scr, &_scr, &_scr_types);
		//MLOG("resp_rc = %d\n", resp_rc);
		//MLOG("i %d , pid %d ts %d n_scr %d scr %f\n", i, pid, ts, n_scr, _scr[0]);
		MedSample s;
		s.id = pid;
		s.time = (int)ts;
		if (n_scr > 0)
			s.prediction.push_back(_scr[0]);
		res.push_back(s);
	}


	AM_API_DisposeRequest(req);
	AM_API_DisposeResponses(resp);

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


	// apply model (+ print top 50 scores)
	model.apply(rep, samples);

	// printing
	vector<MedSample> res1;
	samples.export_to_sample_vec(res1);
	for (int i=0; i<min(50, (int)res1.size()); i++) {
		MLOG("#Res1 :: pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0]);
	}

	// res1 now is the gold standard we compare to : was calculated directly using the infrastructure

	//===============================================================================
	// TEST1: testing internal in_mem in a repository
	//===============================================================================
	AlgoMarker *test_am;

	if (AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, "Pre2D", &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating test algomarker\n");
		return -1;
	}

	MLOG("Name is %s\n", test_am->get_name());

	if ((rc = AM_API_Load(test_am, vm["amconfig"].as<string>().c_str())) != AM_OK_RC) {
		MERR("ERROR: Failed loading algomarker %s with config file %s ERR_CODE: %d\n", test_am->get_name(), vm["amconfig"].as<string>().c_str(), rc);
		return -1;
	}

	MLOG("Algomarker %s was loaded with config file %s\n", test_am->get_name(), test_am->get_config());
	vector<MedSample> res2;
	get_preds_from_algomarker(test_am, vm["rep"].as<string>(), rep, model, samples2, pids, sigs, res2);
	for (int i=0; i<min(50, (int)res1.size()); i++) {
		MLOG("#Res1 :: pid %d time %d pred %f #Res2 pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], res2[i].id, res2[i].time, res2[i].prediction[0]);
	}

	// test results
	int nbad = 0;
	if (res1.size() != res2.size()) {
		MLOG("ERROR:: Didn't get the same number of tests ... %d vs %d\n", res1.size(), res2.size());
	}

	MLOG("Comparing %d scores\n", res1.size());
	for (int i=0; i<res1.size(); i++) {
		if (res1[i].prediction[0] != res2[i].prediction[0]) {
			MLOG("ERROR !!!: #Res1 :: pid %d time %d pred %f #Res2 pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], res2[i].id, res2[i].time, res2[i].prediction[0]);
			nbad++;
		}
	}

	MLOG(">>>>>TEST1: test DLL API batch: ");
	if (nbad == 0) MLOG("PASSED\n"); else MLOG("FAILED\n");

}

