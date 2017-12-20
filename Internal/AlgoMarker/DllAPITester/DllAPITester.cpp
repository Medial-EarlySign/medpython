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
			("direct_test" , "split to a dedicated debug routine")
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
	vector<long long> times(max_vals);
	vector<float> vals(max_vals);

	AM_API_ClearData(am);

	//MLOG("Going over %d pids\n", pids.size());
	for (auto pid : pids)
		for (auto &sig : sigs) {
			rep.uget(pid, sig, usv);
			int nelem = usv.len;
			if (nelem > 0) {
				long long *p_times = &times[0];
				float *p_vals = &vals[0];
				int i_time = 0;
				int i_val = 0;

				if (usv.n_time_channels() > 0) {
					for (int i=0; i<nelem; i++)
						for (int j=0; j<usv.n_time_channels(); j++)
							p_times[i_time++] = (long long)usv.Time(i, j);
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
	vector<long long> _timestamps;
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
	AM_API_CreateResponses(&resp);
	int calc_rc = AM_API_Calculate(am, req, resp);
	//MLOG("After Calculate: rc = %d\n", calc_rc);

	// go over reponses and pack them to a MesSample vector
	int n_resp = AM_API_GetResponsesNum(resp);
	MLOG("Got %d responses\n", n_resp);
	res.clear();
	int n_scr = 0;
	float _scr;
	int pid;
	long long ts;
	char *_scr_type;
	AMResponse *response;
	for (int i=0; i<n_resp; i++) {
		//MLOG("Getting response no. %d\n", i);
		int resp_rc = AM_API_GetResponseAtIndex(resp, i, &response);
		int n_scores;
		AM_API_GetResponseScoresNum(response, &n_scores);
		//int resp_rc = AM_API_GetResponse(resp, i, &pid, &ts, &n_scr, &_scr, &_scr_types);
		//MLOG("resp_rc = %d\n", resp_rc);
		//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);
		
		MedSample s;
		s.id = pid;
		s.time = (int)ts;
		if (resp_rc == AM_OK_RC && n_scores > 0) {
			resp_rc = AM_API_GetResponseScoreByIndex(response, 0, &pid, &ts, &_scr, &_scr_type);
			s.prediction.push_back(_scr);
		}
		else {
			s.prediction.push_back(AM_UNDEFINED_VALUE);
		}
		res.push_back(s);
	}



	// print error messages

	// AM level
	int n_msgs, *msg_codes;
	char **msgs_errs;
	AM_API_GetSharedMessages(resp, &n_msgs, &msg_codes, &msgs_errs);
	for (int i=0; i<n_msgs; i++) {
		MLOG("Shared Message %d : code %d : err: %s\n", n_msgs, msg_codes[i], msgs_errs[i]);
	}

	n_resp = AM_API_GetResponsesNum(resp);
	for (int i=0; i<n_resp; i++) {
		AMResponse *r;
		AM_API_GetResponseAtIndex(resp, i, &r);
		int n_scores;
		AM_API_GetResponseScoresNum(r, &n_scores);
		int p;
		long long t;
		float scr;
		char *scr_t;

		AM_API_GetResponseMessages(r, &n_msgs, &msg_codes, &msgs_errs);
		for (int k=0; k<n_msgs; k++) {
			MLOG("Response %d : Message %d : code %d : err: %s\n", i, k, msg_codes[k], msgs_errs[k]);
		}

		for (int j=0; j<n_scores; j++) {
			AM_API_GetScoreMessages(r, j, &n_msgs, &msg_codes, &msgs_errs);
			for (int k=0; k<n_msgs; k++) {
				MLOG("Response %d : score %d : Message %d : code %d : err: %s\n", i, j, k, msg_codes[k], msgs_errs[k]);
			}
		}
	}

	AM_API_DisposeRequest(req);
	AM_API_DisposeResponses(resp);

	MLOG("Finished getting preds from algomarker\n");
	return 0;
}



//=============================================================================================================================
int debug_me(po::variables_map &vm)
{

	// init AM
	AlgoMarker *test_am;

	if (AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating test algomarker\n");
		return -1;
	}

	MLOG("Name is %s\n", test_am->get_name());

	int load_rc;
	if ((load_rc = AM_API_Load(test_am, vm["amconfig"].as<string>().c_str())) != AM_OK_RC) {
		MERR("ERROR: Failed loading algomarker %s with config file %s ERR_CODE: %d\n", test_am->get_name(), vm["amconfig"].as<string>().c_str(), load_rc);
		return -1;
	}
	MLOG("Algomarker %s was loaded with config file %s\n", test_am->get_name(), test_am->get_config());


	// Load Data
	vector<long long> times ={ 20160101 };
	vector<float> vals ={ 6 };
	AM_API_AddData(test_am, 1, "RBC", (int)times.size(), &times[0], (int)vals.size(), &vals[0]);
	/*vector<float>*/ vals ={ 1960 };
	AM_API_AddData(test_am, 1, "BYEAR", 0, NULL, (int)vals.size(), &vals[0]);
	/*vector<float>*/ vals ={ 1 };
	//AM_API_AddData(test_am, 1, "GENDER", 0, NULL, (int)vals.size(), &vals[0]);

	// Calculate
	char *stypes[] ={ "Raw" };
	vector<int> _pids ={ 1 };
	vector<long long> _timestamps ={ 20160101 };
	AMRequest *req;
	MLOG("Creating Request\n");
	int req_create_rc = AM_API_CreateRequest("test_request", stypes, 1, &_pids[0], &_timestamps[0], (int)_pids.size(), &req);
	if (req == NULL)
		MLOG("ERROR: Got a NULL request !!\n");
	AMResponses *resp;

	// calculate scores
	MLOG("Before Calculate\n");
	AM_API_CreateResponses(&resp);
	int calc_rc = AM_API_Calculate(test_am, req, resp);
	

	// Shared messages
	int n_shared_msgs;
	int *shared_codes;
	char **shared_args;
	AM_API_GetSharedMessages(resp, &n_shared_msgs, &shared_codes, &shared_args);
	MLOG("Shared Messages: %d\n", n_shared_msgs);
	for (int i=0; i<n_shared_msgs; i++) {
		MLOG("Shared message %d : [%d] %s\n", i, shared_codes[i], shared_args[i]);
	}


	//AM_API_DisposeResponses(resp); resp=NULL;
	//AM_API_AddData(test_am, 1, "GENDER", 0, NULL, (int)vals.size(), &vals[0]);
	//AM_API_Calculate(test_am, req, &resp);

	// print result
	int n_resp = AM_API_GetResponsesNum(resp);
	MLOG("Got %d responses\n", n_resp);
	int n_scr = 0;
	float _scr;
	int pid;
	long long ts;
	char *_scr_type;
	AMResponse *response;
	for (int i=0; i<n_resp; i++) {
		MLOG("Getting response no. %d\n", i);

		AM_API_GetResponseAtIndex(resp, i, &response);
		int resp_rc = AM_API_GetResponseScoreByIndex(response, 0, &pid, &ts, &_scr, &_scr_type);
		MLOG("resp_rc = %d\n", resp_rc);
		MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, n_scr, _scr, _scr_type);
	}


	// print error messages

	// AM level
	int n_msgs, *msg_codes;
	char **msgs_errs;
	AM_API_GetSharedMessages(resp, &n_msgs, &msg_codes, &msgs_errs);
	for (int i=0; i<n_msgs; i++) {
		MLOG("Shared Message %d : code %d : err: %s\n", n_msgs, msg_codes[i], msgs_errs[i]);
	}


	// Dispose
	AM_API_DisposeRequest(req);
	AM_API_DisposeResponses(resp);
	AM_API_DisposeAlgoMarker(test_am);

	MLOG("Finished debug_me() test\n");

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


	if (vm.count("direct_test"))
		return debug_me(vm);

	// read model file
	MedModel model;
	if (model.read_from_file(vm["model"].as<string>()) < 0) {
		MERR("FAILED reading model file %s\n", vm["model"].as<string>().c_str());
		return -1;
	}

	unordered_set<string> sigs_set;
	vector<string> sigs;
	model.get_required_signal_names(sigs_set);

	MLOG("Required signals:");
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

	if (AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am) != AM_OK_RC) {
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
	int nbad = 0, n_miss = 0, n_similar = 0;
	if (res1.size() != res2.size()) {
		MLOG("ERROR:: Didn't get the same number of tests ... %d vs %d\n", res1.size(), res2.size());
	}

	MLOG("Comparing %d scores\n", res1.size());
	for (int i=0; i<res1.size(); i++) {

		if (res2[i].prediction[0] == (float)AM_UNDEFINED_VALUE) {
			n_miss++;
		} else if (res1[i].prediction[0] != res2[i].prediction[0]) {
			MLOG("ERROR !!!: #Res1 :: pid %d time %d pred %f #Res2 pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], res2[i].id, res2[i].time, res2[i].prediction[0]);
			nbad++;
		}
		else
			n_similar++;

	}


	MLOG(">>>>>TEST1: test DLL API batch: total %d : n_similar %d : n_bad %d : n_miss %d\n", res1.size(), n_similar, nbad, n_miss);
	if (nbad == 0) MLOG("PASSED\n"); else MLOG("FAILED\n");

}

//
// keep command line:
//
// typical test:
// Linux/Release/DllAPITester --model /nas1/Work/Users/Avi/Diabetes/order/pre2d/runs/partial/pre2d_partial_S6.model --samples test_100k.samples --amconfig /nas1/Work/Users/Avi/AlgoMarkers/pre2d/pre2d.amconfig
