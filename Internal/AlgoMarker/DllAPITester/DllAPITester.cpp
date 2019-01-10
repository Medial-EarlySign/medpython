//
// Test Program to the Dll API
//
// General Plan :
// 
// Compare same data/model/points prediction using the infrastructure directly AND using the DLL.
//

#define AM_DLL_IMPORT

#include <AlgoMarker/AlgoMarker/AlgoMarker.h>

#include <string>
#include <iostream>
#include <boost/program_options.hpp>


#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedIO/MedIO/MedIO.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

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
			("direct_test", "split to a dedicated debug routine")
			("single", "run test in single mode, instead of the default batch")
			("print_msgs", "print algomarker messages when testing batches or single (direct test always prints them)")
			("test_data", po::value<string>()->default_value(""), "test data for --direct_test option")
			("json_data", po::value<string>()->default_value(""), "test json data for --direct_test option")
			("date", po::value<long long>()->default_value(20180101), "test date")
			("egfr_test", "split to a debug routine for the simple egfr algomarker")
			("data_api_test", "split to a test of data api")
			("pid", po::value<int>()->default_value(5000100), "test data_api for this pid, use --data_api_test option")
			("sig", po::value<string>()->default_value("Creatinine"), "test data_api for this signal, use --data_api_test option")

			("kp_test", "split to a dedicated test for kp data")
			("kp_demographic", po::value<string>()->default_value(""), "demographic for --kp_test option: each line: <pid> <byear> <F/M>")
			("kp_data", po::value<string>()->default_value(""), "lab tests for --kp_test option: each line: <pid> <code> <date> <value>")
			("kp_scores", po::value<string>()->default_value(""), "which scores to generate for  --kp_test option: each line start: <pid> <date> ...")
			("kp_codes", po::value<string>()->default_value(""), "lab tests codes for --kp_test option: each line: <code> <name>")
			("kp_fout", po::value<string>()->default_value(""), "output file  for --kp_test option")

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
int get_preds_from_algomarker(AlgoMarker *am, string rep_conf, MedPidRepository &rep, MedModel &model, MedSamples &samples, vector<int> &pids, vector<string> &sigs, vector<MedSample> &res, int print_msgs)
{
	UniversalSigVec usv;

	int max_vals = 100000;
	vector<long long> times(max_vals);
	vector<float> vals(max_vals);

	AM_API_ClearData(am);

	MLOG("=====> now running get_preds_from_algomarker()\n");
	MLOG("Going over %d pids\n", pids.size());
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

	((MedialInfraAlgoMarker *)am)->set_sort(0); // getting rid of cases in which multiple data sets on the same day cause differences and fake failed tests.

	MLOG("After AddData for all batch\n");
	// finish rep loading 
	char *stypes[] ={ "Raw" };
	vector<int> _pids;
	vector<long long> _timestamps;
	vector<MedSample> _vsamp;
	samples.export_to_sample_vec(_vsamp);
	for (auto &s : _vsamp) {
		_pids.push_back(s.id);
		//_timestamps.push_back((long long)s.time*10000 + 1010);
		_timestamps.push_back((long long)s.time);
		//MLOG("pid %d time %lld\n", _pids.back(), _timestamps.back());
	}



	//MLOG("Before CreateRequest\n");
	// prep request
	AMRequest *req;
	int req_create_rc = AM_API_CreateRequest("test_request", stypes, 1, &_pids[0], &_timestamps[0], (int)_pids.size(), &req);
	if (req == NULL)
		MLOG("ERROR: Got a NULL request rc = %d!!\n", req_create_rc);
	AMResponses *resp;

	// calculate scores
	MLOG("Before Calculate\n");
	AM_API_CreateResponses(&resp);
	int calc_rc = AM_API_Calculate(am, req, resp);
	MLOG("After Calculate: rc = %d\n", calc_rc);


	// go over reponses and pack them to a MesSample vector
	int n_resp = AM_API_GetResponsesNum(resp);
	MLOG("Got %d responses\n", n_resp);
	res.clear();
	int pid;
	long long ts;
	char *_scr_type = NULL;
	AMResponse *response;
	for (int i=0; i<n_resp; i++) {
		//MLOG("Getting response no. %d\n", i);
		int resp_rc = AM_API_GetResponseAtIndex(resp, i, &response);
		int n_scores;
		AM_API_GetResponseScoresNum(response, &n_scores);
		//int resp_rc = AM_API_GetResponse(resp, i, &pid, &ts, &n_scr, &_scr, &_scr_type);
		//MLOG("resp_rc = %d\n", resp_rc);
		//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);

		AM_API_GetResponsePoint(response, &pid, &ts);
		MedSample s;
		s.id = pid;
		if (ts > 30000000)
			s.time = (int)(ts/10000);
		else
			s.time = ts;
		if (resp_rc == AM_OK_RC && n_scores > 0) {
			float _scr;
			resp_rc = AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
			//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);
			s.prediction.push_back(_scr);
		}
		else {
			s.prediction.push_back((float)AM_UNDEFINED_VALUE);
		}
		res.push_back(s);
	}


	if (print_msgs) {

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
	}

	AM_API_DisposeRequest(req);
	AM_API_DisposeResponses(resp);

	MLOG("Finished getting preds from algomarker\n");
	return 0;
}


#if 1
//=================================================================================================================
// same test, but running each point in a single mode, rather than batch on whole.
//=================================================================================================================
int get_preds_from_algomarker_single(AlgoMarker *am, string rep_conf, MedPidRepository &rep, MedModel &model, MedSamples &samples, vector<int> &pids, vector<string> &sigs, vector<MedSample> &res, vector<MedSample> &compare_res, int print_msgs)
{
	UniversalSigVec usv;

	int max_vals = 100000;
	vector<long long> times(max_vals);
	vector<float> vals(max_vals);

	AM_API_ClearData(am);

	MLOG("=====> now running get_preds_from_algomarker_single()\n");
	MLOG("Going over %d samples\n", samples.nSamples());
	int n_tested = 0;

	MedTimer timer;
	timer.start();
	for (auto &id : samples.idSamples)
		for (auto &s : id.samples) {

			// adding all data 
			for (auto &sig : sigs) {
				rep.uget(s.id, sig, usv);
				int nelem = usv.len;
				if (nelem > 0) {
					long long *p_times = &times[0];
					float *p_vals = &vals[0];
					int i_time = 0;
					int i_val = 0;

					int nelem_before = 0;

					if (usv.n_time_channels() > 0) {
						for (int i=0; i<nelem; i++)
							for (int j=0; j<usv.n_time_channels(); j++) {
								if (usv.Time(i, j) <= s.time) nelem_before = i+1;
								else break;
								p_times[i_time++] = (long long)usv.Time(i, j);
							}
					}
					else
						p_times = NULL;

					if (usv.n_val_channels() > 0) {
						if (p_times != NULL) nelem = nelem_before;
						for (int i=0; i<nelem; i++)
							for (int j=0; j<usv.n_val_channels(); j++)
								p_vals[i_val++] = usv.Val(i, j);
					}
					else
						p_vals = NULL;

					//MLOG("Adding data: pid %d time %d sig %s n_times %d n_vals %d\n", s.id, s.time, sig.c_str(), i_time, i_val);
					if ((i_val > 0) || (i_time > 0))
						AM_API_AddData(am, s.id, sig.c_str(), i_time, p_times, i_val, p_vals);
				}
			}

			// At this point we can send to the algomarker and ask for a score

			// a small technicality
			((MedialInfraAlgoMarker *)am)->set_sort(0); // getting rid of cases in which multiple data sets on the same day cause differences and fake failed tests.

			// preparing a request
			char *stypes[] ={ "Raw" };
			long long _timestamp = (long long)s.time;

			AMRequest *req;
			int req_create_rc = AM_API_CreateRequest("test_request", stypes, 1, &s.id, &_timestamp, 1, &req);
			if (req == NULL) {
				MLOG("ERROR: Got a NULL request for pid %d time %d rc %d!!\n", s.id, s.time, req_create_rc);
				return -1;
			}

			// create a response
			AMResponses *resp;
			AM_API_CreateResponses(&resp);

			// Calculate
			AM_API_Calculate(am, req, resp);
			//int calc_rc = AM_API_Calculate(am, req, resp);
			//MLOG("after Calculate: calc_rc %d\n", calc_rc);

			int n_resp = AM_API_GetResponsesNum(resp);

			//MLOG("pid %d time %d n_resp %d\n", s.id, s.time, n_resp);
			// get scores
			if (n_resp == 1) {
				AMResponse *response;
				int resp_rc = AM_API_GetResponseAtIndex(resp, 0, &response);
				int n_scores;
				AM_API_GetResponseScoresNum(response, &n_scores);
				if (n_scores == 1) {
					float _scr;
					int pid;
					long long ts;
					char *_scr_type = NULL;
					AM_API_GetResponsePoint(response, &pid, &ts);

					MedSample rs;
					rs.id = pid;
					if (ts > 30000000)
						rs.time = (int)(ts/10000);
					else
						rs.time = ts;

					if (resp_rc == AM_OK_RC && n_scores > 0) {
						resp_rc = AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
						//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);
						rs.prediction.push_back(_scr);
					}
					else {
						rs.prediction.push_back((float)AM_UNDEFINED_VALUE);
					}
					res.push_back(rs);

					//MLOG("pid %d ts %d scr %f %s\n", pid, ts, _scr, _scr_type);
				} 

				//int resp_rc = AM_API_GetResponse(resp, i, &pid, &ts, &n_scr, &_scr, &_scr_type);
				//MLOG("resp_rc = %d\n", resp_rc);

			}
			else {
				MedSample rs = s;
				rs.prediction.clear();
				rs.prediction.push_back((float)AM_UNDEFINED_VALUE);
				res.push_back(rs);
			}

			if (print_msgs) {
				// print error messages
				// AM level
				int n_msgs, *msg_codes;
				char **msgs_errs;
				AM_API_GetSharedMessages(resp, &n_msgs, &msg_codes, &msgs_errs);
				for (int i=0; i<n_msgs; i++) {
					MLOG("pid %d time %d Shared Message %d : code %d : err: %s\n", s.id, s.time, n_msgs, msg_codes[i], msgs_errs[i]);
				}

				n_resp = AM_API_GetResponsesNum(resp);
				for (int i=0; i<n_resp; i++) {
					AMResponse *r;
					AM_API_GetResponseAtIndex(resp, i, &r);
					int n_scores;
					AM_API_GetResponseScoresNum(r, &n_scores);

					AM_API_GetResponseMessages(r, &n_msgs, &msg_codes, &msgs_errs);
					for (int k=0; k<n_msgs; k++) {
						MLOG("pid %d time %d Response %d : Message %d : code %d : err: %s\n", s.id, s.time, i, k, msg_codes[k], msgs_errs[k]);
					}

					for (int j=0; j<n_scores; j++) {
						AM_API_GetScoreMessages(r, j, &n_msgs, &msg_codes, &msgs_errs);
						for (int k=0; k<n_msgs; k++) {
							MLOG("pid %d time %d Response %d : score %d : Message %d : code %d : err: %s\n", s.id, s.time, i, j, k, msg_codes[k], msgs_errs[k]);
						}
					}
				}
			}
			// and now need to dispose responses and request
			AM_API_DisposeRequest(req);
			AM_API_DisposeResponses(resp);

			// clearing data in algomarker
			AM_API_ClearData(am);

			float pred = res.back().prediction[0];
			float compare_pred = compare_res[n_tested].prediction[0];

			if ((pred != (float)AM_UNDEFINED_VALUE) && (pred != compare_pred)) {
				MLOG("ERROR Found: pid %d time %d : pred %f compared to %f ...\n", s.id, s.time, pred, compare_pred);
			}

			n_tested++;
			if ((n_tested % 100) == 0) {
				timer.take_curr_time();
				double dt = timer.diff_sec();
				MLOG("Tested %d samples : time %f sec\n", n_tested, dt);
				dt = (double)n_tested/dt;
				MLOG("%f samples/sec\n", dt);
			}
		}





	MLOG("Finished getting preds from algomarker in a single manner\n");
	return 0;
}


#endif

//=============================================================================================================================
int load_algomarker_from_string(AlgoMarker *am, int pid, const string &sdata)
{

	// example format:
	// Glucose:80.5,123:20100103,20120809;BYEAR:1985;GENDER:1

	vector<string> sigs;
	boost::split(sigs, sdata, boost::is_any_of(";"));

	for (auto &s : sigs) {
		if (s == "") continue;

		vector<string> f;
		boost::split(f, s, boost::is_any_of(":"));
		string sig = f[0];
		// values
		vector<float> vals;
		if (f.size() >= 2) {
			vector<string> svals;
			boost::split(svals, f[1], boost::is_any_of(","));
			for (auto v : svals) vals.push_back(stof(v));
		}

		vector<long long> times;
		if (f.size() >= 3) {
			vector<string> tvals;
			boost::split(tvals, f[2], boost::is_any_of(","));
			for (auto v : tvals) times.push_back(stoll(v));
		}

		MLOG("Adding Data: sig %s :: vals: ", sig.c_str());
		for (auto v : vals) MLOG("%f, ", v);
		MLOG(" times: ");
		for (auto v : times) MLOG("%lld, ", v);
		MLOG("\n");

		long long *pt = NULL;
		float *pv = NULL;
		if (times.size() > 0) pt = &times[0];
		if (vals.size() > 0) pv = &vals[0];
		AM_API_AddData(am, pid, sig.c_str() , (int)times.size(), pt, (int)vals.size(), pv);

	}

	return 0;
}

//=============================================================================================================================
int input_json_to_string(string json_fname, string &sdata)
{
	sdata = "";
	string jdata;
	if (read_file_into_string(json_fname, jdata) < 0) return -1;
	istringstream s(jdata);

	MLOG("jdata is %s\n", jdata.c_str());

	pt::ptree ptree;
	pt::read_json(s, ptree);

	for (auto &p : ptree.get_child("signals")) {
		auto& sig = p.second;
		string sig_name = sig.get<string>("code");
		
		vector<long long> ts;
		vector<float> vals;
		for (auto &data :sig.get_child("data")) {
			for (auto &t : data.second.get_child("timestamp"))
				ts.push_back(t.second.get_value<long long>());
			for (auto &v : data.second.get_child("value"))
				vals.push_back(v.second.get_value<float>());
		}
		sdata += sig_name;
		if (vals.size() > 0) {
			sdata += ":";
			for (auto &v : vals) {
				sdata += to_string(v);
				if (&v != &vals.back()) sdata += ",";
			}
		}
		if (ts.size() > 0) {
			sdata += ":";
			for (auto &t : ts) {
				sdata += to_string(t);
				if (&t != &ts.back()) sdata += ",";
			}
		}
		sdata += ";";
	}

	MLOG("Converted to string: %s\n", sdata.c_str());
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
	string sdata = vm["test_data"].as<string>();
	if (vm["json_data"].as<string>() != "") input_json_to_string(vm["json_data"].as<string>(), sdata);
	load_algomarker_from_string(test_am, 1, sdata);
	
	// Calculate
	char *stypes[] ={ "Raw" };
	vector<int> _pids ={ 1 };
	vector<long long> _timestamps ={ (long long)vm["date"].as<long long>() };
	AMRequest *req;
	MLOG("Creating Request\n");
	int req_create_rc = AM_API_CreateRequest("test_request", stypes, 1, &_pids[0], &_timestamps[0], (int)_pids.size(), &req);
	if (req == NULL)
		MLOG("ERROR: Got a NULL request rc %d!!\n", req_create_rc);
	AMResponses *resp;

	// calculate scores
	MLOG("Before Calculate\n");
	AM_API_CreateResponses(&resp);
	AM_API_Calculate(test_am, req, resp);
	

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
	float _scr;
	int pid;
	long long ts;
	char *_scr_type;
	AMResponse *response;
	for (int i=0; i<n_resp; i++) {
		MLOG("Getting response no. %d\n", i);
		AM_API_GetResponseAtIndex(resp, i, &response);

		n_shared_msgs = 0;
		AM_API_GetResponseMessages(response, &n_shared_msgs, &shared_codes, &shared_args);
		MLOG("Response Messages: %d\n", n_shared_msgs);
		for (int i=0; i<n_shared_msgs; i++) {
			MLOG("Response message %d : [%d] %s\n", i, shared_codes[i], shared_args[i]);
		}

		n_shared_msgs = 0;
		AM_API_GetScoreMessages(response, 0, &n_shared_msgs, &shared_codes, &shared_args);
		MLOG("Score 0 Messages: %d\n", n_shared_msgs);
		for (int i=0; i<n_shared_msgs; i++) {
			MLOG("Score 0 message %d : [%d] %s\n", i, shared_codes[i], shared_args[i]);
		}

		AM_API_GetResponsePoint(response, &pid, &ts);
		int resp_rc = AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
		MLOG("resp_rc = %d\n", resp_rc);
		MLOG("i %d , pid %d ts %lld scr %f\n", i, pid, ts, _scr);
		MLOG("ptr for _scr_type %x\n", _scr_type);
		if (_scr_type != NULL) MLOG("_scr_type %s\n", _scr_type);
		//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, n_scr, _scr, _scr_type);
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

//--------------------------------------------------------------------------------------------------------------------------------
int simple_egfr_test()
{
	// init AM
	AlgoMarker *test_am;

	if (AM_API_Create((int)AM_TYPE_SIMPLE_EXAMPLE_EGFR, &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating test algomarker\n");
		return -1;
	}

	MLOG("Name is %s\n", test_am->get_name());

	int load_rc;
	if ((load_rc = AM_API_Load(test_am, "AUTO") != AM_OK_RC)) {
		MERR("ERROR: Failed loading algomarker %s  , rc: %d\n", test_am->get_name(), load_rc);
		return -1;
	}
	MLOG("Algomarker %s was loaded\n", test_am->get_name());


	// Load Data
	vector<long long> times ={ 20160101 };
	vector<float> vals ={ 2.0 };
	AM_API_AddData(test_am, 1, "Creatinine", (int)times.size(), &times[0], (int)vals.size(), &vals[0]);
	/*vector<float>*/ vals ={ 55 };
	AM_API_AddData(test_am, 1, "Age", 0, NULL, (int)vals.size(), &vals[0]);
	/*vector<float>*/ vals ={ 1 };
	AM_API_AddData(test_am, 1, "GENDER", 0, NULL, (int)vals.size(), &vals[0]);

	// Calculate
	char *stypes[] ={ "Raw" };
	vector<int> _pids ={ 1 };
	vector<long long> _timestamps ={ 20160101 };
	AMRequest *req;
	MLOG("Creating Request\n");
	int req_create_rc = AM_API_CreateRequest("test_request", stypes, 1, &_pids[0], &_timestamps[0], (int)_pids.size(), &req);
	if (req == NULL)
		MLOG("ERROR: Got a NULL request rc %d!!\n", req_create_rc);
	AMResponses *resp;

	// calculate scores
	MLOG("Before Calculate\n");
	AM_API_CreateResponses(&resp);
	AM_API_Calculate(test_am, req, resp);


	// Shared messages
	int n_shared_msgs;
	int *shared_codes;
	char **shared_args;
	AM_API_GetSharedMessages(resp, &n_shared_msgs, &shared_codes, &shared_args);
	MLOG("Shared Messages: %d\n", n_shared_msgs);
	for (int i=0; i<n_shared_msgs; i++) {
		MLOG("Shared message %d : [%d] %s\n", i, shared_codes[i], shared_args[i]);
	}

	// print result
	int n_resp = AM_API_GetResponsesNum(resp);
	MLOG("Got %d responses\n", n_resp);
	float _scr;
	int pid;
	long long ts;
	char *_scr_type;
	AMResponse *response;
	for (int i=0; i<n_resp; i++) {
		MLOG("Getting response no. %d\n", i);

		AM_API_GetResponseAtIndex(resp, i, &response);
		AM_API_GetResponsePoint(response, &pid, &ts);
		int resp_rc = AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
		MLOG("_scr %f _scr_type %s\n", _scr, _scr_type);
		MLOG("resp_rc = %d\n", resp_rc);
		MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);
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

	MLOG("Finished egfr_test()\n");

	return 0;
}


//----------------------------------------------------------------------------------------
int test_data_api(po::variables_map &vm) {

	MedRepository rep;
	RepositoryHandle *rep_h;
	SignalDataHandle *sdh;

	string fname = vm["rep"].as<string>();
	int pid = vm["pid"].as<int>();
	string sig = vm["sig"].as<string>();

	if (rep.read_all(fname, { pid }, { sig }) < 0) return -1;
	
	char *sigs = (char *)sig.c_str();
	if (DATA_API_RepositoryHandle_Create(&rep_h, (char *)fname.c_str(), &pid, 1, &sigs, 1) < 0) return -1;
	if (DATA_API_SignalDataHandle_Create(&sdh) < 0) return -1;

	int len;
	DATA_API_ReadData(rep_h, pid, (char *)sig.c_str(), sdh, &len);

	UniversalSigVec usv;
	rep.uget(pid, sig, usv);

	MLOG("len %d ulen %d\n", len, usv.len);

	if (len != usv.len) { MERR("ERROR: different lengths\n"); }

	for (int i=0; i<len; i++) {

		for (int t=0; t<usv.n_time_channels(); t++) {
			int time;
			DATA_API_GetTime(sdh, i, t, &time);
			MLOG("Time %d,%d : %d , %d\n", i, t, usv.Time(i, t), time);
			if (time != usv.Time(i, t)) MERR("ERROR different times\n");
		}

		for (int v=0; v<usv.n_val_channels(); v++) {
			float val;
			DATA_API_GetVal(sdh, i, v, &val);
			MLOG("Time %d,%d : %f , %f\n", i, v, usv.Val(i, v), val);
			if (val != usv.Val(i, v)) MERR("ERROR different vals\n");
		}

	}

	return 0;
}


//----------------------------------------------------------------------------------------
int test_kp_format(po::variables_map &vm) {

	MLOG("Testing model %s in kp format\n", vm["model"].as<string>().c_str());
	
	// read model file
	MedModel model;
	if (model.read_from_file(vm["model"].as<string>()) < 0) {
		MERR("Could not read model file %s\n", vm["model"].as<string>().c_str());
		return -1;
	}
	MLOG("Read model file %s\n", vm["model"].as<string>().c_str());


	// read demographic file
	vector<vector<string>> raw_demographics;
	if (read_text_file_cols(vm["kp_demographic"].as<string>(), " \t", raw_demographics) < 0) {
		MERR("Could not read demographics file %s\n", vm["kp_demographic"].as<string>().c_str());
		return -1;
	}
	MLOG("Read %d lines from demographics file %s\n", raw_demographics.size(), vm["kp_demographic"].as<string>().c_str());
	unordered_map<string, int> name_to_pid;
	unordered_map<int, string> pid_to_name;

	int curr_id = 1000000;
	for (auto &v : raw_demographics)
		if (v.size() >= 2) {
			if (name_to_pid.find(v[0]) != name_to_pid.end()) {
				MERR("ERROR: Got the same pid %s twice in demographics file.\n", v[0].c_str());
				return -1;
			}
			name_to_pid[v[0]] = curr_id;
			pid_to_name[curr_id] = v[0];
			curr_id++;
		}


	// read data file 
	vector<vector<string>> raw_data;
	if (read_text_file_cols(vm["kp_data"].as<string>(), " \t", raw_data) < 0) {
		MERR("Could not read lab tests data file %s\n", vm["kp_data"].as<string>().c_str());
		return -1;
	}
	MLOG("Read %d lines from data file %s\n", raw_data.size(), vm["kp_data"].as<string>().c_str());

	// read scores file
	vector<vector<string>> raw_scores;
	if (read_text_file_cols(vm["kp_scores"].as<string>(), " \t", raw_scores) < 0) {
		MERR("Could not read scores file %s\n", vm["kp_scores"].as<string>().c_str());
		return -1;
	}
	MLOG("Read %d lines from scores file %s\n", raw_scores.size(), vm["kp_scores"].as<string>().c_str());


	// read code file
	vector<vector<string>> raw_codes;
	if (read_text_file_cols(vm["kp_codes"].as<string>(), " \t", raw_codes) < 0) {
		MERR("Could not read lab codes file %s\n", vm["kp_codes"].as<string>().c_str());
		return -1;
	}
	MLOG("Read %d lines from codes file %s\n", raw_codes.size(), vm["kp_codes"].as<string>().c_str());

	// prepare codes map
	unordered_map<string, string> codes;
	for (auto &v : raw_codes) {
		if (v.size() > 1)
			codes[v[0]] = v[1];
	}

	// prepare MedSamples
	MedSamples samples;
	for (auto &v : raw_scores)
		if (v.size() >= 2) {
			if (name_to_pid.find(v[0]) == name_to_pid.end()) {
				MERR("ERROR: pid %s appears in scores file but not in demographics file\n", v[0].c_str());
				return -1;
			}
			samples.insertRec(name_to_pid[v[0]], stoi(v[1]));
		}
	samples.normalize();
	MLOG("Prepared MedSamples\n");

	// Read (empty) repository 
	MedPidRepository rep;
	if (rep.MedRepository::init(vm["rep"].as<string>()) < 0) {
		MERR("Could not read repository definitions from %s\n", vm["rep"].as<string>().c_str());
		return -1;
	}
	MLOG("Read repository definitions from %s\n", vm["rep"].as<string>().c_str());

	// move to in_mem mode and push all data into it
	rep.switch_to_in_mem_mode();

	// load BYEAR and GENDER
	for (auto &v : raw_demographics) {
		if (v.size() >= 3) {
			int pid = name_to_pid[v[0]];
			float byear = stof(v[1]);
			float gender = (float)1.0;
			if (v[2] == "F") gender = (float)2.0;
			rep.in_mem_rep.insertData(pid, "BYEAR", NULL, &byear, 0, 1);
			rep.in_mem_rep.insertData(pid, "GENDER", NULL, &gender, 0, 1);
		}
	}
	MLOG("Loaded demographics into repository\n");

	// load lab tests
	int n = 0;
	for (auto &v : raw_data) {
		if (v.size() >= 4) {
			if (name_to_pid.find(v[0]) == name_to_pid.end()) {
				MERR("ERROR: pid %s appears in data file but not in demographics file ... ignoring it\n", v[0].c_str());
				continue;
			}
			if (codes.find(v[1]) == codes.end()) {
				MERR("ERROR: code %s in data file is undefined... ignoring it\n", v[1].c_str());
				continue;
			}

			int pid = name_to_pid[v[0]];
			int date = stoi(v[2]);
			float val = stof(v[3]);
			rep.in_mem_rep.insertData(pid, codes[v[1]].c_str(), &date, &val, 1, 1);
			n++;
			if (n % 500000 == 0)
				MLOG("Loaded %d lab tests into in_mem_rep\n", n);
		}
	}
	MLOG("Loaded %d data lines into repository\n", n);


	// sort in_mem , before using it
	rep.in_mem_rep.sortData();
	MLOG("Repository ready\n");

	// apply model
	model.apply(rep, samples);
	MLOG("Model applied\n");

	// write results to output file
	string sout;
	int no = 0;
	int ni = 0;
	for (auto &ids : samples.idSamples) {
		ni++;
		for (auto &s : ids.samples) {
			sout += pid_to_name[s.id] + " " + to_string(s.time) + " " + to_string(s.prediction[0]) + "\n";
			no++;
		}
	}

	
	if (write_string(vm["kp_fout"].as<string>(), sout) < 0) {
		MERR("Could not write results to output file %s\n", vm["kp_fout"].as<string>().c_str());
		return -1;
	}

	MLOG("Wrote %d predictions (for %d different ids) to file %s\n", no, ni, vm["kp_fout"].as<string>().c_str());

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

	if (vm.count("egfr_test"))
		return simple_egfr_test();

	if (vm.count("data_api_test"))
		return test_data_api(vm);

	if (vm.count("kp_test"))
		return test_kp_format(vm);

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
	//sigs ={ "BYEAR","GENDER","Glucose", "BMI" ,"WBC","Triglycerides","ALT","RBC","GENDER","HDL","Na","HbA1C","Weight","Drug" }; // DEBUG
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

#if 0

	//MLOG("SIGNALS IN REP ===========> :: \n");
	//vector<string> print_sigs ={ "Drug" };
	//for (int i=0; i<pids.size(); i++) {
	//	for (int j=0; j< print_sigs.size(); j++)
	//		rep.print_data_vec(pids[i], print_sigs[j]);
	//}

	model.write_feature_matrix("direct.csv");

#endif


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

	int print_msgs = (vm.count("print_msgs")) ? 1 : 0;
	if (vm.count("single"))
		get_preds_from_algomarker_single(test_am, vm["rep"].as<string>(), rep, model, samples2, pids, sigs, res2, res1, print_msgs);
	else
		get_preds_from_algomarker(test_am, vm["rep"].as<string>(), rep, model, samples2, pids, sigs, res2, print_msgs);
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
