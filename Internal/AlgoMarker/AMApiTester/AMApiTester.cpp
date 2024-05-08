//
// Test Program to the Dll API
//
// General Plan :
// 
// Compare same data/model/points prediction using the infrastructure directly AND using the DLL.
//

#define AM_DLL_IMPORT

#include <AlgoMarker/AlgoMarker/AlgoMarker.h>
#include <AlgoMarker/DynAMWrapper/DynAMWrapper.h>
#include <AlgoMarker/CommonTestingTools/CommonTestingTools.h>
#include <AlgoMarker/CommonTestingTools/DataLoader.h>
#include <AlgoMarker/CommonTestingTools/internal_am.h>

#include <string>
#include <iostream>
#include <boost/program_options.hpp>


#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedIO/MedIO/MedIO.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>




#ifdef __linux__ 
#define DEFAULT_AM_LOCATION "${MR_ROOT}/Libs/Internal/AlgoMarker/Linux/Release/libdyn_AlgoMarker.so"
#elif _WIN32
#define DEFAULT_AM_LOCATION "%MR_ROOT%\\Libs\\Internal\\AlgoMarker\\x64\\ReleaseDLL\\AlgoMarker.dll"
#endif

#include <climits>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;
namespace po = boost::program_options;
namespace pt = boost::property_tree;



using namespace CommonTestingTools;

//all program parameters organized in a class 
class testing_context {
public:
	string med_csv_file;
	string am_csv_file;
	bool force_add_data;
	vector<string> ignore_sig;
	string msgs_file;
	ofstream msgs_stream;
	string rep, samples, model;
	string scores_file;
	bool score_to_date_format_is_samples;
	bool test_am = true;
	bool test_med = true;
	string amfile;
	bool egfr_test;
	string amconfig;
	bool print_msgs;
	bool single;
	string med_res_file;
	string am_res_file;
	ofstream json_reqfile_stream;
	string json_reqfile;
	ofstream json_resfile_stream;
	string json_resfile;

	int read_from_var_map(po::variables_map vm) {
		med_csv_file = vm["med_csv_file"].as<string>();
		am_csv_file = vm["am_csv_file"].as<string>();
		force_add_data = vm.count("force_add_data") != 0;
		if (vm["ignore_sig"].as<string>() != "")
			split(ignore_sig, vm["ignore_sig"].as<string>(), boost::is_any_of(","));
		msgs_file = (vm["msgs_file"].as<string>());
		if (msgs_file != "") {
			msgs_stream.open(msgs_file);
			msgs_stream << "msg_type\tpid\tdate\ti\tj\tk\tcode\tmsg_text" << endl;
		}
		rep = vm["rep"].as<string>();
		samples = vm["samples"].as<string>();
		model = vm["model"].as<string>();
		if (vm.count("only_am")) test_med = false;
		if (vm.count("only_med")) test_am = false;
		amfile = vm["amfile"].as<string>();
		egfr_test = (vm.count("egfr_test") != 0);
		amconfig = vm["amconfig"].as<string>();
		print_msgs = (vm.count("print_msgs") != 0);
		single = (vm.count("single") != 0);
		med_res_file = vm["med_res_file"].as<string>();
		am_res_file = vm["am_res_file"].as<string>();
		json_reqfile = vm["json_reqfile"].as<string>();
		json_resfile = vm["json_resfile"].as<string>();
		if (json_reqfile != "") {
			json_reqfile_stream.open(json_reqfile);
		}
		if (json_resfile != "") {
			json_resfile_stream.open(json_resfile);
		}

		return 0;
	}
};


//=========================================================================================================
int read_run_params(int argc, char *argv[], po::variables_map& vm) {
	po::options_description desc("Program options");

	try {
		desc.add_options()
			("help", "Produce help message")
			("rep", po::value<string>()->default_value("/home/Repositories/THIN/thin_mar2017/thin.repository"), "Repository file name")
			("amfile", po::value<string>()->default_value(expandEnvVars(DEFAULT_AM_LOCATION)), "AlgoMarker .so/.dll file")
			("am_res_file", po::value<string>()->default_value(""), "File name to save AlgoMarker API results to")
			("med_res_file", po::value<string>()->default_value(""), "File name to save Medial API results to")
			("med_csv_file", po::value<string>()->default_value(""), "file to write Med API feature matrix after apply")
			("am_csv_file", po::value<string>()->default_value(""), "file to write AM API feature matrix after apply")
			("samples", po::value<string>()->default_value(""), "MedSamples file to use")
			("model", po::value<string>()->default_value(""), "model file to use")
			("amconfig", po::value<string>()->default_value(""), "AlgoMarker configuration file")
			("msgs_file", po::value<string>()->default_value(""), "file to save messages codes to")
			("ignore_sig", po::value<string>()->default_value(""), "Comma-seperated list of signals to ignore, data from these signals will bot be sent to the am")
			("single", "Run test in single mode, instead of the default batch")
			("print_msgs", "Print algomarker messages when testing batches or single (direct test always prints them)")
			("only_am", "Test only the AlgoMarker API with no compare")
			("only_med", "Test only the direct Medial API with no compare")
			("egfr_test", "Test simple egfr algomarker")
			("force_add_data", "Force using the AddData() API call instead of the AddDataStr()")
			("json_reqfile", po::value<string>()->default_value(""), "JSON request file name")
			("json_resfile", po::value<string>()->default_value(""), "JSON result file name")
			;

		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		MLOG("=========================================================================================\n");
		MLOG("Command Line:");
		for (int i = 0; i < argc; i++) MLOG(" %s", argv[i]);
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
	if (vm.count("help")) {
		cerr << desc << "\n";
	}

	return 0;
}

//--------------------------------------------------------------------------------------------------------------------------------
int simple_egfr_test()
{
	// init AM
	AlgoMarker *test_am;

	if (DynAM::AM_API_Create((int)AM_TYPE_SIMPLE_EXAMPLE_EGFR, &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating test algomarker\n");
		return -1;
	}


	int load_rc;
	if ((load_rc = DynAM::AM_API_Load(test_am, "AUTO") != AM_OK_RC)) {
		MERR("ERROR: Failed loading algomarker , rc: %d\n", load_rc);
		return -1;
	}
	MLOG("Algomarker was loaded\n");


	// Load Data
	vector<long long> times = { 20160101 };
	vector<float> vals = { 2.0 };
	DynAM::AM_API_AddData(test_am, 1, "Creatinine", (int)times.size(), &times[0], (int)vals.size(), &vals[0]);
	/*vector<float>*/ vals = { 55 };
	DynAM::AM_API_AddData(test_am, 1, "Age", 0, NULL, (int)vals.size(), &vals[0]);
	/*vector<float>*/ vals = { 1 };
	DynAM::AM_API_AddData(test_am, 1, "GENDER", 0, NULL, (int)vals.size(), &vals[0]);

	// Calculate
	char *stypes[] = { "Raw" };
	vector<int> _pids = { 1 };
	vector<long long> _timestamps = { 20160101 };
	AMRequest *req;
	MLOG("Creating Request\n");
	int req_create_rc = DynAM::AM_API_CreateRequest("test_request", stypes, 1, &_pids[0], &_timestamps[0], (int)_pids.size(), &req);
	if (req == NULL)
		MLOG("ERROR: Got a NULL request rc %d!!\n", req_create_rc);
	AMResponses *resp;

	// calculate scores
	MLOG("Before Calculate\n");
	DynAM::AM_API_CreateResponses(&resp);
	DynAM::AM_API_Calculate(test_am, req, resp);


	// Shared messages
	int n_shared_msgs;
	int *shared_codes;
	char **shared_args;
	DynAM::AM_API_GetSharedMessages(resp, &n_shared_msgs, &shared_codes, &shared_args);
	MLOG("Shared Messages: %d\n", n_shared_msgs);
	for (int i = 0; i < n_shared_msgs; i++) {
		MLOG("Shared message %d : [%d] %s\n", i, shared_codes[i], shared_args[i]);
	}

	// print result
	int n_resp = DynAM::AM_API_GetResponsesNum(resp);
	MLOG("Got %d responses\n", n_resp);
	float _scr;
	int pid;
	long long ts;
	char *_scr_type;
	AMResponse *response;
	for (int i = 0; i < n_resp; i++) {
		MLOG("Getting response no. %d\n", i);

		DynAM::AM_API_GetResponseAtIndex(resp, i, &response);
		DynAM::AM_API_GetResponsePoint(response, &pid, &ts);
		int resp_rc = DynAM::AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
		MLOG("_scr %f _scr_type %s\n", _scr, _scr_type);
		MLOG("resp_rc = %d\n", resp_rc);
		MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);
	}


	// print error messages

	// AM level
	int n_msgs, *msg_codes;
	char **msgs_errs;
	DynAM::AM_API_GetSharedMessages(resp, &n_msgs, &msg_codes, &msgs_errs);
	for (int i = 0; i < n_msgs; i++) {
		MLOG("Shared Message %d : code %d : err: %s\n", n_msgs, msg_codes[i], msgs_errs[i]);
	}


	// Dispose
	DynAM::AM_API_DisposeRequest(req);
	DynAM::AM_API_DisposeResponses(resp);
	DynAM::AM_API_DisposeAlgoMarker(test_am);

	MLOG("Finished egfr_test()\n");

	return 0;
}

vector<MedSample> apply_am_api(testing_context& t_ctx, DataLoader& d) {
	//const string& amconfig, DataLoader& d, bool print_msgs, bool single, const string& am_csv_file,bool force_add_data, ofstream& msgs_stream, vector<string> ignore_sig){
	vector<MedSample> res2;
	AlgoMarker *test_am;

	if (DynAM::AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating test algomarker\n");
		throw runtime_error("ERROR: Failed creating test algomarker\n");
	}

	// put fix here

	if (t_ctx.am_csv_file != "") {
		set_am_matrix(test_am, t_ctx.am_csv_file);
	}

	int rc = 0;
	if ((rc = DynAM::AM_API_Load(test_am, t_ctx.amconfig.c_str())) != AM_OK_RC) {
		MERR("ERROR: Failed loading algomarker with config file %s ERR_CODE: %d\n", t_ctx.amconfig.c_str(), rc);
		throw runtime_error(string("ERROR: Failed loading algomarker with config file ") + t_ctx.amconfig + " ERR_CODE: " + to_string(rc));
	}

	if (t_ctx.single)
		get_preds_from_algomarker_single(test_am, res2, t_ctx.print_msgs, d, t_ctx.force_add_data, t_ctx.msgs_stream, t_ctx.ignore_sig, t_ctx.json_reqfile_stream);
	else
		get_preds_from_algomarker(test_am, res2, t_ctx.print_msgs, d, t_ctx.force_add_data, t_ctx.msgs_stream, t_ctx.ignore_sig);

	return res2;
}

vector<MedSample> apply_med_api(MedPidRepository& rep, MedModel& model, MedSamples& samples, const string& med_csv_file, vector<string> ignore_sig) {

	if (ignore_sig.size() > 0) {
		string ppjson = "{\"pre_processors\":[{\"action_type\":\"rep_processor\",\"rp_type\":\"history_limit\",\"signal\":[";
		ppjson += string("\"") + ignore_sig[0] + "\"";
		for (int i = 1; i < ignore_sig.size(); i++)
			ppjson += string(",\"") + ignore_sig[i] + "\"";
		ppjson += "],\"delete_sig\":\"1\"}]}";
		MLOG("Adding pre_processor = \n'%s'\n", ppjson.c_str());
		model.add_pre_processors_json_string_to_model(ppjson, "");
	}

	// apply model (+ print top 50 scores)
	model.apply(rep, samples);

	if (med_csv_file != "")
		model.write_feature_matrix(med_csv_file);

	/////// REMOVE THIS
	//model.write_feature_matrix("/nas1/Work/Users/Shlomi/apply-program/generated/fmat-apply-program.csv");

	// printing
	vector<MedSample> ret;
	samples.export_to_sample_vec(ret);
	return ret;
}

void compare_results(const vector<MedSample>& res1, const vector<MedSample>& res2) {
	for (int i = 0; i < min(50, (int)res1.size()); i++) {
		MLOG("#Res1 :: pid %d time %d pred %f #Res2 pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], res2[i].id, res2[i].time, res2[i].prediction[0]);
	}

	// test results
	int nbad = 0, n_miss = 0, n_similar = 0;
	if (res1.size() != res2.size()) {
		MLOG("ERROR:: Didn't get the same number of tests ... %d vs %d\n", res1.size(), res2.size());
	}

	MLOG("Comparing %d scores\n", res1.size());
	for (int i = 0; i < res1.size(); i++) {

		if (res2[i].prediction[0] == (float)AM_UNDEFINED_VALUE) {
			n_miss++;
		}
		else if (res1[i].prediction[0] != res2[i].prediction[0]) {
			MLOG("ERROR !!!: #Res1 :: pid %d time %d pred %f #Res2 pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], res2[i].id, res2[i].time, res2[i].prediction[0]);
			nbad++;
		}
		else
			n_similar++;

	}


	MLOG(">>>>>TEST1: test DLL API batch: total %d : n_similar %d : n_bad %d : n_miss %d\n", res1.size(), n_similar, nbad, n_miss);
	if (nbad == 0) MLOG("PASSED\n"); else MLOG("FAILED\n");

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

	if (vm.count("help")) {
		return 0;
	}
	testing_context t_ctx;
	t_ctx.read_from_var_map(vm);

	if (t_ctx.test_am)
		load_am(t_ctx.amfile.c_str());

	if (t_ctx.egfr_test)
		return simple_egfr_test();

	DataLoader d;
	vector<MedSample> res1;
	MedSamples samples, samples2;

	try {
		d.load(t_ctx.rep,
			t_ctx.model,
			t_ctx.samples);

		samples2 = samples = d.samples;
		if (t_ctx.test_med) {
			res1 = apply_med_api(d.rep, d.model, samples, t_ctx.med_csv_file, t_ctx.ignore_sig);
			for (int i = 0; i < min(50, (int)res1.size()); i++) {
				MLOG("#Res1 :: pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0]);
			}
		}
	}
	catch (runtime_error &e) {
		cout << "(EE) Error: " << e.what() << "\n";
		return -1;
	}


	//fake failed
	//res1[3].prediction[0] = 0.1;

	//===============================================================================
	// TEST1: testing internal in_mem in a repository
	//===============================================================================
	vector<MedSample> res2;
	try {
		if (t_ctx.test_am)
			res2 = apply_am_api(t_ctx, d);
	}
	catch (runtime_error &e) {
		return -1;
	}

	if (t_ctx.test_med && t_ctx.med_res_file != "")
		save_sample_vec(res1, t_ctx.med_res_file);
	if (t_ctx.test_am && t_ctx.am_res_file != "")
		save_sample_vec(res2, t_ctx.am_res_file);

	if (t_ctx.test_am && t_ctx.test_med)
		compare_results(res1, res2);

	if (t_ctx.msgs_file != "")
		t_ctx.msgs_stream.close();

	if (t_ctx.json_reqfile_stream.is_open()) {
		t_ctx.json_reqfile_stream.close();
	}
	if (t_ctx.json_resfile_stream.is_open()) {
		t_ctx.json_resfile_stream.close();
	}

	return 0;
}

//
// keep command line:
//
// typical test:
//  ../Linux/Release/SOAPITester --single --print_msgs --rep /home/Repositories/THIN/thin_jun2017/thin.repository --samples ./Build/test/test.samples --model ./Build/test/Partial_All_S6.model --amconfig ./Build/test/pre2d.amconfig
//
// old typical test:
// Linux/Release/DllAPITester --model /nas1/Work/Users/Avi/Diabetes/order/pre2d/runs/partial/pre2d_partial_S6.model --samples test_100k.samples --amconfig /nas1/Work/Users/Avi/AlgoMarkers/pre2d/pre2d.amconfig
//
// ./Linux/Release/AMApiTester --generate_data --generate_data_outfile /tmp/out2.txt --rep /home/Repositories/THIN/thin_final/thin.repository --model /nas1/Products/Pre2D/FrozenVersions/1.0.0.9/pre2d.model --samples /nas1/Work/Users/Avi/GAN/prep_pre2d_mat/pre2d_check_bw.samples
//
