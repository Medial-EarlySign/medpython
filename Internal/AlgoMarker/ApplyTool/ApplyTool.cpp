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
//#include "internal_am.h"
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
	string rep, samples, model, generate_data_outfile, generate_data_cat_prefix;
	bool generate_data;
	bool generate_data_force_cat_prefix;
	bool apply;
	string apply_outfile, apply_repdata, apply_repdata_jsonreq;
	string apply_amconfig;
	string scores_file;
	bool score_to_date_format_is_samples;
	string apply_dates_to_score;
	bool test_am = true;
	bool test_med = true;
	string amfile;
	bool egfr_test;
	bool print_msgs;
	bool single;
	ofstream json_reqfile_stream;
	string json_reqfile;
	bool convert_reqfile_to_data;
	string convert_reqfile_to_data_infile;
	string convert_reqfile_to_data_outfile;
	bool extended_score;

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
		generate_data = (vm.count("generate_data") != 0);
		if (generate_data) {
			if (vm["rep"].as<string>() == "" || vm["samples"].as<string>() == "" || vm["model"].as<string>() == "" || vm["generate_data_outfile"].as<string>() == "")
			{
				std::cerr << "Missing argument, Please specify --rep, --samples, --model, --generate_data_outfile.\n";
				return -1;
			}
			generate_data_outfile = vm["generate_data_outfile"].as<string>();
			generate_data_cat_prefix = vm["generate_data_cat_prefix"].as<string>();
			generate_data_force_cat_prefix = (vm.count("generate_data_force_cat_prefix") != 0);
		}
		apply = (vm.count("apply") != 0);
		apply_outfile = vm["apply_outfile"].as<string>();
		apply_repdata = vm["apply_repdata"].as<string>();
		apply_repdata_jsonreq = vm["apply_repdata_jsonreq"].as<string>();
		apply_amconfig = vm["apply_amconfig"].as<string>();
		apply_dates_to_score = vm["apply_dates_to_score"].as<string>();
		if (apply || (vm.count("apply_amconfig") && apply_amconfig != "")) {
			if (rep == "" ||
				(samples == "" && apply_dates_to_score == "") ||
				model == "" ||
				apply_outfile == "" ||
				(apply_repdata == "" && apply_repdata_jsonreq == ""))
			{
				MERR("Missing arguments, Please specify --rep, --model, --apply_outfile, --apply_repdata, --samples (or --apply_dates_to_score).\n");
				return -1;
			}
			scores_file = vm["samples"].as<string>();
			score_to_date_format_is_samples = true;
			if (vm["apply_dates_to_score"].as<string>() != "") {
				scores_file = vm["apply_dates_to_score"].as<string>();
				score_to_date_format_is_samples = false;
			}
		}
		if (vm.count("only_am")) test_med = false;
		if (vm.count("only_med")) test_am = false;
		amfile = vm["amfile"].as<string>();
		egfr_test = (vm.count("egfr_test") != 0);
		print_msgs = (vm.count("print_msgs") != 0);
		single = (vm.count("single") != 0);
		json_reqfile = vm["json_reqfile"].as<string>();
		if (json_reqfile != "") {
			json_reqfile_stream.open(json_reqfile);
		}
		convert_reqfile_to_data = (vm.count("convert_reqfile_to_data") != 0);
		convert_reqfile_to_data_infile = vm["convert_reqfile_to_data_infile"].as<string>();
		convert_reqfile_to_data_outfile = vm["convert_reqfile_to_data_outfile"].as<string>();

		extended_score = (vm.count("extended_score") == 1);

		if (amfile == "" && apply_amconfig != "") {
			MTHROW_AND_ERR("To apply an AlgoMarker with apply_amconfig option please use --amfile to specify AlgoMarker .so file");
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
			("med_csv_file", po::value<string>()->default_value(""), "file to write Med API feature matrix after apply")
			("am_csv_file", po::value<string>()->default_value(""), "file to write AM API feature matrix after apply")
			("samples", po::value<string>()->default_value(""), "MedSamples file to use")
			("model", po::value<string>()->default_value(""), "model file to use")
			("msgs_file", po::value<string>()->default_value(""), "file to save messages codes to")
			("ignore_sig", po::value<string>()->default_value(""), "Comma-seperated list of signals to ignore, data from these signals will bot be sent to the am")
			("single", "Run test in single mode, instead of the default batch")
			("print_msgs", "Print algomarker messages when testing batches or single (direct test always prints them)")
			("force_add_data", "Force using the AddData() API call instead of the AddDataStr()")
			("generate_data", "Generate a unified repository data file for all the signals a model needs (required options: rep,samples,model)")
			("generate_data_outfile", po::value<string>()->default_value(""), "file to output the Generated unified signal file")
			("generate_data_cat_prefix", po::value<string>()->default_value(""), "If provided, prefer to convert a catogorial channel to a name/setname with given prefix")
			("generate_data_force_cat_prefix", "Ignore signals categories which do not conform to generate_data_cat_prefix")
			("apply", "Apply a model using Medial API, given --model, --rep, --apply_repdata, --samples, --apply_outfile, will write scores to output file")
			("apply_repdata", po::value<string>()->default_value(""), "Unified signal data to be used by apply action")
			("apply_repdata_jsonreq", po::value<string>()->default_value(""), "Same as apply_repdat but using JSON requests files")
			("apply_dates_to_score", po::value<string>()->default_value(""), "File containing a list of tab seperated pid and date to score to beused instead of scores for performing apply")
			("apply_amconfig", po::value<string>()->default_value(""), "Same as --apply but will use the AlgoMarker API and given amconfig")
			("apply_outfile", po::value<string>()->default_value(""), "Output file to save scores from apply")
			("convert_reqfile_to_data", "convert a json requests file to signal data file")
			("convert_reqfile_to_data_infile", po::value<string>()->default_value(""), "json file to load")
			("convert_reqfile_to_data_outfile", po::value<string>()->default_value(""), "data file name to write")
			("json_reqfile", po::value<string>()->default_value(""), "JSON request file name to generate from data being sent to AM")
			("extended_score", "use extended score api")
			;

		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		MLOG("=========================================================================================\n");
		MLOG("Command Line:");
		for (int i = 0; i<argc; i++) MLOG(" %s", argv[i]);
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

//=================================================================================================================

//--------------------------------------------------------------------------------------------------------------------------------

int generate_data(testing_context& t_ctx) {
	DataLoader l;
	l.load(t_ctx.rep, t_ctx.model, t_ctx.samples);
	l.export_required_data(t_ctx.generate_data_outfile, t_ctx.generate_data_cat_prefix, t_ctx.generate_data_force_cat_prefix);
	return 0;
}

vector<MedSample> apply_am_api(testing_context& t_ctx, DataLoader& d) {
	//const string& amconfig, DataLoader& d, bool print_msgs, bool single, const string& am_csv_file,bool force_add_data, ofstream& msgs_stream, vector<string> ignore_sig){
	vector<MedSample> res2;
	AlgoMarker *test_am;
	MLOG("(II) apply_am_api()\n");
	if (DynAM::AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating test algomarker\n");
		throw runtime_error("ERROR: Failed creating test algomarker\n");
	}

	MLOG("(II) AM_API_Create [V]\n");

	// put fix here

	if (t_ctx.am_csv_file != "") {
		set_am_matrix(test_am, t_ctx.am_csv_file);
	}

	int rc = 0;
	if ((rc = DynAM::AM_API_Load(test_am, t_ctx.apply_amconfig.c_str())) != AM_OK_RC) {
		MERR("ERROR: Failed loading algomarker with config file %s ERR_CODE: %d\n", t_ctx.apply_amconfig.c_str(), rc);
		throw runtime_error(string("ERROR: Failed loading algomarker with config file ") + t_ctx.apply_amconfig + " ERR_CODE: " + to_string(rc));
	}

	MLOG("(II) AM_API_Load [V]\n");

	if (t_ctx.single)
		get_preds_from_algomarker_single(test_am, res2, t_ctx.print_msgs, d, t_ctx.force_add_data, t_ctx.msgs_stream, t_ctx.ignore_sig, t_ctx.json_reqfile_stream, t_ctx.extended_score);
	else
		get_preds_from_algomarker(test_am, res2, t_ctx.print_msgs, d, t_ctx.force_add_data, t_ctx.msgs_stream, t_ctx.ignore_sig, t_ctx.extended_score);

	MLOG("(II) get_preds [V]\n");

	return res2;
}

vector<MedSample> apply_med_api(MedPidRepository& rep, MedModel& model, MedSamples& samples, const string& med_csv_file, vector<string> ignore_sig) {

	if (ignore_sig.size()>0) {
		string ppjson = "{\"pre_processors\":[{\"action_type\":\"rep_processor\",\"rp_type\":\"history_limit\",\"signal\":[";
		ppjson += string("\"") + ignore_sig[0] + "\"";
		for (int i = 1; i<ignore_sig.size(); i++)
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

int apply_data(testing_context& t_ctx)
{

	DataLoader l;
	MLOG("(II) Starting apply with:\n(II)   apply_repdata='%s'\n(II)   apply_repdata_jsonreq=%s\n(II)   rep='%s'\n(II)   scores_file='%s' %s\n(II)   model='%s'\n(II)   apply_outfile='%s'\n(II)   apply_amconfig='%s'\n"
		, t_ctx.apply_repdata.c_str(), t_ctx.apply_repdata_jsonreq.c_str(), t_ctx.rep.c_str(), t_ctx.scores_file.c_str(), t_ctx.score_to_date_format_is_samples ? "(samples format)" : "", t_ctx.model.c_str(), t_ctx.apply_outfile.c_str(), t_ctx.apply_amconfig.c_str());
	MLOG("(II) Loading mock repo, model and date for scoring\n");

	if (!t_ctx.score_to_date_format_is_samples) {
		l.load_samples_from_dates_to_score(t_ctx.scores_file);
		l.load(t_ctx.rep, t_ctx.model, "", false);
		MLOG("\n(II) Loading tab seperated pid+dates for scoring from %s\n", t_ctx.scores_file.c_str());
	}
	else {
		MLOG("\n(II) Loading dates for scoring from samples file %s\n", t_ctx.scores_file.c_str());
		l.load(t_ctx.rep, t_ctx.model, t_ctx.scores_file, false);
	}

	if (t_ctx.apply_repdata != "") {
		MLOG("(II) Importing data from '%s'\n", t_ctx.apply_repdata.c_str());
		l.import_required_data(t_ctx.apply_repdata);
	}
	else if (t_ctx.apply_repdata_jsonreq != "") {
		MLOG("(II) Importing json data from '%s'\n", t_ctx.apply_repdata_jsonreq.c_str());
		l.import_json_request_data(t_ctx.apply_repdata_jsonreq);
	}

	if (t_ctx.apply_amconfig == "") {
		MLOG("(II) Starting apply using Medial API\n");
		auto ret = apply_med_api(l.rep, l.model, l.samples, t_ctx.med_csv_file, t_ctx.ignore_sig);
		MLOG("(II) Saving results to %s\n", t_ctx.apply_outfile.c_str());
		save_sample_vec(ret, t_ctx.apply_outfile);
	}
	else {
		MLOG("(II) Starting apply using Algomarker API\n");
		auto ret = apply_am_api(t_ctx, l);
		MLOG("(II) Saving results to %s\n", t_ctx.apply_outfile.c_str());
		save_sample_vec(ret, t_ctx.apply_outfile);
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

	if (vm.count("help")) {
		return 0;
	}
	testing_context t_ctx;
	t_ctx.read_from_var_map(vm);

	if (t_ctx.amfile != "")
		load_am(t_ctx.amfile.c_str());

	if (t_ctx.convert_reqfile_to_data) {
		DataLoader::convert_reqfile_to_data(t_ctx.convert_reqfile_to_data_infile, t_ctx.convert_reqfile_to_data_outfile);
		return 0;
	}

	if (t_ctx.generate_data) {
		return generate_data(t_ctx);
	}
	if (t_ctx.apply || t_ctx.apply_amconfig != "") {
		return apply_data(t_ctx);
	}

	if (t_ctx.msgs_file != "")
		t_ctx.msgs_stream.close();

	if (t_ctx.json_reqfile_stream.is_open()) {
		t_ctx.json_reqfile_stream.close();
	}
//	if (t_ctx.json_resfile_stream.is_open()) {
//		t_ctx.json_resfile_stream.close();
//	}

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
