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
#include "internal_am.h"
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

static const int base_pid = 10000000;

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
	string amconfig;
	bool print_msgs;
	bool single;
	string med_res_file;
	string am_res_file;
	ofstream json_reqfile_stream;
	string json_reqfile;
	ofstream json_resfile_stream;
	string json_resfile;
	bool convert_reqfile_to_data;
	string convert_reqfile_to_data_infile;
	string convert_reqfile_to_data_outfile;

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
				(apply_repdata == "" && apply_repdata_jsonreq == "") )
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
		convert_reqfile_to_data = (vm.count("convert_reqfile_to_data") != 0);
		convert_reqfile_to_data_infile = vm["convert_reqfile_to_data_infile"].as<string>();
		convert_reqfile_to_data_outfile = vm["convert_reqfile_to_data_outfile"].as<string>();
	
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
			("amconfig" , po::value<string>()->default_value(""), "AlgoMarker configuration file")
			("msgs_file", po::value<string>()->default_value(""), "file to save messages codes to")
			("ignore_sig", po::value<string>()->default_value(""), "Comma-seperated list of signals to ignore, data from these signals will bot be sent to the am")
			("single", "Run test in single mode, instead of the default batch")
			("print_msgs", "Print algomarker messages when testing batches or single (direct test always prints them)")
			("only_am", "Test only the AlgoMarker API with no compare")
			("only_med", "Test only the direct Medial API with no compare")
			("egfr_test", "Test simple egfr algomarker")
			("force_add_data","Force using the AddData() API call instead of the AddDataStr()")
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
			("json_reqfile", po::value<string>()->default_value(""), "JSON request file name")
			("json_resfile", po::value<string>()->default_value(""), "JSON result file name")
			;

		po::store(po::parse_command_line(argc, argv, desc), vm);
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
	if (vm.count("help")) {
		cerr << desc << "\n";
	}

	return 0;
}

class DataLoader {
public:
	MedModel model;
	MedSamples samples;
	MedPidRepository rep;
	vector<int> pids;
	vector<string> sigs;
	map<int, MedIdSamples* > pid2samples;
	map<string, vector<map<int, string> > > sig_dict_cached;

	void load(const string& rep_fname, const string& model_fname, const string& samples_fname = "", bool read_signals = true) {
		// read model file
		if (model.read_from_file(model_fname) < 0) {
			MERR("FAILED reading model file %s\n", model_fname.c_str());
			throw runtime_error(string("FAILED reading model file ") + model_fname);
		}

		unordered_set<string> sigs_set;
		model.get_required_signal_names(sigs_set);

		MLOG("Required signals:");
		for (auto &sig : sigs_set) {
			MLOG(" %s", sig.c_str());
			sigs.push_back(sig);
		}
		MLOG("\n");
		if (samples_fname != "") {
			if (samples.read_from_file(samples_fname)) {
				MERR("FAILED reading samples file %s\n", samples_fname.c_str());
				throw runtime_error(string("FAILED reading samples file ") + samples_fname);
			}
		}
		MLOG("\n");
		samples.get_ids(pids);
		if (read_signals) {
			if (rep.read_all(rep_fname, pids, sigs) < 0) {
				MERR("FAILED loading pids and signals from repository %s\n", rep_fname.c_str());
				throw runtime_error(string("FAILED loading pids and signals from repository"));
			}
		}
		else {
			if (rep.MedRepository::init(rep_fname) < 0) {
				MERR("Could not read repository definitions from %s\n", rep_fname.c_str());
				throw runtime_error(string("FAILED MedRepository::init(") + rep_fname + "\")");
			}
		}
		for (auto &id : samples.idSamples)
			pid2samples[id.id] = &id;
	}

	void get_sig_dict_cached(const string& cat_prefix = "", bool force_cat_prefix = false) {
		sig_dict_cached = get_sig_dict(cat_prefix, force_cat_prefix);
	}

	map<string, vector<map<int, string> > > get_sig_dict(const string& cat_prefix = "", bool force_cat_prefix = false) {
		map<string, vector<map<int, string> > > sig_dict;
		for (auto& sig : sigs) {
			vector<map<int, string > > chan_dict;
			int section_id = rep.dict.section_id(sig);
			int sid = rep.sigs.Name2Sid[sig];
			int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
			for (int vchan = 0; vchan < n_vchan; ++vchan) {
				if (rep.sigs.is_categorical_channel(sig, vchan)) {
					map<int, string> new_dict;
					const auto& Id2Names = rep.dict.dict(section_id)->Id2Names;
					const auto& Member2Sets = rep.dict.dict(section_id)->Member2Sets;
					for (const auto& entry : Id2Names) {
						if (boost::starts_with(entry.second[0], cat_prefix)) {
							new_dict[entry.first] = entry.second[0];
							continue;
						}
						string new_ent = entry.second[0];
						if (Member2Sets.count(entry.first) != 0)
							for (const auto& setid : Member2Sets.at(entry.first)) {
								if (Id2Names.count(setid) != 0 && boost::starts_with(Id2Names.at(setid)[0], cat_prefix)) {
									if (!boost::starts_with(new_ent, cat_prefix) || new_ent.length() > Id2Names.at(setid)[0].length())
										new_ent = Id2Names.at(setid)[0];
								}
							}
						if (!force_cat_prefix || boost::starts_with(new_ent, cat_prefix))
							new_dict[entry.first] = new_ent;
					}

					chan_dict.push_back(new_dict);
				}
				else chan_dict.push_back(map<int, string>());
			}
			sig_dict[sig] = chan_dict;
		}
		return sig_dict;
	}

	map<string, vector<map<string, int>* > > get_sig_reverse_dict() {
		map<string, vector<map<string, int >* > > sig_dict;
		MLOG("(II)   Preparing signal reverse dictionary for signals\n");
		for (auto& sig : sigs) {
			//MLOG("(II)   Preparing signal dictionary for signal '%s'\n", sig.c_str());
			vector<map<string, int >* > chan_dict;
			if (rep.sigs.Name2Sid.count(sig) == 0) {
				MERR("no Name2Sid entry for signal '%s'\n", sig.c_str());
				exit(-1);
			}
			int section_id = rep.dict.section_id(sig);
			int sid = rep.sigs.Name2Sid[sig];
			int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
			for (int vchan = 0; vchan < n_vchan; ++vchan) {
				if (rep.sigs.is_categorical_channel(sig, vchan))
				{
					chan_dict.push_back(&(rep.dict.dict(section_id)->Name2Id));
				}
				else {
					chan_dict.push_back(nullptr);
				}
			}
			sig_dict[sig] = chan_dict;
		}
		return sig_dict;
	}

	void export_required_data(const string& fname, const string& cat_prefix, bool force_cat_prefix) {
		ofstream outfile(fname, ios::binary | ios::out);

		MLOG("(II) Preparing dictinaries to export\n", fname.c_str());

		auto sig_dict = get_sig_dict(cat_prefix, force_cat_prefix);

		MLOG("(II) Exporting required data to %s\n", fname.c_str());

		UniversalSigVec usv;

		for (int pid : pids) {
			for (auto &sig : sigs) {
				rep.uget(pid, sig, usv);
				for (int i = 0; i < usv.len; ++i) {
					stringstream outss;
					outss << pid << '\t';
					outss << sig;
					for (int tchan = 0, n_tchan = usv.n_time_channels(); tchan < n_tchan; ++tchan) {
						outss << '\t' << usv.Time(i, tchan);
					}
					bool ignore_line = false;
					for (int vchan = 0, n_vchan = usv.n_val_channels(); vchan < n_vchan; ++vchan) {
						if (sig_dict.at(sig)[vchan].size() == 0)
							outss << '\t' << setprecision(10) << usv.Val(i, vchan);
						else {
							if (sig_dict.at(sig)[vchan].count((int)(usv.Val(i, vchan))) != 0) {
								outss << '\t' << sig_dict.at(sig)[vchan].at((int)(usv.Val(i, vchan)));
							}
							else {
								ignore_line = true;
							}
						}
					}
					if (!ignore_line)
						outfile << outss.str() << '\n';
				}
			}
		}
		outfile.close();
	}

	static void convert_reqfile_to_data(const string& input_json_fname, const string& output_data_fname) {
		ofstream outfile(output_data_fname, ios::binary | ios::out);
		ifstream infile(input_json_fname, ios::binary | ios::in);

		MLOG("(II) Exporting required data to %s\n", output_data_fname.c_str());

		json j;
		infile >> j;

		MLOG("(II) num of requests = %d\n", j.size());

		for (int pid = 0; pid < j.size(); ++pid) {
			json j_req_signals;
			if (j[pid].count("body") != 0)
				j_req_signals = j[pid]["body"]["signals"];
			else if (j[pid].count("signals") != 0)
				j_req_signals = j[pid]["signals"];
			else throw runtime_error("Unrecognized JSON fromat");

			for (const auto& j_sig : j_req_signals)
			{
				string sig = j_sig["code"];
				for (const auto& j_data : j_sig["data"]) {
					outfile << pid + base_pid << '\t';
					outfile << sig;
					for (const auto& j_time : j_data["timestamp"]) {
						outfile << '\t' << j_time;
					}
					for (const auto& j_val : j_data["value"]) {
						if (boost::to_upper_copy(sig) == "GENDER")
							outfile << '\t' << (boost::to_upper_copy(j_val.get<string>()) == "MALE" ? "1" : "2");
						else
							outfile << '\t' << j_val.get<string>();
					}

					outfile << "\n";
				}

			}
		}
		outfile.close();
	}



	void import_required_data(const string& fname) {
		ifstream infile(fname, ios::binary | ios::in);

		auto sig_dict = get_sig_reverse_dict();
		MLOG("(II)   Switching repo to in-mem mode\n");
		string curr_line;
		rep.switch_to_in_mem_mode();

		vector<int> tchan_vec;
		vector<float> vchan_vec;
		tchan_vec.reserve(10);
		vchan_vec.reserve(10);

		MLOG("(II)   reading data in to in-mem repository\n");
		int cur_line = 0;
		while (getline(infile, curr_line)) {
			cur_line++;
			if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
				if (curr_line[curr_line.size() - 1] == '\r')
					curr_line.erase(curr_line.size() - 1);
				vector<string> fields;
				split(fields, curr_line, boost::is_any_of("\t"));
				int fields_i = 0;
				const auto& pid_str = fields[fields_i++];
				int pid = -1;
				try {
					pid = stoi(pid_str);
				}
				catch (...) {
					MERR("failed reading pid, performing stoi(\"%s\") at %s:%d\n", pid_str.c_str(), fname.c_str(), cur_line);
					exit(-1);
				}
				string sig = fields[fields_i++];
				int sid = rep.sigs.Name2Sid[sig];
				int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
				int n_tchan = rep.sigs.Sid2Info[sid].n_time_channels;
				tchan_vec.clear();
				vchan_vec.clear();
				for (int tchan = 0; tchan < n_tchan; ++tchan) {
					const auto& field_str = fields[fields_i++];
					try {
						tchan_vec.push_back(stoi(field_str));
					}
					catch (...) {
						MERR("failed reading time channel #%d, performing stoi(\"%s\") at %s:%d\n", tchan, field_str.c_str(), fname.c_str(), cur_line);
						exit(-1);
					}

				}
				for (int vchan = 0; vchan < n_vchan; ++vchan) {
					if (sig_dict[sig][vchan] == nullptr) {
						const auto& field_str = fields[fields_i++];
						try {
							vchan_vec.push_back(stof(field_str));
						}
						catch (...) {
							MERR("failed reading value channel #%d, performing stof(\"%s\") at %s:%d\n", vchan, field_str.c_str(), fname.c_str(), cur_line);
							exit(-1);
						}
					}
					else
					{
						try {
							vchan_vec.push_back((*(sig_dict.at(sig)[vchan])).at(fields[fields_i++]));
						}
						catch (...) {
							MERR("Error converting sig %s, chan %d, '%s' back to code\n", sig.c_str(), vchan, fields[fields_i - 1].c_str());
							exit(-1);
						}
					}
				}
				rep.in_mem_rep.insertData(pid, sid, tchan_vec.data(), vchan_vec.data(), n_tchan, n_vchan);
			}
		}

		rep.in_mem_rep.sortData();

		infile.close();

		////REMOVE THIS
		//export_required_data("/nas1/Work/Users/Shlomi/apply-program/generated/repdata-re-export-after-import.txt", "ATC_", true);

	}

	void import_json_request_data(const string& fname) {
		ifstream infile(fname, ios::binary | ios::in);

		auto sig_dict = get_sig_reverse_dict();
		MLOG("(II)   Switching repo to in-mem mode\n");
		rep.switch_to_in_mem_mode();

		vector<int> tchan_vec;
		vector<float> vchan_vec;
		vector<string> vchan_vec_actual;
		tchan_vec.reserve(10);
		vchan_vec.reserve(10);

		MLOG("(II)   Started reading json data to in-mem repository\n");
		bool context = false;
		int cur_rec_no = 0;
		for(json j=read_json_array_next_chunk(infile, context); j != nullptr ; j= read_json_array_next_chunk(infile, context)){
			int pid = base_pid + cur_rec_no;
			if (j.count("body"))
				j = j["body"];
			if(j.count("signals")==0 || !j.at("signals").is_array())
				MTHROW_AND_ERR("In file %s, Failed reading req #%d, no signals. request = \n'%s'\n", fname.c_str(), cur_rec_no, j.dump(1).c_str());
			for (auto j_sig : j.at("signals")) {
				string sig = j_sig.at("code");
				int sid = rep.sigs.Name2Sid[sig];
				int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
				int n_tchan = rep.sigs.Sid2Info[sid].n_time_channels;
				tchan_vec.clear();
				vchan_vec.clear();
				vchan_vec_actual.clear();

				if(j_sig.count("data") == 0 || !j_sig.at("data").is_array())
					MTHROW_AND_ERR("In file %s, Failed reading data in req #%d . signal = '%s' json = \n'%s'\n\n", fname.c_str(), cur_rec_no, sig.c_str(), j_sig.dump(1).c_str());

				for (auto d_sig : j_sig.at("data")) {

					if (n_tchan > 0 && (d_sig.count("timestamp") == 0 || (!d_sig.at("timestamp").is_array())))
						MTHROW_AND_ERR("In file %s, Failed reading timestamp in req #%d . signal = '%s' json = \n'%s'\n\n", fname.c_str(), cur_rec_no, sig.c_str(), d_sig.dump(1).c_str());
					for (int tchan = 0; tchan < n_tchan; ++tchan) {
						string field_str = to_string(d_sig.at("timestamp")[tchan].get<int>());
						try {
							tchan_vec.push_back(stoi(field_str));
						}
						catch (...) {
							MERR("failed reading time channel #%d, performing stoi(\"%s\") at %s:%d\n", tchan, field_str.c_str(), fname.c_str(), cur_rec_no);
							exit(-1);
						}

					}
					if (n_vchan > 0 && (d_sig.count("value") == 0 || (!d_sig.at("value").is_array())))
						MTHROW_AND_ERR("In file %s, Failed reading value in req #%d . signal = '%s'\n", fname.c_str(), cur_rec_no, sig.c_str());
					for (int vchan = 0; vchan < n_vchan; ++vchan) {
						string field_str = d_sig.at("value")[vchan].get<string>();
						if (boost::to_upper_copy(sig) == "GENDER") {
							if (boost::to_upper_copy(field_str) == "MALE") field_str = "1";
							else if (boost::to_upper_copy(field_str) == "FEMALE") field_str = "2";
						}
						if (sig_dict[sig][vchan] == nullptr) {
							try {
								vchan_vec.push_back(stof(field_str));
								vchan_vec_actual.push_back(field_str);
							}
							catch (...) {
								MERR("failed reading value channel #%d, performing stof(\"%s\") at %s:%d\n", vchan, field_str.c_str(), fname.c_str(), cur_rec_no);
								exit(-1);
							}
						}
						else
						{
							try {
								vchan_vec.push_back((*(sig_dict.at(sig)[vchan])).at(field_str));
								vchan_vec_actual.push_back(field_str);
							}
							catch (...) {
								MERR("Error converting sig %s, chan %d, '%s' back to code in request #%d\n", sig.c_str(), vchan, field_str.c_str(), cur_rec_no);
								exit(-1);
							}
						}
					}
				}
				/* Write .data repo for testing
				int nelem = 0;
				if (tchan_vec.size() != 0)
					nelem = tchan_vec.size() / n_tchan;
				else nelem = vchan_vec.size() / n_vchan;
				int ti = 0;
				int vi = 0;
				for(int j = 0; j < nelem; j++) {
					MLOG("%d\t%s", pid, sig.c_str());
					for (int i = 0; i < n_tchan; i++)
						MLOG("\t%d", tchan_vec[ti++]);
					for (int i = 0; i < n_vchan; i++)
						MLOG("\t%s", vchan_vec_actual[vi++].c_str());
					MLOG("\n");
				}
				*/
				rep.in_mem_rep.insertData(pid, sid, tchan_vec.data(), vchan_vec.data(), (int)tchan_vec.size(), (int)vchan_vec.size());
			}
			cur_rec_no++;
		}

		rep.in_mem_rep.sortData();

		infile.close();
		MLOG("(II)   Finished loading json data to in-mem repository\n");
	}

	int load_samples_from_dates_to_score(const string& fname)
	{
		// read scores file
		vector<vector<string>> raw_scores;
		if (read_text_file_cols(fname, " \t", raw_scores) < 0) {
			MERR("Could not read scores file %s\n", fname.c_str());
			return -1;
		}
		MLOG("(II) Read %d lines from scores file %s\n", raw_scores.size(), fname.c_str());

		// prepare MedSamples
		for (auto &v : raw_scores)
			if (v.size() >= 2) {
				samples.insertRec(stoi(v[0]), stoi(v[1]));
			}
		samples.normalize();
		MLOG("(II) Prepared MedSamples\n");
		for (auto &id : samples.idSamples)
			pid2samples[id.id] = &id;
		return 0;
	}

    void am_add_data(AlgoMarker *am, int pid, int max_date, bool force_add_data, vector<string> ignore_sig, json& json_out) {
		static bool print_once = false;
		UniversalSigVec usv;
	    int reserve_capacity = 100000;
    	vector<long long> times;
    	vector<float> vals;
		vector<bool> take_nelem;
		charpp_adaptor str_vals;
		times.reserve(reserve_capacity);
		vals.reserve(reserve_capacity);
		str_vals.reserve(reserve_capacity);
		if (!print_once) {
			print_once = true;
			MLOG("(INFO) force_add_data=%d\n", ((int)force_add_data));
			MLOG("(INFO) Will use %s API to insert data\n", (DynAM::so->addr_AM_API_AddDataStr == nullptr || force_add_data) ? "AddData()" : "AddDataStr()");
		}
		string reqId = string("req_") + to_string(pid)+ "_" + to_string(max_date);
		json_out = json({});
		json_out["body"] = {
			{"accountId", "A"},
			{"requestId", reqId.c_str() },
			{"customerId", "Earlysign"},
			{"calculator" , "LC"},
			{"signals",json::array() } 
			};
		json_out["header"] = {
			{"Accept", "application/json"},
			{"Content-Type", "application/json"}
		};
		
		for (auto &sig : sigs) {
			if (std::find(ignore_sig.begin(), ignore_sig.end(), sig) != ignore_sig.end())
				continue;
			json json_sig;
			int sid = rep.sigs.Name2Sid[sig];
			//			int section_id = rep.dict.section_id(sig);
			usv.init(rep.sigs.Sid2Info[sid]);
			rep.uget(pid, sig, usv);
			int nelem = usv.len;
			if (nelem == 0)
				continue;
			vals.clear();
			times.clear();
			take_nelem.resize(nelem);

			if (usv.n_time_channels() <= 0) {
				std::fill(take_nelem.begin(), take_nelem.end(), true);
			}
			else {
				std::fill(take_nelem.begin(), take_nelem.end(), false);
				for (int i = 0; i < nelem; i++) {
					bool take_elem = true;
					for (int j = 0; j < usv.n_time_channels(); j++) {
						if (usv.Time(i, j) > max_date) {
							take_elem = false;
							break;
						}
					}
					if (take_elem) {
						for (int j = 0; j < usv.n_time_channels(); j++) {
							times.push_back((long long)usv.Time(i, j));
						}
						take_nelem[i] = true;
					}
					else {
						if (usv.n_time_channels() == 1)
							break;
					}
				}
			}

			if (DynAM::so->addr_AM_API_AddDataStr == nullptr || force_add_data) {
				vals.clear();
				if (usv.n_val_channels() > 0) {
					for (int i = 0; i < nelem; i++) {
						if (!take_nelem[i])
							continue;
						for (int j = 0; j < usv.n_val_channels(); j++) {
							vals.push_back(usv.Val(i, j));
						}

					}
				}

				if ((times.size() > 0) || (vals.size() > 0)) {
					get_volatile_data_adaptor<long long> p_times(times);
					get_volatile_data_adaptor<float> p_vals(vals);
					DynAM::AM_API_AddData(am, pid, sig.c_str(), (int)times.size(), p_times.get_volatile_data(), (int)vals.size(), p_vals.get_volatile_data());
					json_sig = json_AddData(sig.c_str(), (int)times.size(), p_times.get_volatile_data(), (int)vals.size(), p_vals.get_volatile_data(), usv.n_time_channels(), usv.n_val_channels());
				}
			}
			else {
				str_vals.clear();
				if (usv.n_val_channels() > 0) {
					for (int i = 0; i < nelem; i++) {
						if (!take_nelem[i])
							continue;
						for (int j = 0; j < usv.n_val_channels(); j++) {
							string val = "";
							if (rep.sigs.is_categorical_channel(sid, j)) {
								val = sig_dict_cached.at(sig)[j].at((int)(usv.Val(i, j)));
							}
							else {
								val = precision_float_to_string(usv.Val(i, j));
							}
							str_vals.push_back(val);
						}
					}
				}

				if ((times.size() > 0) || (str_vals.size() > 0)) {
					get_volatile_data_adaptor<long long> p_times(times);
					DynAM::AM_API_AddDataStr(am, pid, sig.c_str(), (int)times.size(), p_times.get_volatile_data(), (int)str_vals.size(), str_vals.get_charpp());
					json_sig = json_AddDataStr(sig.c_str(), (int)times.size(), p_times.get_volatile_data(), (int)vals.size(), str_vals.get_charpp(), usv.n_time_channels(), usv.n_val_channels());
				}
			}
			if(!json_sig.is_null())
				json_out["body"]["signals"].push_back(json_sig);
		}
    }

};

//=================================================================================================================
int get_preds_from_algomarker(AlgoMarker *am, vector<MedSample> &res, bool print_msgs, DataLoader& d, bool force_add_data, ofstream& msgs_stream, vector<string> ignore_sig)
{
	DynAM::AM_API_ClearData(am);

	MLOG("=====> now running get_preds_from_algomarker()\n");
    
	MLOG("Going over %d pids\n", d.pids.size());
	d.get_sig_dict_cached();
	for (auto pid : d.pids) {
		json json_req;
        d.am_add_data(am, pid, INT_MAX, force_add_data, ignore_sig, json_req);
    }

    //ASK_AVI: Is this needed?
	//((MedialInfraAlgoMarker *)am)->set_sort(0); // getting rid of cases in which multiple data sets on the same day cause differences and fake failed tests.

	MLOG("After AddData for all batch\n");
	// finish rep loading 
	char *stypes[] ={ "Raw" };
	vector<int> _pids;
	vector<long long> _timestamps;
	vector<MedSample> _vsamp;
	d.samples.export_to_sample_vec(_vsamp);
	for (auto &s : _vsamp) {
		_pids.push_back(s.id);
		//_timestamps.push_back((long long)s.time*10000 + 1010);
		_timestamps.push_back((long long)s.time);
		//MLOG("pid %d time %lld\n", _pids.back(), _timestamps.back());
	}


	//MLOG("Before CreateRequest\n");
	// prep request
	AMRequest *req;
	int req_create_rc = DynAM::AM_API_CreateRequest("test_request", stypes, 1, &_pids[0], &_timestamps[0], (int)_pids.size(), &req);
	if (req == NULL)
		MLOG("ERROR: Got a NULL request rc = %d!!\n", req_create_rc);
	AMResponses *resp;

	// calculate scores
	MLOG("Before Calculate\n");
	DynAM::AM_API_CreateResponses(&resp);
	int calc_rc = DynAM::AM_API_Calculate(am, req, resp);
	MLOG("After Calculate: rc = %d\n", calc_rc);


	// go over reponses and pack them to a MesSample vector
	int n_resp = DynAM::AM_API_GetResponsesNum(resp);
	MLOG("Got %d responses\n", n_resp);
	res.clear();
	int pid;
	long long ts;
	char *_scr_type = NULL;
	AMResponse *response;
	for (int i=0; i<n_resp; i++) {
		//MLOG("Getting response no. %d\n", i);
		int resp_rc = DynAM::AM_API_GetResponseAtIndex(resp, i, &response);
		int n_scores;
		DynAM::AM_API_GetResponseScoresNum(response, &n_scores);
		//int resp_rc = AM_API_GetResponse(resp, i, &pid, &ts, &n_scr, &_scr, &_scr_type);
		//MLOG("resp_rc = %d\n", resp_rc);
		//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);

		DynAM::AM_API_GetResponsePoint(response, &pid, &ts);
		MedSample s;
		s.id = pid;
		if (ts > 30000000)
			s.time = (int)(ts/10000);
		else
			s.time = ts;
		if (resp_rc == AM_OK_RC && n_scores > 0) {
			float _scr;
			resp_rc = DynAM::AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
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
		DynAM::AM_API_GetSharedMessages(resp, &n_msgs, &msg_codes, &msgs_errs);
		for (int i=0; i<n_msgs; i++) {
			if (msgs_stream.is_open())
				msgs_stream << "SharedMessages\t" << 0 << "\t" << 0 << "\t" << i << "\t" << 0 << "\t" << 0 << "\t" << msg_codes[i] << "\t\"" << msgs_errs[i] << "\"" << endl;
			else
				MLOG("Shared Message %d : code %d : err: %s\n", n_msgs, msg_codes[i], msgs_errs[i]);
		}

		n_resp = DynAM::AM_API_GetResponsesNum(resp);
		for (int i=0; i<n_resp; i++) {
			AMResponse *r;
			DynAM::AM_API_GetResponseAtIndex(resp, i, &r);
			int n_scores;
			DynAM::AM_API_GetResponseScoresNum(r, &n_scores);

			DynAM::AM_API_GetResponseMessages(r, &n_msgs, &msg_codes, &msgs_errs);
			for (int k=0; k<n_msgs; k++) {
				if (msgs_stream.is_open())
					msgs_stream << "ResponseMessages\t" << 0 << "\t" << 0 << "\t" << i << "\t0\t" << k << "\t" << msg_codes[k] << "\t\"" << msgs_errs[k] << "\"" << endl;
				else
					MLOG("Response %d : Message %d : code %d : err: %s\n", i, k, msg_codes[k], msgs_errs[k]);
			}

			for (int j=0; j<n_scores; j++) {
				DynAM::AM_API_GetScoreMessages(r, j, &n_msgs, &msg_codes, &msgs_errs);
				for (int k=0; k<n_msgs; k++) { 
					if (msgs_stream.is_open())
						msgs_stream << "ScoreMessages\t" << 0 << "\t" << 0 << "\t" << i << "\t" << j << "\t" << k << "\t" << msg_codes[k] << "\t\"" << msgs_errs[k] << "\"" << endl;
					else
						MLOG("Response %d : score %d : Message %d : code %d : err: %s\n", i, j, k, msg_codes[k], msgs_errs[k]);
				}
			}
		}
	}

	DynAM::AM_API_DisposeRequest(req);
	DynAM::AM_API_DisposeResponses(resp);

	MLOG("Finished getting preds from algomarker\n");
	return 0;
}


#if 1
//=================================================================================================================
// same test, but running each point in a single mode, rather than batch on whole.
//=================================================================================================================
int get_preds_from_algomarker_single(AlgoMarker *am, vector<MedSample> &res, bool print_msgs, DataLoader& d, bool force_add_data, ofstream& msgs_stream, vector<string> ignore_sig, ofstream& json_reqfile_stream)
{

	DynAM::AM_API_ClearData(am);

	MLOG("=====> now running get_preds_from_algomarker_single()\n");
	MLOG("Going over %d samples\n", d.samples.nSamples());
	int n_tested = 0;

	MedTimer timer;
	d.get_sig_dict_cached();
	timer.start();

	bool first_json_req = true;

	json json_resp_byid;

	for (auto &id : d.samples.idSamples){
		for (auto &s : id.samples) {
			// clearing data in algomarker
			DynAM::AM_API_ClearData(am);

			// adding all data 
			json json_req;
			d.am_add_data(am, s.id, s.time, force_add_data, ignore_sig, json_req);
			if (json_reqfile_stream.is_open()) {
				json_reqfile_stream << (first_json_req ? "[\n" : ",\n");
				json_reqfile_stream << json_req.dump(1) << "\n";
				first_json_req = false;
			}


			// At this point we can send to the algomarker and ask for a score

			// a small technicality
			// ASK_AVI
			//((MedialInfraAlgoMarker *)am)->set_sort(0); // getting rid of cases in which multiple data sets on the same day cause differences and fake failed tests.

			// preparing a request
			char *stypes[] ={ "Raw" };
			long long _timestamp = (long long)s.time;

			AMRequest *req;
			int req_create_rc = DynAM::AM_API_CreateRequest("test_request", stypes, 1, &s.id, &_timestamp, 1, &req);
			if (req == NULL) {
				MLOG("ERROR: Got a NULL request for pid %d time %d rc %d!!\n", s.id, s.time, req_create_rc);
				return -1;
			}

			// create a response
			AMResponses *resp;
			DynAM::AM_API_CreateResponses(&resp);

			// Calculate
			DynAM::AM_API_Calculate(am, req, resp);
			//int calc_rc = AM_API_Calculate(am, req, resp);
			//MLOG("after Calculate: calc_rc %d\n", calc_rc);
			string reqId = string("req_") + to_string(s.id) + "_" + to_string(s.time);
			json_resp_byid[reqId]["messages"] = json::array();
			json_resp_byid[reqId]["result"] = nullptr;

			int n_resp = DynAM::AM_API_GetResponsesNum(resp);

			//MLOG("pid %d time %d n_resp %d\n", s.id, s.time, n_resp);
			// get scores
			if (n_resp == 1) {
				AMResponse *response;
				int resp_rc = DynAM::AM_API_GetResponseAtIndex(resp, 0, &response);
				int n_scores;
				DynAM::AM_API_GetResponseScoresNum(response, &n_scores);
				if (n_scores == 1) {
					float _scr;
					int pid;
					long long ts;
					char *_scr_type = NULL;
					DynAM::AM_API_GetResponsePoint(response, &pid, &ts);
					json_resp_byid[reqId]["requestId"] = string("req_") + to_string(pid) + to_string((int)ts);
					json_resp_byid[reqId]["status"] = 0;
					MedSample rs;
					rs.id = pid;
					if (ts > 30000000)
						rs.time = (int)(ts/10000);
					else
						rs.time = ts;

					if (resp_rc == AM_OK_RC && n_scores > 0) {
						resp_rc = DynAM::AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
						//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);
						rs.prediction.push_back(_scr);
						json_resp_byid[reqId]["result"] = { { "resultType", "Numeric" } };
						json_resp_byid[reqId]["result"]["value"] = _scr;
						json_resp_byid[reqId]["result"]["validTime"] = ts * 1000000;
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
				DynAM::AM_API_GetSharedMessages(resp, &n_msgs, &msg_codes, &msgs_errs);
				for (int i=0; i<n_msgs; i++) {

					if (msgs_stream.is_open())
						msgs_stream << "SharedMessages\t" << s.id << "\t" << s.time << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << msg_codes[i] << "\t\"" << msgs_errs[i] << "\"" << endl;
					else
						MLOG("pid %d time %d Shared Message %d : code %d : err: %s\n", s.id, s.time, n_msgs, msg_codes[i], msgs_errs[i]);
				}

				n_resp = DynAM::AM_API_GetResponsesNum(resp);
				for (int i=0; i<n_resp; i++) {
					AMResponse *r;
					DynAM::AM_API_GetResponseAtIndex(resp, i, &r);
					int n_scores;
					DynAM::AM_API_GetResponseScoresNum(r, &n_scores);

					DynAM::AM_API_GetResponseMessages(r, &n_msgs, &msg_codes, &msgs_errs);
					for (int k=0; k<n_msgs; k++) {
						json json_msg;
						json_msg["code"] = msg_codes[k];
						json_msg["text"] = msgs_errs[k];
						json_msg["status"] = code_to_status_tbl.at(msg_codes[k]);
						json_resp_byid[reqId]["messages"].push_back(json_msg);

						if (msgs_stream.is_open())
							msgs_stream << "ResponseMessages\t" << s.id << "\t" << s.time << "\t" << i << "\t0\t" << k << "\t" << msg_codes[k] << "\t\"" << msgs_errs[k] << "\"" << endl;
						else
							MLOG("pid %d time %d Response %d : Message %d : code %d : err: %s\n", s.id, s.time, i, k, msg_codes[k], msgs_errs[k]);
					}

					for (int j=0; j<n_scores; j++) {
						DynAM::AM_API_GetScoreMessages(r, j, &n_msgs, &msg_codes, &msgs_errs);
						for (int k=0; k<n_msgs; k++) {
							if (msgs_stream.is_open())
								msgs_stream << "ScoreMessages\t" << s.id << "\t" << s.time << "\t" << i << "\t" << j << "\t" << k << "\t" << msg_codes[k] << "\t\"" << msgs_errs[k] << "\"" << endl;
							else
								MLOG("pid %d time %d Response %d : score %d : Message %d : code %d : err: %s\n", s.id, s.time, i, j, k, msg_codes[k], msgs_errs[k]);
						}
					}
				}
			}
			// and now need to dispose responses and request
			DynAM::AM_API_DisposeRequest(req);
			DynAM::AM_API_DisposeResponses(resp);

			// clearing data in algomarker
			DynAM::AM_API_ClearData(am);

			n_tested++;
			if ((n_tested % 100) == 0) {
				timer.take_curr_time();
				double dt = timer.diff_sec();
				MLOG("Tested %d samples : time %f sec\n", n_tested, dt);
				dt = (double)n_tested/dt;
				MLOG("%f samples/sec\n", dt);
			}
		}
   }
	if (json_reqfile_stream.is_open()) {
		json_reqfile_stream  << "]";
	}

	MLOG("Finished getting preds from algomarker in a single manner\n");
	return 0;
}


#endif

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
	vector<long long> times ={ 20160101 };
	vector<float> vals ={ 2.0 };
	DynAM::AM_API_AddData(test_am, 1, "Creatinine", (int)times.size(), &times[0], (int)vals.size(), &vals[0]);
	/*vector<float>*/ vals ={ 55 };
	DynAM::AM_API_AddData(test_am, 1, "Age", 0, NULL, (int)vals.size(), &vals[0]);
	/*vector<float>*/ vals ={ 1 };
	DynAM::AM_API_AddData(test_am, 1, "GENDER", 0, NULL, (int)vals.size(), &vals[0]);

	// Calculate
	char *stypes[] ={ "Raw" };
	vector<int> _pids ={ 1 };
	vector<long long> _timestamps ={ 20160101 };
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
	for (int i=0; i<n_shared_msgs; i++) {
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
	for (int i=0; i<n_resp; i++) {
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
	for (int i=0; i<n_msgs; i++) {
		MLOG("Shared Message %d : code %d : err: %s\n", n_msgs, msg_codes[i], msgs_errs[i]);
	}


	// Dispose
	DynAM::AM_API_DisposeRequest(req);
	DynAM::AM_API_DisposeResponses(resp);
	DynAM::AM_API_DisposeAlgoMarker(test_am);

	MLOG("Finished egfr_test()\n");

	return 0;
}

int generate_data(testing_context& t_ctx) {
	DataLoader l;
	l.load(t_ctx.rep, t_ctx.model, t_ctx.samples);
	l.export_required_data(t_ctx.generate_data_outfile, t_ctx.generate_data_cat_prefix, t_ctx.generate_data_force_cat_prefix);
	return 0;
}

vector<MedSample> apply_am_api(testing_context& t_ctx, DataLoader& d){
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

    int rc=0;
	if ((rc = DynAM::AM_API_Load(test_am, t_ctx.amconfig.c_str())) != AM_OK_RC) {
		MERR("ERROR: Failed loading algomarker with config file %s ERR_CODE: %d\n", t_ctx.amconfig.c_str(), rc);
		throw runtime_error(string("ERROR: Failed loading algomarker with config file ")+ t_ctx.amconfig +" ERR_CODE: "+to_string(rc));
	}

	if (t_ctx.single)
		get_preds_from_algomarker_single(test_am, res2, t_ctx.print_msgs, d, t_ctx.force_add_data, t_ctx.msgs_stream, t_ctx.ignore_sig, t_ctx.json_reqfile_stream);
	else
		get_preds_from_algomarker(test_am, res2, t_ctx.print_msgs, d, t_ctx.force_add_data, t_ctx.msgs_stream, t_ctx.ignore_sig);

    return res2;
}

vector<MedSample> apply_med_api(MedPidRepository& rep, MedModel& model, MedSamples& samples, const string& med_csv_file, vector<string> ignore_sig){
	
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

void compare_results(const vector<MedSample>& res1, const vector<MedSample>& res2){
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

void save_sample_vec(vector<MedSample> sample_vec, const string& fname){
    MedSamples s;
    s.import_from_sample_vec(sample_vec);
    s.write_to_file(fname, 4);
}

int apply_data(testing_context& t_ctx)
{
	
	DataLoader l;
	MLOG("(II) Starting apply with:\n(II)   apply_repdata='%s'\n(II)   apply_repdata_jsonreq=%s\n(II)   rep='%s'\n(II)   scores_file='%s' %s\n(II)   model='%s'\n(II)   apply_outfile='%s'\n(II)   apply_amconfig='%s'\n"
		, t_ctx.apply_repdata.c_str(), t_ctx.apply_repdata_jsonreq.c_str(), t_ctx.rep.c_str(), t_ctx.scores_file.c_str(), t_ctx.score_to_date_format_is_samples ? "(samples format)" : "", t_ctx.model.c_str(), t_ctx.apply_outfile.c_str(), t_ctx.apply_amconfig.c_str());
	MLOG("(II) Loading mock repo, model and date for scoring\n");

	if (!t_ctx.score_to_date_format_is_samples) {
		l.load_samples_from_dates_to_score(t_ctx.scores_file);
		l.load(t_ctx.rep, t_ctx.model,"",false);
		MLOG("\n(II) Loading tab seperated pid+dates for scoring from %s\n", t_ctx.scores_file.c_str());
	}
	else { 
		MLOG("\n(II) Loading dates for scoring from samples file %s\n", t_ctx.scores_file.c_str());
		l.load(t_ctx.rep, t_ctx.model, t_ctx.scores_file,false);
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
	
	if (t_ctx.convert_reqfile_to_data) {
		DataLoader::convert_reqfile_to_data(t_ctx.convert_reqfile_to_data_infile, t_ctx.convert_reqfile_to_data_outfile);
		return 0;
	}

	if (t_ctx.generate_data) {
		return generate_data(t_ctx);
	}
	if (t_ctx.apply|| t_ctx.apply_amconfig != "") {
		return apply_data(t_ctx);
	}
    
	if(t_ctx.test_am)
		load_am(t_ctx.amfile.c_str());

	if (t_ctx.egfr_test)
		return simple_egfr_test();

    DataLoader d;
	vector<MedSample> res1;
	MedSamples samples,samples2;

    try {
        d.load(t_ctx.rep,
			t_ctx.model,
			t_ctx.samples);

        samples2 = samples = d.samples;
		if (t_ctx.test_med) {
			res1 = apply_med_api(d.rep, d.model, samples, t_ctx.med_csv_file, t_ctx.ignore_sig);
			for (int i = 0; i<min(50, (int)res1.size()); i++) {
				MLOG("#Res1 :: pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0]);
			}
		}
    }catch(runtime_error e){
		cout << "(EE) Error: " << e.what() << "\n";
      return -1;
    }
 

    //fake failed
    //res1[3].prediction[0] = 0.1;
   
	//===============================================================================
	// TEST1: testing internal in_mem in a repository
	//===============================================================================
    vector<MedSample> res2;
    try{
		if (t_ctx.test_am)
			res2 = apply_am_api(t_ctx, d);
    }catch(runtime_error e){
      return -1;
    }

    if(t_ctx.test_med && t_ctx.med_res_file!="")
        save_sample_vec(res1, t_ctx.med_res_file);
    if(t_ctx.test_am && t_ctx.am_res_file!="")
        save_sample_vec(res2, t_ctx.am_res_file);

	if(t_ctx.test_am && t_ctx.test_med)
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
