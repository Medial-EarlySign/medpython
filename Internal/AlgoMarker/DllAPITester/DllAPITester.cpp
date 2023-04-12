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
#include <boost/regex.hpp>

#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedIO/MedIO/MedIO.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <json/json.hpp>
#include <AlgoMarker/DynAMWrapper/DynAMWrapper.h>


#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

#define DYN(s) ((DynAM::initialized()) ? DynAM::s : s)

//=========================================================================================================
int read_run_params(int argc, char *argv[], po::variables_map& vm) {
	po::options_description desc("Program options");

	try {
		desc.add_options()
			("help", "produce help message")
			("rep", po::value<string>()->default_value("/home/Repositories/THIN/thin_mar2017/thin.repository"), "repository file name")
			("samples", po::value<string>()->default_value(""), "medsamples file to use")
			("allow_rep_adjustment", po::value<bool>()->default_value(false), "If true will allow fit into repository before applying model")
			("am_res_file", po::value<string>()->default_value(""), "File name to save AlgoMarker API results to")
			("model", po::value<string>()->default_value(""), "model file to use")
			("amconfig", po::value<string>()->default_value(""), "algo marker configuration file")
			("amlib", po::value<string>()->default_value(""), "algo marker .so library")
			("in_jsons", po::value<string>()->default_value(""), "input jsons to read")
			("out_jsons", po::value<string>()->default_value(""), "output jsons file : only for --single mode")
			("accountId", po::value<string>()->default_value("A"), "accountId for output json")
			("calculator", po::value<string>()->default_value("COVID19"), "calculator name for output json")
			("units", po::value<string>()->default_value("BMI,kg/m^2,Weight,kg,Height,cm,Pack_Years,pack*years,Smoking_Intensity,cigs/day,Smoking_Quit_Date,date,Smoking_Duration,years"), "accountId for output json")
			("scoreOnDate", "use to have that field in output jsons")
			("vals_json", "use to get jsons that write values for categorial signals")
			("preds_jsons", po::value<string>()->default_value(""), "output jsons preds")
			("direct_test", "split to a dedicated debug routine")
			("single", "run test in single mode, instead of the default batch")
			("print_msgs", "print algomarker messages when testing batches or single (direct test always prints them)")
			("msgs_file", po::value<string>()->default_value(""), "file to save messages codes to")
			("ignore_sig", po::value<string>()->default_value(""), "Comma-seperated list of signals to ignore, data from these signals will bot be sent to the am")
			("test_data", po::value<string>()->default_value(""), "test data for --direct_test option")
			("json_data", po::value<string>()->default_value(""), "test json data for --direct_test option")
			("signal_categ_regex", po::value<string>()->default_value(""), "path file with tab delimeted of 2 columns: signal + regex")
			("rename_signal", po::value<string>()->default_value(""), "path file with tab delimeted of 2 columns: signal name in rep + target name in model")
			("discovery", po::value<string>()->default_value(""), "path to output discovery")

			("direct_csv", po::value<string>()->default_value(""), "output matrix of the direct run (not via AM)")
			("am_csv", po::value<string>()->default_value(""), "output matrix of the run via AM")

			("dicts_config", po::value<string>()->default_value(""), "configuration file for json dictionary creation")
			("out_json_dict", po::value<string>()->default_value(""), "output json dict to test with AA as an additional load")
			("json_dict", po::value<string>()->default_value(""), "input json dict to test with AA as an additional load")
			("simple_dict", "use to generate simpler dictionaries")

			("json_req_test", "use to do a direct test using a request given as json, load data in req or in in_jsons")
			("jreq", po::value<string>()->default_value(""), "input json request")
			("jresp", po::value<string>()->default_value(""), "output json response")

			("create_jreq", "use to create a json request for given samples")
			("add_data_to_jreq", "use to create a json request for given samples")
			("jreq_defs", po::value<string>()->default_value(""), "defaulted json to start with (if not empty) when creating a request")
			("jreq_out", po::value<string>()->default_value(""), "output json request")


			;


		po::store(po::parse_command_line(argc, argv, desc), vm);
		if (vm.count("help")) {
			cerr << desc << "\n";
			exit(-1);

		}
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

	return 0;
}


//=================================================================================================================
double add_pid_to_am(AlgoMarker *am, MedPidRepository &rep, vector<string> &ignore_sig, int pid, int time, string &sig, vector<long long>& times, vector<float> &vals, char **str_values,
	const map<string, boost::regex> &signal_categ_regex)
{
	// a small technicality
	((MedialInfraAlgoMarker *)am)->set_sort(0); // getting rid of cases in which multiple data sets on the same day cause differences and fake failed tests.

	if (std::find(ignore_sig.begin(), ignore_sig.end(), sig) != ignore_sig.end())
		return 0;

	UniversalSigVec usv;

	int n_time_channels, n_val_channels, *is_categ;
	((MedialInfraAlgoMarker *)am)->get_sig_structure(sig, n_time_channels, n_val_channels, is_categ);

	rep.uget(pid, sig, usv);
	int nelem = usv.len;
	string sv;
	double t_add = 0;
	if (nelem > 0) {
		long long *p_times = &times[0];
		float *p_vals = &vals[0];
		int i_time = 0;
		int i_val = 0;

		int nelem_before = 0;

		if (usv.n_time_channels() > 0) {
			int stop_at = 0;
			int tm_ch_check = 0;
			while (stop_at < nelem && usv.Time(stop_at, tm_ch_check) <= time)
				++stop_at;
			nelem_before = stop_at;
			for (int i = 0; i < stop_at; i++) {
				//for (int j = 0; j < usv.n_time_channels(); j++) {
				for (int j = 0; j < n_time_channels; j++) {
					//MLOG("%s i=%d,j=%d,i_time=%d,stop_at=%d\n", sig.c_str(), i, j, i_time, stop_at);
					p_times[i_time] = min((long long)usv.Time(i, j), (long long)time);
					++i_time;
				}
			}
		}
		else
			p_times = NULL;

		if (usv.n_val_channels() > 0) {
			if (p_times != NULL) nelem = nelem_before;
			for (int i = 0; i < nelem; i++)
				//				for (int j = 0; j < usv.n_val_channels(); j++) {
				for (int j = 0; j < n_val_channels; j++) {
					if (str_values != NULL) {
						if (is_categ[j]) {
							int sec_id = rep.dict.section_id(sig);
							if (sec_id < 0)
								MTHROW_AND_ERR("Unknown section for sig %s\n", sig.c_str());

							const vector<string> &all_code_names = rep.dict.dicts[sec_id].Id2Names[usv.Val(i, j)];
							sv = "";
							if (signal_categ_regex.find(sig) != signal_categ_regex.end()) {
								for (size_t iii = 0; iii < all_code_names.size(); ++iii)
									if (boost::regex_match(all_code_names[iii], signal_categ_regex.at(sig))) {
										sv = all_code_names[iii];
										break;
									}
							}
							if (sv.empty())
								sv = all_code_names[0];
						}
						else
							sv = to_string(usv.Val(i, j));

						sv.copy(str_values[i_val], sv.length());
						str_values[i_val][sv.length()] = 0;
					}
					else
						p_vals[i_val] = usv.Val(i, j);
					++i_val;
				}
		}
		else
			p_vals = NULL;

		//MLOG("Adding data: pid %d time %d sig %s n_times %d n_vals %d\n", pid, time, sig.c_str(), i_time, i_val);
		MedTimer t; t.start();
		if ((i_val > 0) || (i_time > 0)) {
			if (str_values != NULL)
				DYN(AM_API_AddDataStr(am, pid, sig.c_str(), i_time, p_times, i_val, str_values));
			else
				DYN(AM_API_AddData(am, pid, sig.c_str(), i_time, p_times, i_val, p_vals));
		}
		t.take_curr_time();  t_add += t.diff_sec();
	}

	return t_add;
}

//=================================================================================================================
void print_response_msgs(AMResponses *resp, int pid, int time, ofstream &msgs_stream)
{
	// print error messages
	// AM level
	int n_msgs, *msg_codes;
	char **msgs_errs;
	DYN(AM_API_GetSharedMessages(resp, &n_msgs, &msg_codes, &msgs_errs));
	for (int i = 0; i < n_msgs; i++) {
		if (msgs_stream.is_open())
			msgs_stream << "SharedMessages\t" << pid << "\t" << time << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << msg_codes[i] << "\t\"" << msgs_errs[i] << "\"" << endl;
		else
			MLOG("pid %d time %d Shared Message %d : code %d : err: %s\n", pid, time, n_msgs, msg_codes[i], msgs_errs[i]);
	}

	int n_resp = DYN(AM_API_GetResponsesNum(resp));
	for (int i = 0; i < n_resp; i++) {
		AMResponse *r;
		DYN(AM_API_GetResponseAtIndex(resp, i, &r));
		int n_scores;
		DYN(AM_API_GetResponseScoresNum(r, &n_scores));

		DYN(AM_API_GetResponseMessages(r, &n_msgs, &msg_codes, &msgs_errs));
		for (int k = 0; k < n_msgs; k++) {
			if (msgs_stream.is_open())
				msgs_stream << "ResponseMessages\t" << pid << "\t" << time << "\t" << i << "\t0\t" << k << "\t" << msg_codes[k] << "\t\"" << msgs_errs[k] << "\"" << endl;
			else
				MLOG("pid %d time %d Response %d : Message %d : code %d : err: %s\n", pid, time, i, k, msg_codes[k], msgs_errs[k]);
		}

		for (int j = 0; j < n_scores; j++) {
			DYN(AM_API_GetScoreMessages(r, j, &n_msgs, &msg_codes, &msgs_errs));
			for (int k = 0; k < n_msgs; k++) {
				if (msgs_stream.is_open())
					msgs_stream << "ScoreMessages\t" << pid << "\t" << time << "\t" << i << "\t" << j << "\t" << k << "\t" << msg_codes[k] << "\t\"" << msgs_errs[k] << "\"" << endl;
				else
					MLOG("pid %d time %d Response %d : score %d : Message %d : code %d : err: %s\n", pid, time, i, j, k, msg_codes[k], msgs_errs[k]);
			}
		}
	}
}

//=================================================================================================================
int get_response_score_into_sample(AMResponse *response, int resp_rc, MedSample &s)
{
	int n_scores;
	int pid;
	long long ts;
	char *_scr_type = NULL;

	DYN(AM_API_GetResponseScoresNum(response, &n_scores));
	//int resp_rc = AM_API_GetResponse(resp, i, &pid, &ts, &n_scr, &_scr, &_scr_type);
	//MLOG("resp_rc = %d\n", resp_rc);
	//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);

	DYN(AM_API_GetResponsePoint(response, &pid, &ts));
	s.id = pid;
	if (ts > 30000000)
		s.time = (int)(ts / 10000);
	else
		s.time = ts;
	if (resp_rc == AM_OK_RC && n_scores > 0) {
		float _scr;
		resp_rc = DYN(AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type));
		//MLOG("i %d , pid %d ts %d scr %f %s\n", i, pid, ts, _scr, _scr_type);
		s.prediction.push_back(_scr);
	}
	else {
		s.prediction.push_back((float)AM_UNDEFINED_VALUE);
	}
	return n_scores;
}

//=================================================================================================================
int get_preds_from_algomarker(AlgoMarker *am, string rep_conf, MedPidRepository &rep, MedModel &model, MedSamples &samples, vector<int> &pids, vector<string> &sigs, vector<MedSample> &res, int print_msgs, ofstream& msgs_stream, vector<string> ignore_sig)
{
	UniversalSigVec usv;

	int max_vals = 100000;
	vector<long long> times(max_vals);
	vector<float> vals(max_vals);

	int MAX_VALS = 100000, MAX_VAL_LEN = 200;
	vector<char *>str_vals(MAX_VALS);
	vector<char> buf(MAX_VALS*MAX_VAL_LEN);
	for (int i = 0; i < MAX_VALS; i++) { str_vals[i] = &buf[i*MAX_VAL_LEN]; }

	map<string, boost::regex> signal_categ_regex;

	DYN(AM_API_ClearData(am));

	MLOG("=====> now running get_preds_from_algomarker()\n");
	MLOG("Going over %d pids\n", pids.size());
	int pid_cnt = 0;
	for (auto pid : pids) {
		for (auto &sig : sigs)	add_pid_to_am(am, rep, ignore_sig, pid, 20300101, sig, times, vals, &str_vals[0], signal_categ_regex);

		pid_cnt++;
		if (pid_cnt % 1000 == 0)
			MLOG("Loaded %d pids already\n", pid_cnt);
	}

	MLOG("After AddData for all batch\n");
	// finish rep loading 
	char *stypes[] = { "Raw" };
	vector<int> _pids;
	vector<long long> _timestamps;
	vector<MedSample> _vsamp;
	samples.export_to_sample_vec(_vsamp);
	for (auto &s : _vsamp) {
		_pids.push_back(s.id);
		_timestamps.push_back((long long)s.time);
	}

	// prep request
	AMRequest *req;
	int req_create_rc = DYN(AM_API_CreateRequest("test_request", stypes, 1, &_pids[0], &_timestamps[0], (int)_pids.size(), &req));
	if (req == NULL)
		MLOG("ERROR: Got a NULL request rc = %d!!\n", req_create_rc);
	AMResponses *resp;

	// calculate scores
	MLOG("Before Calculate\n");

	DYN(AM_API_CreateResponses(&resp));
	int calc_rc = DYN(AM_API_Calculate(am, req, resp));
	MLOG("After Calculate: rc = %d\n", calc_rc);


	// go over reponses and pack them to a MesSample vector
	int n_resp = DYN(AM_API_GetResponsesNum(resp));
	MLOG("Got %d responses\n", n_resp);
	res.clear();
	AMResponse *response;
	for (int i = 0; i < n_resp; i++) {
		//MLOG("Getting response no. %d\n", i);
		int resp_rc = DYN(AM_API_GetResponseAtIndex(resp, i, &response));
		MedSample s;
		get_response_score_into_sample(response, resp_rc, s);
		res.push_back(s);
	}


	if (print_msgs)	print_response_msgs(resp, 0, 0, msgs_stream);

	DYN(AM_API_DisposeRequest(req));
	DYN(AM_API_DisposeResponses(resp));

	MLOG("Finished getting preds from algomarker\n");
	return 0;
}


//=================================================================================================================
string float_to_string(float val) {
	stringstream ss;
	ss << std::setprecision(10) << val;
	return ss.str();
}

//=================================================================================================================

bool any_regex_matcher(const boost::regex &reg_pat, const vector<string> &nms) {
	bool res = false;
	for (size_t i = 0; i < nms.size() && !res; ++i)
		res = boost::regex_match(nms[i], reg_pat);
	return res;
}

void get_parents(int codeGroup, vector<int> &parents, int max_depth, const boost::regex &reg_pat,
	const map<int, vector<int>> &member2Sets, const map<int, vector<string>> &id_to_names) {
	vector<int> filtered_p;
	//First test codeGroup:
	if (id_to_names.find(codeGroup) == id_to_names.end())
		MTHROW_AND_ERR("get_parents - code %d wasn't found in dict\n", codeGroup);
	const vector<string> &names_cg = id_to_names.at(codeGroup);
	bool pass_regex_filter = any_regex_matcher(reg_pat, names_cg);
	if (pass_regex_filter) {
		filtered_p.push_back(codeGroup);
		parents.swap(filtered_p);
		return;
	}

	vector<int> last_parents = { codeGroup };
	if (last_parents.front() < 0)
		return; //no parents
	parents = {};


	for (size_t k = 0; k < max_depth; ++k) {
		vector<int> new_layer;
		for (int par : last_parents)
			if (member2Sets.find(par) != member2Sets.end()) {
				new_layer.insert(new_layer.end(), member2Sets.at(par).begin(), member2Sets.at(par).end());
			}
		new_layer.swap(last_parents);
		if (last_parents.empty())
			break; //no more parents to loop up

		//test if current layer has legal parents:
		for (int code : last_parents)
		{
			if (id_to_names.find(code) == id_to_names.end())
				MTHROW_AND_ERR("get_parents - code %d wasn't found in dict\n", code);
			const vector<string> &names_ = id_to_names.at(code);
			bool pass_regex_filter = any_regex_matcher(reg_pat, names_);
			if (pass_regex_filter)
				filtered_p.push_back(code);
		}
		if (!filtered_p.empty()) //found
			break;
	}

	parents.swap(filtered_p);

	//uniq:
	unordered_set<int> uniq(parents.begin(), parents.end());
	vector<int> fnal(uniq.begin(), uniq.end());
	parents.swap(fnal);
}

//=================================================================================================================
void add_sig_to_json(AlgoMarker *am, MedPidRepository &rep, vector<string> &ignore_sig,
	int pid, int time, string &sig, unordered_map<string, string> &units, json &json_out,
	const map<string, boost::regex> &signal_to_regex,
	const map<string, string> &rename_signal, int values_flag = 0)
{
	if (std::find(ignore_sig.begin(), ignore_sig.end(), sig) != ignore_sig.end())
		return;

	UniversalSigVec usv;

	int n_time_channels, n_val_channels, *is_categ;
	((MedialInfraAlgoMarker *)am)->get_sig_structure(sig, n_time_channels, n_val_channels, is_categ);

	string physical_signal = sig;
	if (rename_signal.find(sig) != rename_signal.end())
		physical_signal = rename_signal.at(sig);

	int sid = rep.sigs.sid(physical_signal);
	if (sid < 0) {
		MWARN("Signal %s doesn't exists in repository - skipping\n", physical_signal.c_str());
		return;
	}
	rep.uget(pid, sid, usv);

	// find the number of relevant elements to add to json (all those below time)
	int nelem = usv.len;
	int tm_ch_check = 0;
	if (usv.n_time_channels() > 0) {
		nelem = 0;
		while (nelem < usv.len && usv.Time(nelem, tm_ch_check) <= time) nelem++;
	}


	int section_id = rep.dict.section_id(physical_signal);

	const map<int, vector<int>> &member_to_sets = rep.dict.dict(section_id)->Member2Sets;
	const map<int, vector<string>> &id_to_names = rep.dict.dict(section_id)->Id2Names;

	if (nelem > 0) {

		// maybe we need to add unit
		string sunit = "";
		if (units.find(sig) != units.end())
			sunit = units[sig];

		//Search in signal def when empty:
		if (sunit.empty()) {
			bool has_units_in_rep = false;
			if (n_val_channels > 0) {
				sunit = ((MedialInfraAlgoMarker *)am)->get_sig_unit(sig, 0);
				has_units_in_rep = !sunit.empty();
			}
			for (int j = 1; j < n_val_channels; j++) {
				sunit += "," + ((MedialInfraAlgoMarker *)am)->get_sig_unit(sig, j);
				if (!has_units_in_rep)
					has_units_in_rep = !((MedialInfraAlgoMarker *)am)->get_sig_unit(sig, j).empty();
			}
			if (!has_units_in_rep)
				sunit = "";
		}

		// add code to json
		json json_sig = json({ { "code", sig },{ "data", json::array() } });
		if (sunit != "")
			json_sig["unit"] = sunit;

		// push elements and channels
		for (int i = 0; i < nelem; i++) {

			json json_sig_data_item = json({ { "timestamp" , json::array() },{ "value" , json::array() } });
			// add data entry to json

			// time
			for (int j = 0; j < n_time_channels; j++) {
				long long t = min((long long)time, (long long)usv.Time(i, j)); // chopping future time channels to the current time
				json_sig_data_item["timestamp"].push_back(t);
			}

			// vals
			bool no_need_to_add_data = false;
			for (int j = 0; j < n_val_channels; j++) {
				if (values_flag || rep.sigs.is_categorical_channel(sid, j) == 0)
					json_sig_data_item["value"].push_back(float_to_string(usv.Val<float>(i, j)));
				else {
					// using the first name in the dictionary
					int signal_id_val = usv.Val<int>(i, j);
					if (signal_to_regex.find(sig) == signal_to_regex.end() || n_val_channels > 1) {
						string code_txt = rep.dict.dicts[section_id].Id2Names[signal_id_val][0];
						json_sig_data_item["value"].push_back(code_txt);
						if (signal_to_regex.find(sig) != signal_to_regex.end())
							MWARN("No support for signal_to_regex on multiple channels\n");
					}
					else {
						//Currentlly supports only 1 val cahnnel
						//go up in parent till regex:
						const boost::regex &reg_filter = signal_to_regex.at(sig);
						vector<int> codes;
						get_parents(signal_id_val, codes, 1000, reg_filter, member_to_sets, id_to_names);
						//Break into all "codes":
						for (int cd_ : codes)
						{
							json json_sig_data_item_ = json({ { "timestamp" , json::array() },{ "value" , json::array() } });

							for (int j = 0; j < n_time_channels; j++) {
								long long t = min((long long)time, (long long)usv.Time(i, j)); // chopping future time channels to the current time
								json_sig_data_item_["timestamp"].push_back(t);
							}

							const vector<string> &all_code_names = rep.dict.dicts[section_id].Id2Names[cd_];
							string code_txt = "";
							for (size_t iii = 0; iii < all_code_names.size(); ++iii)
								if (boost::regex_match(all_code_names[iii], reg_filter)) {
									code_txt = all_code_names[iii];
									break;
								}
							if (code_txt.empty())
								code_txt = all_code_names[0];
							json_sig_data_item_["value"].push_back(code_txt);

							json_sig["data"].push_back(json_sig_data_item_);
						}
						no_need_to_add_data = true; // the data was added in this loop - skip the push_back of data in the end on vlaue loop
					}
				}
			}
			if (!no_need_to_add_data)
				json_sig["data"].push_back(json_sig_data_item);
		}

		json_out["body"]["signals"].push_back(json_sig);
	}
}

//=================================================================================================================
void init_output_json(po::variables_map &vm, json &json_out, int pid, int time)
{
	string reqId = string("req_") + to_string(pid) + "_" + to_string(time);
	json_out = json({});
	json_out["body"] = {
		{ "accountId", vm["accountId"].as<string>().c_str() },
	{ "requestId", reqId.c_str() },
	{ "calculator" , vm["calculator"].as<string>().c_str() },
	{ "signals",json::array() },
	{ "patient_id" , pid}
	};
	json_out["header"] = {
		{ "Accept", "application/json" },
	{ "Content-Type", "application/json" }
	};

	if (vm.count("scoreOnDate")) json_out["body"] += {"scoreOnDate", time};
}

#if 1
//=================================================================================================================
// same test, but running each point in a single mode, rather than batch on whole.
//=================================================================================================================
int get_preds_from_algomarker_single(po::variables_map &vm, AlgoMarker *am, string rep_conf,
	MedPidRepository &rep, MedModel &model, MedSamples &samples,
	vector<int> &pids, vector<string> &sigs, vector<MedSample> &res, vector<MedSample> &compare_res, ofstream& msgs_stream,
	vector<string> ignore_sig, const map<string, string> &rename_signal)
{
	UniversalSigVec usv;

	int max_vals = 100000;
	vector<long long> times(max_vals);
	vector<float> vals(max_vals);
	int print_msgs = (int)vm.count("print_msgs");
	int write_jsons = (vm["out_jsons"].as<string>() != "");
	unordered_map<string, string> units;
	map<string, boost::regex> signal_categ_regex;

	//Read for vm syntax:
	if (!vm["signal_categ_regex"].as<string>().empty()) {
		//read and prase file:
		ifstream file_reader(vm["signal_categ_regex"].as<string>());
		if (!file_reader.good())
			MTHROW_AND_ERR("Error can't read file %s\n", vm["signal_categ_regex"].as<string>().c_str());
		MLOG("Reading file %s for regex filter\n", vm["signal_categ_regex"].as<string>().c_str());
		string line;
		while (getline(file_reader, line)) {
			boost::trim(line);
			if (line.empty() || line[0] == '#')
				continue;
			vector<string> tokens;
			boost::split(tokens, line, boost::is_any_of("\t"));
			if (tokens.size() != 2)
				MTHROW_AND_ERR("Error bad file format %s. expecting 2 tokens\n",
					vm["signal_categ_regex"].as<string>().c_str());
			signal_categ_regex[tokens[0]] = boost::regex(tokens[1]);
		}
		file_reader.close();
	}

	int MAX_VALS = 100000, MAX_VAL_LEN = 200;
	vector<char *>str_vals(MAX_VALS);
	vector<char> buf(MAX_VALS*MAX_VAL_LEN);
	for (int i = 0; i < MAX_VALS; i++) { str_vals[i] = &buf[i*MAX_VAL_LEN]; }


	double t_all = 0;
	double t_add_data = 0;
	double t_add_data_api = 0;
	double t_calculate = 0;


	DYN(AM_API_ClearData(am));

	MLOG("=====> now running get_preds_from_algomarker_single()\n");
	MLOG("Going over %d samples (DllAPITester)\n", samples.nSamples());
	int n_tested = 0;

	ofstream out_json_f;
	if (write_jsons) {
		out_json_f.open(vm["out_jsons"].as<string>());
		if (!out_json_f.is_open())
			MTHROW_AND_ERR("Can't open output jsons file %s for writing\n", vm["out_jsons"].as<string>().c_str());
		out_json_f << "[" << endl;

		vector<string> f;
		boost::split(f, vm["units"].as<string>(), boost::is_any_of(",;"));
		for (int j = 0; j < f.size(); j += 2)
			units[f[j]] = f[j + 1];
	}
	bool first = true;

	MedTimer timer;	timer.start();
	MedTimer _t_all; _t_all.start();

	int smp_counts = samples.nSamples();
	MedProgress progress("get_preds_from_algomarker_single", smp_counts);
	for (auto &id : samples.idSamples)
		for (auto &s : id.samples) {

			json js;
			if (write_jsons) init_output_json(vm, js, s.id, s.time);
			// adding all data for this sample
			MedTimer _t_add_data; _t_add_data.start();
			for (auto &sig : sigs) {
				t_add_data_api += add_pid_to_am(am, rep, ignore_sig, s.id, s.time, sig, times, vals, &str_vals[0], signal_categ_regex);
				if (write_jsons) add_sig_to_json(am, rep, ignore_sig, s.id, s.time, sig, units, js, signal_categ_regex, rename_signal, (int)vm.count("vals_json"));
			}
			_t_add_data.take_curr_time(); t_add_data += _t_add_data.diff_sec();

			// preparing a request
			char *stypes[] = { "Raw" };
			long long _timestamp = (long long)s.time;

			AMRequest *req;
			int req_create_rc = DYN(AM_API_CreateRequest("test_request", stypes, 1, &s.id, &_timestamp, 1, &req));
			if (req == NULL) {
				MLOG("ERROR: Got a NULL request for pid %d time %d rc %d!!\n", s.id, s.time, req_create_rc);
				return -1;
			}

			// create a response
			AMResponses *resp;
			DYN(AM_API_CreateResponses(&resp));

			MedTimer _t_calculate; _t_calculate.start();
			DYN(AM_API_Calculate(am, req, resp)); // calculate
			_t_calculate.take_curr_time(); t_calculate += _t_calculate.diff_sec();

			int n_resp = DYN(AM_API_GetResponsesNum(resp));

			//MLOG("pid %d time %d n_resp %d\n", s.id, s.time, n_resp);

			// get scores
			if (n_resp == 1) {
				AMResponse *response;
				int resp_rc = DYN(AM_API_GetResponseAtIndex(resp, 0, &response));
				MedSample rs;
				get_response_score_into_sample(response, resp_rc, rs);
				res.push_back(rs);
			}
			else {
				MedSample rs = s;
				rs.prediction.clear();
				rs.prediction.push_back((float)AM_UNDEFINED_VALUE);
				res.push_back(rs);
			}

			if (print_msgs) print_response_msgs(resp, s.id, s.time, msgs_stream);

			// and now need to dispose responses and request
			DYN(AM_API_DisposeRequest(req));
			DYN(AM_API_DisposeResponses(resp));
			DYN(AM_API_ClearData(am)); // clearing data in algomarker


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
				dt = (double)n_tested / dt;
				MLOG("%f samples/sec\n", dt);
			}

			if (write_jsons) {
				if (!first) out_json_f << "," << endl;
				first = false;
				out_json_f << js.dump(1) << endl;
			}
			progress.update();
		}

	_t_all.take_curr_time(); t_all += _t_all.diff_sec();
	MLOG("Time report: t_all %f sec , t_add_data %f sec , t_add_data_api %f sec , t_calculate %f sec\n", t_all, t_add_data, t_add_data_api, t_calculate);


	if (write_jsons) {
		out_json_f << "]" << endl;
		out_json_f.close();
		MLOG("Wrote jsons output file %s\n", vm["out_jsons"].as<string>().c_str());
	}

	MLOG("\nFinished getting preds from algomarker in a single manner\n");
	return 0;
}
#endif

//========================================================================================
void save_sample_vec(vector<MedSample> sample_vec, const string& fname)
{
	MedSamples s;
	s.import_from_sample_vec(sample_vec);
	s.write_to_file(fname);
}

//========================================================================================
int split_file_to_jsons(string fin_name, vector<string> &jsons)
{
	jsons.clear();
	ifstream fin(fin_name);
	if (!fin.is_open()) {
		MTHROW_AND_ERR("Can't open jsons file %s\n", fin_name.c_str());
	}

	string curr_s = "";
	char ch;
	int counter = 0;
	int nch = 0;
	int in_string = 0;
	while (fin >> noskipws >> ch) {
		nch++;
		if (ch == '\"' || ch == '\'') in_string = 1 - in_string;
		if ((!in_string) && ch == '{') counter++;
		if (counter > 0) curr_s += ch;

		if ((!in_string) && ch == '}') counter--;
		if (counter == 0 && curr_s != "") {
			jsons.push_back(curr_s);
			curr_s = "";
			if (jsons.size() > 0 && jsons.size() % 1000 == 0)
				MLOG("Found %d jsons so far\n", jsons.size());
		}
		if (counter < 0) MTHROW_AND_ERR("Mismatch in {} count in file %s\n", fin_name.c_str());
		//MLOG("ch %c %d nch %d counter %d curr_s %d\n", ch, ch, nch, counter, curr_s.length());
	}

	MLOG("Read %d jsons from %s into strings (debug info: counter = %d)\n", jsons.size(), fin_name.c_str(), counter);
	//	MLOG("debug info: counter = %d nch %d len %d curr_s %s \n", counter , nch, curr_s.length(), curr_s.c_str());
	fin.close();
	return 0;
}

//========================================================================================
int add_data_to_am_from_json(AlgoMarker *am, json &js_in, const char *json_str, char **str_values, long long *p_times, int &pid, long long &score_time)
{
	json &js = js_in;
	if (js_in.find("body") != js_in.end())
		js = js_in["body"];

	//MLOG("%s\n", js["body"]["requestId"].dump().c_str());
	string reqId = js["requestId"].get<string>();
	vector<string> fields;
	boost::split(fields, reqId, boost::is_any_of("_"));
	pid = 0;
	score_time = 0;
	if (fields.size() < 3) MERR("Can't parse pid and time from reqId %s\n", reqId.c_str());
	if (fields.size() >= 3) {
		pid = stoi(fields[1]);
		score_time = (long long)stol(fields[2]);
	}

	if (js.find("scoreOnDate") != js.end())
		score_time = js["scoreOnDate"].get<long long>();

	return DYN(AM_API_AddDataByType(am, pid, DATA_JSON_FORMAT, json_str));

	//char str_values[MAX_VALS][MAX_VAL_LEN];
	for (auto &s : js["signals"]) {
		string sig = s["code"].get<string>();
		int n_time_channels, n_val_channels, *is_categ;
		((MedialInfraAlgoMarker *)am)->get_sig_structure(sig, n_time_channels, n_val_channels, is_categ);
		//MLOG("%s %d %d\n", sig.c_str(), n_time_channels, n_val_channels);
		int n_times = 0;
		int n_vals = 0;
		for (auto &d : s["data"]) {
			//MLOG("time ");
			int nt = 0;
			for (auto &t : d["timestamp"]) {
				if (nt < n_time_channels) p_times[n_times++] = t.get<long long>();
				nt++;
				//MLOG("%d ", itime);
			}
			//MLOG("val ");
			int nv = 0;
			for (auto &v : d["value"]) {
				string sv = v.get<string>().c_str();
				if (nv < n_val_channels) {
					sv.copy(str_values[n_vals], sv.length());
					str_values[n_vals][sv.length()] = 0;
					nv++;
					n_vals++;
				}
				//MLOG("%s ", v.get<string>().c_str());
			}
			//MLOG("\n");
		}
		DYN(AM_API_AddDataStr(am, pid, sig.c_str(), n_times, p_times, n_vals, str_values));
	}
}

//========================================================================================
void get_json_examples_preds(po::variables_map vm, AlgoMarker *am)
{
	// a small technicality
	((MedialInfraAlgoMarker *)am)->set_sort(0); // getting rid of cases in which multiple data sets on the same day cause differences and fake failed tests.

	vector<string> jsons;
	split_file_to_jsons(vm["in_jsons"].as<string>(), jsons);
	vector<json> all_jsons(jsons.size());

	for (int i = 0; i < jsons.size(); i++) {
		if (i % 1000 == 0) MLOG("so far parsed as objects %d jsons\n", i);
		all_jsons[i] = json::parse(jsons[i]);
	}

	MLOG("Parsed %d jsons\n", all_jsons.size());
	int MAX_VALS = 10000, MAX_VAL_LEN = 200;
	vector<char *>str_vals(MAX_VALS);
	vector<char> buf(MAX_VALS*MAX_VAL_LEN);
	for (int i = 0; i < MAX_VALS; i++) { str_vals[i] = &buf[i*MAX_VAL_LEN]; }
	vector<long long> times(MAX_VALS);
	int pid;
	long long _timestamp;
	char *stypes[] = { "Raw" };
	vector<MedSample> res;
	ofstream dummy_stream;
	int print_msgs = (int)vm.count("print_msgs");
	int n_tested = 0;

	for (int j = 0; j < jsons.size(); j++) {
		//for (auto &js : all_jsons) {
		json &js = all_jsons[j];
		int add_data_rc = add_data_to_am_from_json(am, js, jsons[j].c_str(), &str_vals[0], &times[0], pid, _timestamp);

		AMRequest *req;
		int req_create_rc = DYN(AM_API_CreateRequest("test_request", stypes, 1, &pid, &_timestamp, 1, &req));
		if (req == NULL) MTHROW_AND_ERR("ERROR: Got a NULL request for pid %d time %d rc %d!!\n", pid, (int)_timestamp, req_create_rc);

		// create a response
		AMResponses *resp;
		DYN(AM_API_CreateResponses(&resp));

		if (add_data_rc > 0) {
			//add error message to mark error when adding data
			AMMessages *msgs = resp->get_shared_messages();
			msgs->insert_message(add_data_rc, "AddData Failed");
		}

		// Calculate
		DYN(AM_API_Calculate(am, req, resp));

		int n_resp = DYN(AM_API_GetResponsesNum(resp));
		//MLOG("pid %d time %d n_resp %d\n", s.id, s.time, n_resp);

		// get scores
		if (n_resp == 1) {
			AMResponse *response;
			int resp_rc = DYN(AM_API_GetResponseAtIndex(resp, 0, &response));
			MedSample rs;
			get_response_score_into_sample(response, resp_rc, rs);
			res.push_back(rs);
		}
		else {
			MedSample rs;
			rs.id = pid;
			rs.time = (int)_timestamp;
			rs.prediction.clear();
			rs.prediction.push_back((float)AM_UNDEFINED_VALUE);
			res.push_back(rs);
		}

		if (print_msgs) print_response_msgs(resp, pid, _timestamp, dummy_stream);

		// and now need to dispose responses and request
		DYN(AM_API_DisposeRequest(req));
		DYN(AM_API_DisposeResponses(resp));

		// clearing data in algomarker
		DYN(AM_API_ClearData(am));

		float pred = res.back().prediction[0];

		MLOG("(%d) pid %d time %d score %f\n", n_tested++, pid, _timestamp, pred);

	}

	if (vm["preds_jsons"].as<string>() != "") {
		MedSamples samples;
		samples.import_from_sample_vec(res);
		samples.write_to_file(vm["preds_jsons"].as<string>());
	}

}


//========================================================================================
void json_req_test(po::variables_map &vm, AlgoMarker *am)
{
	// loading data from in_jsons if provided
	if (vm["in_jsons"].as<string>() != "") {
		string in_jsons;
		read_file_into_string(vm["in_jsons"].as<string>(), in_jsons);
		MLOG("read %d characters from input jsons file %s\n", in_jsons.length(), vm["in_jsons"].as<string>().c_str());
		int load_status = DYN(AM_API_AddDataByType(am, 0, DATA_BATCH_JSON_FORMAT, in_jsons.c_str()));
		MLOG("Added data from %s\n", vm["in_jsons"].as<string>().c_str());
		if (load_status != AM_OK_RC)
			MERR("Error code returned from calling AddDataByType: %d\n", load_status);
	}


	// direct call to CalculateByType
	string sjreq;
	read_file_into_string(vm["jreq"].as<string>(), sjreq);
	char *jreq = (char *)(sjreq.c_str());
	char *jresp;
	MLOG("Before Calculate jreq len %d\n", sjreq.length());
	DYN(AM_API_CalculateByType(am, JSON_REQ_JSON_RESP, jreq, &jresp));
	//((MedialInfraAlgoMarker *)am)->CalculateByType(JSON_REQ_JSON_RESP, jreq, &jresp);	
	MLOG("After Calculate jresp len %d\n", strlen(jresp));

	if (vm["jresp"].as<string>() != "") {
		string s = string(jresp);
		write_string(vm["jresp"].as<string>(), s);
	}

	DYN(AM_API_Dispose(jresp));
}


//========================================================================================
// prep dictionaries for AA
//========================================================================================
void prep_dicts(po::variables_map &vm)
{

	ofstream out_json_f;
	out_json_f.open(vm["out_json_dict"].as<string>());
	if (!out_json_f.is_open())
		MTHROW_AND_ERR("Can't open output json dict file %s for writing\n", vm["out_json_dict"].as<string>().c_str());

	// input config file format:

	// SIGNAL	IN/OUT	code
	vector<vector<string>> res;
	read_text_file_cols(vm["dicts_config"].as<string>(), "\t", res);
	unordered_map<string, unordered_set<string>> in_vals, out_vals;
	for (auto &v : res) {
		//MLOG("%s %s %s\n", v[0].c_str(), v[1].c_str(), v[2].c_str());
		if (v[1] == "IN")
			in_vals[v[0]].insert(v[2]);
		else if (v[1] == "OUT")
			out_vals[v[0]].insert(v[2]);
	}


	MedRepository rep;
	if (rep.init(vm["rep"].as<string>()) < 0) MTHROW_AND_ERR("Can't open rep %s\n", vm["rep"].as<string>().c_str());

	json dj = json({});

	if (vm.count("simple_dict") == 0) dj += { "dictionary", json::array() };

	for (auto &sig : in_vals) {
		MLOG("Working on sig %s\n", sig.first.c_str());
		int section_id = rep.dict.section_id(sig.first);
		// preparing helper lists
		unordered_map<int, string> id2out;
		if (out_vals.find(sig.first) != out_vals.end()) {
			for (auto &v : out_vals[sig.first]) {
				if (rep.dict.dicts[section_id].Name2Id.find(v) == rep.dict.dicts[section_id].Name2Id.end())
					MTHROW_AND_ERR("ERROR: Can't find value %s in dictionary for signal %s\n", v.c_str(), sig.first.c_str());
				int id = rep.dict.dicts[section_id].Name2Id[v];
				id2out[id] = v;
				//MLOG("id2out: %s %d %s\n", sig.first.c_str(), id, v.c_str());
			}
		}

		if (vm.count("simple_dict") == 0) {
			json js = json({});
			js += { "signal", sig.first };
			js += { "signal_map", json::array() };

			for (auto &v : sig.second) {
				if (rep.dict.dicts[section_id].Name2Id.find(v) == rep.dict.dicts[section_id].Name2Id.end())
					MTHROW_AND_ERR("ERROR: Can't find value %s in dictionary for signal %s\n", v.c_str(), sig.first.c_str());
				json jdef = json({});
				jdef += { "def", v };
				jdef += { "sets", json::array() };
				int v_id = rep.dict.dicts[section_id].Name2Id[v];
				if (rep.dict.dicts[section_id].Member2Sets[v_id].size() > 0) {
					for (auto s : rep.dict.dicts[section_id].Member2Sets[v_id])
						if (id2out.find(s) != id2out.end())
							jdef["sets"].push_back(id2out[s]);
				}
				js["signal_map"].push_back(jdef);
			}
			dj["dictionary"].push_back(js);
		}
		else {
			//dj += {sig, {}};
			json jdefs = json({});
			for (auto &v : sig.second) {
				if (rep.dict.dicts[section_id].Name2Id.find(v) == rep.dict.dicts[section_id].Name2Id.end())
					MTHROW_AND_ERR("ERROR: Can't find value %s in dictionary for signal %s\n", v.c_str(), sig.first.c_str());
				jdefs += { v, json::array() };
				int v_id = rep.dict.dicts[section_id].Name2Id[v];
				if (rep.dict.dicts[section_id].Member2Sets[v_id].size() > 0) {
					for (auto s : rep.dict.dicts[section_id].Member2Sets[v_id])
						if (id2out.find(s) != id2out.end())
							jdefs[v].push_back(id2out[s]);
				}
			}
			dj += {sig.first, jdefs};
		}
	}

	out_json_f << dj.dump(1) << endl;
	out_json_f.close();
	MLOG("Wrote jsons output file %s\n", vm["out_json_dict"].as<string>().c_str());

}

//========================================================================================
// loaders and helpers
//========================================================================================

//========================================================================================
void load_amlib(po::variables_map &vm)
{
	// upload dynamic library if needed
	if (vm["amlib"].as<string>() != "") {
		load_am(vm["amlib"].as<string>().c_str());
	}

	if (DynAM::initialized())
		MLOG("Dynamic %s library loaded\n", vm["amlib"].as<string>().c_str());
	else
		MLOG("Dynamic library not loaded\n");

}

//========================================================================================
void initialize_algomarker(po::variables_map &vm, AlgoMarker *&test_am)
{
	// Initialize AlgoMarker
	MLOG("Creating AM\n");

	if (DYN(AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am)) != AM_OK_RC)
		MTHROW_AND_ERR("ERROR: Failed creating test algomarker\n");

	if (vm["am_csv"].as<string>() != "")	((MedialInfraAlgoMarker *)test_am)->set_am_matrix(vm["am_csv"].as<string>());

	// Load
	MLOG("Loading AM\n");
	int rc = DYN(AM_API_Load(test_am, vm["amconfig"].as<string>().c_str()));
	if (rc != AM_OK_RC) MTHROW_AND_ERR("ERROR: Failed loading algomarker %s with config file %s ERR_CODE: %d\n", test_am->get_name(), vm["amconfig"].as<string>().c_str(), rc);

	MLOG("Name is %s\n", test_am->get_name());

	// Additional load of dictionaries (if provided)
	if (vm["json_dict"].as<string>() != "") {
		MLOG("Additional Loading AM\n");
		vector<string> dictionaries;
		boost::split(dictionaries, vm["json_dict"].as<string>(), boost::is_any_of(","));
		for (auto fdict : dictionaries) {
			int rc = DYN(AM_API_AdditionalLoad(test_am, LOAD_DICT_FROM_FILE, fdict.c_str()));
			if (rc != AM_OK_RC) {
				MTHROW_AND_ERR("ERROR: Failed additional loading of dict : algomarker %s with dict file %s ERR_CODE: %d\n", test_am->get_name(), fdict.c_str(), rc);
			}
			else
				MLOG("Loaded dictionary %s\n", fdict.c_str());
		}
	}

	if (vm["discovery"].as<string>() != "") {
		ofstream fw(vm["discovery"].as<string>());
		if (!fw.good())
			MTHROW_AND_ERR("Error can't write to %s\n", vm["discovery"].as<string>().c_str());

		char *res;
		DYN(AM_API_Discovery(test_am, &res));


		if (res != NULL) {
			fw << res << endl;
			delete[] res;
		}
		fw.close();
	}
}

//========================================================================================
void create_jreq(po::variables_map &vm)
{
	json jreq = json({});
	if (vm["jreq_defs"].as<string>() != "") {
		string s;
		if (read_file_into_string(vm["jreq_defs"].as<string>(), s) < 0) MTHROW_AND_ERR("could not open file for jreq_defs\n");
		jreq = json::parse(s);
	}

	MedSamples samples;
	samples.read_from_file(vm["samples"].as<string>());

	jreq.push_back({ "requests" , json::array() });
	for (auto &ids : samples.idSamples)
		for (auto &s : ids.samples) {
			json js = json({});
			js.push_back({ "patient_id", to_string(s.id) });
			js.push_back({ "time" , to_string(s.time) });
			jreq["requests"].push_back(js);
		}

	if (vm["jreq_out"].as<string>() != "") {
		string s = jreq.dump(1);
		write_string(vm["jreq_out"].as<string>(), s);
		MLOG("Wrote request json to: %s\n", vm["jreq_out"].as<string>().c_str());
	}
}

//========================================================================================
// MAIN
//========================================================================================

void load_repo(MedPidRepository &rep, MedModel &model, const vector<int> pids,
	const string &rep_path, bool allow_adjustment, map<string, string> &rename_signal) {
	if (allow_adjustment) {
		if (rep.init(rep_path) < 0)
			MTHROW_AND_ERR("Cannot initialize repository from %s\n", rep_path.c_str());
		model.fit_for_repository(rep);
	}

	// Get Required signals
	vector<string> req_signals;
	model.get_required_signal_names(req_signals);

	for (size_t i = 0; i < req_signals.size(); ++i)
		if (rename_signal.find(req_signals[i]) != rename_signal.end())
			req_signals[i] = rename_signal.at(req_signals[i]);

	// Read Repository
	MLOG("Reading Repository from %s\n", rep_path.c_str());
	if (rep.read_all(rep_path, pids, req_signals) != 0)
		MTHROW_AND_ERR("Read repository from %s failed\n", rep_path.c_str());
}

int main(int argc, char *argv[])
{
	int rc = 0;
	po::variables_map vm;

	// Running Parameters
	MLOG("Reading params\n");
	rc = read_run_params(argc, argv, vm);
	assert(rc >= 0);

	if (vm["out_json_dict"].as<string>() != "") {
		prep_dicts(vm);
		return 0;
	}

	if (vm.count("create_jreq") > 0) {
		create_jreq(vm);
		return 0;
	}

	//--------------------------------------------------------------------
	// preparations for tests
	//--------------------------------------------------------------------
	AlgoMarker *test_am;

	load_amlib(vm);
	initialize_algomarker(vm, test_am);

	if (vm.count("json_req_test") > 0) {
		json_req_test(vm, test_am);
		return 0;
	}


	// from here and on: tests that involve loading a model, and comparing the 
	// results from a repository + model with the repository + algomarker. (score compare test).
	// read model file
	MedModel model;
	if (model.read_from_file(vm["model"].as<string>()) < 0) {
		MERR("FAILED reading model file %s\n", vm["model"].as<string>().c_str());
		return -1;
	}

	// init ignore_sig
	vector<string> ignore_sig;
	if (vm["ignore_sig"].as<string>() != "")
		split(ignore_sig, vm["ignore_sig"].as<string>(), boost::is_any_of(","));

	vector<string> sigs;
	model.get_required_signal_names(sigs);



	// test that AM repo. contains all the required signal
	unordered_set<string> am_signals;
	((MedialInfraAlgoMarker *)test_am)->get_am_rep_signals(am_signals);
	for (auto &sig : sigs)
	{
		if (am_signals.count(sig) == 0)
		{
			MTHROW_AND_ERR("AlgoMarker's repository doesn't contain sig [%s]", sig.c_str());
		}
	}

	if (vm["in_jsons"].as<string>() != "") {
		get_json_examples_preds(vm, test_am);
		return 0;
	}


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
	map<string, string> rename_signal;
	//Read for vm syntax:
	if (!vm["rename_signal"].as<string>().empty()) {
		//read and prase file:
		ifstream file_reader(vm["rename_signal"].as<string>());
		if (!file_reader.good())
			MTHROW_AND_ERR("Error can't read file %s\n", vm["rename_signal"].as<string>().c_str());
		MLOG("Reading file %s for regex filter\n", vm["rename_signal"].as<string>().c_str());
		string line;
		while (getline(file_reader, line)) {
			boost::trim(line);
			if (line.empty() || line[0] == '#')
				continue;
			vector<string> tokens;
			boost::split(tokens, line, boost::is_any_of("\t"));
			if (tokens.size() != 2)
				MTHROW_AND_ERR("Error bad file format %s. expecting 2 tokens\n",
					vm["rename_signal"].as<string>().c_str());
			rename_signal[tokens[1]] = tokens[0];
		}
		file_reader.close();
	}

	//To support rename:
	//model.load_repository(vm["rep"].as<string>(), pids, rep, vm["allow_rep_adjustment"].as<bool>());
	load_repo(rep, model, pids, vm["rep"].as<string>(), vm["allow_rep_adjustment"].as<bool>(), rename_signal);

	//commit rename after reading repo:
	for (const auto &it : rename_signal)
	{
		if (rep.sigs.sid(it.first) < 0) {
			int source_sid = rep.sigs.sid(it.second);
			int currect_section = rep.dict.section_id(it.second);
			if (source_sid < 0)
				MTHROW_AND_ERR("Unknown signal %s\n", it.second.c_str());
			MWARN("Adding signal %s\n", it.first.c_str());
			rep.sigs.Name2Sid[it.first] = source_sid;
			rep.sigs.Sid2Name[source_sid] = it.first;
			rep.sigs.signals_names.push_back(it.first);

			if (currect_section > 0)
				rep.dict.connect_to_section(it.first, currect_section);
		}
	}

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

	if (vm["direct_csv"].as<string>() != "")
		model.write_feature_matrix(vm["direct_csv"].as<string>());


	// printing
	vector<MedSample> res1;
	samples.export_to_sample_vec(res1);
	for (int i = 0; i < min(50, (int)res1.size()); i++) {
		MLOG("#Res1 :: pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0]);
	}

	// res1 now is the gold standard we compare to : was calculated directly using the infrastructure

	//===============================================================================
	// TEST1: testing internal in_mem in a repository
	//===============================================================================


	MLOG("Algomarker %s was loaded with config file %s\n", test_am->get_name(), test_am->get_config());
	vector<MedSample> res2;

	int print_msgs = (vm.count("print_msgs")) ? 1 : 0;
	string msgs_file = (vm["msgs_file"].as<string>());
	ofstream msgs_stream;
	if (msgs_file != "") {
		msgs_stream.open(msgs_file);
		msgs_stream << "msg_type\tpid\tdate\ti\tj\tk\tcode\tmsg_text" << endl;
	}

	if (vm.count("single"))
		get_preds_from_algomarker_single(vm, test_am, vm["rep"].as<string>(), rep, model, samples2, pids, sigs, res2, res1, msgs_stream, ignore_sig, rename_signal);
	else
		get_preds_from_algomarker(test_am, vm["rep"].as<string>(), rep, model, samples2, pids, sigs, res2, print_msgs, msgs_stream, ignore_sig);
	for (int i = 0; i < min(50, (int)res1.size()); i++) {
		MLOG("#Res1 :: pid %d time %d pred %f #Res2 pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0], res2[i].id, res2[i].time, res2[i].prediction[0]);
	}

	AM_API_DisposeAlgoMarker(test_am);

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

	if (msgs_file != "")
		msgs_stream.close();
	if (vm["am_res_file"].as<string>() != "")
		save_sample_vec(res2, vm["am_res_file"].as<string>());
}

//
// keep command line:
//
// typical test:
// Linux/Release/DllAPITester --model /nas1/Work/Users/Avi/Diabetes/order/pre2d/runs/partial/pre2d_partial_S6.model --samples test_100k.samples --amconfig /nas1/Work/Users/Avi/AlgoMarkers/pre2d/pre2d.amconfig
