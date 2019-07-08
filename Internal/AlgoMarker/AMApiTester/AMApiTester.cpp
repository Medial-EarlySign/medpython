//
// Test Program to the Dll API
//
// General Plan :
// 
// Compare same data/model/points prediction using the infrastructure directly AND using the DLL.
//

#define AM_DLL_IMPORT

//#include <AlgoMarker/AlgoMarker/AlgoMarker.h>
#include <AlgoMarker/DynAMWrapper/DynAMWrapper.h>

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

#ifdef __linux__ 
#include <wordexp.h>
#define DEFAULT_AM_LOCATION "${MR_ROOT}/Libs/Internal/AlgoMarker/Linux/Release/libdyn_AlgoMarker.so"
#elif _WIN32
#include "windows.h" 
#define DEFAULT_AM_LOCATION "%MR_ROOT%\\Libs\\Internal\\AlgoMarker\\x64\\ReleaseDLL\\AlgoMarker.dll"
#endif

#include <climits>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

string expandEnvVars(const string &str) {
  string ret = "";
#ifdef __linux__ 
  wordexp_t p;
  char** w;
  wordexp(str.c_str(), &p, 0 );
  w = p.we_wordv;
  for (size_t i=0; i<p.we_wordc;i++ ) ret+= w[i];
  wordfree( &p );
#elif _WIN32
  DWORD max_str_len = 4 * 1024;
  auto buf = new char[max_str_len];
  DWORD req_len = ExpandEnvironmentStrings(str.c_str(), buf, max_str_len);
  if (req_len > max_str_len) {
	  delete buf;
	  buf = new char[req_len];
	  req_len = ExpandEnvironmentStrings(str.c_str(), buf, req_len);
  }
  if (req_len > 0)
	  ret = buf;
  delete buf;
#endif
  return ret;
}

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
			("samples", po::value<string>()->default_value(""), "MedSamples file to use")
			("model", po::value<string>()->default_value(""), "model file to use")
			("amconfig" , po::value<string>()->default_value(""), "AlgoMarker configuration file")
			("single", "Run test in single mode, instead of the default batch")
			("print_msgs", "Print algomarker messages when testing batches or single (direct test always prints them)")
			("egfr_test", "Test simple egfr algomarker")
			("signalsum_test", "Test signalsum algomarker")
			("generate_data", "Generate a unified repository data file for all the signals a model needs (required options: rep,samples,model)")
			("generate_data_outfile", po::value<string>()->default_value(""), "file to output the Generated unified signal file")
			("generate_data_cat_prefix", po::value<string>()->default_value(""), "If provided, prefer to convert a catogorial channel to a name/setname with given prefix")
			("generate_data_force_cat_prefix", "Ignore signals categories which do not conform to generate_data_cat_prefix")
			("apply", "Apply a model using Medial API, given --model, --rep, --apply_repdata, --samples, --apply_outfile, will write scores to output file")
			("apply_repdata", po::value<string>()->default_value(""), "Unified signal data to be used by apply action")
			("apply_dates_to_score", po::value<string>()->default_value(""), "File containing a list of tab seperated pid and date to score to beused instead of scores for performing apply")
			("apply_amconfig", po::value<string>()->default_value(""), "Same as --apply but will use the AlgoMarker API and given amconfig")
			("apply_outfile", po::value<string>()->default_value(""), "Output file to save scores from apply")
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

class DataLoader{
public:
	MedModel model;
	MedSamples samples;
	MedPidRepository rep;
	vector<int> pids;
	vector<string> sigs;
    map<int, MedIdSamples* > pid2samples;

    void load(const string& rep_fname, const string& model_fname, const string& samples_fname="",bool read_signals=true) {
		// read model file
	    if (model.read_from_file(model_fname) < 0) {
		    MERR("FAILED reading model file %s\n", model_fname.c_str());
            throw runtime_error(string("FAILED reading model file ")+model_fname);
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
				throw runtime_error(string("FAILED MedRepository::init(")+rep_fname+"\")");
			}
		}
		for (auto &id : samples.idSamples)
			pid2samples[id.id] = &id;		
    }

	void export_required_data(const string& fname, const string& cat_prefix, bool force_cat_prefix) {
		ofstream outfile(fname, ios::binary | ios::out);
		
		MLOG("(II) Preparing dictinaries to export\n", fname.c_str());

		map<string, vector<map<int, string> > > sig_dict;
		for (auto& sig : sigs) {
			vector<map<int, string > > chan_dict;
			int section_id = rep.dict.section_id(sig);
			int sid = rep.sigs.Name2Sid[sig];
			int n_vchan = rep.sigs.Sid2Info[sid].n_val_channels;
			for (int vchan = 0; vchan < n_vchan; ++vchan) {
				if (rep.sigs.is_categorical_channel(sig, vchan)) {
					map<int, string> new_dict;
					const auto& Id2Name = rep.dict.dict(section_id)->Id2Name;
					const auto& Member2Sets = rep.dict.dict(section_id)->Member2Sets;
					for (const auto& entry : Id2Name) {
						if (boost::starts_with(entry.second, cat_prefix)) {
							new_dict[entry.first] = entry.second;
							continue;
						}
						string new_ent = entry.second;
						if(Member2Sets.count(entry.first) != 0)
						for (const auto& setid : Member2Sets.at(entry.first)) {
							if (Id2Name.count(setid) != 0 && boost::starts_with(Id2Name.at(setid), cat_prefix)) {
								if(!boost::starts_with(new_ent, cat_prefix) || new_ent.length() > Id2Name.at(setid).length())
								new_ent = Id2Name.at(setid);
							}
						}
						if(!force_cat_prefix || boost::starts_with(new_ent, cat_prefix))
							new_dict[entry.first] = new_ent;
					}
					
					chan_dict.push_back(new_dict);
					auto& dict = new_dict; //rep.dict.dict(section_id)->Id2Name;
					ofstream f;
					f.open("/nas1/Work/Users/Shlomi/apply-program/generated/dict.orig.tsv");
					for (const auto& entry : dict) {
						f << entry.first << '\t' << entry.second << '\n';
					}
					f.close();
				}
				else chan_dict.push_back(map<int, string>());
			}
			sig_dict[sig] = chan_dict;
		}
		
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
						if(sig_dict.at(sig)[vchan].size() == 0)
							outss << '\t' << setprecision(10) << usv.Val(i, vchan);
						else {
							if (sig_dict.at(sig)[vchan].count((int)(usv.Val(i, vchan))) != 0) {
								outss << '\t' << sig_dict.at(sig)[vchan].at((int)(usv.Val(i, vchan)));
							}
							else{
								ignore_line = true;
							}
						}
					}
					if(!ignore_line)
						outfile << outss.str() << '\n';
				}
			}
		}
		outfile.close();
	}

	void import_required_data(const string& fname) {
		ifstream infile(fname, ios::binary | ios::in);

		MLOG("(II)   Preparing signal dictionaries\n");

		map<string, vector<map<string, int >* > > sig_dict;
		for (auto& sig : sigs) {
			MLOG("(II)   Preparing signal dictionary for signal '%s'\n", sig.c_str());
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
				for (int vchan = 0 ; vchan < n_vchan; ++vchan) {
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
						catch(...){
							MERR("Error converting sig %s, chan %d, '%s' back to code\n",sig.c_str(), vchan, fields[fields_i-1].c_str());
							/*
							auto& dict = *(sig_dict.at(sig)[vchan]);
							ofstream f;
							f.open("/nas1/Work/Users/Shlomi/apply-program/generated/dict.tsv");
							for (const auto& entry : dict) {
								f << entry.first << '\t' << entry.second << '\n';
							}
							f.close();
							*/
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

	int load_samples_from_dates_to_score(const string& fname)
	{
		// read scores file
		vector<vector<string>> raw_scores;
		if (read_text_file_cols(fname, " \t", raw_scores) < 0) {
			MERR("Could not read scores file %s\n", fname.c_str());
			return -1;
		}
		MLOG("Read %d lines from scores file %s\n", raw_scores.size(), fname.c_str());

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

    void am_add_data(AlgoMarker *am, int pid, int max_date=INT_MAX){
        UniversalSigVec usv;

	    int max_vals = 100000;
    	vector<long long> times(max_vals);
    	vector<float> vals(max_vals);

	    for (auto &sig : sigs) {
		    rep.uget(pid, sig, usv);
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
						    if (usv.Time(i, j) <= max_date) nelem_before = i+1;
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
   
			    if ((i_val > 0) || (i_time > 0))
				    AM_API_AddData(am, pid, sig.c_str(), i_time, p_times, i_val, p_vals);
		    }
	    }
    }

};

//=================================================================================================================
int get_preds_from_algomarker(AlgoMarker *am, vector<MedSample> &res, bool print_msgs, DataLoader& d)
{
	AM_API_ClearData(am);

	MLOG("=====> now running get_preds_from_algomarker()\n");
    
	MLOG("Going over %d pids\n", d.pids.size());
	for (auto pid : d.pids) {
        d.am_add_data(am, pid);
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
int get_preds_from_algomarker_single(AlgoMarker *am, vector<MedSample> &res, bool print_msgs, DataLoader& d)
{

	AM_API_ClearData(am);

	MLOG("=====> now running get_preds_from_algomarker_single()\n");
	MLOG("Going over %d samples\n", d.samples.nSamples());
	int n_tested = 0;

	MedTimer timer;
	timer.start();
	for (auto &id : d.samples.idSamples){
		for (auto &s : id.samples) {

			// adding all data 
			d.am_add_data(am, s.id, s.time);

			// At this point we can send to the algomarker and ask for a score

			// a small technicality
			// ASK_AVI
			//((MedialInfraAlgoMarker *)am)->set_sort(0); // getting rid of cases in which multiple data sets on the same day cause differences and fake failed tests.

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

	MLOG("Finished getting preds from algomarker in a single manner\n");
	return 0;
}


#endif

//--------------------------------------------------------------------------------------------------------------------------------
int simple_egfr_test()
{
	// init AM
	AlgoMarker *test_am;

	if (AM_API_Create((int)AM_TYPE_SIMPLE_EXAMPLE_EGFR, &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating test algomarker\n");
		return -1;
	}


	int load_rc;
	if ((load_rc = AM_API_Load(test_am, "AUTO") != AM_OK_RC)) {
		MERR("ERROR: Failed loading algomarker , rc: %d\n", load_rc);
		return -1;
	}
	MLOG("Algomarker was loaded\n");


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

int signalsum_test()
{
	// init AM
	AlgoMarker *test_am;

	if (AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating SignalSum algomarker\n");
		return -1;
	}


	int load_rc;
	if ((load_rc = AM_API_Load(test_am, "SignalSum_LOG|") != AM_OK_RC)) {
		MERR("ERROR: Failed loading algomarker , rc: %d\n", load_rc);
		return -1;
	}
	MLOG("Algomarker was loaded\n");

	/*
	// Load Data
	vector<long long> times = { 20160101 };
	vector<float> vals = { 2.0 };
	AM_API_AddData(test_am, 1, "Creatinine", (int)times.size(), &times[0], (int)vals.size(), &vals[0]);
	vals = { 55 };
	AM_API_AddData(test_am, 1, "Age", 0, NULL, (int)vals.size(), &vals[0]);
	vals = { 1 };
	AM_API_AddData(test_am, 1, "GENDER", 0, NULL, (int)vals.size(), &vals[0]);

	// Calculate
	char *stypes[] = { "Raw" };
	vector<int> _pids = { 1 };
	vector<long long> _timestamps = { 20160101 };
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
	for (int i = 0; i<n_shared_msgs; i++) {
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
	for (int i = 0; i<n_resp; i++) {
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
	for (int i = 0; i<n_msgs; i++) {
		MLOG("Shared Message %d : code %d : err: %s\n", n_msgs, msg_codes[i], msgs_errs[i]);
	}


	// Dispose
	AM_API_DisposeRequest(req);
	AM_API_DisposeResponses(resp);
	*/
	AM_API_DisposeAlgoMarker(test_am);

	MLOG("Finished signalsum_test()\n");

	return 0;
}

int generate_data(const string& rep_file, const string& samples_file, const string& model_file, const string& output_file, const string& cat_prefix, bool force_cat_prefix) {
	DataLoader l;
	l.load(rep_file, model_file, samples_file);
	l.export_required_data(output_file, cat_prefix, force_cat_prefix);
	return 0;
}

vector<MedSample> apply_am_api(const string& amconfig, DataLoader& d, bool print_msgs, bool single){
	vector<MedSample> res2;
	AlgoMarker *test_am;

	if (AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am) != AM_OK_RC) {
		MERR("ERROR: Failed creating test algomarker\n");
		throw runtime_error("ERROR: Failed creating test algomarker\n");
	}

    int rc=0;
	if ((rc = AM_API_Load(test_am, amconfig.c_str())) != AM_OK_RC) {
		MERR("ERROR: Failed loading algomarker with config file %s ERR_CODE: %d\n", amconfig.c_str(), rc);
		throw runtime_error(string("ERROR: Failed loading algomarker with config file ")+amconfig+" ERR_CODE: "+to_string(rc));
	}

	if (single)
		get_preds_from_algomarker_single(test_am, res2, print_msgs, d);
	else
		get_preds_from_algomarker(test_am, res2, print_msgs, d);

    return res2;
}

vector<MedSample> apply_med_api(MedPidRepository& rep, MedModel& model, MedSamples& samples){
	// apply model (+ print top 50 scores)
	model.apply(rep, samples);

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
    s.write_to_file(fname);
}

int apply_data(const string& repdata_file, const string& mock_rep_file, const string& scores_file, bool score_format_is_samples, const string& model_file, const string& scores_output_file, const string& amconfig_file) {
	
	DataLoader l;
	MLOG("(II) Starting apply with:\n(II)   repdata_file='%s'\n(II)   mock_rep_file='%s'\n(II)   scores_file='%s' %s\n(II)   model_file='%s'\n(II)   scores_output_file='%s'\n(II)   amconfig_file='%s'\n"
		, repdata_file.c_str(), mock_rep_file.c_str(), scores_file.c_str(), score_format_is_samples ? "(samples format)" : "", model_file.c_str(), scores_output_file.c_str(), amconfig_file.c_str());
	MLOG("(II) Loading mock repo, model and date for scoring\n");

	if (!score_format_is_samples) {
		l.load_samples_from_dates_to_score(scores_file);
		l.load(mock_rep_file, model_file,"",false);
		MLOG("\n(II) Loading tab seperated pid+dates for scoring from %s\n", scores_file.c_str());
	}
	else { 
		MLOG("\n(II) Loading dates for scoring from samples file %s\n", scores_file.c_str());
		l.load(mock_rep_file, model_file, scores_file,false); 
	}
	
	//l.rep.switch_to_in_mem_mode();
	MLOG("(II) Importing data from '%s'\n", repdata_file.c_str());
	l.import_required_data(repdata_file);

	if (amconfig_file == "") {
		MLOG("(II) Starting apply using Medial API\n");
		auto ret = apply_med_api(l.rep, l.model, l.samples);
		MLOG("(II) Saving results to %s\n", scores_output_file.c_str());
		save_sample_vec(ret, scores_output_file);
	}
	else {
		MLOG("(II) Starting apply using Algomarker API\n");
		auto ret = apply_am_api(amconfig_file, l, false, false);
		MLOG("(II) Saving results to %s\n", scores_output_file.c_str());
		save_sample_vec(ret, scores_output_file);
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

	if (vm.count("generate_data")) {
		if (vm["rep"].as<string>() == "" || vm["samples"].as<string>() == "" || vm["model"].as<string>() == "" || vm["generate_data_outfile"].as<string>() == "")
		{
			std::cerr << "Missing argument, Please specify --rep, --samples, --model, --generate_data_outfile.\n";
			return -1;
		}
		return generate_data(
			vm["rep"].as<string>(), 
			vm["samples"].as<string>(), 
			vm["model"].as<string>(), 
			vm["generate_data_outfile"].as<string>(), 
			vm["generate_data_cat_prefix"].as<string>(), 
			vm.count("generate_data_force_cat_prefix")!=0);
	}
	if (vm.count("apply") || (vm.count("apply_amconfig") && vm["apply_amconfig"].as<string>() != "")) {
		if (vm["rep"].as<string>() == "" || 
			(vm["samples"].as<string>() == "" && vm["apply_dates_to_score"].as<string>() =="" ) ||
			vm["model"].as<string>() == "" || 
			vm["apply_outfile"].as<string>() == "" ||
			vm["apply_repdata"].as<string>() == "" ) 
		{
			MERR("Missing arguments, Please specify --rep, --model, --apply_outfile, --apply_repdata, --samples (or --apply_dates_to_score).\n");
			return -1;
		}
		string scores_file = vm["samples"].as<string>();
		bool score_to_date_format_is_samples = true;
		if (vm["apply_dates_to_score"].as<string>() != "") {
			scores_file = vm["apply_dates_to_score"].as<string>();
			score_to_date_format_is_samples = false;
		}
		return apply_data(vm["apply_repdata"].as<string>(), vm["rep"].as<string>(), scores_file, score_to_date_format_is_samples, vm["model"].as<string>(), vm["apply_outfile"].as<string>(), vm["apply_amconfig"].as<string>());
	}
    
    load_am(vm["amfile"].as<string>().c_str());

	if (vm.count("egfr_test"))
		return simple_egfr_test();

	if (vm.count("signalsum_test"))
		return signalsum_test();

    DataLoader d;
	vector<MedSample> res1;
	MedSamples samples,samples2;

    try {
        d.load(vm["rep"].as<string>(),
            vm["model"].as<string>(),
            vm["samples"].as<string>());

        samples2 = samples = d.samples;
    
        res1 = apply_med_api(d.rep, d.model,samples);
    }catch(runtime_error e){
      return -1;
    }
 
	for (int i=0; i<min(50, (int)res1.size()); i++) {
		MLOG("#Res1 :: pid %d time %d pred %f\n", res1[i].id, res1[i].time, res1[i].prediction[0]);
	}

    //fake failed
    //res1[3].prediction[0] = 0.1;
   
	//===============================================================================
	// TEST1: testing internal in_mem in a repository
	//===============================================================================
    vector<MedSample> res2;
    try{
        res2 = apply_am_api(vm["amconfig"].as<string>(), d, (vm.count("print_msgs")!=0),(vm.count("single")!=0) );
    }catch(runtime_error e){
      return -1;
    }

    if(vm["med_res_file"].as<string>()!="")
        save_sample_vec(res1, vm["med_res_file"].as<string>());
    if(vm["am_res_file"].as<string>()!="")
        save_sample_vec(res2, vm["am_res_file"].as<string>());

    compare_results(res1, res2);

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