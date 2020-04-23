#include <AlgoMarker/CommonTestingTools/CommonTestingTools.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL


string CommonTestingTools::precision_float_to_string(float val) {
	stringstream ss;
	ss << std::setprecision(10) << val;
	return ss.str();
}

//Expand string with embedded Environment variables in it
string CommonTestingTools::expandEnvVars(const string &str) {
	string ret = "";
#ifdef __linux__ 
	wordexp_t p;
	char** w;
	wordexp(str.c_str(), &p, 0);
	w = p.we_wordv;
	for (size_t i = 0; i < p.we_wordc; i++) ret += w[i];
	wordfree(&p);
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


char** CommonTestingTools::charpp_adaptor::get_charpp() {
	if (this->size() == 0)
		return nullptr;
	size_t charpp_arr_sz = this->size() * sizeof(char*);
	size_t charpp_buf_sz = 0;
	for (auto& str : *this) {
		charpp_buf_sz += str.size() + 1;
	}
	charpp_buf_sz *= sizeof(char);

	charpp_arr = (char**)realloc(charpp_arr, charpp_arr_sz);
	charpp_buf = (char*)realloc(charpp_buf, charpp_buf_sz);

	char** charpp_arr_i = charpp_arr;
	char* charpp_buf_i = charpp_buf;

	for (auto& str : *this)
	{
		*charpp_arr_i = charpp_buf_i;
		charpp_arr_i++;
		for (int i = 0; i < str.size(); ++i) {
			*charpp_buf_i = str[i];
			charpp_buf_i++;
		}
		*charpp_buf_i = '\0';
		charpp_buf_i++;
	}
	return charpp_arr;
}

json CommonTestingTools::read_json_array_next_chunk(ifstream& infile, bool& in_array) {
	char prev_c = '\0';
	char c = '\0';
	bool in_string = false;
	string ret_str = "";
	int block_depth = 0;
	while (infile.get(c)) {
		switch (c) {
		case '"':
			if (!in_string)
				in_string = true;
			else if (prev_c != '\\')
				in_string = false;
			break;
		case '{':
			if (!in_array)
				throw runtime_error("File should be a JSON array containing objects");
			if (!in_string)
				block_depth++;
			break;
		case '}':
			if (!in_array)
				throw runtime_error("Did not expect a '}'");
			if (!in_string)
				block_depth--;
			if (block_depth < 0)
				throw runtime_error("Did not expect a '}'");
			break;
		}
		if (c == '[' && !in_array) {
			in_array = true;
			continue;
		}
		if ((c == ']' || c == ',') && in_array && block_depth == 0)
			break;
		ret_str += c;
		prev_c = c;
	}
	json ret;
	if (std::all_of(ret_str.begin(), ret_str.end(), [](char c) { return c == ' ' || c == '\n' || c == '\t' || c == '\r'; }))
		return ret;
	try {
		ret = json::parse(ret_str);
	}
	catch (...) {
		MERR("Error parsing chunk: \n'%s'\n", ret_str.c_str());
	}
	return ret;
}

json CommonTestingTools::json_AddData(const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values, int n_time_channels, int n_val_channels) {
	//MLOG("json_AddData sig %s nT(%d) %d nV(%d) %d\n", signalName, n_time_channels, TimeStamps_len, n_val_channels, Values_len);
	json json_sig = json({ { "code", signalName },{ "data", json::array() } });
	if (units_tbl.count(signalName) != 0)
		json_sig["unit"] = units_tbl.at(signalName);
	int nelem = 0;
	if (TimeStamps_len != 0)
		nelem = TimeStamps_len / n_time_channels;
	else nelem = Values_len / n_val_channels;
	for (int i = 0; i < nelem; i++) {
		json json_sig_data_item = json({ { "timestamp" , json::array() },{ "value" , json::array() } });

		for (int j = 0; j < n_time_channels; j++) {
			json_sig_data_item["timestamp"].push_back(*TimeStamps);
			TimeStamps++;
		}

		for (int j = 0; j < n_val_channels; j++) {
			json_sig_data_item["value"].push_back(*Values);
			Values++;
		}

		json_sig["data"].push_back(json_sig_data_item);
	}
	return json_sig;
}

json CommonTestingTools::json_AddDataStr(const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values, int n_time_channels, int n_val_channels) {
	//MLOG("json_AddData sig %s nT(%d) %d nV(%d) %d\n", signalName, n_time_channels, TimeStamps_len, n_val_channels, Values_len);
	json json_sig = json({ { "code", signalName },{ "data", json::array() } });
	if (units_tbl.count(signalName) != 0)
		json_sig["unit"] = units_tbl.at(signalName);
	int nelem = 0;
	if (TimeStamps_len != 0)
		nelem = TimeStamps_len / n_time_channels;
	else nelem = Values_len / n_val_channels;
	for (int i = 0; i < nelem; i++) {
		json json_sig_data_item = json({ { "timestamp" , json::array() },{ "value" , json::array() } });

		for (int j = 0; j < n_time_channels; j++) {
			json_sig_data_item["timestamp"].push_back(*TimeStamps);
			TimeStamps++;
		}

		for (int j = 0; j < n_val_channels; j++) {
			json_sig_data_item["value"].push_back(*Values);
			Values++;
		}

		json_sig["data"].push_back(json_sig_data_item);
	}
	return json_sig;
}

int CommonTestingTools::get_preds_from_algomarker(AlgoMarker *am, vector<MedSample> &res, bool print_msgs, DataLoader& d, bool force_add_data, ofstream& msgs_stream, vector<string> ignore_sig, bool extended_score)
{
	DynAM::AM_API_ClearData(am);

	MLOG("=====> now running get_preds_from_algomarker()\n");

	MLOG("Going over %d pids\n", d.pids.size());
	d.get_sig_dict_cached();
	int pid_cnt = 0;
	for (auto pid : d.pids) {
		json json_req;
		d.am_add_data(am, pid, INT_MAX, force_add_data, ignore_sig, json_req);
		pid_cnt++;
		if (pid_cnt % 1000 == 0)
			MLOG("Loaded %d pids...\n", pid_cnt);
	}

	//ASK_AVI: Is this needed?
	//((MedialInfraAlgoMarker *)am)->set_sort(0); // getting rid of cases in which multiple data sets on the same day cause differences and fake failed tests.

	MLOG("After AddData for all batch\n");
	// finish rep loading 
	char *stypes[] = { "Raw" };
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
	//MLOG("Before Calculate\n");
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
	for (int i = 0; i<n_resp; i++) {
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
			s.time = (int)(ts / 10000);
		else
			s.time = ts;
		if (resp_rc == AM_OK_RC && n_scores > 0) {
			float _scr=(float)AM_UNDEFINED_VALUE;
			char *_ext_scr = nullptr;
			resp_rc = DynAM::AM_API_GetResponseScoreByIndex(response, 0, &_scr, &_scr_type);
			if (extended_score) {
				resp_rc = DynAM::AM_API_GetResponseExtendedScoreByIndex(response, 0, &_ext_scr, &_scr_type);
				s.str_attributes["extended_score"] = string(_ext_scr);
			}
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
		for (int i = 0; i<n_msgs; i++) {
			if (msgs_stream.is_open())
				msgs_stream << "SharedMessages\t" << 0 << "\t" << 0 << "\t" << i << "\t" << 0 << "\t" << 0 << "\t" << msg_codes[i] << "\t\"" << msgs_errs[i] << "\"" << endl;
			else
				MLOG("Shared Message %d : code %d : err: %s\n", n_msgs, msg_codes[i], msgs_errs[i]);
		}

		n_resp = DynAM::AM_API_GetResponsesNum(resp);
		for (int i = 0; i<n_resp; i++) {
			AMResponse *r;
			DynAM::AM_API_GetResponseAtIndex(resp, i, &r);
			int n_scores;
			DynAM::AM_API_GetResponseScoresNum(r, &n_scores);

			DynAM::AM_API_GetResponseMessages(r, &n_msgs, &msg_codes, &msgs_errs);
			for (int k = 0; k<n_msgs; k++) {
				if (msgs_stream.is_open())
					msgs_stream << "ResponseMessages\t" << 0 << "\t" << 0 << "\t" << i << "\t0\t" << k << "\t" << msg_codes[k] << "\t\"" << msgs_errs[k] << "\"" << endl;
				else
					MLOG("Response %d : Message %d : code %d : err: %s\n", i, k, msg_codes[k], msgs_errs[k]);
			}

			for (int j = 0; j<n_scores; j++) {
				DynAM::AM_API_GetScoreMessages(r, j, &n_msgs, &msg_codes, &msgs_errs);
				for (int k = 0; k<n_msgs; k++) {
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


//=================================================================================================================
// same test, but running each point in a single mode, rather than batch on whole.
//=================================================================================================================
int CommonTestingTools::get_preds_from_algomarker_single(AlgoMarker *am, vector<MedSample> &res, bool print_msgs, DataLoader& d, bool force_add_data, ofstream& msgs_stream, vector<string> ignore_sig, ofstream& json_reqfile_stream, bool extended_score)
{

	DynAM::AM_API_ClearData(am);

	MLOG("=====> now running get_preds_from_algomarker_single()\n");
	MLOG("Going over %d samples (CommonTestingTools)\n", d.samples.nSamples());
	int n_tested = 0;

	MedTimer timer;
	d.get_sig_dict_cached();
	timer.start();

	bool first_json_req = true;

	json json_resp_byid;

	int counter = 0;
	for (auto &id : d.samples.idSamples) {
		for (auto &s : id.samples) {
			MLOG("===> running sample %d : pid %d time %d\n", counter++, s.id, s.time);
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
			char *stypes[] = { "Raw" };
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
						rs.time = (int)(ts / 10000);
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
				for (int i = 0; i<n_msgs; i++) {

					if (msgs_stream.is_open())
						msgs_stream << "SharedMessages\t" << s.id << "\t" << s.time << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << msg_codes[i] << "\t\"" << msgs_errs[i] << "\"" << endl;
					else
						MLOG("pid %d time %d Shared Message %d : code %d : err: %s\n", s.id, s.time, n_msgs, msg_codes[i], msgs_errs[i]);
				}

				n_resp = DynAM::AM_API_GetResponsesNum(resp);
				for (int i = 0; i<n_resp; i++) {
					AMResponse *r;
					DynAM::AM_API_GetResponseAtIndex(resp, i, &r);
					int n_scores;
					DynAM::AM_API_GetResponseScoresNum(r, &n_scores);

					DynAM::AM_API_GetResponseMessages(r, &n_msgs, &msg_codes, &msgs_errs);
					for (int k = 0; k<n_msgs; k++) {
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

					for (int j = 0; j<n_scores; j++) {
						DynAM::AM_API_GetScoreMessages(r, j, &n_msgs, &msg_codes, &msgs_errs);
						for (int k = 0; k<n_msgs; k++) {
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
				MLOG("==== TIME MEASURE ====> Tested %d samples : time %f sec\n", n_tested, dt);
				dt = (double)n_tested / dt;
				MLOG("==== TIME MEASURE ====> %f samples/sec\n", dt);
			}
		}
	}
	if (json_reqfile_stream.is_open()) {
		json_reqfile_stream << "]";
	}

	MLOG("Finished getting preds from algomarker in a single manner\n");
	return 0;
}

void CommonTestingTools::save_sample_vec(vector<MedSample> sample_vec, const string& fname) {
	MedSamples s;
	s.import_from_sample_vec(sample_vec);
	s.write_to_file(fname, 4);
}


