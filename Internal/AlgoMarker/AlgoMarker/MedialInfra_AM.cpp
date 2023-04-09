#include "AlgoMarker.h"

#include <Logger/Logger/Logger.h>
#include <MedTime/MedTime/MedTime.h>
#include <MedUtils/MedUtils/MedUtils.h>
#include <json/json.hpp>
#include "AlgoMarkerErr.h"

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL


#define PREDICTION_SOURCE_UNKNOWN	0
#define PREDICTION_SOURCE_ATTRIBUTE	1
#define PREDICTION_SOURCE_ATTRIBUTE_AS_JSON	2
#define PREDICTION_SOURCE_JSON	3
#define PREDICTION_SOURCE_PREDICTIONS	4
//#define AM_TIMING_LOGS

class json_req_export {
public:
	string field;
	int type = PREDICTION_SOURCE_UNKNOWN;
	int pred_channel = -1; // relevant only if type is PREDICTION_SOURCE_PREDICTIONS
};

class json_req_info {
public:
	int sample_pid = -1;
	long long sample_time = -1;
	int load_data = 0;
	unordered_map<string, json_req_export> exports;

	int conv_time = -1; // this one is calculated
	int sanity_test_rc = 0; // calculated, keeping eligibility testing result
	int sanity_caught_err = 0;
	vector<InputSanityTesterResult> sanity_res;
	MedSample *res = NULL;
};

// local helper functions (these are CalculateByType helpers)
void add_to_json_array(json &js, const string &key, const string &s_add);
void add_to_json_array(nlohmann::ordered_json &js, const string &key, const string &s_add);
void json_to_char_ptr(json &js, char **jarr);
void json_to_char_ptr(nlohmann::ordered_json &js, char **jarr);
bool json_verify_key(json &js, const string &key, int verify_val_flag, const string &val);
bool json_verify_key(nlohmann::ordered_json &js, const string &key, int verify_val_flag, const string &val);
int json_parse_request(json &jreq, json_req_info &defaults, json_req_info &req_i);

//===========================================================================================================
//===========================================================================================================
// MedialInfraAlgoMarker Implementations ::
// Follows is an implementation of an AlgoMarker , which basically means filling in the:
// Load , Unload, ClearData, AddData and Calculate APIs. ( + private internal functions)
// This specific implementation uses medial internal infrastructure for holding data, models, and getting
// predictions.
//===========================================================================================================
//===========================================================================================================
//-----------------------------------------------------------------------------------
// Load() - reading a config file and initializing repository and model
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Load(const char *config_f)
{
	int rc;

	// read config and check some basic sanities
	rc = read_config(string(config_f));

	if (rc != AM_OK_RC) return rc;

	if (type_in_config_file != "MEDIAL_INFRA")
		return AM_ERROR_LOAD_NON_MATCHING_TYPE;

	if (strlen(get_name()) == 0) {
		MERR("ERROR: Name is missing\n");
		return AM_ERROR_LOAD_BAD_NAME;
	}

	// loading tester file if needed
	if (input_tester_config_file != "") {
		if (ist.read_config(input_tester_config_file) < 0) {
			MERR("ERROR: Could not read testers config file %s\n", input_tester_config_file.c_str());
			return AM_ERROR_LOAD_BAD_TESTERS_FILE;
		}
	}

	// prepare internal ma for work: set name, rep and model
	ma.set_name(get_name());
	ma.set_model_end_stage(model_end_stage);

	try {
		if (ma.init_rep_config(rep_fname.c_str()) < 0)
			return AM_ERROR_LOAD_READ_REP_ERR;
	}
	catch (...) {
		return AM_ERROR_LOAD_READ_REP_ERR;
	}

	ma.set_time_unit_env(get_time_unit());

	try {
		if (ma.init_model_from_file(model_fname.c_str()) < 0)
			return AM_ERROR_LOAD_READ_MODEL_ERR;
		if (ma.model_check_required_signals() < 0)
			return AM_ERROR_LOAD_MISSING_REQ_SIGS;
		//if (ma.init_model_for_apply() < 0)
		//	return AM_ERROR_LOAD_READ_MODEL_ERR;

	}
	catch (...) {
		return AM_ERROR_LOAD_READ_MODEL_ERR;
	}


	ma.data_load_init();
	// That's it. All is ready for data insert and prediction cycles
	is_loaded = true;
	string vers_info = ma.model_version_info();
	if (vers_info.empty())
		vers_info = "Old model without documented version!";
	MLOG("################ LOADED MODEL VERSION INFO: ##############################\n");
	MLOG("%s\n", vers_info.c_str());
	MLOG("##########################################################################\n");
	return AM_OK_RC;
}

//------------------------------------------------------------------------------------------------
// UnLoad() - clears all data, repository and model, making object ready to be deleted and freed
//------------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Unload()
{
	ClearData();
	ma.clear();
	is_loaded = false;
	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// ClearData() - clearing current data inserted inside. 
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::ClearData()
{
	ma.clear_data();
	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// AddData() - adding data for a signal with values and timestamps
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddData(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values)
{
	// At the moment MedialInfraAlgoMarker only loads timestamps given as ints.
	// This may change in the future as needed.
	int *i_times = NULL;
	vector<int> times_int;

	int tu = get_time_unit();
	if (TimeStamps_len > 0) {
		times_int.resize(TimeStamps_len);

		// currently assuming we only work with dates ... will have to change this when we'll move to other units
		for (int i = 0; i < TimeStamps_len; i++) {
			times_int[i] = AMPoint::auto_time_convert(TimeStamps[i], tu);
			if (times_int[i] < 0) {
				MERR("Error in AddData :: patient %d, signals %s, timestamp %lld is ilegal\n",
					patient_id, signalName, TimeStamps[i]);
				return AM_ERROR_ADD_DATA_FAILED;
			}
			//MLOG("time convert %ld to %d\n", TimeStamps[i], times_int[i]);
		}
		i_times = &times_int[0];
	}

	if (ma.data_load_pid_sig(patient_id, signalName, i_times, TimeStamps_len, Values, Values_len) < 0)
		return AM_ERROR_ADD_DATA_FAILED;

	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// AddDatStr() - adding data for a signal with values and timestamps
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddDataStr(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values)
{
	vector<float> converted_Values;
	vector<long long> final_tm;
	final_tm.reserve(TimeStamps_len);
	converted_Values.reserve(Values_len);
	MedRepository &rep = ma.get_rep();

	try {
		string sig = signalName;
		int section_id = rep.dict.section_id(sig);
		int sid = rep.sigs.Name2Sid[sig];
		if (rep.sigs.Sid2Info[sid].n_val_channels > 0) {
			int Values_i = 0;
			int Time_i = 0;
			const auto& category_map = rep.dict.dict(section_id)->Name2Id;
			int n_elem = (int)(Values_len / rep.sigs.Sid2Info[sid].n_val_channels);
			for (int i = 0; i < n_elem; i++) {
				bool skip_val = false;
				int val_start = Values_i;
				for (int j = 0; j < rep.sigs.Sid2Info[sid].n_val_channels; j++) {
					if (rep.sigs.is_categorical_channel(sid, j)) {
						if (category_map.find(Values[Values_i]) == category_map.end()) {
							MWARN("Found undefined code for signal \"%s\" and value \"%s\"\n",
								sig.c_str(), Values[Values_i]);
							(*ma.get_unknown_codes(patient_id))[sig].insert(Values[Values_i]);
							skip_val = true;
						}
						++Values_i;
					}
				}
				if (skip_val) {
					//remove element!
					Values_len -= rep.sigs.Sid2Info[sid].n_val_channels;
					TimeStamps_len -= rep.sigs.Sid2Info[sid].n_time_channels;
				}
				else {
					//All done 
					for (int j = 0; j < rep.sigs.Sid2Info[sid].n_time_channels; j++) {
						final_tm.push_back(TimeStamps[Time_i]);
						++Time_i;
					}
					for (int j = 0; j < rep.sigs.Sid2Info[sid].n_val_channels; j++) {
						float val;
						if (!rep.sigs.is_categorical_channel(sid, j))
							val = stof(Values[Values_i++]);
						else
							val = category_map.at(Values[val_start + j]);
						converted_Values.push_back(val);
					}
				}
			}
		}
	}
	catch (...) {
		MERR("Catched Error MedialInfraAlgoMarker::AddDataStr!!\n");
		return AM_FAIL_RC;
	}

	if (TimeStamps_len > 0 || Values_len > 0)
		return AddData(patient_id, signalName, TimeStamps_len, final_tm.data(), Values_len, converted_Values.data());
	return AM_OK_RC;

}

//-----------------------------------------------------------------------------------
// AddDataByType() : 
// Supporting loading data directly from a json
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddDataByType(int DataType, int patient_id, const char *data)
{
	if (DataType != DATA_JSON_FORMAT && DataType != DATA_BATCH_JSON_FORMAT)
		return AM_ERROR_DATA_UNKNOWN_ADD_DATA_TYPE;

	int ret_code = 0;
	if (DataType == DATA_BATCH_JSON_FORMAT) {
		vector<size_t> j_start, j_len;
		vector<char> cdata;
		get_jsons_locations(data, j_start, j_len);
		for (int j = 0; j < j_start.size(); j++) {
			if (cdata.size() < j_len[j] + 10) cdata.resize(j_len[j] + 10);
			cdata[j_len[j]] = 0;
			strncpy(&cdata[0], &data[j_start[j]], j_len[j]);
			int rc = AddDataByType(DATA_JSON_FORMAT, patient_id, &cdata[0]);
			if (rc != 0)
				ret_code = rc;
		}
		return ret_code;
	}


	json jsdata;

	try {
		jsdata = json::parse(data);
	}
	catch (...) {
		return AM_ERROR_DATA_JSON_PARSE;
	}

	int rc = AddJsonData(patient_id, jsdata);
	if (rc != 0)
		ret_code = rc;
	return ret_code;
}

void process_explainability(nlohmann::ordered_json &jattr, const Explainer_parameters &ex_params) {
	if (jattr.find("explainer_output") != jattr.end())
		for (auto &e : jattr["explainer_output"]) {
			e["contributor_description"] = "Hi";
			if (e.find("contributor_level") != e.end() && ex_params.max_threshold > 0
				&& ex_params.num_groups > 0) {
				float level_bin;
				if (ex_params.use_perc)
					level_bin = (abs(e["contributor_percentage"].get<float>()) / ex_params.max_threshold);
				else
					level_bin = (abs(e["contributor_level"].get<float>()) / ex_params.max_threshold);
				if (level_bin > 1)
					level_bin = 1;
				level_bin *= 100;
				level_bin = level_bin / (100 / ex_params.num_groups);
				if (level_bin > 0)
					level_bin = (int)(level_bin)+1;
				else
					level_bin = int(level_bin);

				if (level_bin > ex_params.num_groups)
					level_bin = ex_params.num_groups;
				level_bin *= (100 / ex_params.num_groups);

				e["contributor_level_group"] = (int)level_bin;
			}
		}
}

//------------------------------------------------------------------------------------------
// Calculate() - after data loading : get a request, get predictions, and pack as responses
//------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Calculate(AMRequest *request, AMResponses *responses)
{
	MWARN("Warning : Calculate is deprecated and will not be supported in the future, please use CalculateByType\n");
#ifdef AM_TIMING_LOGS
	MedTimer timer;
	timer.start();
#endif
	if (sort_needed) {
		if (ma.data_load_end() < 0)
			return AM_FAIL_RC;
	}
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: data_load_end %2.1f milisecond\n", timer.diff_milisec());
#endif

	if (!ma.model_initiated()) {
#ifdef AM_TIMING_LOGS
		timer.start();
#endif
		if (ma.init_model_for_apply() < 0)
			return AM_FAIL_RC;
#ifdef AM_TIMING_LOGS
		timer.take_curr_time();
		MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: init_model_for_apply %2.1f milisecond\n", timer.diff_milisec());
#endif
	}

	if (responses == NULL)
		return AM_FAIL_RC;

#ifdef AM_TIMING_LOGS
	timer.start();
#endif
	AMMessages *shared_msgs = responses->get_shared_messages();
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: get_shared_messages %2.1f milisecond\n", timer.diff_milisec());
#endif

	if (request == NULL) {
		string msg = "Error :: (" + to_string(AM_MSG_NULL_REQUEST) + " ) NULL request in Calculate()";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	//string msg_prefix = "reqId: " + string(request->get_request_id()) + " :: ";
	string msg_prefix = ""; // asked not to put reqId in messages.... (not sure it's a good idea, prev code above in comment)
	responses->set_request_id(request->get_request_id());

#ifdef AM_TIMING_LOGS
	timer.start();
#endif
	for (int i = 0; i < request->get_n_score_types(); i++) {
		char *stype = request->get_score_type(i);
		responses->insert_score_types(&stype, 1);
	}
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: get_score_type %2.1f milisecond\n", timer.diff_milisec());
#endif

	// We now have to prepare samples for the requested points
	// again - we only deal with int times in this class, so we convert the long long stamps to int
#ifdef AM_TIMING_LOGS
	timer.start();
#endif
	ma.clear_samples();
	int n_points = request->get_n_points();
	int tu = get_time_unit();
	vector<int> conv_times;

	for (int i = 0; i < n_points; i++) {
		conv_times.push_back(AMPoint::auto_time_convert(request->get_timestamp(i), tu));
		if ((ma.insert_sample(request->get_pid(i), conv_times.back()) < 0) || (conv_times.back() <= 0)) {
			string msg = msg_prefix + "(" + to_string(AM_MSG_BAD_PREDICTION_POINT) + ") Failed insert prediction point " + to_string(i) + " pid: " + to_string(request->get_pid(i)) + " ts: " + to_string(request->get_timestamp(i));
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	}

	ma.normalize_samples();
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: prepared_samples %2.1f milisecond\n", timer.diff_milisec());

	timer.start();
#endif
	// Checking score types and verify they are supported
	int n_score_types = request->get_n_score_types();
	for (int i = 0; i < n_score_types; i++) {
		if (!IsScoreTypeSupported(request->get_score_type(i))) {
			//string msg = msg_prefix + "(" + to_string(AM_MSG_BAD_SCORE_TYPE) + ") AlgoMarker of type " + string(get_name()) + " does not support score type " + string(request->get_score_type(i));
			string msg = msg_prefix + "AlgoMarker of type " + string(get_name()) + " does not support score type " + string(request->get_score_type(i));
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	}

	// At this stage we need to create a response entry for each of the requested points
	// Then we have to test for eligibility - err the ones that are not eligible
	// And then score all the eligible ones in a single batch.
	vector<int> eligible_pids, eligible_timepoints;
	vector<long long> eligible_ts;
	MedPidRepository &rep = ma.get_rep();
	unordered_map<unsigned long long, vector<long long>> sample2ts; // conversion of each sample to all the ts that were mapped to it.

	int n_bad_scores = 0;
	for (int i = 0; i < n_points; i++) {
		int _pid = request->get_pid(i);
		long long _ts = request->get_timestamp(i);

		// create a response
		AMResponse *res = responses->create_point_response(_pid, _ts);
		//AMResponse *res = responses->get_response_by_point(_pid, (long long)conv_times[i]);
//		if (res == NULL)
//			res = responses->create_point_response(_pid, _ts);
//		res = responses->create_point_response(_pid, (long long)conv_times[i]);

		// test this point for eligibility and add errors if needed
		vector<InputSanityTesterResult> test_res;
		int test_rc = ist.test_if_ok(_pid, (long long)conv_times[i], *ma.get_unknown_codes(_pid), test_res);
		int test_rc2 = ist.test_if_ok(rep, _pid, (long long)conv_times[i], test_res);
		if (test_rc2 < 1)
			test_rc = test_rc2;

		// push messages if there are any
		AMMessages *msgs = res->get_msgs();
		for (auto &tres : test_res) {
			//string msg = msg_prefix + tres.err_msg + " Internal Code: " + to_string(tres.internal_rc);
			string msg = msg_prefix + tres.err_msg; // messages without Internal codes...(prev code in comment above).
			msgs->insert_message(tres.external_rc, msg.c_str());
		}

		if (test_rc <= 0) {
			n_bad_scores++;
		}
		else {
			//MLOG("DEBUG ===> i %d _pid %d conv %d _ts %lld size %d\n", i, _pid, conv_times[i], _ts, eligible_pids.size());
			eligible_pids.push_back(_pid);
			eligible_timepoints.push_back(conv_times[i]);
			eligible_ts.push_back(_ts);
			unsigned long long p = ((unsigned long long)_pid << 32) | (conv_times[i]);
			if (sample2ts.find(p) == sample2ts.end()) sample2ts[p] = vector<long long>();
			sample2ts[p].push_back(_ts);

		}

	}

	int _n_points = (int)eligible_pids.size();

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: tested_eligibility %2.1f milisecond. Has %d samples\n", timer.diff_milisec(), _n_points);

	ma.model_apply_verbose(true);
	timer.start();
#endif

	// Calculating raw scores for eligble points
	vector<float> raw_scores(_n_points, (float)AM_UNDEFINED_VALUE);
	int get_preds_rc = -1;
	try {
		if ((get_preds_rc = ma.get_preds(&eligible_pids[0], &eligible_timepoints[0], &raw_scores[0], _n_points)) < 0) {
			string msg = msg_prefix + "(" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting RAW scores in AlgoMarker " + string(get_name()) + " With return code " + to_string(get_preds_rc);
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	}
	catch (...) {
		string msg = msg_prefix + "(" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting RAW scores in AlgoMarker " + string(get_name()) + " caught a crash. ";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	if (am_matrix != "" && _n_points > 0) { // debug only
		if (first_write)
			ma.write_features_mat(am_matrix);
		else
			ma.add_features_mat(am_matrix);
		first_write = false;
	}

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: get_preds %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif
	// going over scores, and adding them to the right responses
	char **_score_types;
	int _n_score_types;
	responses->get_score_types(&_n_score_types, &_score_types);

	MedSamples *meds = ma.get_samples_ptr();
	Explainer_parameters ex_params;
	ma.get_explainer_params(ex_params);

	for (auto &id_s : meds->idSamples) {
		for (auto &s : id_s.samples) {

			// basic info for current sample
			int c_pid = s.id;
			int c_ts = s.time;
			float c_scr = s.prediction.size() > 0 ? s.prediction[0] : (float)AM_UNDEFINED_VALUE;
			string c_ext_scr = "";
			if (s.str_attributes.size() > 0) {
				nlohmann::ordered_json c_ext_scr_json({});

				for (auto &ex_res_field_name : extended_result_fields) {
					if (s.str_attributes.count(ex_res_field_name)) {
						c_ext_scr_json[ex_res_field_name] = nlohmann::ordered_json::parse(s.str_attributes[ex_res_field_name]);
						process_explainability(c_ext_scr_json[ex_res_field_name], ex_params);
					}
				}
				c_ext_scr = c_ext_scr_json.dump();
			}

			unsigned long long p = ((unsigned long long)c_pid << 32) | (c_ts);

			for (auto ts : sample2ts[p]) {

				// DEBUG
				//for (auto &attr : s.attributes) MLOG("pid %d time %d score %f attr %s %f\n", c_pid, c_ts, c_scr, attr.first.c_str(), attr.second);

				// get the matching response (should be prepared already)
				AMResponse *res = responses->get_response_by_point(c_pid, ts);

				if (res != NULL) {

					// we now test the attribute tests
					vector<InputSanityTesterResult> test_res;
					int test_rc = ist.test_if_ok(s, test_res);

					AMMessages *msgs = res->get_msgs();
					for (auto &tres : test_res) {
						//string msg = msg_prefix + tres.err_msg + " Internal Code: " + to_string(tres.internal_rc);
						string msg = msg_prefix + tres.err_msg; // no Internal Code message (prev code in comment above).
						//MLOG("Inserting attr error to pid %d ts %d : %d : %s\n", c_pid, ts, tres.external_rc, msg.c_str());
						msgs->insert_message(tres.external_rc, msg.c_str());
					}

					if (test_rc <= 0)
						n_bad_scores++;
					else {

						// all is fine, we insert the score into its place
						res->init_scores(_n_score_types);

						for (int j = 0; j < _n_score_types; j++) {

							if (strcmp(_score_types[j], "Raw") == 0) {
								res->set_score(j, c_scr, _score_types[j], c_ext_scr);
							}
							else {
								res->set_score(j, (float)AM_UNDEFINED_VALUE, _score_types[j], "");
								AMScore *am_scr = res->get_am_score(j);
								AMMessages *msgs = am_scr->get_msgs();
								string msg = msg_prefix + "Undefined Score Type: " + string(_score_types[j]);
								msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
							}

						}
					}

				}
			}
		}
	}

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::Calculate :: finished_response %2.1f milisecond\n", timer.diff_milisec());
#endif

	if (n_bad_scores > 0) {
		string msg = msg_prefix + "Failed input tests for " + to_string(n_bad_scores) + " out of " + to_string(n_points) + " scores";
		if (n_bad_scores < n_points) {
			shared_msgs->insert_message(AM_RESPONSES_ELIGIBILITY_ERROR, msg.c_str());
			return AM_OK_RC;
		}

		shared_msgs->insert_message(AM_RESPONSES_ELIGIBILITY_ERROR, msg.c_str());
		return AM_FAIL_RC;
	}

	return AM_OK_RC;
		}


//------------------------------------------------------------------------------------------
// CalculateByType : alllows for a general json in -> json out API with many more options 
//------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::CalculateByType(int CalculateType, char *request, char **response)
{
	if (CalculateType != JSON_REQ_JSON_RESP)
		return AM_FAIL_RC;

#ifdef AM_TIMING_LOGS
	ma.model_apply_verbose(true);
	MedTimer timer;
	timer.start();
#endif
	if (sort_needed) {
		if (ma.data_load_end() < 0)
			return AM_FAIL_RC;
	}
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: data_load_end %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	if (!ma.model_initiated()) {
		if (ma.init_model_for_apply() < 0)
			return AM_FAIL_RC;
	}
#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: init_model_for_apply %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	json jreq;
	nlohmann::ordered_json jresp;

	jresp = nlohmann::ordered_json({ { "type", "response" } });
	try {
		jreq = json::parse(request);
	}
	catch (...) {
		add_to_json_array(jresp, "errors", "ERROR: Could not parse request as a valid json");
		json_to_char_ptr(jresp, response);
		return AM_FAIL_RC;
	}

	// verify the "type" : "request" , and the "request_id" : something fields
	if (!json_verify_key(jreq, "type", 1, "request"))
		add_to_json_array(jresp, "errors", "ERROR: missing type request");
	string request_id;
	if (!json_verify_key(jreq, "request_id", 0, ""))
		add_to_json_array(jresp, "errors", "ERROR: no request_id provided");
	else {
		request_id = jreq["request_id"].get<string>();
		jresp.push_back({ "request_id", request_id });
	}

	if (!json_verify_key(jreq, "requests", 0, ""))
		add_to_json_array(jresp, "errors", "ERROR: missing actual requests in request " + request_id);

	if (json_verify_key(jresp, "errors", 0, "")) { json_to_char_ptr(jresp, response); return AM_FAIL_RC; } // Leave now if there are errors

	// default parameters
	json_req_info defaults;

	vector<json_req_info> sample_reqs;

	//	try {
	json_parse_request(jreq, defaults, defaults);

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: json_parse_request %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	for (auto &jreq_i : jreq["requests"]) {
		json_req_info j_i;
		json_parse_request(jreq_i, defaults, j_i);
		sample_reqs.push_back(j_i);
		if (j_i.load_data && json_verify_key(jreq_i, "data", 0, "")) {
			if (AddJsonData(j_i.sample_pid, jreq_i["data"]) != AM_OK_RC) {
				add_to_json_array(jresp, "errors", "ERROR: error when loading data for patient id " + to_string(j_i.sample_pid));
			}
		}
	}
	//	}

	//	catch (...) {
	//		add_to_json_array(jresp, "errors", "ERROR: error when parsing requests (or loading)");
	//	}
	if (json_verify_key(jresp, "errors", 0, "")) { json_to_char_ptr(jresp, response); return AM_FAIL_RC; } // Leave now if there are errors

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: end load_data %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	// We now convert times and do an initial sanity checks
	// again - we only deal with int times in this class, so we convert the long long stamps to int
	// we also run the eligibility tests, keep the results, and make lists of all eligible points for scoring.
	int n_points = (int)sample_reqs.size();
	int tu = get_time_unit();
	MedPidRepository &rep = ma.get_rep();
	vector<int> eligible_pids, eligible_timepoints;
	unordered_map<unsigned long long, vector<int>> sample2ind; // conversion of each sample to all the ts that were mapped to it.
	int n_failed = 0, n_bad = 0;
	//#pragma omp parallel for if (n_points > 10)
	for (int i = 0; i < n_points; i++) {
		json_req_info &req_i = sample_reqs[i];
		req_i.conv_time = AMPoint::auto_time_convert(req_i.sample_time, tu);
		int ok_time = 1;
		if (tu == MedTime::Date && (req_i.conv_time < 19500000 || req_i.conv_time > 30000000)) ok_time = 0;
		if ((req_i.sample_pid <= 0) || (req_i.conv_time <= 0) || ok_time == 0) {
#pragma omp critical 
			{
				add_to_json_array(jresp, "errors", "ERROR: BAD request patient id or time : failed in inserting pid: " + to_string(req_i.sample_pid) + " ts: " + to_string(req_i.sample_time));
				req_i.sanity_test_rc = -2;
				n_failed++;
			}
		}

		else {
			try {
				req_i.sanity_test_rc = ist.test_if_ok(req_i.sample_pid, (long long)req_i.conv_time, *ma.get_unknown_codes(req_i.sample_pid), req_i.sanity_res);
				int rc_res = ist.test_if_ok(rep, req_i.sample_pid, (long long)req_i.conv_time, req_i.sanity_res);
				if (rc_res < 1)
					req_i.sanity_test_rc = rc_res;
			}
			catch (...) {
				req_i.sanity_caught_err = 1;
			}
		}

#pragma omp critical
		if (req_i.sanity_caught_err == 0 && req_i.sanity_test_rc > 0) {
			unsigned long long p = ((unsigned long long)req_i.sample_pid << 32) | req_i.conv_time;
			if (sample2ind.find(p) == sample2ind.end()) {
				eligible_pids.push_back(req_i.sample_pid);
				eligible_timepoints.push_back(req_i.conv_time);
				sample2ind[p] = vector<int>();
			}
			sample2ind[p].push_back(i);
		}
		else
			n_bad++;
	}

	if (n_failed > 0) { json_to_char_ptr(jresp, response);	return AM_FAIL_RC; }

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: end calc_eligiblilty %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	// at this point in time we are ready to score eligible_pids,eligible_timepoints. We will do that, and later wrap it all up into a single json back.
	int _n_points = (int)eligible_pids.size();
	int get_preds_rc = -1;
	try {
		if ((get_preds_rc = ma.get_preds(&eligible_pids[0], &eligible_timepoints[0], NULL, _n_points)) < 0) {
			add_to_json_array(jresp, "errors", "ERROR: (" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting scores in AlgoMarker " + string(get_name()) + " With return code " + to_string(get_preds_rc));
			json_to_char_ptr(jresp, response);
			return AM_FAIL_RC;
		}
	}
	catch (...) {
		add_to_json_array(jresp, "errors", "ERROR: (" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting scores in AlgoMarker " + string(get_name()) + " caught a crash");
		json_to_char_ptr(jresp, response);
		return AM_FAIL_RC;
	}

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: end get_preds %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	// now we are ready ... we have the results, and we need to put it all into the response json one by one.

	// adding result samples pointers to sample_reqs
	MedSamples *samps = ma.get_samples_ptr();
	for (auto &ids : samps->idSamples)
		for (auto &s : ids.samples) {
			unsigned long long p = ((unsigned long long)s.id << 32) | s.time;
			for (auto j : sample2ind[p])
				sample_reqs[j].res = &s;
		}

	jresp.push_back({ "responses", json::array() });
	for (int i = 0; i < sample_reqs.size(); i++) {
		json_req_info &req_i = sample_reqs[i];
		if (req_i.res != NULL) {


			try {
				int test_rc = ist.test_if_ok(*(req_i.res), req_i.sanity_res);
				if (test_rc == 0) req_i.sanity_test_rc = 0;
			}
			catch (...) {
				req_i.sanity_caught_err = 1;
			}
		}
		Explainer_parameters ex_params;
		ma.get_explainer_params(ex_params);
		//MLOG("=====> Working on i %d pid %d time %d sanity_test_rc %d sanity_caught_err %d\n", i, req_i.sample_pid, req_i.sample_time, req_i.sanity_test_rc, req_i.sanity_caught_err);
		nlohmann::ordered_json js = nlohmann::ordered_json({});

		js.push_back({ "patient_id" , to_string(req_i.sample_pid) });
		js.push_back({ "time" , to_string(req_i.sample_time) });
		if (req_i.sanity_caught_err)
			add_to_json_array(js, "messages", "ERROR: sanity tests crashed");
		if (req_i.sanity_res.size() > 0)
			for (auto &its : req_i.sanity_res) {
				add_to_json_array(js, "messages", "(" + to_string(its.external_rc) + ")" + its.err_msg);
			}
		for (auto &e : req_i.exports) {
			if (e.second.type == PREDICTION_SOURCE_PREDICTIONS) {
				if (req_i.res != NULL && req_i.res->prediction.size() > e.second.pred_channel && req_i.sanity_caught_err == 0 && req_i.sanity_test_rc > 0)
					js.push_back({ e.first, to_string(req_i.res->prediction[e.second.pred_channel]) });
				else
					js.push_back({ e.first, to_string(AM_UNDEFINED_VALUE) });
				if (req_i.res == NULL)
					add_to_json_array(js, "messages", "ERROR: did not get result for field " + e.first + " : " + e.second.field);
				else if (req_i.res->prediction.size() <= e.second.pred_channel || e.second.pred_channel < 0)
					add_to_json_array(js, "messages", "ERROR: prediction channel " + to_string(e.second.pred_channel) + " is illegal");
			}
			else if (e.second.type == PREDICTION_SOURCE_ATTRIBUTE && req_i.res != NULL) {
				if (req_i.res->attributes.find(e.second.field) != req_i.res->attributes.end())
					js.push_back({ e.first, to_string(req_i.res->attributes[e.second.field]) });
				else if (req_i.res->str_attributes.find(e.second.field) != req_i.res->str_attributes.end())
					js.push_back({ e.first, req_i.res->str_attributes[e.second.field] });
			}
			else if (e.second.type == PREDICTION_SOURCE_ATTRIBUTE_AS_JSON && req_i.res != NULL) {
				if (req_i.res->str_attributes.find(e.second.field) != req_i.res->str_attributes.end()) {
					nlohmann::ordered_json jattr;
					try {
						jattr = nlohmann::ordered_json::parse(req_i.res->str_attributes[e.second.field]);
						process_explainability(jattr, ex_params);
						js.push_back({ e.first, jattr });
					}
					catch (...) {
						add_to_json_array(jresp, "messages", "ERROR: could not parse attribute " + e.second.field + " as a valid json");
					}
				}
			}
			else if (e.second.type == PREDICTION_SOURCE_JSON) {
				if (req_i.res->jrec.find(e.second.field) != req_i.res->jrec.end()) {
					auto &jj = req_i.res->jrec[e.second.field];
					js.push_back({ e.first, jj });
				}
			}
		}
		jresp["responses"].push_back(js);
	}

	json_to_char_ptr(jresp, response);

#ifdef AM_TIMING_LOGS
	timer.take_curr_time();
	MLOG("INFO:: MedialInfraAlgoMarker::CalculateByType :: end create_response %2.1f milisecond\n", timer.diff_milisec());
	timer.start();
#endif

	if (am_matrix != "" && _n_points > 0) { // debug only
		if (first_write)
			ma.write_features_mat(am_matrix);
		else
			ma.add_features_mat(am_matrix);
		first_write = false;
	}

	return AM_OK_RC;
	}

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AdditionalLoad(const int LoadType, const char *load)
{
	if (!is_loaded) return AM_ERROR_MUST_BE_LOADED;
	if (LoadType != LOAD_DICT_FROM_FILE && LoadType != LOAD_DICT_FROM_JSON)
		return AM_ERROR_UNKNOWN_LOAD_TYPE;

	json js;

	if (LoadType == LOAD_DICT_FROM_FILE) {
		string sload;
		string f_in(load);
		if (read_file_into_string(f_in, sload) < 0)
			return AM_ERROR_READING_DICT_FILE;
		js = json::parse(sload.c_str());
	}
	else
		js = json::parse(load);

	try {
		ma.add_json_dict(js);
	}
	catch (...) {
		return AM_ERROR_PARSING_JSON_DICT;
	}

	// now that we added the json dictionary, we need to reinitialize the model ! as it needs to prepare potential tables using these new definitions and sets
	//ma.init_model_for_apply();

	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// private internals for class MedialInfraAlgoMarker
//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::read_config(const string &conf_f)
{
	set_config(conf_f.c_str());

	ifstream inf(conf_f);

	if (!inf)
		return AM_ERROR_LOAD_NO_CONFIG_FILE;

	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			if (fields.size() >= 2) {
				if (fields[0] == "TYPE") type_in_config_file = fields[1];
				else if (fields[0] == "REPOSITORY") rep_fname = fields[1];
				else if (fields[0] == "MODEL") model_fname = fields[1];
				else if (fields[0] == "MODEL_END_STAGE") try { model_end_stage = stoi(fields[1]); }
				catch (...) { MTHROW_AND_ERR("Could not parse given value MODEL_END_STAGE='%s'\n", fields[1].c_str()); }
				else if (fields[0] == "EXTENDED_RESULT_FIELDS") split(extended_result_fields, fields[1], boost::is_any_of(";"));
				else if (fields[0] == "INPUT_TESTER_CONFIG") input_tester_config_file = fields[1];
				else if (fields[0] == "NAME")  set_name(fields[1].c_str());
				else if (fields[0] == "TIME_UNIT") {
					set_time_unit(med_time_converter.string_to_type(fields[1].c_str()));
				}
				else if (fields[0] == "DEBUG_MATRIX")  am_matrix = fields[1];
				else if (fields[0] == "AM_UDI_DI")  set_am_udi_di(fields[1].c_str());
				else if (fields[0] == "AM_VERSION")  set_am_version(fields[1].c_str());
				else if (fields[0] == "EXPLAINABILITY_PARAMS") ma.set_explainer_params(fields[1]);
				else if (fields[0] == "TESTER_NAME") {}
				else if (fields[0] == "FILTER") {}
				else MWARN("WRAN: unknown parameter \"%s\". Read and ignored\n", fields[0].c_str());
			}
		}
	}

	string dir = conf_f.substr(0, conf_f.find_last_of("/\\"));
	if (rep_fname != "" && rep_fname[0] != '/' && rep_fname[0] != '\\') {
		// relative path
		rep_fname = dir + "/" + rep_fname;
	}

	if (model_fname != "" && model_fname[0] != '/' && model_fname[0] != '\\') {
		// relative path
		model_fname = dir + "/" + model_fname;
	}

	if (input_tester_config_file == ".") {
		input_tester_config_file = conf_f;  // option to use the general config file as the file to config the tester as well.
	}
	else if (input_tester_config_file != "" && input_tester_config_file[0] != '/' && input_tester_config_file[0] != '\\') {
		// relative path
		input_tester_config_file = dir + "/" + input_tester_config_file;
	}


	return AM_OK_RC;
}

// maximal input of 32GB
#define MAX_POSSIBLE_STRING_LEN ((size_t)1 << 35)
//-----------------------------------------------------------------------------------
void MedialInfraAlgoMarker::get_jsons_locations(const char *data, vector<size_t> &j_start, vector<size_t> &j_len)
{
	j_start.clear();
	j_len.clear();

	size_t j = 0;
	int counter = 0;
	int in_string = 0;
	size_t start = 0;
	size_t len = 0;
	while (data[j] != 0 || j > MAX_POSSIBLE_STRING_LEN) {
		char ch = data[j];
		if (ch == '\"' || ch == '\'') in_string = 1 - in_string;
		if ((!in_string) && ch == '{') {
			if (counter == 0) start = j;
			counter++;
		}
		if (counter > 0) len++;
		if ((!in_string) && ch == '}') counter--;

		if (counter == 0 && len > 0) {
			j_start.push_back(start);
			j_len.push_back(len);
			len = 0;
			if (j_start.size() > 0 && j_start.size() % 1000 == 0)
				MLOG("Found %d jsons so far\n", j_start.size());
		}
		if (counter < 0) MTHROW_AND_ERR("Mismatch in {} count in jsons string\n");

		j++;
	}

	MLOG("Read %d jsons from data string (debug info: counter = %d j = %ld)\n", j_start.size(), counter, j);
}

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddJsonData(int patient_id, json &j_data)
{
	try {
		json &js = j_data;

		// supporting also older style jsons that were embeded in a "body" section
		if (j_data.find("body") != j_data.end())
			js = j_data["body"];

		if (patient_id <= 0) {
			// in this case we take the patient id directly from the json itself
			if (js.find("patient_id") != js.end())
				patient_id = js["patient_id"].get<long long>();
			else if (js.find("pid") != js.end())
				patient_id = js["pid"].get<long long>();
		}

		//MLOG("Loading with pid %d\n", patient_id);

		vector<long long> times;
		int s_data_size = 100000;
		vector<char> sdata(s_data_size);
		vector<int> sinds;
		int curr_s = 0;
		//char str_values[MAX_VALS][MAX_VAL_LEN];
		for (auto &s : js["signals"]) {
			times.clear();
			sinds.clear();
			curr_s = 0;
			string sig = s["code"].get<string>();
			int n_time_channels, n_val_channels, *is_categ;
			get_sig_structure(sig, n_time_channels, n_val_channels, is_categ);
			//MLOG("%s %d %d\n", sig.c_str(), n_time_channels, n_val_channels);
			int n_data = 0;
			for (auto &d : s["data"]) {
				//MLOG("time ");
				int nt = 0;
				for (auto &t : d["timestamp"]) {
					times.push_back(t.get<long long>());
					nt++;
					//MLOG("%d ", itime);
				}
				//MLOG("val ");
				int nv = 0;
				for (auto &v : d["value"]) {
					string sv = v.get<string>().c_str();
					int slen = (int)sv.length();
					//MLOG("val %d : %s len: %d curr_s %d s_data_size %d %d n_val_channels %d\n", nv, sv.c_str(), slen, curr_s, s_data_size, sdata.size(), n_val_channels);
					if (curr_s + 1 + slen > s_data_size) {
						s_data_size *= 2;
						sdata.resize(s_data_size);
					}
					if (nv < n_val_channels) {
						sv.copy(&sdata[curr_s], slen);
						sdata[curr_s + slen] = 0;
						sinds.push_back(curr_s);
						curr_s += slen + 1;
						nv++;
					}
					//char *sp = &sdata[sinds.back()];
					//MLOG("val %d %d %s : %s len: %d curr_s %d s_data_size %d %d\n", sinds.size(), sinds.back(), sp, sv.c_str(), slen, curr_s, s_data_size, sdata.size());
					//MLOG("%s ", v.get<string>().c_str());
				}
				//MLOG("\n");
				n_data++;
			}
			vector<char *> p_str;
			for (auto j : sinds)
				p_str.push_back(&sdata[j]);
			long long *p_times = &times[0];
			int n_times = (int)times.size();
			char **str_values = &p_str[0];
			int n_vals = (int)p_str.size();

			//MLOG("%s n_times %d n_vals %d n_data %d\n", sig.c_str(), n_times, n_vals, n_data);
			//MLOG("times: "); for (int j = 0; j < n_times; j++) MLOG("%d,", p_times[j]); 	MLOG("\nvals: ");
			//for (int j = 0; j < n_vals; j++) MLOG("%s, ", str_values[j]); MLOG("\n");

			if (AddDataStr(patient_id, sig.c_str(), n_times, p_times, n_vals, str_values) != AM_OK_RC) {
				MLOG("Failed AddDataStr() in AddDataByType()\n");
				return AM_FAIL_RC;
			}
		}
	}
	catch (...) {
		return AM_FAIL_RC;
	}

	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
string MedialInfraAlgoMarker::get_lib_code_version() {
	return medial::get_git_version();
}

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Discovery(char **response) {
	nlohmann::ordered_json jresp;
	jresp = nlohmann::ordered_json(
		{
		{ "name", get_name() },
		{ "version", get_am_version() },
		{ "udi_di", get_am_udi_di() },
		{ "model_code_version", ma.model_version_info() },
		{ "lib_code_version", get_lib_code_version() },
		{ "signals", nlohmann::ordered_json::array() }
		}
	);

	nlohmann::ordered_json &json_signals = jresp["signals"];

	//Add signals to json_signals
	if (!ma.model_initiated())
		ma.init_model_for_rep();

	unordered_set<string> sigs;
	get_am_rep_signals(sigs);
	vector<string> req_sigs;
	unordered_map<string, vector<string>> sig_categ, sig_categ_final;
	ma.get_model_signals_info(req_sigs, sig_categ);
	unordered_set<string> req_set(req_sigs.begin(), req_sigs.end());
	//rename sig_categ if needed in some cases!
	for (const auto &it : sig_categ)
	{
		if (req_set.find(it.first) != req_set.end())
			sig_categ_final[it.first] = it.second;
		else {
			//test by section name!
			int sig_section = ma.get_rep().dict.section_id(it.first);
			//retrieve name that exists in SECTION:
			const unordered_set<string> &sec_names = ma.get_rep().dict.dict(sig_section)->section_name;
			vector<string> candidates;
			for (const string &cand : sec_names)
			{
				if (req_set.find(cand) == req_set.end())
					continue;
				candidates.push_back(cand);
			}

			//Add to all:
			if (candidates.empty())
				MWARN("Warn - has used categorical signal %s without mapping\n",
					it.first.c_str());
			else if (candidates.size() > 1)
				MWARN("Warn - has used categorical signal %s with multiple mapping\n",
					it.first.c_str());
			else {
				unordered_set<string> sig_list_c(sig_categ_final[candidates[0]].begin(), sig_categ_final[candidates[0]].end());
				sig_list_c.insert(it.second.begin(), it.second.end());
				vector<string> final_list(sig_list_c.begin(), sig_list_c.end());
				sort(final_list.begin(), final_list.end());
				sig_categ_final[candidates[0]] = move(final_list);
			}
		}
	}

	for (const string &sig_name : req_sigs)
	{
		string sig_nm = sig_name;
		string sig_unit = "";
		string sig_type = "";
		vector<int> categorical_ch;
		vector<string> categorical_vals;
		if (sigs.find(sig_name) != sigs.end()) {
			int sid = ma.get_rep().sigs.sid(sig_name);
			const SignalInfo &si = ma.get_rep().sigs.Sid2Info[sid];
			for (int i = 0; i < si.n_val_channels; ++i) {
				if (i > 0)
					sig_unit += ",";
				sig_unit += si.unit_of_measurement_per_val_channel[i];
			}
			UniversalSigVec usv;
			usv.init(si);
			sig_type = usv.get_signal_generic_spec();
			categorical_ch.resize(si.n_val_channels);
			for (size_t i = 0; i < si.n_val_channels; ++i)
				if (si.is_categorical_per_val_channel[i])
					categorical_ch[i] = 1;

		}
		else {
			sig_nm = sig_nm + "(virtual)";
		}
		if (sig_categ_final.find(sig_name) != sig_categ_final.end())
			categorical_vals = std::move(sig_categ_final.at(sig_name));

		nlohmann::ordered_json sig_js;
		sig_js = {
			{ "code", sig_nm },
			{ "unit", sig_unit },
			{ "type", sig_type },
			{ "categorical_channels", categorical_ch },
			{ "categorical_values", categorical_vals }
		};

		json_signals += sig_js;
	}

	//ma.get_rep().sigs.Sid2Info[1].

	json_to_char_ptr(jresp, response);
	return 0;
}

//-----------------------------------------------------------------------------------
void add_to_json_array(json &js, const string &key, const string &s_add)
{
	if (js.find(key) == js.end())
		js.push_back({ key, json::array() });
	js[key] += s_add;
}
void add_to_json_array(nlohmann::ordered_json &js, const string &key, const string &s_add)
{
	if (js.find(key) == js.end())
		js.push_back({ key, json::array() });
	js[key] += s_add;
}


//-----------------------------------------------------------------------------------
void json_to_char_ptr(json &js, char **jarr)
{
	*jarr = NULL;
	string sj = js.dump(1, '\t');

	*jarr = new char[sj.length() + 1];

	if (*jarr != NULL) {
		(*jarr)[sj.length()] = 0;
		strncpy(*jarr, sj.c_str(), sj.length());
	}
}

void json_to_char_ptr(nlohmann::ordered_json &js, char **jarr)
{
	*jarr = NULL;
	string sj = js.dump(1, '\t');

	*jarr = new char[sj.length() + 1];

	if (*jarr != NULL) {
		(*jarr)[sj.length()] = 0;
		strncpy(*jarr, sj.c_str(), sj.length());
	}
}

//-----------------------------------------------------------------------------------
bool json_verify_key(json &js, const string &key, int verify_val_flag, const string &val)
{
	bool is_in = false;
	if (js.find(key) != js.end()) is_in = true;

	if (is_in && verify_val_flag) {
		if (js[key].get<string>() != val)
			is_in = false;
	}

	return is_in;
}
bool json_verify_key(nlohmann::ordered_json &js, const string &key, int verify_val_flag, const string &val)
{
	bool is_in = false;
	if (js.find(key) != js.end()) is_in = true;

	if (is_in && verify_val_flag) {
		if (js[key].get<string>() != val)
			is_in = false;
	}

	return is_in;
}


//------------------------------------------------------------------------------------------
int json_parse_request(json &jreq, json_req_info &defaults, json_req_info &req_i)
{
	req_i = defaults;
	// read defaults (if exist)
	if (json_verify_key(jreq, "patient_id", 0, "") || json_verify_key(jreq, "pid", 0, "")) {
		if (json_verify_key(jreq, "patient_id", 0, ""))
			req_i.sample_pid = stoi(jreq["patient_id"].get<string>());
		else
			req_i.sample_pid = stoi(jreq["pid"].get<string>());
	}

	if (json_verify_key(jreq, "scoreOnDate", 0, "") || json_verify_key(jreq, "time", 0, "")) {
		if (json_verify_key(jreq, "scoreOnDate", 0, ""))
			req_i.sample_time = stoll(jreq["scoreOnDate"].get<string>());
		else
			req_i.sample_time = stoll(jreq["time"].get<string>());
	}

	if (json_verify_key(jreq, "load", 0, "")) {
		req_i.load_data = stoi(jreq["load"].get<string>());
	}

	if (json_verify_key(jreq, "export", 0, "")) {

		for (auto &jexp : jreq["export"].items()) {

			string name = jexp.key();
			string field = jexp.value().get<string>();

			//MLOG("Working on %s : %s\n", name.c_str(), field.c_str());
			int type = PREDICTION_SOURCE_UNKNOWN;
			int pred_channel = -1;

			vector<string> f;
			boost::split(f, field, boost::is_any_of(" "));
			if (f.size() == 2) {
				if (f[0] == "attr") { type = PREDICTION_SOURCE_ATTRIBUTE; field = f[1]; }
				else if (f[0] == "json_attr") { type = PREDICTION_SOURCE_ATTRIBUTE_AS_JSON; field = f[1]; }
				else if (f[0] == "pred") { type = PREDICTION_SOURCE_PREDICTIONS; field = "pred_" + f[1]; }
				else if (f[0] == "json") { type = PREDICTION_SOURCE_JSON; field = f[1]; }
			}

			if ((type == PREDICTION_SOURCE_UNKNOWN || type == PREDICTION_SOURCE_PREDICTIONS) && (field.length() > 5) && (field.substr(0, 5) == "pred_")) {
				type = PREDICTION_SOURCE_PREDICTIONS;
				pred_channel = stoi(field.substr(5));
			}

			if (type == PREDICTION_SOURCE_UNKNOWN) type = PREDICTION_SOURCE_ATTRIBUTE;

			json_req_export jexport;

			jexport.field = field;
			jexport.pred_channel = pred_channel;
			jexport.type = type;

			req_i.exports[name] = jexport;

		}

	}

	return 0;
}
