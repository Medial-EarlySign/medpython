#include "AlgoMarker.h"

#include <Logger/Logger/Logger.h>
#include <MedTime/MedTime/MedTime.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL



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
		MERR("ERROR: Name is %s\n", get_name());
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
	}
	catch (...) {
		return AM_ERROR_LOAD_READ_MODEL_ERR;
	}


	ma.data_load_init();
	// That's it. All is ready for data insert and prediction cycles
	return AM_OK_RC;
}

//------------------------------------------------------------------------------------------------
// UnLoad() - clears all data, repository and model, making object ready to be deleted and freed
//------------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Unload()
{
	ClearData();
	ma.clear();
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
		for (int i=0; i<TimeStamps_len; i++) {
			times_int[i] = AMPoint::auto_time_convert(TimeStamps[i], tu);
		}
		i_times = &times_int[0];
	}

	if (ma.data_load_pid_sig(patient_id, signalName, i_times, TimeStamps_len, Values, Values_len) < 0)
		return AM_ERROR_ADD_DATA_FAILED;

	return AM_OK_RC;
}

//------------------------------------------------------------------------------------------
// Calculate() - after data loading : get a request, get predictions, and pack as responses
//------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Calculate(AMRequest *request, AMResponses *responses)
{
	if (sort_needed) {
		if (ma.data_load_end() < 0)
		return AM_FAIL_RC;
	}

	if (responses == NULL)
		return AM_FAIL_RC;

	AMMessages *shared_msgs = responses->get_shared_messages();

	if (request == NULL) {
		string msg = "Error :: (" + to_string(AM_MSG_NULL_REQUEST) + " ) NULL request in Calculate()";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	//string msg_prefix = "reqId: " + string(request->get_request_id()) + " :: ";
	string msg_prefix = ""; // asked not to put reqId in messages.... (not sure it's a good idea, prev code above in comment)
	responses->set_request_id(request->get_request_id());

	for (int i=0; i<request->get_n_score_types(); i++) {
		char *stype = request->get_score_type(i);
		responses->insert_score_types(&stype, 1);
	}


	// We now have to prepare samples for the requested points
	// again - we only deal with int times in this class, so we convert the long long stamps to int
	ma.clear_samples();
	int n_points = request->get_n_points();
	int tu = get_time_unit();
	vector<int> conv_times;

	for (int i=0; i<n_points; i++) {
		conv_times.push_back(AMPoint::auto_time_convert(request->get_timestamp(i), tu));
		if ((ma.insert_sample(request->get_pid(i), conv_times.back()) < 0) || (conv_times.back() <= 0)) {
			string msg = msg_prefix + "(" + to_string(AM_MSG_BAD_PREDICTION_POINT) + ") Failed insert prediction point " + to_string(i) + " pid: " + to_string(request->get_pid(i)) + " ts: " + to_string(request->get_timestamp(i));
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	}

	ma.normalize_samples();

	// Checking score types and verify they are supported
	int n_score_types = request->get_n_score_types();
	for (int i=0; i<n_score_types; i++) {
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
	MedRepository &rep = ma.get_rep();
	unordered_map<unsigned long long, vector<long long>> sample2ts; // conversion of each sample to all the ts that were mapped to it.

	int n_bad_scores = 0;
	for (int i=0; i<n_points; i++) {
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
		int test_rc = ist.test_if_ok(rep, _pid, (long long)conv_times[i], test_res);

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

	// Calculating raw scores for eligble points
	vector<float> raw_scores(_n_points, (float)AM_UNDEFINED_VALUE);
	int get_preds_rc;
	if ((get_preds_rc = ma.get_preds(&eligible_pids[0], &eligible_timepoints[0], &raw_scores[0], _n_points)) < 0) {
		string msg = msg_prefix + "(" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting RAW scores in AlgoMarker " + string(get_name()) + " With return code " + to_string(get_preds_rc);
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	if (am_matrix != "")
		ma.write_features_mat(am_matrix); // debug only

	// going over scores, and adding them to the right responses
	char **_score_types;
	int _n_score_types;
	responses->get_score_types(&_n_score_types, &_score_types);

	MedSamples *meds = ma.get_samples_ptr();

	for (auto &id_s : meds->idSamples) {
		for (auto &s : id_s.samples) {

			// basic info for current sample
			int c_pid = s.id;
			int c_ts = s.time;
			float c_scr = s.prediction[0];

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

						for (int j=0; j<_n_score_types; j++) {

							if (strcmp(_score_types[j], "Raw") == 0) {
								res->set_score(j, c_scr, _score_types[j]);
							}
							else {
								res->set_score(j, (float)AM_UNDEFINED_VALUE, _score_types[j]);
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

//-----------------------------------------------------------------------------------
// private internals for class MedialInfraAlgoMarker
//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::read_config(string conf_f)
{
	set_config(conf_f.c_str());

	ifstream inf(conf_f);

	if (!inf)
		return AM_ERROR_LOAD_NO_CONFIG_FILE;

	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size()-1] == '\r')
				curr_line.erase(curr_line.size()-1);

			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			if (fields.size() >= 2) {
				if (fields[0] == "TYPE") type_in_config_file = fields[1];
				else if (fields[0] == "REPOSITORY") rep_fname = fields[1];
				else if (fields[0] == "MODEL") model_fname = fields[1];
				else if (fields[0] == "INPUT_TESTER_CONFIG") input_tester_config_file = fields[1];
				else if (fields[0] == "NAME")  set_name(fields[1].c_str());
				else if (fields[0] == "TIME_UNIT") {
					set_time_unit(med_time_converter.string_to_type(fields[1].c_str()));
				}
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
