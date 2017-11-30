#include "AlgoMarker.h"

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
//-----------------------------------------------------------------------------------
void AMMessages::get_messages(int *n_msgs, int **msgs_codes, char ***msgs_args) 
{
	*n_msgs = get_n_msgs();
	if (*n_msgs > 0) {
		*msgs_codes = &codes[0];
		*msgs_args = &args[0];
	}
	else {
		*msgs_codes = NULL;
		*msgs_args = NULL;
	}
}

//-----------------------------------------------------------------------------------
void AMMessages::insert_message(int code, const char *arg_ch)
{
	string arg = string(arg_ch);
	codes.push_back(code); 
	args_strs.push_back(arg); 
	args.push_back((char *)args_strs.back().c_str()); 
}

//-----------------------------------------------------------------------------------
// if does not exist returns -1.
int AMResponses::get_response_index_by_point(int _pid, long long _timestamp)
{
	pair<int, long long> p(_pid, _timestamp);

	if (point2response_idx.find(p) == point2response_idx.end())
		return -1;
	
	return point2response_idx[p];

}

//-----------------------------------------------------------------------------------
// if does not exist returns NULL
AMResponse *AMResponses::get_response_by_point(int _pid, long long _timestamp)
{
	pair<int,long long> p(_pid,_timestamp);

	if (point2response_idx.find(p) == point2response_idx.end())
		return NULL;

	return &responses[point2response_idx[p]];
}

//-----------------------------------------------------------------------------------
void AMResponses::get_score_types(int *n_score_types, char ***_score_types) 
{
	*n_score_types = (int)score_types.size(); 
	if (n_score_types == 0) 
		*_score_types = NULL; 
	else 
		*_score_types = &score_types[0]; 
}

//-----------------------------------------------------------------------------------
int AMResponses::get_score(int _pid, long long _timestamp, char *_score_type, float *out_score)
{
	pair<int, long long> p(_pid, _timestamp);

	if (point2response_idx.find(p) == point2response_idx.end())
		return AM_FAIL_RC;
	int pidx = point2response_idx[p];

	return get_score_by_type(pidx, _score_type, out_score);
}

//-----------------------------------------------------------------------------------
int AMResponses::get_score_by_type(int index, char *_score_type, float *out_score)
{
	string s = string(_score_type);

	if (index < 0 || index >= get_n_responses())
		return AM_FAIL_RC;
	if (stype2idx.find(s) == stype2idx.end())
		return AM_FAIL_RC;
	int sidx = stype2idx[s];
	char *dummy_type;
	if (responses[index].get_score(sidx, out_score, &dummy_type) != AM_OK_RC) return AM_FAIL_RC;
	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
void AMResponses::insert_score_types(char **_score_type, int n_score_types) {
	for (int i=0; i<n_score_types; i++) {
		string s = string(_score_type[i]);
		score_types_str.push_back(s);
		score_types.push_back((char *)score_types_str.back().c_str());
		stype2idx[s] = (int)score_types.size() - 1;
	}
}

//-----------------------------------------------------------------------------------
AMResponse *AMResponses::create_point_response(int _pid, long long _timestamp)
{
	pair<int, long long> p(_pid, _timestamp);

	AMResponse response;

	response.set_patient_id(_pid);
	response.set_timestamp(_timestamp);
	response.init_scores((int)score_types.size());

	responses.push_back(response);
	
	point2response_idx[p] = (int)responses.size() - 1;

	return &responses.back();
}

//-----------------------------------------------------------------------------------
int AlgoMarker::IsScoreTypeSupported(const char *_stype)
{
	string stype = string(_stype);

	for (auto &s : supported_score_types)
		if (stype == s)
			return 1;
	return 0;
}

//-----------------------------------------------------------------------------------
AlgoMarker *AlgoMarker::make_algomarker(AlgoMarkerType am_type)
{
	if (am_type == AM_TYPE_MEDIAL_INFRA)
		return new MedialInfraAlgoMarker;

	return NULL;
}

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


	// prepare internal ma for work: set name, rep and model
	ma.set_name(get_name());
	
	if (ma.init_rep_config(rep_fname.c_str()) < 0)
		return AM_ERROR_LOAD_READ_REP_ERR;

	if (ma.init_model_from_file(model_fname.c_str()) < 0)
		return AM_ERROR_LOAD_READ_MODEL_ERR;


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

	if (TimeStamps_len > 0) {
		times_int.resize(TimeStamps_len);
		for (int i=0; i<TimeStamps_len; i++)
			times_int[i] = (int)TimeStamps[i];
		i_times = &times_int[0];
	}

	if (ma.data_load_pid_sig(patient_id, signalName, i_times, TimeStamps_len, Values, Values_len) < 0)
		return AM_ERROR_ADD_DATA_FAILED;

	return AM_OK_RC;
}

//------------------------------------------------------------------------------------------
// Calculate() - after data loading : get a request, get predictions, and pack as responses
//------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Calculate(AMRequest *request, AMResponses **responses)
{
	*responses = new AMResponses; // allocating responses, should be disposed by user after usage.

	(*responses)->set_request_id(request->get_request_id());
	for (int i=0; i<request->get_n_score_types(); i++) {
		char *stype = request->get_score_type(i);
		(*responses)->insert_score_types(&stype, 1);
	}

	string msg_prefix = "reqId: " + string(request->get_request_id()) + " :: ";

	AMMessages *shared_msgs = (*responses)->get_shared_messages();

	if (request == NULL) {
		string msg = msg_prefix + "(" + to_string(AM_MSG_NULL_REQUEST) + " ) NULL request in Calculate()";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}

	// We now have to prepare samples for the requested points
	// again - we only deal with int times in this class, so we convert the long long stamps to int
	ma.clear_samples();
	int n_points = request->get_n_points();

	for (int i=0; i<n_points; i++)
		if (ma.insert_sample(request->get_pid(i), (int)request->get_timestamp(i)) < 0) {
			string msg = msg_prefix + "(" + to_string(AM_MSG_BAD_PREDICTION_POINT) + ") Failed insert prediction point " + to_string(i) + " pid: " + to_string(request->get_pid(i)) + " ts: " + to_string(request->get_timestamp(i));
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	ma.normalize_samples();

	// Checking score types and verify they are supported
	int n_score_types = request->get_n_score_types();
	for (int i=0; i<n_score_types; i++) {
		if (!IsScoreTypeSupported(request->get_score_type(i))) {
			string msg = msg_prefix + "(" + to_string(AM_MSG_BAD_SCORE_TYPE) + ") AlgoMarker of type " + string(get_name()) + " does not support score type " + string(request->get_score_type(i));
			shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
			return AM_FAIL_RC;
		}
	}

	// Calculating raw scores
	vector<int> _pids(n_points, -1), _times(n_points, -1);
	vector<float> raw_scores(n_points, (float)AM_UNDEFINED_VALUE);

	int get_preds_rc;
	if ((get_preds_rc = ma.get_raw_preds(&_pids[0], &_times[0], &raw_scores[0])) < 0) {
		string msg = msg_prefix + "(" + to_string(AM_MSG_RAW_SCORES_ERROR) + ") Failed getting RAW scores in AlgoMarker " + string(get_name()) + " With return code " + to_string(get_preds_rc);
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return AM_FAIL_RC;
	}


	// Reporting back and building responses
	// going over raw scores, and for each create a response
	char **_score_types;
	int _n_score_types;
	(*responses)->get_score_types(&_n_score_types, &_score_types);
	for (int i=0; i<n_points; i++) {

		// create a response
		AMResponse *res = (*responses)->create_point_response(_pids[i], (long long)_times[i]);

		//res->set_score_types((*responses)->get_score_type_vec_ptr());
		res->init_scores(_n_score_types);

		for (int j=0; j<_n_score_types; j++) {

			if (strcmp(_score_types[j], "Raw") == 0) {
				res->set_score(j, raw_scores[i], _score_types[j]);
			}
			else
				res->set_score(j, (float)AM_UNDEFINED_VALUE, _score_types[j]);

		}

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

	return AM_OK_RC;
}

//===========================================================================================================
//===========================================================================================================
//===========================================================================================================
// A P I   I M P L E M E N T A T I O N S
//===========================================================================================================
//===========================================================================================================
//===========================================================================================================

//-----------------------------------------------------------------------------------------------------------
// create a new AlgoMarker of type am_type and init its name
//-----------------------------------------------------------------------------------------------------------
int AM_API_Create(int am_type, const char *name, AlgoMarker **new_am)
{
	*new_am = AlgoMarker::make_algomarker((AlgoMarkerType)am_type);

	(*new_am)->set_name(name);

	if (new_am == NULL)
		return AM_FAIL_RC;

	return AM_OK_RC;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// loading AlgoMarker and making it ready to get Requests
//-----------------------------------------------------------------------------------------------------------
int AM_API_Load(AlgoMarker* pAlgoMarker, const char *config_fname)
{
	if (pAlgoMarker == NULL)
		return AM_FAIL_RC;

	return pAlgoMarker->Load(config_fname);
}
//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------
// clearing data from AlgoMarker (recommended at the start and/or end of each query session
//-----------------------------------------------------------------------------------------------------------
int AM_API_ClearData(AlgoMarker* pAlgoMarker)
{
	return pAlgoMarker->ClearData();
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// adding data to an AlgoMarker
// this API allows adding a specific signal, with matching arrays of times and values
//-----------------------------------------------------------------------------------------------------------
int AM_API_AddData(AlgoMarker* pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values)
{
	if (pAlgoMarker == NULL)
		return AM_FAIL_RC;

	return pAlgoMarker->AddData(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values);
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Prepare a Request
// Null RC means failure
// pids and timestamps here are the timepoints to give predictions at
//-----------------------------------------------------------------------------------------------------------
int AM_API_CreateRequest(char *requestId, char **_score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req)
{
	(*new_req) = new AMRequest;

	if ((*new_req) == NULL)
		return AM_FAIL_RC;
	
	(*new_req)->set_request_id(requestId);
	(*new_req)->insert_score_types(_score_types, n_score_types);
	for (int i=0; i<n_points; i++)
		(*new_req)->insert_point(patient_ids[i], time_stamps[i]);

	return AM_OK_RC;
}
//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------
// Get scores for a ready request
//-----------------------------------------------------------------------------------------------------------
int AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses **responses)
{
	if (pAlgoMarker == NULL || request == NULL)
		return AM_FAIL_RC;

	return pAlgoMarker->Calculate(request, responses);
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of AlgoMarker - free all memory 
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeAlgoMarker(AlgoMarker *pAlgoMarker)
{
	if (pAlgoMarker == NULL)
		return;

	pAlgoMarker->Unload();

	delete pAlgoMarker;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of AMRequest - free all memory 
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeRequest(AMRequest *pRequest)
{
	if (pRequest == NULL)
		return;
	delete pRequest;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of responses - free all memory
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeResponses(AMResponses *responses)
{
	if (responses == NULL)
		return;
	delete responses;
}
//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------
// get number of responses (= no. of pid,time result points)
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponsesNum(AMResponses *responses)
{
	if (responses == NULL)
		return 0;
	return responses->get_n_responses();
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get shared msgs. Not a copy - direct pointers, so do not free.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetSharedMessages(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	if (responses == NULL)
		return AM_FAIL_RC;

	AMMessages *shared_m = responses->get_shared_messages();
	shared_m->get_messages(n_msgs, msgs_codes, msgs_args);

	return AM_OK_RC;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get an index of a specific pid,time response, or -1 if it doesn't exist
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseIndex(AMResponses *responses, int _pid, long long _timestamp)
{
	return responses->get_response_index_by_point(_pid, _timestamp);
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get scores for a scpefic response given its index.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponse(AMResponses *responses, int res_index, AMResponse **res)
{
	*res = NULL;
	if (responses == NULL)
		return AM_FAIL_RC;

	if (res_index < 0 || res_index >= responses->get_n_responses())
		return AM_FAIL_RC;

	*res = responses->get_response(res_index);

	return AM_OK_RC;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get number of scores in a response (could contain several score types)
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseScoresNum(AMResponse *response, int *n_scores)
{
	if (response == NULL)
		return AM_FAIL_RC;

	*n_scores = response->get_n_scores();
	return AM_FAIL_RC;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// given a score index , return all we need about it : pid , timestamp, score and score type
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseScoreByIndex(AMResponse *response, int score_index, int *pid, long long *timestamp, float *_score, char **_score_type)
{
	if (response == NULL)
		return AM_FAIL_RC;

	if (score_index < 0 || score_index >= response->get_n_scores())
		return AM_FAIL_RC;

	*pid = response->get_patient_id();
	*timestamp = response->get_timestamp();
	if (response->get_score(score_index, _score, _score_type) != AM_OK_RC)
		return AM_FAIL_RC;

	return AM_OK_RC;	
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get all messages for a specific response given its index
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseMessages(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	if (response == NULL)
		return AM_FAIL_RC;

	if (score_index < 0 || score_index >= response->get_n_scores())
		return AM_FAIL_RC;

	response->get_msgs(score_index)->get_messages(n_msgs, msgs_codes, msgs_args);
	return AM_OK_RC;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get request id . Direct pointer so do not free.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseRequestId(AMResponses *responses, char **requestId)
{
	if (responses == NULL)
		return AM_FAIL_RC;

	*requestId = responses->get_request_id();
	return AM_OK_RC;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get a score using the response index and the score type. RC: fail if something is wrong.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseScoreByType(AMResponses *responses, int res_index, char *_score_type, float *out_score)
{
	if (responses == NULL)
		return AM_FAIL_RC;
	return responses->get_score_by_type(res_index, _score_type, out_score);
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get the nameof an algo marker
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetName(AlgoMarker *pAlgoMarker, char **name)
{
	*name = NULL;
	if (pAlgoMarker == NULL) 
		return AM_FAIL_RC;

	*name = pAlgoMarker->get_name();
	return AM_OK_RC;
}
//-----------------------------------------------------------------------------------------------------------

