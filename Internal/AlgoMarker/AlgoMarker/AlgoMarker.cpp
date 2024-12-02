#include "AlgoMarker.h"

#include <Logger/Logger/Logger.h>
#include <MedTime/MedTime/MedTime.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//-----------------------------------------------------------------------------------
int AMPoint::auto_time_convert(long long ts, int to_type)
{
	long long date_t = 0;
	long long hhmm = 0;

	//MLOG("auto time convert: Date is %d , ts %lld , to_type %d\n", MedTime::Date, ts, to_type);

	if ((ts / (long long)1000000000) == 0) {
		date_t = ts; // yyyymmdd
		hhmm = 0;
	}
	else if (((ts / (long long)100000000000) == 0)) {
		date_t = ts / 100; // yyyymmddhh
		hhmm = 60 * (ts % 100);
	}
	else if (((ts / (long long)10000000000000) == 0)) {
		date_t = ts / 10000; // yyyymmddhhmm
		hhmm = 60 * ((ts % 10000) / 100) + (ts % 100);
	}
	else {
		date_t = ts / 1000000; // yyyymmddhhmmss
		hhmm = 60 * ((ts % 1000000) / 10000) + ((ts % 10000) / 100);
	}

	//MLOG("auto_time_converter: ts %lld to_type %d data_t %lld hhmm %lld\n", ts, to_type, date_t, hhmm);

	if (to_type == MedTime::Date) {
		//MLOG("auto time convert: date_t %d\n", date_t);
		//Ensure valid date:
		int year = int(date_t / 10000);
		if (year < 1900 || year > 3000) {
			//MTHROW_AND_ERR("Error invalid date %lld\n", ts);
			return -1;
		}
		return (int)date_t;
	}

	if (to_type == MedTime::Minutes) {
		int year = int(date_t / 10000);
		if (year < 1900 || year >2100) {
			//MTHROW_AND_ERR("Error invalid timestamp %lld\n", ts);
			return -1;
		}
		int minutes = med_time_converter.convert_date(MedTime::Minutes, (int)date_t);
		return minutes + (int)hhmm;
	}

	return 0;
}

//-----------------------------------------------------------------------------------
void AMMessages::get_messages(int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	if (need_to_update_args) {
		args.clear();
		for (auto &s : args_strs)
			args.push_back((char *)s.c_str());
		need_to_update_args = 0;
	}

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
	need_to_update_args = 1;
	//args.push_back((char *)args_strs.back().c_str()); 
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
	pair<int, long long> p(_pid, _timestamp);

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
	for (int i = 0; i < n_score_types; i++) {
		string s = string(_score_type[i]);
		score_types_str.push_back(s);
		stype2idx[s] = (int)score_types.size() - 1;
	}

	for (int i = 0; i < n_score_types; i++)
		score_types.push_back((char *)score_types_str[i].c_str());

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
	if (am_type == AM_TYPE_SIMPLE_EXAMPLE_EGFR)
		return new SimpleExampleEGFRAlgoMarker;

	return NULL;
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
int AM_API_Create(int am_type, AlgoMarker **new_am)
{
	try {
		*new_am = AlgoMarker::make_algomarker((AlgoMarkerType)am_type);

		if (new_am == NULL)
			return AM_FAIL_RC;

		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// loading AlgoMarker and making it ready to get Requests
//-----------------------------------------------------------------------------------------------------------
int AM_API_Load(AlgoMarker* pAlgoMarker, const char *config_fname)
{
	try {
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		return pAlgoMarker->Load(config_fname);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------
// Additional load options for AlgoMarker
//-----------------------------------------------------------------------------------------------------------
int AM_API_AdditionalLoad(AlgoMarker* pAlgoMarker, const int load_type, const char *load)
{
	try {
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		return pAlgoMarker->AdditionalLoad(load_type, load);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// clearing data from AlgoMarker (recommended at the start and/or end of each query session
//-----------------------------------------------------------------------------------------------------------
int AM_API_ClearData(AlgoMarker* pAlgoMarker)
{
	try {
		return pAlgoMarker->ClearData();
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// adding data to an AlgoMarker
// this API allows adding a specific signal, with matching arrays of times and values
//-----------------------------------------------------------------------------------------------------------
int AM_API_AddData(AlgoMarker* pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values)
{
	try {
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		return pAlgoMarker->AddData(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// adding data to an AlgoMarker using string parameters
// this API allows adding a specific signal, with matching arrays of times and string values which 
// also enables providing the categorial string representation of a value.
//-----------------------------------------------------------------------------------------------------------
int AM_API_AddDataStr(AlgoMarker* pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values)
{
	try {
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		return pAlgoMarker->AddDataStr(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Adding data by type
// DataType signals the format, data contains the actual data.
// Currently MedInfraAlgoMarker implements the DATA_JSON_FORMAT with data given as a json string
//-----------------------------------------------------------------------------------------------------------
int AM_API_AddDataByType(AlgoMarker* pAlgoMarker, const char *data, char **messages)
{
	try {
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		return pAlgoMarker->AddDataByType(data, messages);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Prepare a Request
// Null RC means failure
// pids and timestamps here are the timepoints to give predictions at
//-----------------------------------------------------------------------------------------------------------
int AM_API_CreateRequest(char *requestId, char **_score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req)
{
	try {
		(*new_req) = new AMRequest;

		if ((*new_req) == NULL)
			return AM_FAIL_RC;

		(*new_req)->set_request_id(requestId);
		(*new_req)->insert_score_types(_score_types, n_score_types);
		for (int i = 0; i < n_points; i++)
			(*new_req)->insert_point(patient_ids[i], time_stamps[i]);

		return AM_OK_RC;
	}
	catch (...) {
		(*new_req) = NULL;
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Create a new empty responses object to be later used 
//-----------------------------------------------------------------------------------------------------------
int AM_API_CreateResponses(AMResponses **new_responses)
{
	try {
		(*new_responses) = new AMResponses;

		if ((*new_responses) == NULL)
			return AM_FAIL_RC;

		return AM_OK_RC;
	}
	catch (...) {
		(*new_responses) = NULL;
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------
// Get scores for a ready request
//-----------------------------------------------------------------------------------------------------------
int AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses)
{
	try {
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		return pAlgoMarker->Calculate(request, responses);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Get general Calculate results 
//-----------------------------------------------------------------------------------------------------------
int AM_API_CalculateByType(AlgoMarker *pAlgoMarker, int CalcType, char *request, char **responses)
{
	try {
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		return pAlgoMarker->CalculateByType(CalcType, request, responses);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of AlgoMarker - free all memory 
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeAlgoMarker(AlgoMarker *pAlgoMarker)
{
	try {
		if (pAlgoMarker == NULL)
			return;

		pAlgoMarker->Unload();

		delete pAlgoMarker;
	}
	catch (...) {

	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of AMRequest - free all memory 
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeRequest(AMRequest *pRequest)
{
	try {
		if (pRequest == NULL)
			return;
		delete pRequest;
	}
	catch (...) {

	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of responses - free all memory
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeResponses(AMResponses *responses)
{
	try {
		if (responses == NULL)
			return;
		delete responses;
	}
	catch (...) {

	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of general allocated memory
//-----------------------------------------------------------------------------------------------------------
void AM_API_Dispose(char *data)
{
	try {
		if (data == NULL)
			return;
		delete[] data;
	}
	catch (...) {

	}
}
//-----------------------------------------------------------------------------------------------------------

void AM_API_Discovery(AlgoMarker *pAlgoMarker, char **resp) {
	*resp = NULL;
	try {
		if (pAlgoMarker == NULL)
			return;
		pAlgoMarker->Discovery(resp);
	}
	catch (...) {
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get number of responses (= no. of pid,time result points)
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponsesNum(AMResponses *responses)
{
	try {
		if (responses == NULL)
			return 0;
		return responses->get_n_responses();
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get shared msgs. Not a copy - direct pointers, so do not free.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetSharedMessages(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	try {
		if (responses == NULL)
			return AM_FAIL_RC;

		AMMessages *shared_m = responses->get_shared_messages();
		shared_m->get_messages(n_msgs, msgs_codes, msgs_args);

		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get an index of a specific pid,time response, or -1 if it doesn't exist
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseIndex(AMResponses *responses, int _pid, long long _timestamp)
{
	try {
		return responses->get_response_index_by_point(_pid, _timestamp);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get scores for a scpefic response given its index.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseAtIndex(AMResponses *responses, int res_index, AMResponse **res)
{
	try {
		*res = NULL;
		if (responses == NULL)
			return AM_FAIL_RC;

		if (res_index < 0 || res_index >= responses->get_n_responses())
			return AM_FAIL_RC;

		*res = responses->get_response(res_index);

		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get number of scores in a response (could contain several score types)
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseScoresNum(AMResponse *response, int *n_scores)
{
	try {
		if (response == NULL)
			return AM_FAIL_RC;

		*n_scores = response->get_n_scores();
		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// given a score index , return all we need about it : pid , timestamp, score and score type
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseScoreByIndex(AMResponse *response, int score_index, float *_score, char **_score_type)
{
	try {
		if (response == NULL)
			return AM_FAIL_RC;

		if (score_index < 0 || score_index >= response->get_n_scores())
			return AM_FAIL_RC;

		if (response->get_score(score_index, _score, _score_type) != AM_OK_RC)
			return AM_FAIL_RC;

		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------
// given a score index , return all we need about it : pid , timestamp, score and score type
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseExtendedScoreByIndex(AMResponse *response, int score_index, char **_score, char **_score_type)
{
	try {
		if (response == NULL)
			return AM_FAIL_RC;

		if (score_index < 0 || score_index >= response->get_n_scores())
			return AM_FAIL_RC;

		if (response->get_ext_score(score_index, _score, _score_type) != AM_OK_RC)
			return AM_FAIL_RC;

		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------
// get all messages for a specific response given its index
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseMessages(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	try {
		if (response == NULL)
			return AM_FAIL_RC;

		response->get_msgs()->get_messages(n_msgs, msgs_codes, msgs_args);
		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get all messages for a specific response given its index
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetScoreMessages(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	try {
		if (response == NULL)
			return AM_FAIL_RC;

		if (score_index < 0 || score_index >= response->get_n_scores())
			return AM_FAIL_RC;

		response->get_score_msgs(score_index)->get_messages(n_msgs, msgs_codes, msgs_args);
		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------
// get pid and timestamp of a response
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponsePoint(AMResponse *response, int *pid, long long *timestamp)
{
	try {
		if (response == NULL)
			return AM_FAIL_RC;

		*pid = response->get_patient_id();
		*timestamp = response->get_timestamp();
		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------
// get request id . Direct pointer so do not free.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponsesRequestId(AMResponses *responses, char **requestId)
{
	try {
		if (responses == NULL)
			return AM_FAIL_RC;

		*requestId = responses->get_request_id();
		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get a score using the response index and the score type. RC: fail if something is wrong.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseScoreByType(AMResponses *responses, int res_index, char *_score_type, float *out_score)
{
	try {
		if (responses == NULL)
			return AM_FAIL_RC;
		return responses->get_score_by_type(res_index, _score_type, out_score);
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get the nameof an algo marker
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetName(AlgoMarker *pAlgoMarker, char **name)
{
	try {
		*name = NULL;
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		*name = pAlgoMarker->get_name();
		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}
//-----------------------------------------------------------------------------------------------------------

int AM_API_DebugRepMemory(AlgoMarker *pAlgoMarker, char **resp)
{
	try {
		*resp = NULL;
		if (pAlgoMarker == NULL)
			return AM_FAIL_RC;

		pAlgoMarker->show_rep_data(resp);
		return AM_OK_RC;
	}
	catch (...) {
		return AM_FAIL_RC;
	}
}

//===========================================================================================================
// DATA API Implementation
//===========================================================================================================

//-----------------------------------------------------------------------------------------------------------
// create a Repository Handle and read data into memory, given file name, pids, and signals return: 0: OK -1: failed
//-----------------------------------------------------------------------------------------------------------
int DATA_API_RepositoryHandle_Create(RepositoryHandle **new_rep, char *fname, int *pids, int n_pids, char **sigs, int n_sigs)
{
	(*new_rep) = new RepositoryHandle;

	if ((*new_rep) == NULL) return -1;

	for (int i = 0; i < n_pids; i++)
		(*new_rep)->pids.push_back(pids[i]);

	for (int i = 0; i < n_sigs; i++)
		(*new_rep)->signals.push_back(string(sigs[i]));

	(*new_rep)->fname = string(fname);

	if ((*new_rep)->rep.read_all((*new_rep)->fname, (*new_rep)->pids, (*new_rep)->signals) < 0)
		return -1;

	return 0;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Create a SignalDataHandle : can be reused for each read later, create several if working in parallel
//-----------------------------------------------------------------------------------------------------------
int DATA_API_SignalDataHandle_Create(SignalDataHandle **new_sdh)
{
	(*new_sdh) = new SignalDataHandle;

	if ((*new_sdh) == NULL) return -1;

	return 0;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Read Data into a signal data handle from a repository, returns len=number of elements read
//-----------------------------------------------------------------------------------------------------------
int DATA_API_ReadData(RepositoryHandle *rep_h, int pid, char *sig, SignalDataHandle *sdh, int *len)
{
	rep_h->rep.uget(pid, string(sig), sdh->usv);
	(*len) = sdh->usv.len;

	return 0;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Get a time channel at an index from a loaded SignalDataHandle
//-----------------------------------------------------------------------------------------------------------
int DATA_API_GetTime(SignalDataHandle *sdh, int idx, int time_channel, int *time)
{
	(*time) = sdh->usv.Time(idx, time_channel);
	return 0;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Get a value channel at an index from a loaded SignalDataHandle
//-----------------------------------------------------------------------------------------------------------
int DATA_API_GetVal(SignalDataHandle *sdh, int idx, int val_channel, float *val)
{
	(*val) = sdh->usv.Val(idx, val_channel);
	return 0;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// dispose of SignalDataHandle
//-----------------------------------------------------------------------------------------------------------
void DATA_API_Dispose_SignalDataHandle(SignalDataHandle *sdh)
{
	if (sdh != NULL) delete sdh;
	sdh = NULL;
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// dispose of RepositoryHandle
//-----------------------------------------------------------------------------------------------------------
void DATA_API_Dispose_RepositoryHandle(RepositoryHandle *rep_h)
{
	if (rep_h != NULL) {
		rep_h->rep.clear();
		delete rep_h;
	}
	rep_h = NULL;
}
//-----------------------------------------------------------------------------------------------------------
