#include <stdio.h>
#include <stdlib.h>


#ifdef __linux__ 
#include <dlfcn.h>
#elif _WIN32
#include <windows.h>
#endif // linux/win

#include "DynAMWrapper.h"

DynAM* DynAM::so = nullptr;
std::vector<DynAM> DynAM::sos;


void* load_sym(void* lib_h, const char* sym_name, bool exit_on_fail = true)
{
	printf("Loading %s ... ", sym_name);
#ifdef __linux__ 
	void* ret = dlsym(lib_h, sym_name);
	if (ret == nullptr) {
		char * err = (char*)dlerror();
		printf("Failed: %s\n", err);
#elif _WIN32
	void* ret = GetProcAddress((HMODULE)lib_h, sym_name);
	if (ret == nullptr) {
		printf("Failed\n");
#endif
		if (exit_on_fail)
			exit(0);
	}
	printf("OK\n");
	return ret;
}

void load_am(const char * am_fname) {
	if (DynAM::load(am_fname) < 0)
		exit(0);
}

int DynAM::load(const char * am_fname) {
	printf("Loading %s ... ", am_fname);
#ifdef __linux__ 
	void* lib_handle = dlopen(am_fname, RTLD_NOW); //RTLD_LAZY
#elif _WIN32
	void* lib_handle = (void*)LoadLibrary(am_fname);
#endif // linux/win


	if (lib_handle == NULL) {
#ifdef __linux__ 
		char * err = (char*)dlerror();
		if (err) printf("%s\n", err);
#elif _WIN32
		printf("Failed loading %s\n", am_fname);
#endif	
		return -1;
	}
	sos.push_back(DynAM());
	so = &sos.back();
	printf("OK\n");
	so->addr_AM_API_Create = load_sym(lib_handle, "AM_API_Create");
	so->addr_AM_API_Load = load_sym(lib_handle, "AM_API_Load");
	so->addr_AM_API_AdditionalLoad = load_sym(lib_handle, "AM_API_AdditionalLoad");
	so->addr_AM_API_ClearData = load_sym(lib_handle, "AM_API_ClearData");
	so->addr_AM_API_AddData = load_sym(lib_handle, "AM_API_AddData");
	so->addr_AM_API_AddDataStr = load_sym(lib_handle, "AM_API_AddDataStr", false);
	so->addr_AM_API_AddDataByType = load_sym(lib_handle, "AM_API_AddDataByType", false);
	so->addr_AM_API_CreateRequest = load_sym(lib_handle, "AM_API_CreateRequest");
	so->addr_AM_API_CreateResponses = load_sym(lib_handle, "AM_API_CreateResponses");
	so->addr_AM_API_Calculate = load_sym(lib_handle, "AM_API_Calculate");
	so->addr_AM_API_CalculateByType = load_sym(lib_handle, "AM_API_CalculateByType");
	so->addr_AM_API_GetResponsesNum = load_sym(lib_handle, "AM_API_GetResponsesNum");
	so->addr_AM_API_GetSharedMessages = load_sym(lib_handle, "AM_API_GetSharedMessages");
	so->addr_AM_API_GetResponseIndex = load_sym(lib_handle, "AM_API_GetResponseIndex");
	so->addr_AM_API_GetResponsesRequestId = load_sym(lib_handle, "AM_API_GetResponsesRequestId");
	so->addr_AM_API_GetResponseScoreByType = load_sym(lib_handle, "AM_API_GetResponseScoreByType");
	so->addr_AM_API_GetResponseAtIndex = load_sym(lib_handle, "AM_API_GetResponseAtIndex");
	so->addr_AM_API_GetResponseScoresNum = load_sym(lib_handle, "AM_API_GetResponseScoresNum");
	so->addr_AM_API_GetResponseScoreByIndex = load_sym(lib_handle, "AM_API_GetResponseScoreByIndex");
	so->addr_AM_API_GetResponseExtendedScoreByIndex = load_sym(lib_handle, "AM_API_GetResponseExtendedScoreByIndex");
	so->addr_AM_API_GetResponseMessages = load_sym(lib_handle, "AM_API_GetResponseMessages");
	so->addr_AM_API_GetScoreMessages = load_sym(lib_handle, "AM_API_GetScoreMessages");
	so->addr_AM_API_GetResponsePoint = load_sym(lib_handle, "AM_API_GetResponsePoint");
	so->addr_AM_API_GetName = load_sym(lib_handle, "AM_API_GetName");
	so->addr_AM_API_DisposeAlgoMarker = load_sym(lib_handle, "AM_API_DisposeAlgoMarker");
	so->addr_AM_API_DisposeRequest = load_sym(lib_handle, "AM_API_DisposeRequest");
	so->addr_AM_API_DisposeResponses = load_sym(lib_handle, "AM_API_DisposeResponses");
	so->addr_AM_API_Discovery = load_sym(lib_handle, "AM_API_Discovery", false);
	so->addr_AM_API_Dispose = load_sym(lib_handle, "AM_API_Dispose");
	return (int)sos.size() - 1;
}

int DynAM::AM_API_ClearData(AlgoMarker * pAlgoMarker) {
	return (*((DynAM::t_AM_API_ClearData)DynAM::so->addr_AM_API_ClearData))
		(pAlgoMarker);
}

void DynAM::AM_API_DisposeAlgoMarker(AlgoMarker * pAlgoMarker) {
	(*((DynAM::t_AM_API_DisposeAlgoMarker)DynAM::so->addr_AM_API_DisposeAlgoMarker))
		(pAlgoMarker);
}

void DynAM::AM_API_DisposeRequest(AMRequest *pRequest) {
	(*((DynAM::t_AM_API_DisposeRequest)DynAM::so->addr_AM_API_DisposeRequest))
		(pRequest);
}

void DynAM::AM_API_Dispose(char *data) {
	(*((DynAM::t_AM_API_Dispose)DynAM::so->addr_AM_API_Dispose))
		(data);
}

void DynAM::AM_API_DisposeResponses(AMResponses *responses) {
	(*((DynAM::t_AM_API_DisposeResponses)DynAM::so->addr_AM_API_DisposeResponses))
		(responses);
}

int DynAM::AM_API_GetResponseScoresNum(AMResponse *response, int *n_scores) {
	return (*((DynAM::t_AM_API_GetResponseScoresNum)DynAM::so->addr_AM_API_GetResponseScoresNum))
		(response, n_scores);
}
int DynAM::AM_API_GetName(AlgoMarker * pAlgoMArker, char **name) {
	return (*((DynAM::t_AM_API_GetName)DynAM::so->addr_AM_API_GetName))
		(pAlgoMArker, name);
}
int DynAM::AM_API_GetScoreMessages(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args) {
	return (*((DynAM::t_AM_API_GetScoreMessages)DynAM::so->addr_AM_API_GetScoreMessages))
		(response, score_index, n_msgs, msgs_codes, msgs_args);
}
int DynAM::AM_API_GetResponseMessages(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args) {
	return (*((DynAM::t_AM_API_GetResponseMessages)DynAM::so->addr_AM_API_GetResponseMessages))
		(response, n_msgs, msgs_codes, msgs_args);
}
int DynAM::AM_API_GetResponseScoreByType(AMResponses *responses, int res_index, char *_score_type, float *out_score) {
	return (*((DynAM::t_AM_API_GetResponseScoreByType)DynAM::so->addr_AM_API_GetResponseScoreByType))
		(responses, res_index, _score_type, out_score);
}
int DynAM::AM_API_GetResponseScoreByIndex(AMResponse *response, int score_index, float *score, char **_score_type) {
	return (*((DynAM::t_AM_API_GetResponseScoreByIndex)DynAM::so->addr_AM_API_GetResponseScoreByIndex))
		(response, score_index, score, _score_type);
}

int DynAM::AM_API_GetResponseExtendedScoreByIndex(AMResponse *response, int score_index, char **ext_score, char **_score_type) {
	return (*((DynAM::t_AM_API_GetResponseExtendedScoreByIndex)DynAM::so->addr_AM_API_GetResponseExtendedScoreByIndex))
		(response, score_index, ext_score, _score_type);
}

int DynAM::AM_API_GetResponsePoint(AMResponse *response, int *pid, long long *timestamp) {
	return (*((DynAM::t_AM_API_GetResponsePoint)DynAM::so->addr_AM_API_GetResponsePoint))
		(response, pid, timestamp);
}

int DynAM::AM_API_GetSharedMessages(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args) {
	return (*((DynAM::t_AM_API_GetSharedMessages)DynAM::so->addr_AM_API_GetSharedMessages))
		(responses, n_msgs, msgs_codes, msgs_args);
}

int DynAM::AM_API_GetResponsesNum(AMResponses *responses) {
	return (*((DynAM::t_AM_API_GetResponsesNum)DynAM::so->addr_AM_API_GetResponsesNum))
		(responses);
}
int DynAM::AM_API_GetResponseIndex(AMResponses *responses, int _pid, long long _timestamp) {
	return (*((DynAM::t_AM_API_GetResponseIndex)DynAM::so->addr_AM_API_GetResponseIndex))
		(responses, _pid, _timestamp);
}
int DynAM::AM_API_GetResponseIndex(AMResponses *responses, char **requestId) {
	return (*((DynAM::t_AM_API_GetResponsesRequestId)DynAM::so->addr_AM_API_GetResponsesRequestId))
		(responses, requestId);
}

int DynAM::AM_API_GetResponseAtIndex(AMResponses *responses, int index, AMResponse **response) {
	return (*((DynAM::t_AM_API_GetResponseAtIndex)DynAM::so->addr_AM_API_GetResponseAtIndex))
		(responses, index, response);
}

int DynAM::AM_API_Create(int am_type, AlgoMarker **new_am) {
	return (*((DynAM::t_AM_API_Create)DynAM::so->addr_AM_API_Create))
		(am_type, new_am);
}

int DynAM::AM_API_Load(AlgoMarker * pAlgoMarker, const char *config_fname) {
	return (*((DynAM::t_AM_API_Load)DynAM::so->addr_AM_API_Load))
		(pAlgoMarker, config_fname);
}

int DynAM::AM_API_AdditionalLoad(AlgoMarker * pAlgoMarker, const int load_type, const char *load) {
	return (*((DynAM::t_AM_API_AdditionalLoad)DynAM::so->addr_AM_API_AdditionalLoad))
		(pAlgoMarker, load_type, load);
}

int DynAM::AM_API_AddData(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values) {
	return (*((DynAM::t_AM_API_AddData)DynAM::so->addr_AM_API_AddData))
		(pAlgoMarker, patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values);
}

int DynAM::AM_API_AddDataStr(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values) {
	return (*((DynAM::t_AM_API_AddDataStr)DynAM::so->addr_AM_API_AddDataStr))
		(pAlgoMarker, patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values);
}

int DynAM::AM_API_AddDataByType(AlgoMarker * pAlgoMarker, const char *data, char ** messages) {
	return (*((DynAM::t_AM_API_AddDataByType)DynAM::so->addr_AM_API_AddDataByType))
		(pAlgoMarker, data, messages);
}


int DynAM::AM_API_CreateRequest(char *requestId, char **score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req) {
	return (*((DynAM::t_AM_API_CreateRequest)DynAM::so->addr_AM_API_CreateRequest))
		(requestId, score_types, n_score_types, patient_ids, time_stamps, n_points, new_req);
}

int DynAM::AM_API_CreateResponses(AMResponses **new_responses) {
	return (*((DynAM::t_AM_API_CreateResponses)DynAM::so->addr_AM_API_CreateResponses))
		(new_responses);
}

int DynAM::AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses) {
	return (*((DynAM::t_AM_API_Calculate)DynAM::so->addr_AM_API_Calculate))
		(pAlgoMarker, request, responses);
}

int DynAM::AM_API_CalculateByType(AlgoMarker *pAlgoMarker, int CalcType, char *request, char **responses) {
	return (*((DynAM::t_AM_API_CalculateByType)DynAM::so->addr_AM_API_CalculateByType))
		(pAlgoMarker, CalcType, request, responses);
}

void DynAM::AM_API_Discovery(AlgoMarker *pAlgoMarker, char **resp) {
	if (DynAM::so->addr_AM_API_Discovery != NULL)
		(*((DynAM::t_AM_API_Discovery)DynAM::so->addr_AM_API_Discovery))
		(pAlgoMarker, resp);
	else
		*resp = NULL;
}