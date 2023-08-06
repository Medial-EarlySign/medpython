#ifndef __SODYNWRAPPER_H
#define __SODYNWRAPPER_H

#include <vector>
#include <assert.h>

class AlgoMarker;
class AMRequest;
class AMResponse;
class AMResponses;

class DynAM {
public:
	typedef int(*t_AM_API_Create)(int am_type, AlgoMarker **new_am);
	typedef int(*t_AM_API_Load)(AlgoMarker * pAlgoMarker, const char *config_fname);
	typedef int(*t_AM_API_AdditionalLoad)(AlgoMarker * pAlgoMarker, const int load_type, const char *load);
	typedef int(*t_AM_API_ClearData)(AlgoMarker * pAlgoMarker);
	typedef int(*t_AM_API_AddData)(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
	typedef int(*t_AM_API_AddDataStr)(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values);
	typedef int(*t_AM_API_AddDataByType)(AlgoMarker * pAlgoMarker, const char *data, char ** messages);
	typedef int(*t_AM_API_CreateRequest)(char *requestId, char **score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req);
	typedef int(*t_AM_API_CreateResponses)(AMResponses **);
	typedef int(*t_AM_API_Calculate)(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses);
	typedef int(*t_AM_API_CalculateByType)(AlgoMarker *pAlgoMarker, int CalcType, char *request, char **responses);
	typedef int(*t_AM_API_GetSharedMessages)(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args);
	typedef int(*t_AM_API_GetResponsesNum)(AMResponses *responses);
	typedef int(*t_AM_API_GetResponseIndex)(AMResponses *responses, int _pid, long long _timestamp);
	typedef int(*t_AM_API_GetResponsesRequestId)(AMResponses *responses, char **requestId);
	typedef int(*t_AM_API_GetResponseScoreByType)(AMResponses *responses, int res_index, char *_score_type, float *out_score);
	typedef int(*t_AM_API_GetResponseAtIndex)(AMResponses *responses, int index, AMResponse **response);
	typedef int(*t_AM_API_GetResponseScoresNum)(AMResponse *response, int *n_scores);
	typedef int(*t_AM_API_GetResponseScoreByIndex)(AMResponse *response, int score_index, float *score, char **_score_type);
	typedef int(*t_AM_API_GetResponseExtendedScoreByIndex)(AMResponse *response, int score_index, char **ext_score, char **_score_type);
	typedef int(*t_AM_API_GetResponseMessages)(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args);
	typedef int(*t_AM_API_GetScoreMessages)(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args);
	typedef int(*t_AM_API_GetResponsePoint)(AMResponse *response, int *pid, long long *timestamp);
	typedef int(*t_AM_API_GetName)(AlgoMarker * pAlgoMArker, char **name);
	typedef void(*t_AM_API_DisposeAlgoMarker)(AlgoMarker*);
	typedef void(*t_AM_API_DisposeResponses)(AMResponses*);
	typedef void(*t_AM_API_DisposeRequest)(AMRequest*);
	typedef void(*t_AM_API_Discovery)(AlgoMarker*, char **);
	typedef void(*t_AM_API_Dispose)(char *);
	void *addr_AM_API_Create = nullptr;
	void *addr_AM_API_Load = nullptr;
	void *addr_AM_API_AdditionalLoad = nullptr;
	void *addr_AM_API_ClearData = nullptr;
	void *addr_AM_API_AddData = nullptr;
	void *addr_AM_API_AddDataStr = nullptr;
	void *addr_AM_API_AddDataByType = nullptr;
	void *addr_AM_API_CreateRequest = nullptr;
	void *addr_AM_API_CreateResponses = nullptr;
	void *addr_AM_API_Calculate = nullptr;
	void *addr_AM_API_CalculateByType = nullptr;
	void *addr_AM_API_GetSharedMessages = nullptr;
	void *addr_AM_API_GetResponsesNum = nullptr;
	void *addr_AM_API_GetResponseIndex = nullptr;
	void *addr_AM_API_GetResponsesRequestId = nullptr;
	void *addr_AM_API_GetResponseScoreByType = nullptr;
	void *addr_AM_API_GetResponseAtIndex = nullptr;
	void *addr_AM_API_GetResponseScoresNum = nullptr;
	void *addr_AM_API_GetResponseScoreByIndex = nullptr;
	void *addr_AM_API_GetResponseExtendedScoreByIndex = nullptr;
	void *addr_AM_API_GetResponseMessages = nullptr;
	void *addr_AM_API_GetScoreMessages = nullptr;
	void *addr_AM_API_GetResponsePoint = nullptr;
	void *addr_AM_API_GetName = nullptr;
	void *addr_AM_API_DisposeAlgoMarker = nullptr;
	void *addr_AM_API_DisposeResponses = nullptr;
	void *addr_AM_API_DisposeRequest = nullptr;
	void *addr_AM_API_Discovery = nullptr;
	void *addr_AM_API_Dispose = nullptr;
	// returns index in sos
	static int load(const char * am_fname);
	static DynAM* so;
	static std::vector<DynAM> sos;
	static void set_so_id(int id) { assert(id>=0 && id < sos.size()); so = &sos[id]; };

	static int AM_API_ClearData(AlgoMarker * pAlgoMarker);
	static void AM_API_DisposeAlgoMarker(AlgoMarker * pAlgoMarker);
	static void AM_API_DisposeRequest(AMRequest *pRequest);
	static void AM_API_Dispose(char *data);
	static void AM_API_Discovery(AlgoMarker * pAlgoMarker, char **resp);
	static void AM_API_DisposeResponses(AMResponses *responses);
	static int AM_API_GetResponseScoresNum(AMResponse *response, int *n_scores);
	static int AM_API_GetName(AlgoMarker * pAlgoMArker, char **name);
	static int AM_API_GetScoreMessages(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args);
	static int AM_API_GetResponseMessages(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args);
	static int AM_API_GetResponseScoreByType(AMResponses *responses, int res_index, char *_score_type, float *out_score);
	static int AM_API_GetResponseScoreByIndex(AMResponse *response, int score_index, float *score, char **_score_type);
	static int AM_API_GetResponseExtendedScoreByIndex(AMResponse *response, int score_index, char **ext_score, char **_score_type);
	static int AM_API_GetResponsePoint(AMResponse *response, int *pid, long long *timestamp);
	static int AM_API_GetSharedMessages(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args);
	static int AM_API_GetResponsesNum(AMResponses *responses);
	static int AM_API_GetResponseIndex(AMResponses *responses, int _pid, long long _timestamp);
	static int AM_API_GetResponseIndex(AMResponses *responses, char **requestId);
	static int AM_API_GetResponseAtIndex(AMResponses *responses, int index, AMResponse **response);
	static int AM_API_Create(int am_type, AlgoMarker **new_am);
	static int AM_API_Load(AlgoMarker * pAlgoMarker, const char *config_fname);
	static int AM_API_AdditionalLoad(AlgoMarker * pAlgoMarker, const int load_type, const char *load);
	static int AM_API_AddData(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
	static int AM_API_AddDataStr(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values);
	static int AM_API_AddDataByType(AlgoMarker * pAlgoMarker, const char *data, char ** messages);
	static int AM_API_CreateRequest(char *requestId, char **score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req);
	static int AM_API_CreateResponses(AMResponses **new_responses);
	static int AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses);
	static int AM_API_CalculateByType(AlgoMarker *pAlgoMarker, int CalcType, char *request, char **responses);

	static bool initialized() { return (sos.size() > 0); }
};



void load_am(const char * am_fname);

#endif //__SODYNWRAPPER_H
