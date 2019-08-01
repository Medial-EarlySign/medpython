#ifndef __SODYNWRAPPER_H
#define __SODYNWRAPPER_H

#include "AlgoMarkerFlat.h"
#include <vector>
#include <assert.h>

class so_functions {
public:
	typedef int(*t_AM_API_Create)(int am_type, AlgoMarker **new_am);
	typedef int(*t_AM_API_Load)(AlgoMarker * pAlgoMarker, const char *config_fname);
	typedef int(*t_AM_API_ClearData)(AlgoMarker * pAlgoMarker);
	typedef int(*t_AM_API_AddData)(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
	typedef int(*t_AM_API_AddDataStr)(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values);
	typedef int(*t_AM_API_CreateRequest)(char *requestId, char **score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req);
	typedef int(*t_AM_API_CreateResponses)(AMResponses **);
	typedef int(*t_AM_API_Calculate)(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses);
	typedef int(*t_AM_API_GetSharedMessages)(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args);
	typedef int(*t_AM_API_GetResponsesNum)(AMResponses *responses);
	typedef int(*t_AM_API_GetResponseIndex)(AMResponses *responses, int _pid, long long _timestamp);
	typedef int(*t_AM_API_GetResponsesRequestId)(AMResponses *responses, char **requestId);
	typedef int(*t_AM_API_GetResponseScoreByType)(AMResponses *responses, int res_index, char *_score_type, float *out_score);
	typedef int(*t_AM_API_GetResponseAtIndex)(AMResponses *responses, int index, AMResponse **response);
	typedef int(*t_AM_API_GetResponseScoresNum)(AMResponse *response, int *n_scores);
	typedef int(*t_AM_API_GetResponseScoreByIndex)(AMResponse *response, int score_index, float *score, char **_score_type);
	typedef int(*t_AM_API_GetResponseMessages)(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args);
	typedef int(*t_AM_API_GetScoreMessages)(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args);
	typedef int(*t_AM_API_GetResponsePoint)(AMResponse *response, int *pid, long long *timestamp);
	typedef int(*t_AM_API_GetName)(AlgoMarker * pAlgoMArker, char **name);
	typedef void(*t_AM_API_DisposeAlgoMarker)(AlgoMarker*);
	typedef void(*t_AM_API_DisposeResponses)(AMResponses*);
	typedef void(*t_AM_API_DisposeRequest)(AMRequest*);
	void *AM_API_Create = nullptr;
	void *AM_API_Load = nullptr;
	void *AM_API_ClearData = nullptr;
	void *AM_API_AddData = nullptr;
	void *AM_API_AddDataStr = nullptr;
	void *AM_API_CreateRequest = nullptr;
	void *AM_API_CreateResponses = nullptr;
	void *AM_API_Calculate = nullptr;
	void *AM_API_GetSharedMessages = nullptr;
	void *AM_API_GetResponsesNum = nullptr;
	void *AM_API_GetResponseIndex = nullptr;
	void *AM_API_GetResponsesRequestId = nullptr;
	void *AM_API_GetResponseScoreByType = nullptr;
	void *AM_API_GetResponseAtIndex = nullptr;
	void *AM_API_GetResponseScoresNum = nullptr;
	void *AM_API_GetResponseScoreByIndex = nullptr;
	void *AM_API_GetResponseMessages = nullptr;
	void *AM_API_GetScoreMessages = nullptr;
	void *AM_API_GetResponsePoint = nullptr;
	void *AM_API_GetName = nullptr;
	void *AM_API_DisposeAlgoMarker = nullptr;
	void *AM_API_DisposeResponses = nullptr;
	void *AM_API_DisposeRequest = nullptr;
	// returns index in sos
	static int load(const char * am_fname);
	static so_functions* so;
	static std::vector<so_functions> sos;
	static void set_so_id(int id) { assert(id>=0 && id < sos.size()); so = &sos[id]; };
};



void load_am(const char * am_fname);

#endif //__SODYNWRAPPER_H
