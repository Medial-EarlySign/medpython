#include <stdio.h>
#include <stdlib.h>


#ifdef __linux__ 
#include <dlfcn.h>
#elif _WIN32
#include <windows.h>
#endif // linux/win

#include "DynAMWrapper.h"
//#include "AlgoMarkerFlat.h"

class so_functions {
public:
  typedef int (*t_AM_API_Create)(int am_type, AlgoMarker **new_am);
  typedef int (*t_AM_API_Load)(AlgoMarker * pAlgoMarker, const char *config_fname);
  typedef int (*t_AM_API_ClearData)(AlgoMarker * pAlgoMarker);
  typedef int (*t_AM_API_AddData)(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
  typedef int (*t_AM_API_CreateRequest)(char *requestId, char **score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req);
  typedef int (*t_AM_API_CreateResponses)(AMResponses **);
  typedef int (*t_AM_API_Calculate)(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses);
  typedef int (*t_AM_API_GetSharedMessages)(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args);
  typedef int (*t_AM_API_GetResponsesNum)(AMResponses *responses);
  typedef int (*t_AM_API_GetResponseIndex)(AMResponses *responses, int _pid, long long _timestamp);
  typedef int (*t_AM_API_GetResponsesRequestId)(AMResponses *responses, char **requestId);
  typedef int (*t_AM_API_GetResponseScoreByType)(AMResponses *responses, int res_index, char *_score_type, float *out_score);
  typedef int (*t_AM_API_GetResponseAtIndex)(AMResponses *responses, int index, AMResponse **response);
  typedef int (*t_AM_API_GetResponseScoresNum)(AMResponse *response, int *n_scores);
  typedef int (*t_AM_API_GetResponseScoreByIndex)(AMResponse *response, int score_index, float *score, char **_score_type);
  typedef int (*t_AM_API_GetResponseMessages)(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args);
  typedef int (*t_AM_API_GetScoreMessages)(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args);
  typedef int (*t_AM_API_GetResponsePoint)(AMResponse *response, int *pid, long long *timestamp);
  typedef int (*t_AM_API_GetName)(AlgoMarker * pAlgoMArker, char **name);
  typedef void (*t_AM_API_DisposeAlgoMarker)(AlgoMarker*);
  typedef void (*t_AM_API_DisposeResponses)(AMResponses*);
  typedef void (*t_AM_API_DisposeRequest)(AMRequest*);
  void *AM_API_Create=nullptr;
  void *AM_API_Load=nullptr;
  void *AM_API_ClearData=nullptr;
  void *AM_API_AddData=nullptr;
  void *AM_API_CreateRequest=nullptr;
  void *AM_API_CreateResponses=nullptr;
  void *AM_API_Calculate=nullptr;
  void *AM_API_GetSharedMessages=nullptr;
  void *AM_API_GetResponsesNum=nullptr;
  void *AM_API_GetResponseIndex=nullptr;
  void *AM_API_GetResponsesRequestId=nullptr;
  void *AM_API_GetResponseScoreByType=nullptr;
  void *AM_API_GetResponseAtIndex=nullptr;
  void *AM_API_GetResponseScoresNum=nullptr;
  void *AM_API_GetResponseScoreByIndex=nullptr;
  void *AM_API_GetResponseMessages=nullptr;
  void *AM_API_GetScoreMessages=nullptr;
  void *AM_API_GetResponsePoint=nullptr;
  void *AM_API_GetName=nullptr;
  void *AM_API_DisposeAlgoMarker=nullptr;
  void *AM_API_DisposeResponses=nullptr;
  void *AM_API_DisposeRequest=nullptr;
};


so_functions so;

void* load_sym(void* lib_h, const char* sym_name)
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
    exit(0);
  }
  printf("OK\n");
  return ret;
}

void load_am(const char * am_fname){
  printf("Loading %s ... ",am_fname);
#ifdef __linux__ 
  void* lib_handle = dlopen(am_fname, RTLD_NOW); //RTLD_LAZY
#elif _WIN32
  void* lib_handle = (void*)LoadLibrary(am_fname);
#endif // linux/win


  if(lib_handle == NULL){
#ifdef __linux__ 
	char * err = (char*)dlerror();
    if(err) printf("%s\n",err);
#elif _WIN32
	  printf("Failed loading %s\n", am_fname);
#endif	
    exit(0);
  }
  printf("OK\n");
  so.AM_API_Create = load_sym(lib_handle, "AM_API_Create");
  so.AM_API_Load = load_sym(lib_handle, "AM_API_Load");
  so.AM_API_ClearData = load_sym(lib_handle, "AM_API_ClearData");
  so.AM_API_AddData = load_sym(lib_handle, "AM_API_AddData");
  so.AM_API_CreateRequest = load_sym(lib_handle, "AM_API_CreateRequest");
  so.AM_API_CreateResponses = load_sym(lib_handle, "AM_API_CreateResponses");
  so.AM_API_Calculate = load_sym(lib_handle, "AM_API_Calculate");
  so.AM_API_GetResponsesNum = load_sym(lib_handle, "AM_API_GetResponsesNum");
  so.AM_API_GetSharedMessages = load_sym(lib_handle, "AM_API_GetSharedMessages");
  so.AM_API_GetResponseIndex = load_sym(lib_handle, "AM_API_GetResponseIndex");
  so.AM_API_GetResponsesRequestId = load_sym(lib_handle, "AM_API_GetResponsesRequestId");
  so.AM_API_GetResponseScoreByType = load_sym(lib_handle, "AM_API_GetResponseScoreByType");
  so.AM_API_GetResponseAtIndex = load_sym(lib_handle, "AM_API_GetResponseAtIndex");
  so.AM_API_GetResponseScoresNum = load_sym(lib_handle, "AM_API_GetResponseScoresNum");
  so.AM_API_GetResponseScoreByIndex = load_sym(lib_handle, "AM_API_GetResponseScoreByIndex");
  so.AM_API_GetResponseMessages = load_sym(lib_handle, "AM_API_GetResponseMessages");
  so.AM_API_GetScoreMessages = load_sym(lib_handle, "AM_API_GetScoreMessages");
  so.AM_API_GetResponsePoint = load_sym(lib_handle, "AM_API_GetResponsePoint");
  so.AM_API_GetName = load_sym(lib_handle, "AM_API_GetName");
  so.AM_API_DisposeAlgoMarker = load_sym(lib_handle, "AM_API_DisposeAlgoMarker");
  so.AM_API_DisposeRequest = load_sym(lib_handle, "AM_API_DisposeRequest");
  so.AM_API_DisposeResponses = load_sym(lib_handle, "AM_API_DisposeResponses");
}

int AM_API_ClearData(AlgoMarker * pAlgoMarker){
  return (*((so_functions::t_AM_API_ClearData)so.AM_API_ClearData))
    (pAlgoMarker);
}

void AM_API_DisposeAlgoMarker(AlgoMarker * pAlgoMarker){
  (*((so_functions::t_AM_API_DisposeAlgoMarker)so.AM_API_DisposeAlgoMarker))
    (pAlgoMarker);
}

void AM_API_DisposeRequest(AMRequest *pRequest){
  (*((so_functions::t_AM_API_DisposeRequest)so.AM_API_DisposeRequest))
    (pRequest);
}

void AM_API_DisposeResponses(AMResponses *responses){
  (*((so_functions::t_AM_API_DisposeResponses)so.AM_API_DisposeResponses))
    (responses);
}

int AM_API_GetResponseScoresNum(AMResponse *response, int *n_scores){
  return (*((so_functions::t_AM_API_GetResponseScoresNum)so.AM_API_GetResponseScoresNum))
    (response, n_scores);
}
int AM_API_GetName(AlgoMarker * pAlgoMArker, char **name){
  return (*((so_functions::t_AM_API_GetName)so.AM_API_GetName))
    (pAlgoMArker, name);
}
int AM_API_GetScoreMessages(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args){
  return (*((so_functions::t_AM_API_GetScoreMessages)so.AM_API_GetScoreMessages))
    (response, score_index, n_msgs, msgs_codes, msgs_args);
}
int AM_API_GetResponseMessages(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args){
  return (*((so_functions::t_AM_API_GetResponseMessages)so.AM_API_GetResponseMessages))
    (response, n_msgs, msgs_codes, msgs_args);
}
int AM_API_GetResponseScoreByType(AMResponses *responses,int res_index, char *_score_type, float *out_score){
  return (*((so_functions::t_AM_API_GetResponseScoreByType)so.AM_API_GetResponseScoreByType))
    (responses, res_index, _score_type, out_score);
}
int AM_API_GetResponseScoreByIndex(AMResponse *response, int score_index, float *score, char **_score_type){
  return (*((so_functions::t_AM_API_GetResponseScoreByIndex)so.AM_API_GetResponseScoreByIndex))
    (response,score_index,score,_score_type);
}

int AM_API_GetResponsePoint(AMResponse *response, int *pid, long long *timestamp){
  return (*((so_functions::t_AM_API_GetResponsePoint)so.AM_API_GetResponsePoint))
    (response,pid,timestamp);
}

int AM_API_GetSharedMessages(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args){
  return (*((so_functions::t_AM_API_GetSharedMessages)so.AM_API_GetSharedMessages))
    (responses,n_msgs,msgs_codes,msgs_args);
}

int AM_API_GetResponsesNum(AMResponses *responses){
   return (*((so_functions::t_AM_API_GetResponsesNum)so.AM_API_GetResponsesNum))
     (responses);
}
int AM_API_GetResponseIndex(AMResponses *responses, int _pid, long long _timestamp){
   return (*((so_functions::t_AM_API_GetResponseIndex)so.AM_API_GetResponseIndex))
     (responses,_pid,_timestamp);
}
int AM_API_GetResponseIndex(AMResponses *responses, char **requestId){
   return (*((so_functions::t_AM_API_GetResponsesRequestId)so.AM_API_GetResponsesRequestId))
     (responses, requestId);
}

int AM_API_GetResponseAtIndex(AMResponses *responses, int index, AMResponse **response){
   return (*((so_functions::t_AM_API_GetResponseAtIndex)so.AM_API_GetResponseAtIndex))
     (responses,index,response);
}

int AM_API_Create(int am_type, AlgoMarker **new_am){
  return (*((so_functions::t_AM_API_Create)so.AM_API_Create))
    (am_type, new_am);
}

int AM_API_Load(AlgoMarker * pAlgoMarker, const char *config_fname){
  return (*((so_functions::t_AM_API_Load)so.AM_API_Load))
    (pAlgoMarker, config_fname);
}

int AM_API_AddData(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values){
  return (*((so_functions::t_AM_API_AddData)so.AM_API_AddData))
    (pAlgoMarker, patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values);
}

int AM_API_CreateRequest(char *requestId, char **score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req){
  return (*((so_functions::t_AM_API_CreateRequest)so.AM_API_CreateRequest))
    (requestId,score_types,n_score_types,patient_ids,time_stamps,n_points,new_req);
}

int AM_API_CreateResponses(AMResponses **new_responses){
  return (*((so_functions::t_AM_API_CreateResponses)so.AM_API_CreateResponses))
    (new_responses);
}

int AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses){
  return (*((so_functions::t_AM_API_Calculate)so.AM_API_Calculate))
    (pAlgoMarker, request, responses);
}
