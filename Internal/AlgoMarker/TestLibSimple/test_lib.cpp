#include <cstdio>
#include <cstring>
#include <dlfcn.h>
#include <string>
#include <memory>

#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include <stdexcept>
#include <fstream>

typedef enum {
	AM_TYPE_UNDEFINED = 0,
	AM_TYPE_MEDIAL_INFRA = 1,
	AM_TYPE_SIMPLE_EXAMPLE_EGFR = 2,
} AlgoMarkerType;

#define AM_UNDEFINED_VALUE -9999.99
#define DATA_BATCH_JSON_FORMAT		2002
#define JSON_REQ_JSON_RESP			3001

#define AM_OK_RC									0

// General FAIL RC
#define AM_FAIL_RC									-1

using namespace std;

class AMPoint {
public:
	int pid = -1;
	long long timestamp = -1;

	void set(int _pid, long long _timestamp) { pid = _pid; timestamp = _timestamp; }

	void clear() { pid = -1; timestamp = -1; }
};

class AMRequest {

private:

	// Id for tracking given by user
	string requestId = "";

	// score types asked for
	// Currently supporting : "Raw" 
	vector<string> score_types_str;

	// list of points to give scores at
	vector<AMPoint> points;

public:

	vector<string> verificationConfig;


	// get things
	char *get_request_id() { return (char *)requestId.c_str(); }
	int get_n_score_types() { return (int)score_types_str.size(); }
	char *get_score_type(int index) { if (index >= get_n_score_types()) return NULL; return (char *)score_types_str[index].c_str(); }
	int get_n_points() { return (int)points.size(); }
	AMPoint *get_point(int index) { if (index >= get_n_points()) return NULL;  return &points[index]; }
	int get_pid(int index) { if (index >= get_n_points()) return -1; return points[index].pid; }
	long long get_timestamp(int index) { if (index >= get_n_points()) return -1; return points[index].timestamp; }

	// set things
	void set_request_id(char *req_id) { requestId = string(req_id); }
	void insert_point(int _pid, long long _timestamp) { AMPoint p; p.set(_pid, _timestamp); points.push_back(p); }
	void insert_score_types(char **_score_types, int n_score_types) { for (int i = 0; i<n_score_types; i++) score_types_str.push_back(string(_score_types[i])); }

	// clear
	void clear()
	{
		for (auto element : points)
		{
			element.clear();
		}
		requestId.clear();
		score_types_str.clear();
		points.clear();
	}
};

class AMMessages {

private:

	vector<int> codes;
	vector<string> args_strs;

	vector<char *> args; // for easier c c# export. pointing to strings , so no need to free.
	int need_to_update_args = 0;

public:

	// get things
	int get_n_msgs() { return (int)codes.size(); }
	void get_messages(int *n_msgs, int **msgs_codes, char ***msgs_args);

	// insert
	void insert_message(int code, const char *arg_ch);

	// clear
	void clear() {
		codes.clear();
		args_strs.clear();
		args.clear();
	}

};


class AMScore {
private:
	// no need to release the score type pointer
	char *p_score_type = NULL;
	float score = (float)AM_UNDEFINED_VALUE;
	AMMessages msgs;

public:
	// get things
	void get_score(float *_score, char **_score_type)
	{
		*_score = score;
		*_score_type = p_score_type;

		//char* x = p_score_type;

		/*	if (p_score_type != NULL)
		{
		string x = string(p_score_type);
		*_score_type = (char *)x.c_str();
		}*/
	}
	AMMessages *get_msgs() { return &msgs; }

	// set things
	void set_score_type(char *_score_type) { p_score_type = _score_type; }
	void set_score(float _score) { score = _score; }

	// clear
	void clear()
	{
		msgs.clear();
		p_score_type = NULL;
		score = (float)AM_UNDEFINED_VALUE;
	}
};


class AMResponse {

private:

	// p_score_types just points to the common info in the AMResponses class, no need to free 
	vector<AMScore> scores;

	AMPoint point;

	AMMessages msgs;

public:

	vector<string> verificationConfig;
	map<string, int> scoresIndexs;
	int need_to_update_scoreTypes = 0;

	// get things
	int get_patient_id() { return point.pid; }
	long long get_timestamp() { return point.timestamp; }
	int get_n_scores() { return (int)scores.size(); }
	AMScore *get_am_score(int idx) { if (idx < 0 || idx >= (int)scores.size()) return NULL; return &scores[idx]; }
	int get_score(int idx, float *_score, char **_score_type) {
		if (idx < 0 || idx >= (int)scores.size()) return AM_FAIL_RC;
		scores[idx].get_score(_score, _score_type);
		return AM_OK_RC;
	}
	AMMessages *get_score_msgs(int idx) { if (idx < 0 || idx >= (int)scores.size()) return NULL; return scores[idx].get_msgs(); }
	AMMessages *get_msgs() { return &msgs; }

	// set things
	void set_patient_id(int _patient_id) { point.pid = _patient_id; }
	void set_timestamp(long long _timestamp) { point.timestamp = _timestamp; }
	void set_score(int idx, float _score, char *_score_type) { if (idx >= 0 && idx < (int)scores.size()) scores[idx].set_score(_score); scores[idx].set_score_type(_score_type); }
	void init_scores(int size) { scores.clear(); scores.resize(size); }

	// clear
	void clear()
	{
		for (auto element : scores)
		{
			element.clear();
		}
		scores.clear();
		point.clear();
	}


};

class AMResponses;

class AMResponses {

private:

	string requestId = "";
	string version = "";

	// For each point: pid , time : we hold an AMResponse object that contains all the results on all types for this time point
	// plus all its specific messages
	vector<AMResponse> responses;
	//map<pair<int, long long>, int> point2response_idx;

	// score_types : these are common to all responses
	vector<string> score_types_str;
	vector<char *> score_types;
	//unordered_map<string, int> stype2idx;

	// In here we report messages not specific to a single Response
	AMMessages shared_msgs;
public:

	vector<string> verificationConfig;

	// get things
	int get_n_responses()
	{
		return (int)responses.size();
	}
	AMResponse *get_response(int index) { if (index >= (int)responses.size()) return NULL; return &(responses[index]); }
	int get_response_index_by_point(int _pid, long long _timestamp); // if does not exist returns -1.
																	 //	AMResponse *get_response_by_point(int _pid, long long _timestamp); // if does not exist, return NULL
	void get_score_types(int *n_score_types, char ***_score_types);
	AMMessages *get_shared_messages() { return &shared_msgs; }
	char *get_request_id() { return (char *)requestId.c_str(); }
	char *get_version() { return (char *)version.c_str(); }
	int get_score(int _pid, long long _timestamp, char *_score_type, float *out_score);
	int get_score_by_type(int index, char *_score_type, float *out_score);
	vector<char *> *get_score_type_vec_ptr() { return &score_types; }

	// set things
	void set_request_id(char *request_id) { requestId = string(request_id); }
	void set_version(char *_version) { version = string(_version); }
	void insert_score_types(char **_score_type, int n_score_types);
	AMResponse *create_point_response(int _pid, long long _timestamp);

	// clear
	void clear() {
		requestId.clear();
		version.clear();
		for (auto element : responses)
		{
			element.clear();
		}
		responses.clear();
		//point2response_idx.clear(); 
		score_types_str.clear();
		score_types.clear();
		//stype2idx.clear(); 
		shared_msgs.clear();
	}

};


class AlgoMarker {
private:
	AlgoMarkerType type;
	string name = "";
	string am_udi_di = "";
	string am_version = "";
	string config_fname = "";
	vector<string> supported_score_types;
	int time_unit = 1; // typically Date (for outpatient) or Minutes (for in patients)

public:

	// major APIs
	// When creating a new type of algomarker one needs to inherit from this class, and
	// make sure to implement the following virtual APIs. This will suffice.
	virtual int Load(const char *config_f) { return 0; }
	virtual int Unload() { return 0; }
	virtual int ClearData() { return 0; }
	virtual int AddData(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values) { return 0; }
	virtual int AddDataStr(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values) { return 0; }
	virtual int Calculate(AMRequest *request, AMResponses *responses) { return 0; }

	// Extentions
	virtual int AdditionalLoad(const int LoadType, const char *load) { return 0; } // options for LoadType: LOAD_DICT_FROM_FILE , LOAD_DICT_FROM_JSON
	virtual int AddDataByType(const char *data, char **messages) { return 0; } // options: DATA_JSON_FORMAT
	virtual int CalculateByType(int CalculateType, char *request, char **response) { return 0; } // options: JSON_REQ_JSON_RESP


																								 // check supported score types in the supported_score_types vector
	int IsScoreTypeSupported(const char *_stype);

	// get things
	int get_type() { return (int)type; }
	char *get_name() { return  (char *)name.c_str(); }
	char *get_config() { return (char *)config_fname.c_str(); }
	int get_time_unit() { return time_unit; }
	char *get_am_udi_di() { return  (char *)am_udi_di.c_str(); }
	char *get_am_version() { return  (char *)am_version.c_str(); }

	// set things
	void set_type(int _type) { type = (AlgoMarkerType)_type; }
	void set_name(const char *_name) { name = string(_name); }
	void set_config(const char *_config_f) { config_fname = string(_config_f); }
	void add_supported_stype(const char *stype) { supported_score_types.push_back(string(stype)); }
	void set_time_unit(int tu) { time_unit = tu; }
	void set_am_udi_di(const char *_am_udi_di) { am_udi_di = string(_am_udi_di); }
	void set_am_version(const char *_am_version) { am_version = string(_am_version); }

	// get a new AlgoMarker
	static AlgoMarker *make_algomarker(AlgoMarkerType am_type);

	virtual ~AlgoMarker() { ClearData(); Unload(); };

};


class DynAM {
public:
	typedef int(*t_AM_API_Create)(int am_type, AlgoMarker **new_am);
	typedef int(*t_AM_API_Load)(AlgoMarker * pAlgoMarker, const char *config_fname);
	typedef int(*t_AM_API_AdditionalLoad)(AlgoMarker * pAlgoMarker, const int load_type, const char *load);
	typedef int(*t_AM_API_ClearData)(AlgoMarker * pAlgoMarker);
	typedef int(*t_AM_API_AddData)(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
	typedef int(*t_AM_API_AddDataStr)(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values);
	typedef int(*t_AM_API_AddDataByType)(AlgoMarker * pAlgoMarker, const char *data, char **messages);
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
	void *addr_AM_API_Dispose = nullptr;
	// returns index in sos
	static int load(const char * am_fname);
	static DynAM* so;
	static std::vector<DynAM> sos;
	static void set_so_id(int id) { assert(id >= 0 && id < (int)sos.size()); so = &sos[id]; };

	static int AM_API_ClearData(AlgoMarker * pAlgoMarker);
	static void AM_API_DisposeAlgoMarker(AlgoMarker * pAlgoMarker);
	static void AM_API_DisposeRequest(AMRequest *pRequest);
	static void AM_API_Dispose(char *data);
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
	static int AM_API_AddDataByType(AlgoMarker * pAlgoMarker, const char *data, char **messages);
	static int AM_API_CreateRequest(char *requestId, char **score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req);
	static int AM_API_CreateResponses(AMResponses **new_responses);
	static int AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses);
	static int AM_API_CalculateByType(AlgoMarker *pAlgoMarker, int CalcType, char *request, char **responses);

	static bool initialized() { return (sos.size() > 0); }
};

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
			return NULL;
	}
	printf("OK\n");
	return ret;
	}

void load_am(const char * am_fname) {
	if (DynAM::load(am_fname)<0)
		printf("Error\n");
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
	if (DynAM::so->addr_AM_API_Load == NULL)
		printf("AM_API_Load is NULL\n");
	else
		printf("running AM_API_Load\n");
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

int DynAM::AM_API_AddDataByType(AlgoMarker * pAlgoMarker, const char *data, char **messages) {
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


void initialize_algomarker(const char *amconfig, AlgoMarker *&test_am)
{
	// Load
	printf("Loading AM\n");
	int rc = DynAM::AM_API_Load(test_am, amconfig);
	printf("Loaded\n");
	if (rc != AM_OK_RC) {
		printf("ERROR: Failed loading algomarker %s with config file %s ERR_CODE: %d\n", test_am->get_name(), amconfig, rc);
	}
	printf("Name is %s\n", test_am->get_name());
}

int read_file_into_string(const char *fname, string &data)
{
	ifstream inf(fname);
	if (!inf) {
		printf("MedUtils:MedIO :: read_file_inot_string: Can't open file %s\n", fname);
		return -1;
	}

	inf.seekg(0, std::ios::end);
	size_t size = inf.tellg();
	data.resize(size);
	inf.seekg(0);
	inf.read(&data[0], size);
	return 0;
}


void init_and_load_data(const char *input_json_path, AlgoMarker *am) {
	DynAM::AM_API_ClearData(am);

	string in_jsons;
	char * out_messages;
	if (read_file_into_string(input_json_path, in_jsons) < 0) {
		printf("Error on loading file %s\n", in_jsons.c_str());
		throw logic_error("Error");
	}
	printf("read %zu characters from input jsons file %s\n", in_jsons.length(), input_json_path);
	int load_status = DynAM::AM_API_AddDataByType(am, in_jsons.c_str(), &out_messages);
	if (out_messages != NULL) {
		string msgs = string(out_messages); //New line for each message:
		printf("AddDataByType has messages:\n");
		printf("%s\n", msgs.c_str());
	}
	printf("Added data from %s\n", input_json_path);
	if (load_status != AM_OK_RC)
		printf("Error code returned from calling AddDataByType: %d\n", load_status);
}

int get_preds_from_algomarker_single(AlgoMarker *am,
	const string &sjreq, bool calc_by_type,
	int pred_id, int pred_time, bool print_resp = false)
{

	const char * stypes[] = { "Raw" };
	int pid_id = pred_id;
	long long _timestamp = pred_time;
	char *jreq = (char *)(sjreq.c_str());
	char *jresp;

	if (calc_by_type) {
		DynAM::AM_API_CalculateByType(am, JSON_REQ_JSON_RESP, jreq, &jresp);
		if (print_resp) {
			string s = string(jresp);
			printf("%s\n", s.c_str());
		}
		DynAM::AM_API_Dispose(jresp);
	}
	else {
		AMResponses *resp;
		DynAM::AM_API_CreateResponses(&resp);

		AMRequest *req;
		int req_create_rc = DynAM::AM_API_CreateRequest("test_request", stypes, 1, &pid_id, &_timestamp, 1, &req);
		if (req_create_rc > 0) {
			printf("Faield to create req\n");
			throw logic_error("Error");
		}

		DynAM::AM_API_Calculate(am, req, resp); // calculate

		DynAM::AM_API_GetResponsesNum(resp);

		AMResponse *response;
		int resp_rc = DynAM::AM_API_GetResponseAtIndex(resp, 0, &response);
		//get_response_score_into_sample(response, resp_rc);
		if (resp_rc <0) {
			printf("error in score calc\n");
		}

		if (print_resp) {
			//string s = string(response);
			//MLOG("%s\n", s.c_str());
		}

		DynAM::AM_API_DisposeRequest(req);
		DynAM::AM_API_DisposeResponses(resp);
	}

	return 0;
}

int main(int argc, char *argv[]) {

	if (argc <= 2) {
		printf("Please pass path to lib + amconfig + (optional data_json) (optional data_output)\n");
		return -1;
	}
	char *am_fname = argv[1];
	char *amconfig = argv[2];

	int pid_id = 1;
	int prediction_time = 20210101;

	printf("Loading %s ... ", am_fname);
	load_am(am_fname);

	printf("Creating AM\n");

	AlgoMarker *test_am;
	if (DynAM::AM_API_Create((int)AM_TYPE_MEDIAL_INFRA, &test_am) != AM_OK_RC) {
		printf("ERROR: Failed creating test algomarker\n");
		return -1;
	}
	printf("Created!\n");

	initialize_algomarker(amconfig, test_am);

	string sjreq = "";
	sjreq = "{ \"type\" : \"request\", \"request_id\" : \"my test\", \"export\" : {\"prediction\" : \"pred_0\"}, \"requests\" : [{ \"patient_id\": \"" + to_string(pid_id) +
		"\", \"time\" : \"" + to_string(prediction_time) + "\" }] }";

	if (argc > 3) {
		char *data_json_path = argv[3];
		init_and_load_data(data_json_path, test_am);
	}

	if (argc > 4) {//3 and up
		bool calc_type = stoi(string(argv[4]))>0;

		get_preds_from_algomarker_single(test_am, sjreq, calc_type, pid_id, prediction_time, true);
	}

	printf("Clear data!\n");
	DynAM::AM_API_ClearData(test_am);

	printf("Disposing!\n");
	DynAM::AM_API_DisposeAlgoMarker(test_am);
	printf("Done all!\n");

	return 0;
}

//g++ -Wall --std=c++11 -ldl -march=x86-64 -msse2 -msse3 -msse4 test.cpp -o test_lib
//sudo docker cp ./test_lib 31dbefe0000f:/work