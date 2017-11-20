#pragma once

//===============================================================================
// AlgoMarker.h
//-------------------------------------------------------------------------------
//
// Wrapper for a DLL that will contain :
// (1) All the needed API's to perform as an AlgoMarker
// (2) The option to run full MedProcessTools models with MedRepository inputs.
// (3) All in one single DLL managing it all.
//
//===============================================================================


// AM_DLL_EXPORT is defined only in the matching .cpp file to handle the dll building
// apps just include this h file and hence will work in import mode.


#ifdef _WIN32
#if defined AM_DLL_IMPORT
#define DLL_WORK_MODE __declspec(dllimport)
#else
#define DLL_WORK_MODE __declspec(dllexport)
#endif
#else
#define DLL_WORK_MODE
#endif

#ifdef _WIN32
#pragma warning(disable: 4251)
#endif

//
// includes of Medial Internal Libraries
//
#include "AlgoMarkerInternal.h"
#include "AlgoMarkerErr.h"


#define AM_UNDEFINED_VALUE -9999.99

typedef enum {
	AM_TYPE_UNDEFINED = 0,
	AM_TYPE_MEDIAL_INFRA = 1
} AlgoMarkerType;


//===============================================================================
// Responses and Requests classes
//===============================================================================
extern "C" class DLL_WORK_MODE AMPoint {
public:
	int pid = -1;
	long timestamp = -1;

	void set(int _pid, long _timestamp) { pid = _pid; timestamp = _timestamp; }
	
	void clear() { pid = -1; timestamp = -1; }
};

//-------------------------------------------------------------------------------
extern "C" class DLL_WORK_MODE AMMessages {
	
  private:

	vector<int> codes;
	vector<string> args_strs;

	vector<char *> args; // for easier c c# export. pointing to strings , so no need to free.

  public:

	// get things
	int get_n_msgs() { return (int)codes.size(); }
	void get_messages(int *n_msgs, int **msgs_codes, char ***msgs_args);

	// insert
	void insert_message(int code, const char *arg_ch);

	// clear
	void clear() { codes.clear(); args_strs.clear(); args.clear(); }

};

//-------------------------------------------------------------------------------
extern "C" class DLL_WORK_MODE AMResponse {

private:

	// p_score_types just points to the common info in the AMResponses class, no need to free 
	vector<char *> *p_score_types = NULL;

	// actual scores
	vector<float> scores;

	// In here we report messages not specific to a single Response
	AMMessages msgs;

	AMPoint point;

public:

	// get things
	int get_patient_id() { return point.pid; }
	long get_timestamp() { return point.timestamp; }
	void get_scores(float **_scores, char ***_types) { *_scores = &scores[0];  *_types = &((*p_score_types)[0]); }
	int get_n_scores() { return (int)scores.size(); }
	float get_score(int idx) { return scores[idx]; }
	AMMessages *get_msgs() { return &msgs; }

	// set things
	void set_patient_id(int _patient_id) { point.pid = _patient_id; }
	void set_timestamp(long _timestamp) { point.timestamp = _timestamp; }
	void set_score_types(vector<char *> *_types) { p_score_types = _types; }
	void set_score(int index, float score) { if (index < scores.size()) scores[index] = score; }
	void init_scores(int size) { scores.clear(); scores.resize(size, (float)AM_UNDEFINED_VALUE); }

	// clear
	void clear() { scores.clear(); msgs.clear(); p_score_types=NULL; point.clear(); }


};


//-------------------------------------------------------------------------------
extern "C" class DLL_WORK_MODE AMResponses {

private:

	string requestId = "";
	string version = "";

	// For each point: pid , time : we hold an AMResponse object that contains all the results on all types for this time point
	// plus all its specific messages
	vector<AMResponse> responses;
	map<pair<int,long>, int> point2response_idx;

	// score_types : these are common to all responses
	vector<string> score_types_str;
	vector<char *> score_types;
	unordered_map<string, int> stype2idx;

	// In here we report messages not specific to a single Response
	AMMessages shared_msgs;
public:

	// get things
	int get_n_responses() { return (int)responses.size(); }
	AMResponse *get_response(int index) { if (index >= (int)responses.size()) return NULL; return &(responses[index]); }
	int get_response_index_by_point(int _pid, long _timestamp); // if does not exist returns -1.
	AMResponse *get_response_by_point(int _pid, long _timestamp); // if does not exist, return NULL
	void get_score_types(int *n_score_types, char ***_score_types);
	AMMessages *get_shared_messages() { return &shared_msgs; }
	char *get_request_id() { return (char *)requestId.c_str(); }
	char *get_version() { return (char *)version.c_str(); }
	int get_score(int _pid, long _timestamp, char *_score_type, float *out_score);
	int get_score_by_type(int index, char *_score_type, float *out_score);
	vector<char *> *get_score_type_vec_ptr() { return &score_types; }

	// set things
	void set_request_id(char *request_id) {	requestId = string(request_id); }
	void set_version(char *_version) { version = string(_version); }
	void insert_score_types(char **_score_type, int n_score_types); 
	AMResponse *create_point_response(int _pid, long _timestamp);

	// clear
	void clear() { requestId=""; version=""; responses.clear(); point2response_idx.clear(); score_types_str.clear(); score_types.clear(); stype2idx.clear(); shared_msgs.clear(); }

};

//-------------------------------------------------------------------------------
extern "C" class DLL_WORK_MODE AMRequest {

private:

	// Id for tracking given by user
	string requestId = "";

	// score types asked for
	// Currently supporting : "Raw" 
	vector<string> score_types_str;

	// list of points to give scores at
	vector<AMPoint> points;

public:
	
	// get things
	char *get_request_id() { return (char *)requestId.c_str(); }
	int get_n_score_types() { return (int)score_types_str.size(); }
	char *get_score_type(int index) { if (index >= get_n_score_types()) return NULL; return (char *)score_types_str[index].c_str(); }
	int get_n_points() { return (int)points.size(); }
	AMPoint *get_point(int index) { if (index >= get_n_points()) return NULL;  return &points[index]; }
	int get_pid(int index) { if (index >= get_n_points()) return -1; return points[index].pid; }
	long get_timestamp(int index) { if (index >= get_n_points()) return -1; return points[index].timestamp; }

	// set things
	void set_request_id(char *req_id) { requestId = string(req_id); }
	void insert_point(int _pid, long _timestamp) { AMPoint p; p.set(_pid,_timestamp) ; points.push_back(p); }
	void insert_score_types(char **_score_types, int n_score_types) { for (int i=0; i<n_score_types; i++) score_types_str.push_back(string(_score_types[i])); }

	// clear
	void clear() { requestId=""; score_types_str.clear(); points.clear(); }

};


//===============================================================================
// Base AlgoMarker class
//===============================================================================
extern "C" class DLL_WORK_MODE AlgoMarker {
private:
	AlgoMarkerType type;
	string name = "";
	string config_fname = "";
	vector<string> supported_score_types;

public:

	// major APIs
	// When creating a new type of algomarker one needs to inherit from this class, and
	// make sure to implement the following virtual APIs. This will suffice.
	virtual int Load(const char *config_f) { return 0; }
	virtual int Unload() { return 0; }
	virtual int ClearData() { return 0; }
	virtual int AddData(int patient_id, const char *signalName,  int TimeStamps_len, long* TimeStamps, int Values_len, float* Values) { return 0; }
	virtual int Calculate(AMRequest *request, AMResponses **responses) { return 0; }

	// check supported score types in the supported_score_types vector
	int IsScoreTypeSupported(const char *_stype);

	// get things
	int get_type() { return (int)type; }
	char *get_name() { return  (char *)name.c_str(); }
	char *get_config() { return (char *)config_fname.c_str(); }

	// set things
	void set_type(int _type) { type = (AlgoMarkerType)_type; }
	void set_name(const char *_name) { name = string(_name); }
	void set_config(const char *_config_f) { config_fname = string(_config_f); }
	void add_supported_stype(const char *stype) { supported_score_types.push_back(string(stype)); }

	// get a new AlgoMarker
	static AlgoMarker *make_algomarker(AlgoMarkerType am_type);

};


//===============================================================================
// MedialInfraAlgoMarker - an AlgoMarker that works with Medial infrastructure
//===============================================================================
extern "C" class DLL_WORK_MODE MedialInfraAlgoMarker : public AlgoMarker {

private:
	MedAlgoMarkerInternal ma;

	// some configs
	string type_in_config_file = "";
	string rep_fname = "";
	string model_fname = "";
	//string input_tests_filters = "";

	int read_config(string conf_f);

	//vector<string> supported_score_types ={ "Raw" };


public:
	MedialInfraAlgoMarker() { set_type((int)AM_TYPE_MEDIAL_INFRA); add_supported_stype("Raw");}
	
	int Load(const char *config_f);
	int Unload();
	int ClearData();
	int AddData(int patient_id, const char *signalName, int TimeStamps_len, long* TimeStamps, int Values_len, float* Values);
	int Calculate(AMRequest *request, AMResponses **responses);

};

//========================================================================================
// Actual API to be used via C# :: External users should work ONLY with these !!
// Any change to the below functions must rely only on exposed API of the above classes.
//========================================================================================

// all return codes are defined in AlgoMarkerErr.h

// create a new AlgoMarker of type am_type and init its name
extern "C" DLL_WORK_MODE int AM_API_Create(int am_type, const char *name, AlgoMarker **new_am);

// loading AlgoMarker and making it ready to get Requests
extern "C" DLL_WORK_MODE int AM_API_Load(AlgoMarker* pAlgoMarker, const char *config_fname);

// clearing data from AlgoMarker (recommended at the start and/or end of each query session
extern "C" DLL_WORK_MODE int AM_API_ClearData(AlgoMarker* pAlgoMarker);

// adding data to an AlgoMarker
// this API allows adding a specific signal, with matching arrays of times and values
extern "C" DLL_WORK_MODE int AM_API_AddData(AlgoMarker* pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long* TimeStamps, int Values_len, float* Values);

// Prepare a Request
// Null RC means failure
extern "C" DLL_WORK_MODE int AM_API_CreateRequest(char *requestId, char **score_types, int n_score_types, int *patient_ids, long *time_stamps, int n_points, AMRequest **new_req);

// Get scores for a ready request
extern "C" DLL_WORK_MODE int AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses **responses);

// exploring responses
extern "C" DLL_WORK_MODE int AM_API_GetResponsesNum(AMResponses *responses);
extern "C" DLL_WORK_MODE int AM_API_GetSharedMessages(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args);
extern "C" DLL_WORK_MODE int AM_API_GetResponseIndex(AMResponses *responses, int _pid, long _timestamp);
extern "C" DLL_WORK_MODE int AM_API_GetResponse(AMResponses *responses, int index, int *pid, long *timestamp, int *n_scores, float **scores, char ***_score_types);
extern "C" DLL_WORK_MODE int AM_API_GetResponseMessages(AMResponses *responses, int index, int *n_msgs, int **msgs_codes, char ***msgs_args);
extern "C" DLL_WORK_MODE int AM_API_GetResponseRequestId(AMResponses *responses, char **requestId);
extern "C" DLL_WORK_MODE int AM_API_GetResponseScoreByType(AMResponses *responses, int res_index, char *_score_type, float *out_score);

// Dispose of AlgoMarker - free all memory 
extern "C" DLL_WORK_MODE void AM_API_DisposeAlgoMarker(AlgoMarker *pAlgoMarker);

// Dispose of AMRequest - free all memory 
extern "C" DLL_WORK_MODE void AM_API_DisposeRequest(AMRequest *pRequest);

// Dispose of responses - free all memory
extern "C" DLL_WORK_MODE void AM_API_DisposeResponses(AMResponses *responses);