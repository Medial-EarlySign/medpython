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

// Branched to a separate branch (test)

// AM_DLL_EXPORT is defined only in the matching .cpp file to handle the dll building
// apps just include this h file and hence will work in import mode.

// set next to:
// 1 : and then use Release mode - this is to compile without DLLs (needed for direct tests of performance , etc)
// 0 : and then use ReleaseDLL mode - this is to compile with DLLs
#if 0
#define DLL_WORK_MODE
#else
#ifdef _WIN32
#if defined AM_DLL_IMPORT
#define DLL_WORK_MODE __declspec(dllimport)
#else    // !AM_DLL_IMPORT
#define DLL_WORK_MODE __declspec(dllexport)
#endif   // AM_DLL_IMPORT
#else    // !_WIN32
#if defined AM_DLL_IMPORT
#define DLL_WORK_MODE
#else    // !AM_DLL_IMPORT
#define DLL_WORK_MODE __attribute__ ((visibility ("default")))
#endif   //AM_DLL_IMPORT
#endif   // _WIN32
#endif // 0/1

#ifdef _WIN32
#pragma warning(disable: 4251)
#endif

//
// includes of Medial Internal Libraries
//
#ifndef ALGOMARKER_FLAT_API
#include "AlgoMarkerInternal.h"
#endif // ALGOMARKER_FLAT_API
#include "AlgoMarkerErr.h"


typedef enum {
	AM_TYPE_UNDEFINED = 0,
	AM_TYPE_MEDIAL_INFRA = 1,
	AM_TYPE_SIMPLE_EXAMPLE_EGFR = 2,
} AlgoMarkerType;


// Extension modes options
#define LOAD_DICT_FROM_JSON	1001
#define LOAD_DICT_FROM_FILE	1002

#define DATA_JSON_FORMAT			2001
#define DATA_BATCH_JSON_FORMAT		2002

#define JSON_REQ_JSON_RESP			3001

#ifndef ALGOMARKER_FLAT_API

//===============================================================================
// Responses and Requests classes
//===============================================================================
extern "C" class DLL_WORK_MODE AMPoint {
public:
	int pid = -1;
	long long timestamp = -1;

	void set(int _pid, long long _timestamp) { pid = _pid; timestamp = _timestamp; }

	void clear() { pid = -1; timestamp = -1; }

	// auto time convertor helper function
	static int auto_time_convert(long long ts, int to_format);
};

//-------------------------------------------------------------------------------
extern "C" class DLL_WORK_MODE AMMessages {

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
	void clear() { codes.clear(); args_strs.clear(); args.clear(); }

};

//-------------------------------------------------------------------------------
extern "C" class DLL_WORK_MODE AMScore {
private:
	// no need to release the score type pointer
	char *p_score_type = NULL;
	float score = (float)AM_UNDEFINED_VALUE;
	AMMessages msgs;
	string extended_score;
public:
	// get things
	void get_score(float *_score, char **_score_type) { *_score = score; *_score_type = p_score_type; }
	void get_extended_score(char** _ext_score, char **_score_type) { *_ext_score = &extended_score[0]; *_score_type = p_score_type; }
	AMMessages *get_msgs() { return &msgs; }

	// set things
	void set_score_type(char *_score_type) { p_score_type = _score_type; }
	void set_score(float _score) { score = _score; }
	void set_ext_score(const string& _ext__score) { extended_score = _ext__score; }

	// clear
	void clear() { msgs.clear(); p_score_type = NULL; score = (float)AM_UNDEFINED_VALUE; }
};

//-------------------------------------------------------------------------------
extern "C" class DLL_WORK_MODE AMResponse {

private:

	// p_score_types just points to the common info in the AMResponses class, no need to free 
	vector<AMScore> scores;

	AMPoint point;

	AMMessages msgs;

public:

	// get things
	int get_patient_id() { return point.pid; }
	long long get_timestamp() { return point.timestamp; }
	int get_n_scores() { return (int)scores.size(); }
	AMScore *get_am_score(int idx) { if (idx < 0 || idx >= scores.size()) return NULL; return &scores[idx]; }
	int get_score(int idx, float *_score, char **_score_type) {
		if (idx < 0 || idx >= scores.size()) return AM_FAIL_RC;
		scores[idx].get_score(_score, _score_type);
		return AM_OK_RC;
	}
	int get_ext_score(int idx, char **_ext_score, char **_score_type) {
		if (idx < 0 || idx >= scores.size()) return AM_FAIL_RC;
		scores[idx].get_extended_score(_ext_score, _score_type);
		return AM_OK_RC;
	}
	AMMessages *get_score_msgs(int idx) { if (idx < 0 || idx >= scores.size()) return NULL; return scores[idx].get_msgs(); }
	AMMessages *get_msgs() { return &msgs; }

	// set things
	void set_patient_id(int _patient_id) { point.pid = _patient_id; }
	void set_timestamp(long long _timestamp) { point.timestamp = _timestamp; }
	void set_score(int idx, float _score, char *_score_type, const string& _ext_score) {
		if (idx >= 0 && idx < scores.size()) scores[idx].set_score(_score);
		scores[idx].set_score_type(_score_type);
		scores[idx].set_ext_score(_ext_score);
	}
	void init_scores(int size) { scores.clear(); scores.resize(size); }

	// clear
	void clear() { scores.clear(); point.clear(); }


};

//-------------------------------------------------------------------------------
extern "C" class DLL_WORK_MODE AMResponses {

private:

	string requestId = "";
	string version = "";

	// For each point: pid , time : we hold an AMResponse object that contains all the results on all types for this time point
	// plus all its specific messages
	vector<AMResponse> responses;
	map<pair<int, long long>, int> point2response_idx;

	// score_types : these are common to all responses
	vector<string> score_types_str;
	vector<char *> score_types;
	unordered_map<string, int> stype2idx;

	// In here we report messages not specific to a single Response
	AMMessages shared_msgs;
public:

	AMResponses() { clear(); }
	~AMResponses() { clear(); }

	// get things
	int get_n_responses() { return (int)responses.size(); }
	AMResponse *get_response(int index) { if (index >= (int)responses.size()) return NULL; return &(responses[index]); }
	int get_response_index_by_point(int _pid, long long _timestamp); // if does not exist returns -1.
	AMResponse *get_response_by_point(int _pid, long long _timestamp); // if does not exist, return NULL
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
	void clear() { requestId = ""; version = ""; responses.clear(); point2response_idx.clear(); score_types_str.clear(); score_types.clear(); stype2idx.clear(); shared_msgs.clear(); }

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
	long long get_timestamp(int index) { if (index >= get_n_points()) return -1; return points[index].timestamp; }

	// set things
	void set_request_id(char *req_id) { requestId = string(req_id); }
	void insert_point(int _pid, long long _timestamp) { AMPoint p; p.set(_pid, _timestamp); points.push_back(p); }
	void insert_score_types(char **_score_types, int n_score_types) { for (int i = 0; i < n_score_types; i++) score_types_str.push_back(string(_score_types[i])); }

	// clear
	void clear() { requestId = ""; score_types_str.clear(); points.clear(); }

};


//===============================================================================
// Base AlgoMarker class
//===============================================================================
extern "C" class DLL_WORK_MODE AlgoMarker {
private:
	AlgoMarkerType type;
	string name = "";
	string am_udi_di = "";
	string am_manfactor_date = "";
	string am_version = "";
	string config_fname = "";
	vector<string> supported_score_types;
	int time_unit = MedTime::Date; // typically Date (for outpatient) or Minutes (for in patients)

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
	virtual int AddDataByType(const char *data, char **messages) { *messages = NULL; return 0; }
	virtual int CalculateByType(int CalculateType, char *request, char **response) { *response = NULL; return 0; } // options: JSON_REQ_JSON_RESP

	//Discovery api
	virtual int Discovery(char **response) { *response = NULL; return 0; }

	// check supported score types in the supported_score_types vector
	int IsScoreTypeSupported(const char *_stype);

	// get things
	int get_type() { return (int)type; }
	char *get_name() { return  (char *)name.c_str(); }
	char *get_config() { return (char *)config_fname.c_str(); }
	int get_time_unit() { return time_unit; }
	char *get_am_udi_di() { return  (char *)am_udi_di.c_str(); }
	char *get_manfactor_date() { return (char *)am_manfactor_date.c_str(); }
	char *get_am_version() { return  (char *)am_version.c_str(); }
	virtual void show_rep_data(char **response) { *response = NULL; }

	// set things
	void set_type(int _type) { type = (AlgoMarkerType)_type; }
	void set_name(const char *_name) { name = string(_name); }
	void set_config(const char *_config_f) { config_fname = string(_config_f); }
	void add_supported_stype(const char *stype) { supported_score_types.push_back(string(stype)); }
	void set_time_unit(int tu) { time_unit = tu; }
	void set_am_udi_di(const char *_am_udi_di) { am_udi_di = string(_am_udi_di); }
	void set_manafactur_date(const char *_am_man_date) { am_manfactor_date = string(_am_man_date); }
	void set_am_version(const char *_am_version) { am_version = string(_am_version); }

	// get a new AlgoMarker
	static AlgoMarker *make_algomarker(AlgoMarkerType am_type);

	virtual ~AlgoMarker() { ClearData(); Unload(); };

};


//===============================================================================
// MedialInfraAlgoMarker - an AlgoMarker that works with Medial infrastructure
//===============================================================================
extern "C" class DLL_WORK_MODE MedialInfraAlgoMarker : public AlgoMarker {

private:
	MedAlgoMarkerInternal ma;
	InputSanityTester ist;

	// some configs
	string type_in_config_file = "";
	string rep_fname = "";
	string model_fname = "";
	string input_tester_config_file = "";
	bool allow_rep_adjustments = false;

	int read_config(const string &conf_f);

	//vector<string> supported_score_types ={ "Raw" };

	int sort_needed = 1; // in some debug cases we ommit the sort od data at the end of loading to do that this needs to be 0
	string am_matrix = ""; // for debugging : if not empty will write matrix to given file name
	bool first_write = true; ///< in debug mode - mark first write flag
	int model_end_stage = MED_MDL_END;
	vector<string> extended_result_fields;
	bool is_loaded = false;

	void get_jsons_locations(const char *data, vector<size_t> &j_start, vector<size_t> &j_len); // helper to split given string to jsons within it. Used in batch json mode.
	int AddJsonData(int patient_id, json &j_data, vector<string> &messages, map<pair<int, int>, pair<int, vector<char>>> *data = NULL);
	int rec_AddDataByType(int DataType, const char *data, vector<string> &messages);
	void clear_patients_data(const vector<int> &pids);
	int AddDataStr_data(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values, 
	map<pair<int, int>, pair<int, vector<char>>> *data);
	int AddData_data(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values, 
	map<pair<int, int>, pair<int, vector<char>>> *data);
public:
	MedialInfraAlgoMarker() { set_type((int)AM_TYPE_MEDIAL_INFRA); add_supported_stype("Raw"); }

	int Load(const char *config_f);
	int Unload();
	int ClearData();
	int AddData(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
	int AddDataStr(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values);
	int Calculate(AMRequest *request, AMResponses *responses);
	int AdditionalLoad(const int LoadType, const char *load); // options for LoadType: LOAD_DICT_FROM_FILE , LOAD_DICT_FROM_JSON
	int AddDataByType(const char *data, char **messages);
	int CalculateByType(int CalculateType, char *request, char **response); // options: JSON_REQ_JSON_RESP
	int Discovery(char **response);

	int set_sort(int s) { sort_needed = s; return 0; } // use only for debug modes.
	void set_am_matrix(string s) { am_matrix = s; }
	void get_am_rep_signals(unordered_set<string> &am_sigs) { ma.get_rep_signals(am_sigs); } // returns the available 
	void get_sig_structure(string &sig, int &n_time_channels, int &n_val_channels, int* &is_categ) { ma.get_signal_structure(sig, n_time_channels, n_val_channels, is_categ); }

	string get_sig_unit(const string &sig, int val_channel) { return ma.get_rep().sigs.unit_of_measurement(sig, val_channel); }

	string get_lib_code_version();

	void show_rep_data(char **response);
};

//===============================================================================
// SimpleExampleEGFRAlgoMarker - the simplest example for a different AM
//===============================================================================
extern "C" class DLL_WORK_MODE SimpleExampleEGFRAlgoMarker : public AlgoMarker {

private:

	// inputs for egfr
	float age = -1, gender = -1, creatinine = -1;
	int ethnicity = 0;

	int pid = -1; // this example AM supports only a single pid at a time and no batches

public:
	SimpleExampleEGFRAlgoMarker() { set_type((int)AM_TYPE_MEDIAL_INFRA); add_supported_stype("Raw"); this->set_name("SimpleEGFR"); }

	int Load(const char *config_f) { ClearData(); return AM_OK_RC; }
	int Unload() { return AM_OK_RC; }
	int ClearData() { age = -1; gender = -1; creatinine = -1; return AM_OK_RC; }
	int AddData(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
	int Calculate(AMRequest *request, AMResponses *responses);

};

#endif //ALGOMARKER_FLAT_API


//========================================================================================
// Actual AlgoMarker API to be used via C# :: External users should work ONLY with these !!
// Any change to the below functions must rely only on exposed API of the above classes.
//========================================================================================

// all return codes are defined in AlgoMarkerErr.h

// create a new AlgoMarker of type am_type and init its name
extern "C" DLL_WORK_MODE int AM_API_Create(int am_type, AlgoMarker **new_am);

// loading AlgoMarker and making it ready to get Requests
extern "C" DLL_WORK_MODE int AM_API_Load(AlgoMarker* pAlgoMarker, const char *config_fname);

// Additional load options for AlgoMarker
extern "C" DLL_WORK_MODE int AM_API_AdditionalLoad(AlgoMarker* pAlgoMarker, const int load_type, const char *load);

// clearing data from AlgoMarker (recommended at the start and/or end of each query session
extern "C" DLL_WORK_MODE int AM_API_ClearData(AlgoMarker* pAlgoMarker);

// adding data to an AlgoMarker
// this API allows adding a specific signal, with matching arrays of times and values
extern "C" DLL_WORK_MODE int AM_API_AddData(AlgoMarker* pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
extern "C" DLL_WORK_MODE int AM_API_AddDataStr(AlgoMarker* pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values);

// adding data in a new DataType
extern "C" DLL_WORK_MODE int AM_API_AddDataByType(AlgoMarker* pAlgoMarker, const char *data, char **messages);

// Prepare a Request
// Null RC means failure
extern "C" DLL_WORK_MODE int AM_API_CreateRequest(char *requestId, char **score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req);

// Create a new empty responses object
extern "C" DLL_WORK_MODE int AM_API_CreateResponses(AMResponses **new_responses);

// Get scores for a ready request
extern "C" DLL_WORK_MODE int AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses);

// Get scores in general types
extern "C" DLL_WORK_MODE int AM_API_CalculateByType(AlgoMarker *pAlgoMarker, int CalcType, char *request, char **responses);

// get Responses num
extern "C" DLL_WORK_MODE int AM_API_GetResponsesNum(AMResponses *responses);

// get messages shared for all responses
extern "C" DLL_WORK_MODE int AM_API_GetSharedMessages(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args);

// get a response at a specific index
extern "C" DLL_WORK_MODE int AM_API_GetResponseIndex(AMResponses *responses, int _pid, long long _timestamp);

// get the request id that was used to create a responses object
extern "C" DLL_WORK_MODE int AM_API_GetResponsesRequestId(AMResponses *responses, char **requestId);

// get a score from responses given the response index and the type of the score
extern "C" DLL_WORK_MODE int AM_API_GetResponseScoreByType(AMResponses *responses, int res_index, char *_score_type, float *out_score);

// get a response at a certain index in responses
extern "C" DLL_WORK_MODE int AM_API_GetResponseAtIndex(AMResponses *responses, int index, AMResponse **response);

// get number of different scores in a response
extern "C" DLL_WORK_MODE int AM_API_GetResponseScoresNum(AMResponse *response, int *n_scores);

// get response score at a given score_index, returns pid, ts, score, and score_type
extern "C" DLL_WORK_MODE int AM_API_GetResponseScoreByIndex(AMResponse *response, int score_index, float *score, char **_score_type);

// get response score at a given score_index, returns pid, ts, score, and score_type
extern "C" DLL_WORK_MODE int AM_API_GetResponseExtendedScoreByIndex(AMResponse *response, int score_index, char **ext_score, char **_score_type);

// get messages for a response : messages that are score independent (such as raw eligibility tests)
extern "C" DLL_WORK_MODE int AM_API_GetResponseMessages(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args);

// get score messages for a specific score index
extern "C" DLL_WORK_MODE int AM_API_GetScoreMessages(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args);

// get pid and timepoint of a response 
extern "C" DLL_WORK_MODE int AM_API_GetResponsePoint(AMResponse *response, int *pid, long long *timestamp);

// get the name of an algomarker
extern "C" DLL_WORK_MODE int AM_API_GetName(AlgoMarker *pAlgoMArker, char **name);

// Dispose of AlgoMarker - free all memory 
extern "C" DLL_WORK_MODE void AM_API_DisposeAlgoMarker(AlgoMarker *pAlgoMarker);

// Dispose of AMRequest - free all memory 
extern "C" DLL_WORK_MODE void AM_API_DisposeRequest(AMRequest *pRequest);

// Dispose of responses - free all memory
extern "C" DLL_WORK_MODE void AM_API_DisposeResponses(AMResponses *responses);

extern "C" DLL_WORK_MODE void AM_API_Discovery(AlgoMarker *pAlgoMarker, char **resp);

// Dispose of allocated memory
extern "C" DLL_WORK_MODE void AM_API_Dispose(char *data);

//Show memory:
extern "C" DLL_WORK_MODE int AM_API_DebugRepMemory(AlgoMarker *pAlgoMarker, char **resp);

//========================================================================================
// Follows is a simple API to allow access to data repositories via c#
//========================================================================================
#ifndef ALGOMARKER_FLAT_API

extern "C" class DLL_WORK_MODE RepositoryHandle {
public:
	string fname;
	vector<int> pids;
	vector<string> signals;
	MedRepository rep;
};

extern "C" class DLL_WORK_MODE SignalDataHandle {
public:
	UniversalSigVec usv;
};

// create a Repository Handle and read data into memory, given file name, pids, and signals return: 0: OK -1: failed
extern "C" DLL_WORK_MODE int DATA_API_RepositoryHandle_Create(RepositoryHandle **new_rep, char *fname, int *pids, int n_pids, char **sigs, int n_sigs);

// Create a SignalDataHandle : can be reused for each read later, create several if working in parallel
extern "C" DLL_WORK_MODE int DATA_API_SignalDataHandle_Create(SignalDataHandle **new_sdh);

// Read Data into a signal data handle from a repository, returns len=number of elements read
extern "C" DLL_WORK_MODE int DATA_API_ReadData(RepositoryHandle *rep_h, int pid, char *sig, SignalDataHandle *sdh, int *len);

// Get a time channel at an index from a loaded SignalDataHandle
extern "C" DLL_WORK_MODE int DATA_API_GetTime(SignalDataHandle *sdh, int idx, int time_channel, int *time);

// Get a value channel at an index from a loaded SignalDataHandle
extern "C" DLL_WORK_MODE int DATA_API_GetVal(SignalDataHandle *sdh, int idx, int val_channel, float *val);

// dispose of SignalDataHandle
extern "C" DLL_WORK_MODE void DATA_API_Dispose_SignalDataHandle(SignalDataHandle *sdh);

// dispose of RepositoryHandle
extern "C" DLL_WORK_MODE void DATA_API_Dispose_RepositoryHandle(RepositoryHandle *rep_h);

#endif //ALGOMARKER_FLAT_API
