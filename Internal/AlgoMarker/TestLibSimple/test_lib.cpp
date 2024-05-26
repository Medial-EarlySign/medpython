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

	virtual int Discovery(char **response) { *response = NULL; return 0; }
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
	typedef int(*t_AM_API_CalculateByType)(AlgoMarker *pAlgoMarker, int CalcType, char *request, char **responses);
	typedef int(*t_AM_API_GetName)(AlgoMarker * pAlgoMArker, char **name);
	typedef void(*t_AM_API_DisposeAlgoMarker)(AlgoMarker*);
	typedef void(*t_AM_API_Dispose)(char *);
	typedef void(*t_AM_API_Discovery)(AlgoMarker *pAlgoMarker, char **resp);
	void *addr_AM_API_Create = nullptr;
	void *addr_AM_API_Load = nullptr;
	void *addr_AM_API_AdditionalLoad = nullptr;
	void *addr_AM_API_ClearData = nullptr;
	void *addr_AM_API_AddData = nullptr;
	void *addr_AM_API_AddDataStr = nullptr;
	void *addr_AM_API_AddDataByType = nullptr;
	void *addr_AM_API_CalculateByType = nullptr;
	void *addr_AM_API_GetName = nullptr;
	void *addr_AM_API_DisposeAlgoMarker = nullptr;
	void *addr_AM_API_Dispose = nullptr;
	void *addr_AM_API_Discovery = nullptr;
	// returns index in sos
	static int load(const char * am_fname);
	static DynAM* so;
	static std::vector<DynAM> sos;
	static void set_so_id(int id) { assert(id >= 0 && id < (int)sos.size()); so = &sos[id]; };

	static int AM_API_ClearData(AlgoMarker * pAlgoMarker);
	static void AM_API_DisposeAlgoMarker(AlgoMarker * pAlgoMarker);
	static void AM_API_Dispose(char *data);
	static int AM_API_GetName(AlgoMarker * pAlgoMArker, char **name);
	static int AM_API_Create(int am_type, AlgoMarker **new_am);
	static int AM_API_Load(AlgoMarker * pAlgoMarker, const char *config_fname);
	static int AM_API_AdditionalLoad(AlgoMarker * pAlgoMarker, const int load_type, const char *load);
	static int AM_API_AddData(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
	static int AM_API_AddDataStr(AlgoMarker * pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values);
	static int AM_API_AddDataByType(AlgoMarker * pAlgoMarker, const char *data, char **messages);
	static int AM_API_CalculateByType(AlgoMarker *pAlgoMarker, int CalcType, char *request, char **responses);
	static void AM_API_Discovery(AlgoMarker *pAlgoMarker, char **resp);

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
	if (DynAM::load(am_fname) < 0)
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
	so->addr_AM_API_CalculateByType = load_sym(lib_handle, "AM_API_CalculateByType");
	so->addr_AM_API_GetName = load_sym(lib_handle, "AM_API_GetName");
	so->addr_AM_API_DisposeAlgoMarker = load_sym(lib_handle, "AM_API_DisposeAlgoMarker");
	so->addr_AM_API_Dispose = load_sym(lib_handle, "AM_API_Dispose");
	so->addr_AM_API_Discovery = load_sym(lib_handle, "AM_API_Discovery", false);
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

void DynAM::AM_API_Dispose(char *data) {
	(*((DynAM::t_AM_API_Dispose)DynAM::so->addr_AM_API_Dispose))
		(data);
}

int DynAM::AM_API_GetName(AlgoMarker * pAlgoMArker, char **name) {
	return (*((DynAM::t_AM_API_GetName)DynAM::so->addr_AM_API_GetName))
		(pAlgoMArker, name);
}

void DynAM::AM_API_Discovery(AlgoMarker * pAlgoMArker, char **discovery) {
	if (DynAM::so->addr_AM_API_Load == NULL) {
		printf("AM_API_Discovery is NULL\n");
		return;
	}
	return (*((DynAM::t_AM_API_Discovery)DynAM::so->addr_AM_API_Discovery))
		(pAlgoMArker, discovery);
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


void init_and_load_data(const char *input_json_path, AlgoMarker *am, int &pid) {
	DynAM::AM_API_ClearData(am);

	string in_jsons;
	char *out_messages;
	if (read_file_into_string(input_json_path, in_jsons) < 0) {
		printf("Error on loading file %s\n", in_jsons.c_str());
		throw logic_error("Error");
	}
	int pos_pid = in_jsons.find_first_of("\"patient_id\"");
	if (pos_pid != string::npos) {
		string rest = in_jsons.substr(pos_pid + 13);
		int comma_pos = rest.find_first_of(',');
		if (comma_pos != string::npos) {
			rest = rest.substr(0, comma_pos);
			try {
				pid = stoi(rest);
			}
			catch (...) {
				printf("Failed in fetching patient id from: \"%s\"\n", rest.c_str());
			}
		}
	}
	printf("read %zu characters from input jsons file %s\n", in_jsons.length(), input_json_path);
	int load_status = DynAM::AM_API_AddDataByType(am, in_jsons.c_str(), &out_messages);
	if (out_messages != NULL) {
		string msgs = string(out_messages); //New line for each message:
		printf("AddDataByType has messages:\n");
		printf("%s\n", msgs.c_str());
	}
	DynAM::AM_API_Dispose(out_messages);
	printf("Added data from %s\n", input_json_path);
	if (load_status != AM_OK_RC)
		printf("Error code returned from calling AddDataByType: %d\n", load_status);
}


int main(int argc, char *argv[]) {

	if (argc <= 2) {
		printf("Please pass path to lib + amconfig + (optional data_json)\n");
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

	char *jdiscovery = NULL;
	DynAM::AM_API_Discovery(test_am, &jdiscovery);
	if (jdiscovery != NULL) {
		printf("Got discovery output:\n%s\n", jdiscovery);
		DynAM::AM_API_Dispose(jdiscovery);
	}

	if (argc > 3) {
		char *data_json_path = argv[3];
		init_and_load_data(data_json_path, test_am, pid_id);
	}

	string sjreq = "";
	sjreq = "{ \"type\" : \"request\", \"request_id\" : \"my test\", \"export\" : {\"prediction\" : \"pred_0\"}, \"requests\" : [{ \"patient_id\": \"" + to_string(pid_id) +
		"\", \"time\" : \"" + to_string(prediction_time) + "\" }] }";

	char *jresp = NULL;
	int calc_status = -1;
	if (argc > 3) {//3 and up
		char *jreq = (char *)(sjreq.c_str());
		calc_status = DynAM::AM_API_CalculateByType(test_am, JSON_REQ_JSON_RESP, jreq, &jresp);
	}
	if (jresp != NULL) {
		printf("Got respond with status %d\n%s\n", calc_status, jresp);
		DynAM::AM_API_Dispose(jresp);
	}
	else {
		if (argc > 3)
			printf("Got calcuate with status %d\n", calc_status);
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