#ifndef __MP_PidRepository_H
#define __MP_PidRepository_H

#include "MedPyCommon.h"
#include "MPDictionary.h"

//#define AM_API_FOR_CLIENT

class MPSigVectorAdaptor;
class MPSig;
class MedPidRepository;
class MedConvert;
//class UniversalSigVec;
class MPSigExporter;

class MPPidRepository {
private:
	void get_sig_structure(string &sig, int &n_time_channels, int &n_val_channels, int* &is_categ);
	void AddDataStr(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, char** Values);
	void AddData(int patient_id, const char *signalName, int TimeStamps_len, long long* TimeStamps, int Values_len, float* Values);
	int _load_single_json(void *_js);
	bool data_load_sorted = false;
public:
	MEDPY_IGNORE(MedPidRepository* o);
	MPDictionary dict;

	MPPidRepository();
	~MPPidRepository();

#ifndef AM_API_FOR_CLIENT
	MEDPY_DOC(read_all, "read_all(conf_file_fname, [pids_to_take_array], [list_str_signals_to_take]) -> int\n"
	"returns -1 if fails\n"
	"reading a repository for a group of pids and signals.Empty group means all of it.");
	int read_all(const std::string &conf_fname);
	int read_all(const std::string &conf_fname, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take), const std::vector<std::string> &signals_to_take);

	int read_all_i(const std::string &conf_fname, const std::vector<int> &pids_to_take, const std::vector<int> &signals_to_take);
#endif

	
	MEDPY_DOC(init, "init(conf_file_name) -> int\n"
	"returns -1 if fails");
	int init(const std::string &conf_fname);
	
#ifndef AM_API_FOR_CLIENT
	MEDPY_DOC(loadsig, "loadsig(str_signame) -> int\n" 
	"  load a signal");
	int loadsig(const std::string& signame);

	MEDPY_DOC(loadsig, "loadsig(str_signame) -> int\n"
		"  load a signal for pids");
	int loadsig_pids(const std::string& signame, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take));
#endif

	MEDPY_DOC(sig_id, "list_signals() -> vector<string>\n"
		"  returns all repository signal names");
	std::vector<std::string> list_signals();

	MEDPY_DOC(sig_id, "sig_id(str_signame) -> int\n"
	"  returns signal id number for a given signal name");
	int sig_id(const std::string& signame);

	MEDPY_DOC(sig_type, "sig_type(str_signame) -> int\n"
	"  returns signal type");
	int sig_type(const std::string& signame);

	MEDPY_DOC(sig_description, "sig_description(str_signame) -> string\n"
		"  returns signal string description");
	std::string sig_description(const std::string& signame);

	MEDPY_DOC(is_categorical, "is_categorical(str_signame, int_val_channel) -> bool\n"
		"  returns True is channel is categorical");
	bool is_categorical(const std::string& signame, int val_channel);
	
	MEDPY_DOC(pids, "pids ; property(read) -> list_Int\n"
	"  returns array of pids");
	const std::vector<int>& MEDPY_GET_pids();
	
	MEDPY_DOC(uget, "uget(int_pid, int_sid) -> SigVectorAdaptor\n"
	"  returns a vector of universal signals");
	MPSigVectorAdaptor uget(int pid, int sid);

	MEDPY_DOC(dict_section_id, "dict_section_id(str_secName) -> int\n"
	"  returns section id number for a given section name");
	int dict_section_id(const std::string &secName);

	MEDPY_DOC(dict_name, "dict_name(int_section_id, int_id) -> string\n"
	"  returns name of section + id");
	std::string dict_name(int section_id, int id);
	
	MEDPY_DOC(dict_prep_sets_lookup_table, "dict_prep_sets_lookup_table(int_section_id, list_String set_names) -> BoolVector\n"
	"  returns a look-up-table for given set names");
	std::vector<bool> dict_prep_sets_lookup_table(int section_id, const std::vector<std::string> &set_names);
	
	MEDPY_DOC(get_lut_from_regex, "get get_lut_from_regex to names -> BoolVector\n"
		"  returns a lut  - boolean vector");
	std::vector<bool> get_lut_from_regex(int section_id, const std::string & regex_s);

#ifndef AM_API_FOR_CLIENT
	MEDPY_DOC(export_to_numpy, "export_to_numpy(str_signame) -> SigExporter\n"
	"  Returns the signal data represented as a list of numpy arrays, one for each field");
	MPSigExporter export_to_numpy(string signame, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take), int use_all_pids, int translate_flag, int free_sig, string filter_regex_str);

	MEDPY_DOC(free, "free(signame) -> int\n"
		"  Free the signal data specified by signame");
	int free(string signame);

	MEDPY_DOC_Dyn("void get_sig(const char * sig_name_str, bool translate=true, std::vector<int> * pids=nullptr, bool float32to64=true, bool free_signal=true, const char * regex_str=nullptr, const char * regex_filter=nullptr)");

#endif

	MEDPY_DOC(switch_to_in_mem, "switch_to_in_mem()\n"
		"  Switch to in mem repository mode");
	void switch_to_in_mem();
	MEDPY_DOC(load_from_json, "load_from_json(json_file_path)\n"
		"  Loads patient data into in-mem repository. If patient exists, adds more data");
	void load_from_json(const std::string &json_file_path);
	MEDPY_DOC(load_from_json_str, "load_from_json_str(json_string)\n"
		"  Loads patient data into in-mem repository. If patient exists, adds more data");
	void load_from_json_str(const std::string &json_content);
	MEDPY_DOC(clear, "load_from_json()\n"
		"  Clear repository memory. erase in-memory patient data");
	void clear();

	MEDPY_DOC(finish_load_data, "finish_load_data()\n"
		"  prepare for model apply. not need to call it specifically.");
	void finish_load_data();
};


class MPSig {
	int idx;
	void* o;

public:
	MEDPY_IGNORE(MPSig(void* _o, int index));
	MPSig(const MPSig& other);

	int time(int chan = 0);
	float val(int chan = 0);
	int timeU(int to_time_unit);
	int date(int chan = 0);
	int years(int chan = 0);
	int months(int chan = 0);
	int days(int chan = 0);
	int hours(int chan = 0);
	int minutes(int chan = 0);
};

class MPSigVectorAdaptor {
public:
	MEDPY_IGNORE(void* o);
	MEDPY_IGNORE(MPSigVectorAdaptor());
	MPSigVectorAdaptor(const MPSigVectorAdaptor& other);
	~MPSigVectorAdaptor();
	int __len__();
	MPSig __getitem__(int i);

	int MEDPY_GET_type();
	int MEDPY_GET_n_time_channels();
	int MEDPY_GET_n_val_channels();
	int MEDPY_GET_time_unit();
	int MEDPY_GET_size();
};


class MPMedConvert {
	public:
	MEDPY_IGNORE(MedConvert* o);
	MPMedConvert();
	~MPMedConvert();

	MEDPY_DOC(init, "init_load_params(load_args)");
	void init_load_params(const std::string &load_args);

	MEDPY_DOC(init, "create_rep(conf_fname)");
	void create_rep(const std::string &conf_fname);

	MEDPY_DOC(init, "create_index(conf_fname) - creates reverse index by pid");
	int create_index(std::string &conf_fname);
};

#endif	// !__MP_PidRepository_H