#ifndef __MP_PidRepository_H
#define __MP_PidRepository_H

#include "MedPyCommon.h"
#include "MPDictionary.h"

class MPSigVectorAdaptor;
class MPSig;
class MedPidRepository;
//class UniversalSigVec;
class MPSigExporter;


class MPPidRepository {
public:
	MEDPY_IGNORE(MedPidRepository* o);
	MPDictionary dict;

	MPPidRepository();
	~MPPidRepository();

	MEDPY_DOC(read_all, "read_all(conf_file_fname, [pids_to_take_array], [list_str_signals_to_take]) -> int\n"
	"returns -1 if fails\n"
	"reading a repository for a group of pids and signals.Empty group means all of it.");
	int read_all(const std::string &conf_fname);
	int read_all(const std::string &conf_fname, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take), const std::vector<std::string> &signals_to_take);

	int read_all_i(const std::string &conf_fname, const std::vector<int> &pids_to_take, const std::vector<int> &signals_to_take);


	
	MEDPY_DOC(init, "init(conf_file_name) -> int\n"
	"returns -1 if fails");
	int init(const std::string &conf_fname);
	
	MEDPY_DOC(loadsig, "loadsig(str_signame) -> int\n" 
	"  load a signal");
	int loadsig(const std::string& signame);
	
	MEDPY_DOC(sig_id, "sig_id(str_signame) -> int\n"
	"  returns signal id number for a given signal name");
	int sig_id(const std::string& signame);

	MEDPY_DOC(sig_type, "sig_type(str_signame) -> int\n"
	"  returns signal type id for a given signal name");
	int sig_type(const std::string& signame);
	
	MEDPY_DOC(sig_id, "pids ; property(read) -> list_Int\n"
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

	MEDPY_DOC(export_to_numpy, "export_to_numpy(str_signame) -> SigExporter\n"
	"  Returns the signal data represented as a list of numpy arrays, one for each field");
	MPSigExporter export_to_numpy(string signame, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take), int use_all_pids, int translate_flag, int free_sig);

	MEDPY_DOC(free, "free(signame) -> int\n"
		"  Free the signal data specified by signame");
	int free(string signame);

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



#endif	// !__MP_PidRepository_H