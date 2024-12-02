//
// InfraMed.h - classes to use data repositories and their indexes
//

#ifndef __INFRAMED__H__
#define __INFRAMED__H__
//#define _USE_AS_DLL
#ifdef _USE_AS_DLL 
#ifdef __INFRAMED_DLL
#define DLLEXTERN __declspec(dllexport)
#else
#define DLLEXTERN __declspec(dllimport)
#endif
#else
#define DLLEXTERN   
#endif

#include <Logger/Logger/Logger.h>
#include "MedSparseVec.h"
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <map>
#include <vector>
#include <string>
#include "MedDictionary.h"
#include "MedSignals.h"
#include <fstream>
#include <thread>
#include <mutex>

using namespace std;

#define MED_MAGIC_NUM	0x0123456789abcdef
#define REPOSITORY_FULL_FORMAT 0x0
#define REPOSITORY_STRIPPED_FORMAT 0x1

#define MAX_PID_NUMBER	20000000

#define MAX_SID_NUMBER	100000
#define MAX_SIGNALS		1000

// next is needed in some places, leaving last 16 places for marking some needed states
#define MAX_ID_NUM ((unsigned int)0xfffffff0)

// general agreed upon usefull defines
#define GENDER_BOTH		0
#define GENDER_MALE		1
#define GENDER_FEMALE	2
#define MIN_DATE		19000101
#define MAX_DATE		21000101
#define INFRAMED_MIN_AGE			0
#define INFRAMED_MAX_AGE		150

class MedRepository;
class PidRec;
class InMemRepData;

class MedBufferedFile {
public:
	string name = "";
	int buf_size = 0;
	int buf_len = 0;
	unsigned char *buf = NULL;
	unsigned long long buf_start_in_file = 0;

	static const int default_buf_size = 128 * 1024;

	int open(const string &fname);
	int open(const string &fname, const int bsize);
	void close();
	unsigned long long read(unsigned char *outb, unsigned long long pos, unsigned long long len);

	ifstream *inf = NULL;


};

class IndexElem {
public:
	short file_num;
	unsigned long long pos_in_file;
	void *data;
	int len; // len is always in bytes (to allow for size/alloc calculations)

	IndexElem() { file_num = -1; pos_in_file = 0; data = NULL; len = 0; }
};

class CompactIndexElem {
public:	// 12 bytes per entry instead of 22
	unsigned int pos_add;
	unsigned int data_ptr_add;
	unsigned int len; // could be short however class is always being padded to 4x so using a full int

	CompactIndexElem() { pos_add = 0; data_ptr_add = 0; len = 0; }
};

//extern mutex *index_table_locks; //[MAX_SID_NUMBER];

// next is a fairly memory/speed efficient index to hold positions of a continous signal in a file/memory.
class IndexTable {
public:
	// thread safing
	//mutex m_lock;


	// using MedSparseVec - but also saving mem by getting len with the diff to the next entry
	MedSparseVec<unsigned int> sv;
	int last_len; // length of last element inserted to sv
	int sid;	  // sid for this table

	unsigned long long base;	// base position in a file or in memory
	unsigned int factor;		// size of a single record for the sid


	unsigned char *work_area;	// ptr to data in memory (NULL if there's no data loaded)
	int work_area_allocated;	// signing if we allocated the work_area for data (and need to delete[] it) or not.
	unsigned long long w_size;
	int is_loaded;
	int full_load;				// signing wheather we loaded a full signal or a partial one
	unsigned long long tot_size;
	double tot_size_gb;
	int is_locked;				// will not be freed if locked !!!

	unsigned int acc;			// accumulator for easier insertions

	void init() {
		sid = 0; sv.set_def(0); acc = 0; w_size = 0; tot_size = 0; tot_size_gb = 0; is_loaded = 0; full_load = 0;
		is_locked = 0; work_area_allocated = 0; work_area = NULL;
	}
	IndexTable() { init(); }

	int insert(unsigned int pid, int len) { int rc = insert(pid, acc, len); if (rc >= 0) acc += (unsigned int)len; return rc; }
	int insert(unsigned int pid, unsigned int delta, int len); // { last_len = len; return sv.insert(pid, delta);  }
	int get(const unsigned int pid, unsigned long long &pos, int &len);
	int get_len(const unsigned int pid);
	unsigned long long get_data_size(); // returns the size in bytes required for the actual data accompanying this index table

	void clear(); // also deallocated data if needed

	// locking prevents clearing or reading if already loaded
	void lock(); //{ /*lock_guard<mutex> guard(m_lock);*/ is_locked = 1; }
	void unlock(); // { /*lock_guard<mutex> guard(m_lock);*/ is_locked = 0; }

	// serializations
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// i/o to/from bin file
	int write_to_file(string &fname);
	int read_from_file(string &fname);

	// i/o to read index and data as well
	int read_index_and_data(string &idx_fname, string &data_fname);
	int read_index_and_data(string &idx_fname, string &data_fname, const vector<int> &pids_to_include);

};

class MedIndex {
public:
	int rep_mode;

	vector<string> ifnames;

	void clear() { ifnames.clear(); sid2idx.clear(); pid2idx.clear(); idx_i.clear(); idx.clear(); pids.clear(); signals.clear(); n_pids = 0; n_signals = 0; sid_in.clear(); }

	int read_index(vector<string> &fnames);
	int read_sub_index(vector<string> &fnames, const vector<int> &pids_to_include, const vector<int> &signals_to_include); // empty vector - means all are chosen

	inline void *get_ind_elem(int pid, int sid, unsigned int &len);
	inline IndexElem *get_ind_elem2(int pid, int sid);

	//int write_index(string &fname);

	unsigned long long get_index_max_data_size();
	void set_mem_ptrs_off();

	// next reads all the data parts that match the index in memory
	int read_all_data(unsigned char *&work_area, unsigned long long &wlen, vector<string> &data_fnames);
	//int MedIndex::read_all_data(unsigned char *&work_area, unsigned long long &wlen, vector<string> &data_fnames, vector<int> pids_to_take, vector<int> sids_to_take);
	int read_full_data(unsigned char *&work_area, unsigned long long &wlen, vector<string> &data_fnames);

	// next reads index and data for a single signal using an IndexTable index (mode 3 and up).
	int read_index_table_and_data(int sid, string &idx_fname, string &data_fname, const vector<int> &pids_to_include,
		unsigned char *w_area, unsigned long long &data_size);
	int read_index_table_and_data(int sid, string &idx_fname, string &data_fname, const vector<int> &pids_to_include);
	int update_pids();

	// new modes (2 and up) related
	//int get_idx_file_mode(const string &fname);


	// lists of all pids and signals in the index
	vector<int> pids;
	vector<int> signals;
	map<int, int> sid_in;

	// mode 3 index - (much less variables !!)
	vector<IndexTable> index_table;


	int contains_pid(int pid) { return (pid_idx[pid] >= 0); }

	MedRepository *my_rep;
	int min_pid_num;
	int max_pid_num;

	MedIndex() { min_pid_num = -1; max_pid_num = -1; rep_mode = 0; }

private:
	int mode;

	int get_mode(const string &fname);
	int read_index_mode0(const string &fname, const vector<int> &pids_to_include, const vector<int> &signals_to_include);
	int read_index_mode0_direct(const string &fname, const vector<int> &pids_to_include, const vector<int> &signals_to_include);
	int read_index_mode0_new_direct(const string &fname, int f_factor, const vector<int> &pids_to_include, const vector<int> &signals_to_include);
	int prep_idx_i();
	int prep_idx_i_direct();

	// mode 0 data elements
	int n_pids;
	int n_signals;
	map<int, int> sid2idx;
	map<int, int> pid2idx;
	vector<vector<unsigned int>> idx_i; // first dimension: signals, second: pids
	vector<IndexElem> idx;

	vector<vector<unsigned int>> idx_i_base;
	vector<vector<unsigned char>> idx_i_add;

	vector<vector<CompactIndexElem>> idx_recs;
	vector<IndexElem> idx_recs_base;
	vector<int> i_sid_type_byte_len;
	vector<int> i_sid_factor;

	vector<vector<unsigned int>> last_pid;
	//vector<int> idx_factor;
	//vector<int> idx_fno;


	// direct mapping
	vector<unsigned int> sid_idx;
	vector<unsigned int> pid_idx;
	int min_pid, max_pid;
	unsigned int n_pids_in_index;
	unsigned int n_signals_in_index;
	vector<int> pid_seen;

	vector<unsigned int> idx_pid;
	vector<unsigned int> idx_sid;



};

// InMemRepData is a class holding RAW data to be loaded into a repository 
// it works under the assumption of a small number of pids (1-1000 more or less).
// Large datasets would better be from an already made repository
// The advantages are that no loading from a file is needed and that there are API's to
// load data in and to connect this to a repository.
// API's are more c-style to enable easier work with C# API's
class InMemRepData : public SerializableObject {
public:
	MedRepository * my_rep;
	map<pair<int, int>, pair<int, vector<char>>> data;  // map from a pair of pid,sid to a pair of nvals, data vector

														// init_rep must be called first , as we must know the sigs names/types/etc...
	void init_rep(MedRepository &rep) { my_rep = &rep; }
	//int init_sigs(string &sigs_fname) { return my_rep->sigs.read({ sigs_fname }); }

	// insert data basic API's . These are for a single pid/sid but a vector of elements can be loaded
	// data can be insertred in any order and several times per pid/sid , different pids are supported too.
	int insertData(int pid, int sid, int *time_data, float *val_data, int n_time, int n_val);
	int insertData(int pid, const char *sig, int *time_data, float *val_data, int n_time, int n_val);
	static int insertData_to_buffer(int pid, int sid, int *time_data, float *val_data, int n_time, int n_val, 
	const MedSignals &sigs, map<pair<int, int>, pair<int, vector<char>>> &data);

	/// Erase pid data
	void erase_pid_data(int pid);


	// This sort action MUST be called after inserting all data, otherwise the order of the elements in each pid-sid vector will be the inserting order
	int sortData();
	int sort_pid_sid(int pid, int sid); // helper func - sort single vector

										// clearing 
	void clear() { data.clear(); }

	// a repository get function to use with this type of data
	static void *get_from_buffer(int pid, int sid, int &len,const map<pair<int, int>, pair<int, vector<char>>> &data);
	void *get(int pid, int sid, int &len);

	// debug and prints
	int print_all();
	int print(int pid);
	int print(int pid, int sid);

	// serializer for data
	ADD_SERIALIZATION_FUNCS(data);

};


class DLLEXTERN  MedRepository {
public:
	int rep_mode;
	string rep_files_prefix;

	string desc;
	string config_fname;
	string path;
	string metadata_path;
	vector<string> dictionary_fnames;
	vector<string> signal_fnames;
	vector<string> data_fnames;
	vector<string> index_fnames;
	string	fsignals_to_files;
	int min_pid_num;
	int max_pid_num;
	int time_unit;
	int format;

	MedIndex index;
	MedDictionarySections dict;
	MedSignals sigs;
	vector<int> pids;
	vector<int> all_pids_list;

	InMemRepData in_mem_rep; // for mode of running in in_mem

	//-------------------------------------------------------------------
	// most useful APIs - more overloads below
	//-------------------------------------------------------------------

	// reading a repository for a group of pids and signals. Empty group means all of it.
	int read_all(const string &conf_fname);
	int read_all(const string &conf_fname, const vector<int> &pids_to_take, const vector<string> &signals_to_take);
	int read_all(const string &conf_fname, const vector<int> &pids_to_take, const vector<int> &signals_to_take);

	// reading without any signal (will later be loaded with load()
	int init(const string &conf_fname);

	// getting the data for a pid,signal . Pointer to start returned, len elements inside. If not found NULL and 0 returned.
	inline void *get(int pid, const string &sig_name, int &len);
	inline void *uget(int pid, const string &sig_name, UniversalSigVec &usv);
	//		inline void *get(int pid, int sid, int &len);					// use this variant inside big loops to avoid map from string to int. // default variant
	inline void *get_all_modes(int pid, int sid, int &len);
	inline void *uget(int pid, int sid, UniversalSigVec_legacy &usv);		// Universal vec API, use this inside loops to avoid string map
	inline void *uget(int pid, int sid, GenericSigVec &gsv);
	void * (MedRepository::*get_ptr)(int, int, int&) = &MedRepository::get3;
	inline void *get(int pid, int sid, int &len) { return (this->*get_ptr)(pid, sid, len); }

	//-------------------------------------------------------------------

	inline void *get_in_mem(int pid, int sid, int &len) { return in_mem_rep.get(pid, sid, len); }
	void switch_to_in_mem_mode() { in_mem_rep.clear(); in_mem_rep.my_rep = (MedRepository *)this; get_ptr = &MedRepository::get_in_mem; }
	bool in_mem_mode_active() { return (get_ptr == &MedRepository::get_in_mem); }

	//-------------------------------------------------------------------
	int   read_config(const string &fname);

	int   read_dictionary();
	int   read_dictionary(const string &dname);

	int   read_signals();
	int   read_signals(const string &fname);

	int   read_index();
	int   read_index(const vector<int> &pids_to_take, const vector<int> &signals_to_take);
	int   read_index(vector<string> &fnames);
	int	  read_index(vector<string> &fnames, const vector<int> &pids_to_take, const vector<int> &signals_to_take);
	int	  read_index_tables(const vector<int> &pids_to_take, const vector<int> &signals_to_take);

	// reading the data matching the current index (always read an index BEFORE reading data)
	int read_data();
	int read_data(const string &fname);
	int read_data(vector<string> &fnames);
	//int read_data(string &fname, vector<int> &pids_to_take, vector<int> &signals_to_take);
	void free_data();

	// main initializing routines of repository from disk
	int read_all(const string &conf_fname, const vector<int> &pids_to_take, const vector<int> &signals_to_take, int read_data_flag);

	int read_all(const string &conf_fname, const vector<int> &pids_to_take, const vector<string> &signals_to_take, int read_data_flag);

	// mode 3 and up allows load and free of signals - overloads below
	int load(const vector<int> &sids, vector<int> &pids_to_take);

	// locking a sig in memory: will not be freed with free() functions - overloads below
	int lock_all_sigs();
	int lock(const vector<int> &sids);

	int unlock_all_sigs();
	int unlock(const vector<int> &sids);

	// when freeing: locked signals will NOT be freed, in order to free all unlock all first - overloads below
	int free_all_sigs();
	int free(const vector<int> &sids);

	double bound_gb;
	int free_to_bound() { return free_to_bound(bound_gb); }				// calls free_to_bound with default set memory bound
	int free_to_bound(double _bound_gb); // frees unlocked signals until getting below the bound

	// default is inifinity (1000gb), however - if set, the repository will try to free signals in order to get below this bound
	// it can not free locked in signals, so in some cases the overall memory usage could be higher
	int set_max_mem(double gb_mem) { bound_gb = gb_mem; return 0; }

	MedRepository() { path = ""; metadata_path = ""; work_area = NULL; work_size = 0; fsignals_to_files = ""; index.my_rep = this; min_pid_num = -1; max_pid_num = -1; rep_mode = 0; rep_files_prefix = "rep"; bound_gb = 100.0; sigs.my_repo = this; }
	~MedRepository() {
		//fprintf(stderr, "rep free\n"); fflush(stderr);
		if (work_area) {
			//fprintf(stderr, "~MedRepository before delete[] work_area %d", work_area);
			delete[] work_area;
			//fprintf(stderr, "~MedRepository after delete[] work_area %d", work_area);
			work_area = NULL;
		}
		free_all_sigs();
		//fprintf(stderr, "rep free ended\n"); fflush(stderr);
	};
	//int build_full_format_index();

	void clear();

	// main access routine: gets a pointer to data given pid and sid and its len (in vaiable type units).
	// When there's no data returns NULL and len==0
	inline void *get3(int pid, int sid, int &len);

	// common request: get a list of all values BEFORE (<) a given date
	SDateVal *get_before_date(int pid, int sid, int date, int &len);
	SDateVal *get_before_date(int pid, const string &sig_name, int date, int &len);

	// api for SDateVal (most common), to get ptr to a valur relative to a specific date
	// mode options are "==" "<=" ">=" "<" ">"
	// these return the first one exactly at the same date, or the first one just before (<= or <) etc.
	SDateVal *get_date(int pid, int sid, int date, const string &mode);
	SDateVal *get_date(int pid, const string &sig_name, int date, const string &mode);


	// test if repository has a specific pid or a specific signal
	int contains_pid(int pid) { return index.contains_pid(pid); }
	int contains_sid(int sid);

	// printing routines (mainly for debugging)
	void print_vec_dict(void *data, int len, int pid, int sid);
	void long_print_vec_dict(void *data, int len, int pid, int sid);
	void long_print_vec_dict(void *data, int len, int pid, int sid, int from, int to);
	void long_print_vec_dict(void *data, int len, int pid, int sid, int index);
	void print_channel_helper(int sid, int channel, float val);
	string get_channel_info(int sid, int channel, float val);
	string convert_date(int d, int sid);
	void print_data_vec_dict(int pid, int sid);
	void long_print_data_vec_dict(int pid, int sid);
	void long_print_data_vec_dict(int pid, int sid, int from, int to);
	void print_csv_vec(void * data, int len, int pid, int sid, bool dict_val);
	void convert_pid_sigs(const UniversalSigVec &usv, vector<pair<vector<string>, vector<string>>> &pid_result, const string &sig_name, int sig_id, int limit_count);
	void print_pid_sig(int pid, const string &sig_name, const vector<pair<vector<string>, vector<string>>> &usv);

	// getting all the dates for a pid in which there was at least one of the given signals
	int get_dates_with_signal(int pid, vector<string> &sig_names, vector<int> &dates);

	int get_pids_with_sig(const string &sig_name, vector<int> &in_pids); // getting all pids in rep that have at least once the signal.

	//-----------------------------------------------------------
	// useful overloads
	//-----------------------------------------------------------
	int load(const string &sig_name);
	int load(const int sid);
	int load(const vector<string> &sig_names);
	int load(const vector<int> &sids);
	int load(const string &sig_name, vector<int> &pids_to_take);
	int load(const int sid, vector<int> &pids_to_take);
	int load_pids_sorted(const int sid, vector<int> &pids_to_take);
	int load(const vector<string> &sig_names, vector<int> &pids_to_take);

	int lock(const string &sig_name);
	int lock(const int sid);
	int lock(const vector<string> &sig_names);

	int unlock(const string &sig_name);
	int unlock(const int sid);
	int unlock(const vector<string> &sig_names);

	int free(const string &sig_name);
	int free(const int sid);
	int free(const vector<string> &sig_names);

	int read_pid_list();

	int generate_fnames_for_prefix();

	/// <summary>
	///Option to read additional dictionary - should be called before!! read_all or init. Can be called multiple time for multiple dicts
	/// </summary>
	void load_additional_dict(const string &dict_path);
	/// <summary>
	///Clear additional dict list
	/// </summary>
	void clear_additional_dict();
private:
	int get_data_mode(const string &fname);
	unsigned char *work_area;
	unsigned long long work_size;

	vector<string> _addtional_dict_path;
};


//===============================================================
// simplifying iterations over several signals
class UsvsIterator {
public:
	MedRepository * rep;
	vector<string> sig_names;
	vector<UniversalSigVec *> usvs;

	int init(MedRepository *_rep, const vector<string> &_sig_names, const vector<UniversalSigVec *> &_usvs);
	int read_pid(int pid);
	int read_pid(int pid, const vector<UniversalSigVec *> &_usvs);

private:
	vector<int> sids;
};

//=============================================================================================
// Inline functions
//=============================================================================================
//------------------------------------------------------------
//inline IndexElem *MedIndex::get_ind_elem(int pid, int sid)
//{
//	if (pid2idx.find(pid) == pid2idx.end() || sid2idx.find(sid) == sid2idx.end())
//		return NULL;
//	int i_pid = pid2idx[pid];
//	int i_sid = sid2idx[sid];
//	int ie_i = idx_i[i_pid*n_signals + i_sid];
//
//	if (ie_i < 0)
//		return NULL;
//
//	return &idx[ie_i];
//}
inline IndexElem *MedIndex::get_ind_elem2(int pid, int sid)
{
	unsigned int i_pid = pid_idx[pid];
	unsigned int i_sid = sid_idx[sid];

	//if (i_pid<0 || i_sid<0)
	if (i_sid >= MAX_ID_NUM || i_pid >= MAX_ID_NUM)
		return NULL;

	//	int ie_i = idx_i[i_pid*n_signals + i_sid];
	int ie_i = idx_i[i_sid][i_pid];

	if (ie_i < 0)
		return NULL;

	return &idx[ie_i];
}


inline void *MedIndex::get_ind_elem(int pid, int sid, unsigned int &len)
{
	unsigned int i_pid = pid_idx[pid];
	unsigned int i_sid = sid_idx[sid];

	if (i_sid >= MAX_ID_NUM || i_pid >= MAX_ID_NUM)
		return NULL;

	unsigned char ie_i_add = idx_i_add[i_sid][i_pid];

	if (ie_i_add == 0) {
		//fprintf(stderr, "pid %d sid %d i_pid %d i_sid %d ie_i_add is 0 !!!!\n", pid, sid, i_pid, i_sid);
		return NULL;
	}
	//fprintf(stderr, "pid %d sid %d i_pid %d i_sid %d ie_i_add is not 0 !!!!\n", pid, sid, i_pid, i_sid);
	int i_base = i_pid >> 7;

	int ie_i = idx_i_base[i_sid][i_base] + ie_i_add - 1;

	len = idx_recs[i_sid][ie_i].len;

	unsigned long long d_add = (unsigned long long)idx_recs[i_sid][ie_i].data_ptr_add*(unsigned long long)i_sid_factor[i_sid];

	return (void *)((unsigned long long)idx_recs_base[i_sid].data + d_add);
}

//-----------------------------------------------------------
// returns also len - in units of the relevant sid type (!) 
inline void *MedRepository::get(int pid, const string &sig_name, int &len)
{
	len = 0;
	int sid = sigs.sid(sig_name);
	if (sid < 0)
		return NULL;
	return(get(pid, sid, len));
}


//-----------------------------------------------------------
// returns also len - in units of the relevant sid type (!) 
inline void *MedRepository::get_all_modes(int pid, int sid, int &len)
{
	//return ((get_ptr)(pid, sid, len));

	if (rep_mode >= 3)
		return get3(pid, sid, len);

	void *ie;

	unsigned int blen;

	len = 0;
	ie = index.get_ind_elem(pid, sid, blen);
	if (ie == NULL)
		return NULL;

	len = blen;
	return (ie);
}

//-----------------------------------------------------------
// get for mode 3 using index table 
inline void *MedRepository::get3(int pid, int sid, int &len)
{
	unsigned long long pos;
	index.index_table[sid].get(pid, pos, len);

	if (len == 0)
		return NULL;

	return ((void *)&index.index_table[sid].work_area[pos]);
}


//======================================================================
// Universal APIs
//======================================================================
inline void *MedRepository::uget(int pid, const string &sig_name, UniversalSigVec &usv)
{
	usv.len = 0;
	int sid = sigs.sid(sig_name);
	if (sid < 0)
		return NULL;
	return(uget(pid, sid, usv));
}

inline void *MedRepository::uget(int pid, int sid, UniversalSigVec_legacy &usv)
{
	usv.init(sigs.Sid2Info[sid]);

	usv.data = get(pid, sid, usv.len);
	return usv.data;

}

inline void *MedRepository::uget(int pid, int sid, GenericSigVec &gsv)
{
	gsv.init(sigs.Sid2Info[sid]);

	gsv.data = (char*)get(pid, sid, gsv.len);
	return gsv.data;
}



/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	* \brief signal_hierarchy namespace
	*/
	namespace signal_hierarchy {
		/// \brief filtering hierarchy codes
		string filter_code_hierarchy(const vector<string> &vec, const string &signalHirerchyType);
		/// \brief getting parents in hierarchy codes
		vector<int> parents_code_hierarchy(MedDictionarySections &dict, const string &group,
			const string &signalHirerchyType, int depth = 1, int max_nodes = 10);
		/// \brief getting sons in hierarchy codes
		vector<int> sons_code_hierarchy(MedDictionarySections &dict, const string &group, const string &signalHirerchyType);
		/// \brief getting sons in hierarchy codes - recursively till leaves. max_depth ==0 means no limit
		void sons_code_hierarchy_recursive(MedDictionarySections &dict,
			const string &signalName, const string &code, vector<int> &flat_all_sons, int max_depth = 0);
		/// \brief gets codes
		string get_readcode_code(MedDictionarySections &dict, int id, const string &signalHirerchyType);
	}
	/*!
	* \brief repository namespace
	*/
	namespace repository {
		/// \brief Helper function to calc diff between dates in years
		float DateDiff(int refDate, int dateSample);
		/// \brief Helper function to add days to date
		int DateAdd(int refDate, int daysAdd);
		/// \brief fetching specifc signal code value
		int get_value(MedRepository &rep, int pid, int sigCode);
		/// \brief fetching specifc signal code value
		int get_value(PidRec &rep, int sigCode);

		enum fix_method {
			none = 0,
			drop = 1,
			take_first = 2,
			take_last = 3,
			take_mean = 4,
			take_max = 5,
			take_min = 6
		};
		/// \brief fix contradicting signal values in same time for same patient. return true if changed data
		bool fix_contradictions(UniversalSigVec &s, fix_method method, UniversalSigVec_mem &edited);

		/// \brief fetches the next date from all signals in patientFile by date order.
		/// the signalPointers is array of indexes of each signal. it also advances the right index
		/// returns the signal with the minimal date - "the next date"
		template<class T> int fetch_next_date(vector<T> &patientFile, vector<int> &signalPointers);
		/// \brief sets global time unit by repository config file - has TIME_UNIT field which differ in ICU.
		/// so using it.
		void set_global_time_unit(const string &repository_path);
	}
}

#endif
