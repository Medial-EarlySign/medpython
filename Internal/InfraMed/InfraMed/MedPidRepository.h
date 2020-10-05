//
// MedPidRepository :
// Extending the MedRepository data structure to support creating and using a repository
// organized as patient by patient rather than signal by signal.
// This allows for (Very) fast retrieval of ALL the signals for a certain pid
// Also - this allows for memory efficient use of a repository in cases where a sweep over all (or part of) patients is needed
// This helps when creating feature matrices for a patient.
//
// Tools for handling thread safe reading of data using a predefined amount of memory to use are also available here.
// The general idea is holding a pre allocated cyclic buffer for the data (enough such that there's always free space for the next ones needed by the next threads.)
//

#ifndef __MED_PID_REPOSITORY_
#define __MED_PID_REPOSITORY_

#include "InfraMed.h"
#include <InfraMed/InfraMed/MedSparseVec.h>

#define MAX_PID_DATA_SIZE	10000000
class MedPidRepository;

class PosLen {
  public:
	int pos;
	int len;
	int do_split; // if 1, means that when updating this pos,len we need to move version to a new allocated memory. if 0 , it may move due to other reasons.
	PosLen& operator =(const int a) { pos = (unsigned long long)a; len=a; return *this; }
	bool operator==(const PosLen a) { return (pos == a.pos && len == a.len);}
};

class PidIdxRec {
  public:
	unsigned short fnum;
	unsigned long long pos;
	unsigned int byte_len;
	int idx;

	PidIdxRec& operator =(const int a) { fnum=(unsigned short)a; pos = (unsigned long long)a; byte_len=(unsigned int)a; idx = a; return *this; }

	PidIdxRec() { idx = -1; }
};

class PidRec {
	public:
		int pid;							// pid num of this record
		unsigned char *data;				// pointer to actual record in memory
		unsigned int data_len;				// actual size used (always <= data_size)
		unsigned int data_size;				// max size available
		int is_allocated;					// was the space for data allocated (and hence we need to free it) or not (allocated some other place).
		MedSparseVec<PosLen> sv;			// from serial sid to a pair of pos, len
		MedPidRepository *my_rep;			// needed for the get() method in order to get access to dictionaries
		MedRepository *my_base_rep;			// needed for the get() method in order to get access to dictionaries
		int allow_realloc;					// allow reallocation of data in read if given not enough space

		PidRec() { pid = -1; data = NULL; data_len = 0; data_size = 0; is_allocated = 0; my_rep = NULL; sv.clear(); allow_realloc = 1; }

		// after reading the data to *data this operation is needed to build the sparse vec from (serial) sid to PosLen
		int init_sv();

		// get methods - no need for pid just by signal
		void *get(const string &sig_name, int &len);
		void *get(int sid, int &len);

		// universal API
		UniversalSigVec usv;	// we keep a usv inside, to allow saving of the init() time
		inline void *uget(int sid, UniversalSigVec &_usv) { _usv.init(my_base_rep->sigs.Sid2Info.at(sid)); return (_usv.data = get(sid, _usv.len)); }
		inline void *uget(int sid) { return uget(sid, usv); }
		inline void *uget(const string &sig_name, UniversalSigVec &_usv) { return uget(my_base_rep->sigs.sid(sig_name), _usv); }
		inline void *uget(const string &sig_name) { return uget(sig_name, usv); }


		int init_from_rep(MedRepository *rep, int pid, vector<int> &sids_to_use);

		// memory alloc & free
		void prealloc(unsigned int len);
		int realloc(unsigned int len);
		int resize_data(unsigned int len) { return realloc(len); }
		void set_data_to_buffer() { data = &data_buffer[0]; }
		void free();

	private:
		vector<unsigned char> data_buffer;	// the actual holder of the data when using prealloc/realloc/resize_data
											// this makes for much easier data garbage collection (no need for free)
};


class MedPidRepository : public MedRepository {

  public:

	MedSparseVec<PidIdxRec> pids_idx;		// sparse vec from pid to an index record on files
	vector<MedBufferedFile> in_files;		// keeping open input files


	int init(const string &conf_fname);		// when using MedPidRepository, init it with this API, then use the load() APIs in MedRepository to load full signals.

	// creating the "by pid" index and data files for a range of given pids with at most "jump" pids in each file
	int create(string &rep_fname, int from_pid, int to_pid, int jump);

	// get data size of a pid (0 => pid not in data)
	unsigned int get_data_size(int pid);

	// if data is NULL it will be allocated and data_size will be the allocated size
	// if data is not NULL data_size should contain the max size allowed on input, and on output contains the actual size used
	// general error : -1
	// error due to insufficient data_size (non NULL data) : -2
	int get_pid_rec(int pid, unsigned char *&data, unsigned int &data_size, PidRec &prec);

	// same but simpler API, will use the data and data_len inside the prec instead
	// it is recommended for use with pre allocation of enough space in prec.data when going to reuse the same prec for reads.
	int get_pid_rec(int pid, PidRec &prec);

	~MedPidRepository();

};


//
// PidDynamicRec is a PidRec with additional options to have "versions" for each signal
// These "versions" are initially all pointing to the original signal version.
// A new version can be created out of an existing one by changing it or loading it.
// This allows for example for complex scenarios in cleaning and feature generation in which 
// a version of the signal per time point is needed.
// A PidDynamicRec is a PidRec and hence can be read from a repository using the get_pid_rec() API.
// Another option to initialize a PidDynamicRec is to use init_from_rep , which will load it from repository data that's already in memory.
//
class PidDynamicRec : public PidRec {

public:

	int set_n_versions(int n_ver);				// a version is always a positive (>0) number. 0 is kept as the version number of the original data.
												// this method should always be called AFTER the original version had been read.
	int set_n_versions(vector<int> &time_points); // works with init_from_rep version that uses time_points , see below

	int get_n_versions() { return n_versions; }

	// calling get without version is defined in PidRec and will return the original version
	void *get(const string &sig_name, int version, int &len);
	void *get(int sid, int version, int &len);
	void *get(int sid, int &len) { return PidRec::get(sid, len); }

	// universal API
	inline void *uget(int sid, int version, UniversalSigVec &_usv) { _usv.init(my_base_rep->sigs.Sid2Info.at(sid)); return (_usv.data = get(sid, version, _usv.len)); }
	inline void *uget(int sid, int version) { return uget(sid, version, usv); }
	inline void *uget(const string &sig_name, int version, UniversalSigVec &_usv) { int sid = my_base_rep->sigs.sid(sig_name); if (sid < 0) HMTHROW_AND_ERR("Error: can't find signal %s in repository\n", sig_name.c_str()); return uget(sid, version, _usv); }
	inline void *uget(const string &sig_name, int version) { return uget(sig_name, version, usv); }

	// clearing
	void clear_vers(); // deletes all versions and remains just with the original one.

	// creating and changing versions
	int set_version_data(int sid, int version, void *datap, int len);
	int set_version_off_orig(int sid, int version); // if the version is still pointing the original area, we will make a copy of it outside
	int point_version_to(int sid, int v_src, int v_dst);	// will point version v_dst to the data of version v_src
	int remove(int sid, int version, int idx);			// removing element idx from version
	int remove(int sid, int v_in, int idx, int v_out);	// removing element idx from version v_in and putting it in v_out
	int change(int sid, int version, int idx, void *new_elem);	// changing element idx in version to hold *new_elem
	int change(int sid, int v_in, int idx, void *new_elem, int v_out);	// changing element idx in v_in to *new_elem, and putting it all in v_out
	int update(int sid, int v_in, vector<pair<int, void *>>& changes, vector<int>& removes); // Apply changes and removals
	int update(int sid, int v_in, int val_channel, vector<pair<int, float>>& changes, vector<int>& removes); // Apply val changes and removals, unified variant
	int update(int sid, int v_in, vector<pair<int, vector<float>>>& changes, vector<int>& removes); /// Apply val changes to all channels + removals

	// an API to push a data vector universal signal style , given all its value and time channels
	// when there are several channels (for time or values) , they should be placed one after the other for each member
	// so for example a *_times may have : t_0_ch0, t_0_ch1, t_1_ch0, t_1_ch1, ...
	int set_version_universal_data(int sid, int version, int *_times, float *_vals, int len);

	// test if two versions point to the same place in memory
	int versions_are_the_same(int sid, int v1, int v2) { return ((int)((*get_poslen(sid, v1)) == (*get_poslen(sid, v2)))); }
	int versions_are_the_same(set<int> sids, int v1, int v2) { 
		for (int sid : sids)
			if (!versions_are_the_same(sid, v1, v2))return(0);
		return(1);
	}
	// a few debug helpers
	int print_ver(int sid, int ver);
	int print_all_vers(int sid);
	int print_all(); // print all len > 0 signals in all vers
	int print_sigs(const vector<string> &sigs); // print some sigs we need (debugging tool)
	
	PidDynamicRec() { n_versions = 0; }

	// next are options to init a PidDynamicRec from data that already resides in some part of a MedRepository that is already in memory.

	// most general initialization - just stating how many versions
	int init_from_rep(MedRepository *rep, int pid, vector<int> &sids_to_use, int n_versions);

	// init with time points: given N time points (in the SAME time units of the signals), this will initialize N versions for each signal
	// BUT will also set the len of each version such that all the data points for this signal are <= the matching time point.
	// This is very useful, as now version(i) will be a virtual snapshot of how the data looks if looking only at times <= time_point(i).
	int init_from_rep(MedRepository *rep, int pid, vector<int> &sids_to_use, vector<int> &time_points);


	vector<UniversalSigVec> usvs; /// for efficient usage of a record we may want to keep several usvs initialized and save init time.

private:
	int n_versions;
	unsigned int curr_len;
	MedSparseVec<PosLen> sv_vers;
	PosLen *get_poslen(int sid, int version) { if (version >= n_versions) return NULL; return sv_vers.get((unsigned int)(my_base_rep->sigs.sid2serial[sid])*n_versions+version); }
	void set_poslen(int sid, int version, PosLen pl) { sv_vers[(unsigned int)my_base_rep->sigs.sid2serial[sid]*n_versions+version] = pl; }
};


//
// Dynamic Version iterator allows iteration over versions in two manners -
// 1. For version-dependent operations - forward iterate over all versions
// 2. For verion-independent operations - backward iterate only over versions that are already different
//
// Methods are -
// constructors
// init : Get first version to work on
// next : Get next version to work on
// done : have we analyzed all versions ?
//

class versionIterator {

protected:
	PidDynamicRec *my_rec;
	int iVersion;
	set<int> signalIds;

public:
	versionIterator(PidDynamicRec& _rec, int signalId) { my_rec = &_rec; signalIds = { signalId }; }
	versionIterator(PidDynamicRec& _rec, set<int>& _signalIds) { my_rec = &_rec; signalIds = _signalIds; }

	virtual int init() { return 0; };
	virtual int next() { return 0; };
	virtual bool done() { return true; };
};

class allVersionsIterator : public versionIterator {

public:
	allVersionsIterator(PidDynamicRec& _rec, int signalId) : versionIterator(_rec, signalId) {}
	allVersionsIterator(PidDynamicRec& _rec, set<int>& _signalIds) : versionIterator(_rec, _signalIds) {}

	int init() { iVersion = 0; return iVersion; }
	int next() { return ++iVersion; }
	bool done() { return iVersion >= my_rec->get_n_versions(); }
};

class differentVersionsIterator : public versionIterator {

	int jVersion;

public:
	differentVersionsIterator(PidDynamicRec& _rec, int signalId) : versionIterator(_rec, signalId) {}
	differentVersionsIterator(PidDynamicRec& _rec, set<int>& _signalIds) : versionIterator(_rec, _signalIds) {}

	int init();
	int next();
	bool done() { return iVersion < 0; }

	inline int block_first() { return jVersion+1; }
	inline int block_last() { return iVersion; }
};

#endif

