//
// a few tools to help converting data files into a repository
//

#ifndef __MED_CONVERT_H__
#define __MED_CONVERT_H__

#include "InfraMed.h"
#include "MedDictionary.h"
#include "MedSignals.h"

#include <vector>
#include <string>
#include <fstream>

using namespace std;

#define MAX_COLLECTED_DATA_SIZE 56   // = previous sizeof(collected_data)

class collected_data {
public:
	char buf[MAX_COLLECTED_DATA_SIZE];
	/*
	int type;
	int date;
	long long time ;
	long long time2 ;
	float val;
	int date2;
	short val1;
	short val2;
	short val3;
	short val4;
	long long longVal ;
	float f_val2;
	*/
	// nextis helpers in speeding up parsing
	int serial = -1;
	//int idx = -1;
	void zero() {
		memset(buf, 0, sizeof(buf));
		//type = 0; date = 0; val = 0; date2 = 0; time = 0; time2 = 0; longVal = 0; val1 = 0; val2 = 0; val3 = 0; val4 = 0; f_val2 = 0; 
	}
};


class pid_data {
public:
	int pid;
	vector<vector<collected_data>> raw_data;
};

struct sig_info {
	int sid;
	int serial;
	int fno;
	int type;
};

struct file_stat {
	string fname;
	int id;
	int n_lines;
	int n_relevant_lines;
	int n_parsed_lines;
	int n_bad_format_lines;

	file_stat() { fname = ""; id = -1; n_lines = 0; n_relevant_lines = 0; n_parsed_lines = 0; n_bad_format_lines = 0; }
};

class  MedConvert {
public:
	int mode;						// 0/1 - original mode (currently default) 2 - new mode (data and index file for each signal)
	int safe_mode;					// 0/1 - in safe_mode==1 loading will exit in several inconsistencies
	string	rep_files_prefix;		// general prefix for all files created in mode 2

	string config_fname;

	string path;
	string out_path;
	string code_to_signal_fname;			// format: signal_code, signal_name
	vector<string> dict_fnames;
	vector<string> sig_fnames;
	string signal_to_files_fname;
	string prefixes_fname;
	int relative;				// if 1 - we will put "." in dir in .repository in case of OUTDIR
	int default_time_unit;		/// internal representation for all dates. 
								/// MedTime::Date == 1 means store as YYYYMMDD
								/// MedTime::Minutes == 6 means store as minutes since 01/01/1900

	// next files should be sorted by pid
	string registry_fname;		// format: pid , date, location(string), stage (number)	// tab delimited 
	//string demographic_fname;	// format: pid , birth year, M/F
	vector<string> in_data_fnames;	// format: pid , signal_code, date/time , value
	vector<string> in_strings_data_fnames; // format: pid , signale_code, date/time , 
	vector<string> forced; // signal names of signals that each id MUST have (like GENDER and/or BYEAR , etc), ids without those will not be loaded into the repository
	vector<string> load_only; // loading signals only from this list and leave others as is, usefull for fixes and updates
							  // note that for efficiency it is recommended to have in the data file only the files needed for forced and load_only
							  // the default is an empty load_only, which means load all possible signals in the sig gile.
	vector<int> sids_to_load;

	//running parameters for load:
	int check_for_error_pid_cnt = 10000; ///< after how many pids to check for error. If 0 only at the end
	int read_lines_buffer = 100000; ///< how much lines to read for file
	int test_run_max_pids = 0; ///< If bigger than 1 - will run in dry run till that number of pids
	int allowed_unknown_catgory_cnt = 50; ///< how many unknown categories are allowed
	int allowed_missing_pids_from_forced_cnt = 0; ///< how many pids are allowed to be missing in forced signals. 0 means no limit
	double allowed_missing_pids_from_forced_ratio = 0.05; ///< how many pids are allowed to be missing in forced signals. 0 means no limit
	double max_bad_line_ratio = 0.05; ///< maximal ratio for bad lines in file
	double min_parsed_line_ratio = 0.01; ///< minimal ratio for parsed lines in file
	bool verbose_open_files = false; ///< If true will print when openning files
	bool run_parallel = true; ///< If true will load in parallel
	bool run_parallel_files = false; ///< If true will read files in parallel
	string full_error_file; ///< provide full path to error file

	void init_load_params(const string &init_str);



	// outputs
	string repository_config_fname;
	vector<string> prefix_names;
	vector<string> index_fnames;
	vector<string> data_fnames;
	string description;

	// next are for debug and statistics
	vector<file_stat> fstats;
	map<string, int> missing_forced_signals;

	// internal variables
	unordered_map<string, string> codes2names;
	MedDictionarySections   dict;
	MedSignals sigs;
	map<int, int> sid2fno;
	unordered_map<int, int> sid2serial;
	vector<int> serial2sid;
	vector<sig_info> serial2siginfo;

	vector<int> pid_in_file;

	vector<IndexTable> indexes;

	void clear();

	MedConvert(const string &prefix) { rep_files_prefix = prefix; }
	MedConvert() { rep_files_prefix = "rep"; }

	// main entry points
	int create_rep(const string &config_fname, int _mode) { mode = _mode; return read_all(config_fname); }
	int create_rep(const string &config_fname) { return create_rep(config_fname, 1); }

	// configuration and preparations
	int read_config(const string &fname);
	int read_code_to_signal(const string &fname);
	int read_prefix_names(const string &fname);
	int read_signal_to_files(const string &fname);

	// mode 2 related
	int generate_prefix_names();

	// general prep function
	int read_all(const string &config_fname);

	int n_open_in_files;
	// actually reading data and creating index and data files
	void collect_lines(vector<string> &lines, vector<int> &f_i, int file_i, vector<string> &buffered_lines, int &buffer_pos, ifstream &inf, int file_type, pid_data &curr, int &fpid, file_stat& curr_fstat, map<pair<string, string>, int>&);
	void get_next_signal_all_lines(vector<string> &lines, vector<int> &f_i, pid_data &curr, vector<file_stat> &fstat, map<pair<string, string>, int>&);
	void parse_fields_into_gsv(string &curr_line, vector<string> &fields, int sid, GenericSigVec &cd_sv);
	int create_indexes();
	int create_repository_config();
	int create_signals_config();

	// legacy
	//void get_next_signal(vector<string> &buffered_lines, int &buffer_pos, ifstream &inf, int file_type, pid_data &curr, int &fpid, file_stat& curr_fstat, map<pair<string, string>, int>&);
	//int write_indexes(pid_data &curr);
	//void get_next_signal_new_modes(vector<string> &buffered_lines, int &buffer_pos, ifstream &inf, int file_type, pid_data &curr, int &fpid, file_stat& curr_fstat, map<pair<string, string>, int>&);


	// output files related
	ofstream signals_config_f;
	ofstream repository_config_f;
	vector<ofstream *> index_f;
	vector<ofstream *> data_f;
	vector<unsigned long long> data_f_pos;
	int open_indexes();
	int write_all_indexes(vector<int> &all_pids);
	int write_indexes_new_modes(pid_data &curr);
	int close_indexes();

	// loading subsets
	int prep_sids_to_load();
private:
	/// tests for load error during load. The input flag is indicator for test after finish load all
	void test_for_load_error(const map<pair<string, string>, int> &missing_dict_vals, int n_pids_extracted, bool final_test
		,int prev_total_missings, int &total_missing, const map<string, int> &prev_missing_forced_signals) const;

	ofstream err_log_file;
};

#endif