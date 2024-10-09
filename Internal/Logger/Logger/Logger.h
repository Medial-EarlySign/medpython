/// @file
/// Logger.h - allowing logs with more control
///
/// next are codes for libraries sections using Logger\n
/// Application codes should use LOG_APP as the LOCAL_SECTION code\n
/// each c file using the logger should contain at top something like:\n
/// #include "Logger/Logger/Logger.h"\n
/// #define LOCAL_SECTION LOG_APP\n
/// #define LOCAL_LEVEL	LOG_DEF_LEVEL\n
/// extern MedLogger global_logger;\n
///\n
/// specific sections should use the section code instead of LOG_APP . \n
/// LOCAL_LEVEL can be changed to something else to allow more or less verbal logs.\n
///\n

#ifndef __LOGGER__H__
#define __LOGGER__H__

#include <stdio.h>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdarg.h>
#include <chrono>
#include <stdexcept>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/// user apps section
#define LOG_APP		 0
/// general default section
#define LOG_DEF		 1
// InfraMed sections
#define LOG_INFRA	 2
#define LOG_REP		 3
#define LOG_INDEX	 4
#define LOG_DICT	 5
#define LOG_SIG		 6
#define LOG_CONVERT	 7
// MedUtils sections
#define LOG_MED_UTILS	8
#define LOG_MEDMAT	 9
#define LOG_MEDIO	 10
#define LOG_DATA_STRUCTURES	11
// MedAlgo sections
#define LOG_MEDALGO	12
// MedFeat sections
#define LOG_MEDFEAT	13
// MedStat sections
#define LOG_MEDSTAT	14
// RepCleaner section
#define LOG_REPCLEANER 15
// FtrGenerator section
#define LOG_FTRGNRTR 16
// CrossValidator section
#define LOG_CV 17
// FeatCleaner section
#define LOG_FEATCLEANER 18
// ValueCleaner section
#define LOG_VALCLNR 19
// MedSamples section
#define MED_SAMPLES_CV 20
// FeatsSelector section
#define LOG_FEAT_SELECTOR 21
// SampleFilter section
#define LOG_SMPL_FILTER 22
// SerializableObject section
#define LOG_SRL 23
// MedModel section
#define LOG_MED_MODEL 24
// RepositoryType section
#define LOG_REPTYPE 25

#define MAX_LOG_SECTION	26


/// logs get printed if their given level (which can be different for different code sections)\n
/// is > log level in levels.\n
/// this means that min level never be printed, and max level is always printed\n

#define NO_LOG_LEVEL		0
#define MIN_LOG_LEVEL		0

#define DEBUG_LOG_LEVEL		3

#define LOG_DEF_LEVEL		5

#define MAX_LOG_LEVEL		10
#define VERBOSE_LOG_LEVEL	10


extern vector<string> log_section_to_name;
extern vector<string> log_level_to_name;

class MedLogger {
public:
	vector<int> levels;
	vector<vector<FILE *>> fds; //each section can write to multipal buffers
	vector<vector<string>> file_names; // required to avoid multiple closing of same file
	FILE *out_fd;
	string out_file_name;
	vector<string> format; ///< log format for each section - the message is %s, reset is variables with $: time,section,level

	MedLogger();
	~MedLogger();

	void init_all_levels(int level);
	int init_all_files(const string &fname);
	void init_level(int section, int level);
	void init_file(int section, FILE *of);
	int init_file(int section, const string &fname);
	int add_file(int section, const string &fname);

	void init_out(); //default output (stdout)
	void init_out(FILE *of);
	void init_out(const string &fname);

	int log(int section, int print_level, const char *fmt, ...);
	void out(char *fmt, ...);

	/// <summary>
	/// sets log format. it has the following variables: time,section,level
	/// to use variable use $ before the varaible name. example:
	/// init_format("$time $level $section %s")
	/// </summary>
	void init_format(int section, const string &new_format);

	/// <summary>
	/// sets log format for all sections - refer to init_format for more details
	/// </summary>
	void init_all_formats(int section, const string &new_format);
};

extern MedLogger global_logger;

/// LOG() - all print options : section and level
#define MEDLOG(Section,Level,fmt,...) global_logger.log(Section,Level,fmt, ##__VA_ARGS__)
/// MDBG() - use LOCAL_SECTION, Level is given
#define MDBG(Level,fmt,...)  global_logger.log(LOCAL_SECTION,Level,fmt, ##__VA_ARGS__)
/// MLOG() - use LOCAL_SECTION and LOCAL_LEVEL
#define MLOG(fmt,...) global_logger.log(LOCAL_SECTION, LOCAL_LEVEL, fmt, ##__VA_ARGS__)
/// MLOG_V() - use LOCAL_SECTION and VERBOSE_LOG_LEVEL
#define MLOG_V(fmt,...) global_logger.log(LOCAL_SECTION, VERBOSE_LOG_LEVEL, fmt, ##__VA_ARGS__)
/// MLOG_D() - use LOCAL_SECTION and DEBUG_LOG_LEVEL
#define MLOG_D(fmt,...) global_logger.log(LOCAL_SECTION, DEBUG_LOG_LEVEL, fmt, ##__VA_ARGS__)
/// MERR() - use LOCAL_SECTION , always print
#define MERR(fmt,...) global_logger.log(LOCAL_SECTION, MAX_LOG_LEVEL, fmt, ##__VA_ARGS__)
/// MWARN - use LOCAL_SECTION and one less than MAX level (used for MERR), so that we can easily skip them
#define MWARN(fmt,...) global_logger.log(LOCAL_SECTION, MAX_LOG_LEVEL-1, fmt, ##__VA_ARGS__)

#define MOUT(fmt,...) global_logger.out(fmt, ##__VA_ARGS__)

#define MTHROW_AND_ERR_STR(fmt,...) {char buff[300];snprintf(buff, sizeof(buff), fmt, ##__VA_ARGS__);global_logger.log(LOCAL_SECTION, MAX_LOG_LEVEL, buff); throw runtime_error(string(buff));}

#define MTHROW_AND_ERR(fmt,...) {char buff[300];snprintf(buff, sizeof(buff), fmt, ##__VA_ARGS__);global_logger.log(LOCAL_SECTION,MAX_LOG_LEVEL,"RunTime ERROR: ");global_logger.log(LOCAL_SECTION, MAX_LOG_LEVEL, buff); throw std::runtime_error("Error");}


// next works and compiles also in H files
#define HMTHROW_AND_ERR(fmt,...) {char buff[300];snprintf(buff, sizeof(buff), fmt, ##__VA_ARGS__);global_logger.log(0, MAX_LOG_LEVEL, buff); throw runtime_error(string(buff));}


/**
* MedTimer - a very simple class to allow very easy time measures
*/
class MedTimer {
public:
	string name;
	chrono::high_resolution_clock::time_point t[2];
	unsigned long long diff;

	MedTimer(const string &tname) { name = tname; }
	MedTimer() { name = string(""); }

	void start() { t[0] = chrono::high_resolution_clock::now(); }
	void take_curr_time() { t[1] = chrono::high_resolution_clock::now(); diff = (unsigned long long)(chrono::duration_cast<chrono::microseconds>(t[1] - t[0]).count()); }
	unsigned long long get_clock_micro() {
		auto t_now = chrono::high_resolution_clock::now();
		auto micro = chrono::duration_cast<chrono::microseconds>(t_now.time_since_epoch());
		return (unsigned long long)(micro.count());
	}

	double diff_microsec() { return (double)diff; }
	double diff_milisec() { return (double)diff / 1000.0; }
	double diff_sec() { return (double)diff / 1000000.0; }

};

/**
* class to print progress of long process - multithreaded or not by time interval
*/
class MedProgress {
public:
	chrono::high_resolution_clock::time_point tm_start;
	chrono::high_resolution_clock::time_point tm_prog;
	int progress;
	int max_loop;
	int print_interval;
	string print_title;
	int max_threads;
	bool alway_print_total;

	int print_section;
	int print_level;

	/// <summary>
	/// @param title title to print in each time
	/// @param mprocess_cnt number of actions to preform in total
	/// @param print_inter number of seconds between each print
	/// </summary>
	MedProgress(const string &title, int mprocess_cnt, int print_inter = 30, int max_th = 50);

	/// update when completed process 1 item
	void update();

	/// update action when skiping action - update job counter
	void skip_update();

	/// update action when without checking if needs to print, faster
	void update_no_check();
};

void get_current_time(string &time_str);

#endif

