#define _CRT_SECURE_NO_WARNINGS

//
// Logger.cpp
//

#include <stdio.h>
#include <stdarg.h>
#include <cctype>
#include <thread>
#include <ctime>
#include <time.h>
#include "Logger.h"
#include <unordered_set>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN

#include <chrono>
#include <Windows.h>
#include <winsock2.h>

int gettimeofday(struct timeval* tp, struct timezone* tzp) {
	namespace sc = std::chrono;
	sc::system_clock::duration d = sc::system_clock::now().time_since_epoch();
	sc::seconds s = sc::duration_cast<sc::seconds>(d);
	tp->tv_sec = (long)s.count();
	tp->tv_usec = (long)sc::duration_cast<sc::microseconds>(d - s).count();

	return 0;
}
#else
#include <sys/time.h>
#endif // _WIN32

MedLogger global_logger;

vector<string> log_section_to_name = { "APP", "DEFAULT", "INFRA", "REP", "INDEX", "DICT", "SIG", "CONVERT", "MED_UTILS",
"MEDMAT", "MEDIO", "DATA_STRUCTURES", "MEDALGO", "MEDFEAT", "MEDSTAT", "REPCLEANER", "FTRGNRTR", "CV", "FEATCLEANER",
"VALCLNR", "MED_SAMPLES_CV", "FEAT_SELECTOR", "SMPL_FILTER", "SRL", "MED_MODEL", "REPTYPE" };
vector<string> log_level_to_name = { "min_level", "", "", "DEBUG", "", "INFO", "", "", "", "WARN", "ERROR" };

void get_current_time(string &time_str) {
	if (!time_str.empty())
		return;
	time_t theTime = time(NULL);
	struct tm *now;
	now = localtime(&theTime);
	char buffer[500];
	snprintf(buffer, sizeof(buffer), "[%d-%02d-%02d %02d:%02d:%02d]",
				now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
	time_str = string(buffer);
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_out()
{
	out_fd = stdout;
	out_file_name = "stdout";
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_out(FILE *of)
{
	out_fd = of;
	out_file_name = "__UNKNOWN__";
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_out(const string &fname)
{
	FILE *inf = fopen(fname.c_str(), "w");
	if (inf == NULL) {
		fprintf(stderr, "MedLogger: init_all_files: Can't open file %s - using default stdout instead\n", fname.c_str());
		out_fd = stdout;
	}
	out_fd = inf;
	out_file_name = fname;
}

//-----------------------------------------------------------------------------------------------
MedLogger::MedLogger()
{
	levels.resize(MAX_LOG_SECTION);
	fds.resize(MAX_LOG_SECTION);
	file_names.resize(MAX_LOG_SECTION);
	format.resize(MAX_LOG_SECTION, "%s");

	for (int i = 0; i < MAX_LOG_SECTION; i++) {
		levels[i] = LOG_DEF_LEVEL;
		fds[i].push_back(stderr);
		file_names[i].push_back("stderr");
	}
	init_out();
}

//-----------------------------------------------------------------------------------------------
MedLogger::~MedLogger()
{
	unordered_set<string> closed_files;

	for (int i = 0; i < MAX_LOG_SECTION; i++) {
		for (size_t j = 0; j < fds[i].size(); ++j)
		{
			if (fds[i][j] != NULL)
				fflush(fds[i][j]);

			if (fds[i][j] != NULL && fds[i][j] != stderr && fds[i][j] != stdout) {
				if (closed_files.find(file_names[i][j]) == closed_files.end()) {
					fclose(fds[i][j]);
					closed_files.insert(file_names[i][j]);
				}
				fds[i][j] = NULL;
			}
		}

	}
	if (out_fd != NULL && out_fd != stderr && out_fd != stdout && closed_files.find(out_file_name) == closed_files.end())
		fclose(out_fd);
}

//-----------------------------------------------------------------------------------------------
// Sets all logs to a given file
int MedLogger::init_all_files(const string &fname)
{
	FILE *inf = fopen(fname.c_str(), "w");
	if (inf == NULL) {
		fprintf(stderr, "MedLogger: init_all_files: Can't open file %s\n", fname.c_str());
		return -1;
	}

	for (int i = 0; i < MAX_LOG_SECTION; i++) {
		fds[i].assign(1, inf);
		file_names[i].assign(1, fname);
	}

	return 0;
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_all_levels(int level)
{
	for (int i = 0; i < MAX_LOG_SECTION; i++)
		levels[i] = level;
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_level(int section, int level)
{
	if (section < MAX_LOG_SECTION)
		levels[section] = level;
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_file(int section, FILE *of)
{
	if (section < MAX_LOG_SECTION) {
		if (fds[section].empty())
			fds[section].resize(1);
		fds[section].back() = of;
		file_names[section].back() = "__UNKNOWN__";
	}
}

//-----------------------------------------------------------------------------------------------
int MedLogger::init_file(int section, const string &fname)
{
	FILE *inf;

	if (section >= MAX_LOG_SECTION)
		return -2;

	inf = fopen(fname.c_str(), "w");
	if (inf == NULL) {
		fprintf(stderr, "MedLogger: init_file: Can't open file %s\n", fname.c_str());
		return -1;
	}

	if (fds[section].empty())
		fds[section].resize(1);
	fds[section].back() = inf;
	file_names[section].back() = fname;

	return 0;
}

int MedLogger::add_file(int section, const string &fname) {
	FILE *inf;

	if (section >= MAX_LOG_SECTION)
		return -2;

	inf = fopen(fname.c_str(), "w");
	if (inf == NULL) {
		fprintf(stderr, "MedLogger: init_file: Can't open file %s\n", fname.c_str());
		return -1;
	}
	fds[section].push_back(inf);
	file_names[section].push_back(fname);

	return 0;
}

void init_vars(string &str, int section, int print_level) {
	//search for all $ as sepcial char in format
	size_t idx = str.find("$");
	char buff[100];
	while (idx != string::npos) {
		//read till space:
		int var_idx = 0;
		while (var_idx < sizeof(buff) - 1 && !isspace(str[idx + 1 + var_idx])) {
			buff[var_idx] = tolower(str[idx + 1 + var_idx]);
			++var_idx;
		}
		buff[var_idx] = '\0';
		string var_name = string(buff);
		if (var_name == "timestamp") {
			struct tm *now;

			struct timeval now_ms;
			gettimeofday(&now_ms, NULL);
#if defined(__unix__) || defined(__APPLE__)
			now = localtime(&now_ms.tv_sec);
#else
			//allocate local mem
			struct tm now_m;
			now = &now_m;
			time_t timeT = now_ms.tv_sec;
			localtime_s(now, &timeT);
#endif

			snprintf(buff, sizeof(buff), "[%d-%02d-%02d %02d:%02d:%02d.%03ld]",
				now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec,
				now_ms.tv_usec);

		}
		else if (var_name == "time") {
			time_t theTime = time(NULL);
			struct tm *now;
#if defined(__unix__) || defined(__APPLE__)
			now = localtime(&theTime);
#else
			struct tm now_m;
			now = &now_m;
			localtime_s(now, &theTime);
#endif

			snprintf(buff, sizeof(buff), "[%d-%02d-%02d %02d:%02d:%02d]",
				now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);

		}
		else if (var_name == "level") {
			if (print_level < log_level_to_name.size())
				snprintf(buff, sizeof(buff), "[%s]", log_level_to_name[print_level].c_str());
			else
				snprintf(buff, sizeof(buff), "[Unknown_level]");
		}
		else if (var_name == "section") {
			if (section < log_section_to_name.size())
				snprintf(buff, sizeof(buff), "[%s]", log_section_to_name[section].c_str());
			else
				snprintf(buff, sizeof(buff), "[Unknown_section]");
		}
		else
			throw invalid_argument("Invalid log format: " + str);

		str.replace(idx, var_idx + 1, buff);

		idx = str.find("$");
	}
}

void MedLogger::init_format(int section, const string &new_format) {
	if (section < MAX_LOG_SECTION)
		format[section] = new_format;
}

void MedLogger::init_all_formats(int section, const string &new_format) {
	for (int i = 0; i < MAX_LOG_SECTION; i++)
		format[i] = new_format;
}

//-----------------------------------------------------------------------------------------------
int MedLogger::log(int section, int print_level, const char *fmt, ...)
{
	if (section >= MAX_LOG_SECTION)
		return -2;

	if (print_level < levels[section])
		return 1;

	bool wrote_log = false;
	char buff[5000];
	//prepare format:
	string format_e = format[section];
	init_vars(format_e, section, print_level);
	for (size_t i = 0; i < fds[section].size(); ++i)
	{
		if (fds[section][i] == NULL)
			continue;
		wrote_log = true;
		va_list args;
		va_start(args, fmt);
		//vfprintf(fds[section][i], fmt, args);
		vsnprintf(buff, sizeof(buff), fmt, args);
		va_end(args);
		fprintf(fds[section][i], format_e.c_str(), buff);
		fflush(fds[section][i]);
	}

	if (!wrote_log)
		return 1;
	return 0;
}

//-----------------------------------------------------------------------------------------------
void MedLogger::out(char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	vfprintf(out_fd, fmt, args);
	va_end(args);
	fflush(out_fd);
}

MedProgress::MedProgress(const string &title, int mprocess_cnt, int print_inter, int max_th) {
	progress = 0;
	max_loop = mprocess_cnt;
	print_interval = print_inter;
	print_title = title;
	max_threads = max_th;
	alway_print_total = false;

	print_section = LOG_APP;
	print_level = LOG_DEF_LEVEL;
	tm_prog = chrono::high_resolution_clock::now();
	tm_start = chrono::high_resolution_clock::now();
}

void MedProgress::update() {
#pragma omp atomic
	++progress;

	if (progress == max_loop) {
		double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
			- tm_start).count()) / 1000000.0;
		if (alway_print_total || time_elapsed > print_interval)
			global_logger.log(print_section, print_level, "#%s# :: Done processed all %d. Took %2.1f Seconds in total\n",
				print_title.c_str(), max_loop, time_elapsed);
		return;
	}

	double duration = (unsigned long long)(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
		- tm_prog).count()) / 1000000.0;
	if (duration > print_interval && progress % max_threads == 0) {
#pragma omp critical
		tm_prog = chrono::high_resolution_clock::now();
		double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
			- tm_start).count()) / 1000000.0;
		if (max_loop > 0) {
			double estimate_time = int(double(max_loop - progress) / double(progress) * double(time_elapsed));
			global_logger.log(print_section, print_level, "#%s# :: Processed %d out of %d(%2.2f%%) time elapsed: %2.1f Minutes, "
				"estimate time to finish %2.1f Minutes\n", print_title.c_str(),
				progress, max_loop, 100.0*(progress / float(max_loop)), time_elapsed / 60,
				estimate_time / 60.0);
		}
		else { //unknown job count: print job speed - how many jobs in minutes, and how much time a single job takes
			if (time_elapsed > 0) {
				double jobs_in_minute = double(progress) / (double(time_elapsed) / 60);
				double single_job_time_seconds = -1;
				if (jobs_in_minute > 1)
					single_job_time_seconds = 1 / (jobs_in_minute / 60);
				if (single_job_time_seconds < 0)
					global_logger.log(print_section, print_level, "#%s# :: Processed %d time elapsed: %2.1f Minutes, "
						"processing %2.1f jobs in a minute\n", print_title.c_str(),
						progress, time_elapsed / 60, jobs_in_minute);
				else
					global_logger.log(print_section, print_level, "#%s# :: Processed %d time elapsed: %2.1f Minutes, "
						"processing %2.1f jobs in a minute, single job take %2.4f seconds\n", print_title.c_str(),
						progress, time_elapsed / 60, jobs_in_minute, single_job_time_seconds);
			}
		}
	}
}

void MedProgress::skip_update() {
	if (max_loop > 0)
#pragma omp atomic
		--max_loop;
}

void MedProgress::update_no_check() {
#pragma omp atomic
	++progress;

	if (progress == max_loop) {
		double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
			- tm_start).count()) / 1000000.0;
		if (alway_print_total || time_elapsed > print_interval)
			global_logger.log(print_section, print_level, "#%s# :: Done processed all %d. Took %2.1f Seconds in total\n",
				print_title.c_str(), max_loop, time_elapsed);
	}
}