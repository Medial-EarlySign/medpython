#define _CRT_SECURE_NO_WARNINGS

//
// Logger.cpp
//

#include <stdio.h>
#include <stdarg.h>
#include "Logger.h"

MedLogger global_logger;

//-----------------------------------------------------------------------------------------------
void MedLogger::init_out()
{
	out_fd = stdout;
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_out(FILE *of)
{
	out_fd = of;
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_out(const string &fname)
{
	FILE *inf = fopen(fname.c_str(),"w");
	if (inf == NULL) {
		fprintf(stderr,"MedLogger: init_all_files: Can't open file %s - using default stdout instead\n",fname.c_str());
		out_fd = stdout;
	}
	out_fd = inf;
}

//-----------------------------------------------------------------------------------------------
MedLogger::MedLogger()
{
	levels.resize(MAX_LOG_SECTION);
	fds.resize(MAX_LOG_SECTION);

	for (int i=0; i<MAX_LOG_SECTION; i++) {
		levels[i] = LOG_DEF_LEVEL;
		fds[i] = stderr;
	}
	init_out();
}

//-----------------------------------------------------------------------------------------------
MedLogger::~MedLogger()
{
	for (int i=0; i<MAX_LOG_SECTION; i++) {
		if (fds[i] != NULL)
			fflush(fds[i]);
		if (fds[i] != NULL && fds[i] != stderr && fds[i] != stdout) {
			fclose(fds[i]);
		}
	}
	if (out_fd != NULL && out_fd != stderr && out_fd != stdout)
		fclose(out_fd);
}

//-----------------------------------------------------------------------------------------------
// Sets all logs to a given file
int MedLogger::init_all_files(const string &fname)
{
	FILE *inf = fopen(fname.c_str(),"w");
	if (inf == NULL) {
		fprintf(stderr,"MedLogger: init_all_files: Can't open file %s\n",fname.c_str());
		return -1;
	}

	for (int i=0; i<MAX_LOG_SECTION; i++) {
		fds[i] = inf;
	}

	return 0;
}

//-----------------------------------------------------------------------------------------------
void MedLogger::init_all_levels(int level)
{
	for (int i=0; i<MAX_LOG_SECTION; i++)
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
	if (section < MAX_LOG_SECTION)
		fds[section] = of;
}

//-----------------------------------------------------------------------------------------------
int MedLogger::init_file(int section, const string &fname)
{
	FILE *inf;

	if (section >= MAX_LOG_SECTION)
		return -2;

	inf = fopen(fname.c_str(),"w");
	if (inf == NULL) {
		fprintf(stderr,"MedLogger: init_file: Can't open file %s\n",fname.c_str());
		return -1;
	}

	fds[section] = inf;
	return 0;
}

//-----------------------------------------------------------------------------------------------
int MedLogger::log(int section, int print_level,char *fmt,...)
{
	if (section >= MAX_LOG_SECTION)
		return -2;

	if (print_level < levels[section] || fds[section] == NULL)
		return 1;

	va_list args;
	va_start(args,fmt);
	vfprintf(fds[section],fmt,args);
	va_end(args);
	fflush(fds[section]);

	return 0;
}

//-----------------------------------------------------------------------------------------------
void MedLogger::out(char *fmt,...)
{
	va_list args;
	va_start(args,fmt);
	vfprintf(out_fd,fmt,args);
	va_end(args);
	fflush(out_fd);
}
