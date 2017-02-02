//
// Utils.h - general needed and reused usitilities
//

#ifndef __UTILS__H__
#define __UTILS__H__

#include <string>
#include <vector>
#include <algorithm>

using namespace std;

#define SQUARE(x) ((x)*(x))

// files related (needed locally in order to make InfraMed independent of MedUtils

bool file_exists_IM(const string &fname);
unsigned long long get_file_size_IM(const string &fname);
int add_path_to_name_IM(const string &path, string &fname);
int add_path_to_name_IM(const string &path, vector<string> &fnames);
int copy_file_IM(const string& in_file, const string& out_file) ;
int copy_files_IM(const string &in_path, const string &out_path, vector<string>& fnames) ;
int read_bin_file_IM(string &fname, unsigned char* &data, unsigned long long &size);
int read_bin_file_IM_parallel(string &fname, unsigned char* &data, unsigned long long &size);
int write_bin_file_IM(string &fname, unsigned char* data, unsigned long long size);

// forced to keep a copy of these inside in order NOT to depend on external libraries
// for now assuming an int is enough
// All will be computed from 1/1/1900
//int date_to_days_IM(int date);
//int date_to_hours_IM(int date) { return 24*date_to_days_IM(date); };
//int date_to_minutes_IM(int date) { return 24*60*date_to_days_IM(date); }

#endif