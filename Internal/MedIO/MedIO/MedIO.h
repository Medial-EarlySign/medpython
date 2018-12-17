//
// MedIO : tools to handle read/write from files
// 

#ifndef __MED_IO_LIB_H__
#define __MED_IO_LIB_H__

#include <Logger/Logger/Logger.h>
#include <map>

using namespace::std;

//============================================================
// Files/IO related
//============================================================
bool file_exists(const string &fname);
unsigned long long get_file_size(const string &fname);
int copy_file(const string& in_file, const string& out_file);
int copy_files(const string &in_path, const string &out_path, vector<string>& fnames);

//=============================================================
// paths related
//=============================================================
void add_path_to_name(const string &path, string &fname);
void add_path_to_name(const string &path, vector<string> &fnames);
// General function for platform-independent path names.
int fix_path(const string& in, string& out);

/// Fixes filename and replaces forbidden charecters with replace_char
string fix_filename_chars(string* s, char replace_char = ' ');
//=============================================================
// read/write vectors and arrays
//=============================================================
// write vector of type T to file
template <class T> int write_vector(const string &fname, vector<T> &data);
// read vector of type T from file
template <class T> int read_vector(const string &fname, vector<T> &data);
// read vector of type T from file, starting read at position start_pos
template <class T> int read_vector(const string &fname, unsigned long long start_pos, vector<T> &data);
// write array of type T to file, len is given is units of <T>
template <class T> int write_vector(const string &fname, T *data, unsigned long long len);
// read array of type T from file , len is returned in units of <T>. If data in NULL allocates space, otherwise assumes there's enough allocated space.
template <class T> int read_vector(const string &fname, T *&data, unsigned long long &len);
// read array of type T from file , len is returned in units of <T>, starting read at position start_pos. If data in NULL allocates space, otherwise assumes there's enough allocated space.
template <class T> int read_vector(const string &fname, unsigned long long start_pos, T *&data, unsigned long long &len);

// write string to file (not appending !!)
int write_string(const string &fname, string &data);

// read file into a string
int read_file_into_string(const string &fname, string &data);

// write binary data to file (not appending !!)
int write_binary_data(const string &fname, unsigned char *data, unsigned long long size);
int read_binary_data_alloc(const string &fname, unsigned char *&data, unsigned long long &size);
int read_binary_data_prealloc(const string &fname, unsigned char *&data, unsigned long long &size, unsigned long long max_size);


// read text file, parse it and get the i-th column string
int read_text_file_col(string fname, string ignore_pref, string separators, int col_idx, vector<string> &res);
int read_text_file_cols(string fname, string separators, vector<vector<string>> &res);

//=============================================================
// read/write serialized classes
//=============================================================
// following assume class T has a serialize() and deserialize() methods.
// serialize() should return length (in chars) of the created blob (or -1 if failed)
// deserialize() should return the number of bytes consumed (or -1 if failed)

// write vector of serializable type T to file
template <class T> int write_serialized_vector(const string &fname, vector<T> &data);
// read vector of deserializable type T from file
template <class T> int read_deserialized_vector(const string &fname, vector<T> &data);

char path_sep();

#include "MedIO_imp.h"

#endif