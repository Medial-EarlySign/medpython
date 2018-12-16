//
// Utils.cpp
//
#define _CRT_SECURE_NO_WARNINGS
#include <Logger/Logger/Logger.h>
#include "Utils.h"
#include <sys/stat.h>
#include <fstream>
#define _LARGE_FILE_SOURCE
#define _FILE_OFFSET_BITS 64
#include <cstdio>
#include <boost/filesystem.hpp>

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//-----------------------------------------------------------------------------
bool file_exists_IM(const string &fname)
{
	return boost::filesystem::exists(fname);
}

//----------------------------------------------------------------------------
unsigned long long get_file_size_IM(const string &fname)
{
	ifstream inf;
	unsigned long long size;

	inf.open(fname,ios::in|ios::binary);
	if (!inf)
		return 0;
	inf.seekg(0,inf.end);
	size = inf.tellg();
	inf.close();
	return size;
}

//-----------------------------------------------------------------------------
int add_path_to_name_IM(const string &path, string &fname)
{
	string temp_name ;
	if (path == "" || fname == "" || (fname.size() > 0 && fname[0] == '/'))	// if fname starts with '/' the user meant absolute path, leave it alone
		temp_name = fname ;
	else
		temp_name = path + "/" + fname;

	fname = temp_name;
	return 0;
	//return fix_path(temp_name,fname) ;
}

//-----------------------------------------------------------------------------
int add_path_to_name_IM(const string &path, vector<string> &fnames)
{
	
	string temp_name ;
	for (int i=0; i<fnames.size(); i++) {
		if (path == "" || fnames[i] == "" || (fnames[i].size() > 0 && fnames[i][0] == '/'))	// if fname starts with '/' the user meant absolute path, leave it alone
			temp_name = fnames[i] ;
		else
			temp_name = path + "/" + fnames[i];
		fnames[i] = temp_name;
		//if (fix_path(temp_name,fnames[i]) == -1)
		//	return -1 ;
	}

	return 0 ;
}

//-----------------------------------------------------------------------------
int copy_file_IM(const string& in_file, const string& out_file) {

	MLOG_D("InfraMed: copying %s to %s\n",in_file.c_str(),out_file.c_str()) ;

	ifstream src(in_file,ios::binary) ;
	ofstream dst(out_file,ios::binary) ;

	if (! (src.is_open() && dst.is_open()))
		MTHROW_AND_ERR("failed to copy [%s] to [%s]\n", in_file.c_str(), out_file.c_str());

	dst << src.rdbuf() ;

	src.close() ;
	dst.close() ;
	return 0 ;
}

//-----------------------------------------------------------------------------
int copy_files_IM(const string &in_path, const string &out_path, vector<string>& fnames) {

	if (in_path != out_path) {
		for (unsigned int i=0; i<fnames.size(); i++) {
			if (fnames[i].size() > 0 && fnames[i][0] == '/') {
				MLOG("not copying [%s] as its absolute path\n", fnames[i].c_str());
				continue;
			}
			string temp_in_file = in_path + "/" + fnames[i] ;
			string temp_out_file = out_path + "/" + fnames[i] ;

			string in_file,out_file ;
			//if (fix_path(temp_in_file,in_file) == -1 || fix_path(temp_out_file,out_file) == -1)
			//	return -1 ;
			in_file = temp_in_file;
			out_file = temp_out_file;

			if (copy_file_IM(in_file,out_file) < 0)
				return -1 ;
		}
	}
	return 0 ;
}

//-----------------------------------------------------------------------------
int read_bin_file_IM(string &fname, unsigned char* &data, unsigned long long &size)
{
	ifstream inf;
	FILE *inf_read; // faster reading with FILE * ....

	inf.open(fname, ios::in | ios::binary | ios::ate);

	if (!inf) {
		MERR("read_binary_data_alloc(): can't open file %s for read\n", fname.c_str());
		return -1;
	}

	size = inf.tellg();
	data = new unsigned char[size];
	if (data == NULL)
		MTHROW_AND_ERR("read_bin_file_IM could not allocate %zd bytes\n", size);
	//MLOG("allocated data %d with size = %d\n", data, size);
//	inf.seekg(0, ios::beg);
//	inf.read((char *)data, size);
	inf.close();

	inf_read = fopen(fname.c_str(), "rb");
	if (fread(data, 1, size, inf_read) != size)
		return -1;

	fclose(inf_read);

	return 0;
}

//-----------------------------------------------------------------------------
// testing ideas for parallel read
int read_bin_file_IM_parallel(string &fname, unsigned char* &data, unsigned long long &size)
{
	ifstream inf;
	//FILE *inf_read; // faster reading with FILE * ....

	inf.open(fname, ios::in | ios::binary | ios::ate);

	if (!inf) {
		MERR("read_binary_data_alloc(): can't open file %s for read\n", fname.c_str());
		return -1;
	}

	size = inf.tellg();

	if (size == 0) {
		data = NULL;
		return 0;
	}

	data = new unsigned char[size];
	//	inf.seekg(0, ios::beg);
	//	inf.read((char *)data, size);
	inf.close();


	// now splitting read mission to blocks (def 1MB)
	int k = 23;
	size_t mask = (1<<k)-1;
	size_t n_blocks = size >> k;
	int nerr=0;

#pragma omp parallel for 
	for (long long i=0; i<=(long long)n_blocks; i++) {
		FILE *inf_read;
		inf_read = fopen(fname.c_str(), "rb");
		if (inf_read == NULL)
			MERR("Can't open file %s\n", fname.c_str());
		size_t nread = ((size_t)1<<k);
		if (i==n_blocks)
			nread = size & mask;
#if defined (_MSC_VER) || defined (_WIN32)
		_fseeki64(inf_read, i<<k, SEEK_SET);
#else
		fseeko64(inf_read, i<<k, SEEK_SET);
#endif
		if (nread > 0)
			if (fread(&data[i<<k], 1, nread, inf_read) != nread)
				nerr++;

		fclose(inf_read);
	}

	if (nerr > 0)
		return -1;


	return 0;
}

//-----------------------------------------------------------------------------
int write_bin_file_IM(string &fname, unsigned char* data, unsigned long long size)
{
	ofstream of;

	of.open(fname, ios::out | ios::binary);

	if (!of) {
		MERR("InfraMed Utils:: write_bin_file(): can't open file %s for write\n", fname.c_str());
		return -1;
	}

	of.write((char *)data, size);

	of.close();

	return 0;
}
