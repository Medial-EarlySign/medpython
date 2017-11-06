//
// MedIO
//

#include "MedIO.h"
#include <sys/stat.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <boost/crc.hpp>

#define LOCAL_SECTION LOG_MEDIO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//-----------------------------------------------------------------------------
bool file_exists(const string &fname)
{
	struct stat fstat;
	return (stat(fname.c_str(), &fstat) == 0);
}

//----------------------------------------------------------------------------
unsigned long long get_file_size(const string &fname)
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
void add_path_to_name(const string &path, string &fname)
{
	if (path != "") {
		fname = path + "/" + fname;
	}
}

//-----------------------------------------------------------------------------
void add_path_to_name(const string &path, vector<string> &fnames)
{
	if (path != "") {
		for (int i=0; i<fnames.size(); i++)
			fnames[i] = path + "/" + fnames[i];
	}
}

//-----------------------------------------------------------------------------
int copy_file(const string& in_file, const string& out_file) {

	MLOG_D("InfraMed: copying %s to %s\n",in_file.c_str(),out_file.c_str()) ;

	ifstream src(in_file,ios::binary) ;
	ofstream dst(out_file,ios::binary) ;

	if (! (src.is_open() && dst.is_open()))
		return -1 ;

	dst << src.rdbuf() ;

	src.close() ;
	dst.close() ;
	return 0 ;
}

//-----------------------------------------------------------------------------
int copy_files(const string &in_path, const string &out_path, vector<string>& fnames) {

	if (in_path != out_path) {
		for (unsigned int i=0; i<fnames.size(); i++) {
			string in_file = in_path + "/" + fnames[i] ;
			string out_file = out_path + "/" + fnames[i] ;
			if (copy_file(in_file,out_file) < 0)
				return -1 ;
		}
	}
	return 0 ;
}

//-----------------------------------------------------------------------------
int fix_path(const string& in, string& out) {
	// fprintf(stderr, "Converting path \'%s\'\n", in.c_str());
	// fflush(stderr);
	map<string, string> folders ;
	folders["W"] = "Work" ;
	folders["U"] = "UsersData" ;
	folders["X"] = "Temp" ;
	folders["P"] = "Products" ;
	folders["T"] = "Data" ;
	// fprintf(stderr, "Initialized network drive table\n");

#ifndef _WIN32
	// on Linux, handle first the Windows native format: \\\\server\\Work\..  
	if (in.length() >= 2 && in.substr(0, 2) =="\\\\") {
		// just switching '\' to '/'; works, but adjacent slashes should be unified
		out = in;
		char revSlash = '\\';
		char fwdSlash = '/';
		std::replace(out.begin(), out.end(), revSlash, fwdSlash);
		fprintf(stderr, "Converted path \'%s\' to \'%s\'\n", in.c_str(), out.c_str());
		fflush(stderr);

		return 0;
	}
#endif

	// Work only on "X:/...." or  "/cygdrive/X/...." input strings
	if ((in.length() < 3 || in.substr(1,2) != ":/") && 
		(in.length() < 12 || in.substr(0, 10) != "/cygdrive/" || in.substr(11, 1) != "/")) {
		out = in ;
		return 0 ;
	}

	char driveLetter = (in.substr(1, 2) == ":/") ? in.substr(0, 1)[0] : (char)toupper(in.substr(10, 1)[0]);
	string drive = string(1, driveLetter);

	if (drive=="C" || drive=="D" || drive=="H") {
		out = in ;
		return 0 ;
	}

	int pathPos = (in.substr(1, 2) == ":/") ? 3 : 12;

	if (folders.find(drive) == folders.end()) {
		fprintf(stderr, "Unknown Folder Map %s\n", drive.c_str()) ;
		return -1 ;
	}

#ifdef _WIN32
//	out = "\\\\server\\" + folders[drive];
	out = "\\\\nas1\\" + folders[drive];
#else
	out = "/server/" + folders[drive] ;
#endif

	istringstream in_stream(in.substr(pathPos, in.length() - pathPos)) ;	
	string token ;

	while (getline(in_stream,token,'/'))
#ifdef _WIN32
		out += "\\" + token ;
#else
		out += "/" + token ;
#endif
	fprintf(stderr, "Converted path \'%s\' to \'%s\'\n", in.c_str(), out.c_str());
	fflush(stderr);

	return 0 ;

}

//-----------------------------------------------------------------------------
int write_string(const string &fname, string &data)
{
	ofstream out(fname);
	out << data;
	out.close();
	return(0);
}


//-----------------------------------------------------------------------------
int write_binary_data(const string &fname, unsigned char *data, unsigned long long size)
{
	ofstream of;

	//MLOG("Writing file %s :: size %lld\n", fname.c_str(), size);
	of.open(fname,ios::out | ios::binary);

	if (!of) {
		MERR("write_binary_data(): can't open file %s for write\n",fname.c_str());
		return -1;
	}

	of.write((char *)data,size);

	of.close();

	return 0;
}

//-----------------------------------------------------------------------------
int read_binary_data_alloc(const string &fname, unsigned char *&data, unsigned long long &size)
{
	ifstream inf;

	inf.open(fname,ios::in | ios::binary | ios::ate);

	if (!inf) {
		MERR("read_binary_data_alloc(): can't open file %s for read\n",fname.c_str());
		return -1;
	}

	size = inf.tellg();
	data = new unsigned char [size];
	inf.seekg(0,ios::beg);
	inf.read((char *)data,size);
	
	boost::crc_32_type checksum_agent;
	checksum_agent.process_bytes(data, size);
	MLOG("read_binary_data_alloc [%s] with crc32 [%d]\n", fname.c_str(), checksum_agent.checksum());

	inf.close();
	return 0;
}

//-----------------------------------------------------------------------------
int read_binary_data_prealloc(const string &fname, unsigned char *&data, unsigned long long &size, unsigned long long max_size)
{
	ifstream inf;

	inf.open(fname,ios::in | ios::binary | ios::ate);

	if (!inf) {
		MERR("read_binary_data_prealloc(): can't open file %s for read\n",fname.c_str());
		return -1;
	}

	size = inf.tellg();
	if (size > max_size) {
		MERR("read_binary_data_prealloc(): Not enough space in *data (%ld needed while max is %ld)\n",size,max_size);
		inf.close();
		return -1;
	}

	inf.seekg(0,ios::beg);
	inf.read((char *)data,size);
	inf.close();
	return 0;
}

