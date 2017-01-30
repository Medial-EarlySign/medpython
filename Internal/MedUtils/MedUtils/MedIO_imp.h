// implementation for template functions in MedIO.h

#include <fstream>
//#define LOCAL_SECTION LOG_MEDIO
//#define LOCAL_LEVEL	LOG_DEF_LEVEL
//extern MedLogger global_logger;

//-----------------------------------------------------------------------------
template <class T> int write_vector(const string &fname, vector<T> &data)
{
	ofstream of;

//	MLOG_D("MedIO: write_vector: fname %s data size %d\n",fname.c_str(),(int)data.size());
	of.open(fname,ios::out|ios::binary);
	if (!of)
		return -1;

	if (data.size() > 0) {
		char *d = (char *)(&data[0]);
//		MLOG_D("MedIO: writing %d bytes to file %s\n",sizeof(T)*data.size(),fname.c_str());
		of.write(d,sizeof(T)*data.size());
	}

	of.close();
	return 0;
}

//-----------------------------------------------------------------------------
template <class T> int read_vector(const string &fname, unsigned long long start_pos, vector<T> &data)
{
	ifstream inf;
	unsigned long long size;

	data.clear();
	size = get_file_size(fname);

	if (start_pos >= size)
		return 0;

	inf.open(fname,ios::in|ios::binary);
	if (!inf)
		return -1;

	inf.seekg(start_pos,ios_base::beg);

	size = size - start_pos;
	if (size > 0) {
		data.resize(size/(sizeof(T)));

		char *d = (char *)&data[0];
		inf.read(d,size);
	}
	return 0;
}


//-----------------------------------------------------------------------------
template <class T> int read_vector(const string &fname, vector<T> &data)
{
	return(read_vector(fname,data,0));
}
