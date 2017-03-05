// Object that can be serialized and written/read from file

#ifndef _SERIALIZABLE_OBJECT_H_
#define _SERIALIZABLE_OBJECT_H_

#include "Logger/Logger/Logger.h"
#include "MedUtils/MedUtils/MedIO.h"
#include <cstring>

using namespace std;

class SerializableObject {
public:

	// Virtual serialization
	virtual size_t get_size() { return 0; }
	virtual size_t serialize(unsigned char *blob) { return 0; }
	virtual size_t deserialize(unsigned char *blob) { return 0; }

	// read and deserialize model
	virtual int read_from_file(const string &fname);

	// serialize model and write to file
	virtual int write_to_file(const string &fname);

	// Init from string
	void init_from_string(string init_string);
	virtual int  init(map<string, string>& map) { return 0; }
};


// helpers to common serialization
//namespace MedSerialize {
//
//	// templated simple ones : int, float, long, double, etc...
//	template <class T> size_t get_size(T &elem);
//	template <class T> size_t serialize(unsigned char *blob, T &elem);
//	template <class T> size_t deserialize(unsigned char *blob, T &elem);
//
//	// string is special
//	template<> size_t get_size(string &str);
//	template<> size_t serialize(unsigned char *blob, string &str);
//	template<> size_t deserialize(unsigned char *blob, string &str);
//
//	// vector of type T that has a MedSerialize function
//	template <class T> size_t get_size(vector<T> &v);
//	template <class T> size_t serialize(unsigned char *blob, vector<T> &v);
//	template <class T> size_t deserialize(unsigned char *blob, vector<T> &v);
//
//}

#include "SerializableObject_imp.h"

#endif