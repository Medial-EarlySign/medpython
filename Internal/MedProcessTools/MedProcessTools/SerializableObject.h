// Object that can be serialized and written/read from file

#ifndef _SERIALIZABLE_OBJECT_H_
#define _SERIALIZABLE_OBJECT_H_

#include "Logger/Logger/Logger.h"
#include "MedUtils/MedUtils/MedIO.h"

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

#endif