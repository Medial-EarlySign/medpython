// Object that can be serialized and written/read from file and also initialized from string

#include "SerializableObject.h"
#include "MedProcessUtils.h"
#include <assert.h>

#define LOCAL_SECTION LOG_SRL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

// read and deserialize model
int SerializableObject::read_from_file(const string &fname) {
	unsigned char *blob;
	unsigned long long size;

	if (read_binary_data_alloc(fname, blob, size) < 0) {
		MERR("Error reading model from file %s\n", fname.c_str());
		return -1;
	}

	size_t serSize = deserialize(blob);
	assert(serSize == size);
	if (size > 0) delete[] blob;
	return 0;
}

// serialize model and write to file
int SerializableObject::write_to_file(const string &fname)
{
	unsigned char *blob;
	unsigned long long size;

	size = get_size();
	blob = new unsigned char[size];

	size_t serSize = serialize(blob);

	if (write_binary_data(fname, blob, size) < 0) {
		MERR("Error writing model to file %s\n", fname.c_str());
		return -1;
	}

	if (size > 0) delete[] blob;
	return 0;
}

// Init from string
void SerializableObject::init_from_string(string init_string) {

	map<string, string> map;
	init_map_from_string(init_string, map);
	init(map);
}