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

	int vers = *((int*)blob);
	int s = sizeof(int);
	if (vers != version())
		MTHROW_AND_ERR("deserialization error. code version %d. requested file version %d\n",
			version(), vers);
	unsigned char *blob_without_version = blob + s;

	size_t serSize = deserialize(blob_without_version);
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

	blob = new unsigned char[size + sizeof(int)];
	*((int*)blob) = version(); //save version
	size_t serSize = serialize(blob + sizeof(int));
	size_t final_size = serSize + sizeof(int);
	if (size != serSize)
		MLOG("get_size=%d, acctual_szie=%d\n", size, (int)serSize);
	assert(size == serSize);

	if (write_binary_data(fname, blob, size) < 0) {
		MERR("Error writing model to file %s\n", fname.c_str());
		return -1;
	}

	if (size > 0) delete[] blob;
	return 0;
}

// Init from string
int SerializableObject::init_from_string(string init_string) {

	map<string, string> map;
	if (init_map_from_string(init_string, map) < 0) return -1;
	if (init(map) < 0) return -1;

	return 0;
}