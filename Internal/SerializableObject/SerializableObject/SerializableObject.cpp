// Object that can be serialized and written/read from file and also initialized from string

#include "SerializableObject.h"
#include <assert.h>
#include <boost/crc.hpp>
#include <chrono>
#include <thread>
#include <boost/algorithm/string/replace.hpp>


#ifndef _MSC_VER
#include <unistd.h>
#endif

#define LOCAL_SECTION LOG_SRL
#define LOCAL_LEVEL	LOG_DEF_LEVEL


float med_stof(const string& _Str) {
	try {
		return stof(_Str);
	}
	catch (exception e) {
		MTHROW_AND_ERR("invalid stof argument [%s]\n", _Str.c_str());
	}
}

int med_stoi(const string& _Str) {
	try {
		return stoi(_Str);
	}
	catch (exception e) {
		MTHROW_AND_ERR("invalid stoi argument [%s]\n", _Str.c_str());
	}
}


void SerializableObject::_read_from_file(const string &fname, bool throw_on_version_error) {
	unsigned char *blob;
	unsigned long long final_size;

	if (MedSerialize::read_binary_data_alloc(fname, blob, final_size) < 0)
		MTHROW_AND_ERR("Error reading model from file %s\n", fname.c_str());

	boost::crc_32_type checksum_agent;
	checksum_agent.process_bytes(blob, final_size);
	MLOG("read_from_file [%s] with crc32 [%d] and size [%ld]\n", fname.c_str(), checksum_agent.checksum(), final_size);

	int vers = *((int*)blob);
	if (vers != version()) {
		if (throw_on_version_error) {
			MTHROW_AND_ERR("deserialization error. code version %d. requested file version %d\n",
				version(), vers);
		}
		else {
			MWARN("WARNING: SerializableObject::read_from_file - code version %d. requested file version %d\n",
				version(), vers);
		}
	}
	unsigned char *blob_without_version = blob + sizeof(int);

	size_t serSize = deserialize(blob_without_version);
	if (serSize + sizeof(int) != final_size)
		MTHROW_AND_ERR("final_size=%lld, serSize=%d\n", final_size, (int)serSize);
	if (final_size > 0) delete[] blob;
}

//read unsafe without checking version:
int SerializableObject::read_from_file_unsafe(const string &fname) {
	_read_from_file(fname, false);
	return 0;
}

// read and deserialize model
int SerializableObject::read_from_file(const string &fname) {
	_read_from_file(fname, true);
	return 0;
}

// serialize model and write to file
int SerializableObject::write_to_file(const string &fname)
{
	unsigned char *blob;
	size_t size;

	size = get_size();

	blob = new unsigned char[size + sizeof(int)];
	*((int*)blob) = version(); //save version
	size_t serSize = serialize(blob + sizeof(int));
	if (size != serSize)
		MTHROW_AND_ERR("size=%lld, serSize=%d\n", size, (int)serSize);

	size_t final_size = serSize + sizeof(int);

	boost::crc_32_type checksum_agent;
	checksum_agent.process_bytes(blob, final_size);
	MLOG("write_to_file [%s] with crc32 [%d]\n", fname.c_str(), checksum_agent.checksum());

	if (MedSerialize::write_binary_data(fname, blob, final_size) < 0) {
		MERR("Error writing model to file %s\n", fname.c_str());
		return -1;
	}

	if (size > 0) delete[] blob;
	return 0;
}

// Init from string
int SerializableObject::init_from_string(string init_string) {

	map<string, string> map;
	if (MedSerialize::init_map_from_string(init_string, map) < 0) return -1;
	if (map.size() == 1 && map.begin()->first == "pFile") {
		return init_params_from_file(map.begin()->second);
	}
	if (init(map) < 0) return -1;

	return 0;
}

// Init from file
int SerializableObject::init_params_from_file(string fname)
{
	string data;
	if (MedSerialize::read_file_into_string(fname, data) < 0) return -1;
	boost::replace_all(data, "\n", "");
	return init_from_string(data);
}

