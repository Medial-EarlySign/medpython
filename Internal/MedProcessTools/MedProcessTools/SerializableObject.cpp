// Object that can be serialized and written/read from file and also initialized from string

#include "SerializableObject.h"
#include "MedProcessUtils.h"
#include <assert.h>
#include <boost/crc.hpp>
#include <chrono>
#include <thread>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#define LOCAL_SECTION LOG_SRL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

// read and deserialize model
int SerializableObject::read_from_file(const string &fname) {
	int attempts = 0;
	unsigned char *blob;
	unsigned long long final_size;
	for (;;) {
		try {
			if (read_binary_data_alloc(fname, blob, final_size) < 0) {
				MERR("Error reading model from file %s\n", fname.c_str());
				return -1;
			}

			boost::crc_32_type checksum_agent;
			checksum_agent.process_bytes(blob, final_size);
			MLOG("read_from_file [%s] with crc32 [%d]\n", fname.c_str(), checksum_agent.checksum());

			int vers = *((int*)blob);
			if (vers != version())
				MTHROW_AND_ERR("deserialization error. code version %d. requested file version %d\n",
					version(), vers);
			unsigned char *blob_without_version = blob + sizeof(int);

			size_t serSize = deserialize(blob_without_version);
			if (serSize + sizeof(int) != final_size)
				MTHROW_AND_ERR("final_size=%d, serSize=%d\n", final_size, (int)serSize);
			if (final_size > 0) delete[] blob;
			return 0;
		}
		catch (exception e) {
			if (attempts++ >= 10)
				throw e;
			MWARN("[%d] attempt to read [%s] failed with [%s], retyring...\n", attempts, fname.c_str(), e.what());
			std::this_thread::sleep_for(chrono::seconds(10));
		}
	}
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
		MTHROW_AND_ERR("size=%d, serSize=%d\n", size, (int)serSize);

	size_t final_size = serSize + sizeof(int);

	boost::crc_32_type checksum_agent;
	checksum_agent.process_bytes(blob, final_size);
	MLOG("write_to_file [%s] with crc32 [%d]\n", fname.c_str(), checksum_agent.checksum());

	if (write_binary_data(fname, blob, final_size) < 0) {
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