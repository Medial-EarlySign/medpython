// Object that can be serialized and written/read from file and also initialized from string

#include "SerializableObject.h"
#include <assert.h>
#include <boost/crc.hpp>
#include <chrono>
#include <thread>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>
#include "MedIO/MedIO/MedIO.h"
#include "MedUtils/MedUtils/MedRunPath.h"
#include <cctype>
#include <cmath>


#ifndef _MSC_VER
#include <unistd.h>
#endif

#define LOCAL_SECTION LOG_SRL
#define LOCAL_LEVEL	LOG_DEF_LEVEL


float med_stof(const string& _Str) {
	try {
		return stof(_Str);
	}
	catch (exception &e) {
		MTHROW_AND_ERR("invalid stof argument [%s]\n", _Str.c_str());
	}
}

int med_stoi(const string& _Str) {
	try {
		return stoi(_Str);
	}
	catch (exception &e) {
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
			if (abs(vers - version()) <= 3) {
				MTHROW_AND_ERR("deserialization error of %s from %s. code version %d. requested file version %d\n",
					my_class_name().c_str(), fname.c_str(), version(), vers);
			}
			else {
				MTHROW_AND_ERR("deserialization error of %s from %s. Are you sure this is correct file path? Please check the file path\n",
					my_class_name().c_str(), fname.c_str());
			}
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
		MTHROW_AND_ERR("size=%zu, serSize=%d\n", size, (int)serSize);

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
	if (MedSerialize::init_map_from_string(init_string, map) < 0)
		MTHROW_AND_ERR("Error Init from String %s\n", init_string.c_str());

	if (map.size() == 1 && map.begin()->first == "pFile") {
		int rc = init_params_from_file(map.begin()->second);
		if (rc < 0)
			MTHROW_AND_ERR("Error Init params from file %s\n", map.begin()->second.c_str());
	}

	for (auto &e : map) {
		string val = e.second;
		boost::to_upper(val);
		if (val.compare(0, 5, "FILE:") == 0 || val.compare(0, 5, "LIST:") == 0 ||
			val.compare(0, 9, "LIST_REL:") == 0) {
			string param;
			if (init_param_from_file(e.second, param) < 0)
				MTHROW_AND_ERR("Error Init params from file %s\n", e.second.c_str());
			e.second = param;
		}
	}

	if (init(map) < 0)
		MTHROW_AND_ERR("Error Init from string after convertion to map %s\n", init_string.c_str());

	return 0;
}
int SerializableObject::update_from_string(const string &init_string) {

	map<string, string> mapper;
	if (MedSerialize::init_map_from_string(init_string, mapper) < 0)
		MTHROW_AND_ERR("Error Init from String %s\n", init_string.c_str());

	if (mapper.size() == 1 && mapper.begin()->first == "pFile") {
		int rc = init_params_from_file(mapper.begin()->second);
		if (rc < 0)
			MTHROW_AND_ERR("Error Init params from file %s\n", mapper.begin()->second.c_str());
	}

	for (auto &e : mapper) {
		string val = e.second;
		boost::to_upper(val);
		if (val.compare(0, 5, "FILE:") == 0 || val.compare(0, 5, "LIST:") == 0 ||
			val.compare(0, 9, "LIST_REL:") == 0) {
			string param;
			if (init_param_from_file(e.second, param) < 0)
				MTHROW_AND_ERR("Error Init params from file %s\n", e.second.c_str());
			e.second = param;
		}
	}

	if (update(mapper) < 0)
		MTHROW_AND_ERR("Error Init from string after convertion to map %s\n", init_string.c_str());

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

// Init a specific param from a file
int SerializableObject::init_param_from_file(string file_str, string &param)
{
	string upper_cp = boost::to_upper_copy(file_str);

	// prefix is FILE: as file: is reserved for medmodel json usages
	if (upper_cp.compare(0, 5, "FILE:") == 0) {
		string fname = file_str.substr(5);
		if (MedSerialize::read_file_into_string(fname, param) < 0) return -1;
	}

	if (upper_cp.compare(0, 9, "LIST_REL:") == 0) {
		string fname;
		if (!run_current_path.empty())
			fname = run_current_path + path_sep() + file_str.substr(9);
		else
			fname = file_str.substr(9);
		if (MedSerialize::read_list_into_string(fname, param) < 0) return -1;
	}

	if (upper_cp.compare(0, 5, "LIST:") == 0) {
		string fname = file_str.substr(5);
		if (MedSerialize::read_list_into_string(fname, param) < 0) return -1;
	}

	return 0;
}

string SerializableObject::object_json() const {
	stringstream str;

	str << "{\n\t\"Object\":\"" << my_class_name() << "\",\n\t\"Version\":" << version();
	str << ",\n\t\"data\": ";
	vector<string> field_names;
	serialized_fields_name(field_names);
	//try print field name => to field value recursively:
	str << "{";
	for (int i = 0; i < field_names.size(); ++i) {
		if (i > 0)
			str << ",";
		str << "\n\t\"" << field_names[i] << "\":\"NOT_IMPLEMENTED\"";
	}
	str << "}";
	//close json:
	str << "\n}";

	return str.str();
}

void mes_trim(string &s) {
	auto p_start = s.begin();
	auto p_end = s.end() - 1;
	//trim start
	while (p_start <= p_end && isspace(*p_start))
		++p_start;
	//trim end
	while (p_start <= p_end && isspace(*p_end))
		--p_end;

	s = string(p_start, p_end + 1);
}