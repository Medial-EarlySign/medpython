#ifndef __MP_SigExporter_H
#define __MP_SigExporter_H

#include "MedPyCommon.h"
#include "MPPidRepository.h"

class MPSigExporter_iter;


class MPSigExporter {
	MedPidRepository* o;
	std::vector<std::string> data_keys;
	std::vector<void*> data_column;
	std::vector<int> data_column_nptype;
public:
	std::vector<std::string> keys() { return data_keys; }
	std::string sig_name;
	int sig_id = -1;
	int sig_type = -1;
	int record_count = -1;
	MPSigExporter(MPPidRepository& rep, std::string signame_str) : o(rep.o), sig_name(signame_str) {
		if (rep.loadsig(signame_str) != 0)
			throw runtime_error("could not load signal");
		sig_id = rep.sig_id(sig_name);
		if (sig_id == -1)
			throw runtime_error("bad sig id");
		sig_type = rep.sig_type(sig_name);
		update_record_count();
		get_all_data();
	}
	void update_record_count();
	void get_all_data();
	void clear() {
		for (void* ptr : data_column) {
			if (ptr != nullptr)
				free(ptr);
		}
		vector<string>().swap(data_keys);
		vector<void*>().swap(data_column);
		vector<int>().swap(data_column_nptype);
	}

	void transfer_column(const std::string& key, MEDPY_NP_VARIANT_OUTPUT(void** outarr1, int* outarr1_sz, int* outarr1_npytype));
	void __getitem__(const std::string& key,
		MEDPY_NP_VARIANT_OUTPUT(void** outarr1, int* outarr1_sz, int* outarr1_npytype))
	{
		return transfer_column(key, outarr1, outarr1_sz, outarr1_npytype);
	};
	MPSigExporter_iter __iter__();
};

class MPSigExporter_iter {
	MPSigExporter* obj;
	int iterator;
	std::vector<std::string> keys;
public:
	MPSigExporter_iter(MPSigExporter& o, std::vector<std::string> keys_param) : obj(&o), keys(keys_param) {}
	MPSigExporter_iter(const MPSigExporter_iter& orig) : obj(orig.obj), keys(orig.keys) {}
	//The return type of both string and the NumPy outarr will result in a [str,outarr] list which is good 
	//in this case because that makes it convertible to dict easily.
	std::string next(MEDPY_NP_VARIANT_OUTPUT(void** outarr1, int* outarr1_sz, int* outarr1_npytype)) {
		if (iterator >= keys.size()) {
			obj->clear();
			throw StopIterator();
		}
		string cur_key = keys[iterator];
		//advance:
		iterator++;
		//return values
		obj->transfer_column(cur_key, outarr1, outarr1_sz, outarr1_npytype);
		return cur_key;
	}

	std::string __repr__() {
		return string("\"") + keys[iterator] + string("\"");
	}
};




#endif //__MP_SigExporter_H
