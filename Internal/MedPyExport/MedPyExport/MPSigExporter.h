#ifndef __MP_SigExporter_H
#define __MP_SigExporter_H

#include "MedPyCommon.h"
#include "MPPidRepository.h"
#ifndef AM_API_FOR_CLIENT
class MPSigExporter_iter;


class MPSigExporter {
	MEDPY_IGNORE(MedPidRepository* o);
	std::vector<std::string> data_keys;
	std::vector<void*> data_column;
	std::vector<int> data_column_nptype;
	std::map<std::string, std::vector<int> > raw_val_to_new_val;

	std::map<std::string, std::vector<std::string> > categories;
	int __get_key_id_or_throw(const string& key);
	void gen_cat_dict(const string& field_name, int channel);
	std::vector<int> pids;
	bool translate;
	std::string filter_regex;
	bool free_sig;
public:
	std::vector<std::string> keys() { return data_keys; }
	std::string sig_name;
	int sig_id = -1;
	int sig_type = -1;
	size_t record_count = 0;
	bool record_count_updated = false;
	MPSigExporter(MPPidRepository& rep, std::string signame_str, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take), int use_all_pids, int translate_flag, int free_sig_flag, std::string filter_regex_str);
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
		std::map<std::string, std::vector<std::string> >().swap(categories);
	}

	void transfer_column(const std::string& key, MEDPY_NP_VARIANT_OUTPUT(void** outarr1, unsigned long long* outarr1_sz, int* outarr1_npytype));
	void __getitem__(const std::string& key,
		MEDPY_NP_VARIANT_OUTPUT(void** outarr1, unsigned long long* outarr1_sz, int* outarr1_npytype))
	{
		transfer_column(key, outarr1, outarr1_sz, outarr1_npytype);
	};
	MPSigExporter_iter __iter__();
/*	int __len__() {
		return (int)data_column.size();
	}*/
	
	std::vector<std::string> get_categorical_field_dict(const std::string& field)
	{
		if (categories.count(field) == 0)
			throw runtime_error("MedPy: Not categorical key");		
		return categories[field];
	};

	std::vector<int> get_categorical_field_dict_int(const std::string& field)
	{
		if (raw_val_to_new_val.count(field) == 0)
			throw runtime_error("MedPy: int Not categorical key");
		return raw_val_to_new_val[field];

	};

	std::vector<std::string> get_categorical_fields() {
		std::vector<std::string> ret;
		for (const auto& key : categories) {
			ret.push_back(key.first);
		}
		return ret;
	}


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
	std::string next(MEDPY_NP_VARIANT_OUTPUT(void** outarr1, unsigned long long* outarr1_sz, int* outarr1_npytype)) {
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

#endif


#endif //__MP_SigExporter_H
