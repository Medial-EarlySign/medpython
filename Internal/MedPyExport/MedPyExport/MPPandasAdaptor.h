#ifndef __MP_PandasAdaptor_H
#define __MP_PandasAdaptor_H

#include "MedPyCommon.h"
#ifndef SWIG
#include <iostream>
#endif // SWIG

class MPPandasAdaptor_iter;


class MPPandasAdaptor {

	std::map<std::string, std::map<int, std::string> > categories;
	class column_record {
	public:
		void* data;
		unsigned long long num_elements;
		MED_NPY_TYPE nptype;
		bool is_own;
		column_record() : data(nullptr), num_elements(0), nptype(MED_NPY_TYPES::NPY_NOTYPE), is_own(true) {}
		~column_record() {
			if (data != nullptr && is_own) {
				free(data);
				data = nullptr;
			}
		}
		int data_size() { return nptype.size_in_bytes()*num_elements; };
		void allocate_data() {
			if (data != nullptr) 
				throw runtime_error("MPPandasAdaptor:: attempt to alocate with no free"); 
			if (num_elements !=0 && nptype != MED_NPY_TYPES::NPY_NOTYPE)
				throw runtime_error("MPPandasAdaptor:: attempt to allocate but size and type are not set");
			data = malloc(data_size());
			is_own = true;
		}
	};
	std::map<std::string, column_record> columns;
	column_record& __get_col_or_throw(const string& col_name);
public:
	MPStringStringMapAdaptor type_requirements;
	MEDPY_IGNORE(void set_type_requirement(const string& col_name, const string& ctype_str));
	std::vector<std::string> keys() const { 
		vector<string> ret;
		for (auto it = columns.begin(); it != columns.end(); ++it) ret.push_back(it->first);
		return ret;
	}
	virtual ~MPPandasAdaptor() { clear(); }
	void clear() {
		columns.clear();
	}

	void export_column(const std::string& key, MEDPY_NP_VARIANT_OUTPUT(void** outarr1, unsigned long long* outarr1_sz, int* outarr1_npytype));
	void __getitem__(const std::string& key,
		MEDPY_NP_VARIANT_OUTPUT(void** outarr1, unsigned long long* outarr1_sz, int* outarr1_npytype))
	{
		export_column(key, outarr1, outarr1_sz, outarr1_npytype);
	};
	MPPandasAdaptor_iter __iter__();
	
	MPIntStringMapAdaptor get_categorial_col_dict(const std::string& col_name)
	{
		if (categories.count(col_name) == 0)
			throw runtime_error("MedPy: Not categorial column");		
		return MPIntStringMapAdaptor(&categories[col_name]);
	};
	std::vector<std::string> get_categorical_cols() {
		std::vector<std::string> ret;
		for (const auto& key : categories) {
			ret.push_back(key.first);
		}
		return ret;
	}
	void import_column(const string &col_name, void* IN_ARRAY1, unsigned long long DIM1, int NPYDTC1, bool make_a_copy = false);
#ifndef SWIG
	void push_categorial(const string& col_name, std::vector<int> index_column, std::vector<std::string> categories);
	void push_column(const string &col_name, void* arr, int arr_size, const string& ctype_str, bool make_a_copy = false) {
		import_column(col_name, arr, arr_size, MED_NPY_TYPE::ctypestr_to_npytypeid.at(ctype_str), make_a_copy);
	}
	template<typename T>
	void pull_col_as_vector(const string& col_name, vector<T>& dest) {
		if (columns[col_name].nptype.size_in_bytes() != sizeof(T)) {
			//std::cerr << "bad conversion col:nptype:sz1:sz2 - " << col_name << ":" 
			//	<< (int)columns[col_name].nptype << ":" << columns[col_name].nptype.size_in_bytes() << ":" << sizeof(T) << "\n";
			throw runtime_error("MPPandasAdaptor::Element size not compatible for column and given vector");
		}
		//std::cerr << "col:nptype:npsz:szT:num_el -> " << col_name << ":"
		//	<< (int)columns[col_name].nptype << ":" << columns[col_name].nptype.size_in_bytes() << ":" << sizeof(T) << ":"<< columns[col_name].num_elements <<"\n";
		buf_to_vector((T*)columns[col_name].data, columns[col_name].num_elements, dest);
		columns.erase(col_name);
	}
#endif //SWIG
};

class MPPandasAdaptor_iter {
	MPPandasAdaptor* obj;
	int iterator;
	std::vector<std::string> keys;
public:
	MPPandasAdaptor_iter(MPPandasAdaptor& o, std::vector<std::string> keys_param) : obj(&o), keys(keys_param) {}
	MPPandasAdaptor_iter(const MPPandasAdaptor_iter& orig) : obj(orig.obj), keys(orig.keys) {}
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
		obj->export_column(cur_key, outarr1, outarr1_sz, outarr1_npytype);
		return cur_key;
	}

	std::string __repr__() {
		return string("\"") + keys[iterator] + string("\"");
	}
};




#endif //__MP_PandasAdaptor_H
