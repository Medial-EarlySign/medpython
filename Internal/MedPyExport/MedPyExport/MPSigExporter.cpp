#include "MPSigExporter.h"

#include <time.h>
#include <string>

#include "InfraMed/InfraMed/MedConvert.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/Utils.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedModel.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"


MPSigExporter::MPSigExporter(MPPidRepository& rep, std::string signame_str, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take), int use_all_pids, int translate_flag, int free_sig_flag) 
	: o(rep.o), sig_name(signame_str), translate(translate_flag!=0), free_sig(free_sig_flag!=0) {
	if (rep.loadsig(signame_str) != 0)
		throw runtime_error("could not load signal");
	sig_id = rep.sig_id(sig_name);
	if (sig_id == -1)
		throw runtime_error("bad sig id");
	sig_type = rep.sig_type(sig_name);
	if (use_all_pids)
		this->pids = o->all_pids_list;
	else if (num_pids_to_take == 0)
		this->pids = o->pids;
	else
		buf_to_vector(pids_to_take, num_pids_to_take, this->pids);

	update_record_count();
	get_all_data();
	if (free_sig)
		rep.free(signame_str);
}


static int convert_sv_type_to_npy_val_type(int sv_type){
	switch (sv_type) {
	case GenericSigVec::type_enc::INT32:   return (int)MED_NPY_TYPE::values::NPY_INT;  //int
	case GenericSigVec::type_enc::INT64:   return (int)MED_NPY_TYPE::values::NPY_LONGLONG;  //long long
	case GenericSigVec::type_enc::UINT16:  return (int)MED_NPY_TYPE::values::NPY_USHORT;  //unsigned short
	case GenericSigVec::type_enc::UINT8:   return (int)MED_NPY_TYPE::values::NPY_UBYTE;  //unsigned char
	case GenericSigVec::type_enc::UINT32:  return (int)MED_NPY_TYPE::values::NPY_UINT;  //unsigned int
	case GenericSigVec::type_enc::UINT64:  return (int)MED_NPY_TYPE::values::NPY_ULONGLONG;  //unsigned long long
	case GenericSigVec::type_enc::INT8:    return (int)MED_NPY_TYPE::values::NPY_CHAR;  //char
	case GenericSigVec::type_enc::INT16:   return (int)MED_NPY_TYPE::values::NPY_SHORT;  //short
	case GenericSigVec::type_enc::FLOAT32: return (int)MED_NPY_TYPE::values::NPY_FLOAT;  //float
	case GenericSigVec::type_enc::FLOAT64: return (int)MED_NPY_TYPE::values::NPY_DOUBLE;  //double
	case GenericSigVec::type_enc::FLOAT80: return (int)MED_NPY_TYPE::values::NPY_LONGDOUBLE;  //long double
	}
	return (int)MED_NPY_TYPE::values::NPY_NOTYPE;
}


int MPSigExporter::__get_key_id_or_throw(const string& key) {
	for (int i = 0; i < data_keys.size();++i)
		if (data_keys[i] == key) return i;
	throw runtime_error("Unknown row");
}

void MPSigExporter::gen_cat_dict(const string& field_name, int channel) {
	if (!translate) return;
	int key_index = __get_key_id_or_throw(field_name);
	if (!o->sigs.is_categorical_channel(sig_id, channel))
		return;
	int section_id = o->dict.section_id(sig_name);
	void* arr = data_column[key_index];
	size_t arr_sz = this->record_count;
	int* new_arr = nullptr;
	int arr_npytype = data_column_nptype[key_index];
	if (arr_npytype != (int)MED_NPY_TYPES::NPY_INT)
		new_arr = (int*)malloc(sizeof(int)*arr_sz);
	//std::unordered_set<int> values;
	switch (arr_npytype) {
	case (int)MED_NPY_TYPES::NPY_FLOAT: 
		{float* tarr = (float*)arr; for (size_t i = 0; i < arr_sz; ++i) new_arr[i]=(int)tarr[i]; }
		break;
	case (int)MED_NPY_TYPES::NPY_USHORT: 
		{unsigned short* tarr = (unsigned short*)arr; for (size_t i = 0; i < arr_sz; ++i) new_arr[i] = (int)tarr[i]; }
		break;
	case (int)MED_NPY_TYPES::NPY_LONGLONG: 
		{long long* tarr = (long long*)arr; for (size_t i = 0; i < arr_sz; ++i) new_arr[i] = (int)tarr[i]; }
		break;
	case (int)MED_NPY_TYPES::NPY_SHORT: 
		{short* tarr = (short*)arr; for (size_t i = 0; i < arr_sz; ++i) new_arr[i] = (int)tarr[i]; }
		break;
	case (int)MED_NPY_TYPES::NPY_INT:
		{new_arr = (int*)arr; }
		break;
	default:
		if (new_arr) free(new_arr);
		throw runtime_error("MedPy: categorical value type not supported, we only have values of types float, unsigned short, long long, short");
		break;
	}
	auto& Id2Names = o->dict.dict(section_id)->Id2Names;
	std::unordered_map<int, int> translation_dict;
	translation_dict.reserve(Id2Names.size());
	raw_val_to_new_val[field_name].reserve(Id2Names.size());
	std::vector<std::string> category;
	category.push_back("Undefined Category");  // category[0] , (code 0) is undefined
	for (size_t i = 0; i < arr_sz;i++) {
		int raw_val = new_arr[i];
		if (translation_dict.count(raw_val) == 0) {
			do {
				if (!Id2Names.count(raw_val)) {
					translation_dict[raw_val] = 0;
					break;
				}
				auto& names = Id2Names[raw_val];
				if (names.size() == 0) { 
					translation_dict[raw_val] = 0; 
					break;
				}
				string cat_name = "";
				cat_name = names[0];
				for (int j = 1; j < names.size(); j++)
					cat_name += string("|") + names[j];
				category.push_back(cat_name);
				translation_dict[raw_val] = category.size() - 1;
				raw_val_to_new_val[field_name].push_back(raw_val);
			} while (0);
		}
		new_arr[i] = translation_dict[raw_val];
	}
	if (arr_npytype != (int)MED_NPY_TYPES::NPY_INT) {
		free(data_column[key_index]);
		data_column[key_index] = new_arr;
		data_column_nptype[key_index] = (int)MED_NPY_TYPES::NPY_INT;
	}

	categories[field_name].swap(category);
};

void MPSigExporter::get_all_data() {

	if (!this->record_count_updated)
		update_record_count();

	switch (this->sig_type) {

		//Export SDateVal

	case SigType::T_DateVal:
	{
		data_keys = vector<string>({ "pid","date","val" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		int* date_vec = (int*)malloc(sizeof(int)*this->record_count);;
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);;
		if (pid_vec == nullptr)
			throw runtime_error(string("Failed allocating memory of size ") + to_string(sizeof(int)*this->record_count) + " bytes for pid column record_count = "+to_string(this->record_count));

		if (date_vec == nullptr)
			throw runtime_error(string("Failed allocating memory of size ") + to_string(sizeof(int)*this->record_count) + " for date channel");

		if (val_vec == nullptr)
			throw runtime_error(string("Failed allocating memory of size ") + to_string(sizeof(float)*this->record_count) + " for value channel");

		int len;
		SDateVal *sdv = nullptr;
		size_t cur_row = 0;
		for (size_t j = 0; j < this->pids.size();j++) {
			int pid = this->pids[j];
			sdv = (SDateVal *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				date_vec[cur_row] = sdv[i].date;
				val_vec[cur_row] = sdv[i].val;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(date_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		gen_cat_dict("val", 0);
	}
	break;

	//export SVal

	case SigType::T_Value:
	{
		data_keys = vector<string>({ "pid","val" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);;

		int len;
		SVal *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SVal *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				val_vec[cur_row] = sdv[i].val;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		gen_cat_dict("val", 0);
	}
	break;

	//Export STimeVal

	case SigType::T_TimeVal:
	{
		data_keys = vector<string>({ "pid","time","val" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		long long* time_vec = (long long*)malloc(sizeof(long long)*this->record_count);;
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);;

		int len;
		STimeVal *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (STimeVal *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				time_vec[cur_row] = sdv[i].time;
				val_vec[cur_row] = sdv[i].val;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(time_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_LONGLONG);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		gen_cat_dict("val", 0);
	}
	break;

	//Export SDateRangeVal

	case SigType::T_DateRangeVal:
	{
		data_keys = vector<string>({ "pid","date_start","date_end","val" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		int* date_start_vec = (int*)malloc(sizeof(int)*this->record_count);;
		int* date_end_vec = (int*)malloc(sizeof(int)*this->record_count);;
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);;

		int len;
		SDateRangeVal *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SDateRangeVal *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				date_start_vec[cur_row] = sdv[i].date_start;
				date_end_vec[cur_row] = sdv[i].date_end;
				val_vec[cur_row] = sdv[i].val;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(date_start_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(date_end_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		gen_cat_dict("val", 0);
	}
	break;

	// Export STimeRangeVal

	case SigType::T_TimeRangeVal:
	{
		data_keys = vector<string>({ "pid","time_start","time_end","val" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		long long* time_start_vec = (long long*)malloc(sizeof(long long)*this->record_count);;
		long long* time_end_vec = (long long*)malloc(sizeof(long long)*this->record_count);;
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);;

		int len;
		STimeRangeVal *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (STimeRangeVal *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				time_start_vec[cur_row] = sdv[i].time_start;
				time_end_vec[cur_row] = sdv[i].time_end;
				val_vec[cur_row] = sdv[i].val;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(time_start_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_LONGLONG);
		data_column.push_back(time_end_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_LONGLONG);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		gen_cat_dict("val", 0);
	}
	break;

	// Export STimeStamp

	case SigType::T_TimeStamp:
	{
		data_keys = vector<string>({ "pid","time" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		long long* time_vec = (long long*)malloc(sizeof(long long)*this->record_count);;

		int len;
		STimeStamp *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (STimeStamp *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				time_vec[cur_row] = sdv[i].time;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(time_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_LONGLONG);
	}
	break;

	//Export SDateVal2

	case SigType::T_DateVal2:
	{
		data_keys = vector<string>({ "pid","date","val","val2" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		int* date_vec = (int*)malloc(sizeof(int)*this->record_count);;
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);
		unsigned short* val2_vec = (unsigned short*)malloc(sizeof(unsigned short)*this->record_count);;

		int len;
		SDateVal2 *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SDateVal2 *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				date_vec[cur_row] = sdv[i].date;
				val_vec[cur_row] = sdv[i].val;
				val2_vec[cur_row] = sdv[i].val2;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(date_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		data_column.push_back(val2_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_USHORT);
		gen_cat_dict("val", 0);
		gen_cat_dict("val2", 1);
	}
	break;

	//Export STimeLongVal

	case SigType::T_TimeLongVal:
	{
		data_keys = vector<string>({ "pid","time","val" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		long long* time_vec = (long long*)malloc(sizeof(long long)*this->record_count);;
		long long* val_vec = (long long*)malloc(sizeof(long long)*this->record_count);;

		int len;
		STimeLongVal *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (STimeLongVal *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				time_vec[cur_row] = sdv[i].time;
				val_vec[cur_row] = sdv[i].val;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(time_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_LONGLONG);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_LONGLONG);
		gen_cat_dict("val", 0);
	}
	break;

	//Export SDateShort2

	case SigType::T_DateShort2:
	{
		data_keys = vector<string>({ "pid","date","val1","val2" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		int* date_vec = (int*)malloc(sizeof(int)*this->record_count);;
		short* val1_vec = (short*)malloc(sizeof(short)*this->record_count);
		short* val2_vec = (short*)malloc(sizeof(short)*this->record_count);;

		int len;
		SDateShort2 *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SDateShort2 *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				date_vec[cur_row] = sdv[i].date;
				val1_vec[cur_row] = sdv[i].val1;
				val2_vec[cur_row] = sdv[i].val2;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(date_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val1_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_SHORT);
		data_column.push_back(val2_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_SHORT);
		gen_cat_dict("val1", 0);
		gen_cat_dict("val2", 1);
	}
	break;

	//Export SValShort2

	case SigType::T_ValShort2:
	{
		data_keys = vector<string>({ "pid","val1","val2" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		short* val1_vec = (short*)malloc(sizeof(short)*this->record_count);
		short* val2_vec = (short*)malloc(sizeof(short)*this->record_count);;

		int len;
		SValShort2 *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SValShort2 *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				val1_vec[cur_row] = sdv[i].val1;
				val2_vec[cur_row] = sdv[i].val2;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val1_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_SHORT);
		data_column.push_back(val2_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_SHORT);
		gen_cat_dict("val1", 0);
		gen_cat_dict("val2", 1);
	}
	break;

	//Export SValShort4

	case SigType::T_ValShort4:
	{
		data_keys = vector<string>({ "pid","val1","val2","val3","val4" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		short* val1_vec = (short*)malloc(sizeof(short)*this->record_count);
		short* val2_vec = (short*)malloc(sizeof(short)*this->record_count);;
		short* val3_vec = (short*)malloc(sizeof(short)*this->record_count);;
		short* val4_vec = (short*)malloc(sizeof(short)*this->record_count);;

		int len;
		SValShort4 *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SValShort4 *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				val1_vec[cur_row] = sdv[i].val1;
				val2_vec[cur_row] = sdv[i].val2;
				val3_vec[cur_row] = sdv[i].val3;
				val4_vec[cur_row] = sdv[i].val4;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val1_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_SHORT);
		data_column.push_back(val2_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_SHORT);
		data_column.push_back(val3_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_SHORT);
		data_column.push_back(val4_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_SHORT);
		gen_cat_dict("val1", 0);
		gen_cat_dict("val2", 1);
		gen_cat_dict("val3", 2);
		gen_cat_dict("val4", 3);
	}
	break;

	//Export SCompactDateVal

	case SigType::T_CompactDateVal:
	{
		data_keys = vector<string>({ "pid","compact_date","val" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		unsigned short* compact_date_vec = (unsigned short*)malloc(sizeof(unsigned short)*this->record_count);;
		unsigned short* val_vec = (unsigned short*)malloc(sizeof(unsigned short)*this->record_count);;

		int len;
		SCompactDateVal *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SCompactDateVal *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				compact_date_vec[cur_row] = sdv[i].compact_date;
				val_vec[cur_row] = sdv[i].val;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(compact_date_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_USHORT);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_USHORT);
		gen_cat_dict("val", 0);
	}
	break;

	//Export SDateRangeVal

	case SigType::T_DateRangeVal2:
	{
		data_keys = vector<string>({ "pid","date_start","date_end","val","val2" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		int* date_start_vec = (int*)malloc(sizeof(int)*this->record_count);
		int* date_end_vec = (int*)malloc(sizeof(int)*this->record_count);
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);
		float* val2_vec = (float*)malloc(sizeof(float)*this->record_count);

		int len;
		SDateRangeVal2 *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SDateRangeVal2 *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				date_start_vec[cur_row] = sdv[i].date_start;
				date_end_vec[cur_row] = sdv[i].date_end;
				val_vec[cur_row] = sdv[i].val;
				val2_vec[cur_row] = sdv[i].val2;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(date_start_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(date_end_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		data_column.push_back(val2_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		gen_cat_dict("val", 0);
		gen_cat_dict("val2", 1);
	}
	break;

	//Export SDateFloat2

	case SigType::T_DateFloat2:
	{
		data_keys = vector<string>({ "pid","date","val","val2" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		int* date_vec = (int*)malloc(sizeof(int)*this->record_count);
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);
		float* val2_vec = (float*)malloc(sizeof(float)*this->record_count);

		int len;
		SDateFloat2 *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (SDateFloat2 *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				date_vec[cur_row] = sdv[i].date;
				val_vec[cur_row] = sdv[i].val;
				val2_vec[cur_row] = sdv[i].val2;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(date_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(val_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		data_column.push_back(val2_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
		gen_cat_dict("val", 0);
		gen_cat_dict("val2", 1);
	}
	break;
	
	case SigType::T_TimeShort4:
	{
		data_keys = vector<string>({ "pid","time","val1","val2","val3","val4" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		long long* time_vec = (long long*)malloc(sizeof(long long)*this->record_count);
		unsigned short* val1_vec = (unsigned short*)malloc(sizeof(unsigned short)*this->record_count);
		unsigned short* val2_vec = (unsigned short*)malloc(sizeof(unsigned short)*this->record_count);
		unsigned short* val3_vec = (unsigned short*)malloc(sizeof(unsigned short)*this->record_count);
		unsigned short* val4_vec = (unsigned short*)malloc(sizeof(unsigned short)*this->record_count);

		int len;
		STimeShort4 *sdv = nullptr;
		size_t cur_row = 0;
		for (int pid : this->pids) {
			sdv = (STimeShort4 *)o->get(pid, this->sig_id, len);
			if (len == 0)
				continue;
			for (int i = 0; i < len; i++) {
				pid_vec[cur_row] = pid;
				time_vec[cur_row] = sdv[i].time;
				val1_vec[cur_row] = sdv[i].val1;
				val2_vec[cur_row] = sdv[i].val2;
				val3_vec[cur_row] = sdv[i].val3;
				val4_vec[cur_row] = sdv[i].val4;
				cur_row++;
			}
		}
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		data_column.push_back(time_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_LONGLONG);
		data_column.push_back(val1_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_USHORT);
		data_column.push_back(val2_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_USHORT);
		data_column.push_back(val3_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_USHORT);
		data_column.push_back(val4_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_USHORT);
		gen_cat_dict("val1", 0);
		gen_cat_dict("val2", 1);
		gen_cat_dict("val3", 2);
		gen_cat_dict("val4", 3);
	}
	break;

	case SigType::T_Generic:
	{
		typedef struct {
			string data_key;
			char* buf = nullptr; 
			int buf_bytes_len = 0; 
			int gsv_type = 0;
			int gsv_type_offset = 0;
			int npy_type = 0; 
			int gsv_type_bytes_len = 0; 
			int gsv_chan_num = -1;
			int rec_count = 0;
			bool is_timechan = false;
		} chan_info;
		data_keys = vector<string>();
		vector<chan_info> data_vec;
		UniversalSigVec sv;
		
		//o->usv_init(, sv);
		sv.init_from_repo(*o, this->sig_id);

		for (int tchan = 0; tchan < sv.n_time; tchan++) {
			chan_info ci;
			ci.data_key = string("time") + to_string(tchan);
			ci.rec_count = this->record_count;
			ci.gsv_type = sv.time_channel_types[tchan];
			ci.gsv_type_offset = sv.time_channel_offsets[tchan];
			ci.gsv_type_bytes_len = GenericSigVec::type_enc::bytes_len(ci.gsv_type);
			ci.buf_bytes_len = ci.rec_count * ci.gsv_type_bytes_len;
			ci.npy_type = convert_sv_type_to_npy_val_type(ci.gsv_type);
			ci.buf = (char*)malloc(ci.buf_bytes_len);
			if(ci.buf == nullptr)
				throw runtime_error(string("Failed allocating memory of size ")+to_string(ci.buf_bytes_len)+" for time channel "+to_string(tchan));
			ci.gsv_chan_num = tchan;
			ci.is_timechan = true;
			data_vec.push_back(ci);
		}

		for (int vchan = 0; vchan < sv.n_val; vchan++) {			
			chan_info ci;
			ci.data_key = string("val") + to_string(vchan);
			ci.rec_count = this->record_count;
			ci.gsv_type = sv.val_channel_types[vchan];
			ci.gsv_type_offset = sv.val_channel_offsets[vchan];
			ci.gsv_type_bytes_len = GenericSigVec::type_enc::bytes_len(ci.gsv_type);
			ci.buf_bytes_len = ci.rec_count * ci.gsv_type_bytes_len;
			ci.npy_type = convert_sv_type_to_npy_val_type(ci.gsv_type);
			ci.buf = (char*)malloc(ci.buf_bytes_len);
			if (ci.buf == nullptr)
				throw runtime_error(string("Failed allocating memory of size ") + to_string(ci.buf_bytes_len) + " for value channel " + to_string(vchan));
			ci.gsv_chan_num = vchan;
			ci.is_timechan = false;
			data_vec.push_back(ci);
		}
		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		if (pid_vec == nullptr)
			throw runtime_error(string("Failed allocating memory of size ") + to_string(sizeof(int)*this->record_count) + " for pid column");

		
		size_t cur_row = 0;
		for (int pid : this->pids) {
			o->uget(pid, this->sig_id, sv);
			const char* data = (const char*)sv.data;
			if (sv.len == 0)
				continue;
			for (int i = 0; i < sv.len; i++) {
				pid_vec[cur_row] = pid;
				for(auto& ci : data_vec){
					memcpy(&ci.buf[cur_row*ci.gsv_type_bytes_len], &data[ci.gsv_type_offset], ci.gsv_type_bytes_len);
				}
				data += sv.struct_size;
				cur_row++;
			}
		}
		data_keys.push_back("pid");
		data_column.push_back(pid_vec);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
		
		for (auto& ci : data_vec) {
			data_keys.push_back(ci.data_key);
			data_column.push_back(ci.buf);
			data_column_nptype.push_back(ci.npy_type);
		}

		for (auto& ci : data_vec) {
			if (!ci.is_timechan)
				gen_cat_dict(ci.data_key, ci.gsv_chan_num);
		}
	}
	break;


	default:
		throw runtime_error(string("MedPy: sig type not supported: ")+to_string(this->sig_type));
		break;
	}
}

void MPSigExporter::update_record_count() {
	int rec_len;
	size_t total_len = 0;
	if (this->sig_id == -1 || this->sig_type == -1) {
		this->record_count = 0;
		return;
	}
	for (int pid : this->pids)
	{
		o->get(pid, this->sig_id, rec_len);
		total_len += rec_len;
	}
	record_count_updated = true;
	this->record_count = total_len;
}

void MPSigExporter::transfer_column(const std::string& key,
	MEDPY_NP_VARIANT_OUTPUT(void** outarr1, unsigned long long* outarr1_sz, int* outarr1_npytype))
{
	int key_index = __get_key_id_or_throw(key);
	*outarr1 = data_column[key_index];
	*outarr1_sz = this->record_count;
	*outarr1_npytype = data_column_nptype[key_index];
	data_column[key_index] = nullptr;
	data_column_nptype[key_index] = (int)MED_NPY_TYPES::NPY_NOTYPE;
	data_keys[key_index] = "";

}

MPSigExporter_iter MPSigExporter::__iter__() { return MPSigExporter_iter(*this, this->data_keys); };
