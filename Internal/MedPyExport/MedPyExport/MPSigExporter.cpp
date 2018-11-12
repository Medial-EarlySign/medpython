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

MPSigExporter MPPidRepository::export_to_numpy(string signame) {
	return MPSigExporter(*this, signame);
}


int MPSigExporter::__get_key_id_or_throw(const string& key) {
	for (int i = 0; i < data_keys.size();++i)
		if (data_keys[i] == key) return i;
	throw runtime_error("Unknown row");
}

void MPSigExporter::gen_cat_dict(const string& field_name, int channel) {
	int key_index = __get_key_id_or_throw(field_name);
	if (!o->sigs.is_categorical_channel(sig_id, channel))
		return;
	int section_id = o->dict.section_id(sig_name);
	void* arr = data_column[key_index];
	int arr_sz = this->record_count;
	int arr_npytype = data_column_nptype[key_index];
	std::unordered_set<int> values;
	switch (arr_npytype) {
	case (int)MED_NPY_TYPES::NPY_FLOAT: 
		{float* tarr = (float*)arr; for (int i = 0; i < arr_sz; ++i) values.insert((int)tarr[i]); }
		break;
	case (int)MED_NPY_TYPES::NPY_USHORT: 
		{unsigned short* tarr = (unsigned short*)arr; for (int i = 0; i < arr_sz; ++i) values.insert((int)tarr[i]); }
		break;
	case (int)MED_NPY_TYPES::NPY_LONGLONG: 
		{long long* tarr = (long long*)arr; for (int i = 0; i < arr_sz; ++i) values.insert((int)tarr[i]); }
		break;
	case (int)MED_NPY_TYPES::NPY_SHORT: 
		{short* tarr = (short*)arr; for (int i = 0; i < arr_sz; ++i) values.insert((int)tarr[i]); }
		break;
	default:
		throw runtime_error("MedPy: categorical value type not supported, we only have values of types float, unsigned short, long long, short");
		break;
	}
	auto& Id2Names = o->dict.dict(section_id)->Id2Names;
	std::map<int, std::string> cat_dict;
	for (int raw_val : values) {
		if (!Id2Names.count(raw_val)) continue;
		auto& names = Id2Names[raw_val];
		if (names.size() == 0) { cat_dict[raw_val] = ""; continue; }
		cat_dict[raw_val] = names[0];
		for (int j = 1; j < names.size() && j < 3; j++)
			cat_dict[raw_val] += string("|") + names[j];
	}
	categories[field_name] = cat_dict;
};

void MPSigExporter::get_all_data() {

	if (this->record_count <= 0)
		update_record_count();

	switch (this->sig_type) {

		//Export SDateVal

	case SigType::T_DateVal:
	{
		data_keys = vector<string>({ "pid","date","val" });

		int* pid_vec = (int*)malloc(sizeof(int)*this->record_count);
		int* date_vec = (int*)malloc(sizeof(int)*this->record_count);;
		float* val_vec = (float*)malloc(sizeof(float)*this->record_count);;

		int len;
		SDateVal *sdv = nullptr;
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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
		int cur_row = 0;
		for (int pid : o->all_pids_list) {
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

	default:
		throw runtime_error("MedPy: sig type not supported");
		break;
	}
}

void MPSigExporter::update_record_count() {
	int rec_len;
	int total_len = 0;
	if (this->sig_id == -1 || this->sig_type == -1) {
		this->record_count = 0;
		return;
	}
	for (int pid : o->all_pids_list)
	{
		o->get(pid, this->sig_id, rec_len);
		total_len += rec_len;
	}
	this->record_count = total_len;
}

void MPSigExporter::transfer_column(const std::string& key,
	MEDPY_NP_VARIANT_OUTPUT(void** outarr1, int* outarr1_sz, int* outarr1_npytype))
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
