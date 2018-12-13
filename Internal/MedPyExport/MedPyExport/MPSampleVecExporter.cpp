#include "MPSampleVecExporter.h"

#include <time.h>
#include <string>
#include <unordered_set>

#include "InfraMed/InfraMed/MedConvert.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/Utils.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedModel.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"


MPSampleVecExporter::MPSampleVecExporter(MPSampleVectorAdaptor& rep) : o(rep.o) {
	update_record_count();
	get_all_data();
}


int MPSampleVecExporter::__get_key_id_or_throw(const string& key) {
	for (int i = 0; i < data_keys.size();++i)
		if (data_keys[i] == key) return i;
	throw runtime_error("Unknown row");
}

void MPSampleVecExporter::gen_cat_dict(const string& field_name, int channel) {
	/*
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
	*/
};

void MPSampleVecExporter::get_all_data() {
	update_record_count();
	data_keys = vector<string>({ "id","split","time","outcome","outcomeTime" });
	
	int max_predictions = 0;
	vector<string> attribute_names;
	vector<string> str_attribute_names;

	{
		unordered_set<string> attribute_names_set;
		unordered_set<string> str_attribute_names_set;

		for (auto& samp : *o) {
			if (samp.prediction.size() > max_predictions)
				max_predictions = (int)samp.prediction.size();
			if (samp.attributes.size() != 0)
				for (auto& entry : samp.attributes)
					attribute_names_set.insert(entry.first);
			if (samp.str_attributes.size() != 0)
				for (auto& entry : samp.str_attributes)
					str_attribute_names_set.insert(entry.first);
		}
		copy(attribute_names_set.begin(), attribute_names_set.end(), attribute_names.end());
		copy(str_attribute_names_set.begin(), str_attribute_names_set.end(), attribute_names.end());
		//for (auto& s : attribute_names_set) attribute_names.push_back(s);
		//for (auto& s : str_attribute_names_set) str_attribute_names.push_back(s);
	}

	int* id_vec = (int*)malloc(sizeof(int)*this->record_count);
	int* split_vec = (int*)malloc(sizeof(int)*this->record_count);
	int* time_vec = (int*)malloc(sizeof(int)*this->record_count);
	float* outcome_vec = (float*)malloc(sizeof(float)*this->record_count);
	int* outcome_time_vec = (int*)malloc(sizeof(int)*this->record_count);
	
	vector<float*> pred_vecs;
	for (int i = 0; i < max_predictions; i++) {
		pred_vecs.push_back((float*)malloc(sizeof(float)*this->record_count));
		data_keys.push_back((string)"pred_"+to_string(i));
	}

	int attr_num = (int)attribute_names.size();
	vector<float*> attr_vecs;
	for (int i = 0; i < attr_num; i++) {
		attr_vecs.push_back((float*)malloc(sizeof(float)*this->record_count));
		data_keys.push_back((string)"attr_" + attribute_names[i]);
	}

	int cur_row = 0;
	for (auto& samp : *o) {
		id_vec[cur_row] = samp.id;
		split_vec[cur_row] = samp.split;
		time_vec[cur_row] = samp.time;
		outcome_vec[cur_row] = samp.outcome;
		outcome_time_vec[cur_row] = samp.outcomeTime;
		if (max_predictions > 0)
			for (int i = 0; i < max_predictions; i++)
				pred_vecs[i][cur_row] = i < samp.prediction.size() ? samp.prediction[i] : -1;
		if (attr_num > 0)
			for (int i = 0; i < attr_num; ++i)
				attr_vecs[i][cur_row] = GetOrDefault(samp.attributes, attribute_names[i], -1.0f);
				//string attr_name = attribute_names[i];
				//attr_vecs[i][cur_row] = i < samp.attributes.count(attr_name) ? samp.attributes[attr_name] : -1;
				
		cur_row++;
	}
	
	data_column.push_back(id_vec);
	data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
	data_column.push_back(split_vec);
	data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
	data_column.push_back(time_vec);
	data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);
	data_column.push_back(outcome_vec);
	data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
	data_column.push_back(outcome_time_vec);
	data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_INT);

	for (int i = 0; i < max_predictions; i++) {
		data_column.push_back(pred_vecs[i]);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
	}

	for (int i = 0; i < attr_num; i++) {
		data_column.push_back(attr_vecs[i]);
		data_column_nptype.push_back((int)MED_NPY_TYPES::NPY_FLOAT);
	}


	//gen_cat_dict("val", 0);
}

void MPSampleVecExporter::update_record_count() {
	this->record_count = (int)o->size();
}

void MPSampleVecExporter::transfer_column(const std::string& key,
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

MPSampleVecExporter_iter MPSampleVecExporter::__iter__() { return MPSampleVecExporter_iter(*this, this->data_keys); };
