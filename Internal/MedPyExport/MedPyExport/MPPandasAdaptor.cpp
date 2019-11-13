#include "MPPandasAdaptor.h"

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
void MPPandasAdaptor::set_type_requirement(const string& col_name, const string& ctype_str){
	if (MED_NPY_TYPE::ctypestr_to_dtypestr.count(ctype_str) == 0)
		throw runtime_error(string("unrecognized C type for type requirement")+ ctype_str);
	type_requirements.__setitem__(col_name, MED_NPY_TYPE::ctypestr_to_dtypestr.at(ctype_str));
}

MPPandasAdaptor::column_record& MPPandasAdaptor::__get_col_or_throw(const string& col_name) {
	if(columns.count(col_name) == 0)
		throw runtime_error("Key error - No row by that key");
	return columns.at(col_name);
}

void MPPandasAdaptor::push_categorial(const string& col_name, std::vector<int> index_column, std::vector<std::string> categories_vec) {
	if (columns.count(col_name) > 0)
		throw runtime_error("PandasAdaptor: already have a column by this name");
	auto& col = columns[col_name];
	col.nptype = MED_NPY_TYPES::NPY_INT;
	vector_to_buf(index_column, (int**)(&(col.data)), &(col.num_elements));
	col.is_own = true;
	auto& cat_col = categories[col_name];
	for (int i = 0; i < categories_vec.size(); ++i)
		cat_col[i] = categories_vec[i];
};

void MPPandasAdaptor::export_column(const std::string& key,
	MEDPY_NP_VARIANT_OUTPUT(void** outarr1, unsigned long long* outarr1_sz, int* outarr1_npytype))
{
	auto& col = __get_col_or_throw(key);
	if (col.is_own) {
		*outarr1 = col.data;
		col.is_own = false;
		col.data = nullptr;
	}
	else {
		*outarr1 = malloc(col.data_size());
		memcpy(*outarr1, col.data, col.data_size());
	}
	*outarr1_sz = col.num_elements;
	*outarr1_npytype = (int)col.nptype;
	columns.erase(key);
}

void MPPandasAdaptor::import_column(const string &col_name, void* IN_ARRAY1, unsigned long long DIM1, int NPYDTC1, bool make_a_copy) {
	if (columns.count(col_name) > 0)
		throw runtime_error("PandasAdaptor: already have a column by this name");
	auto& col = columns[col_name];
	col.nptype = NPYDTC1;	
	col.num_elements = DIM1;
	if (make_a_copy) {
		col.allocate_data();
		memcpy(col.data, IN_ARRAY1, col.data_size());
		col.is_own = true;
	}
	else {
		col.data = IN_ARRAY1;
		col.is_own = false;
	}
}


MPPandasAdaptor_iter MPPandasAdaptor::__iter__() { return MPPandasAdaptor_iter(*this, this->keys()); };
