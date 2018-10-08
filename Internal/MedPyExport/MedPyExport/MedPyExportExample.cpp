#include "MedPyExportExample.h"

#include <time.h>
#include <string>


void MedPyExportExample::numpy_vec_in_out(MEDPY_NP_INPLACE(double* vec, int m)) {
	int  i;
	for (i = 0; i<m; i++) {
		vec[i] = 2 * vec[i];
	}
}

void MedPyExportExample::numpy_vec_out(MEDPY_NP_OUTPUT(double** array, int* n), int ar_len)
{
	//int ar_len = 10;
	double* arr = (double*)malloc(sizeof(double)*ar_len);
	int i;

	*array = arr;
	*n = ar_len;

	for (i = 0; i<*n; i++)
		(*array)[i] = i;
}
int MedPyExportExample::numpy_vec_input(MEDPY_NP_INPUT(double* vec2, int nn)) {
	double sum = 0;
	for (int i = 0; i<nn; i++)
		sum += vec2[i];
	return (int)sum;
}

void MedPyExportExample::numpy_vec_out2(MEDPY_NP_OUTPUT(double** array11, int* n11), 
	MEDPY_NP_OUTPUT(double** array22, int* n22), int ar_len)
{
	{
		double* arr = (double*)malloc(sizeof(double)*ar_len);
		int i;

		*array11 = arr;
		*n11 = ar_len;

		for (i = 0; i < *n11; i++)
			(*array11)[i] = i;
	}
	{
		ar_len *= 2;
		double* arr = (double*)malloc(sizeof(double)*ar_len);
		int i;

		*array22 = arr;
		*n22 = ar_len;

		for (i = 0; i < *n22; i++)
			(*array22)[i] = i;
	}
}

MedPyExportExample2::MedPyExportExample2(int size) {
	obj = new objType();
	obj->resize(size, false);
}

int MedPyExportExample2::__getitem__(int i) {
	return (*obj)[i];
}

void MedPyExportExample2::__setitem__(int i, int val) {
	(*obj)[i] = (val != 0);
}

MedPyExportExample2_iterator MedPyExportExample2::__iter__() {
	return MedPyExportExample2_iterator(*this, 0, (int)obj->size());
}

/*
bool MedPyExportExample2_iterator::isEnd() {
return (iterator >= obj->instance().size());
}*/


MedPyExportExample3::MedPyExportExample3(MedPyExportExample2& o) {
	count = 0;
	for (auto b : o.instance())
		if (b)
			count++;
}


int MedPyExportExample4::MEDPY_GET_x() { return _x; }

void MedPyExportExample4::MEDPY_SET_x(int new_x) { _x = new_x; }

int MedPyExportExample4::MEDPY_GET_ro_x() { return _x; }

void MedPyExportExample4::MEDPY_SET_wo_x(int new_x) { _x = new_x; }

void MedPyExportExample5::getitem(const std::string& key, MEDPY_NP_VARIANT_OUTPUT(void** var_arr, int* var_arr_sz, int* var_arr_type)) {
	*var_arr_sz = 0;
	if (key == "i")
	{
		*var_arr = (void*)malloc(sizeof(int) * 10);
		*var_arr_sz = 10;
		*var_arr_type = (int)MED_NPY_TYPES::NPY_INT;
		for (int i = 0; i < 10; i++)
			((*(int**)var_arr))[i] = i * 5;
	}
	else if (key == "d")
	{
		*var_arr = (void*)malloc(sizeof(double) * 20);
		*var_arr_sz = 20;
		*var_arr_type = (int)MED_NPY_TYPES::NPY_DOUBLE;
		for (int i = 0; i < 20; i++)
			((*(double**)var_arr))[i] = i * 2.5;
	}
	else if (key == "f")
	{
		*var_arr = (void*)malloc(sizeof(float) * 15);
		*var_arr_sz = 15;
		*var_arr_type = (int)MED_NPY_TYPES::NPY_FLOAT;
		for (int i = 0; i < 15; i++)
			((*(float**)var_arr))[i] = i * 0.33333f;
	}
	else if (key == "n")
	{
		*var_arr = nullptr;
	}
}



int MPExample6::strvec_arg(std::vector<std::string> strarr) {
	return (int)strarr.size();
}

int MPExample6::strvec_arg2(std::vector<std::string>& strarr) {
	strarr.pop_back();
	return (int)strarr.size();
}
int MPExample6::strvec_arg3(vector<std::string> strarr) {
	return (int)strarr.size();
}
int MPExample6::strvec_arg4(std::vector<string> strarr) {
	return (int)strarr.size();
}
int MPExample6::strvec_arg5(vector<string> strarr) {
	return (int)strarr.size();
}
double MPExample6::dbl_arg(double darg) { return darg + 0.5; }


