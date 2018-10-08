#ifndef __MED_PY_EXPORT_EXAMPLE_H
#define __MED_PY_EXPORT_EXAMPLE_H

#include "MedPyCommon.h"

/*
*  Never use const ref ('const x&') passing , only &x
*
*/


/*
*
*  Example 1 : numpy array for input and output
*
*/


class MedPyExportExample {
public:
	void numpy_vec_in_out(MEDPY_NP_INPLACE(double* vec, int m));

	void numpy_vec_out(MEDPY_NP_OUTPUT(double** array, int* n), int ar_len);

	int numpy_vec_input(MEDPY_NP_INPUT(double* vec2, int nn));

	void numpy_vec_out2(MEDPY_NP_OUTPUT(double** array11, int* n11), MEDPY_NP_OUTPUT(double** array22, int* n22), int ar_len);

	int property = 555;
};

/*
*
*  Example 2 - a class with index access and iterator
*
*/

class MedPyExportExample2_iterator;

class MedPyExportExample2 MEDPY_IGNORE(:public TypeWrapper<vector<bool>>) {
public:
	MedPyExportExample2(int size);

	//Index access :
	int __getitem__(int i);
	void __setitem__(int i, int val);

	//iterator
	MedPyExportExample2_iterator __iter__();
};

/* iterator class for MedPyExportExample2 */

class MedPyExportExample2_iterator MEDPY_IGNORE(:public IteratorWrapper<MedPyExportExample2, int, int>) {
public:
	MedPyExportExample2_iterator(MedPyExportExample2& o, int begin, int end) : IteratorWrapper(o, begin, end) {}

	int derefrence() { return obj->__getitem__(iterator); }
	bool isEnd() { return iterator == end_iter; }

	int next() { if (this->isEnd()) throw StopIterator(); auto ret = this->derefrence(); ++iterator;  return ret; }

	std::string __repr__() {
		if (isEnd()) return string("[END]");
		return string("[") + std::to_string(iterator) + string("]=") + std::to_string(this->derefrence());
	}
};

/*
* Example 3 : Another object as a parameter
*/

class MedPyExportExample3 {
public:
	MedPyExportExample3(MedPyExportExample2& o);
	int count;
};

/*
* Getter/setters to become properties
*
*/


class MedPyExportExample4 {
	int _x = 500;
public:
    // property 'x'
	int MEDPY_GET_x();
	void MEDPY_SET_x(int new_x);

	// readonly property 'ro_x'
	int MEDPY_GET_ro_x();

	// writeonly property 'wo_x'
	void MEDPY_SET_wo_x(int new_x);

};

/*
 * returning variant type of array.
 * In this example the returned value can be anything we choose. the buffer is in a void* ptr and the type
 * is set by the programmer using the 3rd argument, i.e. var_arr_type below and the value from MED_NPY_TYPES,
 * defined in MedPyCommon.h
 * when key == "i" it will return numpy array of ints
 * == "d" - array of doubles
 * == "f" - array of floats
 * == "n" - return nothing - None
 */


class MedPyExportExample5 {
public:
	void getitem(const std::string& key, MEDPY_NP_VARIANT_OUTPUT(void** var_arr, int* var_arr_sz, int* var_arr_type));
};



/*
* vector of string example
*
*/

class MPExample6 {
public:
	double dbl_arg(double darg);
	int strvec_arg(std::vector<std::string> strarr);
	int strvec_arg2(std::vector<std::string>& strarr);
	int strvec_arg3(vector<std::string> strarr);
	int strvec_arg4(std::vector<string> strarr);
	int strvec_arg5(vector<string> strarr);
};



#endif // !__MED_PY_EXPORT_EXAMPLE_H

