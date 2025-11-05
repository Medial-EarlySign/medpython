#ifndef __MED_PY_COMMON_H
#define __MED_PY_COMMON_H

//Macros for numPy parameters - these define nothing

// Input arrays are defined as arrays of data that are passed into a routine but are not altered in - 
// place or returned to the user.The Python input array is therefore allowed to be almost any Python
// sequence(such as a list) that can be converted to the requested type of array.The input array signatures are
// ( DATA_TYPE* IN_ARRAY1, int DIM1 )
// ( int DIM1, DATA_TYPE* IN_ARRAY1 )
#define MEDPY_NP_INPUT(x1,x2) x1,x2

// In - place arrays are defined as arrays that are modified in - place.The input values may or may not be
// used, but the values at the time the function returns are significant.The provided Python argument must
// therefore be a NumPy array of the required type.The in - place signatures are
// (DATA_TYPE* INPLACE_ARRAY1, int DIM1)
// (int DIM1, DATA_TYPE* INPLACE_ARRAY1)
#define MEDPY_NP_INPLACE(x1,x2) x1,x2

// Memory Managed Argout View Arrays
// Argout arrays are arrays that appear in the input arguments in C, but are in fact output arrays.
// This pattern occurs often when there is more than one output variable and the single return argument 
// is therefore not sufficient. In Python, the conventional way to return multiple arguments is to pack 
// them into a sequence (tuple, list, etc.) and return the sequence. This is what the argout typemaps do.
// If a wrapped function that uses these argout typemaps has more than one return argument, they are 
// packed into a tuple or list, depending on the version of Python. The Python user does not pass these 
// arrays in, they simply get returned. Additinal parameter in the c++ code sets length of the allocated
// memory. The array should be allocated using malloc() call. Puthon will issue the free() call when it's 
// time to dispose the associated python object using the GC or python's del call.
// The argout signatures are
// (DATA_TYPE** ARGOUTVIEWM_ARRAY1, DIM_TYPE* DIM1)
// (DIM_TYPE* DIM1, DATA_TYPE** ARGOUTVIEWM_ARRAY1)
#define MEDPY_NP_OUTPUT(x1,x2) x1,x2

// Same as MEDPY_NP_OUTPUT but with variant array type - array datatype is set by 3rd argument which
// is a pointer to int, that int should point to a value from enum MED_NPY_TYPES defined below.
// Common types are: 
//   MED_NPY_TYPES::NPY_CHAR
//   MED_NPY_TYPES::NPY_FLOAT
//   MED_NPY_TYPES::NPY_DOUBLE
//   MED_NPY_TYPES::NPY_INT
//   MED_NPY_TYPES::NPY_LONGLONG
//   MED_NPY_TYPES::NPY_DATETIME
//   ... etc. look at MED_NPY_TYPES for more options
//
#define MEDPY_NP_VARIANT_OUTPUT(x1,x2,x3) x1,x2,x3


// Macro to attach docstrings
#ifdef SWIG
#define MEDPY_DOC(function_or_class_name, docstring) %feature("autodoc", docstring) function_or_class_name
#else 
#define MEDPY_DOC(function_or_class_name, docstring)
#endif

#define MEDPY_DOC_Dyn(text)

#define GET_MACRO(_0,_1,_2,_3,_4,_5,_6,NAME,...) NAME

#ifndef SWIG
//#define MEDPY_IGNORE(subj1) subj1

#define __MEDPY_IGNORE6(_1,_2,_3,_4,_5,_6) _1,_2,_3,_4,_5,_6
#define __MEDPY_IGNORE5(_1,_2,_3,_4,_5) _1,_2,_3,_4,_5
#define __MEDPY_IGNORE4(_1,_2,_3,_4) _1,_2,_3,_4
#define __MEDPY_IGNORE3(_1,_2,_3) _1,_2,_3
#define __MEDPY_IGNORE2(_1,_2) _1,_2
#define __MEDPY_IGNORE1(_1) _1
#define __MEDPY_IGNORE0()

#define MEDPY_IGNORE(...) GET_MACRO(_0, ##__VA_ARGS__,__MEDPY_IGNORE6, __MEDPY_IGNORE5,__MEDPY_IGNORE4, __MEDPY_IGNORE3,__MEDPY_IGNORE2,__MEDPY_IGNORE1,__MEDPY_IGNORE0)(__VA_ARGS__)
#else
#define MEDPY_IGNORE(...)
#define COMMA
#endif

#ifndef SWIG

#include <string>
#include <cstring>		// for std::memcpy()
#include <vector>
#include <stdexcept>
#include <boost/format.hpp>
#include <unordered_set>
#include <map>

using std::vector;
using std::string;
using std::runtime_error;

template<typename X>
class TypeWrapper {
public:
	X* obj = nullptr;
	typedef X objType;
	objType& instance() { return *obj; }
	virtual ~TypeWrapper() { if (obj) delete obj; obj = nullptr; }
};
/*
 * O - itearatable object type
 * I - iterator type
 * E - element type
 */
template<typename O,typename E,typename I>
class IteratorWrapper {
public:
	O* obj;
	I iterator;
	I end_iter;
	IteratorWrapper(O& o, I begin_iter_prm, I end_iter_prm) : obj(&o), iterator(begin_iter_prm), end_iter(end_iter_prm) {}
	IteratorWrapper(const IteratorWrapper& orig) : obj(orig.obj), iterator(orig.iterator), end_iter(orig.end_iter) {}

	virtual E next() = 0;

	virtual ~IteratorWrapper() { obj = nullptr; }
};

class MED_NPY_TYPE {
public:
	static const std::map<std::string,std::string> ctypestr_to_dtypestr;
	static const std::map<std::string, int> ctypestr_to_npytypeid;
	enum class values: int {
		NPY_BOOL = 0,
		NPY_BYTE, NPY_UBYTE,
		NPY_SHORT, NPY_USHORT,
		NPY_INT, NPY_UINT,
		NPY_LONG, NPY_ULONG,
		NPY_LONGLONG, NPY_ULONGLONG,
		NPY_FLOAT, NPY_DOUBLE, NPY_LONGDOUBLE,
		NPY_CFLOAT, NPY_CDOUBLE, NPY_CLONGDOUBLE,
		NPY_OBJECT = 17,
		NPY_STRING, NPY_UNICODE,
		NPY_VOID,
		NPY_DATETIME, NPY_TIMEDELTA, NPY_HALF,
		NPY_NTYPES,
		NPY_NOTYPE,
		NPY_CHAR,
		NPY_USERDEF = 256,
		NPY_NTYPES_ABI_COMPATIBLE = 21
	};
	MED_NPY_TYPE(int val) : value(val) {}
	MED_NPY_TYPE(MED_NPY_TYPE::values val) : value((int)val) {}
	int value = (int)values::NPY_NOTYPE;
	MED_NPY_TYPE& operator=(MED_NPY_TYPE o) { value = o.value; return *this; }
	MED_NPY_TYPE& operator=(int o) { value = o; return *this; }
	bool operator==(MED_NPY_TYPE o) const { return value == o.value; }
	bool operator!=(MED_NPY_TYPE o) const { return value != o.value; }
	bool operator==(int o) const { return value == o; }
	bool operator!=(int o) const { return value != o; }
	operator int() const { return value; }
	static const int sizes[];
	int size_in_bytes() {
		if (value > (int)values::NPY_CLONGDOUBLE) return -1;
		return sizes[value];
	}
};
typedef MED_NPY_TYPE::values MED_NPY_TYPES;


// "convert" a vector data buffer to a simple buffer allocated with malloc and a buffer size
template<typename S, typename T>
inline void vector_to_buf(vector<S> vec, T** buf, unsigned long long* buf_size) {
	*buf = (T*)malloc(sizeof(T)*vec.size());
	*buf_size = vec.size();
	std::memcpy(*buf, vec.data(), sizeof(T)*vec.size());
	//
	// TODO: 
	// Use this somethinglike  https://stackoverflow.com/questions/3667580/can-i-detach-a-stdvectorchar-from-the-data-it-contains
	// Detach the original buffer using Noalloc allocator and steal the original buffer with  no buffer copy operations.
	//
}

//specialization for bool/char conversion
template<>
inline void vector_to_buf<bool, char>(vector<bool> vec, char** buf, unsigned long long* buf_size) {
	*buf = (char*)malloc(sizeof(char)*vec.size());
	*buf_size = vec.size();
	for (int i = 0; i < vec.size(); ++i) {
		(*buf)[i] = (vec[i] ? 1 : 0);
	}
}

// "convert" a vector data buffer to a simple buffer allocated with malloc and a buffer size
template<typename S, typename T>
inline void buf_to_vector(S* buf, unsigned long long buf_size, vector<T>& vec) {
	vec.resize(buf_size);
	std::memcpy(vec.data(), buf, sizeof(S)*vec.size());
}

//specialization for bool/char conversion
template<>
inline void buf_to_vector<char, bool>(char* buf, unsigned long long buf_size, vector<bool>& vec) {
	vec.resize(buf_size);
	for (unsigned long long i = 0; i < buf_size; ++i) {
		vec[i] = (buf[i] != 0);
	}
}

template <template<class, class, class...> class C, typename K, typename V, typename... Args>
V GetOrDefault(const C<K, V, Args...>& m, K const& key, const V & defval)
{
	typename C<K, V, Args...>::const_iterator it = m.find(key);
	if (it == m.end())
		return defval;
	return it->second;
}

#pragma warning( disable : 4290 )

int val_or_exception(int val, const string& exception_str = "", int expected_val = 0);

#endif // !SWIG

class StopIterator {};

class MPIntIntMapAdaptor {
	bool o_owned = true;
	std::map<int, int>* o;
public:
	MPIntIntMapAdaptor();
	MEDPY_IGNORE(MPIntIntMapAdaptor(const MPIntIntMapAdaptor& other));
	MEDPY_IGNORE(MPIntIntMapAdaptor(std::map<int, int>* ptr));
	~MPIntIntMapAdaptor();
	int __len__();
	int __getitem__(int i);
	void __setitem__(int i, int val);
	void keys(MEDPY_NP_OUTPUT(int** intkeys_out_buf, unsigned long long* intkeys_out_buf_len));
	MEDPY_IGNORE(MPIntIntMapAdaptor& operator=(const MPIntIntMapAdaptor& other));
};

class MPStringStringMapAdaptor {
	bool o_owned = true;
	std::map<std::string, std::string>* o;
public:
	MPStringStringMapAdaptor();
	MEDPY_IGNORE(MPStringStringMapAdaptor(const MPStringStringMapAdaptor& other));
	MEDPY_IGNORE(MPStringStringMapAdaptor(std::map<std::string, std::string>* ptr));
	~MPStringStringMapAdaptor();
	int __len__();
	std::string __getitem__(const std::string& i);
	void __setitem__(const std::string& i, const std::string& val);
	std::vector<std::string> keys();
	MEDPY_IGNORE(MPStringStringMapAdaptor& operator=(const MPStringStringMapAdaptor& other));
};

class MPStringFloatMapAdaptor {
	bool o_owned = true;
public:
	MEDPY_IGNORE(std::map<std::string, float>* o);
	MPStringFloatMapAdaptor();
	MEDPY_IGNORE(MPStringFloatMapAdaptor(const MPStringFloatMapAdaptor& other));
	MEDPY_IGNORE(MPStringFloatMapAdaptor(std::map<std::string, float>* ptr));
	~MPStringFloatMapAdaptor();
	int __len__();
	float __getitem__(const std::string& i);
	void __setitem__(const std::string& i, float val);
	std::vector<std::string> keys();
	MEDPY_IGNORE(MPStringFloatMapAdaptor& operator=(const MPStringFloatMapAdaptor& other));
};


class MPStringVecFloatMapAdaptor {
	bool o_owned = true;
	std::map<std::string, std::vector<float> >* o;
public:
	MPStringVecFloatMapAdaptor();
	MEDPY_IGNORE(MPStringVecFloatMapAdaptor(const MPStringVecFloatMapAdaptor& other));
	MEDPY_IGNORE(MPStringVecFloatMapAdaptor(std::map<std::string, std::vector<float> >* ptr));
	~MPStringVecFloatMapAdaptor();
	int __len__();
	void __getitem__(std::string key, MEDPY_NP_OUTPUT(float** float_out_buf, unsigned long long* float_out_buf_len));
	void __setitem__(std::string key, MEDPY_NP_INPUT(float* float_in_buf, unsigned long long float_in_buf_len));
	std::vector<std::string> keys();
	MEDPY_IGNORE(MPStringVecFloatMapAdaptor& operator=(const MPStringVecFloatMapAdaptor& other));
};

class MPIntPairIntIntMapAdaptor {
	bool o_owned = true;
	std::map<int, std::pair<int, int> >* o;
public:
	MPIntPairIntIntMapAdaptor();
	MEDPY_IGNORE(MPIntPairIntIntMapAdaptor(const MPIntPairIntIntMapAdaptor& other));
	MEDPY_IGNORE(MPIntPairIntIntMapAdaptor(std::map<int, std::pair<int, int> >* ptr));
	~MPIntPairIntIntMapAdaptor();
	int __len__();
	void __getitem__(int key, MEDPY_NP_OUTPUT(int** int_out_buf, unsigned long long* int_out_buf_len));
	void __setitem__(int key, MEDPY_NP_INPUT(int* int_in_buf, unsigned long long int_in_buf_len));
	void keys(MEDPY_NP_OUTPUT(int** intkeys_out_buf, unsigned long long* intkeys_out_buf_len));
	MEDPY_IGNORE(MPIntPairIntIntMapAdaptor& operator=(const MPIntPairIntIntMapAdaptor& other));
};


class MPStringUOSetStringMapAdaptor {
	bool o_owned = true;
	std::map<std::string, std::unordered_set<std::string> >* o;
public:
	MPStringUOSetStringMapAdaptor();
	MEDPY_IGNORE(MPStringUOSetStringMapAdaptor(const MPStringUOSetStringMapAdaptor& other));
	MEDPY_IGNORE(MPStringUOSetStringMapAdaptor(std::map<std::string, std::unordered_set<std::string> >* ptr));
	~MPStringUOSetStringMapAdaptor();
	int __len__();
	std::vector<std::string> __getitem__(std::string key);
	void __setitem__(std::string key, std::vector<std::string> val);
	std::vector<std::string> keys();
	MEDPY_IGNORE(MPStringUOSetStringMapAdaptor& operator=(const MPStringUOSetStringMapAdaptor& other));
};

class MPIntStringMapAdaptor {
	bool o_owned = true;
	std::map<int, string>* o;
public:
	MPIntStringMapAdaptor();
	MEDPY_IGNORE(MPIntStringMapAdaptor(const MPIntStringMapAdaptor& other));
	MEDPY_IGNORE(MPIntStringMapAdaptor(std::map<int, string>* ptr));
	~MPIntStringMapAdaptor();
	int __len__();
	std::string __getitem__(int i);
	void __setitem__(int i, const string& val);
	std::vector<int> keys();
	MEDPY_IGNORE(MPIntStringMapAdaptor& operator=(const MPIntStringMapAdaptor& other));
};


class MPIntVecIntMapAdaptor {
	bool o_owned = true;
public:
	std::map<int, std::vector<int> >* o;
	MPIntVecIntMapAdaptor();
	MEDPY_IGNORE(MPIntVecIntMapAdaptor(const MPIntVecIntMapAdaptor& other));
	MEDPY_IGNORE(MPIntVecIntMapAdaptor(std::map<int, std::vector<int> >* ptr));
	~MPIntVecIntMapAdaptor();
	int __len__();
	void __getitem__(int key, MEDPY_NP_OUTPUT(int** int_out_buf, unsigned long long* int_out_buf_len));
	void __setitem__(int key, MEDPY_NP_INPUT(int* int_in_buf, unsigned long long int_in_buf_len));
	std::vector<int> keys();
	MEDPY_IGNORE(MPIntVecIntMapAdaptor& operator=(const MPIntVecIntMapAdaptor& other));
};

void logger_use_stdout();

#endif // !__MED_PY_COMMON_H
