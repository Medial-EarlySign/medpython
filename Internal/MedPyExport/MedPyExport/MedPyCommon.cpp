#include "MedPyCommon.h"

/* https ://www.numpy.org/devdocs/user/basics.types.html
np.dtype		C++ type			dtype('..') Desc
np.bool			bool 				bool		Boolean(True or False) stored as a byte
np.byte			signed char 		int8		Platform - defined
np.ubyte		unsigned char 		uint8		Platform - defined
np.short		short 				int16		Platform - defined
np.ushort		unsigned short 		uint16		Platform - defined
np.intc			int 				int32		Platform - defined
np.uintc		unsigned int 		uint32		Platform - defined
np.int_			long 				int64		Platform - defined
np.uint			unsigned long 		uint64		Platform - defined
np.longlong		long long 			int64		Platform - defined
np.ulonglong	unsigned long long	uint64		Platform - defined
np.half/np.float16 	  				float16		Half precision float: sign bit, 5 bits exponent, 10 bits mantissa
np.single		float 				float32		Platform - defined single precision float : typically sign bit, 8 bits exponent, 23 bits mantissa
np.double		double 				float64		Platform - defined double precision float : typically sign bit, 11 bits exponent, 52 bits mantissa.
np.longdouble	long double			float128	Platform - defined extended - precision float
np.csingle		float complex 		complex64	Complex number, represented by two single - precision floats(real and imaginary components)
np.cdouble		double complex 		complex128	Complex number, represented by two double - precision floats(real and imaginary components).
np.clongdouble	long double complex complex256	Complex number, represented by two extended - precision floats(real and imaginary components).
*/
const std::map<std::string, std::string> MED_NPY_TYPE::ctypestr_to_dtypestr = {
{"bool","bool"},
{"char","int8" },
{"signed char","int8"},
{"unsigned char","uint8"},
{"short","int16"},
{"unsigned short","uint16"},
{"int","int32"},
{"unsigned int","uint32"},
{"long","int64"},
{"unsigned long","uint64"},
{"long long","int64"},
{"unsigned long long","uint64"},
{"float","float32"},
{"double","float64"},
{"long double","float128"},
{"float complex","complex64"},
{"double complex","complex128"},
{"long double complex","complex256"}
};

const std::map<std::string, int> MED_NPY_TYPE::ctypestr_to_npytypeid = {
	{ "bool",(int)MED_NPY_TYPES::NPY_BOOL },
	{ "char",(int)MED_NPY_TYPES::NPY_BYTE },
	{ "signed char",(int)MED_NPY_TYPES::NPY_BYTE },
	{ "unsigned char",(int)MED_NPY_TYPES::NPY_UBYTE },
	{ "short",(int)MED_NPY_TYPES::NPY_SHORT },
	{ "unsigned short",(int)MED_NPY_TYPES::NPY_USHORT },
	{ "int",(int)MED_NPY_TYPES::NPY_INT },
	{ "unsigned int",(int)MED_NPY_TYPES::NPY_UINT },
	{ "long",(int)MED_NPY_TYPES::NPY_LONG },
	{ "unsigned long",(int)MED_NPY_TYPES::NPY_ULONG },
	{ "long long",(int)MED_NPY_TYPES::NPY_LONGLONG },
	{ "unsigned long long",(int)MED_NPY_TYPES::NPY_ULONGLONG },
	{ "float",(int)MED_NPY_TYPES::NPY_FLOAT },
	{ "double",(int)MED_NPY_TYPES::NPY_DOUBLE },
	{ "long double",(int)MED_NPY_TYPES::NPY_CLONGDOUBLE },
	{ "float complex",(int)MED_NPY_TYPES::NPY_CFLOAT },
	{ "double complex",(int)MED_NPY_TYPES::NPY_CDOUBLE },
	{ "long double complex",(int)MED_NPY_TYPES::NPY_CLONGDOUBLE }
};

 const int MED_NPY_TYPE::sizes[] {
	sizeof(unsigned char),									//NPY_BOOL
	sizeof(char),sizeof(unsigned char),						//NPY_BYTE, NPY_UBYTE,
	sizeof(short),sizeof(unsigned short),					//NPY_SHORT, NPY_USHORT,
	sizeof(int),sizeof(unsigned int),						//NPY_INT, NPY_UINT,
	sizeof(long int),sizeof(unsigned long int),				//NPY_LONG, NPY_ULONG,
	sizeof(long long int),sizeof(unsigned long long int),	//NPY_LONGLONG, NPY_ULONGLONG,
	sizeof(float),sizeof(double),sizeof(long double),		//NPY_FLOAT, NPY_DOUBLE, NPY_LONGDOUBLE,
	sizeof(float) * 2,sizeof(double) * 2,sizeof(long double) * 2	//NPY_CFLOAT, NPY_CDOUBLE, NPY_CLONGDOUBLE,
																	//-1, //NPY_OBJECT
																	//-1,-1, //NPY_STRING, NPY_UNICODE,
																	//-1, //NPY_VOID
																	//-1,-1,-1 //NPY_DATETIME, NPY_TIMEDELTA, NPY_HALF,
																	//-1, //NPY_NTYPES
																	//-1, //NPY_NOTYPE
																	//-1, //NPY_CHAR
};

 /******************************************************************************************************************************/

MPIntIntMapAdaptor::MPIntIntMapAdaptor() { o = new std::map<int, int>(); };
MPIntIntMapAdaptor::MPIntIntMapAdaptor(const MPIntIntMapAdaptor& other)
{
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<int, int>();
		*o = *other.o;
	}
};

MPIntIntMapAdaptor::MPIntIntMapAdaptor(std::map<int, int>* ptr) { o_owned = false; o = ptr; };
MPIntIntMapAdaptor::~MPIntIntMapAdaptor() { if (o_owned) delete o; };
int MPIntIntMapAdaptor::__len__() { return (int)o->size(); };
int MPIntIntMapAdaptor::__getitem__(int i) { return o->operator[](i); };
void MPIntIntMapAdaptor::__setitem__(int i, int val) { o->operator[](i) = val; };
void MPIntIntMapAdaptor::keys(MEDPY_NP_OUTPUT(int** intkeys_out_buf, unsigned long long* intkeys_out_buf_len)) 
{
	vector<int> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	vector_to_buf(ret, intkeys_out_buf, intkeys_out_buf_len);
};

MPIntIntMapAdaptor& MPIntIntMapAdaptor::operator=(const MPIntIntMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<int, int>();
		*o = *(other.o);
	}
	return *this;
}

/******************************************************************************************************************************/
MPStringStringMapAdaptor::MPStringStringMapAdaptor() { o = new std::map<std::string, std::string>(); };
MPStringStringMapAdaptor::MPStringStringMapAdaptor(const MPStringStringMapAdaptor& other)
{
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, std::string>();
		*o = *other.o;
	}
};

MPStringStringMapAdaptor::MPStringStringMapAdaptor(std::map<std::string, std::string>* ptr) { o_owned = false; o = ptr; };
MPStringStringMapAdaptor::~MPStringStringMapAdaptor() { if (o_owned) delete o; };
int MPStringStringMapAdaptor::__len__() { return (int)o->size(); };
std::string MPStringStringMapAdaptor::__getitem__(const std::string& i) { return o->operator[](i); };
void MPStringStringMapAdaptor::__setitem__(const std::string& i, const std::string& val) { o->operator[](i) = val; };
std::vector<std::string> MPStringStringMapAdaptor::keys()
{
	vector<string> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};

MPStringStringMapAdaptor& MPStringStringMapAdaptor::operator=(const MPStringStringMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, std::string>();
		*o = *(other.o);
	}
	return *this;
}
/******************************************************************************************************************************/
MPStringFloatMapAdaptor::MPStringFloatMapAdaptor() { o = new std::map<std::string, float>(); }
MPStringFloatMapAdaptor::MPStringFloatMapAdaptor(const MPStringFloatMapAdaptor& other)
{
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, float>();
		*o = *other.o;
	}
};
MPStringFloatMapAdaptor::MPStringFloatMapAdaptor(std::map<std::string, float>* ptr) { o_owned = false; o = ptr; };
MPStringFloatMapAdaptor::~MPStringFloatMapAdaptor() { if (o_owned) delete o; };
int MPStringFloatMapAdaptor::__len__() { return (int)o->size(); };
float MPStringFloatMapAdaptor::__getitem__(const std::string& i) { return o->operator[](i); };
void MPStringFloatMapAdaptor::__setitem__(const std::string& i, float val) { o->operator[](i) = val; };
std::vector<std::string> MPStringFloatMapAdaptor::keys()
{
	std::vector<std::string> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};

MPStringFloatMapAdaptor& MPStringFloatMapAdaptor::operator=(const MPStringFloatMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, float>();
		*o = *(other.o);
	}
	return *this;
}

/******************************************************************************************************************************/
MPStringVecFloatMapAdaptor::MPStringVecFloatMapAdaptor() { o = new std::map<std::string, std::vector<float> >(); };
MPStringVecFloatMapAdaptor::MPStringVecFloatMapAdaptor(const MPStringVecFloatMapAdaptor& other) {
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, std::vector<float> >();
		*o = *other.o;
	}
};
MPStringVecFloatMapAdaptor::MPStringVecFloatMapAdaptor(std::map<std::string, std::vector<float> >* ptr) { o_owned = false; o = ptr; };
MPStringVecFloatMapAdaptor::~MPStringVecFloatMapAdaptor() { if (o_owned) delete o; };
int MPStringVecFloatMapAdaptor::__len__() { return (int)o->size(); };
void MPStringVecFloatMapAdaptor::__getitem__(std::string key, MEDPY_NP_OUTPUT(float** float_out_buf, unsigned long long* float_out_buf_len)) {
	auto& fvec = o->operator[](key);
	vector_to_buf(fvec, float_out_buf, float_out_buf_len);
};
void MPStringVecFloatMapAdaptor::__setitem__(std::string key, MEDPY_NP_INPUT(float* float_in_buf, unsigned long long float_in_buf_len)) {
	vector<float> fvec;
	buf_to_vector(float_in_buf, float_in_buf_len, fvec);
	o->operator[](key) = fvec;
};
std::vector<std::string> MPStringVecFloatMapAdaptor::keys() {
	vector<string> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};

MPStringVecFloatMapAdaptor& MPStringVecFloatMapAdaptor::operator=(const MPStringVecFloatMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, std::vector<float> >();
		*o = *(other.o);
	}
	return *this;
}




MPIntPairIntIntMapAdaptor::MPIntPairIntIntMapAdaptor() { o = new std::map<int, std::pair<int, int> >(); };
MPIntPairIntIntMapAdaptor::MPIntPairIntIntMapAdaptor(const MPIntPairIntIntMapAdaptor& other) {
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<int, std::pair<int, int> >();
		*o = *other.o;
	}
}
MPIntPairIntIntMapAdaptor::MPIntPairIntIntMapAdaptor(std::map<int, std::pair<int, int> >* ptr) { o_owned = false; o = ptr; };
MPIntPairIntIntMapAdaptor::~MPIntPairIntIntMapAdaptor() { if (o_owned) delete o; };
int MPIntPairIntIntMapAdaptor::__len__() { return (int)o->size(); };

void MPIntPairIntIntMapAdaptor::__getitem__(int key, MEDPY_NP_OUTPUT(int** int_out_buf, unsigned long long* int_out_buf_len)) {
	if (o->count(key) < 1)
		throw runtime_error(string("Key Error:")+std::to_string(key));
	*int_out_buf = (int*)malloc(sizeof(int)*2);
	*int_out_buf_len = 2;
	auto& v = o->at(key);
	(*int_out_buf)[0] = v.first;
	(*int_out_buf)[1] = v.second;
}

void MPIntPairIntIntMapAdaptor::__setitem__(int key, MEDPY_NP_INPUT(int* int_in_buf, unsigned long long int_in_buf_len)) {
	if (int_in_buf_len <= 1 || int_in_buf == nullptr)
		throw runtime_error("map value type is a 2 item array");
	auto& v = o->operator[](key);
	v.first = int_in_buf[0];
	v.second = int_in_buf[1];
}

void MPIntPairIntIntMapAdaptor::keys(MEDPY_NP_OUTPUT(int** intkeys_out_buf, unsigned long long* intkeys_out_buf_len))
{
	vector<int> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	vector_to_buf(ret, intkeys_out_buf, intkeys_out_buf_len);
};

MPIntPairIntIntMapAdaptor& MPIntPairIntIntMapAdaptor::operator=(const MPIntPairIntIntMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<int, std::pair<int, int> >();
		*o = *(other.o);
	}
	return *this;
}


MPStringUOSetStringMapAdaptor::MPStringUOSetStringMapAdaptor() { o = new std::map<std::string, std::unordered_set<std::string> >(); };
MPStringUOSetStringMapAdaptor::MPStringUOSetStringMapAdaptor(const MPStringUOSetStringMapAdaptor& other) {
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, std::unordered_set<std::string> >();
		*o = *other.o;
	}
};
MPStringUOSetStringMapAdaptor::MPStringUOSetStringMapAdaptor(std::map<std::string, std::unordered_set<std::string> >* ptr) { o_owned = false; o = ptr; };
MPStringUOSetStringMapAdaptor::~MPStringUOSetStringMapAdaptor() { if (o_owned) delete o; };
int MPStringUOSetStringMapAdaptor::__len__() { return (int)o->size(); };

std::vector<std::string> MPStringUOSetStringMapAdaptor::__getitem__(std::string key) {
	vector<string> ret;
	for (auto& s : o->operator[](key)) ret.push_back(s);
	return ret;
};

void MPStringUOSetStringMapAdaptor::__setitem__(std::string key, std::vector<std::string> val) {
	o->operator[](key).clear();
	for (auto& s : val) o->operator[](key).insert(s);
};

std::vector<std::string> MPStringUOSetStringMapAdaptor::keys() {
	vector<string> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};

MPStringUOSetStringMapAdaptor& MPStringUOSetStringMapAdaptor::operator=(const MPStringUOSetStringMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, std::unordered_set<std::string> >();
		*o = *(other.o);
	}
	return *this;
}




/************************************************************************************/

MPIntVecIntMapAdaptor::MPIntVecIntMapAdaptor() { o = new std::map<int, std::vector<int> >(); };
MPIntVecIntMapAdaptor::MPIntVecIntMapAdaptor(const MPIntVecIntMapAdaptor& other) {
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<int, std::vector<int> >();
		*o = *other.o;
	}
};

MPIntVecIntMapAdaptor::MPIntVecIntMapAdaptor(std::map<int, std::vector<int> >* ptr) { o_owned = false; o = ptr; };
MPIntVecIntMapAdaptor::~MPIntVecIntMapAdaptor() { if (o_owned) delete o; };
int MPIntVecIntMapAdaptor::__len__() { return (int)o->size(); };
void MPIntVecIntMapAdaptor::__getitem__(int key, MEDPY_NP_OUTPUT(int** int_out_buf, unsigned long long* int_out_buf_len)) {
	auto& fvec = o->operator[](key);
	vector_to_buf(fvec, int_out_buf, int_out_buf_len);
};
void MPIntVecIntMapAdaptor::__setitem__(int key, MEDPY_NP_INPUT(int* int_in_buf, unsigned long long int_in_buf_len)) {
	vector<int> fvec;
	buf_to_vector(int_in_buf, int_in_buf_len, fvec);
	o->operator[](key) = fvec;
};
std::vector<int> MPIntVecIntMapAdaptor::keys() {
	vector<int> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};

MPIntVecIntMapAdaptor& MPIntVecIntMapAdaptor::operator=(const MPIntVecIntMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<int, std::vector<int> >();
		*o = *(other.o);
	}
	return *this;
}

/*************************************************/

MPIntStringMapAdaptor::MPIntStringMapAdaptor() { o = new std::map<int, std::string>(); };
MPIntStringMapAdaptor::MPIntStringMapAdaptor(const MPIntStringMapAdaptor& other) {
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<int, std::string>();
		*o = *other.o;
	}
}
MPIntStringMapAdaptor::MPIntStringMapAdaptor(std::map<int, string>* ptr) { o_owned = false; o = ptr; };
MPIntStringMapAdaptor::~MPIntStringMapAdaptor() { if (o_owned) delete o; };
int MPIntStringMapAdaptor::__len__() { return (int)o->size(); };
std::string MPIntStringMapAdaptor::__getitem__(int i) { return o->operator[](i); };
void MPIntStringMapAdaptor::__setitem__(int i, const string& val) { o->operator[](i) = val; };
std::vector<int> MPIntStringMapAdaptor::keys()
{
	vector<int> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};

MPIntStringMapAdaptor& MPIntStringMapAdaptor::operator=(const MPIntStringMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<int, std::string>();
		*o = *(other.o);
	}
	return *this;
}

#include <Logger/Logger/Logger.h>

void logger_use_stdout() {
	for (int sect = 0; sect < MAX_LOG_SECTION; ++sect)
		global_logger.init_file(sect, stdout);
}

int val_or_exception(int val, const string& exception_str , int expected_val) {
	string what_msg = exception_str;
	if (exception_str == "")
		what_msg = "Operation failed";
	if (val != expected_val)
		throw runtime_error(what_msg + " (return value=" + std::to_string(val) + ")");
	return val;
}
