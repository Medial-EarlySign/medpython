#include "MedPyCommon.h"

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
void MPIntIntMapAdaptor::__setitem__(int i, int val) { o->insert(o->begin(), std::pair<int, int>(i, val)); };
void MPIntIntMapAdaptor::keys(MEDPY_NP_OUTPUT(int** intkeys_out_buf, int* intkeys_out_buf_len)) 
{
	vector<int> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	vector_to_buf(ret, intkeys_out_buf, intkeys_out_buf_len);
};


MPIntStringMapAdaptor::MPIntStringMapAdaptor() { o = new std::map<int, std::string>(); };
MPIntStringMapAdaptor::MPIntStringMapAdaptor(const MPIntStringMapAdaptor& other) {
	o_owned = other.o_owned;
	if (!other.o_owned) {
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
void MPIntStringMapAdaptor::__setitem__(int i, const string& val) { o->insert(o->begin(), std::pair<int, std::string>(i, val)); };
void MPIntStringMapAdaptor::keys(MEDPY_NP_OUTPUT(int** intkeys_out_buf, int* intkeys_out_buf_len))
{
	vector<int> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	vector_to_buf(ret, intkeys_out_buf, intkeys_out_buf_len);
};




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
void MPStringVecFloatMapAdaptor::__getitem__(std::string key, MEDPY_NP_OUTPUT(float** float_out_buf, int* float_out_buf_len)) {
	auto& fvec = o->operator[](key);
	vector_to_buf(fvec, float_out_buf, float_out_buf_len);
};
void MPStringVecFloatMapAdaptor::__setitem__(std::string key, MEDPY_NP_INPUT(float* float_in_buf, int float_in_buf_len)) { 
	vector<float> fvec;
	buf_to_vector(float_in_buf, float_in_buf_len, fvec);
	o->insert(o->begin(), std::pair<string, vector<float> >(key, fvec)); 
};
std::vector<std::string> MPStringVecFloatMapAdaptor::keys() {
	vector<string> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};

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

void MPIntPairIntIntMapAdaptor::__getitem__(int key, MEDPY_NP_OUTPUT(int** int_out_buf, int* int_out_buf_len)) {
	vector<int> ret;
	ret.push_back(o->operator[](key).first);
	ret.push_back(o->operator[](key).second);
	vector_to_buf(ret, int_out_buf, int_out_buf_len);
}

void MPIntPairIntIntMapAdaptor::__setitem__(int key, MEDPY_NP_INPUT(int* int_in_buf, int int_in_buf_len)) {
	vector<int> ret;
	if (int_in_buf_len <= 1 || int_in_buf == nullptr)
		throw runtime_error("map value type is a 2 item array");
	ret.push_back(o->operator[](key).first);
	ret.push_back(o->operator[](key).second);
}

void MPIntPairIntIntMapAdaptor::keys(MEDPY_NP_OUTPUT(int** intkeys_out_buf, int* intkeys_out_buf_len))
{
	vector<int> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	vector_to_buf(ret, intkeys_out_buf, intkeys_out_buf_len);
};



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
	for (auto& s : val) o->operator[](key).insert(s);
};

std::vector<std::string> MPStringUOSetStringMapAdaptor::keys() {
	vector<string> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};



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
MPIntVecIntMapAdaptor::MPIntVecIntMapAdaptor(std::map<int, std::vector<int> >* ptr) { 
	o_owned = false; o = ptr; 
};

MPIntVecIntMapAdaptor::~MPIntVecIntMapAdaptor() { if (o_owned) delete o; };
int MPIntVecIntMapAdaptor::__len__() { return (int)o->size(); };
void MPIntVecIntMapAdaptor::__getitem__(int key, MEDPY_NP_OUTPUT(int** int_out_buf, int* int_out_buf_len)) {
	auto& fvec = o->operator[](key);
	vector_to_buf(fvec, int_out_buf, int_out_buf_len);
};
void MPIntVecIntMapAdaptor::__setitem__(int key, MEDPY_NP_INPUT(int* int_in_buf, int int_in_buf_len)) {
	vector<int> fvec;
	buf_to_vector(int_in_buf, int_in_buf_len, fvec);
	o->insert(o->begin(), std::pair<int, vector<int> >(key, fvec));
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

