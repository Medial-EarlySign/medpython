#ifndef _SERIALIZABLE_OBJECT_H_
#define _SERIALIZABLE_OBJECT_H_

#include "Logger/Logger/Logger.h"
#include "MedUtils/MedUtils/MedIO.h"

#include <cstring>
#include <unordered_set>
#include <set>
#include <unordered_map>

using namespace std;

/** @file
* An Abstract class that can be serialized and written/read from file
*/
class SerializableObject {
public:
	///Relevant for serializations. if changing serialization, increase version number for the 
	///implementing class
	virtual int version() { return  0; }

	// next adds an option to add some actions before a serialization is done
	virtual void pre_serialization() {};

	// Virtual serialization
	virtual size_t get_size() { return 0; } ///<Gets bytes sizes for serializations
	virtual size_t serialize(unsigned char *blob) { return 0; } ///<Serialiazing object to blob memory. return number ob bytes wrote to memory
	virtual size_t deserialize(unsigned char *blob) { return 0; } ///<Deserialiazing blob to object. returns number of bytes read

	// APIs for vectors
	size_t serialize(vector<unsigned char> &blob) { size_t size = get_size(); blob.resize(size); return serialize(&blob[0]); }
	size_t deserialize(vector<unsigned char> &blob) { return deserialize(&blob[0]); }


	//template <class T> void copy_object(T* dst) { 
	//	vector<unsigned int> blob; 
	//	serialize(blob);
	//	dst->deserialize(blob);
	//}

	/// read and deserialize model
	virtual int read_from_file(const string &fname);

	/// serialize model and write to file
	virtual int write_to_file(const string &fname);

	/// Init from string
	int init_from_string(string init_string);
	virtual int init(map<string, string>& map) { return 0; } ///<Virtual to init object from parsed fields
};

/*! @def MEDSERIALIZE_SUPPORT(Type)
* To Join the MedSerialize Wagon :\n
* (1) include this h file, in your h file\n
* (2) implement the get_size, serialize and deserialize functions for your class, you can use MedSerialize functions for that\n
* (3) add the following macro for your class\n
*/
#define MEDSERIALIZE_SUPPORT(Type)																					\
namespace MedSerialize {																							\
	template<> inline size_t get_size<Type>(Type &elem) { return elem.get_size(); }									\
	template<> inline size_t serialize<Type>(unsigned char *blob, Type &elem) { return elem.serialize(blob); }		\
	template<> inline size_t deserialize<Type>(unsigned char *blob, Type &elem) { return elem.deserialize(blob); }	\
}

/*! @def ADD_SERIALIZATION_FUNCS(...)
* To add automatic serialization to your class you can use the following macro with the list of the \n
* variables to serialize inside your class. They all should be MedSerialize supported\n
*/
//size_t serialize(unsigned char *blob) { pre_serialization(); return MedSerialize::serialize(blob, __VA_ARGS__); }		\

#define ADD_SERIALIZATION_FUNCS(...)																\
	size_t get_size() { return MedSerialize::get_size(__VA_ARGS__); }								\
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob,  __VA_ARGS__);}		\
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, __VA_ARGS__); }

#include "SerializableObject_imp.h"

#endif
