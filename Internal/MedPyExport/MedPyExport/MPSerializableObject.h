#ifndef __MED__MPSerializableObject__H__
#define __MED__MPSerializableObject__H__

#include "MedPyCommon.h"

class SerializableObject;

class MPSerializableObject {
public:
	MEDPY_IGNORE(SerializableObject* o);
	MPSerializableObject();
	MEDPY_IGNORE(MPSerializableObject(SerializableObject* from_obj));
	~MPSerializableObject();

	///Relevant for serializations. if changing serialization, increase version number for the 
	///implementing class
	int version();

	/// For better handling of serializations it is highly recommended that each SerializableObject inheriting class will 
	/// implement the next method. It is a must when one needs to support the new_polymorphic method.
	/// One can simply add the macro ADD_CLASS_NAME(class) as a public member in the class.
	string my_class_name();

	/// The names of the serialized fields. Can be helpful for example to print object fields that can be initialized.
	/// can also print helpfull message for init if we keep initialization with the same names.
	//void serialized_fields_name(vector<string> &field_names);

	// next adds an option to add some actions before and after a serialization is done
	//void pre_serialization();
	//void post_deserialization();

	MEDPY_DOC(get_size, "get_size() -> int\n"
		"    Gets bytes sizes for serializations");
	int get_size();

					   // Virtual serialization
	//size_t serialize(unsigned char *blob) { return 0; } ///<Serialiazing object to blob memory. return number ob bytes wrote to memory
	//size_t deserialize(unsigned char *blob) { return 0; } ///<Deserialiazing blob to object. returns number of bytes read

																  // APIs for vectors
	//size_t serialize_vec(vector<unsigned char> &blob) { size_t size = get_size(); blob.resize(size); return serialize(&blob[0]); }
	//size_t deserialize_vec(vector<unsigned char> &blob) { return deserialize(&blob[0]); }
	//virtual size_t serialize(vector<unsigned char> &blob) { return serialize_vec(blob); }
	//virtual size_t deserialize(vector<unsigned char> &blob) { return deserialize_vec(blob); }

	MEDPY_DOC(read_from_file, "read_from_file(fname) -> int\n"
		"    read and deserialize object");
	int read_from_file(const string fname);

	MEDPY_DOC(write_to_file, "write_to_file(fname) -> int\n"
		"    serialize object and write to file");
	int write_to_file(const string fname);

	MEDPY_DOC(read_from_file_unsafe, "read_from_file_unsafe(fname) -> int\n"
		"    read and deserialize model without checking version number - unsafe read");
	int read_from_file_unsafe(const string fname);

	/// Init from string
	//int init_from_string(string init_string);
	//int init_params_from_file(string init_file);
	//int init_param_from_file(string file_str, string &param);
	//virtual int init(map<string, string>& map) { return 0; } ///<Virtual to init object from parsed fields

};



#endif // !__MED__MPSerializableObject__H__
