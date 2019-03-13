#include "MPSerializableObject.h"
#include <SerializableObject/SerializableObject/SerializableObject.h>

MPSerializableObject::MPSerializableObject() { o = nullptr; }
MPSerializableObject::MPSerializableObject(SerializableObject* from_obj) { o = from_obj; }
MPSerializableObject::~MPSerializableObject() {/*do nothing*/}
int MPSerializableObject::version() { if (o==nullptr) throw runtime_error("ERROR: SerializableObject can only be created by a .asSerializable() method"); return o->version(); }
string MPSerializableObject::my_class_name() { if (o == nullptr) throw runtime_error("ERROR: SerializableObject can only be created by a .asSerializable() method"); return o->my_class_name(); }
int MPSerializableObject::get_size() { if (o == nullptr) throw runtime_error("ERROR: SerializableObject can only be created by a .asSerializable() method"); return (int)o->get_size(); }
int MPSerializableObject::read_from_file(const string fname) { if (o == nullptr) throw runtime_error("ERROR: SerializableObject can only be created by a .asSerializable() method"); return val_or_exception(o->read_from_file(fname),string("Error: SerializableObject::read_from_file() failed reading file ")+ fname); }
int MPSerializableObject::write_to_file(const string fname) { if (o == nullptr) throw runtime_error("ERROR: SerializableObject can only be created by a .asSerializable() method"); return val_or_exception(o->write_to_file(fname), string("Error: SerializableObject::write_to_file() failed writing file ") + fname); }
int MPSerializableObject::read_from_file_unsafe(const string fname){ if (o == nullptr) throw runtime_error("ERROR: SerializableObject can only be created by a .asSerializable() method"); return o->read_from_file_unsafe(fname); }
