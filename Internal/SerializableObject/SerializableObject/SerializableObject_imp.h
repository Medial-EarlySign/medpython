#ifndef __SERIALIZABLE_OBJECT_IMP_LIB_H__
#define __SERIALIZABLE_OBJECT_IMP_LIB_H__

#include <type_traits>
#include <typeinfo>

#if __GNUC__
#if __cplusplus <= 201402L
template< bool B, class T = void >
using enable_if_t = typename enable_if<B, T>::type;
#endif
#endif // !_

namespace MedSerialize {

	// Helper function to parse names and clean spaces, tabs, etc
	static void get_list_names(char *list_of_args, vector<string> &names) {
		// get names
		string s(list_of_args);
		boost::replace_all(s, " ", "");
		boost::replace_all(s, "\t", "");
		boost::replace_all(s, "\n", "");
		boost::split(names, s, boost::is_any_of(","));
	}

	//========================================================================================
	// FUNCTION FROM BELOW THAT IMPLEMENT AN STD CONTAINER (map, vector, etc) 
	// MUST BE PRE DECLARED HERE BEFORE THE IMPLEMENTATION
	// TO MAKE SURE THE CORRECT RECURSIVE RESOLVING IS DONE
	// (simple, string & variadic need not - and will create a compilation error)
	//========================================================================================

	template <class T>  size_t get_size(T *&v);
	template <class T> size_t serialize(unsigned char *blob, T *&v);
	template <class T> size_t deserialize(unsigned char *blob, T *&v);
	template <class T> string object_json(T *&v);
	template <class T> string object_json(const T *&v);
	template <class T> string object_json(T * const &v);
	template <class T>  size_t get_size(vector<T> &v);
	template <class T> size_t serialize(unsigned char *blob, vector<T> &v);
	template <class T> size_t deserialize(unsigned char *blob, vector<T> &v);
	template <class T> string object_json(const vector<const T> &v);
	template <class T> string object_json(const vector<T> &v);
	template <class T> string object_json(vector<const T> &v);
	template <class T, class S> size_t get_size(pair<T, S> &v);
	template <class T, class S> size_t serialize(unsigned char *blob, pair<T, S> &v);
	template <class T, class S> size_t deserialize(unsigned char *blob, pair<T, S> &v);
	template <class T, class S> string object_json(pair<const T, const S> &v);
	template <class T, class S> string object_json(const pair<T, S> &v);
	template <class T, class S> string object_json(const pair<const T, const S> &v);
	template <class T, class S> size_t get_size(map<T, S> &v);
	template <class T, class S> size_t serialize(unsigned char *blob, map<T, S> &v);
	template <class T, class S> size_t deserialize(unsigned char *blob, map<T, S> &v);
	template <class T, class S> string object_json(map<const T, const S> &v);
	template <class T, class S> string object_json(const map<const T, const S> &v);
	template <class T, class S> string object_json(const map<T, S> &v);
	template <class T, class S> size_t get_size(unordered_map<T, S> &v);
	template <class T, class S> size_t serialize(unsigned char *blob, unordered_map<T, S> &v);
	template <class T, class S> size_t deserialize(unsigned char *blob, unordered_map<T, S> &v);
	template <class T, class S> string object_json(unordered_map<const T, const S> &v);
	template <class T, class S> string object_json(const unordered_map<const T, const S> &v);
	template <class T, class S> string object_json(const unordered_map<T, S> &v);
	template <class T> size_t get_size(unordered_set<T> &v);
	template <class T> size_t serialize(unsigned char *blob, unordered_set<T> &v);
	template <class T> size_t deserialize(unsigned char *blob, unordered_set<T> &v);
	template <class T> string object_json(unordered_set<const T> &v);
	template <class T> string object_json(const unordered_set<const T> &v);
	template <class T> string object_json(const unordered_set<T> &v);
	template <class T> size_t get_size(set<T> &v);
	template <class T> size_t serialize(unsigned char *blob, set<T> &v);
	template <class T> size_t deserialize(unsigned char *blob, set<T> &v);
	template <class T> string object_json(set<const T> &v);
	template <class T> string object_json(const set<const T> &v);
	template <class T> string object_json(const set<T> &v);
	template <class T> size_t get_size(unique_ptr<T> &v);
	template <class T> size_t serialize(unsigned char *blob, unique_ptr<T> &v);
	template <class T> size_t deserialize(unsigned char *blob, unique_ptr<T> &v);
	template <class T> string object_json(unique_ptr<T> &v);
	template <class T> string object_json(unique_ptr<const T> &v);

	//========================================================================================
	// IMPLEMANTATIONS 
	//========================================================================================

	//.........................................................................................
	// templated simple ones : int, float, long, double, etc...
	//.........................................................................................
	template <class T> size_t get_size(T &elem)
	{
		//cout << "inside simple getsize with type " << typeid(T).name() << endl;
		return sizeof(T);
	}

	//.........................................................................................
	template <class T> size_t serialize(unsigned char *blob, T &elem)
	{
		//cout << "inside simple serialize with type " << typeid(T).name() << " with sizeof " << sizeof(T) << endl;
		memcpy(blob, &elem, sizeof(T));
		//unsigned char *pelem = (unsigned char *)(&elem);
		//SRL_LOG("simple serialize: %02x %02x %02x %02x\n", blob[0], blob[1], blob[2], blob[3]);
		return sizeof(T);
	}

	//.........................................................................................
	template <class T> size_t deserialize(unsigned char *blob, T &elem)
	{
		//cerr << "inside simple deserialize with type " << typeid(T).name() << " with sizeof " << sizeof(T) << endl;
		memcpy(&elem, blob, sizeof(T));
		//unsigned char *pelem = (unsigned char *)(&elem);
		//SRL_LOG("simple deserialize: %02x %02x %02x %02x\n", blob[0], blob[1], blob[2], blob[3]);
		return sizeof(T);
	}

	template<class T>
	string object_json_spec(T& v, std::true_type /* is_enum */, std::false_type) {
		stringstream str;
		//it's object:
		string enum_nm = typeid(v).name();
		//remove number prefix:
		enum_nm = boost::regex_replace(enum_nm, boost::regex("^([0-9]*|enum )"), "");
		str << enum_nm << "::" << (int)v;
		return str.str();
	}
	template<class T>
	string object_json_spec(T& v, std::false_type, std::true_type /* SerializableObject */) {
		stringstream str;
		str << v.object_json();
		return str.str();
	}
	template<class T>
	string object_json_spec(T& v, std::false_type, std::false_type) {
		// neither
		stringstream str;
		str << "UNSUPPORTED::" << typeid(v).name();
		return str.str();
	}

	template<class T>
	string object_json(T& v) {
		return object_json_spec(v, std::is_enum<T>{}, std::is_base_of<SerializableObject, T>{});
	}

	//.........................................................................................
	// T * case : will ONLY Work for classes implementing the new_polymorphic, and my_class_name methods !!
	// Also assumes a SINGLE new (not an array).
	//.........................................................................................
	template <class T> size_t get_size(T *&elem)
	{
		//cerr << "get size of T *\n";
		size_t size = 0;

		if (elem == NULL) {
			// account for "NULL" string when pointer is null
			string s = "NULL";
			size += MedSerialize::get_size(s);
		}
		else {
			string s = elem->my_class_name(); // help compiler
											  //cerr << "get size of T * (2) : class" << s << "\n";
			size += MedSerialize::get_size(s);		// account for class name
			size += MedSerialize::get_size((*elem)); // account for class serialization
		}
		//cerr << "get size of T * (3) : size" << size << "\n";
		return size;
	}

	//.........................................................................................
	template <class T> size_t serialize(unsigned char *blob, T *&elem)
	{
		//cerr << "Serializing T * for " << (elem != NULL) ? elem->my_class_name().c_str() : string("NULL").c_str() << "\n";
		size_t pos = 0;

		// serializing name of class, for polymorphic support
		// we will detect NULL pointers by simply writing "NULL" into the serialization
		string s = "NULL";
		if (elem != NULL)
			s = elem->my_class_name(); // help compiler
		pos += MedSerialize::serialize<string>(blob + pos, s);

		// serializing the actual class
		//cerr << "Serializing T * (2) : class name is " << s << " pos " << pos << "\n";
		if (elem != NULL)
			pos += MedSerialize::serialize(blob + pos, (*elem));

		//cerr << "Serializing T * (4) : pos " << pos << "\n";
		return pos;
	}

	//.........................................................................................
	template <class T> size_t deserialize(unsigned char *blob, T *&elem)
	{
		//cerr << "DeSerializing T * \n";

		size_t pos = 0;

		// deserialize name of class
		string cname;
		pos += MedSerialize::deserialize<string>(blob + pos, cname);

		//cerr << "Deserializing T * (2) : Got class name " << cname << " pos " << pos << "\n";

		if (cname == "NULL") {
			elem = NULL;
		}
		else {
			// heart of matters: doing the right new operation
			T dummy; // we need access to the new_polymorphic method
			elem = (T *)dummy.new_polymorphic(cname);
			//cerr << "Deserializing T * (3) : elem is " << elem << "\n";
			if (elem == NULL) {
				elem = new T;
			}

			//cerr << "Deserializing T * (4) : elem is " << elem << "\n";

			// now we are ready to deserialize T or its derived
			pos += MedSerialize::deserialize(blob + pos, (*elem));
		}

		//cerr << "Deserializing T * (4) : pos " << pos << "\n";
		return pos;
	}

	template <class T> string object_json(T *&v) {
		stringstream str;
		//it's object:
		if (v != NULL) {
			str << MedSerialize::object_json((const T&)*v);
		}
		else
			str << "NULL";
		return str.str();
	}
	template <class T> string object_json(const T *&v) {
		stringstream str;
		if (v != NULL) {
			//it's object:
			str << MedSerialize::object_json((const T&)*v);
		}
		else
			str << "NULL";
		return str.str();
	}
	template <class T> string object_json(T * const &v) {
		stringstream str;
		if (v != NULL) {
			//it's object:
			str << MedSerialize::object_json((const T&)*v);
		}
		else
			str << "NULL";
		return str.str();
	}



	template <class T> size_t get_size(unique_ptr<T> &elem)
	{
		//cerr << "get size of T *\n";
		size_t size = 0;

		if (elem.get() == NULL) {
			// account for "NULL" string when pointer is null
			string s = "NULL";
			size += MedSerialize::get_size(s);
		}
		else {
			string s = elem->my_class_name(); // help compiler
											  //cerr << "get size of T * (2) : class" << s << "\n";
			size += MedSerialize::get_size(s);		// account for class name
			size += MedSerialize::get_size((*elem)); // account for class serialization
		}
		//cerr << "get size of T * (3) : size" << size << "\n";
		return size;
	}

	//.........................................................................................
	template <class T> size_t serialize(unsigned char *blob, unique_ptr<T> &elem)
	{
		//cerr << "Serializing T * for " << (elem != NULL) ? elem->my_class_name().c_str() : string("NULL").c_str() << "\n";
		size_t pos = 0;

		// serializing name of class, for polymorphic support
		// we will detect NULL pointers by simply writing "NULL" into the serialization
		string s = "NULL";
		if (elem.get() != NULL)
			s = elem->my_class_name(); // help compiler
		pos += MedSerialize::serialize<string>(blob + pos, s);

		// serializing the actual class
		//cerr << "Serializing T * (2) : class name is " << s << " pos " << pos << "\n";
		if (elem.get() != NULL)
			pos += MedSerialize::serialize(blob + pos, (*elem));

		//cerr << "Serializing T * (4) : pos " << pos << "\n";
		return pos;
	}

	//.........................................................................................
	template <class T> size_t deserialize(unsigned char *blob, unique_ptr<T> &elem)
	{
		//cerr << "DeSerializing T * \n";

		size_t pos = 0;

		// deserialize name of class
		string cname;
		pos += MedSerialize::deserialize<string>(blob + pos, cname);

		//cerr << "Deserializing T * (2) : Got class name " << cname << " pos " << pos << "\n";

		if (cname == "NULL") {
			elem = NULL;
		}
		else {
			// heart of matters: doing the right new operation
			T dummy; // we need access to the new_polymorphic method
			elem = unique_ptr<T>((T *)dummy.new_polymorphic(cname));
			//cerr << "Deserializing T * (3) : elem is " << elem << "\n";
			if (elem.get() == NULL) {
				elem = unique_ptr<T>(new T);
			}

			//cerr << "Deserializing T * (4) : elem is " << elem << "\n";

			// now we are ready to deserialize T or its derived
			pos += MedSerialize::deserialize(blob + pos, (*elem.get()));
		}

		//cerr << "Deserializing T * (4) : pos " << pos << "\n";
		return pos;
	}

	template <class T> string object_json(unique_ptr<T> &v) {
		stringstream str;
		//it's object:
		if (v.get() != NULL) {
			str << MedSerialize::object_json((const T&)*v);
		}
		else
			str << "NULL";
		return str.str();
	}
	template <class T> string object_json(unique_ptr<const T> &v) {
		stringstream str;
		//it's object:
		if (v.get() != NULL) {
			str << MedSerialize::object_json((const T&)*v);
		}
		else
			str << "NULL";
		return str.str();
	}


	//.........................................................................................
	// string is special
	//.........................................................................................
	template<> inline size_t get_size<string>(string &str) {
		size_t size = 0;
		size += sizeof(size_t); // length of string
		size += str.length() + 1;
		return size;
	}

	//.........................................................................................
	template<> inline size_t serialize<string>(unsigned char *blob, string &str)
	{
		size_t pos = 0;
		size_t len = str.length();
		//fprintf(stderr, "string serialize(%d) %s\n", len, str.c_str());
		memcpy(blob, &len, sizeof(size_t)); pos += sizeof(size_t);
		memcpy(blob + pos, str.c_str(), len); pos += len;
		blob[pos] = 0; pos++;
		//fprintf(stderr, "string serialize(%d) %s\n", len, &blob[sizeof(size_t)]);
		return pos;
	}

	//.........................................................................................
	template<> inline size_t deserialize<string>(unsigned char *blob, string &str)
	{
		//fprintf(stderr, "string deserialize\n");
		size_t pos = 0;
		size_t len;
		memcpy(&len, blob, sizeof(size_t)); pos += sizeof(size_t);
		//fprintf(stderr, "string deserialize pos %d len %d\n", pos, len);
		string new_s((char *)&blob[pos]);
		str = new_s;
		//fprintf(stderr, "string deserialize pos %d :: %s\n", pos, str.c_str());
		pos += len + 1;
		return pos;
	}

	//.........................................................................................
	// vector of type T that has a MedSerialize function
	//.........................................................................................
	template<class T> size_t get_size(vector<T> &v)
	{
		//cout << "inside vector getsize with type " << typeid(T).name() << endl;
		size_t size = 0, len = v.size();
		size += MedSerialize::get_size<size_t>(len);
		for (T &elem : v)
			size += MedSerialize::get_size(elem);

		return size;
	}

	//.........................................................................................
	template <class T> size_t serialize(unsigned char *blob, vector<T> &v)
	{
		//fprintf(stderr, "vector serialize\n");
		//cout << "inside vector serialize with type " << typeid(T).name() << endl;
		size_t pos = 0, len = v.size();
		pos += MedSerialize::serialize<size_t>(blob + pos, len);
		if (len > 0)
			for (T &elem : v)
				pos += MedSerialize::serialize(blob + pos, elem);
		return pos;
	}

	//.........................................................................................
	template <class T> size_t deserialize(unsigned char *blob, vector<T> &v)
	{
		//fprintf(stderr, "vector deserialize\n");
		//cout << "inside vector deserialize with type " << typeid(T).name() << endl;
		size_t pos = 0, len;
		pos += MedSerialize::deserialize<size_t>(blob + pos, len);
		if (len != v.size()) v.clear();
		if (len > 0) {
			v.resize(len);
			for (T &elem : v)
				pos += MedSerialize::deserialize(blob + pos, elem);
		}
		return pos;
	}

	template <class T> string object_json(const vector<const T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		for (size_t i = 0; i < v.size(); ++i)
		{
			if (i > 0)
				str << ",";
			str << MedSerialize::object_json(v[i]);
		}
		str << "]";
		return str.str();
	}
	template <class T> string object_json(const vector<T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		for (size_t i = 0; i < v.size(); ++i)
		{
			if (i > 0)
				str << ",";
			str << MedSerialize::object_json(v[i]);
		}
		str << "]";
		return str.str();
	}
	template <class T> string object_json(vector<const T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		for (size_t i = 0; i < v.size(); ++i)
		{
			if (i > 0)
				str << ",";
			str << MedSerialize::object_json(v[i]);
		}
		str << "]";
		return str.str();
	}


	//.........................................................................................
	// set<T> with a MedSerialize function
	//.........................................................................................
	template <class T> size_t get_size(set<T> &v)
	{
		size_t size = 0, len = v.size();
		size += MedSerialize::get_size<size_t>(len);

		for (typename set<T>::iterator it = v.begin(); it != v.end(); ++it)
			size += MedSerialize::get_size(*it);

		return size;
	}

	//.........................................................................................
	template <class T> size_t serialize(unsigned char *blob, set<T> &v)
	{
		//fprintf(stderr, "map serialize\n");
		size_t pos = 0, len = v.size();
		pos += MedSerialize::serialize<size_t>(blob + pos, len);

		if (len > 0) {
			for (typename set<T>::iterator it = v.begin(); it != v.end(); ++it)
				pos += MedSerialize::serialize(blob + pos, (T &)(*it));
		}

		return pos;
	}

	//.........................................................................................
	template <class T> size_t deserialize(unsigned char *blob, set<T> &v)
	{
		size_t pos = 0, len;
		pos += MedSerialize::deserialize<size_t>(blob + pos, len);
		v.clear();
		T elem;
		if (len > 0) {
			for (int i = 0; i < len; i++) {
				pos += MedSerialize::deserialize(blob + pos, elem);
				v.insert(elem);
			}
		}
		return pos;
	}

	template <class T> string object_json(set<const T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (typename set<T>::iterator it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			str << MedSerialize::object_json((const T &)(*it));
			start = false;
		}
		str << "]";
		return str.str();
	}
	template <class T> string object_json(const set<const T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (typename set<T>::iterator it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			str << MedSerialize::object_json((const T &)(*it));
			start = false;
		}
		str << "]";
		return str.str();
	}
	template <class T> string object_json(const set<T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (typename set<T>::iterator it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			str << MedSerialize::object_json((T &)(*it));
			start = false;
		}
		str << "]";
		return str.str();
	}


	//.........................................................................................
	// unordered_set<T> with a MedSerialize function
	//.........................................................................................
	template <class T> size_t get_size(unordered_set<T> &v)
	{
		size_t size = 0, len = v.size();
		size += MedSerialize::get_size<size_t>(len);
		for (T elem : v)
			size += MedSerialize::get_size(elem);

		return size;
	}

	//.........................................................................................
	template <class T> size_t serialize(unsigned char *blob, unordered_set<T> &v)
	{
		//fprintf(stderr, "map serialize\n");
		size_t pos = 0, len = v.size();
		pos += MedSerialize::serialize<size_t>(blob + pos, len);
		if (len > 0)
			for (T elem : v) {
				pos += MedSerialize::serialize(blob + pos, elem);
			}
		return pos;
	}

	//.........................................................................................
	template <class T> size_t deserialize(unsigned char *blob, unordered_set<T> &v)
	{
		size_t pos = 0, len;
		pos += MedSerialize::deserialize<size_t>(blob + pos, len);
		v.clear();
		T elem;
		if (len > 0) {
			for (int i = 0; i < len; i++) {
				pos += MedSerialize::deserialize(blob + pos, elem);
				v.insert(elem);
			}
		}
		return pos;
	}

	template <class T> string object_json(unordered_set<const T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (const T &e : v)
		{
			if (!start)
				str << ",";
			str << MedSerialize::object_json(e);
			start = false;
		}
		str << "]";
		return str.str();
	}
	template <class T> string object_json(const unordered_set<const T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (const T &e : v)
		{
			if (!start)
				str << ",";
			str << MedSerialize::object_json(e);
			start = false;
		}
		str << "]";
		return str.str();
	}
	template <class T> string object_json(const unordered_set<T> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (const T &e : v)
		{
			if (!start)
				str << ",";
			str << MedSerialize::object_json(e);
			start = false;
		}
		str << "]";
		return str.str();
	}



	//.........................................................................................
	// pair<T,S> both with a MedSerialize function
	//.........................................................................................
	template <class T, class S> size_t get_size(pair<T, S> &v)
	{
		size_t size = 0;
		T *t = (T *)&v.first;
		S *s = (S *)&v.second;
		size += MedSerialize::get_size((*t));
		size += MedSerialize::get_size((*s));

		return size;
	}

	//.........................................................................................
	template <class T, class S> size_t serialize(unsigned char *blob, pair<T, S> &v)
	{
		//fprintf(stderr, "map serialize\n");
		size_t pos = 0;
		T *t = (T *)&v.first;
		S *s = (S *)&v.second;
		pos += MedSerialize::serialize(blob + pos, (*t));
		pos += MedSerialize::serialize(blob + pos, (*s));

		return pos;
	}

	//.........................................................................................
	template <class T, class S> size_t deserialize(unsigned char *blob, pair<T, S> &v)
	{
		//fprintf(stderr, "map deserialize\n");
		size_t pos = 0;
		T t;
		S s;
		pos += MedSerialize::deserialize(blob + pos, t);
		pos += MedSerialize::deserialize(blob + pos, s);
		v.first = t;
		v.second = s;

		return pos;
	}

	template <class T, class S> string object_json(const pair<const T, const S> &v) {
		stringstream str;
		//it's object: - TODO: print type name
		str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
		str << "[";

		str << MedSerialize::object_json(v.first);
		str << "," << MedSerialize::object_json(v.second);

		str << "]}";
		return str.str();
	}
	template <class T, class S> string object_json(pair<const T, const S> &v) {
		stringstream str;
		//it's object: - TODO: print type name
		str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
		str << "[";

		str << MedSerialize::object_json(v.first);
		str << "," << MedSerialize::object_json(v.second);

		str << "]}";
		return str.str();
	}
	template <class T, class S> string object_json(const pair<T, S> &v) {
		stringstream str;
		//it's object: - TODO: print type name
		str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
		str << "[";

		str << MedSerialize::object_json(v.first);
		str << "," << MedSerialize::object_json(v.second);

		str << "]}";
		return str.str();
	}

	//.........................................................................................
	// map<T,S> both with a MedSerialize function
	//.........................................................................................
	template <class T, class S> size_t get_size(map<T, S> &v)
	{
		size_t size = 0, len = v.size();
		size += MedSerialize::get_size<size_t>(len);
		for (auto &elem : v) {
			T *t = (T *)&elem.first;
			S *s = (S *)&elem.second;
			size += MedSerialize::get_size((*t));
			size += MedSerialize::get_size((*s));
		}
		return size;
	}

	//.........................................................................................
	template <class T, class S> size_t serialize(unsigned char *blob, map<T, S> &v)
	{
		//fprintf(stderr, "map serialize\n");
		size_t pos = 0, len = v.size();
		pos += MedSerialize::serialize<size_t>(blob + pos, len);
		if (len > 0)
			for (auto &elem : v) {
				T *t = (T *)&elem.first;
				S *s = (S *)&elem.second;
				pos += MedSerialize::serialize(blob + pos, (*t));
				pos += MedSerialize::serialize(blob + pos, (*s));
			}
		return pos;
	}

	//.........................................................................................
	template <class T, class S> size_t deserialize(unsigned char *blob, map<T, S> &v)
	{
		//fprintf(stderr, "map deserialize\n");
		size_t pos = 0, len;
		pos += MedSerialize::deserialize(blob + pos, len);
		v.clear();
		T t;
		S s;
		if (len > 0) {
			for (int i = 0; i < len; i++) {
				pos += MedSerialize::deserialize(blob + pos, t);
				pos += MedSerialize::deserialize(blob + pos, s);
				v[t] = s;
			}
		}
		return pos;
	}

	template <class T, class S> string object_json(map<const T, const S> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (auto it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			//print pair:

			str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
			str << "[";

			str << MedSerialize::object_json((const T&)it->first);
			str << "," << MedSerialize::object_json((const S &)it->second);

			str << "]}";
			start = false;

		}
		str << "]";
		return str.str();
	}
	template <class T, class S> string object_json(const map<const T, const S> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (auto it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			//print pair:

			str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
			str << "[";

			str << MedSerialize::object_json((const T&)it->first);
			str << "," << MedSerialize::object_json((const S &)it->second);

			str << "]}";
			start = false;

		}
		str << "]";
		return str.str();
	}
	template <class T, class S> string object_json(const map<T, S> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (auto it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			//print pair:

			str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
			str << "[";

			str << MedSerialize::object_json((T&)it->first);
			str << "," << MedSerialize::object_json((S &)it->second);

			str << "]}";
			start = false;

		}
		str << "]";
		return str.str();
	}


	//.........................................................................................
	// unordered_map<T,S> both with a MedSerialize function
	//.........................................................................................
	template <class T, class S> size_t get_size(unordered_map<T, S> &v)
	{
		size_t size = 0, len = v.size();
		size += MedSerialize::get_size<size_t>(len);

		for (typename unordered_map<T, S>::iterator it = v.begin(); it != v.end(); ++it) {
			size += MedSerialize::get_size((T &)(it->first));
			size += MedSerialize::get_size((S &)(it->second));
		}

		return size;
	}

	//.........................................................................................
	template <class T, class S> size_t serialize(unsigned char *blob, unordered_map<T, S> &v)
	{
		size_t pos = 0, len = v.size();
		pos += MedSerialize::serialize<size_t>(blob + pos, len);

		if (len > 0) {
			for (typename unordered_map<T, S>::iterator it = v.begin(); it != v.end(); ++it) {
				pos += MedSerialize::serialize(blob + pos, (T &)(it->first));
				pos += MedSerialize::serialize(blob + pos, (S &)(it->second));
			}
		}
		return pos;
	}

	//.........................................................................................
	template <class T, class S> size_t deserialize(unsigned char *blob, unordered_map<T, S> &v)
	{
		size_t pos = 0, len;
		pos += MedSerialize::deserialize(blob + pos, len);
		v.clear();
		T t;
		S s;
		if (len > 0) {
			for (int i = 0; i < len; i++) {
				pos += MedSerialize::deserialize(blob + pos, t);
				pos += MedSerialize::deserialize(blob + pos, s);
				v[t] = s;
			}
		}
		return pos;
	}

	template <class T, class S> string object_json(unordered_map<const T, const S> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (auto it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			//print pair:

			str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
			str << "[";

			str << MedSerialize::object_json((const T &)it->first);
			str << "," << MedSerialize::object_json((const S &)it->second);

			str << "]}";
			start = false;

		}
		str << "]";
		return str.str();
	}
	template <class T, class S> string object_json(const unordered_map<const T, const S> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (auto it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			//print pair:

			str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
			str << "[";

			str << MedSerialize::object_json((const T &)it->first);
			str << "," << MedSerialize::object_json((const S &)it->second);

			str << "]}";
			start = false;

		}
		str << "]";
		return str.str();
	}
	template <class T, class S> string object_json(const unordered_map<T, S> &v) {
		stringstream str;
		//it's object:
		str << "[";
		bool start = true;
		for (auto it = v.begin(); it != v.end(); ++it)
		{
			if (!start)
				str << ",";
			//print pair:

			str << "{ \"Object\":\"pair<" << typeid(T).name() << ", " << typeid(S).name() << ">\", ";
			str << "[";

			str << MedSerialize::object_json((T &)it->first);
			str << "," << MedSerialize::object_json((S &)it->second);

			str << "]}";
			start = false;

		}
		str << "]";
		return str.str();
	}


	//.........................................................................................
	// Variadic Wrappers to call several elements is a single line
	//.........................................................................................
	template <class T, class... Ts> size_t get_size(T& elem, Ts&... args)
	{
		size_t size = 0;

		size += MedSerialize::get_size(elem);
		size += MedSerialize::get_size(args...);

		return size;

	}

	//.........................................................................................
	template <class T, class... Ts> size_t serialize(unsigned char *blob, T& elem, Ts&... args)
	{
		size_t pos = 0;

		pos += MedSerialize::serialize(blob, elem);
		pos += MedSerialize::serialize(blob + pos, args...);

		return pos;

	}

	//.........................................................................................
	template <typename T, typename... Ts> size_t deserialize(unsigned char *blob, T& elem, Ts&... args)
	{
		size_t pos = 0;

		pos += MedSerialize::deserialize(blob, elem);
		pos += MedSerialize::deserialize(blob + pos, args...);

		return pos;

	}

	template <class T, class... Ts> string object_json(T& elem, Ts&... args) {
		stringstream str;
		//it's object:
		str << MedSerialize::object_json(elem);
		str << "," << MedSerialize::object_json(args...);

		return str.str();
	}


	//====================================================================================================
	// Next is the implementation of the serialization/deserialization process
	//====================================================================================================
	// last stage serializer :  the last element : size, name, serialization
	template<class T> size_t serializer(unsigned char *blob, int counter, vector<string> &names, T &elem)
	{
		//SRL_LOG("last_serializaer: blob %x counter %d name %s\n", blob, counter, names[counter].c_str());

		// keep size_t place for size
		size_t pos = 0;
		pos += MedSerialize::serialize(blob + pos, pos);

		//SRL_LOG("last_serializer(1) :  blob %x size %d name %s\n", blob, pos, names[counter].c_str());
		// serialize string : the name of the variable
		pos += MedSerialize::serialize(blob + pos, names[counter]);

		//SRL_LOG("last_serializer(2) :  blob %x size %d name %s\n", blob, pos, names[counter].c_str());
		// serialize elem
		pos += MedSerialize::serialize(blob + pos, elem);

		//SRL_LOG("last_serializer(3) :  blob %x size %d name %s\n", blob, pos, names[counter].c_str());

		// going back to 0 and writing the size
		MedSerialize::serialize(blob, pos);

		// return serialized size
		return pos;
	}


	// mid level serializer : looping through the elements to serialize
	template<class T, class... Ts> size_t serializer(unsigned char *blob, int counter, vector<string> &names, T &elem, Ts&...args)
	{
		//SRL_LOG("mid_serializer: name %d is %s\n", counter, names[counter].c_str());

		size_t pos = 0;
		// serialize elem
		pos += MedSerialize::serializer(blob, counter, names, elem);

		//SRL_LOG("mid(1) name %s blob %x pos %d\n", names[counter].c_str(), blob, pos);

		// recursive call to the next elements
		pos += MedSerialize::serializer(blob + pos, counter + 1, names, args...);

		//SRL_LOG("mid(2) name %s blob %x pos %d\n", names[counter].c_str(), blob, pos);

		// return serialized size
		return pos;
	}

	template<typename T, enable_if_t<std::is_base_of<SerializableObject, T>::value, int> = 0> string get_name() {
		unique_ptr<T> cl = unique_ptr<T>(new T); //Not using object
		string cl_name = cl->my_class_name();

		return cl_name;
	}

	template<typename T, enable_if_t<!std::is_base_of<SerializableObject, T>::value, int> = 0> string get_name() {
		string cl_name = "primitive_type";

		return cl_name;
	}

	// last stage deserializer :  the last element : size, name, serialization
	template<class T> size_t deserializer(unsigned char *blob, int counter, vector<string> &names, map <string, size_t> &name2pos, T &elem)
	{
		string name = names[counter];
		//cerr << "last_deserializer: name " << counter << " is " << name << "\n";

		size_t pos = 0;
		// deserialize only if appears in map, otherwise warn
		if (name2pos.find(name) != name2pos.end()) {
			//cerr << "last_deserializer: name " << name << " name2pos " << name2pos[name] << "\n";
			pos += MedSerialize::deserialize(blob + name2pos[name], elem);
			//unsigned char *p = blob + name2pos[name];
			//SRL_LOG("=====> size is %02x %02x %02x %02x\n", p[0], p[1], p[2], p[3]);
		}
		else {
			string cl_name = get_name<T>();
			cerr << "WARNING: In \"" << cl_name << "\" element " << name << " not serialized... will be deserialized to its default\n";
		}

		return pos;
	}


	// mid level deserializer : looping through the elements to serialize
	template<class T, class... Ts> size_t deserializer(unsigned char *blob, int counter, vector<string> &names, map<string, size_t> &name2pos, T &elem, Ts&...args)
	{
		string name = names[counter];
		//cerr << "mid_deserializer: name " << counter << " is " << name << "\n";

		size_t pos = 0;

		// deserialize elem , note we don't use blob+pos , this is due to the usage of name2pos to get exact positions
		pos += MedSerialize::deserializer(blob, counter, names, name2pos, elem);

		// recursive call to the next elements
		pos += MedSerialize::deserializer(blob, counter + 1, names, name2pos, args...);

		// return serialized size
		return pos;
	}


	// last stage get_sizer :  the last element : size, name, serialization
	template<class T> size_t get_sizer(int counter, vector<string> &names, T &elem)
	{
		string name = names[counter];
		//cerr << "last_serializer: name " << counter << " is " << name << "\n";

		size_t size = 0;

		size += MedSerialize::get_size(size); // account for size at start
		size += MedSerialize::get_size(name); // account for string name
		size += MedSerialize::get_size(elem); // account for actual element serialization

		return size;
	}


	// mid level get_sizer : looping through the elements to serialize
	template<class T, class... Ts> size_t get_sizer(int counter, vector<string> &names, T &elem, Ts&...args)
	{
		string name = names[counter];
		//cerr << "mid_get_sizer: name " << counter << " is " << name << "\n";

		size_t size = 0;

		// deserialize elem
		size += MedSerialize::get_sizer(counter, names, elem);

		// recursive call to the next elements
		size += MedSerialize::get_sizer(counter + 1, names, args...);

		return size;
	}


	// get_size before serialization , for top levels
	template<class... Ts> size_t get_size_top(char *list_of_args, Ts&...args) {

		//cerr << "get_size_top\n";
		size_t size = 0;

		// get names
		vector<string> names;
		MedSerialize::get_list_names(list_of_args, names);

		size += MedSerialize::get_size(size); // count for the serialization size kept at the beginning
		size += MedSerialize::get_sizer(0, names, args...);

		//cerr << "get_size_top returned size is " << size << "\n";
		return size;
	}

	// preparing for serializer chain of calls: getting names out of the args list.
	// adding and filling in the total size at start (will be needed in deserialize_top)
	template<class... Ts> size_t serialize_top(unsigned char *blob, char *list_of_args, Ts&...args) {

		//cerr << "serialize_top\n";

		// get names
		vector<string> names;
		MedSerialize::get_list_names(list_of_args, names);

		// save place for size
		size_t pos = 0;
		pos += MedSerialize::serialize(blob + pos, pos);

		// serialize (recursively)
		pos += MedSerialize::serializer(blob + pos, 0, names, args...);

		// write back the size at pos 0
		//cerr << "stop: writing pos " << pos << " at start\n";
		MedSerialize::serialize(blob, pos);

		// return the size
		return pos;
	}

	// preparing for deserializer: (1) get names (2) prepare a map of the names + positions inside blob at this level
	template<class... Ts> size_t deserialize_top(unsigned char *blob, char *list_of_args, Ts&...args) {

		//SRL_LOG("deserialize_top\n");

		// get names
		vector<string> names;
		MedSerialize::get_list_names(list_of_args, names);

		// get total size
		size_t tot_size = 0, pos = 0;
		pos += MedSerialize::deserialize(blob, tot_size);
		//SRL_LOG("dtop(1) pos %d tot_size %d\n", pos, tot_size);

		// prepare the map
		map<string, size_t> name2pos;
		while (pos < tot_size) {

			// read curr_size
			size_t curr_size;
			size_t orig_curr_pos = pos;
			//SRL_LOG("dtop(2) pos %d\n", pos);
			pos += MedSerialize::deserialize(blob + pos, curr_size);

			// read name
			string name;
			//SRL_LOG("dtop(3) pos %d curr_size %d\n", pos, curr_size);
			pos += MedSerialize::deserialize(blob + pos, name);

			// add to map
			//SRL_LOG("dtop(4) pos %d curr_size %d name %s\n", pos, curr_size, name.c_str());
			name2pos[name] = pos;

			// advance pos to next element
			pos = orig_curr_pos + curr_size;
			//cerr << "dtop(5) orig_pos " << orig_curr_pos << " pos " << pos << " curr_size " << curr_size << "\n";

		}

#if 0
		for (auto &e : name2pos)
			SRL_LOG("dtop(5.5) name2pos [ %s ]  =  %d\n", e.first.c_str(), e.second);
#endif
		// ready for deserializer (note that all name2pos addresses are relative to the original blob)
		MedSerialize::deserializer(blob, 0, names, name2pos, args...);

		// return size of deserialized area (taking it from tot_size (!) , and not from the actual deserialized part,
		// in case there were elements in the blob that were not asked for.)

		//cerr << "dtop(6) tot_size " << tot_size << "\n";
		return tot_size;
	}

	template<class T> string object_json_rec(int counter, vector<string> &names, T &elem) {
		stringstream str;
		//need to print only names[counter] elemet. the value is stored in elem

		str << "\n\t\"" << names[counter] << "\": ";
		str << MedSerialize::object_json(elem);

		return str.str();
	}

	template<class T, class... Ts> string object_json_rec(int counter, vector<string> &names, T &elem, Ts&...args) {
		stringstream str;

		//print current and get next:
		str << MedSerialize::object_json_rec(counter, names, elem);

		if (counter + 1 < names.size())
			str << ",";

		str << MedSerialize::object_json_rec(counter + 1, names, args...);

		return str.str();
	}

	template<class... Ts> string object_json_start(const string &cls_name, int vers, char *list_of_args, Ts&...args) {
		vector<string> names;
		MedSerialize::get_list_names(list_of_args, names);

		stringstream str;
		//str << "{\n\t\"Object\":\"" << cls_name << "\",";
		str << "{\n\t\"Object\":\"" << cls_name << "\",\n\t\"Version\":" << vers << ",";

		/*str << "{";
		for (int i = 0; i < names.size(); ++i) {
			if (i > 0)
				str << ",";
			str << "\n\t\"" << names[i] << "\":\"NOT_IMPLEMENTED\"";
		}
		str << "}";*/

		string res = MedSerialize::object_json_rec(0, names, args...);
		str << res << "\n}";

		return str.str();
	}

	//primitive types
	template<> inline string object_json<const int>(const int &v) { return to_string(v); }
	template<> inline string object_json<const unsigned int>(const unsigned int &v) { return to_string(v); }
	template<> inline string object_json<const float>(const float &v) { return to_string(v); }
	template<> inline string object_json<const double>(const double &v) { return to_string(v); }
	template<> inline string object_json<const char>(const char &v) { return to_string(v); }
	template<> inline string object_json<const unsigned char>(const unsigned char &v) { return to_string(v); }
	template<> inline string object_json<const long>(const long &v) { return to_string(v); }
	template<> inline string object_json<const long long>(const long long &v) { return to_string(v); }
	template<> inline string object_json<const unsigned long long>(const unsigned long long &v) { return to_string(v); }
	template<> inline string object_json<const bool>(const bool &v) { return to_string(v); }
	template<> inline string object_json<int>(int &v) { return to_string(v); }
	template<> inline string object_json<unsigned int>(unsigned int &v) { return to_string(v); }
	template<> inline string object_json<float>(float &v) { return to_string(v); }
	template<> inline string object_json<double>(double &v) { return to_string(v); }
	template<> inline string object_json<char>(char &v) { return to_string(v); }
	template<> inline string object_json<unsigned char>(unsigned char &v) { return to_string(v); }
	template<> inline string object_json<long>(long &v) { return to_string(v); }
	template<> inline string object_json<long long>(long long &v) { return to_string(v); }
	template<> inline string object_json<unsigned long long>(unsigned long long &v) { return to_string(v); }
	template<> inline string object_json<bool>(bool &v) { return to_string(v); }
	template<> inline string object_json<string>(string &str) { return "\"" + str + "\""; }
	template<> inline string object_json<const string>(const string &str) { return "\"" + str + "\""; }
	inline string object_json(int v) { return to_string(v); }
	inline string object_json(unsigned int v) { return to_string(v); }
	inline string object_json(float v) { return to_string(v); }
	inline string object_json(double v) { return to_string(v); }
	inline string object_json(char v) { return to_string(v); }
	inline string object_json(unsigned char v) { return to_string(v); }
	inline string object_json(long v) { return to_string(v); }
	inline string object_json(long long v) { return to_string(v); }
	inline string object_json(unsigned long long v) { return to_string(v); }
	inline string object_json(bool v) { return to_string(v); }


}

void mes_trim(string &s);

// A few IO helpers copied here in order to make SerializableObject a PURE h file implementation
namespace MedSerialize {

	inline int read_binary_data_alloc(const string &fname, unsigned char *&data, unsigned long long &size)
	{
		ifstream inf;

		inf.open(fname, ios::in | ios::binary | ios::ate);

		if (!inf) {
			SRL_ERR("read_binary_data_alloc(): can't open file %s for read\n", fname.c_str());
			return -1;
		}

		size = inf.tellg();
		data = new unsigned char[size];
		inf.seekg(0, ios::beg);
		inf.read((char *)data, size);

		boost::crc_32_type checksum_agent;
		checksum_agent.process_bytes(data, size);
		SRL_LOG("read_binary_data_alloc [%s] with crc32 [%d]\n", fname.c_str(), checksum_agent.checksum());

		inf.close();
		return 0;
	}

	//-----------------------------------------------------------------------------
	inline int write_binary_data(const string &fname, unsigned char *data, unsigned long long size)
	{
		ofstream of;

		//MLOG("Writing file %s :: size %lld\n", fname.c_str(), size);
		of.open(fname, ios::out | ios::binary);

		if (!of) {
			SRL_ERR("write_binary_data(): can't open file %s for write\n", fname.c_str());
			return -1;
		}

		of.write((char *)data, size);

		of.close();

		return 0;
	}

	// Initialization Utility
	static int initialization_text_to_map(const string& text, map<string, string>& init_map) {

		if (text == "")
			return 0;

		//MLOG("INPUT TEXT: %s\n", text.c_str());
		// dealing with {}
		// whenever there's a v={S} where S is any string (that may also include {}) we want the map to put S for v ...
		// follows is (hence) an ugly code to parse that
		// but this adds the ability to pass parameters for an embedded element within our current (say a model that holds parameters for other models)
		//

		vector<size_t> start_pos;
		vector<size_t> end_pos;
		vector<pair<size_t, size_t>> from_to;

		// find all positions of "={"
		size_t pos = text.find("={", 0);
		while (pos != string::npos) {
			start_pos.push_back(pos);
			pos = text.find("={", pos + 1);
		}

		// find all positions of "}"
		pos = text.find("}", 0);
		while (pos != string::npos) {
			end_pos.push_back(pos);
			pos = text.find("}", pos + 1);
		}

		// treating nesting 
		if (start_pos.size() > 0 && end_pos.size() > 0) {

			int i = 0, j = 0, stack = 0, stack_first = -1, stack_last = -1;

			while (j < end_pos.size()) {
				if (i < (int)start_pos.size() && start_pos[i] < end_pos[j]) {
					if (stack_first < 0) stack_first = (int)start_pos[i];
					stack++;
					i++;
				}
				else {
					if (stack == 0) {
						SRL_ERR("ERROR: Unmatched {} in string %s\n", text.c_str());
						return -1;
					}
					stack--;
					if (stack == 0) stack_last = (int)end_pos[j];
					j++;
				}

				if (stack == 0) {
					from_to.push_back(pair<size_t, size_t>(stack_first, stack_last));
					stack_first = -1;
				}
			}

			for (auto &ft : from_to) {
				SRL_LOG_D("found substring: %d-%d : %s\n", ft.first, ft.second, text.substr(ft.first, ft.second - ft.first + 1).c_str());
			}

		}


		// replacing {} areas with other strings to allow for correct parsing, and then returning them
		string new_text = "";
		map<string, string> replacers;
		if (from_to.size() == 0) new_text = text;
		else {
			new_text = text.substr(0, from_to[0].first + 1); // up to the first '='
			int j;
			for (j = 0; j < from_to.size(); j++) {
				string name = "REPLACE_ME_LATER_NUMBER_" + to_string(j);
				string replacer = text.substr(from_to[j].first + 2, from_to[j].second - from_to[j].first - 2);
				SRL_LOG_D("replacer %d : %s -> %s\n", j, name.c_str(), replacer.c_str());
				new_text += name;
				replacers[name] = replacer;
				if (j < from_to.size() - 1)
					new_text += text.substr(from_to[j].second + 1, from_to[j + 1].first - from_to[j].second);
			}
			new_text += text.substr(from_to[j - 1].second + 1, text.length() - from_to[j - 1].second);
			SRL_LOG_D("new_text is %s\n", new_text.c_str());

		}

		// TBD


		// get "Name = value" fields
		vector<string> fields;
		boost::split(fields, new_text, boost::is_any_of(";"));

		// get name + value
		vector<string> sub_fields;

		for (string& field : fields) {
			if (field.size() == 0)
				continue;

			boost::split(sub_fields, field, boost::is_any_of("="));
			if (sub_fields.size() != 2) {
				SRL_ERR("Cannot parse \'%s\' from \'%s\'\n", field.c_str(), text.c_str());
				return -1;
			}
			mes_trim(sub_fields[0]);
			mes_trim(sub_fields[1]);
			init_map[sub_fields[0]] = sub_fields[1];
		}

		for (auto &el : init_map) {
			if (el.second.compare(0, 24, "REPLACE_ME_LATER_NUMBER_") == 0) {
				init_map[el.first] = replacers[el.second];
			}
		}

		return 0;
	}
	//..............................................................................
	static int init_map_from_string(string text, map<string, string>& init_map) {

		if (text == "") return 0;

		// parse text of the format "Name = Value ; Name = Value ; ..."

		// remove white spaces
		text.erase(remove_if(text.begin(), text.end(), ::isspace), text.end());

		if (MedSerialize::initialization_text_to_map(text, init_map) == -1)
			return -1;

		//	for (auto rec : init_map)
		//		MLOG("Initializing with \'%s\' = \'%s\'\n", rec.first.c_str(), rec.second.c_str());


		return 0;
	}

	static int read_file_into_string(const string &fname, string &data)
	{
		ifstream inf(fname);
		if (!inf) {
			SRL_ERR("MedSerialize::read_file_into_string: Can't open file %s\n", fname.c_str());
			return -1;
		}

		data = "";
		string curr_line;
		while (getline(inf, curr_line)) {
			if ((curr_line.size() > 1) && (curr_line[0] != '#')) { // ignore empty lines, ignore comment lines

				// get rid of leading spaces, trailing spaced, and shrink inner spaces to a single one, get rid of tabs and end of line (win or linux)
				string fixed_spaces = boost::regex_replace(curr_line, boost::regex("^ +| +$|( ) +|\r|\n|\t+"), string("$1"));
				data += fixed_spaces;
			}
		}

		return 0;
	}


	// similar in concept to read_file_into_string, but gets the strings in the file and adds a comma "," between them
	static int read_list_into_string(const string &fname, string &data)
	{
		ifstream inf(fname);
		if (!inf) {
			SRL_ERR("MedSerialize::read_file_into_string: Can't open file %s\n", fname.c_str());
			return -1;
		}

		data = "";
		string curr_line;
		while (getline(inf, curr_line)) {
			//SRL_LOG("read_list: curr_line: %s\n", curr_line.c_str());
			if ((curr_line.size() > 1) && (curr_line[0] != '#')) { // ignore empty lines, ignore comment lines
				// move all tabs to spaces
				string fixed_spaces = boost::regex_replace(fixed_spaces, boost::regex("\t+"), " ");

				// get rid of leading spaced, ending spaces, \r
				fixed_spaces = boost::regex_replace(curr_line, boost::regex("^ +| +$|\r|\n"), "");

				fixed_spaces += "\n"; // re-adding eol, in case it was missing

				// change all internal spaces and \n to comma
				fixed_spaces = boost::regex_replace(fixed_spaces, boost::regex(" +|\n"), ",");

				// make sure there are no adjacent commas
				fixed_spaces = boost::regex_replace(fixed_spaces, boost::regex(",,+"), ",");

				// add 
				data += fixed_spaces;
				//SRL_LOG("read_list: data: %s\n", data.c_str());
			}
		}
		// could happen that last char is comma, we fix it.
		if (data.back() == ',') data.pop_back();

		return 0;
	}
}


#endif
