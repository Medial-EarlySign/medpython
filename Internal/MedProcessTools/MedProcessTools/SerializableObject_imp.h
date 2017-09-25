//#ifndef __SERIALIZABLE_OBJECT_IMP_H__
//#define __SERIALIZABLE_OBJECT_IMP_H__

namespace MedSerialize {

	//========================================================================================
	// FUNCTION FROM BELOW THAT IMPLEMENT AN STD CONTAINER (map, vector, etc) 
	// MUST BE PRE DECLARED HERE BEFORE THE IMPLEMENTATION
	// TO MAKE SURE THE CORRECT RECURSIVE RESOLVING IS DONE
	// (simple, string & variadic need not - and will create a compilation error)
	//========================================================================================

	template<class T>  size_t get_size(vector<T> &v);
	template <class T> size_t serialize(unsigned char *blob, vector<T> &v);
	template <class T> size_t deserialize(unsigned char *blob, vector<T> &v);
	template <class T, class S> size_t get_size(map<T, S> &v);
	template <class T, class S> size_t serialize(unsigned char *blob, map<T, S> &v);
	template <class T, class S> size_t deserialize(unsigned char *blob, map<T, S> &v);
	template <class T, class S> size_t get_size(unordered_map<T, S> &v);
	template <class T, class S> size_t serialize(unsigned char *blob, unordered_map<T, S> &v);
	template <class T, class S> size_t deserialize(unsigned char *blob, unordered_map<T, S> &v);
	template <class T> size_t get_size(unordered_set<T> &v);
	template <class T> size_t serialize(unsigned char *blob, unordered_set<T> &v);
	template <class T> size_t deserialize(unsigned char *blob, unordered_set<T> &v);
	template <class T> size_t get_size(set<T> &v);
	template <class T> size_t serialize(unsigned char *blob, set<T> &v);
	template <class T> size_t deserialize(unsigned char *blob, set<T> &v);

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
		//cout << "inside simple serialize with type " << typeid(T).name() << endl;
		memcpy(blob, &elem, sizeof(T));
		return sizeof(T);
	}

	//.........................................................................................
	template <class T> size_t deserialize(unsigned char *blob, T &elem)
	{
		//cout << "inside simple deserialize with type " << typeid(T).name() << endl;
		memcpy(&elem, blob, sizeof(T));
		return sizeof(T);
	}
	
	//.........................................................................................
	// string is special
	//.........................................................................................
	template<> inline size_t get_size<string>(string &str) { 
		size_t size = 0;
		size += sizeof(size_t); // length of string
		size += str.length()+1;
		return size;
	}

	//.........................................................................................
	template<> inline size_t serialize<string>(unsigned char *blob, string &str)
	{
		size_t pos = 0;
		size_t len = str.length();
		//fprintf(stderr, "string serialize(%d) %s\n", len, str.c_str());
		memcpy(blob, &len, sizeof(size_t)); pos += sizeof(size_t);
		memcpy(blob+pos, str.c_str(), len); pos += len;
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
		pos += len+1;
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
				pos += MedSerialize::serialize(blob+pos, elem);
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
				pos += MedSerialize::deserialize(blob+pos, elem);
		}
		return pos;
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


	//.........................................................................................
	// map<T,S> both with a MedSerialize function
	//.........................................................................................
	template <class T, class S> size_t get_size(map<T,S> &v)
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
		int i = 0;
		if (len > 0)
			for (auto &elem : v) {
				T *t = (T *)&elem.first;
				S *s = (S *)&elem.second;
				pos += MedSerialize::serialize(blob+pos, (*t));
				pos += MedSerialize::serialize(blob+pos, (*s));
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
			for (int i=0; i<len; i++) {
				pos += MedSerialize::deserialize(blob+pos, t);
				pos += MedSerialize::deserialize(blob+pos, s);
				v[t] = s;
			}
		}
		return pos;
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
			for (int i = 0; i<len; i++) {
				pos += MedSerialize::deserialize(blob + pos, t);
				pos += MedSerialize::deserialize(blob + pos, s);
				v[t] = s;
			}
		}
		return pos;
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

}


//#endif