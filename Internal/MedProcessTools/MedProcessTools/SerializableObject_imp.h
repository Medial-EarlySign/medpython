//#ifndef __SERIALIZABLE_OBJECT_IMP_H__
//#define __SERIALIZABLE_OBJECT_IMP_H__

namespace MedSerialize {


	//.........................................................................................
	// templated simple ones : int, float, long, double, etc...
	//.........................................................................................
	template <class T> size_t get_size(T &elem)
	{
		return sizeof(T);
	}

	//.........................................................................................
	template <class T> size_t serialize(unsigned char *blob, T &elem)
	{
		memcpy(blob, &elem, sizeof(T));
		return sizeof(T);
	}

	//.........................................................................................
	template <class T> size_t deserialize(unsigned char *blob, T &elem)
	{
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
		memcpy(blob, &len, sizeof(size_t)); pos += sizeof(size_t);
		memcpy(blob+pos, str.c_str(), len); pos += len;
		blob[pos] = 0; pos++;
		return pos;
	}

	//.........................................................................................
	template<> inline size_t deserialize<string>(unsigned char *blob, string &str)
	{
		//fprintf(stderr,"string deserialize\n");
		size_t pos = 0;
		size_t len;
		memcpy(&len, blob, sizeof(size_t)); pos += sizeof(size_t);
		string new_s((char *)blob+pos);
		str = new_s;
		pos += len+1;
		return pos;
	}

	//.........................................................................................
	// vector of type T that has a MedSerialize function
	//.........................................................................................
	template<class T>  size_t get_size(vector<T> &v)
	{
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
		size_t pos = 0, len;
		pos += MedSerialize::deserialize<size_t>(blob + pos, len);
		v.clear();
		if (len > 0) {
			v.resize(len);
			for (T &elem : v)
				pos += MedSerialize::deserialize(blob+pos, elem);
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
	// Wrappers to call several elements is a single line
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