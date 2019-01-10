// MedSparseVec
// ------------
//
// A useful data structure to hold a vector of elements of type T when the vector is sparse.
// If the vector V is in the range [a,b] (b-a+1 length) we would like to be able to insert, and query 
//
// We would like our memory to be as low as possible. 
// This implementation will add 1.5bits for each possible range.
// So if the range (always an unsigned int) of the min,max values of keys for V is [a,b], and N = b-a+1 ,
// and there are actually only M<N T elements, we will use memory of size: 1.5N/8 + M*sizeof(T) bytes.
//
// since we want to save memory this is better than a simple array only if:
//   N*sizeof(T) > 1.5N/8 + M*sizeof(T)
// or:
// M < N - 1.5N/(8*sizeof(T)) =(for sizeof(T)=4) N - 1.5N/32 = 0.953N , 
//
// so...in practice even if we keep just an int of memory as an item, we start saving memory even if the table is 0.95 full !!!
//
// other options are to use map<int,T> but map is using lots of memory per key, something like 32 bytes....
// hence map will be better only for much sparser cases:
//
// M*(sizeof(T)+32) < 1.5N/8 + M*sizeof(T)
// or:
// M < 1.5N/(8*32) = 3N/512 = 0.005N !! 
//
// Another option to compare to is keep a vector of the the pair <uint, T> (maybe even sorted)
// this costs M*(sizeof(T)+4) and is VERY slow in lookup
// however even here we get
// M*(sizeof(T)+4) < 1.5N/8 + M*sizeof(T)
// or:
// M < 1.5N/32 = 0.046N .... so in sparsness of the 5%-95% we will be better even than just KEEPING the data AS IS !!!!
// 
// so for sparsness in the range of 0.5% to 95% this data structure will save memory.
//
// Another drawback on this data structure: 
// - in order to get a fast insert time, it is best to:
//   (1) set the range [a,b] in advance.
//   (2) insert the elements with the keys SORTED, such that if we enter V[i] after V[j] then j > i (or appeared before).
//       failing to do so will return an error (or a terrible insertion time).
// 
// so the actual situation to start with when using this package is :
// (1) have a vector of size M of keys K[i] , which are unsigned int and sorted.
// (2) have a vector of size M of elements E[i] of type T, such that the key for E[i] is K[i].
// 
// OR:
// (1) make sure that whenever a new pair {k,e} is inserted then k is >= of all the keys that were before or appeared already (in which case we replace its 
// ... value with the new e.
//
// 

#ifndef __MED__SPARSE__VEC__
#define __MED__SPARSE__VEC__

#include <vector>
#include <immintrin.h>
#include <nmmintrin.h>
#include <cstring>

using namespace std;

#define MED_SPARSE_VEC_MAGIC_NUM 0x0102030405060708

template <class T> class MedSparseVec {

  public:

	unsigned int min_val;
	unsigned int max_val;
	int max_set;
	unsigned int max_key;
	T	def_val;
	
	vector<unsigned int> counts;
	vector<unsigned long long> is_in_bit;
	vector<T> data;

	void set_min(int _min) { min_val= _min; }
	void set_max(int _max) { max_val= _max; max_set = 1; }
	MedSparseVec() { min_val = 0; max_val = 0; max_set = 0; max_key = 0; init();  }
	MedSparseVec(int _min) { min_val = _min; max_val = 0; max_set = 0; max_key = 0; init();  }
	MedSparseVec(int _min, int _max) { min_val = _min; max_set = 1; max_val = _max; max_key = 0; init(); }
	void set_def(const T val) { def_val = val; data[0] = def_val; }

	inline T operator[] (const unsigned int key) const { return (*get(key)); }

	inline T &operator[] (const unsigned int key) { return (*get(key)); }
	
	void reserve(unsigned int size) {data.reserve(size);}
	
	//------------------------------------------------------
	// init major arrays
	//------------------------------------------------------
	void init() {
		counts.resize(2, 0);
		is_in_bit.resize(1, 0);
		data.clear();
//		def_val = (T)0;
		def_val = 0;
		data.push_back(def_val);
	}


	//------------------------------------------------------
	// clear a used sparse vec, going back to all defaults
	void clear() {
		counts.assign(2, 0);
		is_in_bit.assign(1, 0ULL);
		data.clear();
		//def_val = (T)0;
		def_val = 0;
		data.push_back(def_val);
		min_val = 0; max_val = 0; max_set = 0; max_key = 0;
	}

	//------------------------------------------------------
	// get index of a data item - 0 is not found
	//------------------------------------------------------
	inline unsigned int get_ind(const unsigned int key)
	{

		if (key < min_val || key > max_key || (max_set && key>max_val))
			return 0; // not in range

		unsigned int mkey = key - min_val;
		unsigned int i_count = mkey >> 6;
		unsigned int i_bit = 63 - (mkey & 0x3f);
		unsigned long long kbits = is_in_bit[i_count]>>i_bit;
		//unsigned long long i_mask = (((unsigned long long)1)<<(i_bit));

		//fprintf(stderr, "key=%d mkey=%d i_count=%d %d i_bit=%d kbits=%d ind= %d\n", key, mkey, i_count, counts[i_count], i_bit, kbits, -1); fflush(stderr);
		//if ((is_in_bit[i_count]&i_mask) == 0)
		if ((kbits & 0x1) == 0)
			return 0; // no value inserted

//		unsigned int ind = counts[i_count] + (int)_mm_popcnt_u64(is_in_bit[i_count]>>i_bit);
		unsigned int ind = counts[i_count] + (int)_mm_popcnt_u64(kbits);

		//fprintf(stderr, "key=%d mkey=%d i_count=%d %d i_bit=%d i_mask=%llx %llx pos = %d val = %d\n", key, mkey, i_count, counts[i_count], i_bit, i_mask, is_in_bit[i_count]>>i_bit, pos, data[pos]); fflush(stderr);

		return ind;
	}

	//------------------------------------------------------
	// insert a new element
	//------------------------------------------------------
	int insert(const unsigned int key, const T elem)
	{
		// first we check we are allowed to insert it
		if (key<min_val || (max_set && key>max_val)) return -1;

		// get the bit indexes of key
		unsigned int mkey = key - min_val;
		unsigned int i_count = mkey >> 6;
		unsigned int i_bit = 63 - (mkey & 0x3f);
		unsigned long long i_mask = (((unsigned long long)1)<<(i_bit));

//		fprintf(stderr, "key=%d mkey=%d i_count=%d i_bit=%d i_mask=%llx max_key=%d is_in_bit=%llx\n", key, mkey, i_count, i_bit, i_mask,max_key, is_in_bit[i_count]);
		if (key <= max_key) {
			if (is_in_bit[i_count]&i_mask) {
				unsigned int pos = counts[i_count] + (int)_mm_popcnt_u64(is_in_bit[i_count]>>i_bit);
				data[pos] = elem;
				return 0;
			}
			else
				if (key < max_key)
					return -2; // elements were not inserted in the correct sorted order
				// in key == max_key with bit 0 we simply have to insert the value
		}

		// new max key
		unsigned int last_cnt = 0;
		if (counts.size() > 0)
			last_cnt = counts.back();
		counts.resize(i_count+2, last_cnt);
		is_in_bit.resize(i_count+1, 0);
		is_in_bit[i_count] |= i_mask;
		counts[i_count+1] = counts[i_count] + (int)_mm_popcnt_u64(is_in_bit[i_count]);
		data.push_back(elem);
		max_key = key;
		return 0;
	}

	//---------------------------------------------------------------
	// get - NULL is returned for a key not inside.
	//----------------------------------------------------------------
	inline T *get(unsigned int key) {
		return &data[get_ind(key)];
	}


	//---------------------------------------------------------------
	// get_all_keys
	//---------------------------------------------------------------
	int get_all_keys(vector<unsigned int> &keys) {

		keys.resize(data.size()-1);
		unsigned int j=0;

		//unsigned int base = min_val;
		//for (int i=0; i<counts.size(); i++) {
		//	unsigned long long bits = is_in_bit[i];
		//	if (bits != 0)
		//		for (int k=63; k>=0; k--) {
		//			unsigned long long kbits = bits>>k;
		//			if (kbits & 0x1) {
		//				keys[j++] = base + 63 - k;
		//			}
		//		}
		//	base += 64;
		//}


		for (unsigned int i=min_val; i<=max_key; i++)
			if (get_ind(i) > 0)
				keys[j++] = i;

		return 0;
	}

	//---------------------------------------------------------------
	// get_all_intersected key:
	// gets a uniq list of in_keys
	// outputs:
	// keys - the keys in the list that are also in in_keys
	// inds - the indexes for these keys (in a vector of the same size)
	//---------------------------------------------------------------
	int get_all_intersected_keys(const vector<int> &in_keys, vector<int> &keys, vector<int> &inds) {
		
		keys.resize(in_keys.size());
		inds.resize(in_keys.size());
		int i_size = 0;

		for (int i=0; i<in_keys.size(); i++) {
			int ind = get_ind(in_keys[i]);
			if (ind > 0) {
				keys[i_size] = in_keys[i];
				inds[i_size] = ind;
				i_size++;
			}
		}
/*
		for (int i=0; i<in_keys.size(); i++) {
			unsigned int curr_key = in_keys[i];
			//fprintf(stderr, "working on curr_key %d ind %d min_val %d max_key %d\n", curr_key, get_ind(curr_key), min_val, max_key);
			if (curr_key >= min_val && curr_key <= max_key) {
				unsigned int mkey = curr_key - min_val;
				int i_count = mkey>>6;
				int k = 63 - (mkey & 0x3f);
				unsigned long long bits = is_in_bit[i_count];
				if (bits) {
					unsigned long long kbits = bits >> k;
					//fprintf(stderr, "mkey %d i_count %d k %d bits %d\n", mkey, i_count, k, bits);
					if (kbits &0x1) {
						keys[i_size] = curr_key;
						inds[i_size] = counts[i_count] + (int)_mm_popcnt_u64(kbits);
						i_size++;
					}
				}
			}
		}
*/
		keys.resize(i_size);
		inds.resize(i_size);

		//fprintf(stderr, "i_size = %d/%d %d %d \n", i_size, in_keys.size(), keys.size(), inds.size()); fflush(stderr);
		return 0;
	}
	//---------------------------------------------------------------
	// Serializations
	//---------------------------------------------------------------
	size_t get_size() {
		size_t size = 0;

		size += sizeof(unsigned long long); // len of serialization
		size += sizeof(unsigned long long); // magic number recognizer
		size += sizeof(unsigned int); // min_val
		size += sizeof(unsigned int); // max_val
		size += sizeof(int);			// max_set
		size += sizeof(unsigned int); // max_key
		size += sizeof(T);			// def_val
		size += sizeof(unsigned int); // len counts
		size += sizeof(unsigned int) * counts.size(); // counts vector
		size += sizeof(unsigned int); // len is_in_bit
		size += sizeof(unsigned long long) * is_in_bit.size(); // is_in_bit vector
		size += sizeof(unsigned int);	// len data
		size += sizeof(T) * data.size();

		return size;

	}

	//---------------------------------------------------------------
	size_t serialize(unsigned char *blob) {

		unsigned char *curr = blob;

		curr += sizeof(unsigned long long); // bypassing len - will place it at the end.
		((unsigned long long *)curr)[0] = (unsigned long long)MED_SPARSE_VEC_MAGIC_NUM; curr+= sizeof(unsigned long long);
		((unsigned int *)curr)[0] = min_val; curr+= sizeof(unsigned int);
		((unsigned int *)curr)[0] = max_val; curr+= sizeof(unsigned int);
		((int *)curr)[0] = max_set; curr+= sizeof(int);
		((unsigned int *)curr)[0] = max_key; curr+= sizeof(unsigned int);
		((T *)curr)[0] = def_val; curr+= sizeof(T);
		((unsigned int *)curr)[0] = (unsigned int)counts.size(); curr+= sizeof(unsigned int);
		for (int i=0; i<counts.size(); i++) {
			((unsigned int *)curr)[0] = counts[i]; curr+= sizeof(unsigned int);
		}
		((unsigned int *)curr)[0] = (unsigned int)is_in_bit.size(); curr+= sizeof(unsigned int);
		for (int i=0; i<is_in_bit.size(); i++) {
			((unsigned long long *)curr)[0] = is_in_bit[i]; curr+= sizeof(unsigned long long);
		}
		((unsigned int *)curr)[0] = (unsigned int)data.size(); curr+= sizeof(unsigned int);
		for (int i=0; i<data.size(); i++) {
			((T *)curr)[0] = data[i]; curr+= sizeof(T);
		}
		unsigned long long len = (unsigned long long)(curr-blob);
		((unsigned long long *)blob)[0] = len;
		return len;

	}
	
	//---------------------------------------------------------------
	size_t deserialize(unsigned char *blob) {

		unsigned char *curr = blob;

		counts.clear();
		is_in_bit.clear();
		data.clear();
		unsigned long long serialize_len = ((unsigned long long *)curr)[0]; curr += sizeof(unsigned long long);
		unsigned long long magic_num = ((unsigned long long *)curr)[0]; curr += sizeof(unsigned long long);
		if (magic_num != (unsigned long long)MED_SPARSE_VEC_MAGIC_NUM) {
			fprintf(stderr, "ERROR: Sparse Vec Magic Num wrong : can't deserialize() (%llx)\n", magic_num);
			return 0;
		}
		min_val = ((unsigned int *)curr)[0]; curr += sizeof(unsigned int);
		max_val = ((unsigned int *)curr)[0]; curr += sizeof(unsigned int);
		max_set = ((int *)curr)[0]; curr += sizeof(int);
		max_key = ((unsigned int *)curr)[0]; curr += sizeof(unsigned int);
		def_val = ((T *)curr)[0]; curr += sizeof(T);
		unsigned int len_counts = ((unsigned int *)curr)[0]; curr += sizeof(unsigned int);
		counts.resize(len_counts);
		memcpy(&counts[0], curr, len_counts*sizeof(unsigned int));
		curr += sizeof(unsigned int) * len_counts;

		unsigned int len_is_in_bit = ((unsigned int *)curr)[0]; curr += sizeof(unsigned int);
		is_in_bit.resize(len_is_in_bit);
		memcpy(&is_in_bit[0], curr, len_is_in_bit*sizeof(unsigned long long));
		curr += sizeof(unsigned long long)*len_is_in_bit;

		unsigned int len_data = ((unsigned int *)curr)[0]; curr += sizeof(unsigned int);
		data.resize(len_data);
		memcpy(&data[0], curr, len_data*sizeof(T));
		curr += sizeof(T) * len_data;
		
		size_t len = curr - blob;
		if (len != serialize_len) {
			fprintf(stderr, "ERROR: Sparse Vec serialize len not matching decalred one: %zu != %llu\n", len, serialize_len);
			return 0;
		}
		return len;
	}

};

#endif
