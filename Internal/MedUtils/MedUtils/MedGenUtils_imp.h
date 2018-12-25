//
#ifndef __MED_GEN_UTILS_IMP_H__
#define __MED_GEN_UTILS_IMP_H__

//................................................................................................
template <class T> void get_zero_inds(T *v, int len, vector<int> &inds)
{
	inds.clear();
	for (int i=0; i<len; i++)
		if (v[i] == (T)0)
			inds.push_back(i);
}

//................................................................................................
template <class T> void get_nonzero_inds(T *v, int len, vector<int> &inds)
{
	inds.clear();
	for (int i=0; i<len; i++)
		if (v[i] != (T)0)
			inds.push_back(i);
}

//................................................................................................
template <class T> void get_vec_from_vecvec(vector<vector<T>> &v_in, vector<T> &v_out)
{
	v_out.clear();

	for (int i=0; i<v_in.size(); i++)
		for (int j=0; j<v_in[i].size(); j++)
			v_out.push_back(v_in[i][j]);
}

//................................................................................................
// gets number of different values in a vector
template <class T> int get_vec_ndiff_vals(vector<T> &v)
{

	if (v.size() == 0) return 0;

	map<T, int> m;

	for (int i=0; i<v.size(); i++) {
		if (m.find(v[i]) == m.end())
			m[v[i]] = 1;
	}

	return ((int)m.size());
}



// given two sorted vectors, where <in> has unique entries, find where the values of <search> would fit in <in>
// if a searched entry is found exactily in <in> it will be considered to fit after the equal entry in <in>
template <typename T> int find_sorted_vec_in_sorted_vec(const vector<T> &search, const vector<T> &in, vector<size_t>& indices) {
	indices.clear();

	if (in.empty())
		return -1;

	if (search.empty())
		return 0;

	indices.resize(search.size());

	//at this point search.size() >= 1
	
	//search_start, search_end - range of indices in <search> to handle
	int search_start = (int)(search.size()); 
	int search_end = 0;

	//handle first search entries that fit at the beginning of <in>
	for (int j = 0; j < search.size(); ++j) {
		if (search[j] < in[0]) {
			indices[j] = 0;
		}
		else {
			search_start = j;
			break;
		}
	}

	for (int j = (int)(search.size()) - 1; j >= search_start; --j) {
		if (search[j] >= in.back()) {
			indices[j] = in.size();
		}
		else {
			search_end = j;
			break;
		}
	}

	//note that for search_start =< j <= search_end we have search[j] in [in[0],in.back())
	
	if (search_start > search_end)
		return 0;	

	int start = 0;

	for (int j = search_start; j <= search_end; ++j) {
		T e = search[j];

		int end = (int)(in.size()) - 1;

		//we will always have an invariant: e is in [in[start],in[end])

		while (1) {			
			if (start >= end - 1) {
				//cout << indices.size() << " " << j << " " << start << " " << end << " " << search_start  << " " << search_end << endl;
				indices[j] = (size_t)end;
				break;
			}

			int mid = (start + end) / 2;

			if (e >= in[mid]) {
				start = mid;				
			}
			else {
				end = mid;
			}
		}
	}

	return 0;
}


//................................................................................................
// generates an arithmetic sequence
template<typename T> int sequence(T start, T finish, T step, vector<T>& seq, bool isForward) {

	if (step <= 0)
		return -1;

	seq.clear();
	seq.reserve(1 + (size_t)((finish - start) / step));

	if (isForward) {
		T cur = start;

		while (cur <= finish) {
			seq.push_back(cur);
			cur += step;
		}
	}
	else {//backward
		T cur = finish;

		while (cur >= start) {
			seq.push_back(cur);
			cur -= step;
		}

		reverse(seq.begin(), seq.end());
	}

	seq.shrink_to_fit();
	return 0;
}

//................................................................................................
//comparators for pairs - by first and by second elements

template<typename S, typename T> struct ComparePairBySecond {
	bool operator()(const pair<S, T>& left, const pair<S, T>& right) {
		return (left.second < right.second);
	};
};

template<typename S, typename T> struct ComparePairByFirst {
	bool operator()(const pair<S, T>& left, const pair<S, T>& right) {
		return (left.first < right.first);
	};
};

#endif