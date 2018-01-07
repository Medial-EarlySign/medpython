//

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



