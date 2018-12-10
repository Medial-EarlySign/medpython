// Implementation of template functions for MedPerformance.


// Init from array
template <typename T, typename S> MedClassifierPerformance::MedClassifierPerformance(T *preds, S *labels, int n) {
	load(preds,labels,n) ;
}

// Load from arrays
template <typename T, typename S> void MedClassifierPerformance::load(T *_preds, S *_labels, int n) {
	
	preds.resize(1) ;
	preds[0].resize(n) ;
	for (int i=0; i<n; i++) {
		pair<float,float> current((float) _preds[i],(float) _labels[i]) ;
		preds[0][i] = current ;
	}

	init() ;
	ShuffleSort() ;
	getPerformanceValues() ;

	PerformancePointers.resize(preds.size()) ;
}
	
// Init from object
template <typename T> MedClassifierPerformance::MedClassifierPerformance(T& inObj) {
	load(inObj) ;
}

// Load from object
template <typename T> void MedClassifierPerformance::load(T& inObj) {

	_load(inObj) ;	
	post_load();
}

//iterative algorithm for calculating Kendall Tau
//v1,v2 should be two equal length vectors, which may be of any type that can be cast to double (float, int, etc.)
//is01Vec1, is01Vec2 - are v1, v2 0/1 vectors, respectively. It is the user's responsibility that this is correct
//returns the KendallTau vector, or < -1.0 on bad input
//in general, ties are considered as not contributing to the score, while the tied pairs are counted.
//linear-time algorithm used if both vectors are 0/1. Otherwise, the algorithm is nlog(n), where n is the v1/2.size().
//If not both are 0/1 vectors, random noise is used to break ties.
template <typename T, typename S>
double kendallTau(const vector<T>& _v1, const vector<S>& _v2, bool is01Vec1 = false, bool is01Vec2 = false) {

	if (_v1.size() != _v2.size())
		return -2.0;

	if (_v1.size() <= 1)
		return -3.0;

	vector<double> v1(_v1.begin(), _v1.end());
	vector<double> v2(_v2.begin(), _v2.end());

	int len = (int)v1.size();

	long long nPairs = (long long)len * (long long)(len - 1) / 2;

	//case 1 - both vectors are 0/1. in this case, (n11*n00 - n10*n01) / nPairs

	if (is01Vec1 && is01Vec2) {
		long long n10 = 0, n01 = 0, n11 = 0;

		for (int i = 0; i < len; ++i) {
			if (v1[i] == 1) {
				if (v2[i] == 1)
					n11++;
				else
					n10++;
			}
			else {//v1[i] == 0
				if (v2[i] == 1)
					n01++;
			}
		}

		long long n00 = (long long)len - n11 - n10 - n01;

		return (double)(n11*n00 - n01*n10) / (double)nPairs;
	}

	//random mechanism used in resolving ties
	default_random_engine gen;
	unsigned int seed = 13;
	gen.seed(seed);
	uniform_real_distribution<double> dist(-0.25, 0.25); //can only change order for equal integers, not different ones

	long long nBadPairs = 0; //pairs that violate order. 

							 //case 2 - one vector is 0/1. WLOG it is v2 (otherwise we swap)
	if (is01Vec1 && !is01Vec2) {
		swap(v1, v2);
		swap(is01Vec1, is01Vec2);
	}

	if (!is01Vec1 && is01Vec2) {
		//ties are resolved fairly. we add [-0.25,0.25] random noise, after replacing the values with (tied) ranks		
		map<double, int> m;
		for (auto x : v1)
			m[x] = 1;//dummy value

		int i = 0;
		for (auto& x : m)
			x.second = i++;//translate value to rank. Thus distinct values are separated by at least 1

		for (auto& x : v1)
			x = (double)(m[x]) + dist(gen);

		vector<pair<double, double> > p; p.reserve(len);

		for (int i = 0; i < len; ++i)
			p.push_back(pair<double, double>(v1[i], v2[i]));

		//sort by the first vector
		sort(p.begin(), p.end());

		//the number of reversals is the rank sum of the n0 0-elements in v[2] minus n0*(n0-1)/2
		int n0 = 0;
		for (int i = 0; i < len; ++i) {
			if (p[i].second == 0) {
				nBadPairs += i;
				n0++;
			}
		}

		nBadPairs -= ((long long)n0*(long long)(n0 - 1) / 2);

		return (double)((long long)n0*(long long)(len - n0) - 2 * nBadPairs) / (double)nPairs;
	}

	//at this point, neither v1 nor v2 is 0/1

	//ties in both v1 and v2 are resolved fairly. we add [-0.25,0.25] random noise, after replacing the values with (tied) ranks

	map<double, int> m;
	for (auto x : v1)
		m[x] = 1;//dummy value

	for (auto x : v2)
		m[x] = 1;

	int i = 0;
	for (auto& x : m)
		x.second = i++;//translate value to rank

	for (auto& x : v1)
		x = (double)(m[x]) + dist(gen);

	for (auto& x : v2)
		x = (double)(m[x]) + dist(gen);

	//fast algorithm: 
	//basic idea: if we pair the values of the two vectors, sort the pairs by one coordinate, then replace its values with 0,2,..len-1, 
	//and then sort by the second, then we are reduced to counting reversals of the indices (integers) in the first coordinate. 
	//this is done recursively as follows:
	//we start with a range of size len of indices from 0 to len-1. Locate the subset with the values above len/2 (S1). Denote the other subset by S0. 
	//Reversals are either "inter" (one element in each set) or "intra" (two elements inside the same set). We count the "inter" reversals, and then
	//move S0 to the first half of the range and S1 to the second half without changing their relative internal order.
	//We subtract a constant from the values in each subrange (half) so that the values are again 0 - subrange length, and we can call the function recursively.
	//the recursive calls will get the "intra" reversals, and together we get the total reversals. 
	//the heart if the calculation is realizing the the "inter" count is a simple function of the rank sum of S0 elements inside the range; if their positions are 
	//i_0, ..., i_k-1, then the inter count is sum_{j=0}^{k-1} (i_j - j).
	//The recursion is simulated by adding each new subrange to a "to do" stack, instead of a recursive call stack. We continue until all subranges are processed.
	//Note that counting "intra" reversals in a range is independent of all other ranges, so the order of processing the subranges is immaterial.

	vector<pair<double, double> > p; p.reserve(len);

	for (int i = 0; i < len; ++i)
		p.push_back(pair<double, double>(v1[i], v2[i]));

	//sort by the first vector
	sort(p.begin(), p.end(), ComparePairByFirst<double, double>());

	//translate it to rank
	for (int i = 0; i < len; ++i)
		p[i].first = i;

	//sort by the second vector. We will then need only count reversals in the first vector
	//which has values [0,len)
	sort(p.begin(), p.end(), ComparePairBySecond<double, double>());

	//we now want to count the number of reversals in the first coordinate
	vector<int> res; res.reserve(len);

	for (int i = 0; i < len; ++i)
		res.push_back((int)(p[i].first));

	//the algorithm requires two vectors as workspaces
	vector<int> ws[2];
	ws[0].reserve(len);
	ws[1].reserve(len);

	//use a stack to emulate recursion calls iteratively
	//each pair is <start point, length> for a range inside res
	//which requires processing
	vector<pair<int, int> > stack; stack.reserve(len);
	stack.push_back(pair<int, int>(0, len)); //initially, the whole range requires processing

											 //while there are sub-ranges of res to process
	while (!stack.empty()) {
		//	for (int i = 0; i < stack.size(); ++i)
		//		cout << "(" << stack[i].first << ", " << stack[i].second << ") ";
		//	cout << endl;

		//pop the last element
		pair<int, int> inds = stack.back();
		int n = inds.second;
		stack.pop_back();
		ws[0].clear(); ws[1].clear();

		//we now handle the range indicated by inds
		if (n <= 1)
			continue;

		//the lowest n/2 elements and the highest n - n/2 elements are separated
		//the violations *between* the two sets are a function of the rank sum of the former set

		//go over elements in the range (in the res vector, from indices inds.first to inds.first + inds.second-1)
		for (int i = 0; i < n; ++i) {
			int ind = (int)inds.first + i;
			if (res[ind] >= n / 2)
				ws[1].push_back(res[ind]);
			else {
				nBadPairs += i;
				ws[0].push_back(res[ind]);
			}
		}
		nBadPairs -= ((long long)(ws[0].size()) * (long long)(ws[0].size() - 1) / 2);

		//copy back to res
		for (int i = 0; i < n / 2; ++i)
			res[inds.first + i] = ws[0][i];

		stack.push_back(pair<int, int>(inds.first, n / 2));

		for (int i = n / 2; i < n; ++i)
			res[inds.first + i] = ws[1][i - n / 2] - n / 2;

		stack.push_back(pair<int, int>(inds.first + n / 2, n - n / 2));
	}

	return 1.0 - 2.0 * (double)nBadPairs / (double)nPairs;
}


