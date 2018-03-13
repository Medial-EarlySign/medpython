#define _CRT_SECURE_NO_WARNINGS
#include "MedUtils.h"

#define LOCAL_SECTION LOG_MED_UTILS
#define LOCAL_LEVEL	LOG_DEF_LEVEL


// Pearson Correlation: A pure C function
double  get_corr_pearson(float *v1, float *v2, int len)
{
	double sx,sy,sxy,sxx,syy,n;

	sx = sy = sxy = sxx = syy = 0;

	for (int i=0; i<len; i++) {
		sx += v1[i];
		sy += v2[i];
		sxx += v1[i]*v1[i];
		syy += v2[i]*v2[i];
		sxy += v1[i]*v2[i];
	}

	n = (double)len;

	sx /= n;
	sy /= n;
	sxx /= n;
	syy /= n;
	sxy /= n;

	double c1 = sxy - sx*sy;
	double c2 = sxx - sx*sx;
	double c3 = syy - sy*sy;

	double epsilon = 1e-8 ;
	if (c2 < epsilon || c3 < epsilon)
		return 0;

	return (float)(c1/(sqrt(c2)*sqrt(c3)));
}

// Mutual information of binned-vectors
int get_mutual_information(vector<int>& x, vector<int>& y, double& mi, int &n) {

	if (x.size() != y.size()) {
		MERR("Size mismatch. Quitting\n") ;
		return -1 ;
	}

	map<int,int> x_counts,y_counts ;
	map<pair<int,int>,int> xy_counts ;
	n = 0 ;

	for (unsigned int i=0; i<x.size(); i++) {
		if (x[i]!=-1 && y[i]!=-1) {
			x_counts[x[i]]++ ;
			y_counts[y[i]]++ ;
			xy_counts[pair<int,int>(x[i],y[i])]++ ;
			n++ ;
		}
	}

	if (n < 2) {
		MLOG_V("Not enough common non-missing entries for mutual information.\n") ;
		mi = -1 ;
	} else {
		mi = get_mutual_information(x_counts,n,y_counts,n,xy_counts,n) ;
	}

	return 0 ;
}

// Counts from binned vectors
void get_counts(vector<int>& x, map<int,int>& counts, int& n) {

	n=0 ;
	counts.clear() ;
	for (unsigned int i=0; i<x.size(); i++) {
		if (x[i] != -1) {
			counts[x[i]] ++ ;
			n++ ;
		}
	}
}

int get_co_counts(vector<int>& x, vector<int>& y, map<pair<int,int>,int>& counts, int& n) {

	if (x.size() != y.size()) {
		MERR("Size mismatch (%d vs %d). Quitting\n",x.size(),y.size()) ;
		return -1 ;
	}

	n=0 ;
	counts.clear() ;
	for (unsigned int i=0; i<x.size(); i++) {
		if (x[i] != -1 && y[i] != -1) {
			counts[pair<int,int>(x[i],y[i])] ++ ;
			n++ ;
		}
	}

	return 0 ;
}

// Mutual information from counts
double get_mutual_information(map<int,int>& x_count, int nx, map<int,int>& y_count, int ny,  map<pair<int,int>,int>& xy_count, int n) {

	double mi = 0 ;
	for (auto it = xy_count.begin(); it != xy_count.end(); it++) {
		double p = (it->second + 0.0)/n ;
		double px = (x_count[it->first.first] + 0.0)/nx ;
		double py = (y_count[it->first.second] + 0.0)/ny ;

		mi += p * log(p/px/py)/log(2.0) ;
	}

	return mi ;
}

// Get moments of a vector
int get_moments(vector<float>& values, vector<float>& wgts, float missing_value, float& mean, float&sd, bool do_missing) {

	return get_moments(&(values[0]), &(wgts[0]), (int)values.size(), missing_value, mean, sd, do_missing);
}

int get_moments(float *values, float* wgts, int size, float missing_value, float& mean, float&sd, bool do_missing) {
	double  n = 0;
	double s = 0;

	for (int i = 0; i < size; i++) {
		if (!do_missing || values[i] != missing_value) {
			n += wgts[i];
			s += wgts[i] * values[i];
		}
	}
	if (n == 0) {
		mean = 0;
		sd = 1.0;
		return 0;
	}

	mean = (float)(s / n);

	s = 0.0;
	for (int i = 0; i < size; i++) {
		if (!do_missing || values[i] != missing_value)
			s += wgts[i] * (values[i]-mean)*(values[i]-mean);
	}

	if (n > 1)
		sd = (float)sqrt((s / n));
	else
		sd = (float) 1.0;

	if (sd == 0) {
		MWARN("get_moments for all-zeros vector, fixing SD from 0.0 to 1.0\n");
		sd = 1.0;
	}
	return (int)n;
}

template<class T> string medial::print::print_obj(T obj, string format) {
	//return to_string((round(num * 1000) / 1000));
	char res[50];
	snprintf(res, sizeof(res), format.c_str(), obj);
	return string(res);
}
template string medial::print::print_obj<float>(float obj, string format);
template string medial::print::print_obj<double>(double obj, string format);
template<class T> void medial::process::prctils(const vector<T> &x, const vector<double> &prc, vector<T> &res) {
	if (x.size() == 0)
		throw invalid_argument("x is empty");

	vector<T> cp(x);
	T *data = cp.data();
	sort(cp.begin(), cp.end());

	res.resize((int)prc.size());
	for (size_t i = 0; i < res.size(); ++i)
	{
		double pos = prc[i] * (x.size() - 1);
		int pos_a = (int)pos;
		T r = data[pos_a];
		res[i] = r;
		if (pos_a + 1 < x.size())
			res[i] = T(res[i] * (1 - (pos - pos_a)) + (pos - pos_a) * data[pos_a + 1]);
	}
}
template void medial::process::prctils<float>(const vector<float> &x, const vector<double> &prc, vector<float> &res);
template void medial::process::prctils<double>(const vector<double> &x, const vector<double> &prc, vector<double> &res);

template<class T> void medial::print::print_vec(const vector<T> &vec, const string &title, const string &format) {
	if (vec.empty()) {
		MLOG("%s: EMPTY\n", title.c_str());
		return;
	}
	string bf = print_obj(vec[0], format);
	for (size_t i = 1; i < vec.size(); ++i)
		bf += ", " + print_obj(vec[i], format);

	MLOG("%s: [%s]\n", title.c_str(), bf.c_str());
}
template void medial::print::print_vec<double>(const vector<double> &vec, const string &title, const string &format);
template void medial::print::print_vec<float>(const vector<float> &vec, const string &title, const string &format);

template<class T> void medial::print::print_hist_vec(const vector<T> &vec, const string &title, const string &format) {
	if (vec.size() <= 10)
		print_vec(vec, title, format);
	vector<T> prcs;
	map<T, int> uniq_vals;
	for (size_t i = 0; i < vec.size(); ++i)
		++uniq_vals[vec[i]];

	if (uniq_vals.size() > 10) {
		medial::process::prctils(vec, { 0, 0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9 ,1 }, prcs);
		string bf = print_obj(prcs[0], format);
		for (size_t i = 1; i < prcs.size(); ++i)
			bf += ", " + print_obj(prcs[i], format);
		MLOG("%s: HISTOGRAM[%s]\n", title.c_str(), bf.c_str());
	}
	else {
		auto ii = uniq_vals.begin();
		string bf = print_obj(ii->first, format) + ":" +
			print_obj(100 * ii->second / double(vec.size()), "%2.2f");
		++ii;
		for (; ii != uniq_vals.end(); ++ii)
			bf += ", " + print_obj(ii->first, format) + ":" +
			print_obj(100 * ii->second / double(vec.size()), "%2.2f");
		MLOG("%s: VALUES[%s]\n", title.c_str(), bf.c_str());
	}
}
template void medial::print::print_hist_vec<double>(const vector<double> &vec, const string &title, const string &format);
template void medial::print::print_hist_vec<float>(const vector<float> &vec, const string &title, const string &format);

template<typename T> int medial::process::binary_search_index(const T *begin, const T *end, T val) {
	int maxSize = (int)(end - begin) + 1;
	int mid = int((maxSize - 1) / 2);
	if (maxSize <= 2) {
		if (*begin == val) {
			return 0;
		}
		else if (*end == val) {
			return 1;
		}
		else {
			return -1;
		}
	}

	if (begin[mid] == val) {
		//return first Index if there are more then one
		if (mid > 0 && begin[mid - 1] == val)
		{
			if (*begin == val) {
				return 0;
			}
			else {
				return binary_search_index(begin, begin + mid + 1, val);
			}
		}
		else {
			return mid;
		}
	}
	else {
		if (val > begin[mid]) {
			int p = binary_search_index(begin + mid, end, val);
			if (p > 0) {
				return mid + p;
			}
			else {
				return -1;
			}
		}
		else {
			return binary_search_index(begin, begin + mid, val);
		}
	}
}
template<typename T> int medial::process::binary_search_position(const T *begin, const T *end, T val, bool reversed) {
	int maxSize = (int)(end - begin) + 1;
	int mid = int((maxSize - 1) / 2);
	if (maxSize <= 2) {
		if (!reversed) {
			if (val < *begin) {
				return 0;
			}
			else if (val < *end) {
				return 1;
			}
			else {
				return maxSize;
			}
		}
		else {
			if (val >= *begin) {
				return 0;
			}
			else if (val >= *end) {
				return 1;
			}
			else {
				return maxSize;
			}
		}
	}

	if (!reversed) {
		if (val <= begin[mid]) {
			return binary_search_position(begin, begin + mid, val, reversed);
		}
		else {
			return mid + binary_search_position(begin + mid, end, val, reversed);
		}
	}
	else {
		if (val >= begin[mid]) {
			return binary_search_position(begin, begin + mid, val, reversed);
		}
		else {
			return mid + binary_search_position(begin + mid, end, val, reversed);
		}
	}
}
template int medial::process::binary_search_position<int>(const int *begin, const int *end, int val, bool reversed);
template int medial::process::binary_search_position<double>(const double *begin, const double *end, double val, bool reversed);
template int medial::process::binary_search_position<float>(const float *begin, const float *end, float val, bool reversed);

template int medial::process::binary_search_index(const float *begin, const float *end, float val);
template int medial::process::binary_search_index(const int *begin, const int *end, int val);
template int medial::process::binary_search_index(const double *begin, const double *end, double val);
template int medial::process::binary_search_index(const string *begin, const string *end, string val);