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
int get_moments(vector<float>& values, vector<float>& wgts, float missing_value, float& mean, float&sd) {

	int n = 0;
	double s = 0;

	for (unsigned int i = 0; i < values.size(); i++) {
		if (values[i] != missing_value) {
			n++;
			s += values[i];
		}
	}
	if (n == 0) {
		mean = 0;
		sd = 1.0;
		return 0;
	}

	mean = (float)(s / n);

	s = 0.0;
	for (unsigned int i = 0; i < values.size(); i++) {
		if (values[i] != missing_value)
			s += (values[i]-mean)*(values[i]-mean);
	}

	if (n > 1)
		sd = (float)sqrt((s / n));
	else
		sd = (float) 1.0;

	if (sd == 0) {
		MWARN("get_moments for all-zeros vector, fixing SD from 0.0 to 1.0\n");
		sd = 1.0;
	}
	return n;
}
