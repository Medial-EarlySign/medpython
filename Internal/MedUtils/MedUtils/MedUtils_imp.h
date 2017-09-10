//
// templated code for MedUtils class, included after class definition
//

#include <math.h>

// Correlation of vectors. Missing value = MED_MAT_MISSING_VAL.
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n) {
	return get_corr(_x,_y, corr, n, (float) MED_MAT_MISSING_VALUE,CORR_PEARSON) ;
}

template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, float missing_value, int& n) {
	return get_corr(_x,_y,corr,n,missing_value,CORR_PEARSON) ;
}

template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n, MedCorrelationType type) {
	return get_corr(_x,_y,corr,n,(float) MED_MAT_MISSING_VALUE,type) ;
}

template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n, float missing_value, MedCorrelationType type) {
		
	if (type >= CORR_LAST) {
		MEDLOG(LOG_MED_UTILS,MAX_LOG_LEVEL,"Unknown correlation type %d\n",type) ;
		return -1 ;
	}

	if (_x.size() != _y.size()) {
		MEDLOG(LOG_MED_UTILS,MAX_LOG_LEVEL,"Size mismatch. Quitting\n") ;
		return -1 ;
	}

	vector<double> x,y ;
	for (unsigned int i=0; i<_x.size(); i++) {
		if (_x[i] != missing_value && _y[i] != missing_value) {
			x.push_back((double)_x[i]) ;
			y.push_back((double)_y[i]) ;
		}
	}

	n = (int) x.size() ;
	if (n < 2) {
		MEDLOG(LOG_MED_UTILS,VERBOSE_LOG_LEVEL,"Not enough common non-missing entries for correlation.\n") ;
		corr = -2.0 ;
	} else {
		if (type == CORR_PEARSON) 
			corr = get_corr_pearson(x,y) ;
	}

	return 0 ;
}

// Find Pearson correlation of two vectors. NOTE: A pure C version is in MedUtils.cpp
template <class S> double get_corr_pearson(vector<S>& x, vector<S>& y) {

	double sx = 0 ;
	double sy = 0 ;
	int n = (int) x.size() ;

	for (unsigned int i=0; i<x.size(); i++) {
		sx += x[i] ;
		sy += y[i] ;
	}

	double meanx = sx/n ;
	double meany = sy/n ;

	double sxx = 0 ;
	double sxy = 0 ;
	double syy = 0 ;

	for (unsigned int j=0; j<x.size(); j++) {
		sxx += (x[j]-meanx)*(x[j]-meanx) ;
		syy += (y[j]-meany)*(y[j]-meany) ;
		sxy += (x[j]-meanx)*(y[j]-meany) ;
	}

	double epsilon = 1e-8;

	if (sxx<epsilon) sxx=1 ;
	if (syy<epsilon) syy=1 ;


	double r2 = sxy/sqrt(sxx*syy) ;
	return r2 ;
}

// Find Pearson correlation of two vectors. NOTE: A pure C version is in MedUtils.cpp
template <class S> double get_squared_dist(vector<S>& x, vector<S>& y) {

	double sxy2 = 0;
	int n = (int)x.size();

	for (unsigned int i=0; i<x.size(); i++)
		sxy2 += (double)(x[i] - y[i])*(x[i] - y[i]);

	if (n > 0)
		return sxy2/(double)n;

	return 0;
}

template <class S> double get_abs_avg_dist(vector<S>& x, vector<S>& y) {
	double dabs = 0;
	int n = (int)x.size();

	for (unsigned int i=0; i<x.size(); i++)
		dabs += abs((double)(x[i]-y[i]));

	if (n > 0)
		return dabs/(double)n;

	return 0;
}

template <class S> double get_abs_relative_avg_dist(vector<S>& x, vector<S>& y) {
	double dabs_rel = 0;
	int n = (int)x.size();

	for (unsigned int i=0; i<x.size(); i++)
		if (x[i] != 0) {
			dabs_rel += abs((double)(x[i]-y[i]))/abs((double)x[i]);
		}

	if (n > 0)
		return dabs_rel/(double)n;

	return 0;
}

template <class S> double get_vecs_accuracy(vector<S>& x, vector<S>& y, double epsilon) {
	double n_similar = 0;
	int n = (int)x.size();

	for (unsigned int i=0; i<x.size(); i++)
		if (abs((double)(x[i]-y[i])) <= epsilon)
			n_similar++;

	if (n > 0)
		return n_similar/(double)n;

	return 0;
}


// Discretization of values
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins) {
	return discretize(x,binned_x,nbins,max_bins,MED_MAT_MISSING_VALUE,BIN_EQUISIZE) ;
}

template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, float missing_value) {
	return discretize(x,binned_x,nbins,max_bins,missing_value,BIN_EQUISIZE) ;
}

template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, MedBinningType binning) {
	return discretize(x,binned_x,nbins,max_bins,MED_MAT_MISSING_VALUE,binning) ;
}

template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins,  int max_bins, float missing_value, MedBinningType binning) {

	binned_x.clear();
	binned_x.reserve(x.size());

	if (binning >= BIN_LAST) {
		MEDLOG(LOG_MED_UTILS,MAX_LOG_LEVEL,"Unknown binning type %d\n",binning) ;
		return -1 ;
	}

	map<float,int> x_values ;
	for (unsigned int i=0; i<x.size(); i++) {
		if (x[i] != missing_value)
			x_values[(float)x[i]]++ ;
	}

	nbins = (int) x_values.size() ;
	map<float,int> x_index ;

	if (nbins <= max_bins) { // Leave as is
		int idx = 0 ;
		for (auto it=x_values.begin(); it!=x_values.end(); it++)
			x_index[it->first] = idx++ ;
		assert(idx == nbins) ;
	} else { // Need to combine values into bins
		if (binning == BIN_EQUIDIST) {
			float min_val = x_values.begin()->first ;
			float max_val = x_values.rbegin()->first ;
			float bin_size = (max_val - min_val)/max_bins ;

			for (auto it=x_values.begin(); it!=x_values.end(); it++)
				x_index[it->first] = (int) ((it->first - min_val)/bin_size) ;
		} else if (binning == BIN_EQUISIZE) {
			int tot = 0 ;
			for (auto it=x_values.begin(); it!=x_values.end(); it++)
				tot += it->second ;
			int bin_size = tot/max_bins ;

			tot = 0 ;
			for (auto it=x_values.begin(); it!=x_values.end(); it++) {
				x_index[it->first] = tot/bin_size ;
				tot += it->second ;
			}
		}

		nbins = max_bins ;
	}


	for (unsigned int i=0; i<x.size(); i++) {
		if (x[i] == missing_value)
			binned_x.push_back(-1) ;
		else
			binned_x.push_back(x_index[x[i]]) ;
	}

	return 0 ;
}

//// Get moments of a vector
//int get_moments(vector<float>& values, vector<float>& wgts, float missing_value, float& mean, float&sd) {
//
//	int n = 0;
//	double s = 0;
//
//	for (unsigned int i = 0; i < values.size(); i++) {
//		if (values[i] != missing_value) {
//			n++;
//			s += values[i];
//		}
//	}
//
//	mean = (float)(s / n);
//
//	s = 0.0;
//	for (unsigned int i = 0; i < values.size(); i++) {
//		if (values[i] != missing_value)
//			s += (values[i]-mean)*(values[i]-mean);
//	}
//
//	if (n > 1)
//		sd = (float)sqrt((s / n));
//	else
//		sd = (float) 1.0;
//
//	return n;
//}
/*
template <class S> int get_moments(vector<S>& values, vector<float>& wgts, S missing_value, S& mean, S& sd) {

	int n = 0;
	double s = 0;

	for (unsigned int i = 0; i < values.size(); i++) {
		if (values[i] != missing_value) {
			n++;
			s += values[i];
		}
	}

	mean = (S)(s / n);

	s = 0.0;
	for (unsigned int i = 0; i < values.size(); i++) {
		if (values[i] != missing_value) 
			s += (values[i]-mean)*(values[i]-mean);
	}

	if (n > 1)
		sd = (S) sqrt((s / n));
	else
		sd = (S) 1.0;

	return n;
}
*/

/*
template <class S> int get_moments(vector<S>& values, vector<float>& wgts, S missing_value, S& mean, S& sd) {

	int n = 0;
	double s = 0, s2 = 0;

	for (unsigned int i = 0; i < values.size(); i++) {
		if (values[i] != missing_value) {
			n++;
			s += values[i];
			s2 += values[i] * values[i];
		}
	}

	mean = (S) (s / n);
	if (n > 1)
		sd = (S) sqrt((s2 - mean*s / n) / (n - 1));
	else
		sd = (S) 1.0;

	return n;
}
*/	