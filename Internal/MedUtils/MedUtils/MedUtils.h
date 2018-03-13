//
// MedUtils: Collection of several needed utilities (Files IO, Matrix data type, Data structures, etc)
//

// This is a simple h file to include all library parts

#ifndef __MED_UTILS_H__
#define __MED_UTILS_H__

#include "MedGenUtils.h"
#include "MedIO.h"
#include "MedMat.h"
#include "MedMedical.h"
#include "MedDataStructures.h"
#include "Logger/Logger/Logger.h"
#include "assert.h"
#include "MedPlot.h"

enum MedCorrelationType {
	CORR_PEARSON,
	CORR_LAST
} ;

enum MedBinningType {
	BIN_EQUIDIST,
	BIN_EQUISIZE,
	BIN_LAST
} ;

// Correlation of vectors
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n) ;
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n, float missing_value) ;
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n, MedCorrelationType type) ;
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n, float missing_value, MedCorrelationType type) ;


// Pearson correlation of two vectors
template <class S> double get_corr_pearson(vector<S>& x, vector<S>& y);
double get_corr_pearson(float *v1, float *v2, int len); // A COPY OF THE TEMPLATED VERSION FOR C INTERFACING

template <class S> double get_squared_dist(vector<S>& x, vector<S>& y);
template <class S> double get_abs_avg_dist(vector<S>& x, vector<S>& y);
template <class S> double get_abs_relative_avg_dist(vector<S>& x, vector<S>& y);
template <class S> double get_vecs_accuracy(vector<S>& x, vector<S>& y, double epsilon);

// Discretization
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins) ;
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, float missing_value) ;
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, MedBinningType binning) ;
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins,  int max_bins, float missing_value, MedBinningType binning) ;

// Counts from binned vectors
void get_counts(vector<int>& x, map<int,int>& counts, int& n) ;
int get_co_counts(vector<int>& x, vector<int>& y, map<pair<int,int>,int>& counts, int& n) ;

// Mutual Information
double get_mutual_information(map<int,int>& x_count, int nx, map<int,int>& y_count, int ny,  map<pair<int,int>,int>& xy_count, int n) ;
int get_mutual_information(vector<int>& x, vector<int>& y, double& mi, int &n) ;
double get_mutual_information(map<int,int>& x_count, int nx, map<int,int>& y_count, int ny,  map<pair<int,int>,int>& xy_count, int n) ;

// Moments
//template <class S> int get_moments(vector<S>& values, vector<float>& wgts, S missing_value, S& mean, S& sd);
int get_moments(vector<float>& values, vector<float>& wgts, float missing_value, float& mean, float& sd, bool do_missing=true);
int get_moments(float *values,float *wgts, int n, float missing_value, float& mean, float& sd, bool do_missing = true);

namespace medial {
	namespace print {
		template<class T> string print_obj(T obj, string format);
		template<class T> void print_vec(const vector<T> &vec, const string &title, const string &format);
		template<class T> void print_hist_vec(const vector<T> &vec, const string &title, const string &format);
	}
	namespace process {
		template<class T> void prctils(const vector<T> &x, const vector<double> &prc, vector<T> &res);
		template<typename T> int binary_search_index(const T *begin, const T *end, T val);
		template<typename T> int binary_search_position(const T *begin, const T *end, T val, bool reversed);
	}
}

#include "MedUtils_imp.h"

#endif