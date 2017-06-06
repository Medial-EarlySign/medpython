//
// MedStat - Statistics utilities :
//			 2. MedPerformance used for analyzing performance
//

#ifndef _MED_PERFORMANCE_T_H_
#define _MED_PERFORMANCE_T_H_

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include "assert.h"
#include "math.h"

#include <vector>
#include <map>
#include <string>
#include <algorithm>

#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedFeat/MedFeat/MedFeat.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"

#include "string.h"

#define MED_CLEANER_MAX_Z 15
#define MED_CLEANER_EPSILON 0.0001

using namespace std ;

// Performance class : Calculation and Comparison of performance parameters.
enum PerformanceCompareMode {
	PRF_COMPARE_FULL = 1,
	PRF_COMPARE_ALL = 2,
	PRF_COMPARE_SPLITS = 3,
	PRF_COMPARE_FULL_AND_PART_SPLITS = 4,
	PRF_COMPARE_LAST = 5,
} ;


class Measurement {
public:
	string setParam ;
	string queriedParam ;
	float setValue ;

	Measurement(const string& sParam, const string& qParam, float sValue) ;
	Measurement(const string& qParam, float sValue) ;
	Measurement(const string& qParam) ;
	Measurement() ;

	inline bool operator<(const Measurement& otherMeasuemrent) const {
		return ((queriedParam < otherMeasuemrent.queriedParam) ||
			    (queriedParam == otherMeasuemrent.queriedParam && setValue < otherMeasuemrent.setValue) || 
				(queriedParam == otherMeasuemrent.queriedParam && setValue == otherMeasuemrent.setValue && setParam < otherMeasuemrent.setParam)) ;

	}
} ;


class MedClassifierPerformance {
public:

	// Prediction data (first elemnt = full data; next elements = splits)
	vector<vector<pair<float, float> > > preds ;

	// counters
	vector<int> npos,nneg ;
	vector<vector<int> > tps,fps ;
	// Performance
	vector<map<string, vector<float> > > PerformanceValues ;
	map<Measurement,vector<float> > MeasurementValues ;

	// Comparison
	Measurement comparePoint ;
	PerformanceCompareMode compareMode ;
	float partialCompareRatio ;

	// Helper - locations of queried PerformancePointer
	vector<map<pair<string,float> , pair<int,int> > > PerformancePointers ;
	map<string,int> compareDirections ;

	// Init
	void init() ;

	// Load
	void _load(vector<pair<float,float> >& in_preds) ;
	void _load(vector<vector<pair<float,float> > >& in_split_preds) ;
	void _load(MedFeaturesData& predictor_data) ;
	void _load(MedSamples& inSamples);
	void load_preds_on_train(MedFeaturesData& predictor_data);
	void post_load();
	template <typename T> void load(T& object) ;
	template <typename T, typename S> void load(T *preds, S *labels, int n) ;

	// Helper functions for loading
	void SplitsToComplete() ;
	void ShuffleSort() ;
	void Count() ;
	void getPerformanceValues() ;

	// Constructors
	MedClassifierPerformance()  ;
	template <typename T> MedClassifierPerformance(T& object) ;
	template <typename T, typename S> MedClassifierPerformance(T *preds, S *labels, int n) ;

	// Queries
	vector<float> operator() (Measurement& inMeasurement) ;
	// Parameter at point determined by another parameters (e.g. PPV at Specificity = 0.99 is GetPerformanceParam("PPV","Spec",0.99,outPPV). setParams = (Score,Sens,Spec), queriedParams = (Score,Sens,Spec,PPV,NPV,OR)
	int GetPerformanceParam(const string& setParam, const string& queriedParam, float setValue) ;
	int GetPerformanceParam(Measurement& inMeasurement) ;
	// General performance parameter, with optional value (e.g. AUC = GetPerformanceParam("AUC",outAuc) or GetPerformanceParam("AUC",1.0,outAUC). Partial AUC = GetPerformanceParam("AUC",0.2,partAUC)
	int GetPerformanceParam(const string& queriedParam, float setValue) ;
	int GetPerformanceParam(const string& queriedParam) ;
	// Performance Graph
	int GetPrformanceGraph(const string& xParam, const string& yParam, vector<vector<float> >& x, vector<vector<float> >& y) ;

	// Comparison
	int compare(MedClassifierPerformance& other) ;

	// Helpers for queries
	int getPerformancePointer(pair<string,float>& set, int index) ; 
	int getPointer(const string& param, float value, int index, int direction) ;
	void getAUC(float maxFPR, vector<float>& qValues) ;
	float getAUC(float maxFPR, int index) ;

	int getPerformanceValues(pair<string,float>& set, const string &queriedParam, int index, vector<float>& queriedValues) ;

private:
	struct _PredsCompare {
		bool operator()(const pair<float, float>& left, const pair<float, float>& right) {
			return (left.first > right.first) ;
		} ;
	} ;
} ;

// general useful routines
// Pearson Correlation // DO NOT REMOVE !!! WE USE THESE
float get_pearson_corr(float *v1, float *v2, int len);
float get_pearson_corr(vector<float> &v1, vector<float> &v2); // {return get_pearson_corr(VEC_DATA(v1),VEC_DATA(v2),(int)v1.size());}
// Pearson Correaltion after removing missing values. return number of values left in n.
float get_pearson_corr(vector<float> &v1, vector<float> &v2, int &n, float missing_val) ;


// n x m contigency table expected and chi2 score
double get_chi2_n_x_m(vector<int> &cnts, int n, int m, vector<double> &exp);
double get_chi2_n_x_m(vector<int> &cnts, int n, int m);

// AUC
float get_preds_auc(vector<float> &preds, vector<float> &y);
float get_preds_auc_q(const vector<float> &preds, const vector<float> &y);

// returns 4 counts a,b,c,d for the top 'size' percent of the data. a,b,c,d are:
//              0    1
// >= size ::   b    a
//  < size ::   d    c
int get_preds_perf_cnts(vector<float> &preds, vector<float> &y, vector<float> &size, int direction, vector<vector<int>> &cnts);
int cnts_to_perf(vector<int> &cnt, float &sens, float &spec, float &ppv, float &rr);

// multi category helpers
int multicateg_get_max_pred(vector<float> &probs, int nsamples, int ncateg, vector<float> &max_pred);
int multicateg_get_avg_pred(vector<float> &probs, int nsamples, int ncateg, vector<float> &avg_pred);
int multicateg_get_error_rate(vector<float> &probs, vector<float> &y, int nsamples, int ncateg, float &err_rate, float &rms, float &avg_rms);

int get_quantized_breakdown(vector<float> &preds, vector<float> &y, vector<float> &bounds, MedMat<int> &counts);
void print_quantized_breakdown(MedMat<int> &cnt, vector<float> &bounds);


//void get_mean_std(float *v, int len, float &mean, float &std);
//void get_mean_std(vector<float> &v, float &mean, float &std);

#include "MedPerformance_imp.h"

#endif