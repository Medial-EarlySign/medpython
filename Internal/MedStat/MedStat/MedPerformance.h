//
// MedPerformance used for analyzing performance
//

#ifndef _MED_PERFORMANCE_T_H_
#define _MED_PERFORMANCE_T_H_

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include <assert.h>
#include <math.h>

#include <vector>
#include <map>
#include <string>
#include <algorithm>

#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <MedUtils/MedUtils/MedGenUtils.h>
#include <string.h>

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

	Measurement(const string& qParam, const string& sParam, float sValue) ; // Look for value of qParams when sParam=sValue
	Measurement(const string& qParam, float sValue) ; // Look at value of qParam parameterized by sValue (e.g. "auc",0.2 for partial AUC)
	Measurement(const string& qParam) ; // Look at value of qParam (e.g. "auc")
	Measurement() ;

	inline bool operator<(const Measurement& otherMeasuemrent) const {
		return ((queriedParam < otherMeasuemrent.queriedParam) ||
			    (queriedParam == otherMeasuemrent.queriedParam && setValue < otherMeasuemrent.setValue) || 
				(queriedParam == otherMeasuemrent.queriedParam && setValue == otherMeasuemrent.setValue && setParam < otherMeasuemrent.setParam)) ;

	}

	static void get_from_name(string& fileName, Measurement& msrs);
	static void read_from_file(string& fileName, vector<Measurement>& msrs);
	string name();
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
	void _load(MedSamples& inSamples);
	void _load(MedFeatures& ftrs);
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

#include "MedPerformance_imp.h"

#endif