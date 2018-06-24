// Cross Validator = an object for performing cross validation = learning set filterer + model + test set filterer

#ifndef _MED_CV_H_
#define _MED_CV_H_

#include <InfraMed/InfraMed/InfraMed.h>
#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/RepProcess.h>
#include <MedProcessTools/MedProcessTools/FeatureProcess.h>
#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedProcessTools/MedProcessTools/SampleFilter.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>

#define DEFAULT_CV_NTHREADS 8

// TBD : splits object should be integrated here

//.......................................................................................
//.......................................................................................
// An object for performing cross-validation
//.......................................................................................
//.......................................................................................

class CrossValidator {
public:

	// Threading
	int nthreads;

	// Learning set filter(s)
	vector <SampleFilter *> learning_set_filters;

	// Test set filter(s)
	vector< SampleFilter *> test_set_filters;

	// Model 
	MedModel *model;

	// Constructor/Destructor
	CrossValidator() { nthreads = DEFAULT_CV_NTHREADS; };
	CrossValidator(MedModel *_model) { nthreads = DEFAULT_CV_NTHREADS; model = _model; };
	~CrossValidator() {};

	// Add
	void add_learning_set_filter(SampleFilter *filter) { learning_set_filters.push_back(filter); };
	void add_test_set_filter(SampleFilter *filter) { test_set_filters.push_back(filter); };

	// Run
	int doCV(MedPidRepository& rep, MedSamples& samples, int nfolds, MedSamples& outSamples);
	int doCV(MedPidRepository& rep, MedSamples& samples, int nfolds, map<int, int>& folds, MedSamples& outSamples);
};

#endif