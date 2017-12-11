#ifndef _SAMPLE_FILTER_H_
#define _SAMPLE_FILTER_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedValueCleaner.h"

#define DEFAULT_SMPL_FLTR_NTHREADS 8

//.......................................................................................
//.......................................................................................
// SampleFilter is the parent class for creating a list of samples to work on from a 
//		another list and/or a repository
// Basic functionalities:
//		learn : learn the filtering parameters from a given MedSamples + optional
//				MedRepository (this method can be empty)
//		filter : generate a new MedSample from a given MedSamples + optional MedRepository
//				This method must be implemented for each inheriting class
//		Variants for learn and filter :
//				- no repository required
//				- in-place filtering
//.......................................................................................
//.......................................................................................
// Define types of sample filter
typedef enum {
	SMPL_FILTER_TRN,
	SMPL_FILTER_TST,
	SMPL_FILTER_OUTLIERS,
	SMPL_FILTER_MATCH,
	SMPL_FILTER_REQ_SIGNAL,
	SMPL_FILTER_BASIC,
	SMPL_FILTER_LAST
} SampleFilterTypes;

class SampleFilter : public SerializableObject {
public:

	// Threading 
	int nthreads;

	// Type
	SampleFilterTypes filter_type;

	// Constructor/Destructor
	SampleFilter() { nthreads = DEFAULT_SMPL_FLTR_NTHREADS; };
	~SampleFilter() {};

	// Init from name or type. optinally with a parameters string
	static SampleFilter *make_filter(string name);
	static SampleFilter *make_filter(SampleFilterTypes type);
	static SampleFilter *make_filter(string name, string params);
	static SampleFilter *make_filter(SampleFilterTypes type, string params);

	// initialize : from object/string/defaults.
	// Should be implemented for inheriting classes that have parameters
	virtual int init(void *params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	// Learn with Repository
	// _learn should be implemented for inheriting classes that learn parameters
	virtual int _learn(MedRepository& rep, MedSamples& samples) { return _learn(samples); }
	int learn(MedRepository& rep, MedSamples& samples) { return _learn(rep, samples); }

	// Learn without Repository
	// _learn should be implemented for inheriting classes that learn parameters
	virtual int _learn(MedSamples& samples) { return 0; }
	int learn(MedSamples& samples) { return _learn(samples); }

	// Filter with Repository
	// _filter must be implemented for inheriting classes that use repository for filtering
	virtual int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {return _filter(inSamples,outSamples) ;}
	int filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) { return _filter(rep, inSamples, outSamples); }
	int filter(MedRepository& rep, MedSamples& samples);

	// Filter without Repository
	// _filter must be implemented for all inheriting classes
	virtual int _filter(MedSamples& inSamples, MedSamples& outSamples) { return 0; }
	int filter(MedSamples& inSamples, MedSamples& outSamples) { return _filter(inSamples, outSamples); }
	int filter(MedSamples& samples);

	// Serialization (including type)
	size_t get_filter_size();
	size_t filter_serialize(unsigned char *blob);
};

// Utilities
SampleFilterTypes sample_filter_name_to_type(const string& filter_name);

//.......................................................................................
//.......................................................................................
// Training set filter
//	- take all controls samples (outcome=0) and all cases before outcomeTime
//.......................................................................................
//.......................................................................................
class BasicTrainFilter : public SampleFilter {
public:

	BasicTrainFilter() : SampleFilter() { filter_type = SMPL_FILTER_TRN; };
	~BasicTrainFilter() {};

	// Filter
	int _filter(MedSamples& inSamples, MedSamples& outSamples);
};

//.......................................................................................
//.......................................................................................
// Test set filter
//	- dummy filter - take everything
//.......................................................................................
//.......................................................................................
class BasicTestFilter : public SampleFilter {
public:

	BasicTestFilter() : SampleFilter() { filter_type = SMPL_FILTER_TST; };
	~BasicTestFilter() {};

	// Filter
	int _filter(MedSamples& inSamples, MedSamples& outSamples);
};

//.......................................................................................
//.......................................................................................
// Outliers filter
//	- A filter that remove samples with outlier-outcomes (suitable for regression)
//	  Outliers detection is done using MedValueCleaner's methods (through inheritance)
//.......................................................................................
//.......................................................................................
#define SMPL_FLTR_TRIMMING_SD_NUM 7
#define SMPL_FLTR_REMOVING_SD_NUM 7

class OutlierSampleFilter : public SampleFilter, public MedValueCleaner {
public:

	OutlierSampleFilter() : SampleFilter() { init_defaults(); };
	~OutlierSampleFilter() {};

	// Filter
	int _filter(MedSamples& inSamples, MedSamples& outSamples);

	// Learning : check outlier-detection method and call appropriate learner (iterative/quantile)
	int _learn(MedSamples& samples);
	// Learning : learn outliers using MedValueCleaner's iterative approximation of moments
	int iterativeLearn(MedSamples& samples);
	// Learning : learn outliers using MedValueCleaner's quantile appeoximation of moments
	int quantileLearn(MedSamples& samples);
	// Helper for learning - extract all outcomes from samples.
	void get_values(MedSamples& samples, vector<float>& values);

	// Init from map
	int init(map<string, string>& mapper) { init_defaults();  return MedValueCleaner::init(mapper); }
	// Init defualts
	void init_defaults() {
		filter_type = SMPL_FILTER_OUTLIERS;
		params.trimming_sd_num = SMPL_FLTR_TRIMMING_SD_NUM; params.removing_sd_num = SMPL_FLTR_REMOVING_SD_NUM;
		params.take_log = 0;
		params.doTrim = false;
		params.doRemove = true;
		params.type = VAL_CLNR_ITERATIVE;
		params.missing_value = MED_MAT_MISSING_VALUE;
	};

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

//.......................................................................................
//.......................................................................................
// Matching filter 
//	- match cases and controls according to a set of matching strata
//  - matching creiterion can be
//		- value of signal within window
//		- age 
//		- sample time
//
//	- Each sample is assigned to a bin according to the vector of strata.
//	  The optimal case/control sampling ratio is then found, when each control that is
//	  removed costs 1 point and each case costs eventToControlPriceRatio points. The 
//	  maximal control to case ratio is set by matchMaxRatio. One ratio is decided, sampling
//	  is perfored randomly within each bin.
//.......................................................................................
//.......................................................................................

typedef enum {
	SMPL_MATCH_SIGNAL, // Match by value of signal
	SMPL_MATCH_AGE,	   // Match by age
	SMPL_MATCH_TIME,   // Match by time
	SMPL_MATCH_LAST
} SampleMatchingType;

class matchingParams {
public:

	// Match by ?
	SampleMatchingType match_type;

	// Matching details
	string signalName; // For matching by signal
	int timeWindow, windowTimeUnit ; // For matching by signal
	int matchingTimeUnit; // For matching by time
	float resolution ;

	// Helpers (for matching by signal)
	int signalId;
	bool isTimeDependent; // is the signal time-dependent (e.g. hemoglobin) or not (e.g. byear)
	int signalTimeUnit;

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

class MatchingSampleFilter : public SampleFilter {
public:

	vector<matchingParams> matchingStrata; // Matching parameters

	float eventToControlPriceRatio = 100.0; // Cost of removing case relative to removing control
	float matchMaxRatio = 10.0; // Maximal control/case ratio
	int verbose = 0; // control level of debug printing

	// helpers
	int samplesTimeUnit; // Time unit of samples
	int byearId, ageId; // signal-id for byear (if given) or age (if directly given)

	// Constructor/Destructor
	MatchingSampleFilter() { init_defaults(); };
	~MatchingSampleFilter() {};

	// Initialization
	int init(map<string, string>& mapper);
	void init_defaults() { filter_type = SMPL_FILTER_MATCH; };

	// Add a matching stratum defined by a string
	int addMatchingStrata(string& init_string);

	// Utilities
	// Check if repository is needed for matching (strata includes signal/age)
	bool isRepRequired();
	// Check if age is needed for matching
	bool isAgeRequired();
	// Indexing of a single sample according to strata
	int getSampleSignature(MedSample& sample, MedRepository& rep, string& signature);
	// add indexing of a single sample according to a single stratum to sample's index
	int addToSampleSignature(MedSample& sample, matchingParams& stratum, MedRepository& rep, string& signature);
	// initialize values of requried helpers
	int initHelpers(MedSamples& inSamples, MedRepository& rep);
	// Get all signals required  for matching
	int get_required_signals(vector<string> req_sigs);

	// Filter with repository
	int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
	// Filter without repository (return -1 if repository is required)
	int _filter(MedSamples& inSamples, MedSamples& outSamples);

	// Utility for Matching optimization : find optimal case/control ratio. 
	// the price of giving up 1 control is 1.0, the price of giving up 1 event is w 
	float get_pairing_ratio(map<string, pair<int, int>> cnts, float w);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};


//.......................................................................................
//.......................................................................................
// Required Signal Filter 
//	- Keep only samples with a required signal appearing in a time-window.
//	- OBSOLETE - REPLACED BY BasicSampleFilter. KEPT HERE FOR BACKWARD COMPETABILITY
//.......................................................................................
//.......................................................................................

class RequiredSignalFilter : public SampleFilter {
public:

	string signalName; // Required signal
	int timeWindow, windowTimeUnit; // Time before sample-time

									// Constructor/Destructor
	RequiredSignalFilter() { init_defaults(); };
	~RequiredSignalFilter() {};

	int init(map<string, string>& mapper);
	void init_defaults();

	// Filter
	int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
	int _filter(MedSamples& inSamples, MedSamples& outSamples);


	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};


//.......................................................................................
//.......................................................................................
// A general filter to allow the following basics:
// (1) min and max time of outcomeTime
// (2) option to allow for a signal to be in some given range (if there is a sample in some given window)
// (3) option to force N samples of a specific signal within the given range (in the given time window)
//.......................................................................................
//.......................................................................................

// Parameters for BasicFilter
struct BasicFilteringParams : public SerializableObject {
	// Name of signal to filter by
	string sig_name;
	// Time window for deciding on filtering
	int win_from = 0;
	int win_to = (int)(1<<30);
	// Allowed values range for signal
	float min_val = -1e10;
	float max_val = 1e10;
	// Required number of instances of signal within time window
	int min_Nvals = 1;
	// Signal parameters
	int time_channel = 0;
	int val_channel = 0;

	// Initialization from string
	int init_from_string(const string &init_str);
	// Test filtering criteria. Return 1 if passing and 0 otherwise
	int test_filter(MedSample &sample, MedRepository &rep, int win_time_unit);

	ADD_SERIALIZATION_FUNCS(sig_name, time_channel, val_channel, win_from, win_to, min_val, max_val, min_Nvals);

private:
	int sig_id = -1; // uninitialized until first usage
};				


class BasicSampleFilter : public SampleFilter {
public:

	// filtering parameters
	int min_sample_time = 0;
	int max_sample_time = (int)(1<<30); // these should always be given in the samples' time-unit
	vector<BasicFilteringParams> bfilters;
	int winsTimeUnit = MedTime::Days;

	// next is initialized with init string
	vector<string> req_sigs; // useful to load the repository needed for this filter

	// init from mapped string
	int init(map<string, string>& mapper);

	int get_req_signals(vector<string> &reqs);

	// Filter
	int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
	int _filter(MedSamples& inSamples, MedSamples& outSamples); // relevant only if bfilters is empty

	ADD_SERIALIZATION_FUNCS(min_sample_time, max_sample_time, bfilters, winsTimeUnit);
};

//
// Sanity Simple Filter helps making sanity tests on input data
//
// The basic tests optional are:
// (1) test that the signal actually exist in name (in the signals list in the repository)
// (2) within a given window: minimal number of tests
// (3) within a given window: maximal number of outliers
// (4) count outliers within a given window
// 
struct SanitySimpleFilter : public SerializableObject {
	// Name of signal to filter by
	string sig_name;
	// Time window for deciding on filtering
	int win_from = 0;
	int win_to = (int)(1<<30);
	// Allowed values range for signal
	float min_val = -1e10;
	float max_val = 1e10;
	// Required number of instances of signal within time window
	int min_Nvals = 0;
	int max_Nvals = -1; // -1 signs to not test this

	int max_outliers = -1; // -1 means don't do the max_outliers test

	int win_time_unit = MedTime::Days;

	// Signal parameters
	int time_channel = 0;
	int val_channel = 0;


	// Initialization from string
	int init_from_string(const string &init_str);
	// Test filtering criteria. Return 1 if passing and 0 otherwise
	int test_filter(MedSample &sample, MedRepository &rep) {
		int nvals, noutliers; return test_filter(sample, rep, nvals, noutliers);
	}
	int test_filter(MedSample &sample, MedRepository &rep, int &nvals, int &noutliers);

	// test_filter return codes
	const static int Passed = 0;
	const static int Failed = 1; // General fail due to reasons different than the following
	const static int Signal_Not_Valid = 2;
	const static int Failed_Min_Nvals = 3;
	const static int Failed_Max_Nvals = 4;
	const static int Failed_Outliers = 5;



	ADD_SERIALIZATION_FUNCS(sig_name, time_channel, val_channel, win_from, win_to, min_val, max_val, min_Nvals, max_Nvals, max_outliers, win_time_unit);

private:
	int sig_id = -1; // uninitialized until first usage
};

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(SanitySimpleFilter);
MEDSERIALIZE_SUPPORT(BasicFilteringParams);
MEDSERIALIZE_SUPPORT(BasicSampleFilter);

#endif
