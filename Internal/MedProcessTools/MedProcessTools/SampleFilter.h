#ifndef _SAMPLE_FILTER_H_
#define _SAMPLE_FILTER_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedValueCleaner.h"

//.......................................................................................
//.......................................................................................
// Sample Filter creates a list of samples to work on from a repository or another list
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

#define DEFAULT_SMPL_FLTR_NTHREADS 8

class SampleFilter : public SerializableObject {
public:

	// Threading
	int nthreads;

	// Type
	SampleFilterTypes filter_type;

	// Constructor/Destructor
	SampleFilter() { nthreads = DEFAULT_SMPL_FLTR_NTHREADS; };
	~SampleFilter() {};

	// Init
	static SampleFilter *make_filter(string name);
	static SampleFilter *make_filter(SampleFilterTypes type);
	static SampleFilter *make_filter(string name, string params);
	static SampleFilter *make_filter(SampleFilterTypes type, string params);

	virtual int init(void *params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	// Learn
	// With Repository
	virtual int _learn(MedRepository& rep, MedSamples& samples) { return _learn(samples); }
	int learn(MedRepository& rep, MedSamples& samples) { return _learn(rep, samples); }

	// Without Repository
	virtual int _learn(MedSamples& samples) { return 0; }
	int learn(MedSamples& samples) { return _learn(samples); }

	// Filter
	// With Repository
	virtual int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {return _filter(inSamples,outSamples) ;}
	int filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) { return _filter(rep, inSamples, outSamples); }
	int filter(MedRepository& rep, MedSamples& samples);

	// Without Repository
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
// A filter that remove samples with outlier-outcomes (suitable for regression)
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

	// Learn
	int _learn(MedSamples& samples);
	int iterativeLearn(MedSamples& samples);
	int quantileLearn(MedSamples& samples);
	void get_values(MedSamples& samples, vector<float>& values);

	// Init
	int init(map<string, string>& mapper) { init_defaults();  return MedValueCleaner::init(mapper); }
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
// A matching filter
//.......................................................................................
//.......................................................................................

typedef enum {
	SMPL_MATCH_SIGNAL,
	SMPL_MATCH_AGE,
	SMPL_MATCH_TIME,
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
	bool isTimeDependent;
	int signalTimeUnit;

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

class MatchingSampleFilter : public SampleFilter {
public:

	vector<matchingParams> matchingStrata;

	float eventToCasePriceRatio = 100.0;
	float matchMaxRatio = 10.0;

	// helpers
	int samplesTimeUnit;
	int byearId, ageId;

	// Constructor/Destructor
	MatchingSampleFilter() { init_defaults(); };
	~MatchingSampleFilter() {};

	int init(map<string, string>& mapper);
	void init_defaults() { filter_type = SMPL_FILTER_MATCH; };
	int addMatchingStrata(string& init_string);

	// Utilities
	bool isRepRequired();
	bool isAgeRequired();
	int getSampleSignature(MedSample& sample, MedRepository& rep, string& signature);
	int addToSampleSignature(MedSample& sample, matchingParams& stratum, MedRepository& rep, string& signature);
	int initHelpers(MedSamples& inSamples, MedRepository& rep);
	int get_required_signals(vector<string> req_sigs);

	// Filter
	int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
	int _filter(MedSamples& inSamples, MedSamples& outSamples);

	// Matching optimization
	float get_pairing_ratio(map<string, pair<int, int>> cnts, float w);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};


//.......................................................................................
//.......................................................................................
// A filter on required signal
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

struct BasicFilteringParams : public SerializableObject {
	string sig_name;
	int win_from = 0;
	int win_to = (int)(1<<30);
	float min_val = -1e10;
	float max_val = 1e10;
	int min_Nvals = 1;
	int time_channel = 0;
	int val_channel = 0;

	int init_from_string(const string &init_str);
	int test_filter(MedSample &sample, MedRepository &rep, int win_time_unit);

	ADD_SERIALIZATION_FUNCS(sig_name, time_channel, val_channel, win_from, win_to, min_val, max_val, min_Nvals);

private:
	int sig_id = -1; // uninitialized until first usage
};

class BasicSampleFilter : public SampleFilter {
public:

	// params
	int min_sample_time = 0;
	int max_sample_time = (int)(1<<30); // these should always be given in the unit time appearing in the sample file
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

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(BasicFilteringParams);
MEDSERIALIZE_SUPPORT(BasicSampleFilter);

#endif
