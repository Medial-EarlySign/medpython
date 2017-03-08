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

	virtual int init(void *cleaner_params) { return 0; };
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
	virtual void _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {_filter(inSamples,outSamples) ;}
	void filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) { _filter(rep, inSamples, outSamples); }
	void filter(MedRepository& rep, MedSamples& samples);

	// Without Repository
	virtual void _filter(MedSamples& inSamples, MedSamples& outSamples) { return; }
	void filter(MedSamples& inSamples, MedSamples& outSamples) { _filter(inSamples, outSamples); }
	void filter(MedSamples& samples);

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
	void _filter(MedSamples& inSamples, MedSamples& outSamples);
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
	void _filter(MedSamples& inSamples, MedSamples& outSamples);
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
	void _filter(MedSamples& inSamples, MedSamples& outSamples);

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

#endif
