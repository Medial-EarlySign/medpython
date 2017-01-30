#ifndef _SAMPLE_FILTER_H_
#define _SAMPLE_FILTER_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "InfraMed/InfraMed/MedPidRepository.h"

//.......................................................................................
//.......................................................................................
// Sample Filter creates a list of samples to work on from a repository or another list
//.......................................................................................
//.......................................................................................
// Define types of sample filter
typedef enum {
	SMPL_FILTER_TRN,
	SMPL_FILTER_TST,
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

	// Filter
	virtual void Filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {return ;}
	void filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) { Filter(rep, inSamples, outSamples); }
	void filter(MedRepository& rep, MedSamples& samples);

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
	void Filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
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
	void Filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
};

#endif
