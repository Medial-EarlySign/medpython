#include "SampleFilter.h"
#include "Logger/Logger/Logger.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/MedPidRepository.h"

#define LOCAL_SECTION LOG_SMPL_FILTER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// SampleFilter
//=======================================================================================
// Filter Types
SampleFilterTypes sample_filter_name_to_type(const string& filter_name) {

	if (filter_name == "train")
		return SMPL_FILTER_TRN;
	else if (filter_name == "test")
		return SMPL_FILTER_TST;
	else if (filter_name == "outliers")
		return SMPL_FILTER_OUTLIERS;
	else
		return SMPL_FILTER_LAST;
}

// Initialization
//.......................................................................................
SampleFilter* SampleFilter::make_filter(string filter_name) {

	return make_filter(sample_filter_name_to_type(filter_name));
}

//.......................................................................................
SampleFilter * SampleFilter::make_filter(string filter_name, string init_string) {

	return make_filter(sample_filter_name_to_type(filter_name), init_string);
}

//.......................................................................................
SampleFilter * SampleFilter::make_filter(SampleFilterTypes filter_type) {

	if (filter_type == SMPL_FILTER_TRN)
		return new BasicTrainFilter;
	else if (filter_type == SMPL_FILTER_TST)
		return new BasicTestFilter;
	else if (filter_type = SMPL_FILTER_OUTLIERS)
		return new OutlierSampleFilter;
	else
		return NULL;

}

//.......................................................................................
SampleFilter * SampleFilter::make_filter(SampleFilterTypes filter_type, string init_string) {

	SampleFilter *newSampleFilter = make_filter(filter_type);
	newSampleFilter->init_from_string(init_string);
	return newSampleFilter;
}

//.......................................................................................

// (De)Serialize
//.......................................................................................
size_t SampleFilter::get_filter_size() {
	return sizeof(filter_type) + get_size();
}

//.......................................................................................
size_t SampleFilter::filter_serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &filter_type, sizeof(SampleFilterTypes)); ptr += sizeof(SampleFilterTypes);
	ptr += serialize(blob + ptr);

	return ptr;
}

// Filter
//.......................................................................................
void SampleFilter::filter(MedSamples& samples) {

	MedSamples out_samples;
	filter(samples, out_samples);

	samples = out_samples;
}

// Filter
//.......................................................................................
void SampleFilter::filter(MedRepository& rep, MedSamples& samples) {

	MedSamples out_samples;
	filter(rep, samples, out_samples);

	samples = out_samples;
}


//=======================================================================================
// BasicTrainFilter
//=======================================================================================
// Filter
//.......................................................................................
void BasicTrainFilter::_filter(MedSamples& inSamples,MedSamples& outSamples) {

	outSamples.time_unit = inSamples.time_unit;

	// Take only samples before outcome
	for (MedIdSamples& idSamples : inSamples.idSamples) {

		MedIdSamples outIdSamples;
		outIdSamples.id = idSamples.id;

		for (MedSample& sample : idSamples.samples) {
			// Negative or pre-outcome
			if (sample.outcome == 0 || sample.outcomeTime > sample.time)
				outIdSamples.samples.push_back(sample);
		}
		
		if (outIdSamples.samples.size() > 0)
			outSamples.idSamples.push_back(outIdSamples);
	}
}

//=======================================================================================
// BasicTestFilter
//=======================================================================================
// Filter
//.......................................................................................
void BasicTestFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {

	// Take them all
	outSamples = inSamples;
}

//=======================================================================================
// OutlierSampleFilter
//=======================================================================================
// Filter
//.......................................................................................
void OutlierSampleFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {

	outSamples.time_unit = inSamples.time_unit;

	// Filter by value of outcome
	for (MedIdSamples& idSample : inSamples.idSamples) {
		MedIdSamples newIdSample;
		newIdSample.id = idSample.id;

		for (MedSample& sample : idSample.samples) {
			if (sample.outcome >= removeMin - NUMERICAL_CORRECTION_EPS && sample.outcome <= removeMax + NUMERICAL_CORRECTION_EPS)
				newIdSample.samples.push_back(sample);
		}

		if (newIdSample.samples.size() > 0)
			outSamples.idSamples.push_back(newIdSample);
	}

}

// Learn
//.......................................................................................
int OutlierSampleFilter::_learn(MedSamples& samples) {

	if (params.type == VAL_CLNR_ITERATIVE)
		return iterativeLearn(samples);
	else if (params.type == VAL_CLNR_QUANTILE)
		return quantileLearn(samples);
	else {
		MERR("Unknown cleaning type %d\n", params.type);
		return -1;
	}
}

// Learn
//.......................................................................................
int OutlierSampleFilter::iterativeLearn(MedSamples& samples) {

	// Get all values
	vector<float> values;
	get_values(samples,values);

	return get_iterative_min_max(values);
}

// Learn
//.......................................................................................
int OutlierSampleFilter::quantileLearn(MedSamples& samples) {

	// Get all values
	vector<float> values;
	get_values(samples,values);

	return get_quantile_min_max(values);
}

// Utility for learning
//.......................................................................................
void OutlierSampleFilter::get_values(MedSamples& samples, vector<float>& values) {

	for (MedIdSamples& idSample : samples.idSamples) {
		for (MedSample& sample : idSample.samples)
			values.push_back(sample.outcome);
	}
}

// (De)Serialization
//.......................................................................................
size_t OutlierSampleFilter::get_size() {

	size_t size = 0;

	size += sizeof(int); // int take_log
	size += 2 * sizeof(float); // float removeMax, removeMin;

	return size;
}

//.......................................................................................
size_t OutlierSampleFilter::serialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(blob + ptr, &params.take_log, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &removeMax, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &removeMin, sizeof(float)); ptr += sizeof(float);

	return ptr;
}

//.......................................................................................
size_t OutlierSampleFilter::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(&params.take_log, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&removeMax, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&removeMin, blob + ptr, sizeof(float)); ptr += sizeof(float);

	return ptr;
}
