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

//.......................................................................................
// Filter
void SampleFilter::filter(MedRepository& rep, MedSamples& samples) {

	MedSamples out_samples;
	filter(rep, samples, out_samples);

	samples = out_samples;
}


//=======================================================================================
// BasicTrainFilter
//=======================================================================================
void BasicTrainFilter::Filter(MedRepository& rep, MedSamples& inSamples,MedSamples& outSamples) {

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
void BasicTestFilter::Filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {

	// Take them all
	outSamples = inSamples;
}
