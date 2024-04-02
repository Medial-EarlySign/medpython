#include "MPSampleFilter.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"

const int MPSampleFilter::SMPL_FILTER_TRN_ = SMPL_FILTER_TRN;
const int MPSampleFilter::SMPL_FILTER_TST_ = SMPL_FILTER_TST;
const int MPSampleFilter::SMPL_FILTER_OUTLIERS_ = SMPL_FILTER_OUTLIERS;
const int MPSampleFilter::SMPL_FILTER_MATCH_ = SMPL_FILTER_MATCH;
const int MPSampleFilter::SMPL_FILTER_REQ_SIGNAL_ = SMPL_FILTER_REQ_SIGNAL;
const int MPSampleFilter::SMPL_FILTER_BASIC_ = SMPL_FILTER_BASIC;
const int MPSampleFilter::SMPL_FILTER_LAST_ = SMPL_FILTER_LAST;


MPSampleFilter::MPSampleFilter() { o = nullptr; }
MPSampleFilter::MPSampleFilter(const MPSampleFilter& other) { *o = *(other.o); };
MPSampleFilter::MPSampleFilter(SampleFilter* ptr) { o = ptr; };
MPSampleFilter::~MPSampleFilter() { delete o; };
MPSampleFilter& MPSampleFilter::operator=(const MPSampleFilter& other) {
	*o = *(other.o);
	return *this;
};

MPSampleFilter MPSampleFilter::from_name(string name) { MPSampleFilter ret; ret.o = SampleFilter::make_filter(name); return ret; };
MPSampleFilter MPSampleFilter::from_type(int type) { MPSampleFilter ret; ret.o = SampleFilter::make_filter((SampleFilterTypes)type); return ret; };
MPSampleFilter MPSampleFilter::from_name_params(string name, string params) { MPSampleFilter ret; ret.o = SampleFilter::make_filter(name, params); return ret; };
MPSampleFilter MPSampleFilter::from_type_params(int type, string params) { MPSampleFilter ret; ret.o = SampleFilter::make_filter((SampleFilterTypes)type, params); return ret; };

int MPSampleFilter::filter(MPPidRepository& rep, MPSamples& samples) {
	o->filter(*(rep.o), *(samples.o));
	return 0;
}

MPSerializableObject MPSampleFilter::asSerializable() { return MPSerializableObject(o); }