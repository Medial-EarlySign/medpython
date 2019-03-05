#ifndef __MED__MPSAMPLEFILTER__H__
#define __MED__MPSAMPLEFILTER__H__

#include "MedPyCommon.h"
#include "MPPidRepository.h"
#include "MPSamples.h"

class SampleFilter;

class MPSampleFilter {
private:
	//always  owned on this type
	//bool o_owned = true;
	MEDPY_IGNORE(SampleFilter * o);
public:
	static const int SMPL_FILTER_TRN_;
	static const int SMPL_FILTER_TST_;
	static const int SMPL_FILTER_OUTLIERS_;
	static const int SMPL_FILTER_MATCH_;
	static const int SMPL_FILTER_REQ_SIGNAL_;
	static const int SMPL_FILTER_BASIC_;
	static const int SMPL_FILTER_LAST_;

	static MPSampleFilter from_name(string name);
	static MPSampleFilter from_type(int type);
	static MPSampleFilter from_name_params(string name, string params);
	static MPSampleFilter from_type_params(int type, string params);

	MEDPY_IGNORE(MPSampleFilter());
	MEDPY_IGNORE(MPSampleFilter(const MPSampleFilter& other));
	MEDPY_IGNORE(MPSampleFilter(SampleFilter* ptr));
	~MPSampleFilter();
	MEDPY_IGNORE(MPSampleFilter& operator=(const MPSampleFilter& other));

	/// <summary> in-place filtering with repository </summary>
	MEDPY_DOC(filter, "filter(rep, samples) -> int\n"
		"in-place filtering with repository. returns -1 if fails.");
	int filter(MPPidRepository& rep, MPSamples& samples);

	MPSerializableObject asSerializable();
	/*

	// Initialization
	/// <summary> initialize from a params object :  Should be implemented for inheriting classes that have parameters </summary>
	virtual int init(void *params) { return 0; };
	/// <summary> initialize from a map :  Should be implemented for inheriting classes that have parameters </summary>
	virtual int init(map<string, string>& mapper) { return 0; };
	/// <summary> initialize to default values :  Should be implemented for inheriting classes that have parameters </summary>
	virtual void init_defaults() {};

	// Learning : Actually learn 
	/// <summary> learn with repository : Should be implemented for inheriting classes that learn parameters using Repository information </summary>
	virtual int _learn(MedRepository& rep, MedSamples& samples) { return _learn(samples); }
	/// <summary> learn without repository : Should be implemented for inheriting classes that learn parameters</summary>
	virtual int _learn(MedSamples& samples) { return 0; }

	// Learning : Envelopes (Here because of probelsm with overload + inheritance)
	/// <summary> learn with repository  </summary>
	virtual int learn(MedRepository& rep, MedSamples& samples) { return _learn(rep, samples); }
	/// <summary> learn without repository </summary>
	virtual int learn(MedSamples& samples) { return _learn(samples); }

	// Filtering
	/// <summary> filter with repository : Should be implemented for inheriting classes that filter using Repository information </summary>
	virtual int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) { return _filter(inSamples, outSamples); }
	/// <summary> _filter without repository : Should be implemented for all inheriting classes </summary>
	virtual int _filter(MedSamples& inSamples, MedSamples& outSamples) = 0;

	// Filtering : Envelopes (Here because of probelsm with overload + inheritance)
	/// <summary> filter with repository </summary>
	virtual int filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) { return _filter(rep, inSamples, outSamples); }

	/// <summary> filter without repository : Should be implemented for all inheriting classes </summary>
	virtual int filter(MedSamples& inSamples, MedSamples& outSamples) { return _filter(inSamples, outSamples); }
	/// <summary> in-place filtering without repository </summary>
	int filter(MedSamples& samples);

	/// <summary>  Get all signals required for filtering : Should be implemented for inheriting classes that filter using Repository information </summary>
	virtual void get_required_signals(vector<string>& req_sigs) { return; }
	*/
}; 
#endif //__MED__MPSAMPLEFILTER__H__