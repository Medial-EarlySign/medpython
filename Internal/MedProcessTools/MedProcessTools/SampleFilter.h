#ifndef _SAMPLE_FILTER_H_
#define _SAMPLE_FILTER_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedValueCleaner.h"
#include <MedMat/MedMat/MedMat.h>

#define DEFAULT_SMPL_FLTR_NTHREADS 8

//.......................................................................................
/**  SampleFilter is the parent class for creating a list of samples to work on from a
*		another list and/or a repository <br>
* Basic functionalities: <br>
*		learn : learn the filtering parameters from a given MedSamples + optional
*				MedRepository (this method can be empty) <br>
*		filter : generate a new MedSample from a given MedSamples + optional MedRepository
*				This method must be implemented for each inheriting class <br>
*		Variants for learn and filter : <br>
*			- no repository required <br>
*			- in-place filtering <br>
*/
//.......................................................................................

//.......................................................................................
/** Define types of sample filte </summary>
*/
//.......................................................................................
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

	// Type
	SampleFilterTypes filter_type; ///< The type of the filter

	// create a new sample filter
	/// <summary> create a new sample filter from name </summary>
	static SampleFilter *make_filter(string name); 
	/// <summary> create a new sample filter from type </summary>
	static SampleFilter *make_filter(SampleFilterTypes type); 
	/// <summary> create a new sample filter from name and a parameters string</summary>
	static SampleFilter *make_filter(string name, string params); 
	///  create a new sample filter from type and a parameters string </summary>
	static SampleFilter *make_filter(SampleFilterTypes type, string params); 

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
	virtual int learn(MedRepository& rep, MedSamples& samples) { return _learn(rep,samples); }
	/// <summary> learn without repository </summary>
	virtual int learn(MedSamples& samples) { return _learn(samples); }

	// Filtering
	/// <summary> filter with repository : Should be implemented for inheriting classes that filter using Repository information </summary>
	virtual int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {return _filter(inSamples,outSamples) ;}
	/// <summary> _filter without repository : Should be implemented for all inheriting classes </summary>
	virtual int _filter(MedSamples& inSamples, MedSamples& outSamples) = 0;

	// Filtering : Envelopes (Here because of probelsm with overload + inheritance)
	/// <summary> filter with repository </summary>
	virtual int filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) { return _filter(rep, inSamples, outSamples); }
	/// <summary> in-place filtering with repository </summary>
	int filter(MedRepository& rep, MedSamples& samples);
	/// <summary> filter without repository : Should be implemented for all inheriting classes </summary>
	virtual int filter(MedSamples& inSamples, MedSamples& outSamples) { return _filter(inSamples, outSamples); }
	/// <summary> in-place filtering without repository </summary>
	int filter(MedSamples& samples);

	/// <summary>  Get all signals required for filtering : Should be implemented for inheriting classes that filter using Repository information </summary>
	virtual void get_required_signals(vector<string>& req_sigs) {return; }

	// Serialization (including type)
	ADD_CLASS_NAME(SampleFilter)
	ADD_SERIALIZATION_FUNCS(filter_type)
	virtual void *new_polymorphic(string derived_class_name);

	/// <summary> get size of filter + filter_type </summary>
	size_t get_filter_size();
	/// <summary> seialize filter + filter_type </summary>
	size_t filter_serialize(unsigned char *blob);
};

// Utilities
/// <summary> get SampleFilterType from name </summary>
SampleFilterTypes sample_filter_name_to_type(const string& filter_name);

//.......................................................................................
/** Training set filter
*	 take all controls samples (outcome=0) and all cases before outcomeTime
*/
//.......................................................................................
class BasicTrainFilter : public SampleFilter {
public:

	/// <summary> constructor </summary>
	BasicTrainFilter() { filter_type = SMPL_FILTER_TRN; };

	/// <summary> Filter without repository </summary>
	int _filter(MedSamples& inSamples, MedSamples& outSamples);
	ADD_CLASS_NAME(BasicTrainFilter)
	ADD_SERIALIZATION_FUNCS(filter_type)
};

//.......................................................................................
/** Test set filter
*	- dummy filter - take everything
*/
//.......................................................................................
class BasicTestFilter : public SampleFilter {
public:

	/// <summary> constructor </summary>
	BasicTestFilter() { filter_type = SMPL_FILTER_TST; };
	~BasicTestFilter() {};

	/// <summary> Filter without repository </summary>
	int _filter(MedSamples& inSamples, MedSamples& outSamples);

	ADD_CLASS_NAME(BasicTestFilter)
	ADD_SERIALIZATION_FUNCS(filter_type)
};

#define SMPL_FLTR_TRIMMING_SD_NUM 7
#define SMPL_FLTR_REMOVING_SD_NUM 7
//.......................................................................................
/** Outliers filter
*	- A filter that remove samples with outlier-outcomes (suitable for regression) <br>
*	  Outliers detection is done using MedValueCleaner's methods (through inheritance)
*/
//.......................................................................................
class OutlierSampleFilter : public SampleFilter, public MedValueCleaner {
public:

	/// <summary> constructor </summary>
	OutlierSampleFilter() : SampleFilter() { init_defaults(); };

	/// <summary> Learning : check outlier-detection method and call appropriate learner (iterative/quantile) </summary>
	int _learn(MedSamples& samples);
	/// <summary> Learning : learn outliers using MedValueCleaner's iterative approximation of moments </summary>
	int iterativeLearn(MedSamples& samples);
	/// <summary> Learning : learn outliers using MedValueCleaner's quantile appeoximation of moments </summary>
	int quantileLearn(MedSamples& samples);
	/// <summary> Helper for learning - extract all outcomes from samples  </summary>
	void get_values(MedSamples& samples, vector<float>& values);

	/// <summary> Filter without repository </summary>
	int _filter(MedSamples& inSamples, MedSamples& outSamples);

	/// <summary> init from map </summary>
	int init(map<string, string>& mapper) { return MedValueCleaner::init(mapper); }
	/// <summary> init to defaults </summary>
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
	ADD_CLASS_NAME(OutlierSampleFilter)
	ADD_SERIALIZATION_FUNCS(filter_type, params.take_log, removeMax, removeMin)
};

//.......................................................................................
/** Define types of matching criteria </summary>
*/
//.......................................................................................
typedef enum {
	SMPL_MATCH_SIGNAL, ///< Match by value of signal
	SMPL_MATCH_AGE,	   ///< Match by age
	SMPL_MATCH_TIME,   ///< Match by time
	SMPL_MATCH_FTR,   ///< Match by value of feature
	SMPL_MATCH_LAST
} SampleMatchingType;

//.......................................................................................
/** MatchingParams defines a specific matching criterion
*/
//.......................................................................................
class matchingParams : public SerializableObject {
public:

	SampleMatchingType match_type; ///< matching criterion

	// Matching details
	string signalName; ///< signal name for matching by signal
	string featureName,resolvedFeatureName; ///< feature name for matching by feature
	int timeWindow, windowTimeUnit ; ///< time-window info For matching by signal
	int matchingTimeUnit; ///< time-unit for matching by time
	float resolution ; ///< binnning resolution

	/// Helpers (for matching by signal)
	int signalId; ///< matching signal id
	bool isTimeDependent; ///< flag: is the signal time-dependent (e.g. hemoglobin) or not (e.g. byear)
	int signalTimeUnit; ///< matching signal time-unit

	// Serialization
	// Serialization
	ADD_CLASS_NAME(matchingParams)
	ADD_SERIALIZATION_FUNCS(match_type, signalName, featureName, timeWindow, matchingTimeUnit, resolution)
};

//.......................................................................................
/** Matching filter <br>
*	- match cases and controls according to a set of matching strata <br>
*   - matching creiterion can be <br>
*		- value of signal within window <br>
*		- age  <br>
*		- sample time <br>
* <br>
*	- Each sample is assigned to a bin according to the vector of strata.
*	  The optimal case/control sampling ratio is then found, when each control that is
*	  removed costs 1 point and each case costs eventToControlPriceRatio points. 
*     Once ratio is decided, sampling is perfored randomly within each bin.
*/
//.......................................................................................
class MatchingSampleFilter : public SampleFilter {
public:

	vector<matchingParams> matchingStrata; ///< Matching parameters

	float eventToControlPriceRatio = 100.0; ///< Cost of removing case relative to removing control
	int min_group_size = 5; ///< minimal group size to take - smaller than that, will drop
	float maxControlToEventRatio = -1.0; ///< maximal allowed control/case ratio
	int verbose = 0; ///< control level of debug printing
	float match_to_prior = -1; ///<If given (0-1) will ignore price ratio and will match to this prior

	// helpers
	int samplesTimeUnit; ///< Time unit of samples
	int bdateId; ///< signal-id for byear

	/// <summary> Constructor </summary>
	MatchingSampleFilter() { init_defaults(); };

	// Initialization
	/// <summary> init from map </summary>
	int init(map<string, string>& mapper);
	/// <summary> init to defaults </summary>
	void init_defaults() { filter_type = SMPL_FILTER_MATCH; };

	/// <summary> Add a matching stratum defined by a string </summary>
	int addMatchingStrata(string& init_string);

	// Utilities
	/// <summary> Check if repository is needed for matching (e.g. strata includes signal/age) </summary>
	bool isRepRequired();
	/// <summary> Check if age is needed for matching </summary>
	bool isAgeRequired();
	/// <summary> Indexing of a single sample according to strata </summary>
	/// <returns> 0 upon success, -1 upon finding an illegal strata definition </returns>
	int getSampleSignature(MedSample& sample, MedFeatures& features, int i, MedRepository& rep, string& signature);
	/// <summary> add indexing of a single sample according to a single stratum to sample's index </summary>
	/// <returns> 0 upon success, -1 upon finding an illegal strata definition </returns>
	int addToSampleSignature(MedSample& sample, matchingParams& stratum, MedFeatures& features, int i, MedRepository& rep, string& signature);
	/// <summary> initialize values of requried helpers </summary>
	/// <returns> 0 upon success, -1 if any of the required signals does not appear in the dictionary or any of the required features is not given </returns>
	int initHelpers(MedSamples& inSamples, MedFeatures& features, MedRepository& rep);
	/// <summary>  Get all signals required for matching </summary>
	void get_required_signals(vector<string>& req_sigs);

	/// <summary> Filter with repository </summary>
	int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
	/// <summary> Filter without repository </summary>
	/// <returns>  -1 if repository is required, 0 othereise </returns>
	int _filter(MedSamples& inSamples, MedSamples& outSamples);

	/// <summary> Filter with matrix + repository
	int _filter(MedFeatures& features, MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
	/// <summary> Filter with matrix 
	/// <returns>  -1 if repository is required, 0 othereise </returns>
	int _filter(MedFeatures& features, MedSamples& inSamples, MedSamples& outSamples);

	// Serialization
	ADD_CLASS_NAME(MatchingSampleFilter)
	ADD_SERIALIZATION_FUNCS(filter_type, matchingStrata, eventToControlPriceRatio, maxControlToEventRatio, min_group_size, match_to_prior)
};

//.......................................................................................
/** Required Signal Filter <br>
*	- Keep only samples with a required signal appearing in a time-window. <br>
*	- OBSOLETE - REPLACED BY BasicSampleFilter. KEPT HERE FOR BACKWARD COMPETABILITY
*/
//.......................................................................................
class RequiredSignalFilter : public SampleFilter {
public:

	string signalName; ///< Required signal
	int timeWindow, windowTimeUnit; ///< Time before sample-time (and time-unit)

	/// <summary> Constructor </summary>
	RequiredSignalFilter() { init_defaults(); };

	/// <summary> init from map </summary>
	int init(map<string, string>& mapper);
	/// <summary> init to defaults </summary>
	void init_defaults();

	// Filter
	/// <summary> Filter with repository </summary>
	int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
	/// <summary> Filter without repository </summary>
	/// <returns>  -1 if repository is required, 0 othereise </returns>
	int _filter(MedSamples& inSamples, MedSamples& outSamples);

	// Serialization
	ADD_CLASS_NAME(RequiredSignalFilter)
	ADD_SERIALIZATION_FUNCS(filter_type, signalName, timeWindow, windowTimeUnit)

};

//.......................................................................................
/** BasicFilteringParams defines filtering parameters for BasicFilter with helpers
*/
//.......................................................................................
struct BasicFilteringParams : public SerializableObject {
	string sig_name; ///< Name of signal to filter by	
	int win_from = 0; ///< Time window for deciding on filtering - start
	int win_to = (int)(1<<30); ///< Time window for deciding on filtering - end
	float min_val = -1e10; ///< Allowed values range for signal - minimum
	float max_val = 1e10;///< Allowed values range for signal - maximum
	int min_Nvals = 1; ;///< Required number of instances of signal within time window
	int time_channel = 0; ///< signal time-channel to consider
	int val_channel = 0; ///< signal value channel to consider

	/// <summary> Initialization from string </summary>
	int init_from_string(const string &init_str);
	/// <summary> Test filtering criteria </summary>
	/// <returns> 1 if passing and 0 otherwise </returns>
	int test_filter(MedSample &sample, MedRepository &rep, int win_time_unit);

	ADD_CLASS_NAME(BasicFilteringParams)
	ADD_SERIALIZATION_FUNCS(sig_name, time_channel, val_channel, win_from, win_to, min_val, max_val, min_Nvals)

private:
	int sig_id = -1; ///< signal-id : uninitialized until first usage
	int use_byear= 0; ///< indicator that required signal is age that should be treated differently
};				

//...........................................................................................................
/** BasicSimpleFilter is a general filter to allow the following basics: <br>
* (1) min and max time of outcomeTime <br>
* (2) option to allow for a signal to be in some given range (if there is a sample in some given window) <br>
* (3) option to force N samples of a specific signal within the given range (in the given time window) <br>
*/
//...........................................................................................................
class BasicSampleFilter : public SampleFilter {
public:

	// filtering parameters
	int min_sample_time = 0; ///< minimal allowed time (should always be given in the samples' time-unit)
	int max_sample_time = (int)(1<<30); ///< maximal allowed time (should always be given in the samples' time-unit)
	vector<BasicFilteringParams> bfilters; ///< vector of filters to apply
	int min_bfilter = -1; ///< -1: force each bfilter to pass , other n : at least n bfilters must pass
	int winsTimeUnit = MedTime::Days; ///< time unit to be used 

	// next is initialized with init string
	vector<string> req_sigs; ///< useful to load the repository needed for this filter

	///<summary> init from mapped string </summary>
	int init(map<string, string>& mapper);
	/// <summary>  Get all signals required for filtering </summary>
	void get_required_signals(vector<string>& req_sigs);

	// Filter
	/// <summary> Filter with repository </summary>
	int _filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples);
	/// <summary> Filter without repository </summary>
	/// <returns>  -1 if repository is required, 0 othereise (when bfilters is empty) </returns>
	int _filter(MedSamples& inSamples, MedSamples& outSamples);

	ADD_CLASS_NAME(BasicSampleFilter)
	ADD_SERIALIZATION_FUNCS(filter_type, min_sample_time, max_sample_time, bfilters, min_bfilter, winsTimeUnit)
};

//...........................................................................................................
/** SanitySimpleFilter helps making sanity tests on input data <br>
* <br>
* The basic tests optional are: <br>
* (1) test that the signal actually exist in name (in the signals list in the repository) <br>
* (2) within a given window: minimal number of tests <br>
* (3) within a given window: maximal number of outliers <br>
* (4) count outliers within a given window <br>
*/
//...........................................................................................................
class SanitySimpleFilter : public SerializableObject {
public:
	string sig_name; ///< Name of signal to filter by	
	int win_from = 0; ///< Time window for deciding on filtering - start
	int win_to = (int)(1 << 30); ///< Time window for deciding on filtering - end
	float min_val = -1e10; ///< Allowed values range for signal - minimum
	float max_val = 1e10;///< Allowed values range for signal - maximum
	int min_Nvals = -1; ;///< Required number of instances of signal within time window
	int max_Nvals = -1; ;///< Maximal allowed number of instances of signal within time window
	int time_channel = 0; ///< signal time-channel to consider
	int val_channel = 0; ///< signal value channel to consider
	int max_outliers = -1; ///< maximla allowed number of outliers. -1 means don't do the max_outliers test
	int win_time_unit = MedTime::Days; ///< time unit to be used 
	int samples_time_unit = MedTime::Date; ///< time unit to be used 
	unordered_set<float> allowed_values; /// list of allowed values for the signal
	int values_in_dictionary = 0; /// flag: if 1: make sure all given values are valid - that is are in the signal dictionary.
	int min_left = -1; /// test the min number of instances left that are not outliers



	///<summary> init from mapped string </summary>
	int init_from_string(const string &init_str);
	///<summary> Test filtering criteria </summary>
	///<returns> 1 if passing and 0 otherwise </returns>
	int test_filter(MedSample &sample, MedRepository &rep) {
		int nvals, noutliers; return test_filter(sample, rep, nvals, noutliers);
	}
	///<summary> Test filtering criteria </summary>
	///<returns> 1 if passing and 0 otherwise </returns>
	int test_filter(MedSample &sample, MedRepository &rep, int &nvals, int &noutliers);

	// test_filter return codes
	const static int Passed = 0;
	const static int Failed = 1; // General fail due to reasons different than the following
	const static int Signal_Not_Valid = 2;
	const static int Failed_Min_Nvals = 3;
	const static int Failed_Max_Nvals = 4;
	const static int Failed_Outliers = 5;
	const static int Failed_Age = 6;
	const static int Failed_Age_No_Byear = 7;
	const static int Failed_Allowed_Values = 8;
	const static int Failed_Dictionary_Test = 9;
	const static int Failed_Not_Enough_Non_Outliers_Left = 10;

	ADD_CLASS_NAME(SanitySimpleFilter)
	ADD_SERIALIZATION_FUNCS(sig_name, time_channel, val_channel, win_from, win_to, min_val, max_val, min_Nvals, max_Nvals, allowed_values, values_in_dictionary, max_outliers, win_time_unit)

private:
	int sig_id = -1; ///< signal-id : uninitialized until first usage (0 is kept for the age special case)
	int section_id = -1; /// uninitialized section_id
	int bdate_id = -1;
	bool used_byear = false; ///< If true will use BYEAR when no BDATE
};

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(SampleFilter)
MEDSERIALIZE_SUPPORT(SanitySimpleFilter)
MEDSERIALIZE_SUPPORT(BasicFilteringParams)
MEDSERIALIZE_SUPPORT(BasicSampleFilter)
MEDSERIALIZE_SUPPORT(MatchingSampleFilter)
MEDSERIALIZE_SUPPORT(BasicTrainFilter)
MEDSERIALIZE_SUPPORT(BasicTestFilter)
MEDSERIALIZE_SUPPORT(OutlierSampleFilter)
MEDSERIALIZE_SUPPORT(matchingParams)
MEDSERIALIZE_SUPPORT(RequiredSignalFilter)

#endif
