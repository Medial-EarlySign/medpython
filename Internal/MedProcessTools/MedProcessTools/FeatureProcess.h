#ifndef _FTR_PROCESS_H_
#define _FTR_PROCESS_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedProcessTools/MedProcessTools/MedValueCleaner.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"
#include "MedStat/MedStat/MedPerformance.h"
#include "MedStat/MedStat/MedStat.h"
#include <unordered_set>

#define DEFAULT_FEAT_CLNR_NTHREADS 24

/** @enum
* Rep Processors types enum
*/
typedef enum {
	FTR_PROCESS_MULTI, ///<"multi_processor" or "multi" to create MultiFeatureProcessor
	FTR_PROCESS_BASIC_OUTLIER_CLEANER, ///<"basic_outlier_cleaner" or "basic_cleaner" to create FeatureBasicOutlierCleaner
	FTR_PROCESS_NORMALIZER, ///<"normalizer" to create FeatureNormalizer
	FTR_PROCESS_IMPUTER, ///<"imputer" to create FeatureImputer
	FTR_PROCESS_DO_CALC, ///<"do_calc" to create DoCalcFeatProcessor
	FTR_PROCESS_UNIVARIATE_SELECTOR, ///<"univariate_selector" to create UnivariateFeatureSelector
	FTR_PROCESSOR_MRMR_SELECTOR, ///<"mrmr" or "mrmr_selector" to create MRMRFeatureSelector
	FTR_PROCESSOR_LASSO_SELECTOR, ///<"lasso" to create LassoSelector
	FTR_PROCESSOR_TAGS_SELECTOR, ///<"tags_selector" to create TagFeatureSelector
	FTR_PROCESSOR_IMPORTANCE_SELECTOR, ///<"importance_selector" to create ImportanceFeatureSelector
	FTR_PROCESS_REMOVE_DGNRT_FTRS, ///<"remove_deg" to create DgnrtFeatureRemvoer
	FTR_PROCESS_ITERATIVE_IMPUTER, ///<"iterative_imputer" to create IterativeImputer
	FTR_PROCESS_ENCODER_PCA, ///<"pca" to create FeaturePCA
	FTR_PROCESS_ONE_HOT, ///< make one-hot features from a given feature
	FTR_PROCESS_LAST
} FeatureProcessorTypes;

/** @file
* A virtual class of processes on MedFeatures;
* E.g. Cleaning
*/
class FeatureProcessor : public SerializableObject {
public:

	/// Feature name ( + name as appears in MedFeatures) ;
	string feature_name = "unset_feature_name";
	string resolved_feature_name;

	// Type
	FeatureProcessorTypes processor_type = FTR_PROCESS_LAST;

	// Threading
	int learn_nthreads, clean_nthreads;

	// Constructor/Destructor
	FeatureProcessor() { init_defaults(); };
	virtual ~FeatureProcessor() { clear(); };
	virtual void clear() {};
	void init_defaults() { learn_nthreads = DEFAULT_FEAT_CLNR_NTHREADS;  clean_nthreads = DEFAULT_FEAT_CLNR_NTHREADS; processor_type = FTR_PROCESS_LAST; };

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *processor; }

	// Virtual Set Feature Name
	virtual void set_feature_name(const string& feature_name) { this->feature_name = feature_name; }
	virtual string get_feature_name() { return this->feature_name; }
	virtual void get_feature_names(vector<string> & feature_names) { feature_names.clear(); feature_names.push_back(feature_name); };

	// Learn cleaning model
	virtual int Learn(MedFeatures& features, unordered_set<int>& ids) { return 0; }

	/// <summary>
	/// PostProcess of MedFeatures - on all ids. stores information to post process
	/// new features. calls virtual function "Learn" for the specific implementation
	/// </summary>
	/// <returns>
	/// 0 if succesfull, otherwise errorcode -1
	/// </returns>
	int learn(MedFeatures& features);
	int learn(MedFeatures& features, unordered_set<int>& ids) { return Learn(features, ids); }

	// Apply cleaning model
	virtual int Apply(MedFeatures& features, unordered_set<int>& ids) { return 0; }

	/// <summary>
	/// PostProcess of MedFeatures - on all ids. apply the post process on the
	/// new features. calls virtaul function "Apply" for the specific implementation
	/// </summary>
	/// <returns>
	/// 0 if succesfull, otherwise errorcode -1
	/// </returns>
	int apply(MedFeatures& features);
	int apply(MedFeatures& features, unordered_set<int>& ids) { return Apply(features, ids); }

	// Init
	static FeatureProcessor *make_processor(string processor_name);
	static FeatureProcessor *make_processor(FeatureProcessorTypes type);
	static FeatureProcessor *make_processor(string processor_name, string params);
	static FeatureProcessor *make_processor(FeatureProcessorTypes type, string params);

	virtual int init(void *processor_params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };

	/// Filter according to a subset of features
	virtual int filter(unordered_set<string>& features) { return (features.find(feature_name) == features.end()) ? 0 : 1; };

	/// Utility : get corresponding name in MedFeatures
	string resolve_feature_name(MedFeatures& features, string substr);

	/// allows testing if this feature processor is a selector
	virtual bool is_selector() { return false; }

	// Serialization (including type)
	ADD_CLASS_NAME(FeatureProcessor)
	ADD_SERIALIZATION_FUNCS(feature_name, resolved_feature_name, processor_type)
	void *new_polymorphic(string derived_class_name);

	size_t get_processor_size();
	size_t processor_serialize(unsigned char *blob);


	// debug prints
	virtual void dprint(const string &pref, int rp_flag);

};

// Utilities
FeatureProcessorTypes feature_processor_name_to_type(const string& cleaner_name);

/**
* A Processor which contains a vector of simpler processors
* Useful for applying same cleaners on a set of features, for example
*
* To Use this selector specify <b>"multi"</b> or <b>multi_processor</b> in the fp_type
*/
class MultiFeatureProcessor : public FeatureProcessor {
public:

	// For generating processors only at learning, we need type + init_string
	FeatureProcessorTypes members_type;
	string init_string;
	int duplicate;
	string tag;

	// Processors (if empty, will be generated upon learning for all featuers)
	vector<FeatureProcessor *> processors;

	// Constructor/Destructor
	MultiFeatureProcessor() { init_defaults(); };
	~MultiFeatureProcessor() { clear(); };

	void init_defaults() { processor_type = FTR_PROCESS_MULTI; duplicate = 0; members_type = FTR_PROCESS_LAST; init_string = ""; tag = ""; };

	void clear();

	/// The parsed fields from init command.
	/// @snippet FeatureProcess.cpp MultiFeatureProcessor::init
	int init(map<string, string>& mapper);

	// Copy
	virtual void copy(FeatureProcessor *processor);

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	virtual void get_feature_names(vector<string>& feature_names);

	// Add processors
	void add_processors_set(FeatureProcessorTypes type, vector<string>& features);
	void add_processors_set(FeatureProcessorTypes type, vector<string>& features, string init_string);

	// Filter according to a subset of features
	int filter(unordered_set<string>& features);

	// debug print
	void dprint(const string &pref, int fp_flag);

	// Serialization
	ADD_CLASS_NAME(MultiFeatureProcessor)
	ADD_SERIALIZATION_FUNCS(processor_type, members_type, init_string, duplicate, tag, processors)
};

#define DEF_FTR_TRIMMING_SD_NUM 7
#define DEF_FTR_REMOVING_SD_NUM 14
/**
* A simple cleaner considering each value of a certain feature separatley
*
* To Use this selector specify <b>"basic_outlier_cleaner"</b> or <b>basic_cleaner</b> in the fp_type
*/
class FeatureBasicOutlierCleaner : public FeatureProcessor, public MedValueCleaner {
public:

	// Constructor
	FeatureBasicOutlierCleaner() : FeatureProcessor() { init_defaults(); }
	FeatureBasicOutlierCleaner(string& feature_name) : FeatureProcessor() { set_feature_name(feature_name); init_defaults(); }
	FeatureBasicOutlierCleaner(string& feature_name, string init_string) : FeatureProcessor() { set_feature_name(feature_name);  init_defaults();  init_from_string(init_string); }
	FeatureBasicOutlierCleaner(string& feature_name, ValueCleanerParams *_params) : FeatureProcessor() { set_feature_name(feature_name);  MedValueCleaner::init(_params); }

	void init_defaults() {
		processor_type = FTR_PROCESS_BASIC_OUTLIER_CLEANER;
		params.missing_value = MED_MAT_MISSING_VALUE;
		params.trimming_sd_num = DEF_FTR_TRIMMING_SD_NUM;
		params.removing_sd_num = DEF_FTR_REMOVING_SD_NUM;
		params.nbrs_sd_num = 0;
		params.take_log = 0;
		params.doTrim = params.doRemove = true;
		params.type = VAL_CLNR_ITERATIVE;

	};

	// Init
	int init(void *processor_params) { return MedValueCleaner::init(processor_params); };
	/// The parsed fields from init command.
	/// @snippet MedValueCleaner.cpp MedValueCleaner::init
	int init(map<string, string>& mapper);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<FeatureBasicOutlierCleaner *>(processor)); }

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);
	int iterativeLearn(MedFeatures& features, unordered_set<int>& ids);
	int quantileLearn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	ADD_CLASS_NAME(FeatureBasicOutlierCleaner)
	ADD_SERIALIZATION_FUNCS(processor_type, feature_name, resolved_feature_name, params.doTrim, params.doRemove, trimMax, trimMin, removeMax, removeMin)

};

/**
* Feature Normalizer
*
* To Use this selector specify <b>"normalizer"</b> in the fp_type
*/
class FeatureNormalizer : public FeatureProcessor {
public:

	/// Missing Value
	float missing_value;

	/// Normalize Standard Deviation
	bool normalizeSd;

	/// Fill missing values with mean 
	bool fillMissing;

	/// Moments
	float mean, sd;

	/// Utility : maximum number of samples to take for moments calculations
	int max_samples = 10000;

	// Constructor
	FeatureNormalizer() : FeatureProcessor() { init_defaults(); }
	FeatureNormalizer(const  string& feature_name) : FeatureProcessor() { init_defaults(); set_feature_name(feature_name); }
	FeatureNormalizer(const  string& feature_name, string init_string) : FeatureProcessor() { init_from_string(init_string);  set_feature_name(feature_name); }

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	/// The parsed fields from init command.
	/// @snippet FeatureProcess.cpp FeatureNormalizer::init
	int init(map<string, string>& mapper);
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; normalizeSd = true; fillMissing = false; processor_type = FTR_PROCESS_NORMALIZER; };

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<FeatureNormalizer *>(processor)); }

	// Serialization
	ADD_CLASS_NAME(FeatureNormalizer)
	ADD_SERIALIZATION_FUNCS(processor_type, feature_name, resolved_feature_name, mean, sd, normalizeSd, fillMissing)

};

//.......................................................................................
//.......................................................................................
// Feature Imputer
//.......................................................................................
//.......................................................................................

typedef enum {
	IMPUTE_MMNT_MEAN,
	IMPUTE_MMNT_MEDIAN,
	IMPUTE_MMNT_COMMON,
	IMPUTE_MMNT_SAMPLE,
	IMPUTE_MMNT_LAST
} imputeMomentTypes;

class featureStrata :SerializableObject {
public:
	string name;
	float resolution, min, max;
	int nValues;

	featureStrata() {};
	featureStrata(string& _name, float _resolution, float _min, float _max) { name = _name; resolution = _resolution; min = _min; max = _max; }

	void SetNValues() { nValues = ((int)(max / resolution) - (int)(min / resolution) + 1); }
	// returns the correct strata for a value. 
	// E.g. if "strata": "Age,40,80,5" 42 will return 0, the first bin
	int getIndex(float value, float missing_val) {
		if (value == missing_val)
			return nValues / 2;
		else {
			if (value >= max)
				return nValues - 1;
			else if (value <= min)
				return 0;
			else
				return ((int)(value / resolution) - (int)(min / resolution));
		}
	}

	// Serialization
	ADD_CLASS_NAME(featureStrata)
	ADD_SERIALIZATION_FUNCS(name, resolution, min, max, nValues)

};
/// When building startas on a set of several features, we build a cartesian product of their combinations:
/// e.g. when "strata": "Age,40,80,5:Gender,1,2,1"
/// starta [Age] factor [1] starta [Gender] factor[9]
/// for a total of [18] stratas
class featureSetStrata :SerializableObject {
public:
	vector<featureStrata> stratas;
	vector<int> factors;

	size_t nStratas() { return stratas.size(); }

	void getFactors() {

		if (stratas.size() == 0)
			return;

		factors.resize(stratas.size());

		for (auto& strata : stratas)
			strata.SetNValues();

		factors[0] = 1;
		for (int i = 1; i < stratas.size(); i++)
			factors[i] = factors[i - 1] * stratas[i - 1].nValues;
	}

	int nValues() {
		if (stratas.size() == 0)
			return 1;
		else
			return factors.back() * stratas.back().nValues;
	}

	// Serialization
	ADD_CLASS_NAME(featureSetStrata)
	ADD_SERIALIZATION_FUNCS(stratas, factors)
};

/**
* Feature Imputer to complete missing values
*
* To Use this selector specify <b>"imputer"</b> in the fp_type
*/
class FeatureImputer : public FeatureProcessor {
public:

	// Missing Value
	float missing_value;
	bool verbose; ///< If true will print how many missing value were in each feature
	bool verbose_learn; ///< If true will call print after learn

	// Strata for setting moment
	featureSetStrata imputerStrata;

	// minimum samples required for learning
	int min_samples = 50;

	// Moment
	imputeMomentTypes moment_type;
	float default_moment;
	vector<float> moments;
	// for sampling-imputation
	vector < pair<float, float> > default_histogram;
	vector < vector<pair<float, float> > > histograms;

	vector<int> strata_sizes;

	/// Utility : maximum number of samples to take for moments calculations
	int max_samples = 100000;

	// Constructor
	FeatureImputer() : FeatureProcessor() { init_defaults(); }
	FeatureImputer(const  string& feature_name) : FeatureProcessor() { init_defaults(); set_feature_name(feature_name); }
	FeatureImputer(const  string& feature_name, string init_string) : FeatureProcessor() { init_from_string(init_string);  set_feature_name(feature_name); }

	// Add stratifier
	void addStrata(string& init_string);
	void addStrata(featureStrata& strata) { imputerStrata.stratas.push_back(strata); }
	void addStrata(string& name, float resolution, float min, float max) { imputerStrata.stratas.push_back(featureStrata(name, resolution, min, max)); }

	/// The parsed fields from init command.
	/// @snippet FeatureProcess.cpp FeatureImputer::init
	int init(map<string, string>& mapper);
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; moment_type = IMPUTE_MMNT_MEAN;  processor_type = FTR_PROCESS_IMPUTER; verbose = true; verbose_learn = false; };
	imputeMomentTypes getMomentType(string& entry);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<FeatureImputer *>(processor)); }

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	int version() { return  1; }
	ADD_CLASS_NAME(FeatureImputer)
	ADD_SERIALIZATION_FUNCS(processor_type, feature_name, resolved_feature_name, missing_value, imputerStrata, moment_type, moments, histograms, strata_sizes, default_moment, default_histogram)

	/// debug and print
	void print();

};

/**
* Feature Selector abstract class
*/
class FeatureSelector : public FeatureProcessor {
public:

	/// Missing Value
	float missing_value = (float)MED_MAT_MISSING_VALUE;

	/// Required Features
	unordered_set<string> required;

	/// Selected Features (ordered)
	vector<string> selected;

	/// Target number to select (if 0, ignored)
	int numToSelect = 0;

	// Constructor
	FeatureSelector() : FeatureProcessor() { missing_value = MED_MAT_MISSING_VALUE; numToSelect = 0; }

	/// Find set of selected features- Calls _learn function, and may be overrided directly
	virtual int Learn(MedFeatures& features, unordered_set<int>& ids);

	/// Apply selection
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	bool is_selctor() { return true; }

	ADD_CLASS_NAME(FeatureSelector)
	ADD_SERIALIZATION_FUNCS(processor_type, missing_value, required, selected, numToSelect)

private:
	/// Find set of selected features
	virtual int _learn(MedFeatures& features, unordered_set<int>& ids) { return 0; }
};

/**
* Feature Selector : lasso
*
* To Use this selector specify <b>"lasso"</b> in the fp_type
*/
class LassoSelector : public FeatureSelector {
public:
	LassoSelector() : FeatureSelector() { init_defaults(); }

	/// The parsed fields from init command.
	/// @snippet FeatureSelector.cpp LassoSelector::init
	int init(map<string, string>& mapper);

	/// Initial lambda
	float initMaxLambda = (float)0.005;

	/// Features less controled in the selection stage (set labmda -> lambda*lambdaRatio)
	float lambdaRatio = (float)0.1;
	vector<string> lax_lasso_features;

	int nthreads = 12;

	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; processor_type = FTR_PROCESSOR_LASSO_SELECTOR; };

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<LassoSelector *>(processor)); }

	// Serialization
	ADD_CLASS_NAME(LassoSelector)
	ADD_SERIALIZATION_FUNCS(processor_type, initMaxLambda, nthreads, missing_value, required, selected, numToSelect)

private:
	// Find set of selected features
	int _learn(MedFeatures& features, unordered_set<int>& ids);
};

/**
* Feature Selector : Remove Degenerate features
*
* To Use this selector specify <b>"remove_deg"</b> in the fp_type
*/
class DgnrtFeatureRemvoer : public FeatureSelector {
public:

	// Percantage covered by single value to define as degenerate
	float percentage = 1.0F;

	// Constructor
	DgnrtFeatureRemvoer() : FeatureSelector() { init_defaults(); }

	/// The parsed fields from init command.
	/// @snippet FeatureSelector.cpp DgnrtFeatureRemvoer::init
	int init(map<string, string>& mapper);
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; processor_type = FTR_PROCESS_REMOVE_DGNRT_FTRS; numToSelect = 0; };

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<DgnrtFeatureRemvoer *>(processor)); }

	// Serialization
	ADD_CLASS_NAME(DgnrtFeatureRemvoer)
	ADD_SERIALIZATION_FUNCS(processor_type, percentage, missing_value, selected)

private:
	// Find set of selected features
	int _learn(MedFeatures& features, unordered_set<int>& ids);
};

typedef enum {
	UNIV_SLCT_PRSN = 0,
	UNIV_SLCT_MI = 1,
	UNIV_SLCT_DCORR = 2,
	UNIV_SLCT_LAST
} UnivariateSelectionMethod;

class univariateSelectionParams : SerializableObject {
public:
	UnivariateSelectionMethod method;
	float minStat;

	/// for mutual information
	int nBins = 10;
	MedBinningType binMethod = BIN_EQUIDIST;

	/// for correlation
	int takeSquare = 0;

	/// for samples distance correlation
	float pDistance;

	/// Utility : maximum number of samples to take for moments calculations
	int max_samples = 10000;

	UnivariateSelectionMethod get_method(string name) {

		boost::algorithm::to_lower(name);
		if (name == "pearson")
			return UNIV_SLCT_PRSN;
		else if (name == "mi" || name == "mutual_information" || name == "mutualinformation")
			return UNIV_SLCT_MI;
		else if (name == "dcorr" || name == "dist_corr" || name == "distcorr")
			return UNIV_SLCT_DCORR;
		else
			return UNIV_SLCT_LAST;
	}

	MedBinningType get_binning_method(string name) {

		boost::algorithm::to_lower(name);
		if (name == "equi_dist")
			return BIN_EQUIDIST;
		else if (name == "equi_size")
			return BIN_EQUISIZE;
		else
			return BIN_LAST;
	}

	ADD_CLASS_NAME(univariateSelectionParams)
	ADD_SERIALIZATION_FUNCS(method, minStat, nBins, binMethod, takeSquare, pDistance, max_samples)
};

/**
* Feature Selector : Univariate
*
* To Use this selector specify <b>"univariate_selector"</b> in the fp_type
*/
class UnivariateFeatureSelector : public FeatureSelector {
public:

	/// Selection Params
	univariateSelectionParams params;

	// Constructor
	UnivariateFeatureSelector() : FeatureSelector() { init_defaults(); }

	/// The parsed fields from init command.
	/// @snippet FeatureSelector.cpp UnivariateFeatureSelector::init
	int init(map<string, string>& mapper);
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; processor_type = FTR_PROCESS_UNIVARIATE_SELECTOR;  params.method = UNIV_SLCT_PRSN;  params.minStat = 0.05F; };

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<UnivariateFeatureSelector *>(processor)); }

	// Serialization
	ADD_CLASS_NAME(UnivariateFeatureSelector)
	ADD_SERIALIZATION_FUNCS(processor_type, params, missing_value, required, selected, numToSelect)

private:
	// Scores 
	int getAbsPearsonCorrs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats);
	int getMIs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats);
	int getDistCorrs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats);
	// Find set of selected features
	int _learn(MedFeatures& features, unordered_set<int>& ids);
};

typedef enum {
	MRMR_MAX = 0,
	MRMR_MEAN = 1,
	MRMR_LAST
} MRMRPenaltyMethod;

/**
* Feature Selector : MRMR
*
* To Use this selector specify <b>"mrmr"</b> or <b>"mrmr_selector"</b> in the fp_type
*/
class MRMRFeatureSelector : public FeatureSelector {
public:
	/// Selection Params
	univariateSelectionParams params;
	float penalty;
	MRMRPenaltyMethod penaltyMethod;

	// Constructor
	MRMRFeatureSelector() : FeatureSelector() { init_defaults(); }

	/// The parsed fields from init command.
	/// @snippet FeatureSelector.cpp MRMRFeatureSelector::init
	int init(map<string, string>& mapper);
	void init_defaults();
	MRMRPenaltyMethod get_penalty_method(string _method);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<MRMRFeatureSelector *>(processor)); }

	// Serialization
	ADD_CLASS_NAME(MRMRFeatureSelector)
	ADD_SERIALIZATION_FUNCS(processor_type, params, penalty, penaltyMethod, missing_value, required, selected, numToSelect)

private:
	// Find set of selected features
	int _learn(MedFeatures& features, unordered_set<int>& ids);
private:
	// Scores 
	int fillStatsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index);
	int fillAbsPearsonCorrsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index);
	int fillMIsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index);
	int fillDistCorrsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index);
};

//.......................................................................................
//.......................................................................................
// Additional implimentations in other h files
//.......................................................................................
//.......................................................................................
#include "IterativeImputer.h"

//.......................................................................................
//.......................................................................................
// Utilities
//.......................................................................................
//.......................................................................................

#define DEF_MAX_SAMPLE 1000
void get_all_values(MedFeatures& features, string& signalName, unordered_set<int>& ids, vector<float>& values, int max_sample = DEF_MAX_SAMPLE);
void get_all_outcomes(MedFeatures& features, unordered_set<int>& ids, vector<float>& values, int max_sample = DEF_MAX_SAMPLE);
void smearBins(vector<int>& bins, int nBins, int reqNbins);

/************************************************************************************//**
* TagFeatureSelector - selector which leave us only with the selected "tags" given as
* param (if empty do nothing) and removes removed_tags (if empty do nothing)
* note that you can use regex notation to specify the tags
* To Use this selector specify <b>"tags_selector"</b> in the fp_type
****************************************************************************************/
class TagFeatureSelector : public FeatureSelector {
public:
	int verbose = 0;
	vector<string> selected_tags; ///< the selected tags
	vector<string> removed_tags; ///< tags to remove
	// Constructor
	TagFeatureSelector() : FeatureSelector() { init_defaults(); }

	/// The parsed fields from init command.
	/// @snippet FeatureSelector.cpp TagFeatureSelector::init
	int init(map<string, string>& mapper);
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; processor_type = FTR_PROCESSOR_TAGS_SELECTOR; };

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<TagFeatureSelector *>(processor)); }

	// Serialization
	ADD_CLASS_NAME(TagFeatureSelector)
	ADD_SERIALIZATION_FUNCS(processor_type, selected_tags, selected)
private:
	// Find set of selected features
	int _learn(MedFeatures& features, unordered_set<int>& ids);
};

/** ImportanceFeatureSelector - selector which uses feature importance method for sepcific
* model to rank the feature importance and select them
*
* To Use this selector specify <b>"importance_selector"</b> in the fp_type
*/
class ImportanceFeatureSelector : public FeatureSelector {
public:
	string predictor; ///<the predictor type - same as in the json file: qrf,lightgbm...
	string predictor_params; ///<the predictor parameters
	string importance_params; ///<additional importance parameters for the feature importance
	float minStat;///<minimal threshold score to select the feature
	bool verbose; ///<print all feature importance
	// Constructor
	ImportanceFeatureSelector() : FeatureSelector() { init_defaults(); }

	/// The parsed fields from init command.
	/// @snippet FeatureSelector.cpp ImportanceFeatureSelector::init
	int init(map<string, string>& mapper);
	virtual void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; processor_type = FTR_PROCESSOR_IMPORTANCE_SELECTOR; };

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<ImportanceFeatureSelector *>(processor)); }


	// Serialization
	ADD_CLASS_NAME(ImportanceFeatureSelector)
	ADD_SERIALIZATION_FUNCS(processor_type, predictor, predictor_params, importance_params, minStat, selected)

private:
	// Find set of selected features
	int _learn(MedFeatures& features, unordered_set<int>& ids);
};

/**
* FeatureEncoder - General class for encoding features - PCA, autoencoder...
*/
class FeatureEncoder : public FeatureProcessor {
public:

	/// generated names
	vector<string> names;

	// Constructor
	FeatureEncoder() : FeatureProcessor() { init_defaults(); }

	/// Generate set of selected features - calls _learn
	virtual int Learn(MedFeatures& features, unordered_set<int>& ids);

	/// Apply selection - calls _apply
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	ADD_CLASS_NAME(FeatureEncoder)
	ADD_SERIALIZATION_FUNCS(processor_type, names)

private:
	/// Specific learner of the encoder
	virtual int _learn(MedFeatures& features, unordered_set<int>& ids) { return 0; }
	/// Specific apply of the encoder
	virtual int _apply(MedFeatures& features, unordered_set<int>& ids) { return 0; }
};

/**
* PCA Parameters class
*/
class FeaturePCAParams : SerializableObject {
public:
	int pca_top; ///<Max Number of PCA Components to take
	float pca_cutoff;///<PCA variance threshold to stop
	int subsample_count;///<subsample in the pca rows to speed up

	ADD_CLASS_NAME(FeaturePCAParams)
	ADD_SERIALIZATION_FUNCS(pca_top, pca_cutoff)
};

/**
* FeaturePCA - PCA encoder
*
* To Use this selector specify <b>"pca"</b> in the fp_type
*/
class FeaturePCA :public FeatureEncoder {
public:
	///PCA parameters
	FeaturePCAParams params;

	// Constructor
	FeaturePCA() : FeatureEncoder() { init_defaults(); }

	/// The parsed fields from init command.
	/// @snippet FeatureEncoder.cpp FeaturePCA::init
	int init(map<string, string>& mapper);
	void init_defaults() { processor_type = FTR_PROCESS_ENCODER_PCA;  params.pca_cutoff = 0; params.pca_top = 100; };

	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<FeaturePCA *>(processor)); }

	ADD_CLASS_NAME(FeaturePCA)
	ADD_SERIALIZATION_FUNCS(processor_type, names, params, selected_indexes, W)

private:
	MedMat<float> W;
	vector<int> selected_indexes;

	int _learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply selection
	int _apply(MedFeatures& features, unordered_set<int>& ids);
};

/**
* OneHotFeatProcessor:
*
* Create one-hot index features from a given feature
*/

class OneHotFeatProcessor :public FeatureProcessor {
public:

	string index_feature_prefix = ""; ///< prefix of index features (names are prefix_value)
	string other_feature_name = ""; ///< name of 'other' feature (if needed)
	string removed_feature_name = ""; ///< name of feature to remove (if needed)
	bool rem_origin = true; ///< if true, remove original feature after creating indeices
	bool add_other = false; ///< if true, add an extra feature for values not in learning-set
	bool allow_other = false; ///< if true, values in test, but not in learning-set are allowed
	bool remove_last = false; ///< if true, remove the feature corresponding to the last value to avoid linear dependency
	int max_values = 32; ///< maximal allowed number of different values

	map<float, string> value2feature;

	// Constructor
	OneHotFeatProcessor() { init_defaults(); }

	/// The parsed fields from init command.
	/// @snippet FeatureProcessor.cpp OneHotFeatProcessor::init
	int init(map<string, string>& mapper);

	void init_defaults() { processor_type = FTR_PROCESS_ONE_HOT; }
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<OneHotFeatProcessor *>(processor)); }

	ADD_CLASS_NAME(OneHotFeatProcessor)
	ADD_SERIALIZATION_FUNCS(processor_type, feature_name, index_feature_prefix, other_feature_name, removed_feature_name, rem_origin, add_other, remove_last, allow_other, value2feature)
private:
	int Learn(MedFeatures& features, unordered_set<int>& ids);
	int Apply(MedFeatures& features, unordered_set<int>& ids);
};


//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(FeatureProcessor)
MEDSERIALIZE_SUPPORT(MultiFeatureProcessor)
MEDSERIALIZE_SUPPORT(FeatureBasicOutlierCleaner)
MEDSERIALIZE_SUPPORT(FeatureNormalizer)
MEDSERIALIZE_SUPPORT(featureStrata)
MEDSERIALIZE_SUPPORT(featureSetStrata)
MEDSERIALIZE_SUPPORT(FeatureImputer)
MEDSERIALIZE_SUPPORT(FeatureSelector)
MEDSERIALIZE_SUPPORT(LassoSelector)
MEDSERIALIZE_SUPPORT(DgnrtFeatureRemvoer)
MEDSERIALIZE_SUPPORT(UnivariateFeatureSelector)
MEDSERIALIZE_SUPPORT(MRMRFeatureSelector)
MEDSERIALIZE_SUPPORT(FeatureEncoder)
MEDSERIALIZE_SUPPORT(FeaturePCAParams)
MEDSERIALIZE_SUPPORT(FeaturePCA)
MEDSERIALIZE_SUPPORT(TagFeatureSelector)
MEDSERIALIZE_SUPPORT(ImportanceFeatureSelector)
MEDSERIALIZE_SUPPORT(OneHotFeatProcessor)
MEDSERIALIZE_SUPPORT(univariateSelectionParams)
#endif
