#ifndef _FTR_PROCESS_H_
#define _FTR_PROCESS_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedProcessTools/MedProcessTools/MedValueCleaner.h"
#include "MedStat/MedStat/MedPerformance.h"
#include "MedStat/MedStat/MedStat.h"
#include <unordered_set>

#define DEFAULT_FEAT_CLNR_NTHREADS 24

//.......................................................................................
//.......................................................................................
// A virtual class of processes on MedFeatures;
// E.g. Cleaning
//.......................................................................................
//.......................................................................................

typedef enum {
	FTR_PROCESS_MULTI,
	FTR_PROCESS_BASIC_OUTLIER_CLEANER,
	FTR_PROCESS_NORMALIZER,
	FTR_PROCESS_IMPUTER,
	FTR_PROCESS_DO_CALC,
	FTR_PROCESS_UNIVARIATE_SELECTOR,
	FTR_PROCESSOR_MRMR_SELECTOR,
	FTR_PROCESS_LAST
} FeatureProcessorTypes;

class FeatureProcessor : public SerializableObject {
public:

	// Feature name ( + name as appears in MedFeatures) ;
	string feature_name;
	string resolved_feature_name;

	// Type
	FeatureProcessorTypes processor_type;

	// Threading
	int learn_nthreads, clean_nthreads;

	// Constructor/Destructor
	FeatureProcessor() { learn_nthreads = DEFAULT_FEAT_CLNR_NTHREADS;  clean_nthreads = DEFAULT_FEAT_CLNR_NTHREADS; };
	~FeatureProcessor() {};

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *processor;}

	// Virtual Set Feature Name
	virtual void set_feature_name(const string& feature_name) { this->feature_name = feature_name; }
	virtual string get_feature_name() { return this->feature_name; }
	virtual void get_feature_names(vector<string> & feature_names) { feature_names.clear(); feature_names.push_back(feature_name); };

	// Learn cleaning model
	virtual int Learn(MedFeatures& features, unordered_set<int>& ids) { return 0; }

	int learn(MedFeatures& features);
	int learn(MedFeatures& features, unordered_set<int>& ids) { return Learn(features, ids); }

	// Apply cleaning model
	virtual int Apply(MedFeatures& features, unordered_set<int>& ids) { return 0; }

	int apply(MedFeatures& features) ;
	int apply(MedFeatures& features, unordered_set<int>& ids) { return Apply(features, ids); }

	// Init
	static FeatureProcessor *make_processor(string processor_name);
	static FeatureProcessor *make_processor(FeatureProcessorTypes type);
	static FeatureProcessor *make_processor(string processor_name, string params);
	static FeatureProcessor *make_processor(FeatureProcessorTypes type, string params);

	virtual int init(void *processor_params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	// Filter according to a subset of features
	virtual int filter(unordered_set<string>& features) { return (features.find(feature_name) == features.end()) ? 0 : 1; };

	// Utility : get corresponding name in MedFeatures
	string resolve_feature_name(MedFeatures& features, string substr);

	// Serialization (including type)
	size_t get_processor_size();
	size_t processor_serialize(unsigned char *blob);
};

// Utilities
FeatureProcessorTypes feature_processor_name_to_type(const string& cleaner_name);

//.......................................................................................
//.......................................................................................
// A Processor which contains a vector of simpler processors
// Useful for applying same cleaners on a set of features, for example
//.......................................................................................
//.......................................................................................

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
	MultiFeatureProcessor() { processor_type = FTR_PROCESS_MULTI; duplicate = 0; };
	~MultiFeatureProcessor() {};

	// Init
	int init(map<string, string>& mapper);

	// Copy
	virtual void copy(FeatureProcessor *processor) ;

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

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

//.......................................................................................
//.......................................................................................
// A simple cleaner considering each value of a certain feature separatley
//.......................................................................................
//.......................................................................................

#define DEF_FTR_TRIMMING_SD_NUM 7
#define DEF_FTR_REMOVING_SD_NUM 14

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
	int init(void *processor_params) {return MedValueCleaner::init(processor_params);};
	int init(map<string, string>& mapper) { init_defaults();  return MedValueCleaner::init(mapper); };

	// Copy
	virtual void copy(FeatureProcessor *processor) {*this = *(dynamic_cast<FeatureBasicOutlierCleaner *>(processor));}

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);
	int iterativeLearn(MedFeatures& features, unordered_set<int>& ids);
	int quantileLearn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

};

//.......................................................................................
//.......................................................................................
// Feature Normalizer
//.......................................................................................
//.......................................................................................

class FeatureNormalizer : public FeatureProcessor {
public:

	// Missing Value
	float missing_value;

	// Normalize Standard Deviation
	bool normalizeSd;

	// Fill missing values with mean 
	bool fillMissing;

	// Moments
	float mean, sd;

	// Constructor
	FeatureNormalizer() : FeatureProcessor() { init_defaults(); }
	FeatureNormalizer(const  string& feature_name) : FeatureProcessor() { init_defaults(); set_feature_name(feature_name); }
	FeatureNormalizer(const  string& feature_name, string init_string) : FeatureProcessor() { init_from_string(init_string);  set_feature_name(feature_name); }

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Init
	int init(map<string, string>& mapper) ;
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; normalizeSd = true; fillMissing = false; processor_type = FTR_PROCESS_NORMALIZER;};

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<FeatureNormalizer *>(processor)); }

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

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
	IMPUTE_MMNT_LAST
} imputeMomentTypes;

class featureStrata {
public:
	string name;
	float resolution, min, max;
	int nValues;

	featureStrata() {};
	featureStrata(string& _name, float _resolution, float _min, float _max) { name = _name; resolution = _resolution; min = _min; max = _max; }

	void SetNValues() { nValues = ((int)(max / resolution) - (int)(min / resolution) + 1);}
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
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

} ;

class featureSetStrata {
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
		for (int i = 2; i < stratas.size(); i++)
			factors[i] = factors[i - 1] * stratas[i - 1].nValues;
	}

	int nValues() { 
		if (stratas.size() == 0)
			return 1;
		else
			return factors.back() * stratas.back().nValues; 
	}

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

class FeatureImputer : public FeatureProcessor {
public:

	// Missing Value
	float missing_value;

	// Strata for setting moment
	featureSetStrata imputerStrata;

	// Moment
	imputeMomentTypes moment_type;
	vector<float> moments;

	// Constructor
	FeatureImputer() : FeatureProcessor() { init_defaults(); }
	FeatureImputer(const  string& feature_name) : FeatureProcessor() { init_defaults(); set_feature_name(feature_name); }
	FeatureImputer(const  string& feature_name, string init_string) : FeatureProcessor() { init_from_string(init_string);  set_feature_name(feature_name);}

	// Add stratirfier
	void addStrata(string& init_string);
	void addStrata(featureStrata& strata) { imputerStrata.stratas.push_back(strata); }
	void addStrata(string& name, float resolution, float min, float max) { imputerStrata.stratas.push_back(featureStrata(name, resolution, min, max)); }

	// Init
	int init(map<string, string>& mapper);
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; moment_type = IMPUTE_MMNT_MEAN;  processor_type = FTR_PROCESS_IMPUTER; };
	imputeMomentTypes getMomentType(string& entry);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<FeatureImputer *>(processor)); }

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// debug and print
	void print();

};

//.......................................................................................
//.......................................................................................
// Feature Selector
//.......................................................................................
//.......................................................................................

class FeatureSelector : public FeatureProcessor {
public:

	// Missing Value
	float missing_value;

	// Reauired Features
	unordered_set<string> required;

	// Selected Features (ordered)
	vector<string> selected;

	// Target number to select (if 0, ignored)
	int numToSelect;

	// Constructor
	FeatureSelector() : FeatureProcessor() {}

	// Find set of selected features
	virtual int Learn(MedFeatures& features, unordered_set<int>& ids);
	virtual int _learn(MedFeatures& features, unordered_set<int>& ids) { return 0; }

	// Apply selection
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

//.......................................................................................
//.......................................................................................
// Feature Selector : Univariate
//.......................................................................................
//.......................................................................................

typedef enum {
	UNIV_SLCT_PRSN = 0,
	UNIV_SLCT_MI = 1,
	UNIV_SLCT_DCORR = 2,
	UNIV_SLCT_LAST
} UnivariateSelectionMethod;

class univariateSelectionParams {
public:
	UnivariateSelectionMethod method;
	float minStat;

	// for mutual information
	int nBins;
	MedBinningType binMethod = BIN_EQUIDIST;

	// for samples distance correlation
	float pDistance; 

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
};

class UnivariateFeatureSelector : public FeatureSelector {
public:

	// Selection Params
	univariateSelectionParams params;

	// Constructor
	UnivariateFeatureSelector() : FeatureSelector() { init_defaults();}

	// Find set of selected features
	virtual int _learn(MedFeatures& features, unordered_set<int>& ids);

	// Init
	int init(map<string, string>& mapper);
	virtual void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; processor_type = FTR_PROCESS_UNIVARIATE_SELECTOR;  params.method = UNIV_SLCT_PRSN; numToSelect = 0; params.minStat = 0.05F;};

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<UnivariateFeatureSelector *>(processor)); }

	// Scores 
	int getAbsPearsonCorrs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats);
	int getMIs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats);
	int getDistCorrs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats);
};

//.......................................................................................
//.......................................................................................
// Feature Selector : MRMR
//.......................................................................................
//.......................................................................................

typedef enum {
	MRMR_MAX = 0,
	MRMR_MEAN = 1,
	MRMR_LAST
} MRMRPenaltyMethod;

class MRMRFeatureSelector : public FeatureSelector {
public:

	// Selection Params
	univariateSelectionParams params;
	float penalty;
	MRMRPenaltyMethod penaltyMethod;

	// Constructor
	MRMRFeatureSelector() : FeatureSelector() { init_defaults(); }

	// Find set of selected features
	virtual int _learn(MedFeatures& features, unordered_set<int>& ids);

	// Init
	int init(map<string, string>& mapper);
	virtual void init_defaults(); 
	MRMRPenaltyMethod get_penalty_method(string _method);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<MRMRFeatureSelector *>(processor)); }

	// Scores 
	int fillStatsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index);
	int fillAbsPearsonCorrsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index);
	int fillMIsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index);
	int fillDistCorrsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats,int index);
};

//.......................................................................................
//.......................................................................................
// Utilities
//.......................................................................................
//.......................................................................................

void get_all_values(MedFeatures& features, string& signalName, unordered_set<int>& ids, vector<float>& values);
void get_all_outcomes(MedFeatures& features, unordered_set<int>& ids, vector<float>& values);
void smearBins(vector<int>& bins, int nBins, int reqNbins);

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(MultiFeatureProcessor)
MEDSERIALIZE_SUPPORT(FeatureBasicOutlierCleaner)
MEDSERIALIZE_SUPPORT(FeatureNormalizer)
MEDSERIALIZE_SUPPORT(featureStrata)
MEDSERIALIZE_SUPPORT(featureSetStrata)
MEDSERIALIZE_SUPPORT(FeatureImputer)


#endif
