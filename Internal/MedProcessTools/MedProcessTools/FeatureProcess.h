#ifndef _FTR_PROCESS_H_
#define _FTR_PROCESS_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedProcessTools/MedProcessTools/MedValueCleaner.h"
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
	FTR_PROCESS_LAST
} FeatureProcessorTypes;

class FeatureProcessor : public SerializableObject  {
public:

	// Type
	FeatureProcessorTypes processor_type;

	// Threading
	int learn_nthreads, clean_nthreads;

	// Constructor/Destructor
	FeatureProcessor() { learn_nthreads = DEFAULT_FEAT_CLNR_NTHREADS;  clean_nthreads = DEFAULT_FEAT_CLNR_NTHREADS; };
	~FeatureProcessor() {};

	// Virtual Set Feature Name
	virtual void set_name(const string& name) { return; }

	// Learn cleaning model
	virtual int Learn(MedFeatures& features, unordered_set<int>& ids) { return 0; }

	int learn(MedFeatures& features);
	int learn(MedFeatures& features, unordered_set<int>& ids) { return Learn(features, ids); }

	// Apply cleaning model
	virtual int Apply(MedFeatures& features, unordered_set<int>& ids) { return 0; }

	int apply(MedFeatures& features) ;
	int apply(MedFeatures& features, unordered_set<int>& ids) { return Apply(features, ids); }

	// Init
	static FeatureProcessor *make_processor(string name);
	static FeatureProcessor *make_processor(FeatureProcessorTypes type);
	static FeatureProcessor *make_processor(string name, string params);
	static FeatureProcessor *make_processor(FeatureProcessorTypes type, string params);

	virtual int init(void *processor_params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

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
	// Cleaners
	vector<FeatureProcessor *> processors;

	// Constructor/Destructor
	MultiFeatureProcessor() { processor_type = FTR_PROCESS_MULTI; };
	~MultiFeatureProcessor() {};

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);


	// Add processors
	void add_processors_set(FeatureProcessorTypes type, vector<string>& features);
	void add_processors_set(FeatureProcessorTypes type, vector<string>& features, string init_string);

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

	// Name
	string feature_name;

	// Constructor
	FeatureBasicOutlierCleaner() : FeatureProcessor() { init_defaults(); }
	FeatureBasicOutlierCleaner(string& name) : FeatureProcessor() { feature_name = name;  init_defaults(); }
	FeatureBasicOutlierCleaner(string& name, string init_string) : FeatureProcessor() { feature_name = name;  init_defaults();  init_from_string(init_string); }
	FeatureBasicOutlierCleaner(string& name, ValueCleanerParams *_params) : FeatureProcessor() { feature_name = name;  MedValueCleaner::init(_params); }

	void init_defaults() {
		processor_type = FTR_PROCESS_BASIC_OUTLIER_CLEANER;
		params.missing_value = MED_MAT_MISSING_VALUE; 
		params.trimming_sd_num = DEF_FTR_TRIMMING_SD_NUM; 
		params.removing_sd_num = DEF_FTR_REMOVING_SD_NUM;
		params.take_log = 0; 
		params.doTrim = params.doRemove = true;
		params.type = VAL_CLNR_ITERATIVE;

	};

	// Set Feature Name
	void set_name(const string& name) { feature_name = name; }

	// Init
	int init(void *processor_params) {return MedValueCleaner::init(processor_params);};
	int init(map<string, string>& mapper) { init_defaults();  return MedValueCleaner::init(mapper); };

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

	// Name
	string feature_name;

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
	FeatureNormalizer(const  string& name) : FeatureProcessor() { init_defaults(); feature_name = name; }
	FeatureNormalizer(const  string& name, string init_string) : FeatureProcessor() { init_from_string(init_string);  feature_name = name; }

	// Set Feature Name
	void set_name(const string& name) { feature_name = name; }

	// Learn cleaning model
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Apply cleaning model
	int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Init
	int init(map<string, string>& mapper) ;
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; normalizeSd = false; fillMissing = false; processor_type = FTR_PROCESS_NORMALIZER;};

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

	// Name
	string feature_name;

	// Missing Value
	float missing_value;

	// Strata for setting moment
	featureSetStrata imputerStrata;

	// Moment
	imputeMomentTypes moment_type;
	vector<float> moments;

	// Constructor
	FeatureImputer() : FeatureProcessor() { init_defaults(); }
	FeatureImputer(const  string& name) : FeatureProcessor() { init_defaults(); feature_name = name; }
	FeatureImputer(const  string& name, string init_string) : FeatureProcessor() { init_from_string(init_string);  feature_name = name; }

	// Set Feature Name
	void set_name(const string& name) { feature_name = name; }

	// Add stratirfier
	void addStrata(string& init_string);
	void addStrata(featureStrata& strata) { imputerStrata.stratas.push_back(strata); }
	void addStrata(string& name, float resolution, float min, float max) { imputerStrata.stratas.push_back(featureStrata(name, resolution, min, max)); }

	// Init
	int init(map<string, string>& mapper);
	void init_defaults() { missing_value = MED_MAT_MISSING_VALUE; moment_type = IMPUTE_MMNT_MEAN;  processor_type = FTR_PROCESS_IMPUTER; };

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
// Utilities
//.......................................................................................
//.......................................................................................

void get_all_values(MedFeatures& features, string& signalName, unordered_set<int>& ids, vector<float>& values);


#endif
