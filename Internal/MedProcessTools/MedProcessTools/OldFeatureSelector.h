#ifndef _FEAT_SELECTOR_H_
#define _FEAT_SELECTOR_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/FeatureGenerator.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"

#define DEFAULT_FEAT_SLCTR_NTHREADS 8


//.......................................................................................
//.......................................................................................
// A virtual class of general selector of features
//.......................................................................................
//.......................................................................................

// Define types of features selector
typedef enum {
	FTR_SLCTR_ALL,
	FTR_SLCTR_LAST
} FeatureSelectorTypes;

class FeatureSelector : public SerializableObject {
public:

	// Type
	FeatureSelectorTypes selector_type;

	// Threading
	int nthreads;

	// Constructor/Destructor
	FeatureSelector() { nthreads = DEFAULT_FEAT_SLCTR_NTHREADS; };
	~FeatureSelector() {};

	// Init
	static FeatureSelector *make_selector(string name);
	static FeatureSelector *make_selector(FeatureSelectorTypes type);
	static FeatureSelector *make_selector(string name, string params);
	static FeatureSelector *make_selector(FeatureSelectorTypes type, string params);

	virtual int init(void *cleaner_params) { return 0; };
	virtual int init(map<string, string>& mapper) { return 0; };
	virtual void init_defaults() {};

	// Filter
	virtual void select(vector<FeatureGenerator *> generators, vector<FeatureGenerator *> selected) {return;}
	virtual void select(vector<FeatureGenerator *> generators) { return; }

	// Serialization (including type)
	size_t get_selector_size();
	size_t selector_serialize(unsigned char *blob);
};

// Utilities
FeatureSelectorTypes feature_selector_name_to_type(const string& selector_name);

//.......................................................................................
//.......................................................................................
// DummyFeatsSelector : Take them all
//.......................................................................................
//.......................................................................................
class DummyFeatsSelector : public FeatureSelector {
public:
	void select(vector<FeatureGenerator *> generators, vector<FeatureGenerator *> selected) {selected = generators;}
	void select(vector<FeatureGenerator *> generators) { return; };
};

#endif
