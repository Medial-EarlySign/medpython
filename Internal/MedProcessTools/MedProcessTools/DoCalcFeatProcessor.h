#pragma once
#include "FeatureProcess.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"

//.......................................................................................
//.......................................................................................
// User defined calculations on other features.
// NOTE: it is the user responsibility to put these generators after their source features generators
//.......................................................................................
//.......................................................................................

class DoCalcFeatProcessor : public FeatureProcessor {
public:
	int serial_id;
	// target_feature_name as specified by the user, will be decorated for uniqueness and extra information
	string raw_target_feature_name = "";

	// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_source_feature_names;

	vector<string> source_feature_names;

	// user function selector (e.g. sum, ratio)
	string calc_type;

	// when a source_feature == missing_value, the calculation would also be missing_value
	float missing_value;

	// for sum
	vector<float> weights;

	// Functions
	DoCalcFeatProcessor() : FeatureProcessor() { serial_id = ++MedFeatures::global_serial_id_cnt; init_defaults(); }
	~DoCalcFeatProcessor() {};

	void init_defaults();

	// init_from_string
	int init(map<string, string>& mapper);

	virtual int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Specific Functions
	void sum(vector<float*> p_sources, float *p_out, int n_samples);
	void chads2(vector<float*> p_sources, float *p_out, int n_samples, int vasc_flag, int max_flag);
	void has_bled(vector<float*> p_sources, float *p_out, int n_samples, int max_flag);
	void fragile(vector<float*> p_sources, float *p_out, int n_samples);

	// Single Input Functions
	void _log(vector<float*> p_sources, float *p_out, int n_samples);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<DoCalcFeatProcessor *>(processor)); }

	// Serialization
	ADD_SERIALIZATION_FUNCS(processor_type, serial_id, raw_target_feature_name, feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights)

	// Naming
	void set_feature_name(const string& feature_name);

private:
	virtual void resolve_feature_names(MedFeatures &features);
};