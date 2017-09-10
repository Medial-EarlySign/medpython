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

	DoCalcFeatProcessor() : FeatureProcessor() { serial_id = ++MedFeatures::global_serial_id_cnt; init_defaults(); }
	~DoCalcFeatProcessor() {};

	virtual void set_feature_name(const string& feature_name) { 
		if (feature_name.substr(0, 4) == "FTR_") 
			this->feature_name = feature_name; // already uniq
		else this->feature_name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + feature_name;
	}

	void init_defaults();

	// init_from_string
	int init(map<string, string>& mapper);

	virtual int Apply(MedFeatures& features, unordered_set<int>& ids);

	// Specific Functions
	void sum(vector<float*> p_sources, float *p_out, int n_samples);
	void chads2(vector<float*> p_sources, float *p_out, int n_samples, int vasc_flag, int max_flag);
	void has_bled(vector<float*> p_sources, float *p_out, int n_samples, int max_flag);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<DoCalcFeatProcessor *>(processor)); }

	// Serialization
	ADD_SERIALIZATION_FUNCS(processor_type, serial_id, raw_target_feature_name, feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights)

private:
	virtual void resolve_feature_names(MedFeatures &features);
};