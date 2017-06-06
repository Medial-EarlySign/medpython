#pragma once
#include "FeatureProcess.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"

//.......................................................................................
//.......................................................................................
// User defined calculations on other features.
// NOTE: it is the user responsibility to put these generators after their source features generators
//.......................................................................................
//.......................................................................................

class DoCalcFeatGenerator : public FeatureProcessor {
public:
	int serial_id;
	// target_feature_name as specified by the user, will be decorated for uniqueness
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

	DoCalcFeatGenerator() : FeatureProcessor() { serial_id = ++MedFeatures::global_serial_id_cnt; init_defaults(); }
	~DoCalcFeatGenerator() {};

	virtual void set_feature_name(const string& feature_name) { 
		if (feature_name.substr(0, 4) == "FTR_") 
			this->feature_name = feature_name; // already uniq
		else this->feature_name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + feature_name;
	}

	void init_defaults();

	// init_from_string
	int init(map<string, string>& mapper);

	virtual int Apply(MedFeatures& features, unordered_set<int>& ids);

	float sum(vector<float*> p_sources, int offset);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<DoCalcFeatGenerator *>(processor)); }

	// Serialization
	size_t get_size() { return MedSerialize::get_size(processor_type, serial_id, raw_target_feature_name, feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, serial_id, processor_type, raw_target_feature_name, feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights); }
	size_t deserialize(unsigned char *blob) { return MedSerialize::deserialize(blob, serial_id, processor_type, raw_target_feature_name, feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights); }

private:
	virtual void resolve_feature_names(MedFeatures &features);
};