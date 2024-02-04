#ifndef __DO_CALC_FEAT_PROCESSOR_H__
#define __DO_CALC_FEAT_PROCESSOR_H__

#include <MedProcessTools/MedProcessTools/FeatureProcess.h>
#include <MedProcessTools/MedProcessTools/MedFeatures.h>

/**
* User defined calculations on other features.\n
* NOTE: it is the user responsibility to put these generators after their source features generators\n
*/
class DoCalcFeatProcessor : public FeatureProcessor {
public:
	int serial_id;
	/// target_feature_name as specified by the user, will be decorated for uniqueness and extra information
	string raw_target_feature_name = "";

	/// source_feature_names as specified by the user, will be resolved to decorated names
	vector<string> raw_source_feature_names;

	vector<string> source_feature_names;

	/// general-purpose parameters which can be used by the calc
	vector<string> parameters;

	/// user function selector (e.g. sum, ratio)
	string calc_type;

	/// when a source_feature == missing_value, the calculation would also be missing_value
	float missing_value;

	/// for sum
	vector<float> weights;

	/// Tags - for defining labels or groups. may be used later for filtering for example
	vector<string> tags;

	// Functions
	DoCalcFeatProcessor() : FeatureProcessor() { serial_id = ++MedFeatures::global_serial_id_cnt; init_defaults(); }
	~DoCalcFeatProcessor() {};

	void init_defaults();

	/// The parsed fields from init command.
	/// @snippet DoCalcFeatProcessor.cpp DoCalcFeatProcessor::init
	int init(map<string, string>& mapper);

	int _apply(MedFeatures& features, unordered_set<int>& ids);

	// Specific Functions
	void sum(vector<float*> p_sources, float *p_out, int n_samples);
	void chads2(vector<float*> p_sources, float *p_out, int n_samples, int vasc_flag, int max_flag);
	void has_bled(vector<float*> p_sources, float *p_out, int n_samples, int max_flag);
	void fragile(vector<float*> p_sources, float *p_out, int n_samples);
	void framingham_chd(vector<float*> p_sources, float *p_out, int n_samples);
	void do_boolean_condition(vector<float*> p_sources, float *p_out, int n_samples);
	void do_boolean_condition_ignore_missing(vector<float*> p_sources, float *p_out, int n_samples);
	void do_not(vector<float*> p_sources, float *p_out, int n_samples);

	// Single Input Functions
	void _log(vector<float*> p_sources, float *p_out, int n_samples);
	void do_threshold(vector<float*> p_sources, float *p_out, int n_samples);

	/// check if a set of features is affected by the current processor
	bool are_features_affected(unordered_set<string>& out_req_features);

	/// update sets of required as input according to set required as output to processor
	void update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features);

	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<DoCalcFeatProcessor *>(processor)); }

	// Serialization
	ADD_CLASS_NAME(DoCalcFeatProcessor)
	ADD_SERIALIZATION_FUNCS(processor_type, serial_id, raw_target_feature_name, feature_name, calc_type, missing_value, raw_source_feature_names, source_feature_names, weights, parameters, tags)

private:
	virtual void resolve_feature_names(MedFeatures &features);
	void prepare_feature(MedFeatures& features, int samples_size) const;

	void max(vector<float*> p_sources, float *p_out, int n_samples) const;
	void min(vector<float*> p_sources, float *p_out, int n_samples) const;
};

MEDSERIALIZE_SUPPORT(DoCalcFeatProcessor);

#endif