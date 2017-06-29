#include "DoCalcFeatProcessor.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void DoCalcFeatProcessor::init_defaults() {
	processor_type = FTR_PROCESS_DO_CALC;
	missing_value = MED_MAT_MISSING_VALUE;
	calc_type = "calc_type_not_set";
}

void DoCalcFeatProcessor::resolve_feature_names(MedFeatures &features) {
	this->source_feature_names.clear();

	for (string name : raw_source_feature_names) {
		string real_feature_name = resolve_feature_name(features, name);
		this->source_feature_names.push_back(real_feature_name);
	}
}

int DoCalcFeatProcessor::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "name")
			raw_target_feature_name = entry.second;
		else if (field == "calc_type")
			calc_type = entry.second;
		else if (field == "source_feature_names")
			split(raw_source_feature_names, entry.second, boost::is_any_of(","));
		else if (field == "weights") {
			weights.clear();
			vector<string> vals;
			split(vals, entry.second, boost::is_any_of(","));
			for (string& s: vals)
				weights.push_back(stof(s));
		}
		else if (field != "fp_type")
			MLOG("Unknown parameter \'%s\' for DoCalcFeatProcessor\n", field.c_str());
	}
	if (weights.size() > 0 && weights.size() != raw_source_feature_names.size())
		MTHROW_AND_ERR("DoCalcFeatProcessor got [%d] weights != [%d] source_feature_names", (int)weights.size(), (int)raw_source_feature_names.size());
	
	set_feature_name(raw_target_feature_name);
	return 0;
}

int DoCalcFeatProcessor::Apply(MedFeatures& features, unordered_set<int>& ids) {
	
	// Prepare new Feature
	int samples_size = (int)features.samples.size();
	features.data[feature_name].clear();
	features.data[feature_name].resize(samples_size);
	float *p_out = &(features.data[feature_name][0]);
	
	// Attributes
	features.attributes[feature_name].normalized = false;
	features.attributes[feature_name].imputed = true;

	// Get Source Features
	resolve_feature_names(features);
	vector<float*> p_sources;
	for (string source : source_feature_names) {
		assert(features.data.find(source) != features.data.end());
		p_sources.push_back(&(features.data[source][0]));
	}
	

	// Do your stuff
	if (calc_type == "sum")
		sum(p_sources, p_out, samples_size);
	else if (calc_type == "min_chads2")
		chads2(p_sources, p_out, samples_size, 0);
	else if (calc_type == "max_chads2")
		chads2(p_sources, p_out, samples_size, 1);
	else
		MTHROW_AND_ERR("CalcFeatGenerator got an unknown calc_type: [%s]", calc_type.c_str());
	
	return 0;
}

// Out = Sum(In) or Sum(In*W)
void DoCalcFeatProcessor::sum(vector<float*> p_sources, float *p_out, int n_samples) {

	for (int i = 0; i < n_samples; i++) {
		float res = 0.0;

		int cnt = 0;
		for (float* p : p_sources) {
			if (p[i] == missing_value) {
				res = missing_value;
				break;
			}
			else if (weights.size() > 0)
				res += p[i] * weights[cnt++];
			else
				res += p[i];
		}
		p_out[i] = res;
	}

	return;
}

// Chads2 Scores: ASSUME order of given names is : age,Diabetes Registry, Hyper-Tenstion Registry, Stroke/TIA indicator, CHF indicator
void DoCalcFeatProcessor::chads2(vector<float*> p_sources, float *p_out, int n_samples, int max_flag) {
	
	for (int i = 0; i < n_samples; i++) {

		float chads2 = 0;
		// Age
		if (p_sources[0][i] >= 75)
			chads2++;

		// Diabetes
		if (p_sources[1][i] == 2)
			chads2++;
		else if (p_sources[1][i] == missing_value && max_flag)
			chads2++;

		// HyperTension
		if (p_sources[2][i] == 1)
			chads2++;
		else if (p_sources[2][i] == missing_value && max_flag)
			chads2++;

		// S2 : Prior Stroke or TIA or thromboembolism
		if (p_sources[3][i] > 0)
			chads2 += 2;
		else if (p_sources[3][i] == missing_value && max_flag)
			chads2 += 2;

		// CHF
		if (p_sources[4][i] > 0)
			chads2++;
		else if (p_sources[4][i] == missing_value && max_flag)
			chads2++;

		p_out[i] = chads2;
	}

	return;
		
}
