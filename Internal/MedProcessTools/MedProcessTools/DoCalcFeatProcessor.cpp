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
		chads2(p_sources, p_out, samples_size, 0, 0);
	else if (calc_type == "max_chads2")
		chads2(p_sources, p_out, samples_size, 0, 1);
	else if (calc_type == "min_chads2_vasc")
		chads2(p_sources, p_out, samples_size, 1, 0);
	else if (calc_type == "max_chads2_vasc")
		chads2(p_sources, p_out, samples_size, 1, 1);
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

// Chads2 Scores: ASSUME order of given data is : age,Diabetes Registry, Hyper-Tenstion Registry, Stroke/TIA indicator, CHF indicator, (optinal: Sex, Vasc indicator)
void DoCalcFeatProcessor::chads2(vector<float*> p_sources, float *p_out, int n_samples, int vasc_flag, int max_flag) {
	
	for (int i = 0; i < n_samples; i++) {

		float chads2 = 0;
		// Age
		if (p_sources[0][i] >= 75) {
			if (vasc_flag == 1)
				chads2 += 2;
			else
				chads2++;
		}
		else if (p_sources[0][i] >= 65 && vasc_flag)
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

		if (vasc_flag) {
			// Sex
			if (p_sources[5][i] == 2)
				chads2++;

			// Vasc
			if (p_sources[6][i] > 0)
				chads2++;
			else if (p_sources[6][i] == missing_value && max_flag)
				chads2++;
		}
		

		p_out[i] = chads2;
	}

	return;
		
}

// HAS-BLED Scores: ASSUME order of given data is : systolic-blood-pressure,dialysis-info,transplant-info,creatinine,cirrhosis-info,bilirubin,AST,ALP,ALKP,Stroke indicator,Past bleeding indicator,
//													unstable-INR indicator ,age,drugs-indicator, alcohol/drug consumption
void DoCalcFeatProcessor::has_bled(vector<float*> p_sources, float *p_out, int n_samples, int max_flag) {

	for (int i = 0; i < n_samples; i++) {

		float score = 0;
		// Blood-pressure (0 : Systolic blood pressure)
		if (p_sources[0][i] == missing_value && max_flag)
			score++;
		if (p_sources[0][i] > 160)
			score++;

		// Kidney (1: Dialysis, 2: Transplanct, 3: Creatinine)
		if (p_sources[3][i] == missing_value && max_flag)
			score++;
		else if (p_sources[1][i] == 1 || p_sources[2][i] == 1 || p_sources[3][i] > 2.26)
			score++;

		// Liver (4: Cirrhosis, 5: Bilirubin, 6: AST, 7: ALT, 8: ALKP)
		if ((p_sources[5][i] == missing_value || p_sources[6][i] == missing_value || p_sources[7][i] == missing_value || p_sources[8][i] == missing_value) && max_flag)
			score++;
		else if (p_sources[4][i] == 1 || p_sources[5][i] > 2 * 1.9 || p_sources[6][i] > 3 * 40 || p_sources[7][i] > 3 * 56 || p_sources[8][i] > 3 * 147)
			score++;

		//Stroke History (9 : Indicator) 
		if (p_sources[9][i] == 1)
			score++;

		// Prior Major Bleeding event (10: Indicator)
		if (p_sources[10][i] == 1)
			score++;

		// Unstable INR (11: Indicator)
		if (p_sources[11][i] == missing_value && max_flag)
			score++;
		if (p_sources[11][i] == 1)
			score++;

		// Age (12)
		if (p_sources[12][i] > 65)
			score++;

		// Drugs (13 : Anti-platelets indicator, 14: NSAID indicator)
		if (p_sources[13][i] == 1 || p_sources[14][i] == 1)
			score++;

		// Alcohol (14 : Amount)
		if (p_sources[14][i] == missing_value && max_flag)
			score++;
		else if (p_sources[14][i] > 8)
			score++;

		p_out[i] = score;
	}

	return;

}