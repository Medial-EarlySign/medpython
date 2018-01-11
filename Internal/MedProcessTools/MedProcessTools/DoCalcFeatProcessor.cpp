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

// Set single-input calculation in/out names
void DoCalcFeatProcessor::set_feature_name(const string& feature_name) {

	if (raw_target_feature_name != "")
		MTHROW_AND_ERR("DoCalcFeatProcessor cannot set feature name for %s on %s . Already given: %s", calc_type.c_str(), feature_name.c_str(), raw_target_feature_name.c_str());

	raw_source_feature_names.assign(1,feature_name);
	this->feature_name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + calc_type + "_" + feature_name;
}

int DoCalcFeatProcessor::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [DoCalcFeatProcessor::init]
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
			for (string& s : vals)
				weights.push_back(stof(s));
		}
		else if (field != "fp_type")
			MLOG("Unknown parameter \'%s\' for DoCalcFeatProcessor\n", field.c_str());
		//! [DoCalcFeatProcessor::init]
	}
	if (weights.size() > 0 && weights.size() != raw_source_feature_names.size())
		MTHROW_AND_ERR("DoCalcFeatProcessor got [%d] weights != [%d] source_feature_names", (int)weights.size(), (int)raw_source_feature_names.size());

	// Default lists of source features : See examples nin Config_Exacmple
	if (raw_source_feature_names.empty()) {
		if (calc_type == "fragile") {
			raw_source_feature_names.push_back("Weight.last.win_180_360");
			raw_source_feature_names.push_back("Weight.min.win_0_180");
			raw_source_feature_names.push_back("BMI.min.win_0_180");
			raw_source_feature_names.push_back("BMI.max.win_0_180");
			raw_source_feature_names.push_back("CRP.max.win_0_180");
			raw_source_feature_names.push_back("Albumin.min.win_0_180");
			raw_source_feature_names.push_back("WBC.min.win_0_180");
			raw_source_feature_names.push_back("WBC.max.win_0_180");
			raw_source_feature_names.push_back("Gender");
			raw_source_feature_names.push_back("Hemoglobin.min.win_0_180");
			raw_source_feature_names.push_back("Current_Smoker");
			raw_source_feature_names.push_back("category_set_RC_Weight");
			raw_source_feature_names.push_back("category_set_RC_FallAdmission");
			raw_source_feature_names.push_back("category_set_RC_Weakness");
			raw_source_feature_names.push_back("category_set_RC_Vision");
			raw_source_feature_names.push_back("category_set_RC_Dyspnea");
			raw_source_feature_names.push_back("category_set_RC_Fatigue");
			raw_source_feature_names.push_back("category_set_RC_Chronic_Pain");
			raw_source_feature_names.push_back("category_set_RC_Urinery_Inconsistence");
			raw_source_feature_names.push_back("category_set_RC_Depression");
			raw_source_feature_names.push_back("category_set_RC_Coginition_Problems");
			raw_source_feature_names.push_back("category_set_RC_Social");
		}
		else if (calc_type == "min_chads2" || calc_type == "max_chads2") { 
			raw_source_feature_names = { "Age", "DM_Registry", "HT_Registry", "strokeIndicator", "chfIndicator" };
		}
		else if (calc_type == "min_chads2_vasc" || calc_type == "max_chads_vasc2") {
			raw_source_feature_names = { "Age", "DM_Registry", "HT_Registry", "strokeIndicator", "chfIndicator","Gender","vascIndicator" };
		}
		else if (calc_type == "framingham_chd") {
			raw_source_feature_names = { "Gender","Age", "DM_Registry", "Current_Smoker", "BP.last.win_0_1095", "BP.last.win_0_1095.t0v1","Cholesterol.last.win_0_1095","LDL.last.win_0_1095"};
		}

	}

	// Set name
	if (raw_target_feature_name == "") {
		string inFeatures = boost::join(raw_source_feature_names, "_");
		feature_name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + calc_type + "_" + inFeatures;
	}
	else if (raw_target_feature_name.substr(0, 4) == "FTR_")
		feature_name = raw_target_feature_name;
	else
		feature_name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + raw_target_feature_name;

	return 0;
}

int DoCalcFeatProcessor::Apply(MedFeatures& features, unordered_set<int>& ids) {
	
	// Get Source Features
	resolve_feature_names(features);

	// Prepare new Feature
	int samples_size = (int)features.samples.size();
	features.data[feature_name].clear();
	features.data[feature_name].resize(samples_size);
	float *p_out = &(features.data[feature_name][0]);

	// Attributes
	features.attributes[feature_name].normalized = false;
	features.attributes[feature_name].imputed = true;

	// Prepare
	vector<float*> p_sources;
	for (string source : source_feature_names) {
		assert(features.data.find(source) != features.data.end());
		p_sources.push_back(&(features.data[source][0]));
	}

	// Do your stuff
	if (calc_type == "sum")
		sum(p_sources, p_out, samples_size);
	else if (calc_type == "framingham_chd")
		framingham_chd(p_sources, p_out, samples_size);
	else if (calc_type == "min_chads2")
		chads2(p_sources, p_out, samples_size, 0, 0);
	else if (calc_type == "max_chads2")
		chads2(p_sources, p_out, samples_size, 0, 1);
	else if (calc_type == "min_chads2_vasc")
		chads2(p_sources, p_out, samples_size, 1, 0);
	else if (calc_type == "max_chads2_vasc")
		chads2(p_sources, p_out, samples_size, 1, 1);
	else if (calc_type == "fragile")
		fragile(p_sources, p_out, samples_size);	
	else if (calc_type == "log")
		_log(p_sources, p_out, samples_size);
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


void DoCalcFeatProcessor::framingham_chd(vector<float*> p_sources, float *p_out, int n_samples) {
	// www.framinghamheartstudy.org/risk-functions/cardiovascular-disease/10-year-risk.php
	for (int i = 0; i < n_samples; i++) {
		float res = 0.0;
		float res1 = 0.0;

		float gender = p_sources[0][i];
		float Age = p_sources[1][i];
		float DM_Registry = p_sources[2][i];
		float Current_Smoker = p_sources[3][i];
		float BP_dia = p_sources[4][i];
		float BP_sys = p_sources[5][i];
		float chol = p_sources[6][i];
		float hdl = p_sources[7][i];
		float BP_drug = p_sources[8][i];

		if (DM_Registry== missing_value || Current_Smoker == missing_value || BP_dia == missing_value || BP_sys == missing_value || chol == missing_value || hdl == missing_value)
			MTHROW_AND_ERR("CalcFeatGenerator framingham_chd found missing value need to use imputer");
		//Men

		if (1 == 1) {
			float sum_beta = 0.0;
			if (gender == 1) {
				sum_beta += log(Age)*3.06117;
				sum_beta += log(chol)*1.12370;
				sum_beta += log(hdl)*-0.93263;
				
				if (BP_drug==0) 
					sum_beta += log(BP_sys)*1.93303;
				else 
					sum_beta += log(BP_sys)*1.99881;

				sum_beta += Current_Smoker*0.65451;
				sum_beta += DM_Registry*0.57367;
				res = 1 - pow(0.88936, exp(sum_beta - 23.9802));

			}
			//Women
			else {
				sum_beta += log(Age)*2.32888;
				sum_beta += log(chol)*1.20904;
				sum_beta += log(hdl)*-0.70833;
				
				if (BP_drug == 0)
					sum_beta += log(BP_sys)*2.76157;
				else 
					sum_beta += log(BP_sys)*2.82263;

				sum_beta += Current_Smoker*0.52873;
				sum_beta += DM_Registry*0.69154;
				res = 1 - pow(0.95012, exp(sum_beta - 26.1931));
			}
		}


		p_out[i] = res;
	}
}





void DoCalcFeatProcessor::fragile(vector<float*> p_sources, float *p_out, int n_samples) {
	for (int i = 0; i < n_samples; i++) {
		float res = 0.0;
		//Weight
		res += ((p_sources[0][i] - p_sources[1][i] > 6.8) ||
			(p_sources[2][i] <= 18.5) || (p_sources[3][i] >= 30)
			|| (p_sources[11][i] > 0));
		//Fall Admission
		res += p_sources[12][i] > 0;
		//Weakness Admission
		res += p_sources[13][i] > 0;
		//Vision:
		res += p_sources[14][i] > 0;
		//Dyspnea:
		res += p_sources[15][i] > 0;
		//Fatigue:
		res += p_sources[16][i] > 0;
		//Chronic Pain:
		res += p_sources[17][i] > 0;
		//Urinery Inconsistence:
		res += p_sources[18][i] > 0;
		//Depression:
		res += p_sources[19][i] > 0;
		//Coginition Problems:
		res += p_sources[20][i] > 0;
		//Social:
		res += p_sources[21][i] > 0;
		//Smoking - current smoker:
		res += p_sources[10][i] > 0;
		//CRP:
		res += p_sources[4][i] >= 5;
		//Albumin:
		res += p_sources[5][i] <= 3.6;
		//WBC:
		res += p_sources[6][i] <= 3.2 || p_sources[7][i] >= 9.8;
		//Hemoglobin:
		int gen = (int) p_sources[8][i];
		if (gen == GENDER_MALE)
			res += p_sources[9][i] < 12;
		else
			res += p_sources[9][i] < 13.7;

		p_out[i] = res;
	}
}

// Out = Log(In)
void DoCalcFeatProcessor::_log(vector<float*> p_sources, float *p_out, int n_samples) {

	float *p = p_sources[0];
	for (int i = 0; i < n_samples; i++) {

		if (p[i] == missing_value || p[i] <= 0)
			p_out[i] = missing_value;
		else
			p_out[i] = log(p[i]);
	}

	return;
}

