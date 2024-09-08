#include "DoCalcFeatProcessor.h"
#include <omp.h>
#include <boost/lexical_cast.hpp>

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
		//! [DoCalcFeatProcessor::init]
		if (field == "name")
			raw_target_feature_name = entry.second;
		else if (field == "calc_type")
			calc_type = entry.second;
		else if (field == "source_feature_names")
			split(raw_source_feature_names, entry.second, boost::is_any_of(","));
		else if (field == "parameters")
			split(parameters, entry.second, boost::is_any_of(","));
		else if (field == "tags")
			split(tags, entry.second, boost::is_any_of(","));
		else if (field == "missing_value")
			missing_value = med_stof(entry.second);
		else if (field == "weights") {
			weights.clear();
			vector<string> vals;
			split(vals, entry.second, boost::is_any_of(","));
			for (string& s : vals)
				weights.push_back(stof(s));
		}
		else if (field != "fp_type")
			MLOG("Unknown parameter [%s] for DoCalcFeatProcessor\n", field.c_str());
		//! [DoCalcFeatProcessor::init]
	}
	if (weights.size() > 0 && weights.size() != raw_source_feature_names.size())
		MTHROW_AND_ERR("DoCalcFeatProcessor got [%d] weights != [%d] source_feature_names\n", (int)weights.size(), (int)raw_source_feature_names.size());

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
			raw_source_feature_names = { "Gender","Age", "DM_Registry", "Current_Smoker", "BP.last.win_0_1095", "BP.last.win_0_1095.t0v1",
				"Cholesterol.last.win_0_1095","HDL.last.win_0_1095", "Drug.category_set_hypertension_drugs.win_0_1095" };
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

void DoCalcFeatProcessor::prepare_feature(MedFeatures& features, int samples_size) const {
#pragma omp critical
	{
		features.data[feature_name].clear();
		features.data[feature_name].resize(samples_size);
		// Attributes
		features.attributes[feature_name].normalized = false;
		features.attributes[feature_name].imputed = true;

		for (const string& tag : tags)
			features.tags[feature_name].insert(tag);

	}
}

int DoCalcFeatProcessor::_apply(MedFeatures& features, unordered_set<int>& ids) {
	// Get Source Features
	resolve_feature_names(features);

	// Prepare new Feature
	int samples_size = (int)features.samples.size();
	prepare_feature(features, samples_size);

	float *p_out = &(features.data[feature_name][0]);

	// Prepare
	vector<float*> p_sources;
	for (const string &source : source_feature_names) {
		if (features.data.find(source) == features.data.end())
			MTHROW_AND_ERR("could not find feature [%s]\n", source.c_str());
		p_sources.push_back(&(features.data[source][0]));
	}

	// Do your stuff
	if (calc_type == "sum")
		sum(p_sources, p_out, samples_size);
	else if (calc_type == "max")
		max(p_sources, p_out, samples_size);
	else if (calc_type == "min")
		min(p_sources, p_out, samples_size);
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
	else if (calc_type == "threshold")
		do_threshold(p_sources, p_out, samples_size);
	else if (calc_type == "and")
		do_boolean_condition(p_sources, p_out, samples_size);
	else if (calc_type == "or")
		do_boolean_condition(p_sources, p_out, samples_size);
	else if (calc_type == "and_ignore_missing")
		do_boolean_condition_ignore_missing(p_sources, p_out, samples_size);
	else if (calc_type == "or_ignore_missing")
		do_boolean_condition_ignore_missing(p_sources, p_out, samples_size);
	else if (calc_type == "not")
		do_not(p_sources, p_out, samples_size);
	else
		MTHROW_AND_ERR("CalcFeatGenerator got an unknown calc_type: [%s]", calc_type.c_str());

	return 0;
}

/// check if a set of features is affected by the current processor
//.......................................................................................
bool DoCalcFeatProcessor::are_features_affected(unordered_set<string>& out_req_features) {

	// If empty = all features are required
	if (out_req_features.empty())
		return true;

	// Otherwise - check  generated features
	if (out_req_features.find(feature_name) != out_req_features.end())
		return true;

	return false;
}

/// update sets of required as input according to set required as output to processor
//.......................................................................................
void DoCalcFeatProcessor::update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features) {

	// If empty, keep as is
	if (out_req_features.empty())
		in_req_features.clear();
	else {
		in_req_features = out_req_features;
		// If active, add originals 
		if (are_features_affected(out_req_features)) {
			for (string& ftr : source_feature_names)
				in_req_features.insert(ftr);
		}
	}
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

void DoCalcFeatProcessor::min(vector<float*> p_sources, float *p_out, int n_samples) const {
	for (int i = 0; i < n_samples; i++) {
		float res = missing_value;
		for (float* p : p_sources)
			if (p[i] != missing_value && (res == missing_value || p[i] < res))
				res = p[i];
		p_out[i] = res;
	}
}

void DoCalcFeatProcessor::max(vector<float*> p_sources, float *p_out, int n_samples) const {
	for (int i = 0; i < n_samples; i++) {
		float res = missing_value;
		for (float* p : p_sources)
			if (p[i] != missing_value && (res == missing_value || p[i] > res))
				res = p[i];
		p_out[i] = res;
	}
}


void DoCalcFeatProcessor::do_threshold(vector<float*> p_sources, float *p_out, int n_samples) {
	//MLOG("DoCalcFeatProcessor::do_threshold start\n");
	if (p_sources.size() != 1)
		MTHROW_AND_ERR("do_threshold expects 1 source_feature_names, got [%d]\n", (int)p_sources.size());
	float *p = p_sources[0];
	if (parameters.size() != 2)
		MTHROW_AND_ERR("do_threshold expects 2 parameters, got [%d]\n", (int)parameters.size());
	if (parameters[0] != ">" && parameters[0] != "<" && parameters[0] != ">=" && parameters[0] != "<=")
		MTHROW_AND_ERR("do_threshold expects the first parameter to be one of [>,<,>=,<=], got [%s]\n", parameters[0].c_str());
	float thresold;
	try
	{
		thresold = boost::lexical_cast<float>(parameters[1]);
	}
	catch (const boost::bad_lexical_cast &) {
		MTHROW_AND_ERR("do_threshold expects the second parameter to be a float threshold, got [%s]\n", parameters[1].c_str());
	}

	for (int i = 0; i < n_samples; i++) {
		float res;
		if (p[i] == missing_value)
			res = missing_value;
		else if (parameters[0] == ">")
			res = p[i] > thresold;
		else if (parameters[0] == "<")
			res = p[i] < thresold;
		else if (parameters[0] == ">=")
			res = p[i] >= thresold;
		else if (parameters[0] == "<=")
			res = p[i] <= thresold;
		else MTHROW_AND_ERR("do_threshold expects the first parameter to be one of [>,<,>=,<=], got [%s]\n", parameters[0].c_str());
		p_out[i] = res;
	}
	//MLOG("DoCalcFeatProcessor::do_threshold end\n");
	return;
}

void DoCalcFeatProcessor::do_boolean_condition(vector<float*> p_sources, float *p_out, int n_samples) {
	//MLOG("DoCalcFeatProcessor::do_boolean_condition start\n");
	if (p_sources.size() < 1)
		MTHROW_AND_ERR("[%s] expects at least 1 source_feature_names, got [%d]\n", calc_type.c_str(), (int)p_sources.size());

	if (calc_type == "and") {
		for (int i = 0; i < n_samples; i++) {
			int res = 1;
			bool any_missing = false;
			for (int j = 0; j < p_sources.size(); j++) {
				if (p_sources[j][i] != missing_value) {
					res &= (p_sources[j][i] != 0.0);
				}
				else
					any_missing = true;
			}
			if (any_missing)
			{
				if (!res)
					p_out[i] = float(0);
				else
					p_out[i] = missing_value;
			}
			else
				p_out[i] = (float)res;
		}
	}
	else if (calc_type == "or") {
		for (int i = 0; i < n_samples; i++) {
			int res = 0;
			bool any_missing = false;
			for (int j = 0; j < p_sources.size(); j++) {
				if (p_sources[j][i] != missing_value) {
					res |= (p_sources[j][i] != 0.0);
				}
				else
					any_missing = true;
			}
			if (any_missing)
			{
				if (res)
					p_out[i] = float(1);
				else
					p_out[i] = missing_value;
			}
			else
				p_out[i] = (float)res;
		}
	}
	else MTHROW_AND_ERR("do_boolean_condition expects the first parameter to be one of [and,or], got [%s]\n", calc_type.c_str());
	return;
}

void DoCalcFeatProcessor::do_boolean_condition_ignore_missing(vector<float*> p_sources, float *p_out, int n_samples) {
	MLOG("DoCalcFeatProcessor::do_boolean_condition_ignore_missing start\n");
	if (p_sources.size() < 1)
		MTHROW_AND_ERR("[%s] expects at least 1 source_feature_names, got [%d]\n", calc_type.c_str(), (int)p_sources.size());

	if (calc_type == "and_ignore_missing") {
		for (int i = 0; i < n_samples; i++) {
			int res = 1;
			for (int j = 0; j < p_sources.size(); j++) {
				if (p_sources[j][i] != missing_value) {
					res &= (p_sources[j][i] != 0.0);
				}
			}

			p_out[i] = (float)res;
		}
	}
	else if (calc_type == "or_ignore_missing") {
		for (int i = 0; i < n_samples; i++) {
			int res = 0;
			for (int j = 0; j < p_sources.size(); j++) {
				if (p_sources[j][i] != missing_value) {
					res |= (p_sources[j][i] != 0.0);
				}
			}
			p_out[i] = (float)res;
		}
	}
	else MTHROW_AND_ERR("do_boolean_condition_ignore_missing expects the first parameter to be one of [and_ignore_missing,or_ignore_missing], got [%s]\n", calc_type.c_str());
	return;
}

void DoCalcFeatProcessor::do_not(vector<float*> p_sources, float *p_out, int n_samples) {
	MLOG("DoCalcFeatProcessor::do_not start\n");
	if (p_sources.size() != 1)
		MTHROW_AND_ERR("do_not expects 1 source_feature_names, got [%d]\n", (int)p_sources.size());
	float *p = p_sources[0];

	for (int i = 0; i < n_samples; i++) {
		if (p[i] == missing_value)
			p_out[i] = missing_value;
		else
			p_out[i] = (p[i] == 0.0);
	}
	MLOG("DoCalcFeatProcessor::do_not end\n");
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
		double res = 0.0;

		float gender = p_sources[0][i];
		float Age = p_sources[1][i];
		float DM_Registry = p_sources[2][i];
		float Current_Smoker = p_sources[3][i];
		float BP_dia = p_sources[4][i];
		float BP_sys = p_sources[5][i];
		float chol = p_sources[6][i];
		float hdl = p_sources[7][i];
		float BP_drug = p_sources[8][i];

		if (DM_Registry == missing_value || Current_Smoker == missing_value || BP_dia == missing_value || BP_sys == missing_value || chol == missing_value || hdl == missing_value)
			MTHROW_AND_ERR("CalcFeatGenerator framingham_chd found missing value need to use imputer");
		//Men

		double sum_beta = 0.0;
		if (gender == 1) {
			sum_beta += log(Age)*3.06117F;
			sum_beta += log(chol)*1.12370F;
			sum_beta += log(hdl)*-0.93263F;

			if (BP_drug == 0)
				sum_beta += log(BP_sys)*1.93303F;
			else
				sum_beta += log(BP_sys)*1.99881F;

			sum_beta += Current_Smoker * 0.65451F;
			sum_beta += DM_Registry * 0.57367F;
			res = 1 - pow(0.88936F, exp(sum_beta - 23.9802F));

		}
		//Women
		else {
			sum_beta += log(Age)*2.32888F;
			sum_beta += log(chol)*1.20904F;
			sum_beta += log(hdl)*-0.70833F;

			if (BP_drug == 0)
				sum_beta += log(BP_sys)*2.76157F;
			else
				sum_beta += log(BP_sys)*2.82263F;

			sum_beta += Current_Smoker * 0.52873F;
			sum_beta += DM_Registry * 0.69154F;
			res = 1 - pow(0.95012F, exp(sum_beta - 26.1931F));
		}



		p_out[i] = (float)res;
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
		int gen = (int)p_sources[8][i];
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

