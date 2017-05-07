#include "DoCalcFeatGenerator.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void DoCalcFeatGenerator::init_defaults() {
	generator_type = FTR_GEN_DO_CALC;
	missing_value = MED_MAT_MISSING_VALUE;
	calc_type = "calc_type_not_set";
};

void DoCalcFeatGenerator::init(MedFeatures &features) {
	FeatureGenerator::init(features);
	for (string substr : raw_source_feature_names) {
		string real_feature_name = "";
		for (auto candidate : features.attributes)
			if (candidate.first.find(substr) != string::npos) {
				if (real_feature_name != "")
					throw runtime_error(string("source_feature_name [") + substr + "] matches both [" +
						real_feature_name + "] and [" + candidate.first + "], can not generate [" + raw_target_feature_name + "]");
				real_feature_name = candidate.first;
			}
		if (real_feature_name == "") {
			string err = string("source_feature_name [") + substr + "] does not match any feature, can not generate [" +
				raw_target_feature_name + "]. Tried matching to these features:\n";
			for (auto candidate : features.attributes)
				err += candidate.first + "\n";
			throw runtime_error(err);
		}
		source_feature_names.push_back(real_feature_name);
	}
}

void DoCalcFeatGenerator::set_names() {
	names.clear();
	string name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + raw_target_feature_name;
	names.push_back(name);
}

int DoCalcFeatGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "target_feature_name")
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
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for DoCalcFeatGenerator\n", field.c_str());
	}
	if (weights.size() > 0 && weights.size() != raw_source_feature_names.size())
		MTHROW_AND_ERR("DoCalcFeatGenerator got [%d] weights != [%d] source_feature_names", (int)weights.size(), (int)raw_source_feature_names.size());
	set_names();
	return 0;
}

int DoCalcFeatGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	string& target = names[0];
	float *p_target = &(features.data[target][index]);
	vector<float*> p_sources;
	for (string source : source_feature_names) {
		assert(features.data.find(source) != features.data.end());
		p_sources.push_back(&(features.data[source][index]));
	}
	for (int i = 0; i < num; i++) {	
		if (calc_type == "sum")
			p_target[i] = sum(p_sources, i);
		else MTHROW_AND_ERR("CalcFeatGenerator got an unknown calc_type: [%s]", calc_type.c_str());
	}
	return 0;
}

float DoCalcFeatGenerator::sum(vector<float*> p_sources, int offset) {
	float res = 0.0;
	
	int cnt = 0;
	for (float* p : p_sources)
		if (p[offset] == missing_value)
			return missing_value;
		else if (weights.size() > 0)
			res += (*p) * weights[cnt++];
		else res += (*p);
		
	return res;
}
