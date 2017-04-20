#include "SmokingGenerator.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace boost;

void SmokingGenerator::set_names() {
	names.clear();
	unordered_set<string> legal_features({ "Current_Smoker", "Ex_Smoker", "Smok_Years_Since_Quitting", "Smoking_Years", "Smok_Pack_Years" });
	if (raw_feature_names.size() == 0)
		MTHROW_AND_ERR("SmokingGenerator got no smoking_features");
	for (string s : raw_feature_names) {
		if (legal_features.find(s) == legal_features.end())
			MTHROW_AND_ERR("SmokingGenerator does not know how to generate [%s]",s.c_str());
		names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);
	}
}

int SmokingGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "smoking_features")
			split(raw_feature_names, entry.second, boost::is_any_of(","));
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for SmokingGenerator\n", field.c_str());
	}
	set_names();
	return 0;
}

int SmokingGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {
	int missing = -1;

	for (int i = 0; i < num; i++) {
		int current_smoker, ex_smoker;
		int years_since_quitting, smoking_years;
		float pack_years;

		int len;
		bool never_smoked = true;

		string sname = "SMOKING_ENRICHED";
		SValShort4 *smx_status = (SValShort4 *)rec.get(sname, 0, len);
		if (len > 0)
			never_smoked = (smx_status[0].val1 == -1);
		assert(len <= 1);

		if (len == 0) { // No Data
			current_smoker = ex_smoker = (int)missing;
			years_since_quitting = smoking_years = (int)missing;
			pack_years = (float)missing;
		}
		else if (never_smoked) { // Non Smoker
			current_smoker = ex_smoker = 0;
			years_since_quitting = 100;
			smoking_years = 0;
			pack_years = 0.0;
		}
		else { // (Ex)Smoker
			int start_year = smx_status[0].val1;
			int end_year = smx_status[0].val2;
			int target_year = (int)(med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time) / 10000);
			if (target_year < end_year) {
				// still in smoking period
				smoking_years = target_year - start_year;
				years_since_quitting = 0;
				current_smoker = 1;
			}
			else {
				// maybe done smoking
				current_smoker = smx_status[0].val4;
				smoking_years = end_year - start_year; // we are merciful
				if (!current_smoker)
					years_since_quitting = target_year - end_year;
				else
					years_since_quitting = 0;
			}
			pack_years = ((float)smx_status[0].val3 / 20) * smoking_years;
			ex_smoker = 1 - current_smoker;
		}
		for (int j = 0; j < raw_feature_names.size(); j++) {			
			if (raw_feature_names[j] == "Current_Smoker")
				features.data[names[j]][index + i] = (float)current_smoker;
			else if (raw_feature_names[j] == "Ex_Smoker")
				features.data[names[j]][index + i] = (float)ex_smoker;
			else if (raw_feature_names[j] == "Smok_Years_Since_Quitting")
				features.data[names[j]][index + i] = (float)years_since_quitting;
			else if (raw_feature_names[j] == "Smoking_Years")
				features.data[names[j]][index + i] = (float)smoking_years;
			else if (raw_feature_names[j] == "Smok_Pack_Years")
				features.data[names[j]][index + i] = (float)pack_years;
			else MTHROW_AND_ERR("unknown feature name [%s]", raw_feature_names[j].c_str());
		}
	}

	return 0;
}
