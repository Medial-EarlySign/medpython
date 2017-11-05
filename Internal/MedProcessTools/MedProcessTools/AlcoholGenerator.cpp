#include "AlcoholGenerator.h"
#include <boost/algorithm/string/predicate.hpp>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace boost;

void generateAlcoholRangeSignal(SDateVal2* rawSignal, SDateRangeVal *outRangeSignal) {

}

void AlcoholGenerator::set_names() {
	names.clear();
	unordered_set<string> legal_features({ "Current_Drinker", "Drinking_Quantity", "PLM_Drinking_Level" });
	if (raw_feature_names.size() == 0)
		MTHROW_AND_ERR("AlcoholGenerator got no alcohol_features");
	for (string s : raw_feature_names) {
		if (legal_features.find(s) == legal_features.end())
			MTHROW_AND_ERR("AlcoholGenerator does not know how to generate [%s]", s.c_str());
		names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);
	}
}

int AlcoholGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "alcohol_features")
			boost::split(raw_feature_names, entry.second, boost::is_any_of(","));
		else if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for AlcoholGenerator\n", field.c_str());
	}
	set_names();
	return 0;
}

int AlcoholGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	int missing = MED_MAT_MISSING_VALUE;
	
	for (int i = 0; i < num; i++) {
		int len;
		string sname = "ALCOHOL";
		SDateShort2 *alcohol_status = (SDateShort2 *)rec.get(sname, 0, len);
		int current_drinker, drinking_level, plm_drinking_level;

		if (len == 0) { // No Data
			current_drinker = drinking_level = plm_drinking_level = missing;
		}
		else {
			/*
PLM_Drinking_Level	per day	per week
0	not at all	`0
1	less than 1 	`1-4 drinks
2	1-2 drinks	`5-17 drinks
3	3-6 drinks	`18-45 drinks
4	7-9 drinks	`46-65 drinks
5	9 or more	`66-99 drinks
			*/
			current_drinker = alcohol_status[len - 1].val1;
			drinking_level = alcohol_status[len - 1].val2;
			if (current_drinker == 1 && drinking_level == -1)
				drinking_level = missing, plm_drinking_level = 2;
			else if (current_drinker == 0 || drinking_level <= 0)
				plm_drinking_level = 0;
			else if (drinking_level >= 1 && drinking_level <= 4)
				plm_drinking_level = 1;
			else if (drinking_level >= 5 && drinking_level <= 17)
				plm_drinking_level = 2;
			else if (drinking_level >= 18 && drinking_level <= 45)
				plm_drinking_level = 3;
			else if (drinking_level >= 46 && drinking_level <= 65)
				plm_drinking_level = 4;
			else 
				plm_drinking_level = 5;
		}

		for (int j = 0; j < names.size(); j++) {
			if (algorithm::ends_with(names[j], "Current_Drinker"))
				features.data[names[j]][index + i] = (float)current_drinker;
			else if (algorithm::ends_with(names[j], "Drinking_Quantity"))
				features.data[names[j]][index + i] = (float)drinking_level;
			else if (algorithm::ends_with(names[j], "PLM_Drinking_Level"))
				features.data[names[j]][index + i] = (float)plm_drinking_level;
			else MTHROW_AND_ERR("unknown feature name [%s]", names[j].c_str());
		}
	}

	return 0;
}
