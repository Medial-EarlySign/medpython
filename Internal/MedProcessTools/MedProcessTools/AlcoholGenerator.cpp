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
	for (string s : legal_features) {
		names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);
	}
}

int AlcoholGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "tags")
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
			current_drinker = alcohol_status[len - 1].val1;
			drinking_level = alcohol_status[len - 1].val2;
			if (drinking_level == -1)
				plm_drinking_level = drinking_level = missing;
			else if (current_drinker == 0 || drinking_level <= 0)
				plm_drinking_level = 0;
			else if (drinking_level <= 1)
				plm_drinking_level = 1;
			else if (drinking_level <= 2)
				plm_drinking_level = 2;
			else if (drinking_level <= 6)
				plm_drinking_level = 3;
			else if (drinking_level <= 9)
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
