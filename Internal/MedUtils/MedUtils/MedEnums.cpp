#include "MedEnums.h"
#include "Logger/Logger/Logger.h"
#include "MedUtils.h"

#define LOCAL_SECTION LOG_INFRA
#define LOCAL_LEVEL	LOG_DEF_LEVEL

vector<string> ConflictMode_to_name = { "all", "drop", "max", "last", "bitwise_max" };
vector<string> TimeWindow_to_name = { "before_end", "before_start" ,"after_start", "within", "all" };

ConflictMode ConflictMode_name_to_type(const string& ConflictMode_name) {
	for (int i = 0; i < ConflictMode_to_name.size(); ++i)
		if (ConflictMode_to_name[i] == ConflictMode_name) {
			return ConflictMode(i);
		}
	MTHROW_AND_ERR("Error in SamplingMode_name_to_type - Unsupported \"%s\". options are: %s\n",
		ConflictMode_name.c_str(), medial::io::get_list(ConflictMode_to_name).c_str());
}

TimeWindowMode TimeWindow_name_to_type(const string& TimeWindow_name) {
	for (int i = 0; i < TimeWindow_to_name.size(); ++i)
		if (TimeWindow_to_name[i] == TimeWindow_name) {
			return TimeWindowMode(i);
		}
	MTHROW_AND_ERR("Error in SamplingMode_name_to_type - Unsupported \"%s\". options are: %s\n",
		TimeWindow_name.c_str(), medial::io::get_list(TimeWindow_to_name).c_str());
}