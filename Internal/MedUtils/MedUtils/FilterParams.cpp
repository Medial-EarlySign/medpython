#include "FilterParams.h"
#include <Logger/Logger/Logger.h>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL LOG_DEF_LEVEL

FilterParams::FilterParams() {
	min_age = 0;
	max_age = -1;
	min_time = 0;
	max_time = -1;
	interaction_mode = TimeWindowMode::Within;
}

int FilterParams::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "min_age")
			min_age = med_stoi(it->second);
		else if (it->first == "max_age")
			max_age = med_stoi(it->second);
		else if (it->first == "min_time")
			min_time = med_stoi(it->second);
		else if (it->first == "max_time")
			max_time = med_stoi(it->second);
		else if (it->first == "interaction_mode")
			interaction_mode = TimeWindow_name_to_type(it->second);
		else
			MTHROW_AND_ERR("Error in FilterParams::init - Unsupported param %s\n", it->first.c_str())
	}
	return 0;
}