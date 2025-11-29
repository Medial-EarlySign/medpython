#ifndef __LABEL_PARAMS_H__
#define __LABEL_PARAMS_H__

#include "MedEnums.h"
#include <unordered_map>
#include <map>
#include <Logger/Logger/Logger.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <array>

using namespace std;
/**
* A warpper class for initializing rules for time window interaction
*/
class TimeWindowInteraction {
private:
	bool has_default_mode; ///< If true has default value
	bool has_default_range; ///< If true has default value range
	TimeWindowMode default_modes[2]; ///< the deafult value
	pair<float, float> default_intersection_range; ///< the default range

public:
	unordered_map<float, std::array<TimeWindowMode,2>> interaction_mode; ///< map from label value to interaction rules
	unordered_map<float, pair<float, float>> intersection_range_condition; ///< intersection rate range with registry condition

	TimeWindowInteraction() {
		has_default_mode = false;
		has_default_range = false;
		default_modes[0] = TimeWindowMode::All_;
		default_modes[1] = TimeWindowMode::All_;
		default_intersection_range.first = 0; default_intersection_range.second = 0;
	}

	TimeWindowMode *operator[] (float x) {
		if (interaction_mode.find(x) != interaction_mode.end() || !has_default_mode)
			return interaction_mode[x].data();
		//has_default_mode is true and not exist
		return default_modes;
	}

	bool find(const float x) const {
		return has_default_mode || interaction_mode.find(x) != interaction_mode.end();
	}

	const TimeWindowMode *at(float x) const {
		if (interaction_mode.find(x) != interaction_mode.end() || !has_default_mode)
			return interaction_mode.at(x).data();
		return default_modes;
	}

	void set_default(TimeWindowMode defaults_m[2]) {
		if (has_default_mode)
			HMTHROW_AND_ERR("Error - TimeWindowInteraction has already default\n");
		has_default_mode = true;
		default_modes[0] = defaults_m[0];
		default_modes[1] = defaults_m[1];
	}

	void set_default_range(float min_range, float max_range) {
		default_intersection_range.first = min_range;
		default_intersection_range.second = max_range;
		has_default_range = true;
	}

	bool get_inresection_range_cond(float x, float &min_range, float &max_range) const {
		if (intersection_range_condition.find(x) != intersection_range_condition.end()) {
			min_range = intersection_range_condition.at(x).first;
			max_range = intersection_range_condition.at(x).second;
			return true;
		}
		else if (has_default_range) {
			min_range = default_intersection_range.first;
			max_range = default_intersection_range.second;
			return true;
		}
		return false;
	}

	void reset_for_init() {
		has_default_mode = false;
		has_default_range = false;
		interaction_mode.clear();
		intersection_range_condition.clear();
	}

	void init_from_string(const string &init);
};

/**
* Parameters for lableing strategy on MedRegistry for given time window
*/
class LabelParams : public SerializableObject {
public:
	TimeWindowInteraction label_interaction_mode; ///< the label interaction definition
	TimeWindowInteraction censor_interaction_mode; ///< the label interaction definition
	int time_from; ///< time window from - for labeling
	int time_to; ///< time window to - for labeling
	ConflictMode conflict_method; ///< resolving conflicts
	int censor_time_from; ///< time window from - for censor registry
	int censor_time_to; ///< time window to - for censor registry
	bool treat_0_class_as_other_classes; ///< used to define the way outcome time is set to 0 class. in case control should be false, inmulticlass true
	LabelParams();

	/// init function
	int init(map<string, string>& map);
};

#endif