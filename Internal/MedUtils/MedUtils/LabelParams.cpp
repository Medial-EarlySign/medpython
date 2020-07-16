#include "LabelParams.h"
#include <MedMat/MedMat/MedMatConstants.h>

#define LOCAL_SECTION LOG_INFRA

void TimeWindowInteraction::init_from_string(const string &init) {
	reset_for_init();
	vector<string> tokens;
	boost::split(tokens, init, boost::is_any_of("|"));
	for (size_t i = 0; i < tokens.size(); ++i)
	{
		vector<string> tokens_inner, tokens_rules, intersection_tokens;
		//Format of tokens[i] is: "label:start,end"
		boost::split(tokens_inner, tokens[i], boost::is_any_of(":"));
		if (tokens_inner.size() != 2)
			MTHROW_AND_ERR("Error in medial::sampling::init_time_window_mode - reading token \"%s\" and missing"
				" \":\". format should be label:start,end(,num-num as optional)\n", tokens[i].c_str());
		const string &label = tokens_inner[0];
		boost::split(tokens_rules, tokens_inner[1], boost::is_any_of(","));
		if (tokens_rules.size() != 2 && tokens_rules.size() != 3)
			MTHROW_AND_ERR("Error in medial::sampling::init_time_window_mode - reading token \"%s\" and missing"
				" \",\". format should be start,end. full_token = \"%s\"\n", tokens_inner[1].c_str(), tokens[i].c_str());
		if (label == "all") {
			//mode
			TimeWindowMode temp_mode[2];
			temp_mode[0] = TimeWindow_name_to_type(tokens_rules[0]);
			temp_mode[1] = TimeWindow_name_to_type(tokens_rules[1]);

			set_default(temp_mode);
		}
		else {
			(*this)[med_stof(label)][0] = TimeWindow_name_to_type(tokens_rules[0]);
			(*this)[med_stof(label)][1] = TimeWindow_name_to_type(tokens_rules[1]);
		}
		if (tokens_rules.size() == 3) {
			//aditional args for intersection:
			boost::split(intersection_tokens, tokens_rules[2], boost::is_any_of("-"));
			if (intersection_tokens.size() != 2)
				MTHROW_AND_ERR("Error in medial::sampling::init_time_window_mode - reading token \"%s\" and missing"
					" \",\". format should be number-number. full_token = \"%s\"\n",
					tokens_rules[2].c_str(), tokens[i].c_str());
			if (label != "all") {
				intersection_range_condition[med_stof(label)].first = med_stof(intersection_tokens[0]);
				intersection_range_condition[med_stof(label)].second = med_stof(intersection_tokens[1]);
			}
			else
				set_default_range(med_stof(intersection_tokens[0]), med_stof(intersection_tokens[1]));
		}
	}


}

LabelParams::LabelParams() {
	label_interaction_mode.init_from_string("0:all,before_end|1:before_start,after_start");
	censor_interaction_mode.init_from_string("all:within,all");
	conflict_method = ConflictMode::Drop;
	time_from = 0; 
	time_to = 0;
	censor_time_from = MED_MAT_MISSING_VALUE;
	censor_time_to = MED_MAT_MISSING_VALUE;
	treat_0_class_as_other_classes = false;
}

int LabelParams::init(map<string, string>& map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "time_from")
			time_from = med_stoi(it->second);
		else if (it->first == "time_to")
			time_to = med_stoi(it->second);
		else if (it->first == "conflict_method")
			conflict_method = ConflictMode_name_to_type(it->second);
		else if (it->first == "label_interaction_mode")
			label_interaction_mode.init_from_string(it->second);
		else if (it->first == "censor_interaction_mode")
			censor_interaction_mode.init_from_string(it->second);
		else if (it->first == "censor_time_from")
			censor_time_from = med_stoi(it->second);
		else if (it->first == "censor_time_to")
			censor_time_to = med_stoi(it->second);
		else if (it->first == "treat_0_class_as_other_classes")
			treat_0_class_as_other_classes = med_stoi(it->second);

		else
			MTHROW_AND_ERR("Error in LabelParams::init - unsupported param %s\n", it->first.c_str());
	}
	return 0;
}
