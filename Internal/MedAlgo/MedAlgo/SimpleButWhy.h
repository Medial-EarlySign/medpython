#pragma once
#pragma once

#include <SerializableObject/SerializableObject/SerializableObject.h>

//
// SimpleButWhy :
//
// Contains the SimpleButWhyParams amd SimpleBuyWhyRes objects .
// Some Predictors implement the get_simple_but_why_preds() API
// Which gets the x matrix and the SimpleButWhyParams and returns the predictions and a SimpleButWhyRes Object for each prediction.
//

class SimpleButWhyParams : SerializableObject {

public:
	string method = "simple_default";
	int report_single_features = 1;
	int report_signal_groups = 1;
	int report_positive = 1;
	int report_negative = 1;
	int top_features = 5;
	int top_signal_groups = 5;

	int init_features_names_auto = 1;
	vector<string> feature_names;

	int init_groups_auto = 1;
	map<string, vector<string>> group2names;
	map<string, vector<int>> group2cols;

	// given standard feature names , build signal based groups
	void get_auto_grouping(vector<string> &_feature_names);

	float do_above_score = 0;
	float do_below_score = 0;

	ADD_CLASS_NAME(SimpleButWhyParams)
	ADD_SERIALIZATION_FUNCS(method, report_single_features, report_signal_groups, report_positive, report_negative, top_features, top_signal_groups,
							init_features_names_auto, feature_names, init_groups_auto, group2names, group2cols)
};


class SimpleButWhyRes : SerializableObject {

public:
	vector<pair<string, float>> top_positive_features;
	vector<pair<string, float>> top_negative_features;
	vector<pair<string, float>> top_positive_groups;
	vector<pair<string, float>> top_negative_groups;

	void get_from_contribs(SimpleButWhyParams &bw_params, vector<float> &contribs);

	ADD_CLASS_NAME(SimpleButWhyRes)
	ADD_SERIALIZATION_FUNCS(top_positive_features, top_negative_features, top_positive_groups, top_negative_groups)
};

MEDSERIALIZE_SUPPORT(SimpleButWhyParams)
MEDSERIALIZE_SUPPORT(SimpleButWhyRes)
