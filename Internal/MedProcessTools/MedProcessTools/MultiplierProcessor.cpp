#include <boost/regex.hpp>
#include "FeatureProcess.h"

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int MultiplierProcessor::init(map<string, string>& mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		string val = boost::trim_copy(entry.second);
		//! [MultiplierProcessor::init]
		if (field == "selected_tags") { if (!val.empty()) boost::split(selected_tags, val, boost::is_any_of(",")); }
		else if (field == "multiplier_name") multiplier_name = val;
		else if (field == "divide") divide = med_stoi(val) > 0;
		else if (field == "verbose") verbose = med_stoi(val) > 0;
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknown parameter \'%s\' for MultiplierProcessor\n", field.c_str());
		//! [MultiplierProcessor::init]
	}
	if (multiplier_name.empty())
		MTHROW_AND_ERR("Error MultiplierProcessor::init - must provide multiplier_name");
	return 0;
}

void MultiplierProcessor::dprint(const string &pref, int fp_flag) {
	if (fp_flag > 0) {
		string tags_str = medial::io::get_list(selected_tags);
		MLOG("%s :: MultiplierProcessor :: selected_tags(%zu)=[%s], feature_multipler=%s, divide=%d\n",
			pref.c_str(), selected_tags.size(), tags_str.c_str(), multiplier_name.c_str(), int(divide));
	}
}

int MultiplierProcessor::_apply(MedFeatures& features, unordered_set<int>& ids) {
	unordered_set<string> s(selected_tags.begin(), selected_tags.end());
	vector<string> all_names;
	features.get_feature_names(all_names);
	string resolved_multi = all_names[find_in_feature_names(all_names, multiplier_name)];
	const vector<float> &multiplier_vec = features.data.at(resolved_multi);

	vector<string> touch_cnt;
	for (auto it = features.tags.begin(); it != features.tags.end(); ++it) {
		if (it->first == resolved_multi)
			continue;
		const unordered_set<string> &feature_tags = it->second;
		bool found_match = s.empty(); //if empty - select all
		for (const string &candidate_tag : feature_tags)
		{
			for (const string& substring : s) {
				boost::regex regi(substring);
				if (boost::regex_match(candidate_tag, regi)) {
					found_match = true;
					break;
				}
			}
			if (found_match)
				break;
		}

		if (found_match) {
			touch_cnt.push_back(it->first);
			//do the multiplication or division:
			vector<float> &vec = features.data.at(it->first);
			for (size_t i = 0; i < vec.size(); ++i)
			{
				if (vec[i] != MED_MAT_MISSING_VALUE) {
					if (!divide)
						vec[i] *= multiplier_vec[i];
					else {
						if (multiplier_vec[i] != 0)
							vec[i] /= multiplier_vec[i];
						else
							vec[i] = MED_MAT_MISSING_VALUE;
					}
				}
			}
		}
	}

	if (verbose && !touch_cnt.empty()) {
		//MLOG("INFO: multiplier using %s touched %zu features: [%s]\n", resolved_multi.c_str(),
		//	touch_cnt.size(), medial::io::get_list(touch_cnt).c_str());
	}
	return 0;
}