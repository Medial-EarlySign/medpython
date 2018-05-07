#include "MedModel.h"
#include "MedProcessUtils.h"
#include <omp.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/find.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <string>
#include "StripComments.h"

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#define CHECK_CRC 0

using namespace boost::property_tree;

void MedModel::parse_action(basic_ptree<string, string>& action, vector<vector<string>>& all_action_attrs, int& duplicate, ptree& root, const string& fname) {
	all_action_attrs.clear();
	for (ptree::value_type &attr : action) {
		string attr_name = attr.first;
		string single_attr_value = attr.second.data();
		if (attr_name == "action_type")
			continue;
		if (attr_name == "duplicate") {
			boost::algorithm::to_lower(single_attr_value);
			if (single_attr_value == "no" || single_attr_value == "n" || single_attr_value == "0")
				duplicate = 0;
			else if (single_attr_value == "yes" || single_attr_value == "y" || single_attr_value == "1")
				duplicate = 1;
			else MTHROW_AND_ERR("unknown value for duplicate [%s]\n", single_attr_value.c_str());				
		}
		else {
			vector<string> current_attr_values;
			if (single_attr_value.length() > 0) {				
				if (boost::starts_with(single_attr_value, "file:")) {
					//e.g. "signal": "file:my_list.txt" - file can be relative
					vector<string> my_list;
					string small_file = single_attr_value.substr(5);
					fill_list_from_file(make_absolute_path(fname, small_file), my_list);
					for (string s : my_list)
						current_attr_values.push_back(parse_key_val(attr_name, s));
				}
				else if (boost::starts_with(single_attr_value, "ref:")) {
					auto my_ref = root.get_child(single_attr_value.substr(4));
					for (auto &r : my_ref)
						//e.g. "signal": "ref:signals"
						current_attr_values.push_back(parse_key_val(attr_name, r.second.data()));
				}
				else
					// e.g. "fg_type": "gender"
					current_attr_values.push_back(parse_key_val(attr_name, single_attr_value));
			}
			else
				//e.g. "type": ["last", "slope"]
				for (ptree::value_type &attr_value : attr.second)
					current_attr_values.push_back(parse_key_val(attr_name, attr_value.second.data()));
			all_action_attrs.push_back(current_attr_values);
		}
	}
}

void MedModel::init_from_json_file_with_alterations(const string &fname, vector<string>& alterations) {
	string json_contents = file_to_string(0, fname, alterations);
	istringstream no_comments_stream(json_contents);

	int i = 5;
	MLOG("init model from json file [%s], stripping comments and displaying first %d lines:\n", fname.c_str(), i);
	string my_line;
	while (i-- > 0 && getline(no_comments_stream, my_line))
		MLOG("%s\n", my_line.c_str());
	no_comments_stream.clear();
	no_comments_stream.seekg(0);

	ptree pt;
	try {
		read_json(no_comments_stream, pt);
	}
	catch (json_parser_error e) {
		no_comments_stream.clear();
		no_comments_stream.seekg(0);
		MLOG("json parsing error [%s] at line %d\n", e.message().c_str(), e.line());
		for (int i = 1; getline(no_comments_stream, my_line); i++) {
			if (abs((int)e.line() - i) < 3)
				MLOG("%d\t%s\n", i, my_line.c_str());
		}
		MTHROW_AND_ERR("json parsing error [%s] at line %d\n", e.message().c_str(), e.line());
	}
	this->model_json_version = pt.get<int>("model_json_version", model_json_version);
	MLOG("\nmodel_json_version [%d]\n", model_json_version);
	if (model_json_version <= 1) {
		init_from_json_file_with_alterations_version_1(fname, alterations);
		return;
	}
	string ser = pt.get<string>("serialize_learning_set", to_string(this->serialize_learning_set).c_str());	
	this->serialize_learning_set = stoi(ser);
	int rp_set = 0, fp_set = 0;
	for (auto &p : pt.get_child("model_actions")) {
		vector<vector<string>> all_action_attrs;	
		auto& action = p.second;
		string action_type = action.get<string>("action_type").c_str();
		if (action_type == "rp_set" || action_type == "fp_set") {
			int process_set;
			if (action_type == "rp_set") process_set = rp_set++;
			else process_set = fp_set++;			
			int num_members = (int)action.get_child("members").size();
			int num_actions = 0;
			string first_action_added = "";
			for (auto &member : action.get_child("members")) {
				int duplicate = 0;
				parse_action(member.second, all_action_attrs, duplicate, pt, fname);
				if (duplicate == 1 && num_members != 1)
					MTHROW_AND_ERR("duplicate is currently supported only for sets with a single action. [%s] set %d has %d members, please separate it to multiple sets\n",
						action_type.c_str(), process_set, num_members);
				vector<string> all_combinations;
				concatAllCombinations(all_action_attrs, 0, "", all_combinations);
				if (all_combinations.empty())
					MTHROW_AND_ERR("set %d expanded to 0 combinations! did you put an empty list inside a []?!\n");
				if (duplicate == 1 && all_combinations.size() != 1)
					MTHROW_AND_ERR("duplicate is currently supported only for sets with a single action. [%s] set %d has one member which expanded to %d actions\n",
						action_type.c_str(), process_set, (int)all_combinations.size());
				for (string c : all_combinations)
					add_process_to_set(process_set, duplicate, c);
				num_actions += (int)all_combinations.size();	
				if (first_action_added == "")
					first_action_added = all_combinations[0];
			}			
			MLOG("added %d actions to [%s] set %d, first of which was [%s]\n", num_actions, action_type.c_str(), process_set, first_action_added.c_str());
		}
		else if (action_type == "rep_processor" || action_type == "feat_generator" || action_type == "feat_processor") {
			int process_set; string set_name = "";
			if (action_type == "rep_processor") {
				process_set = rp_set++;
				set_name = "rp_set";
			}
			else if (action_type == "feat_processor") {
				process_set = fp_set++;
				set_name = "fp_set";
			}
			else if (action_type == "feat_generator") {
				process_set = 0;
				set_name = "fg_set";
			}
			else MTHROW_AND_ERR("unknown action_type [%s]\n", action_type.c_str());
			int duplicate = 0;
			parse_action(action, all_action_attrs, duplicate, pt, fname);
			if (duplicate == 1)
				MTHROW_AND_ERR("duplicate action requested and not inside a set!");
			vector<string> all_combinations;
			concatAllCombinations(all_action_attrs, 0, "", all_combinations);
			if (all_combinations.empty())
				MTHROW_AND_ERR("set %d expanded to 0 combinations! did you put an empty list inside a []?!\n");
			if (all_combinations.size() > 1 && (action_type == "rep_processor" || action_type == "feat_processor"))
				MTHROW_AND_ERR("action_type [%s] expanded to %d combinations, which is possible only inside a set! first instance is [%s]\n",
					action_type.c_str(), (int)all_combinations.size(), all_combinations[0].c_str());			
			for (string c : all_combinations)
				add_process_to_set(process_set, duplicate, c);
			MLOG("added %d actions to [%s] set %d, first of which was [%s]\n", all_combinations.size(), set_name.c_str(), process_set, all_combinations[0].c_str());
		} else MTHROW_AND_ERR("unknown action_type [%s]\n", action_type.c_str());
	}
	if (pt.count("predictor") > 0) {
		auto my_pred = pt.get_child("predictor");
		auto my_pred_params = pt.get_child("predictor_params");
		set_predictor(my_pred.data(), my_pred_params.data());
	}
	else MWARN("NOTE: no [predictor] node found in file\n");
}
