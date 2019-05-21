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
#include <boost/optional/optional.hpp>

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#define CHECK_CRC 0

using namespace boost::property_tree;

void parse_my_json_to_pt(istringstream &no_comments_stream, ptree &pt);

void parse_my_json_to_pt(string &str, ptree &pt) {
	istringstream i_str(str);
	parse_my_json_to_pt(i_str, pt);
}

void parse_my_json_to_pt(istringstream &no_comments_stream, ptree &pt)
{
	string my_line;
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
		MTHROW_AND_ERR("json parsing error [%s] at line %lu\n", e.message().c_str(), e.line());
	}
}


void get_prefix_suffix_from_tokens(const string& single_attr_value, string& small_file, string& ref_node, string& prefix, string& suffix) {
	vector<string> tokens;
	boost::split(tokens, single_attr_value, boost::is_any_of(";"));
	for (string token : tokens) {
		if (boost::starts_with(token, "prefix:"))
			prefix = token.substr(7);
		else if (boost::starts_with(token, "suffix:"))
			suffix = token.substr(7);
		else if (boost::starts_with(token, "file:"))
			small_file = tokens[0].substr(5);
		else if (boost::starts_with(token, "file_rel:"))
			small_file = tokens[0].substr(9);
		else if (boost::starts_with(token, "path_rel:"))
			small_file = tokens[0].substr(9);
		else if (boost::starts_with(token, "ref:"))
			ref_node = tokens[0].substr(4);
		else MTHROW_AND_ERR("dont know how to handle token [%s]\n", token.c_str());
	}
}
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
				string small_file = "", ref_node = "", prefix = "", suffix = "";
				if (boost::starts_with(single_attr_value, "file:")) {
					//e.g. "signal": "file:my_list.txt;prefix:ppp;suffix:sss" - file can be relative					
					get_prefix_suffix_from_tokens(single_attr_value, small_file, ref_node, prefix, suffix);
					vector<string> my_list;
					fill_list_from_file(make_absolute_path(fname, small_file), my_list);
					for (string s : my_list)
						current_attr_values.push_back(prefix + parse_key_val(attr_name, s) + suffix);
				}
				else if (boost::starts_with(single_attr_value, "file_rel:")) { //wih relative paths
					//e.g. "signal": "file:my_list.txt;prefix:ppp;suffix:sss" - file can be relative					
					get_prefix_suffix_from_tokens(single_attr_value, small_file, ref_node, prefix, suffix);
					vector<string> my_list;
					fill_list_from_file(make_absolute_path(fname, small_file, true), my_list);
					for (string s : my_list)
						current_attr_values.push_back(prefix + parse_key_val(attr_name, s) + suffix);
				}
				else if (boost::starts_with(single_attr_value, "path_rel:")) { //wih relative paths
																			   //e.g. "signal": "file:my_list.txt;prefix:ppp;suffix:sss" - file can be relative					
					get_prefix_suffix_from_tokens(single_attr_value, small_file, ref_node, prefix, suffix);

					string abs_path = make_absolute_path(fname, small_file, true);
					current_attr_values.push_back(prefix + parse_key_val(attr_name, abs_path) + suffix);
				}
				else if (boost::starts_with(single_attr_value, "ref:")) {
					get_prefix_suffix_from_tokens(single_attr_value, small_file, ref_node, prefix, suffix);
					auto my_ref = root.get_child(ref_node);
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
			if (current_attr_values.empty())
				MTHROW_AND_ERR("[%s] has an empty val [%s]\n", attr_name.c_str(), single_attr_value.c_str());
			all_action_attrs.push_back(current_attr_values);
		}
	}
}

void MedModel::init_from_json_file_with_alterations(const string &fname, vector<string>& alterations) {
	run_current_path = boost::filesystem::path(fname).parent_path().string();
	string json_contents = json_file_to_string(0, fname, alterations, "", true);

	if (init_from_json_string(json_contents, fname) == 1)
		init_from_json_file_with_alterations_version_1(fname, alterations);
}

int MedModel::init_from_json_string(string& json_contents, const string& fname) {
	istringstream no_comments_stream(json_contents);
	MLOG("MedModel:: init model from json file [%s]:\n", fname.c_str());

	ptree pt;
	parse_my_json_to_pt(no_comments_stream, pt);
	this->model_json_version = pt.get<int>("model_json_version", model_json_version);
	MLOG_D("\nmodel_json_version [%d]\n", model_json_version);
	if (model_json_version <= 1)
		return 1;

	boost::optional<int> v = pt.get_optional<int>("generate_masks_for_features");
	if (v)
		this->generate_masks_for_features = v.get();

	//MLOG("debug=====> :: generate_masks_for_features %d\n", generate_masks_for_features);

	string ser = pt.get<string>("serialize_learning_set", to_string(this->serialize_learning_set).c_str());
	this->serialize_learning_set = stoi(ser);
	int rp_set = 0, fp_set = 0, pp_set = 0;
	for (auto &p : pt.get_child("model_actions")) {
		vector<vector<string>> all_action_attrs;
		auto& action = p.second;

		string action_type = action.get<string>("action_type").c_str();
		if (boost::starts_with(action_type, "change_path:")) {
			//change json base_path fo relative paths to work:
			string new_path = boost::replace_all_copy(action_type, "change_path:", "");
			run_current_path = new_path;
			MLOG_D("Changed base path to %s\n", new_path.c_str());
		}
		else if (action_type == "rp_set" || action_type == "fp_set") {
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
					MTHROW_AND_ERR("[%s] set %d expanded to 0 combinations! did you put an empty list inside a []?!\n", action_type.c_str(), process_set);
				if (duplicate == 1 && all_combinations.size() != 1)
					MTHROW_AND_ERR("duplicate is currently supported only for sets with a single action. [%s] set %d has one member which expanded to %d actions\n",
						action_type.c_str(), process_set, (int)all_combinations.size());
				for (string c : all_combinations)
					add_process_to_set(process_set, duplicate, c);
				num_actions += (int)all_combinations.size();
				if (first_action_added == "")
					first_action_added = all_combinations[0];
			}
			MLOG_D("added %d actions to [%s] set %d, first of which was [%s]\n", num_actions, action_type.c_str(), process_set, first_action_added.c_str());
		}
		else if (action_type == "rep_processor" || action_type == "feat_generator" || action_type == "feat_processor" || action_type == "post_processor") {
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
			else if (action_type == "post_processor") {
				process_set = pp_set++;
				set_name = "pp_set";
			}
			else MTHROW_AND_ERR("unknown action_type [%s]\n", action_type.c_str());
			int duplicate = 0;
			parse_action(action, all_action_attrs, duplicate, pt, fname);
			if (duplicate == 1)
				MTHROW_AND_ERR("duplicate action requested and not inside a set!");
			vector<string> all_combinations;
			concatAllCombinations(all_action_attrs, 0, "", all_combinations);
			if (all_combinations.empty())
				MTHROW_AND_ERR("set %d expanded to 0 combinations! did you put an empty list inside a []?!\n", process_set);
			if (all_combinations.size() > 1 && (action_type == "rep_processor" || action_type == "feat_processor"))
				MTHROW_AND_ERR("action_type [%s] expanded to %d combinations, which is possible only inside a set! first instance is [%s]\n",
					action_type.c_str(), (int)all_combinations.size(), all_combinations[0].c_str());
			for (string c : all_combinations)
				add_process_to_set(process_set, duplicate, c);
			MLOG_D("added %d actions to [%s] set %d, first of which was [%s]\n", all_combinations.size(), set_name.c_str(), process_set, all_combinations[0].c_str());
		}
		else MTHROW_AND_ERR("unknown action_type [%s]\n", action_type.c_str());
	}
	if (pt.count("predictor") > 0) {
		auto my_pred = pt.get_child("predictor");
		auto my_pred_params = pt.get_child("predictor_params");
		set_predictor(my_pred.data(), my_pred_params.data());
	}
	else MWARN("NOTE: no [predictor] node found in file\n");

	return 0;
}

//-----------------------------------------------------------------------------------------------------
// next option gets a separate json with just pre_processors inside it and adds them to rep_processors
// at the moment only direct serial processing is allowed, as there's no learn in pre_processors
// and we anyway parallelize on the pids level
//
// format for pre_processors is very much like the format for rep_processors
//
// format in general is :
// { "pre_processors" : [ {"rp_type": "history_limit", ...} , ... ] }
//
// use "" and a file name to start from a file, or a string and empty or given fname to start from a string
//-----------------------------------------------------------------------------------------------------
void MedModel::add_pre_processors_json_string_to_model(string in_json, string fname)
{
	string json_contents = in_json;
	if (json_contents == "") {
		vector<string> dummy_alts;
		json_contents = json_file_to_string(0, fname, dummy_alts);
	}
	ptree pt;
	parse_my_json_to_pt(json_contents, pt);

	size_t n = 0;
	for (auto &p : pt.get_child("pre_processors")) {
		vector<vector<string>> all_action_attrs;
		auto& action = p.second;
		//string action_type = action.get<string>("action_type").c_str();
		int duplicate = 0;
		parse_action(action, all_action_attrs, duplicate, pt, fname);
		if (duplicate == 1)
			MTHROW_AND_ERR("duplicate action requested and not inside a set!");
		vector<string> all_combinations;
		concatAllCombinations(all_action_attrs, 0, "", all_combinations);
		if (all_combinations.empty())
			MTHROW_AND_ERR("pre processor expanded to 0 combinations! did you put an empty list inside a []?!\n");

		for (int idx = 0; idx < all_combinations.size(); idx++) {
			string c = all_combinations[idx];
			MLOG("Adding pre_processor: %s\n", c.c_str());
			insert_rep_processor(c, idx);
		}
		MLOG("added %d pre processors, first of which was [%s]\n", all_combinations.size(), all_combinations[0].c_str());
		n += all_combinations.size();
	}
	MLOG("Succesfully added %d pre_processors\n", n);
}

//-----------------------------------------------------------------------------------------------------
// Same as above, for adding post-processors for the existing ones. Returing the number of PP's added
//-----------------------------------------------------------------------------------------------------
int MedModel::add_post_processors_json_string_to_model(string in_json, string fname)
{
	string json_contents = in_json;
	if (json_contents == "") {
		vector<string> dummy_alts;
		json_contents = json_file_to_string(0, fname, dummy_alts);
	}
	ptree pt;
	parse_my_json_to_pt(json_contents, pt);

	size_t n = 0;
	for (auto &p : pt.get_child("post_processors")) {
		vector<vector<string>> all_action_attrs;
		auto& action = p.second;
		//string action_type = action.get<string>("action_type").c_str();
		int duplicate = 0;
		parse_action(action, all_action_attrs, duplicate, pt, fname);
		if (duplicate == 1)
			MTHROW_AND_ERR("duplicate action requested and not inside a set!");
		vector<string> all_combinations;
		concatAllCombinations(all_action_attrs, 0, "", all_combinations);
		if (all_combinations.empty())
			MTHROW_AND_ERR("post processor expanded to 0 combinations! did you put an empty list inside a []?!\n");

		for (int idx = 0; idx < all_combinations.size(); idx++) {
			string c = all_combinations[idx];
			MLOG("Adding post_processor: %s\n", c.c_str());
			PostProcessor *post_proc = PostProcessor::create_processor(c);
			post_processors.push_back(post_proc);
		}
		MLOG("added %d post processors, first of which was [%s]\n", all_combinations.size(), all_combinations[0].c_str());
		n += all_combinations.size();
	}
	MLOG("Succesfully added %d post_processors\n", n);

	return n;
}
