#include "FeatureGenExtractTable.h"
#include <fstream>

using namespace std;

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

static unordered_map<string, int> mapper_Rule_Type = {
	{"SET", Rule_Type::SET },
{ "AGE_RANGE", Rule_Type::AGE_RANGE }
};
Rule_Type get_rule_type(const string &op) {
	if (mapper_Rule_Type.find(op) == mapper_Rule_Type.end()) {
		string opts = medial::io::get_list(mapper_Rule_Type, ",");
		MTHROW_AND_ERR("Error can't find option %s, available options are: [%s]\n",
			op.c_str(), opts.c_str());
	}
	return Rule_Type(mapper_Rule_Type.at(op));
}

void KeyRule::parse_rule() {
	parsed_min_range = MED_MAT_MISSING_VALUE;
	parsed_max_range = MED_MAT_MISSING_VALUE;
	if (type == Rule_Type::AGE_RANGE) {
		//parse age rule and test it: - supports "+" "-" with partial range, and just number
		if (rule_value == "All") {} //no limits, take all
		else if (rule_value.find('+') != string::npos)
			parsed_min_range = med_stof(rule_value.substr(0, rule_value.find('+')));
		else if (rule_value.find('-') != string::npos) {
			vector<string> tokens;
			boost::split(tokens, rule_value, boost::is_any_of("-"));
			if (tokens.size() != 2)
				MTHROW_AND_ERR("Error AGE_RANGE type should contain one \"-\" mark and got: \"%s\"\n",
					rule_value.c_str());
			boost::trim(tokens[0]); boost::trim(tokens[1]);
			if (!tokens[0].empty())
				parsed_min_range = med_stof(tokens[0]);
			if (!tokens[1].empty())
				parsed_max_range = med_stof(tokens[1]);
		}
		else { //just a number
			parsed_min_range = med_stof(rule_value);
			parsed_max_range = parsed_min_range;
		}
	}
}

KeyRule::KeyRule(const string &signal, const string &type_str, const string &value) {
	rep_signal = signal;
	type = get_rule_type(type_str);
	rule_value = value;

	parse_rule();
}

void KeyRule::init_lut(MedDictionarySections& dict) {
	if (type == Rule_Type::SET) {
		int section_id = dict.section_id(rep_signal);
		vector<string> set_vals(1);
		set_vals[0] = rule_value;
		dict.prep_sets_lookup_table(section_id, set_vals, lut);
	}
}

void FeatureGenExtractTable::read_rule_table_files() {
	vector<string> key_rows, val_rows;
	vector<string> table_headers, table_data_lines;

	if (rules_config_file.empty())
		MTHROW_AND_ERR("Error must provide rules_config_file in FeatureGenExtractTable\n");
	if (table_file.empty())
		MTHROW_AND_ERR("Error must provide table_file in FeatureGenExtractTable\n");
	//read rules file into key_rows, val_rows
	ifstream rules_fr(rules_config_file);
	string line;
	if (!rules_fr.good())
		MTHROW_AND_ERR("Error FeatureGenExtractTable::read_file - can't read rules file %s\n", rules_config_file.c_str());
	while (getline(rules_fr, line)) {
		boost::trim(line);
		if (line.empty() || line[0] == '#')
			continue;
		if (boost::starts_with(line, "KEY"))
			key_rows.push_back(line);
		else if (boost::starts_with(line, "VALUE"))
			val_rows.push_back(line);
		else {
			rules_fr.close();
			MTHROW_AND_ERR("Error in reading rule in %s. got\n%s\n", rules_config_file.c_str(), line.c_str());
		}
	}
	rules_fr.close();

	//read table file into table_headers, table_data_lines
	ifstream table_fr(table_file);
	if (!table_fr.good())
		MTHROW_AND_ERR("Error FeatureGenExtractTable::read_file - can't read table file %s\n", table_file.c_str());
	while (getline(table_fr, line)) {
		boost::trim(line);
		if (line.empty() || line[0] == '#')
			continue;

		if (table_headers.empty()) //read header
			boost::split(table_headers, line, boost::is_any_of("\t"));
		else //read data
			table_data_lines.push_back(line);
	}
	table_fr.close();

	//parse files using: (key_rows, val_rows, table_headers, table_data_lines)
	unordered_map<string, int> table_header_idx;
	for (int i = 0; i < table_headers.size(); ++i)
		table_header_idx[table_headers[i]] = i;

	//index the mapping of rules into source table column indexes
	unordered_map<string, int> map_rule_signal_to_tbl_key_ind;
	for (const string &rule_line : key_rows)
	{
		vector<string> tokens_rule;
		boost::split(tokens_rule, rule_line, boost::is_any_of("\t"));
		if (tokens_rule.size() != 4)
			MTHROW_AND_ERR("Error FeatureGenExtractTable::read_file - bad file format(%s) in line:\n%s\n",
				rules_config_file.c_str(), rule_line.c_str());
		const string &table_col_name = tokens_rule[1];
		const string &target_sig = tokens_rule[2];

		if (table_header_idx.find(table_col_name) == table_header_idx.end())
			MTHROW_AND_ERR("Error FeatureGenExtractTable::read_file - in rules file(%s) - no match for table col rule %s\nAvailable columns in table:\n%s\n",
				rules_config_file.c_str(), table_col_name.c_str(), medial::io::get_list(table_headers).c_str());

		map_rule_signal_to_tbl_key_ind[target_sig] = table_header_idx.at(table_col_name);
	}

	//update fetched values:
	extracted_names.clear();
	for (const string &val_line : val_rows) {
		vector<string> tokens_val;
		boost::split(tokens_val, val_line, boost::is_any_of("\t"));
		if (tokens_val.size() != 3)
			MTHROW_AND_ERR("Error FeatureGenExtractTable::read_file - bad file format(%s) in rules line:\n%s\n"
				"Should have 3 tokens in VALUE row\n",
				rules_config_file.c_str(), val_line.c_str());
		const string &table_val_col = tokens_val[2];
		extracted_names.push_back(table_val_col);
	}

	//parse each table line to rules:
	key_rules.clear();
	for (const string &table_line : table_data_lines)
	{
		vector<string> tokens_tbl;
		boost::split(tokens_tbl, table_line, boost::is_any_of("\t"));
		if (tokens_tbl.size() > table_headers.size())
			MTHROW_AND_ERR("Error FeatureGenExtractTable::read_file - bad file format (%s) in table line:\n%s\n"
				"Has %zu tokens and in header have %zu tokens\n",
				table_file.c_str(), table_line.c_str(), tokens_tbl.size(), table_headers.size());
		//if smaller - the last are empty values:
		if (tokens_tbl.size() < table_headers.size())
			tokens_tbl.resize(table_headers.size()); //the lasts are empty values

		MapRules full_join_rule;
		//parse table line by rules - construct KeyRule for each row:
		for (const string &rule_line : key_rows)
		{
			vector<string> tokens_rule;
			boost::split(tokens_rule, rule_line, boost::is_any_of("\t"));
			if (tokens_rule.size() != 4)
				MTHROW_AND_ERR("Error FeatureGenExtractTable::read_file - bad file format(%s) in line:\n%s\n",
					rules_config_file.c_str(), rule_line.c_str());
			const string &target_sig = tokens_rule[2];
			const string &rule_type = tokens_rule[3];

			int col_num = map_rule_signal_to_tbl_key_ind.at(target_sig);
			const string &table_val = tokens_tbl[col_num];

			KeyRule kr(target_sig, rule_type, table_val);
			full_join_rule.rules.push_back(kr);
		}

		//add values:
		for (const string &val_line : val_rows) {
			vector<string> tokens_val;
			boost::split(tokens_val, val_line, boost::is_any_of("\t"));
			const string &table_val_col = tokens_val[1];
			if (table_header_idx.find(table_val_col) == table_header_idx.end())
				MTHROW_AND_ERR("Error FeatureGenExtractTable::read_file - can't find rule value column %s in table\n"
					"Available table headers: [%s]\n",
					table_val_col.c_str(), medial::io::get_list(table_headers).c_str());
			int col_idx = table_header_idx.at(table_val_col);
			const string &token_before_parse = tokens_tbl[col_idx];
			float table_val = missing_val;
			if (!token_before_parse.empty())
				table_val = med_stof(token_before_parse);
			full_join_rule.values.push_back(table_val);
		}

		//add rule:
		key_rules.push_back(full_join_rule);
	}

	if (reverse_rule_order)
		reverse(key_rules.begin(), key_rules.end());
	MLOG("Read %zu rules from %s\n", key_rules.size(), table_file.c_str());
}

int FeatureGenExtractTable::init(map<string, string>& mapper) {
	for (const auto &it : mapper)
	{
		//! [FeatureGenExtractTable::init]
		if (it.first == "rules_config_file")
			rules_config_file = it.second;
		else if (it.first == "table_file")
			table_file = it.second;
		else if (it.first == "table_nice_name")
			table_nice_name = it.second;
		else if (it.first == "reverse_rule_order")
			reverse_rule_order = med_stoi(it.second) > 0;
		else if (it.first == "missing_val")
			missing_val = med_stof(it.second);
		else if (it.first == "fg_type" || it.first == "tags") {}
		else
			MTHROW_AND_ERR("Error FeatureGenExtractTable::init - unsupported argument %s\n",
				it.first.c_str());
		//! [FeatureGenExtractTable::init]
	}
	read_rule_table_files();
	if (table_nice_name.empty())
		table_nice_name = table_file;

	req_signals.clear();
	req_signal_ids.clear();
	if (!key_rules.empty()) {
		const vector<KeyRule> &all_rules = key_rules.front().rules;
		//fetch AND condition signal rules
		for (const KeyRule &kr : all_rules)
		{
			if (kr.type == Rule_Type::SET)
				req_signals.push_back(kr.rep_signal);

			else if (kr.type == Rule_Type::AGE_RANGE)
				//requires BDATE - if more complicated (other feature generator - use this feature generator to retrieve req_features)
				req_signals.push_back("BDATE");
			else
				MTHROW_AND_ERR("Error FeatureGenExtractTable::set_signal_ids - not impelmented for rule type %d\n",
					int(kr.type));
		}
	}

	return 0;
}

void FeatureGenExtractTable::set_signal_ids(MedSignals& sigs) {
	//set req_signal from rules:
	req_signals.clear();
	req_signal_ids.clear();
	if (!key_rules.empty()) {
		const vector<KeyRule> &all_rules = key_rules.front().rules;
		//fetch AND condition signal rules
		for (const KeyRule &kr : all_rules)
		{
			if (kr.type == Rule_Type::SET) {
				req_signals.push_back(kr.rep_signal);
				int sid = sigs.sid(kr.rep_signal);
				if (sid < 0)
					MTHROW_AND_ERR("Error FeatureGenExtractTable::set_signal_ids - can't find signal %s\n",
						kr.rep_signal.c_str());
				req_signal_ids.push_back(sid);
			}
			else if (kr.type == Rule_Type::AGE_RANGE) {
				//requires BDATE - if more complicated (other feature generator - use this feature generator to retrieve req_features)
				req_signals.push_back("BDATE");
				int sid = sigs.sid("BDATE");
				if (sid < 0)
					MTHROW_AND_ERR("Error FeatureGenExtractTable::set_signal_ids - can't find BDATE\n");
				req_signal_ids.push_back(sid);
			}
			else
				MTHROW_AND_ERR("Error FeatureGenExtractTable::set_signal_ids - not impelmented for rule type %d\n",
					int(kr.type));
		}
	}

}

void FeatureGenExtractTable::init_tables(MedDictionarySections& dict) {
	for (size_t i = 0; i < key_rules.size(); ++i)
		for (size_t j = 0; j < key_rules[i].rules.size(); ++j)
			key_rules[i].rules[j].init_lut(dict);
	missing_values_cnt = 0;
}

void FeatureGenExtractTable::get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const {
	if (!key_rules.empty()) {
		const vector<KeyRule> &all_rules = key_rules.front().rules;
		for (size_t i = 0; i < all_rules.size(); ++i)
			if (all_rules[i].type == Rule_Type::SET)
				signal_categories_in_use[req_signals[i]].push_back(all_rules[i].rule_value);
	}
}

void FeatureGenExtractTable::set_names() {
	//convert extracted_names => names. maybe wants to use FTR_ or something
	names = extracted_names;
}

bool KeyRule::test_rule(float val) const {

	switch (type)
	{
	case Rule_Type::SET:
		return lut[(int)val] > 0;
	case Rule_Type::AGE_RANGE:
		return (parsed_min_range == MED_MAT_MISSING_VALUE || val >= parsed_min_range)
			&& (parsed_max_range == MED_MAT_MISSING_VALUE || val <= parsed_max_range);
	default:
		MTHROW_AND_ERR("Error KeyRule::test_rule - Unsupported type %d\n", (int)type);;
	}
}

bool MapRules::join(const vector<float> &join_vals) const {
	//check each condition
	bool can_join = true;
	for (size_t i = 0; i < join_vals.size() && can_join; ++i)
		can_join = rules[i].test_rule(join_vals[i]);
	return can_join;
}

int FeatureGenExtractTable::_generate(PidDynamicRec& in_rep, MedFeatures& features,
	int index, int num, vector<float *> &_p_data) {
	if (key_rules.empty())
		MTHROW_AND_ERR("Error FeatureGenExtractTable::_generate no rules\n");

	const vector<KeyRule> &rules_cond = key_rules.front().rules;
	int join_n = (int)rules_cond.size();
	MedSample *p_samples = &(features.samples[index]);
	vector<UniversalSigVec> usvs(req_signals.size());
	vector<float> join_vals(join_n);

	for (int i = 0; i < num; ++i) {
		//get input signals:
		for (size_t k = 0; k < usvs.size(); ++k)
			in_rep.uget(req_signal_ids[k], i, usvs[k]);

		//update values in join_vals to calculate join values in sample:
		int ind_sig = 0;
		for (size_t k = 0; k < join_n; ++k)
		{
			int byear, bdate;
			float set_val;
			switch (rules_cond[k].type)
			{
			case Rule_Type::SET:
				set_val = usvs[ind_sig].Val(0); //currently only channel of value. and no time window - features lie gender, county, race that has no time channels
				join_vals[k] = set_val;
				++ind_sig;
				break;
			case Rule_Type::AGE_RANGE:
				//calculate Age: has byear
				if (usvs[ind_sig].n_val_channels() > 0)
					bdate = usvs[ind_sig].Val(0);
				else
					bdate = usvs[ind_sig].Time(0);
				byear = int(bdate / 10000);
				set_val = float(med_time_converter.convert_times(features.time_unit, MedTime::Date, p_samples[i].time) / 10000) - byear;
				join_vals[k] = set_val;
				++ind_sig;
				break;
			default:
				MTHROW_AND_ERR("Error FeatureGenExtractTable::_generate - Unsupported type %d\n", (int)rules_cond[k].type);
			}
		}

		//find rule to match with join_vals:
		int rule_idx = -1;
		for (int k = 0; k < key_rules.size() && rule_idx < 0; ++k)
			if (key_rules[k].join(join_vals))
				rule_idx = k;

		if (rule_idx < 0)
#pragma omp atomic
			++missing_values_cnt;
		//Fill all feature values:
		for (size_t k = 0; k < names.size(); ++k)
		{
			float *p_feat = _p_data[k] + index;
			if (rule_idx < 0)
				p_feat[i] = missing_val;
			else
				p_feat[i] = key_rules[rule_idx].values[k];
		}
	}

	return 0;
}

int FeatureGenExtractTable::filter_features(unordered_set<string>& validFeatures) {
	//go pass throuhh names and keep only those who are inside:
	vector<string> final_list, final_list_ext;
	vector<int> keep_idx;
	for (int i = 0; i < names.size(); ++i)
		if (validFeatures.find(names[i]) != validFeatures.end()) {
			final_list.push_back(names[i]);
			final_list_ext.push_back(extracted_names[i]);
			keep_idx.push_back(i);
		}
	names = move(final_list);
	extracted_names = move(final_list_ext);
	//extracted_names is not important to update - but it's not time consuming and for the whole validty of the object

	//update rules to keep only those features:
	for (size_t i = 0; i < key_rules.size(); ++i)
	{
		vector<float> keep_vals(keep_idx.size());
		for (size_t j = 0; j < keep_idx.size(); ++j)
			keep_vals[j] = key_rules[i].values[keep_idx[j]];
		key_rules[i].values = move(keep_vals);
	}

	return (int)names.size();
}

void FeatureGenExtractTable::prepare(MedFeatures &features, MedPidRepository& rep, MedSamples& samples) {
	FeatureGenerator::prepare(features, rep, samples);
	missing_values_cnt = 0;
}

void FeatureGenExtractTable::make_summary() {
	if (missing_values_cnt > 0)
		MLOG("FeatureGenExtractTable :: has %d missing samples to join with table %s\n",
			missing_values_cnt, table_nice_name.c_str());
}