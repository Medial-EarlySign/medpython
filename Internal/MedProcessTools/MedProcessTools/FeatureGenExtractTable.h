#ifndef _FTR_FEATURE_GEN_EXTRACT_TABLE_H_
#define _FTR_FEATURE_GEN_EXTRACT_TABLE_H_
#include "FeatureGenerator.h"

typedef enum
{
	SET,
	AGE_RANGE
} Rule_Type;

/**
* A class that represent a simple rule to calculate when "joining" with patient sample
*/
class KeyRule : public SerializableObject {
private:
	float parsed_min_range;
	float parsed_max_range;
	vector<char> lut; //lut for the specific rule

	void parse_rule();
public:
	string rep_signal; ///< the repository target signal
	//maybe in the future support more complex rule to accept feature_generator init to generate feature and than
	// use the _generate of this feature on each sample to calculate the range - not for now!
	Rule_Type type; ///< the rule type
	string rule_value; ///< the rule value to match

	KeyRule() {};

	KeyRule(const string &signal, const string &type_str, const string &value);

	/// function to store lut
	void init_lut(MedDictionarySections& dict);

	/// tests the value condition with the rule
	bool test_rule(float val) const;

	ADD_CLASS_NAME(KeyRule)
		ADD_SERIALIZATION_FUNCS(rep_signal, type, rule_value, parsed_min_range, parsed_max_range)
};

/**
* A class that represents set of rules with AND condition to join with each sample.
* It has vector of match values.
* If the patient has multiple choises - the first one is selected (it's order from more more specific to less specific)
*/
class MapRules : public SerializableObject {
public:
	vector<KeyRule> rules; ///< the rules with AND condition. the rules are prioritories by order the more specific to less specific
	vector<float> values; ///< the matched values for each "join". the names are in different file

	/// test for full join condition with the rules
	bool join(const vector<float> &join_vals) const;

	ADD_CLASS_NAME(MapRules)
		ADD_SERIALIZATION_FUNCS(rules, values)
};

/**
* A Feature generatoor that join population properties like Age, Gender, Race with extrenal table with
* numbers for those population (for example death rate, cancer and more...)
*
*/
class FeatureGenExtractTable : public FeatureGenerator {
private:
	vector<MapRules> key_rules; ///< each row is join combination of rules with and
	vector<string> extracted_names; ///< the extracted names of the values

	/// <summary>
	/// @param table_files - a path to the external table with the numbers and data
	/// @param rules_config_file - a path to rules file (TAB delimeted). Each row either starts with "KEY" or "VALUE"
	/// KEY row - define rule in AND condition to join with: table_column_name, repository_signal_name, rule_type(SET or AGE_RANGE)
	/// Value row - define value to fetch: name_of_column_in_table_file, name_of_feature
	/// @param reverse_table_order - if true will reverse table rules
	/// </summary>
	void read_rule_table_files();
	int missing_values_cnt;
public:
	string rules_config_file; ///< path to rules config file
	string table_file; ///< path to table file with numbers
	bool reverse_rule_order; ///< if true will reverse the rule order for priority
	string table_nice_name; ///< nice name for printing
	//object uses also req_signals, req_signals_ids - calculates from key_rules when set_signal_ids is called

	FeatureGenExtractTable() {
		generator_type = FTR_GEN_EXTRACT_TBL;
		missing_val = MED_MAT_MISSING_VALUE;
		missing_values_cnt = 0;
		table_nice_name = "";
		table_file = "";
		rules_config_file = "";
		reverse_rule_order = false;
	}

	void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<FeatureGenExtractTable *>(generator)); }

	/// The parsed fields from init command.
	/// @snippet FeatureGenExtractTable.cpp FeatureGenExtractTable::init
	int init(map<string, string>& mapper);

	void set_signal_ids(MedSignals& sigs);
	void init_tables(MedDictionarySections& dict);
	void set_names();

	void prepare(MedFeatures &features, MedPidRepository& rep, MedSamples& samples);

	int _generate(PidDynamicRec& in_rep, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	int filter_features(unordered_set<string>& validFeatures);

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

	void make_summary();

	ADD_CLASS_NAME(FeatureGenExtractTable)
		ADD_SERIALIZATION_FUNCS(generator_type, req_signals, names, missing_val, tags, iGenerateWeights, key_rules, extracted_names, table_nice_name)
};

MEDSERIALIZE_SUPPORT(KeyRule)
MEDSERIALIZE_SUPPORT(MapRules)
MEDSERIALIZE_SUPPORT(FeatureGenExtractTable)

#endif
