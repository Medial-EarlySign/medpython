// Various utilities used in MedProcessTools

#include "MedProcessUtils.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include "Logger/Logger/Logger.h"
#include "MedProcessTools/MedProcessTools/RepProcess.h"
#include "MedProcessTools/MedProcessTools/FeatureGenerator.h"

#define LOCAL_SECTION LOG_MED_UTILS
#define LOCAL_LEVEL	LOG_DEF_LEVEL

char signalName_c[MAX_NAME_LEN + 1];

//..............................................................................
int init_map_from_string(string text, map<string, string>& init_map) {

	if (text == "") return 0;

	// parse text of the format "Name = Value ; Name = Value ; ..."

	// remove white spaces
	text.erase(remove_if(text.begin(), text.end(), ::isspace), text.end());

	if (initialization_text_to_map(text, init_map) == -1)
		return -1;

	//	for (auto rec : init_map)
	//		MLOG("Initializing with \'%s\' = \'%s\'\n", rec.first.c_str(), rec.second.c_str());


	return 0;
}

//..............................................................................
int init_dvec(string& in, vector<int>& out) {

	vector<string> vals;
	split(vals, in, boost::is_any_of(","));

	out.resize(vals.size());
	for (unsigned int i = 0; i < vals.size(); i++)
		out[i] = stoi(vals[i]);

	return 0;
}

//..............................................................................
void get_single_val_from_init_string(string init_s, string attr, string &val_s)
{
	val_s = "";
	if (attr == "") return;

	init_s.erase(remove_if(init_s.begin(), init_s.end(), ::isspace), init_s.end());
	vector<string> fields;
	split(fields, init_s, boost::is_any_of(";="));
	for (int i = 0; i < fields.size(); i++) {
		if (fields[i] == attr) {
			val_s = fields[++i];
			return;
		}
	}

}

//..............................................................................
string int_to_string_digits(int i, int ndigits)
{
	string s;

	int imax = (int)pow(10, ndigits) - 1;
	if (i > imax)
		s = to_string(i);
	else {

		s = to_string(i + imax + 1);
		s.erase(0, 1);

	}

	return s;
}

// Create a required signal names by back propograting : First find what's required by
// the feature generators, and then find add signals required by the rep_porcessors that
// are required ....
//.......................................................................................
void get_all_required_signal_names(unordered_set<string>& signalNames, vector<RepProcessor *>& rep_processors, int position, vector<FeatureGenerator *>& generators) {

	// Collect from generators
	for (auto& generator : generators)
		generator->get_required_signal_names(signalNames);

	// Collect from processors itertively
	for (int i = (int)rep_processors.size() - 1; i > position; i--)
		rep_processors[i]->get_required_signal_names(signalNames, signalNames);
}

void get_all_required_signal_ids(unordered_set<int>& signalIds, vector<RepProcessor *>& rep_processors, int position, vector<FeatureGenerator *>& generators) {


	// Collect from generators
	for (auto& generator : generators)
		generator->get_required_signal_ids(signalIds);

	// Collect from processors itertively
	for (int i = (int)rep_processors.size() - 1; i > position; i--)
		rep_processors[i]->get_required_signal_ids(signalIds, signalIds);

}

void handle_required_signals(vector<RepProcessor *>& processors, vector<FeatureGenerator *>& generators, unordered_set<int>& extra_req_signal_ids,
	vector<int>& all_req_signal_ids_v, vector<unordered_set<int> >& current_required_signal_ids) {

	// Get All required signals
	unordered_set<int> all_req_signal_ids;
	all_req_signal_ids.insert(extra_req_signal_ids.begin(), extra_req_signal_ids.end());
	get_all_required_signal_ids(all_req_signal_ids, processors, -1, generators);

	// Feed to vector
	for (int signalId : all_req_signal_ids)
		all_req_signal_ids_v.push_back(signalId);

	// Get required signals iteratively
	for (size_t i = 0; i < processors.size(); i++) {
		current_required_signal_ids[i].insert(extra_req_signal_ids.begin(), extra_req_signal_ids.end());
		get_all_required_signal_ids(current_required_signal_ids[i], processors, (int)i, generators);
	}
}

// Handle feature names
//.......................................................................................
int find_in_feature_names(vector<string>& names, string& substr) {

	int index = -1;
	for (int i = 0; i < names.size(); i++) {
		if (names[i].find(substr) != string::npos) {
			if (index != -1)
				MTHROW_AND_ERR("%s\n", (string("Got source_feature_name [") + substr + "] which matches both [" + names[i] + "] and [" + names[index] + "]").c_str());
			index = i;
			if (names[i] == substr)
				return index;
			if (names[i].substr(0, 4) == "FTR_")
				if (names[i].substr(names[i].find(".") + 1) == substr)
					return index;
		}
	}

	if (index == -1) {
		string err = string("Got source_feature_name [") + substr + "] which does not match any feature (did you forget to set duplicate=1?). Tried matching to these features:\n";
		for (auto candidate : names)
			err += candidate + "\n";
		if (err.size() >= 300) //out of buffer
			MWARN("%s\n", err.c_str());
		MTHROW_AND_ERR("%s\n", err.c_str());
	}
	else
		return index;
}