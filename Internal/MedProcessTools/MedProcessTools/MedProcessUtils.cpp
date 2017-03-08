// Various utilities used in MedProcessTools

#include "MedProcessUtils.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include "Logger/Logger/Logger.h"

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
	if (attr=="") return;

	init_s.erase(remove_if(init_s.begin(), init_s.end(), ::isspace), init_s.end());
	vector<string> fields;
	split(fields, init_s, boost::is_any_of(";="));
	for (int i=0; i<fields.size(); i++) {
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
