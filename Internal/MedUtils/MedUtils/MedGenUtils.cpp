#include "MedGenUtils.h"
#include "Logger/Logger/Logger.h"
#include <algorithm>

#define LOCAL_SECTION LOG_MED_UTILS
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

void set_rand_seed(int seed)
{
	if (seed != -1)
		globalRNG::srand(seed);
	else {
		MedTimer t;
		unsigned long long curr = t.get_clock_micro();
		unsigned int s = (unsigned int)(curr & 0xffffffff);
		globalRNG::srand(s);
	}
}

void get_rand_vector_no_repetitions(vector<int> &v, int N, int size)
{
	vector<int> w(N);

	for (int i=0; i<N; i++)
		w[i] = i;

	random_shuffle(w.begin(),w.end());

	v.resize(size);

	for (int i=0; i<size; i++)
		v[i] = w[i];
}

void get_rand_vector_with_repetitions(vector<int> &v, int N, int size)
{
	v.resize(size);
	for (int i=0; i<size; i++)
		v[i] = rand_N(N);
}


void get_rand_splits(vector<int> &split, int nsplits, int size)
{
	split.resize(size);
	for (int i=0; i<size; i++)
		split[i] = i % nsplits;
	random_shuffle(split.begin(), split.end());
}

void categorize_vec(vector<float> &in, vector<float> &bounds, vector<float> &out)
{
	out.resize(in.size());

	for (int i=0; i<in.size(); i++) {
		int j = 0;
		while (j<bounds.size() && in[i]>bounds[j+1]) j++;
		out[i] = (float)j;
	}
}

void get_probs_vec(vector<float> &v)
{
	float sum = 0;

	for (int i=0; i<v.size(); i++) {
		sum += v[i];
	}

	if (sum > 0) {
		for (int i=0; i<v.size(); i++)
			v[i] = v[i]/sum;
	}
}

// Initialization Utility
int initialization_text_to_map(const string& text, map<string, string>& init_map) {

	if (text == "")
		return 0;

	//MLOG("INPUT TEXT: %s\n", text.c_str());
	// dealing with {}
	// whenever there's a v={S} where S is any string (that may also include {}) we want the map to put S for v ...
	// follows is (hence) an ugly code to parse that
	// but this adds the ability to pass parameters for an embedded element within our current (say a model that holds parameters for other models)
	//

	vector<size_t> start_pos;
	vector<size_t> end_pos;
	vector<pair<size_t, size_t>> from_to;

	// find all positions of "={"
	size_t pos = text.find("={", 0);
	while (pos != string::npos) {
		start_pos.push_back(pos);
		pos = text.find("={", pos+1);
	}

	// find all positions of "}"
	pos = text.find("}", 0);
	while (pos != string::npos) {
		end_pos.push_back(pos);
		pos = text.find("}", pos+1);
	}

	// treating nesting 
	if (start_pos.size()>0 && end_pos.size()>0) {

		int i = 0, j = 0, stack = 0, stack_first=-1, stack_last=-1;

		while (j<end_pos.size()) {
			if (i<(int)start_pos.size() && start_pos[i] < end_pos[j]) {
				if (stack_first < 0) stack_first = (int) start_pos[i];
				stack++;
				i++;
			}
			else {
				if (stack == 0) {
					MERR("ERROR: Unmatched {} in string %s\n", text.c_str());
					return -1;
				}
				stack--;
				if (stack == 0) stack_last = (int)end_pos[j];
				j++;
			}

			if (stack == 0) {
				from_to.push_back(pair<size_t, size_t>(stack_first, stack_last));
				stack_first = -1;
			}
		}

		for (auto &ft : from_to) {
			MLOG("found substring: %d-%d : %s\n", ft.first, ft.second, text.substr(ft.first, ft.second-ft.first+1).c_str());
		}

	}


	// replacing {} areas with other strings to allow for correct parsing, and then returning them
	string new_text = "";
	map<string, string> replacers;
	if (from_to.size() == 0) new_text = text;
	else {
		new_text = text.substr(0, from_to[0].first+1); // up to the first '='
		int j;
		for (j=0; j<from_to.size(); j++) {
			string name = "REPLACE_ME_LATER_NUMBER_"+to_string(j);
			string replacer = text.substr(from_to[j].first+2, from_to[j].second-from_to[j].first-2);
			MLOG("replacer %d : %s -> %s\n", j, name.c_str(), replacer.c_str());
			new_text += name;
			replacers[name] = replacer;
			if (j<from_to.size()-1)
				new_text += text.substr(from_to[j].second+1, from_to[j+1].first-from_to[j].second);
		}
		new_text += text.substr(from_to[j-1].second+1, text.length()-from_to[j-1].second);
		MLOG("new_text is %s\n", new_text.c_str());

	}

	// TBD


	// get "Name = value" fields
	vector<string> fields;
	boost::split(fields, new_text, boost::is_any_of(";"));

	// get name + value
	vector<string> sub_fields;

	for (string& field : fields) {
		if (field.size() == 0)
			continue;

		boost::split(sub_fields, field, boost::is_any_of("="));
		if (sub_fields.size() != 2) {
			fprintf(stderr, "Cannot parse \'%s\' from \'%s\'\n", field.c_str(), text.c_str());
			return -1;
		}
		init_map[sub_fields[0]] = sub_fields[1];
	}

	for (auto &el : init_map) {
		if (el.second.compare(0, 24, "REPLACE_ME_LATER_NUMBER_") == 0) {
			init_map[el.first] = replacers[el.second];
		}
	}

	return 0;
}

bool is_windows_os(void) 
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32)
	return true;
#else
	return false;
#endif
}
