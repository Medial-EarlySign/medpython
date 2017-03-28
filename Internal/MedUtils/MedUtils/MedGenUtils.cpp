#include "MedGenUtils.h"
#include "Logger/Logger/Logger.h"
#include <algorithm>


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

	// get "Name = value" fields
	vector<string> fields;
	boost::split(fields, text, boost::is_any_of(";"));

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

	return 0;
}
