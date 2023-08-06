// Various utilities used in MedProcessTools

#ifndef _MED_PROCESS_UTILS_H_
#define _MED_PROCESS_UTILS_H_

#define MAX_NAME_LEN 256
#define PID_REC_SIZE	 5000000
#define DYNAMIC_REC_SIZE 20000000

#include <stdlib.h>
#include <string>
#include <map>
#include <vector>
#include <unordered_set>

class RepProcessor;
class FeatureGenerator;

using namespace std;

int init_dvec(string& in, vector<int>& out);
void get_single_val_from_init_string(string init_s, string attr, string &val_s);

string int_to_string_digits(int i, int ndigits);

// Required signals
void get_all_required_signal_names(unordered_set<string>& signalNames, const vector<RepProcessor *>& rep_processors, int position, vector<FeatureGenerator *>& generators);
void get_all_required_signal_ids(unordered_set<int>& signalIds, vector<RepProcessor *>& rep_processors, int position, vector<FeatureGenerator *>& generators);
void handle_required_signals(vector<RepProcessor *>& processors, vector<FeatureGenerator *>& generators, unordered_set<int>& extra_req_signal_ids,
	vector<int>& all_req_signal_ids_v, vector<unordered_set<int> >& current_required_signal_ids);

// Feature names
int find_in_feature_names(const vector<string>& names,const string& substr, bool throw_on_error = true);

#endif