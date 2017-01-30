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

using namespace std;

int init_map_from_string(string text, map<string, string>& map);
int init_dvec(string& in, vector<int>& out);

#endif