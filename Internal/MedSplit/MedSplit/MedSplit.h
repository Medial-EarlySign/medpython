#pragma once

#include <vector>
#include <map>
#include <string>

using namespace std;

class MedSplit {
public:
	vector<int> pids;
	vector<int> split;
	int nsplits;
	map<int, int> pid2split;

	MedSplit() { clear(); }
	void clear() { pids.clear(); split.clear(); pid2split.clear(); nsplits = 0; }
	int read_from_file(const string &fname);
	int write_to_file(const string &fname);

	int create_random(vector<int> &in_pids, int nsplits);

};
