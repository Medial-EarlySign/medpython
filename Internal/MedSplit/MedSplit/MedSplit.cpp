#include "MedSplit.h"

#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <Logger/Logger/Logger.h>
#include <MedUtils/MedUtils/MedGlobalRNG.h>
#include <MedUtils/MedUtils/MedGenUtils.h>
#include <fstream>
#include <algorithm>
#include <thread>

#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

using namespace std;

//================================================================================================
// MedSplit
//================================================================================================
//................................................................................................
int MedSplit::read_from_file(const string &fname)
{
	clear();

	ifstream inf(fname);

	if (!inf) {
		MERR("MedSplit: can't open file %s for read\n", fname.c_str());
		return -1;
	}
	bool found_nsplits = false;
	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {

			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			vector<string> fields;
			boost::split(fields, curr_line, boost::is_any_of(" \t"));

			if (fields.size() >= 2) {

				if (!found_nsplits && fields[0] == "NSPLITS") {
					nsplits = stoi(fields[1]);
					found_nsplits = true;
				}
				else {
					int pid = stoi(fields[0]);
					int isplit = stoi(fields[1]);

					pids.push_back(pid);
					split.push_back(isplit);
					pid2split[pid] = isplit;
				}

			}
		}
	}
	if (!found_nsplits) {
		MERR("MedSplit: file [%s] does not contain a line with NSPLITS\n", fname.c_str());
		return -1;
	}

	inf.close();
	return 0;
}

//................................................................................................
int MedSplit::write_to_file(const string &fname)
{
	ofstream of;

	of.open(fname, ios::out);

	if (!of) {
		MERR("MedSplit: can't open file %s for write\n", fname.c_str());
		return -1;
	}

	of << "NSPLITS\t" << nsplits << "\n";
	for (int i = 0; i<pids.size(); i++) {
		of << pids[i] << "\t" << split[i] << "\n";
	}

	of.close();
	return 0;
}

//................................................................................................
int MedSplit::create_random(vector<int> &in_pids, int _nsplits)
{
	clear();

	nsplits = _nsplits;
	pids = in_pids;
	get_rand_splits(split, nsplits, (int)pids.size());

	return 0;
}
