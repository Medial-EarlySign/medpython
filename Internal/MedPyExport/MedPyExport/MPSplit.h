#ifndef __MED__MPSPLIT__H__
#define __MED__MPSPLIT__H__

#include "MedPyCommon.h"

class MedSplit;

class MPSplit {
public:
	MEDPY_IGNORE(MedSplit* o);
	MPSplit();
	~MPSplit();
	//vector<int> pids;
	//vector<int> split;
	
	int MEDPY_GET_nsplits();

	MPIntIntMapAdaptor MEDPY_GET_pid2split();
	//map<int, int> pid2split;

	//MedSplit() { clear(); }
	void clear();
	int read_from_file(const string &fname);
	int write_to_file(const string &fname);

	//int create_random(vector<int> &in_pids, int nsplits);

};

#endif // !__MED__MPSPLIT__H__
