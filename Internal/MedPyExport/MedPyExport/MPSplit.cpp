#include "MPSplit.h"
//#include "MedFeat/MedFeat/MedOutcome.h"
#include "MedSplit/MedSplit/MedSplit.h"



MPSplit::MPSplit() { o = new MedSplit(); };
MPSplit::~MPSplit() { delete o; };

int MPSplit::MEDPY_GET_nsplits() { return o->nsplits; };
void MPSplit::clear() { o->clear(); };
int MPSplit::read_from_file(const string &fname) { return o->read_from_file(fname); };
int MPSplit::write_to_file(const string &fname) { return o->write_to_file(fname); };

MPIntIntMapAdaptor MPSplit::MEDPY_GET_pid2split() {
	return MPIntIntMapAdaptor(&(o->pid2split));
}