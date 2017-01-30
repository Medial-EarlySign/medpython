#ifndef _PID_CLN_REC_H_
#define _PID_CLN_REC_H_

#include "InfraMed/InfraMed/MedPidRepository.h"

class PidClnRec {
public:
	PidClnRec() {};
	~PidClnRec() {};

	// Initialize
	void init_all(PidRec in_rep, int ntime_points);

	// Initialize signal from PidRec
	void init_signal(PidRec& in_rep, int signalId);

	// Initialize from time-point
	void init_signal(int origIdx, int newIdx, int signalId);

	// Change value
	void remove(int signalId, int clnTimeIndex, int valueIndex);
	void change(int signalId, int clnTimeIndex, int valueIndex, float newValue);

	// Get
	void *get(string &sig_name, int idx, int &len);
	void *get(int sid, int idx, int &len);
};

#endif