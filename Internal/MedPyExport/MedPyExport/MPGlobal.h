#ifndef __MED__MPGLOBAL__H__
#define __MED__MPGLOBAL__H__

#include "MedPyCommon.h"
#include "MPLogger.h"

class MPRNG {
public:
	int rand();
	int rand30();
	void srand(int val);
	int max();
};

class MPGlobalClass {
public:
	int MEDPY_GET_default_time_unit();
	void MEDPY_SET_default_time_unit(int new_val);
	int MEDPY_GET_default_windows_time_unit();
	void MEDPY_SET_default_windows_time_unit(int new_val);
	MPRNG RNG;
	MPLogger logger;
};

#endif //!__MED__MPGLOBAL__H__
