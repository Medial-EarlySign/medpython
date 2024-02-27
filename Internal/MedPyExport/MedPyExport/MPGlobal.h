#ifndef __MED__MPGLOBAL__H__
#define __MED__MPGLOBAL__H__

#include "MedPyCommon.h"
#include "MPLogger.h"
#include "MPSamples.h"
#include "MPFeatures.h"
#include "MPPredictor.h"

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
	string version_info;

	MPGlobalClass();
};

class MPCommonLib {
public:
	static void print_auc_performance(MPSamples &samples, MEDPY_NP_INPUT(int* folds, unsigned long long num_folds), string outFile);
	static void shuffleMatrix(MPFeatures& matrix);
	static void read_predictor_from_file(MPPredictor& pred, string predictorFile);
};


#endif //!__MED__MPGLOBAL__H__
