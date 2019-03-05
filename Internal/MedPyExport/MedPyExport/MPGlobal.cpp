#include "MPGlobal.h"
#include "MedTime/MedTime/MedTime.h"
#include "MedUtils/MedUtils/MedGlobalRNG.h"


int MPGlobalClass::MEDPY_GET_default_time_unit() { return global_default_time_unit; }
void MPGlobalClass::MEDPY_SET_default_time_unit(int new_val) { global_default_time_unit = new_val; }
int MPGlobalClass::MEDPY_GET_default_windows_time_unit() { return global_default_windows_time_unit; }
void MPGlobalClass::MEDPY_SET_default_windows_time_unit(int new_val) { global_default_windows_time_unit = new_val; }

int MPRNG::rand() { return globalRNG::rand(); };
int MPRNG::rand30() { return globalRNG::rand30(); };
void MPRNG::srand(int val) { globalRNG::srand(val); };
int MPRNG::max() { return globalRNG::max(); };


#define _CRT_SECURE_NO_WARNINGS
#define _CRT_RAND_S

#include "../Tools/MedProcessUtils/CommonLib/commonHeader.h"

void MPCommonLib::print_auc_performance(MPSamples &samples, int nfolds, string outFile) {
	::print_auc_performance(*(samples.o), nfolds, outFile);
}