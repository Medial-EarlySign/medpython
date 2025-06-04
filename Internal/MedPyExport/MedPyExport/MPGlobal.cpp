#include "MPGlobal.h"
#include "MedTime/MedTime/MedTime.h"
#include "MedUtils/MedUtils/MedGlobalRNG.h"
#include <MedUtils/MedUtils/MedUtils.h>


int MPGlobalClass::MEDPY_GET_default_time_unit() { return global_default_time_unit; }
void MPGlobalClass::MEDPY_SET_default_time_unit(int new_val) { global_default_time_unit = new_val; }
int MPGlobalClass::MEDPY_GET_default_windows_time_unit() { return global_default_windows_time_unit; }
void MPGlobalClass::MEDPY_SET_default_windows_time_unit(int new_val) { global_default_windows_time_unit = new_val; }

MPGlobalClass::MPGlobalClass() {
	 version_info = medial::get_git_version();
}

int MPRNG::rand() { return globalRNG::rand(); };
int MPRNG::rand30() { return globalRNG::rand30(); };
void MPRNG::srand(int val) { globalRNG::srand(val); };
int MPRNG::max() { return globalRNG::max(); };


#define _CRT_SECURE_NO_WARNINGS
#define _CRT_RAND_S

#include "CommonLib/CommonLib/commonHeader.h"

void MPCommonLib::print_auc_performance(MPSamples &samples, MEDPY_NP_INPUT(int* folds, unsigned long long num_folds), string outFile) {
	vector<int> folds_vec;
	buf_to_vector(folds, num_folds, folds_vec);
	::print_auc_performance(*(samples.o), folds_vec, outFile);
}
void MPCommonLib::shuffleMatrix(MPFeatures& matrix) {
	::shuffleMatrix(*(matrix.o));
}

void MPCommonLib::read_predictor_from_file(MPPredictor& pred, string predictorFile) {
	::read_predictor_from_file(pred.o, predictorFile);
}