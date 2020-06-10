#ifndef _FTR__PROCESS_DUPLICATE_H_
#define _FTR__PROCESS_DUPLICATE_H_

#include "FeatureProcess.h"

using namespace std;

/**
* Duplicates the samples in Apply only - can be used for multiple imputations to calculate CI adn more
*/
class DuplicateProcessor : public FeatureProcessor {
public:
	int resample_cnt; ///< how much to resample

	DuplicateProcessor() : FeatureProcessor() { init_defaults(); }

	// Copy
	//void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<PredictorImputer *>(processor)); }

	void init_defaults();
	/// The parsed fields from init command.
	/// @snippet DuplicateProcessor.cpp DuplicateProcessor::init
	int init(map<string, string>& mapper);

	// Apply cleaning model
	int _apply(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	ADD_CLASS_NAME(DuplicateProcessor)
		ADD_SERIALIZATION_FUNCS(processor_type, resample_cnt)
};

MEDSERIALIZE_SUPPORT(DuplicateProcessor)

#endif