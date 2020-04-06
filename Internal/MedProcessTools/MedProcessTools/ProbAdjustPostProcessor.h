#ifndef __PROBADJUSTPOSTPROCESSOR_H__
#define __PROBADJUSTPOSTPROCESSOR_H__

#include <MedAlgo/MedAlgo/MedAlgo.h>
#include "PostProcessor.h"
#include <Logger/Logger/Logger.h>

class ProbAdjustPostProcessor : public PostProcessor {
public:
	vector<string> names;
	vector<string> resolvedNames;
	vector<int> min, max, factors;
	vector<float> probs;
	vector<float> odds ; ///< over all odds. learn if not given

	/// Parameters
	string priorsFile;

	// Functions
	ProbAdjustPostProcessor() { processor_type = PostProcessorTypes::FTR_POSTPROCESS_ADJUST; };

	~ProbAdjustPostProcessor() {};

	/// Global init for general args in all explainers
	int init(map<string, string> &mapper);

	///Learns from predictor and train_matrix (PostProcessor API)
	void Learn(const MedFeatures &matrix);
	void Apply(MedFeatures &matrix) const;

	// Helper functions
	void readPriors();
	void getOdds(const MedFeatures &matrix);

	ADD_CLASS_NAME(ProbAdjustPostProcessor)
	ADD_SERIALIZATION_FUNCS(processor_type, names,resolvedNames,min,max,factors,probs, priorsFile, odds)
};

MEDSERIALIZE_SUPPORT(ProbAdjustPostProcessor)

#endif