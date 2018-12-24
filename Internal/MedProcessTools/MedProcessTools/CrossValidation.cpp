#include "CrossValidator.h"
//#include <MedStat/MedStat/MedPerformance.h>
#include <MedUtils/MedUtils/MedGlobalRNG.h>

#define LOCAL_SECTION LOG_CV
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// CrossValidator
//=======================================================================================

//.......................................................................................
int CrossValidator::doCV(MedPidRepository& rep, MedSamples& samples, int nfolds, MedSamples& outSamples) {

	// Create folds
	map<int, int> folds;
	for (auto& sample : samples.idSamples)
		folds[sample.id] = globalRNG::rand() % nfolds;

	// CV
	return doCV(rep, samples, nfolds, folds, outSamples);

}

//.......................................................................................
int CrossValidator::doCV(MedPidRepository& rep, MedSamples& samples, int nfolds, map<int, int>& folds, MedSamples& outSamples) {

	// CV
	for (int ifold = 0; ifold < nfolds; ifold++) {

		// Learning and test set
		MedSamples learningSamples, testSamples;
		for (auto& sample : samples.idSamples) {
			if (folds[sample.id] != ifold)
				learningSamples.idSamples.push_back(sample);
			else {
				testSamples.idSamples.push_back(sample);
				testSamples.idSamples.back().split = ifold;
				for (auto& sample : testSamples.idSamples.back().samples)
					sample.split = ifold;
			}
		}

		// Filter
		for (int ifilter = 0; ifilter < learning_set_filters.size(); ifilter++) 
			learning_set_filters[ifilter]->filter(rep, learningSamples);

		for (int ifilter = 0; ifilter < test_set_filters.size(); ifilter++)
			test_set_filters[ifilter]->filter(rep,testSamples);

		// Learn Model
		fprintf(stderr, "Fold %d - learning :: learn set %d ids : test set %d ids\n", ifold, (int)learningSamples.idSamples.size(), (int)testSamples.idSamples.size());
		if (model->learn(rep, &learningSamples) < 0) {
			MERR("Learning model for fold %d failed\n", ifold);
			return -1;
		}

		// Apply
		fprintf(stderr, "Fold %d - applying\n", ifold);
		if (model->apply(rep, testSamples) < 0) {
			MERR("Applyinh model for for %d failed\n", ifold);
			return -1;
		}

		// Append with
		outSamples.append(testSamples);
	}

	return 0;
}