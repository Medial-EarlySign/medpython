#ifndef TRAIN_WITH_MISSING_PROCESSOR_H__
#define	TRAIN_WITH_MISSING_PROCESSOR_H__

#include "FeatureProcess.h"

/**
* TrainMissingProcessor:
* Add missing values to the train matrix for the train process.
* Should be first feature_processor before imputations/normalization if exists.
*/
class TrainMissingProcessor : public FeatureProcessor {
public:
	vector<string> selected_tags; ///< the selected tags to activeate on
	float missing_value; ///< missing value 

	//grouping of imputattion - for example handle imputations by signal (or other groups):
	string grouping; ///< grouping file or "BY_SIGNAL" keyword to group by signal or "BY_SIGNAL_CATEG" - for category signal to split by values (aggreagates time windows) or "BY_SIGNAL_CATEG_TREND" - also splitby TRENDS

	int add_new_data; ///< how many new data data points to add for train according to sample masks
	bool sample_masks_with_repeats; ///< Whether or not to sample masks with repeats
	bool uniform_rand; ///< it True will sample masks uniformlly
	bool use_shuffle; ///< if not sampling uniformlly, If true will use shuffle (to speed up runtime)
	int subsample_train; ///< if not zero will use this to subsample original train sampels to this number
	int limit_mask_size; ///< if set will limit mask size in the train - maximal number of missing values

	bool verbose; ///< print verbose

	TrainMissingProcessor() : FeatureProcessor() { init_defaults(); }
	// Copy
	virtual void copy(FeatureProcessor *processor) { *this = *(dynamic_cast<TrainMissingProcessor *>(processor)); }

	/// The parsed fields from init command.
	/// @snippet MultiplierProcessor.cpp MultiplierProcessor::init
	int init(map<string, string>& mapper);
	void init_defaults();

	//print function
	void dprint(const string &pref, int fp_flag);

	// Apply do nothing in apply
	int _apply(MedFeatures& features, unordered_set<int>& ids) { return 0; }
	int Learn(MedFeatures& features, unordered_set<int>& ids);

	// Serialization
	ADD_CLASS_NAME(TrainMissingProcessor)
		ADD_SERIALIZATION_FUNCS(processor_type, selected_tags, missing_value, add_new_data, sample_masks_with_repeats,
			uniform_rand, use_shuffle, subsample_train, limit_mask_size, grouping, groupNames, group2Inds, verbose)
private:
	vector<vector<int>> group2Inds;
	vector<string> groupNames;
};

MEDSERIALIZE_SUPPORT(TrainMissingProcessor)

#endif // !TRAIN_WITH_MISSING_PROCESSOR_H__

