#ifndef __BIN_SPLITTER_H__
#define __BIN_SPLITTER_H__
#include <map>
#include <vector>
#include <limits.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>

using namespace std;

/**
* a bin split using optimizer
*/
class BinSplitOptimizer {
public:
	vector<int> SplitToBins(const vector<float> &vec, const vector<float> &y, size_t kBins,
		size_t min_samples, vector<float> &histAvg, vector<float> &partitionValues, bool &has_error);
private:
	bool evaluate(const vector<int> &indexes, float &totScore);

	float *_sortedArr;
	int _minSamples;
	map<float, int> _histElements;
	map<float, float> _histAvg;

};

enum BinSplitMethod {
	SameValueWidth, ///< simple merge by splitting the range with fixed width. called "same_width"
	PartitaionMover, ///< using the partition moving to get best clusters. called "partition_mover"
	IterativeMerge, ///< merging iterativaly clusters. called "iterative_merge"
	DynamicSplit, ///< splits based on QRF on single feature - dynamic bin size. called "dynamic_split"
};

/**
* A specific settings for binning feature
*/
class BinSettings : public SerializableObject {
public:
	int min_bin_count; ///<minimal count of cases+controls to create bin for feature
	double min_res_value; ///< minimal distance from each feature value between bins. if 0 will not use
	int binCnt; ///< the bin Count for spliting, if 0 will not use
	float min_value_cutoff; ///< a minimal value trim cutoff
	float max_value_cutoff; ///< a maximal value trim cutoff
	bool weighted; ///< if true in fixed width will average the value in each bin

	BinSplitMethod split_method; /// the split method, please reffer to BinSplitMethod

	int init(map<string, string>& map);

	static const unordered_map<int, string> name_to_method;
	static BinSplitMethod bin_method_name_to_type(const string& bin_method);

	BinSettings() {
		min_bin_count = 1;
		min_res_value = 0;
		binCnt = 10;
		min_value_cutoff = (float)INT_MIN;
		max_value_cutoff = (float)INT_MAX;
		weighted = false;
	}

	ADD_CLASS_NAME(BinSettings)
	ADD_SERIALIZATION_FUNCS(min_bin_count, min_res_value, binCnt, min_value_cutoff, max_value_cutoff, 
		split_method, weighted)
};

MEDSERIALIZE_SUPPORT(BinSettings)
/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	*  \brief process namespace
	*/
	namespace process {
		/// <summary>
		/// splits feature to bin using setting
		/// @param setting the settings of split
		/// @param feature the feature vector values
		/// @param sel_indexes the indexes to take from feature. if empty will take all feature vector 
		/// @param y labels if we have. some binning methods uses the labels for better split
		/// </summary>
		/// <returns>
		///it updates feature to the binned values
		/// </returns>
		void split_feature_to_bins(const BinSettings &setting, vector<float> &feature,
			const vector<int> &sel_indexes, vector<float> &y);
		/// <summary>
		/// normalize feature to be between [0-1] but also change the distribution of values to
		/// be uniform
		/// </summary>
		void normalize_feature_to_uniform(const BinSettings &setting, vector<float> &feature);
	}
}

#endif