#ifndef __BIN_SPLITTER_H__
#define __BIN_SPLITTER_H__
#include <map>
#include "PredictiveModel.h"
using namespace std;

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
	SameValueWidth, ///< simple merge by splitting the range with fixed width
	PartitaionMover, ///< using the partition moving to get best clusters
	IterativeMerge, ///< merging iterativaly clusters
	DynamicSplit, ///< splits based on QRF on single feature - dynamic bin size
};

/**
* A specific settings for binning feature
*/
class BinSettings {
public:
	int min_bin_count; ///<minimal count of cases+controls to create bin for feature
	double min_res_value; ///< minimal distance from each feature value between bins. if 0 will not use
	int binCnt; ///< the bin Count for spliting, if 0 will not use

	BinSplitMethod split_method; /// the split method, please reffer to BinSplitMethod
};

namespace medial {
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
	}
}

#endif