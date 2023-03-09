//
// templated code for MedUtils class, included after class definition
//
#ifndef __MED_UTILS_IMP__H_
#define __MED_UTILS_IMP__H_
#include <math.h>

// Discretization of values
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins) {
	return discretize(x, binned_x, nbins, max_bins, MED_MAT_MISSING_VALUE, BIN_EQUISIZE);
}

template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, float missing_value) {
	return discretize(x, binned_x, nbins, max_bins, missing_value, BIN_EQUISIZE);
}

template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, MedBinningType binning) {
	return discretize(x, binned_x, nbins, max_bins, MED_MAT_MISSING_VALUE, binning);
}

template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, float missing_value, MedBinningType binning) {
	if (max_bins <= 0)
		HMTHROW_AND_ERR("Error max_bins should be > 0\n");

	binned_x.clear();
	binned_x.reserve(x.size());

	if (binning >= BIN_LAST) {
		MEDLOG(LOG_MED_UTILS, MAX_LOG_LEVEL, "Unknown binning type %d\n", binning);
		return -1;
	}

	map<float, int> x_values;
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] != missing_value)
			x_values[(float)x[i]]++;
	}

	nbins = (int)x_values.size();
	map<float, int> x_index;

	if (nbins <= max_bins) { // Leave as is
		int idx = 0;
		for (auto it = x_values.begin(); it != x_values.end(); it++)
			x_index[it->first] = idx++;
		assert(idx == nbins);
	}
	else { // Need to combine values into bins
		if (binning == BIN_EQUIDIST) {
			float min_val = x_values.begin()->first;
			float max_val = x_values.rbegin()->first;
			float bin_size = (max_val - min_val) / max_bins;
			if (bin_size > 0)
				for (auto it = x_values.begin(); it != x_values.end(); it++)
					x_index[it->first] = (int)((it->first - min_val) / bin_size);
		}
		else if (binning == BIN_EQUISIZE) {
			int tot = 0;
			for (auto it = x_values.begin(); it != x_values.end(); it++)
				tot += it->second;
			int bin_size = tot / max_bins;

			if (bin_size > 0) {
				tot = 0;
				for (auto it = x_values.begin(); it != x_values.end(); it++) {
					x_index[it->first] = tot / bin_size;
					tot += it->second;
				}
			}
		}

		nbins = max_bins;
	}


	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] == missing_value)
			binned_x.push_back(-1);
		else
			binned_x.push_back(x_index[x[i]]);
	}

	return 0;
}
#endif