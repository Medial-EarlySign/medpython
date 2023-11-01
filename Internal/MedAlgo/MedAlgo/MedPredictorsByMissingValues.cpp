#include "MedPredictorsByMissingValues.h"

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL LOG_DEF_LEVEL

int MedPredictorsByMissingValues::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "predictor_type")
			predictor_type = it->second;
		else if (it->first == "predictor_params")
			predictor_params = it->second;
		else if (it->first == "masks_params")
			masks_params = it->second;
		else
			MTHROW_AND_ERR("Unsupported argument \"%s\" for MedPredictorsByMissingValues\n", it->first.c_str());
	}
	return 0;
}

int MedPredictorsByMissingValues::learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts) {
	vector<string> masks;
	boost::split(masks, masks_params, boost::is_any_of("|"));

	int p = pow(2, masks.size()); // p - number of models to train

	vector<int> empty_is_all;
	vector<int> cols;
	int len_cols;
	vector<int> cols_to_keep;
	vector<string> mask;
	string sig;
	MedMat<float> x_copy;

	if (masks.size() == 0) {
		MTHROW_AND_ERR("The predictor needs at least one mask");
	}

	// for every model we find the relevant masks by binary coding of model number
	for (int i = 0; i < p; i++) {
		int i1 = i;
		for (int j = masks.size(); j > 0; j--) {
			if (i1 >= pow(2, j-1)) {
				// mask j is part of model i
				boost::split(mask, masks[j-1], boost::is_any_of(","));
				if (mask.size() == 0) {
					MTHROW_AND_ERR("Check the masks - one of the masks is empty");
				}
				len_cols = cols.size();
				// for every signal in the mask, we find all features that should be removed
				for (int m = 0; m < mask.size(); m++) {
					sig = '.' + mask[m] + '.';
					for (int s = 0; s < x.signals.size(); s++) {
						bool isFound = x.signals[s].find(sig) != string::npos;
						if (isFound) {
							cols.insert(cols.end(), s);
						}
					}
				}
				if (len_cols == cols.size()) {
					MTHROW_AND_ERR("Mask \"%s\" removed no features", masks[j-1].c_str());
				}
				i1 = i1 - pow(2, j-1);
			}
		}
		// we keep only the columns that are not removed
		for (int s = 0; s < x.signals.size(); s++)
			if (find(cols.begin(), cols.end(), s) == cols.end()) {
				cols_to_keep.insert(cols_to_keep.end(), s);
			}
		x_copy = x;
		x_copy.get_sub_mat(empty_is_all, cols_to_keep); // retrun sub mat with all the rows and just cols_to_keep columns
		cols_to_keep.clear();
		cols.clear();

		// train model i
		predictors.push_back(MedPredictor::make_predictor(predictor_type, predictor_params));
		predictors[i]->learn(x_copy, y, wgts);
		// x_copy.write_to_csv_file("/nas1/Work/Users/Eitan/LungWithMask/predictions/tmp" + to_string(i) + ".csv"); // write to csv
	}
	return 0;
}
int MedPredictorsByMissingValues::predict(MedMat<float> &x, vector<float> &preds) const {
	throw logic_error("please implement");
}