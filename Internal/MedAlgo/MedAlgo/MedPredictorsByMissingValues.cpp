#include "MedPredictorsByMissingValues.h"

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL LOG_DEF_LEVEL

int MedPredictorsByMissingValues::init(map<string, string>& mapper) {
	throw logic_error("please implement");
	//! [MedPredictorsByMissingValues::init]
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "predictor_type")
			predictor_type = it->second;
		else if (it->first == "predictor_params")
			predictor_params = it->second;
		else
			MTHROW_AND_ERR("Unsupported argument \"%s\" for MedPredictorsByMissingValues\n", it->first.c_str());
	}
	//! [MedPredictorsByMissingValues::init]
}

int MedPredictorsByMissingValues::learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts) {
	throw logic_error("please implement");
	//x.write_to_csv_file() - write to csv
	//x(i,j); - access sample i in feature j
	//x.signals - name of features
	//x.get_sub_mat() - return subset of matrix

	//predictors.push_back(MedPredictor::make_predictor(predictor_type,predictor_params));
	//predictors[0]->learn(x, y, wgts);
}
int MedPredictorsByMissingValues::predict(MedMat<float> &x, vector<float> &preds) const {
	throw logic_error("please implement");
}