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

int MedPredictorsByMissingValues::Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
	throw logic_error("please implement");
}
int MedPredictorsByMissingValues::Predict(float *x, float *&preds, int nsamples, int nftrs) const {
	throw logic_error("please implement");
}