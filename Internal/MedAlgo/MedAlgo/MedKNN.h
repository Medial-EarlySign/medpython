#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>


//======================================================================================
// KNN
//======================================================================================
typedef enum {
	KNN_DIST_MEAN,
	KNN_1_DIST,
	KNN_WEIGHTEDLS,
	KNN_AVG_LAST
} knnAveraging;

typedef enum {
	KNN_L1,
	KNN_L2,
	KNN_METRIC_LAST
}knnMetric;

struct MedKNNParams : public SerializableObject {

	int k;
	knnAveraging knnAv;
	knnMetric knnMetr;

	ADD_CLASS_NAME(MedKNNParams)
		ADD_SERIALIZATION_FUNCS(k, knnAv, knnMetr)
};

class MedKNN : public MedPredictor {
public:
	// Model
	int nsamples;
	int nftrs;
	vector<float> x;
	vector<float> y;
	vector<float> w;


	// Parameters
	MedKNNParams params;


	// Function
	MedKNN();
	MedKNN(void *params);
	MedKNN(MedKNNParams& params);
	/// The parsed fields from init command.
	/// @snippet MedKNN.cpp MedKNN::init
	virtual int set_params(map<string, string>& mapper);
	int init(void *params);
	knnAveraging get_knn_averaging(string name);
	knnMetric get_knn_metric(string name);

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_CLASS_NAME(MedKNN)
		ADD_SERIALIZATION_FUNCS(classifier_type, params, nsamples, nftrs, x, y, w)
};

MEDSERIALIZE_SUPPORT(MedKNNParams)
MEDSERIALIZE_SUPPORT(MedKNN)
