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

struct MedKNNParams {

	int k;
	knnAveraging knnAv;
	knnMetric knnMetr;
};

class MedKNN : public MedPredictor {
public:
	// Model
	int nsamples;
	int nftrs;
	float *x;
	float *y;
	float *w;


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
	~MedKNN();
	knnAveraging get_knn_averaging(string name);
	knnMetric get_knn_metric(string name);

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_CLASS_NAME(MedKNN)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};


MEDSERIALIZE_SUPPORT(MedKNN)
