#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>

//======================================================================================
// GBM: C++ version of GBM from R.
//======================================================================================
#define MED_GBM_DEF_SHRINKAGE 0.1
#define MED_GBM_DEF_BAG_P 0.5
#define MED_GBM_DEF_NTREES 100
#define MED_GBM_DEF_DEPTH 6
#define MED_GBM_DEF_TAKE_ALL_POS false
#define MED_GBM_DEF_MIN_OBS_IN_NODE 10

typedef gbm_parameters MedGBMParams;

class MedGBM : public MedPredictor {
public:
	/// Loss function
	GBM_LossFunctions loss_function;

	/// Alpha for quantile loss function
	double alpha_quantile;

	/// Model 
	full_gbm_learn_info_t gbm_model;

	/// Parameters
	MedGBMParams params;

	/// Predicting on subset of trees
	int predict_ntrees;

	// Function
	MedGBM();
	MedGBM(void *params);
	MedGBM(MedGBMParams& params);
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedGBM.cpp MedGBM::init
	virtual int set_params(map<string, string>& mapper);
	void init_defaults();
	GBM_LossFunctions get_loss_function(string name);
	~MedGBM();

	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);

	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	ADD_CLASS_NAME(MedGBM)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Print
	void print(FILE *fp, const string& prefix, int level = 0) const;
};

// Initialization of parameters
void init_default_gbm_params(MedGBMParams& _params);

MEDSERIALIZE_SUPPORT(MedGBM)
