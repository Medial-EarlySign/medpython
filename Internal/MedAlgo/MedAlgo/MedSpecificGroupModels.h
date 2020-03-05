#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>

///wrapper for MedPredictor for certian groups - routes the input to correct model group.
///for example may be used to train specific model for each age group
class MedSpecificGroupModels : public MedPredictor {
public:
	// Model

	int nsamples;
	int nftrs;
	/*double **x;
	double **y;
	float *w;
	*/

	// Function
	MedSpecificGroupModels();
	~MedSpecificGroupModels();

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;
	MedSpecificGroupModels *clone() const;

	ADD_CLASS_NAME(MedSpecificGroupModels)
		ADD_SERIALIZATION_FUNCS(classifier_type, featNum, feat_ths, predictors, model_features, features_count)

		//	void print(FILE *fp, const string& prefix) ;
		// Parameters
		void set_predictors(const vector<MedPredictor *> &predictors); //for each group index
	void set_group_selection(int featNum, const vector<float> &feat_ths);
	MedPredictor *get_model(int ind);
	int model_cnt() const;
private:
	vector<MedPredictor *> predictors;
	int featNum;
	vector<float> feat_ths;
	int selectPredictor(const float *x) const; //retrieve predictor index
};


MEDSERIALIZE_SUPPORT(MedSpecificGroupModels)
