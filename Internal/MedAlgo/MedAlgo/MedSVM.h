#pragma once

#include <MedAlgo/MedAlgo/MedAlgo.h>

class MedSvm : public MedPredictor {
public:
	// Model

	struct svm_parameter params;
	struct svm_model *model;
	/*double **x;
	double **y;
	float *w;
	*/

	// Function
	MedSvm();
	MedSvm(void *params);
	MedSvm(struct svm_parameter &params);
	~MedSvm();

	void init_defaults();
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedSvm.cpp MedSvm::init
	virtual int set_params(map<string, string>& mapper);
	int init(struct svm_parameter &params);

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_CLASS_NAME(MedSvm)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

private:


};
