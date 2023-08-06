#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>

//======================================================================================
// BackProp 
//======================================================================================
typedef enum { SIGMOID, RELU, LINEAR }neuronFunT;

typedef struct {
	int layerIndex;
	neuronFunT neuronFunction;
	int x;
	int y;
	int  firstWeight, lastWeight;
	int firstSource, lastSource;
	double value;
	double error;
	double delta;
}neuronStruct;

//================================================================
typedef struct {
	neuronStruct *neuron;


	int *source;
	double *weight;
	int numLayers;
	int numNeurons;
	int numInputs, numOutputs;
	int numWeights, numSource;

}netStruct;

struct MedBPParams {

	int numLayers;

	int numIterations;
	double alpha; ///< learning rate
	double beta;///< parameter of logistic function
};

class MedBP : public MedPredictor {
public:
	// Model

	int nsamples;
	int nftrs;
	/*double **x;
	double **y;
	float *w;
	*/

	// Function
	MedBP();
	MedBP(MedBPParams& params);
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedBP.cpp MedBP::init
	virtual int set_params(map<string, string>& mapper);
	~MedBP();

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_CLASS_NAME(MedBP)
		size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Parameters
private:
	MedBPParams params;

	netStruct network;
};
MEDSERIALIZE_SUPPORT(MedBP)