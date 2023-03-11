#define _CRT_SECURE_NO_WARNINGS
#include "MedUtils/MedUtils/MedGlobalRNG.h"
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedAlgo/MedAlgo/MedBP.h>

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

// Backprop Model
//Using NeuralNet.lib
MedBP::MedBP() {
	classifier_type = MODEL_BP;
	transpose_for_learn = false;
	transpose_for_predict = false;

	MedBPParams inParams;
	inParams.numLayers = 3;
	inParams.alpha = 1.e-3;
	inParams.numIterations = 1000;
	inParams.beta = 1.;

	init(&inParams);

}


MedBP::MedBP(MedBPParams& _in_params) {

	classifier_type = MODEL_BP;
	transpose_for_learn = false;
	transpose_for_predict = false;

	init((void *)&_in_params);
}
/*
MedBP::MedBP(void *_in_params) {

	classifier_type = MODEL_KNN ;
	transpose_for_learn = false ;
	transpose_for_predict = false ;



	init(_in_params);
}
*/

//..............................................................................
int MedBP::set_params(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [MedBP::init]
		if (field == "numLayers") params.numLayers = med_stoi(entry.second);
		else if (field == "numIterations") params.numIterations = med_stoi(entry.second);
		else if (field == "alpha") params.alpha = stod(entry.second);
		else if (field == "beta") params.beta = stod(entry.second);
		else MLOG("Unknonw parameter \'%s\' for BP\n", field.c_str());
		//! [MedBP::init]
	}

	return 0;
}

int MedBP::init(void *_in_params)
{

	MedBPParams in_params = *(MedBPParams *)_in_params;

	params = in_params;
	assert(in_params.numLayers > 0);
	assert(in_params.alpha >= 0 && in_params.beta > 0 && in_params.numIterations >= 0);//0 build a net and dont train
	/*x=y=NULL;
	w=NULL;*/
	nftrs = nsamples = 0;
	transpose_for_learn = false;
	transpose_for_predict = false;
	normalize_for_learn = true;
	normalize_for_predict = true;
	normalize_y_for_learn = false;

	memset(&network, 0, sizeof(network));
	return 0;
}


//prototypes for the working functions
#define SQR(x) ((x)*(x))



netStruct  buildNet(int numLayers, int shared[], int xd[], int yd[], int zd[], int size2[][2], int skip2[][2]);
netStruct  buildNet(int numLayers, int shared[], int xd[], int yd[], int zd[], int size2[][2], int skip2[][2], neuronFunT layerNeuronFun[]);// an overload that specifies for each later the neurons type
void  feedForward(netStruct myNet, double *inputs, double *outputs, double beta);
double trainOnce(netStruct *myNet, double *inputs, double * desiredOutputs, double *actualOutputs, double alpha, double beta);
double trainOnceLog(netStruct *myNet, double *inputs, double * desiredOutputs, double *actualOutputs, double alpha, double beta, int sync, double *weightChange);
double trainBatch(netStruct *myNet, double **inputs, double **outputs, int numSamples, double alpha, double beta, int sync);
double trainBatch(netStruct *myNet, double *inputs, double *outputs, int numSamples, double alpha, double beta, int sync);//overload for continous input matrix
void printNet(netStruct thisNet);
netStruct  loadNetwork(char *fname);
void saveNetwork(netStruct thisNet, char *fname);
void saveReceptiveFields(netStruct myNet, char *fname);
double trainRP(netStruct *myNet, double **inputs, double **outputs, int numSamples, double alpha, double beta, double *errorArray, double *sumErrorArray);
void randperm(int n, int perm[]);
void destroyNet(netStruct& net);
void serializeNetwork(void *ptr, netStruct network, size_t *netSize);
netStruct deserializeNetwork(unsigned char *blob);
size_t networkGetSize(netStruct network); // size of serialized network


int MedBP::Learn(float *_x, float *_y, const float *_w/*ignored*/, int _nsamples, int _nftrs)
{
	int(*size2)[2];
	int(*skip2)[2];
	double thisError = 0;
	int *xd, *yd, *zd;
	int *shared;
	double **x, **y;

	xd = (int *)malloc(params.numLayers * sizeof(int));
	yd = (int *)malloc(params.numLayers * sizeof(int));
	zd = (int *)malloc(params.numLayers * sizeof(int));
	assert(xd&&yd&&zd);
	nsamples = _nsamples;
	nftrs = _nftrs;
	xd[0] = nftrs;
	for (int layer = 0; layer < params.numLayers; layer++) {
		yd[layer] = zd[layer] = 1;
		if (layer > 0) xd[layer] = (int)(2 * sqrt(nftrs));
		if (layer == params.numLayers - 1) xd[layer] = 1;
	}
	// make a mtrix of data  in double format
	x = (double**)malloc(sizeof(double *)*nsamples);
	y = (double**)malloc(sizeof(double *)*nsamples);
	assert(x&&y);
	for (int row = 0; row < nsamples; row++) {
		x[row] = (double *)malloc(sizeof(double)*nftrs);
		y[row] = (double *)malloc(sizeof(double)); // single output
		assert(x[row] && y[row]);
		y[row][0] = _y[row] > 0;
		for (int col = 0; col < nftrs; col++)x[row][col] = _x[col + nftrs * row];
	}


	shared = (int *)malloc(sizeof(int)*params.numLayers);


	assert(shared);
	for (int i = 0; i < params.numLayers; i++)
		shared[i] = 0;
	size2 = (int(*)[2])malloc(2 * sizeof(int)*params.numLayers);
	skip2 = (int(*)[2])malloc(2 * sizeof(int)*params.numLayers);
	assert(size2);
	assert(skip2);
	for (int i = 1; i < params.numLayers; i++) {
		skip2[i][0] = 0;
		size2[i][0] = xd[i - 1];// we do not have receptive fields here
		skip2[i][1] = 0;
		size2[i][1] = 1;
	}

	destroyNet(network);
	network = buildNet(params.numLayers, shared, xd, yd, zd, size2, skip2);
	for (int i = 0; i < params.numIterations; i++) {
		thisError = trainBatch(&network, x, y, nsamples, params.alpha, params.beta, 0);

		if (i == 0)
			fprintf(stderr, "end %d  %lf\n", i, thisError);
		fflush(stderr);
	}
	fprintf(stderr, "end   %lf\n", thisError);
	for (int row = 0; row < nsamples; row++) {
		if (x[row])free(x[row]);
		x[row] = 0;
		if (y[row])free(y[row]);
		y[row] = 0;
	}
	free(x); x = 0;
	free(y); y = 0;
	free(skip2); skip2 = 0;
	free(size2); size2 = 0;
	free(xd); xd = 0;
	free(yd); yd = 0;
	free(zd); zd = 0;
	free(shared); shared = 0;
	return(0);
}



int MedBP::Predict(float *xPred, float *&preds, int pred_samples, int _nftrs) const {
	assert(preds);
	assert(_nftrs == nftrs);

	// OK, lets go ...
	//fprintf(stderr,"Running BackProp : numLayers = %d , Data = (%d + %d) x %d\n",params.numLayers,nsamples,nsamples,nftrs) ;




	double *input;
	double result;
	input = (double *)malloc(sizeof(double)*nftrs);
	assert(input);
	for (int i = 0; i < pred_samples; i++) {
		for (int ii = 0; ii < nftrs; ii++)input[ii] = xPred[nftrs*i + ii];
		feedForward(network, input, &result, params.beta);
		preds[i] = (float)result;
	}


	free(input);

	return(0);
}


size_t MedBP::get_size() {

	return(sizeof(*this) + MedSerialize::get_size(classifier_type) + networkGetSize(network));
}

size_t MedBP::serialize(unsigned char *blob) {
	//assumes blob already assined to get_size()
	size_t ptr = 0;

	ptr += MedSerialize::serialize(blob, classifier_type);
	size_t advance;
	memcpy(blob + ptr, this, advance = sizeof(*this)); ptr += advance;

	size_t netSize;
	serializeNetwork(blob + ptr, network, &netSize);
	ptr += netSize;
	return (ptr);
}

size_t MedBP::deserialize(unsigned char *blob) {

	size_t ptr = 0;
	ptr += MedSerialize::deserialize(blob, classifier_type);
	memcpy(this, blob + ptr, sizeof(*this));
	this->network = deserializeNetwork(blob + ptr + sizeof(*this));
	return ptr + sizeof(*this);
}






//#include "neuralNetLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "MedAlgo/MedAlgo/MedAlgo.h"


#define EPS (1.e-3) /* bind accuracy to avoid singularity */
#define ALPHA (0.01)// RELU parameter

netStruct  buildNet(int numLayers, int shared[], int xd[], int yd[], int zd[], int size2[][2], int skip2[][2])
{
	netStruct myNet;
	int i, numWeights;
	int neuronPointer = 0, weightPointer = 0, sourcePointer = 0;
	int layerIndex, xNeuron, yNeuron, zNeuron;
	int xw, yw;
	myNet.numNeurons = 1;// bias
	myNet.numInputs = xd[0] * yd[0];
	myNet.numOutputs = xd[numLayers - 1] * yd[numLayers - 1];
	for (i = 0; i < numLayers; i++)
		myNet.numNeurons += xd[i] * yd[i] * zd[i];

	numWeights = myNet.numNeurons*myNet.numNeurons;
	// allocations
	myNet.neuron = (neuronStruct *)malloc(sizeof(neuronStruct)*(myNet.numNeurons + 1));

	myNet.weight = (double *)malloc(sizeof(double)*numWeights);
	myNet.source = (int *)malloc(sizeof(int)*numWeights);

	myNet.numLayers = numLayers;

	//bias
	myNet.neuron[0].value = 1;
	myNet.neuron[0].layerIndex = -1;
	myNet.neuron[0].x = myNet.neuron[0].y = 0;
	myNet.neuron[neuronPointer].firstWeight = 0;
	myNet.neuron[neuronPointer].lastWeight = 0;
	myNet.neuron[neuronPointer].firstSource = 0;
	myNet.neuron[neuronPointer].lastSource = 0;
	neuronPointer++;
	myNet.numLayers = numLayers;


	for (layerIndex = 0; layerIndex < numLayers; layerIndex++) {
		for (zNeuron = 0; zNeuron < zd[layerIndex]; zNeuron++) {//no sharing along z
			int layerWeightStart = weightPointer; //return to it when weights are shared
			for (xNeuron = 0; xNeuron < xd[layerIndex]; xNeuron++)
				for (yNeuron = 0; yNeuron < yd[layerIndex]; yNeuron++) {
					if (shared[layerIndex])
						weightPointer = layerWeightStart;//shared weights use the same weights but keep different sources

					myNet.neuron[neuronPointer].x = xNeuron;
					myNet.neuron[neuronPointer].y = yNeuron;
					myNet.neuron[neuronPointer].firstWeight = weightPointer;
					myNet.neuron[neuronPointer].firstSource = sourcePointer;
					myNet.neuron[neuronPointer].layerIndex = layerIndex;
					myNet.neuron[neuronPointer].neuronFunction = SIGMOID;
					if (layerIndex == 0) {
						myNet.neuron[neuronPointer].lastWeight = weightPointer;
						myNet.neuron[neuronPointer].lastSource = sourcePointer;
						neuronPointer++;
						continue;
					}
					//set the bias weight;
					myNet.source[sourcePointer++] = 0;
					unsigned int urand;
					urand = globalRNG::rand();
					myNet.weight[weightPointer++] = (double)urand / UINT_MAX * 2 - 1;

					for (xw = skip2[layerIndex][0] * xNeuron; xw < skip2[layerIndex][0] * xNeuron + size2[layerIndex][0]; xw++)
						for (yw = skip2[layerIndex][1] * yNeuron; yw < skip2[layerIndex][1] * yNeuron + size2[layerIndex][1]; yw++)
							for (int nIndex = 0; nIndex < neuronPointer; nIndex++)
								if (myNet.neuron[nIndex].layerIndex == layerIndex - 1)
									if (myNet.neuron[nIndex].x == xw)
										if (myNet.neuron[nIndex].y == yw) {
											myNet.source[sourcePointer] = nIndex;

											sourcePointer++;
											urand = globalRNG::rand();
											myNet.weight[weightPointer++] = (double)urand / UINT_MAX * 2 - 1;

										}
					myNet.neuron[neuronPointer].lastSource = sourcePointer;
					myNet.neuron[neuronPointer].lastWeight = weightPointer;
					neuronPointer++;
				}
		}//neuron


	}//layer

	myNet.numWeights = weightPointer;
	myNet.numSource = sourcePointer;
	myNet.weight = (double *)realloc(myNet.weight, weightPointer * sizeof(*(myNet.weight)));
	myNet.source = (int *)realloc(myNet.source, sourcePointer * sizeof(int));
	return(myNet);
}
double trainBatch(netStruct *myNet, double *inputs, double *outputs, int numSamples, double alpha, double beta, int sync) {
	int numInputs = myNet->numInputs;
	int numOutputs = myNet->numOutputs;
	double **inputP = (double **)malloc(numSamples * sizeof(double *));
	assert(inputP);
	double **outputP = (double **)malloc(numSamples * sizeof(double*));
	assert(outputP);


	for (int k = 0; k < numSamples; k++) {
		inputP[k] = inputs + numInputs * k;
		outputP[k] = outputs + numOutputs * k;
	}
	return(
		trainBatch(myNet, inputP, outputP, numSamples, alpha, beta, sync)
		);



}

netStruct  buildNet(int numLayers, int shared[], int xd[], int yd[], int zd[], int size2[][2], int skip2[][2], neuronFunT layerNeuronFun[]) // an overload that specifies for each later the neurons type
{
	netStruct myNet = buildNet(numLayers, shared, xd, yd, zd, size2, skip2);
	for (int i = 0; i < myNet.numNeurons; i++)
		myNet.neuron[i].neuronFunction = layerNeuronFun[myNet.neuron[i].layerIndex];

	return(myNet);

}
void  feedForward(netStruct myNet, double *inputs, double *outputs, double beta)
// allocation of outputs is done by caller.
{


	myNet.neuron[0].value = 1;
	for (int i = 1; i <= myNet.numInputs; i++)
		myNet.neuron[i].value = inputs[i - 1];

	for (int i = myNet.numInputs + 1; i < myNet.numNeurons; i++) {
		myNet.neuron[i].value = 0;
		int weightIndex = myNet.neuron[i].firstWeight;
		for (int source = myNet.neuron[i].firstSource; source < myNet.neuron[i].lastSource; source++)
			myNet.neuron[i].value += myNet.neuron[myNet.source[source]].value*myNet.weight[weightIndex++];
		switch (myNet.neuron[i].neuronFunction) {
		case SIGMOID:	if ((-beta * myNet.neuron[i].value) > -log(EPS))
			myNet.neuron[i].value = EPS * 2 - 1;
						else if ((-beta * myNet.neuron[i].value) < log(EPS))
							myNet.neuron[i].value = -EPS * 2 + 1;
						else
							myNet.neuron[i].value = (1. / (1 + exp(-beta * myNet.neuron[i].value))) * 2 - 1;

			break;
		case RELU: if (myNet.neuron[i].value < 0) myNet.neuron[i].value *= ALPHA;// leaky RELU positive stays the same. negative multiplied by alpha
			break;
		case LINEAR: break;// in linear value stays as is
		}

	}

	for (int i = 0; i < myNet.numOutputs; i++)
		outputs[i] = myNet.neuron[myNet.numNeurons - myNet.numOutputs + i].value;

}
#define COMPENSATION (1)  // prefer positive outputs because they are rare
double trainOnce(netStruct *myNet, double *inputs, double * desiredOutputs, double *actualOutputs, double alpha, double beta)
{
	int neuronIndex;
	double sumError = 0;
	feedForward(*myNet, inputs, actualOutputs, beta);
	//zero deltas
	for (int i = 0; i < myNet->numNeurons; i++)
		myNet->neuron[i].delta = 0;
	//set the delta at output layer and move it backwards

	for (int i = 0; i < myNet->numOutputs; i++) {
		neuronIndex = myNet->numNeurons - myNet->numOutputs + i;
		myNet->neuron[neuronIndex].error = (actualOutputs[i] - desiredOutputs[i]);
		if (myNet->neuron[neuronIndex].error < 0)myNet->neuron[neuronIndex].error *= COMPENSATION;
		sumError += SQR(myNet->neuron[neuronIndex].error);
		switch (myNet->neuron[neuronIndex].neuronFunction) {
		case SIGMOID:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].error*beta * 2 * (myNet->neuron[neuronIndex].value + 1)*(1 - myNet->neuron[neuronIndex].value);
			break;
		case RELU:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].error;
			if (myNet->neuron[neuronIndex].value < 0)myNet->neuron[neuronIndex].delta *= ALPHA;
			break;
		case LINEAR:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].error;
			break;

		}

		int weightIndex = myNet->neuron[neuronIndex].firstWeight;
		for (int source = myNet->neuron[neuronIndex].firstSource; source < myNet->neuron[neuronIndex].lastSource; source++)
			myNet->neuron[myNet->source[source]].delta += myNet->weight[weightIndex++] * myNet->neuron[neuronIndex].delta;


	}


	//set the delta at non outputlayer and move it backwards
	for (neuronIndex = myNet->numNeurons - 1 - myNet->numOutputs; neuronIndex >= 0; neuronIndex--) {
		switch (myNet->neuron[neuronIndex].neuronFunction) {
		case SIGMOID:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].delta*beta * 2 * (myNet->neuron[neuronIndex].value + 1)*(1 - myNet->neuron[neuronIndex].value);
			break;
		case RELU:
			if (myNet->neuron[neuronIndex].value < 0)myNet->neuron[neuronIndex].delta *= ALPHA;
			break;
		case LINEAR:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].error;
			break;
		}
		int weightIndex = myNet->neuron[neuronIndex].firstWeight;
		for (int source = myNet->neuron[neuronIndex].firstSource; source < myNet->neuron[neuronIndex].lastSource; source++)
			myNet->neuron[myNet->source[source]].delta += myNet->weight[weightIndex++] * myNet->neuron[neuronIndex].delta;
	}

	//adjust weights
	for (neuronIndex = myNet->numNeurons - 1; neuronIndex >= 0; neuronIndex--) {
		int weightIndex = myNet->neuron[neuronIndex].firstWeight;
		for (int source = myNet->neuron[neuronIndex].firstSource; source < myNet->neuron[neuronIndex].lastSource; source++)
			myNet->weight[weightIndex++] -= alpha * myNet->neuron[neuronIndex].delta*myNet->neuron[myNet->source[source]].value;
	}

	return(sumError);
}

double trainOnce(netStruct *myNet, double *inputs, double * desiredOutputs, double *actualOutputs, double alpha, double beta, int sync, double *weightChange)
// use log winner prob as target error
{
	int neuronIndex;
	double sumError = 0;
	feedForward(*myNet, inputs, actualOutputs, beta);
	//zero deltas
	for (int i = 0; i < myNet->numNeurons; i++)
		myNet->neuron[i].delta = 0;
	//set the delta at output layer and move it backwards



	for (int i = 0; i < myNet->numOutputs; i++) {

		neuronIndex = myNet->numNeurons - myNet->numOutputs + i;
		myNet->neuron[neuronIndex].error = (actualOutputs[i] - desiredOutputs[i]);
		if (myNet->neuron[neuronIndex].error < 0)myNet->neuron[neuronIndex].error *= COMPENSATION;
		sumError += SQR(myNet->neuron[neuronIndex].error);

		switch (myNet->neuron[neuronIndex].neuronFunction) {
		case SIGMOID:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].error*beta * 2 * (myNet->neuron[neuronIndex].value + 1)*(1 - myNet->neuron[neuronIndex].value);
			break;
		case RELU:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].error;
			if (myNet->neuron[neuronIndex].value < 0)myNet->neuron[neuronIndex].delta *= ALPHA;
			break;
		case LINEAR:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].error;
			break;
		}


		int weightIndex = myNet->neuron[neuronIndex].firstWeight;
		for (int sourceI = myNet->neuron[neuronIndex].firstSource; sourceI < myNet->neuron[neuronIndex].lastSource; sourceI++)
			myNet->neuron[myNet->source[sourceI]].delta += myNet->weight[weightIndex++] * myNet->neuron[neuronIndex].delta;

	}
	//set the delta at non outputlayer and move it backwards
	for (neuronIndex = myNet->numNeurons - 1 - myNet->numOutputs; neuronIndex >= 0; neuronIndex--) {
		switch (myNet->neuron[neuronIndex].neuronFunction) {
		case SIGMOID:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].delta*beta * 2 * (myNet->neuron[neuronIndex].value + 1)*(1 - myNet->neuron[neuronIndex].value);
			break;
		case RELU:
			if (myNet->neuron[neuronIndex].value < 0)myNet->neuron[neuronIndex].delta *= ALPHA;
			break;
		case LINEAR:
			myNet->neuron[neuronIndex].delta = myNet->neuron[neuronIndex].error;
			break;
		}
		int weightIndex = myNet->neuron[neuronIndex].firstWeight;
		for (int sourceI = myNet->neuron[neuronIndex].firstSource; sourceI < myNet->neuron[neuronIndex].lastSource; sourceI++) {



			myNet->neuron[myNet->source[sourceI]].delta += myNet->weight[weightIndex++] * myNet->neuron[neuronIndex].delta;
		}

	}

	//adjust weights
	for (neuronIndex = myNet->numNeurons - 1; neuronIndex >= 0; neuronIndex--) {
		int weightIndex = myNet->neuron[neuronIndex].firstWeight;
		for (int sourceI = myNet->neuron[neuronIndex].firstSource; sourceI < myNet->neuron[neuronIndex].lastSource; sourceI++) {
			if (sync)weightChange[weightIndex] -= alpha * myNet->neuron[neuronIndex].delta*myNet->neuron[myNet->source[sourceI]].value;
			else myNet->weight[weightIndex] -= alpha * myNet->neuron[neuronIndex].delta*myNet->neuron[myNet->source[sourceI]].value;
			weightIndex++;
		}
	}

	return(sumError);
}

double trainBatch(netStruct *myNet, double **inputs, double **outputs, int numSamples, double alpha, double beta, int sync)
{
	double sumError = 0;
	double * weightChange = NULL;
	double *tempOutputs = (double *)malloc(sizeof(double)*myNet->numOutputs);
	if (sync)weightChange = (double *)calloc(sizeof(double), myNet->numWeights);
	int *rperm = (int *)malloc(sizeof(int)*numSamples);
	//randperm(numSamples,rperm);  
	for (int k = 0; k < numSamples; k++)rperm[k] = rand() % numSamples;// deaw with repetitions

	for (int i = 0; i < numSamples; i++) {

		int k = rperm[i];
		//k=i;   //DEBUG
		sumError += trainOnce(myNet, inputs[k], outputs[k], tempOutputs, alpha, beta, sync, weightChange);
		//sumError+=trainOnce(myNet,inputs[k],outputs[k],tempOutputs,alpha,beta);
		//if(!(i%100))fprintf(stderr,"%d %f\r",i,sumError);
	}
	if (sync) {
		for (int weightIndex = 0; weightIndex < myNet->numWeights; weightIndex++)
			myNet->weight[weightIndex] += weightChange[weightIndex] / numSamples;
		free(weightChange);
	}

	free(tempOutputs);
	free(rperm);

	free(outputs);
	outputs = NULL;
	free(inputs);
	inputs = NULL;

	return(sumError / numSamples);
}


double trainRP(netStruct *myNet, double **inputs, double **outputs, int numSamples, double alpha, double beta, double *errorArray, double *sumErrorArray)
{
	// take each sample at a probability proportoional to its lst known error
	double sumError = 0;
	double thisError;
	int ii;

	double *tempOutputs = (double *)malloc(sizeof(double)*myNet->numOutputs);
	for (int i = 0; i < numSamples; i++) {
		unsigned int urand;
		urand = globalRNG::rand();
		double  random01 = (double)urand / globalRNG::max();
		double acc = 0;
		for (ii = 0; ii < numSamples - 1; ii++) {
			acc += errorArray[ii];
			if (acc > *sumErrorArray*random01)break;
		}

		sumError += thisError = -trainOnce(myNet, inputs[ii], outputs[ii], tempOutputs, alpha, beta, 0, NULL);
		*sumErrorArray += thisError - errorArray[ii];
		errorArray[ii] = thisError;
		if (!(i % 100))fprintf(stderr, "%d %f\r", i, *sumErrorArray);
	}

	free(tempOutputs);
	return(*sumErrorArray / numSamples);
}
void printNet(netStruct thisNet)
{
	printf("inputs: %d  outputs: %d\n", thisNet.numInputs, thisNet.numOutputs);

	for (int neuronIndex = 0; neuronIndex < thisNet.numNeurons; neuronIndex++) {
		printf("neuron %d value: %lf\n", neuronIndex, thisNet.neuron[neuronIndex].value);
		int weightIndex = thisNet.neuron[neuronIndex].firstWeight;
		int source = thisNet.neuron[neuronIndex].firstSource;
		for (int i = 0; i < thisNet.neuron[neuronIndex].lastSource - thisNet.neuron[neuronIndex].firstSource; i++)
			printf("source: %d weight: %lf\n", thisNet.source[source++], thisNet.weight[weightIndex++]);
	}
}
void saveReceptiveFields(netStruct myNet, char *fname)
// save if file fname the receptive field weights of all layer1 neurons.
{
	FILE *fi = fopen(fname, "wb");
	assert(fi != NULL);

	for (int k = 0; k < myNet.numNeurons; k++)
		if (myNet.neuron[k].layerIndex == 1)
			for (int kk = myNet.neuron[k].firstWeight; kk < myNet.neuron[k].lastWeight; kk++)
				fprintf(fi, "%f\n", myNet.weight[kk]);

	fclose(fi);
	fi = NULL;
}




void saveNetwork(netStruct thisNet, char *fname)
{
	FILE *fi = fopen(fname, "wb");
	assert(fi != NULL);

	fwrite(&thisNet, 1, sizeof(thisNet), fi);
	fwrite(thisNet.neuron, sizeof(thisNet.neuron[0]), thisNet.numNeurons, fi);
	fwrite(thisNet.weight, sizeof(thisNet.weight[0]), thisNet.numWeights, fi);
	fwrite(thisNet.source, sizeof(thisNet.source[0]), thisNet.numSource, fi);
	fclose(fi);
}
void serializeNetwork(void *ptr, netStruct thisNet, size_t *size)
// ptr is pre alloced for the networkGetSize;
{
	*size = 0;
	memcpy((char *)ptr + *size, &thisNet, sizeof(thisNet));
	*size += sizeof(thisNet);
	memcpy((char *)ptr + *size, thisNet.neuron, sizeof(thisNet.neuron[0])*thisNet.numNeurons);
	*size += sizeof(thisNet.neuron[0])*thisNet.numNeurons;
	memcpy((char *)ptr + *size, thisNet.weight, sizeof(thisNet.weight[0])*thisNet.numWeights);
	*size += sizeof(thisNet.weight[0])*thisNet.numWeights;
	memcpy((char *)ptr + *size, thisNet.source, sizeof(thisNet.source[0])*thisNet.numSource);
	*size += sizeof(thisNet.source[0])*thisNet.numSource;
}
size_t networkGetSize(netStruct thisNet)
{
	size_t size = 0;
	size += sizeof(thisNet);
	size += sizeof(thisNet.neuron[0])*thisNet.numNeurons;
	size += sizeof(thisNet.weight[0])*thisNet.numWeights;
	size += sizeof(thisNet.source[0])*thisNet.numSource;
	return(size);
}

netStruct  loadNetwork(char *fname)
{
	// malloc is used. Ude free when getting rid of network
	netStruct loadedNet;
	FILE *fi = fopen(fname, "rb");
	assert(fi != NULL);

	fread(&loadedNet, sizeof(loadedNet), 1, fi);

	loadedNet.neuron = (neuronStruct *)malloc(sizeof(loadedNet.neuron[0])*loadedNet.numNeurons);
	loadedNet.weight = (double *)malloc(sizeof(loadedNet.weight[0])*loadedNet.numWeights);
	loadedNet.source = (int *)malloc(sizeof(loadedNet.source[0])*loadedNet.numSource);
	fread(loadedNet.neuron, sizeof(loadedNet.neuron[0]), loadedNet.numNeurons, fi);
	fread(loadedNet.weight, sizeof(loadedNet.weight[0]), loadedNet.numWeights, fi);
	fread(loadedNet.source, sizeof(loadedNet.source[0]), loadedNet.numSource, fi);
	fclose(fi);
	return(loadedNet);
}
netStruct deserializeNetwork(unsigned char *blob)
{
	netStruct thisNet;
	size_t advance = 0;
	memcpy(&thisNet, blob, advance = sizeof(thisNet));
	blob += advance;
	thisNet.neuron = (neuronStruct *)malloc(sizeof(thisNet.neuron[0])*thisNet.numNeurons);
	thisNet.weight = (double *)malloc(sizeof(thisNet.weight[0])*thisNet.numWeights);
	thisNet.source = (int *)malloc(sizeof(thisNet.source[0])*thisNet.numSource);
	assert(thisNet.neuron&&thisNet.weight&&thisNet.source);
	memcpy(thisNet.neuron, blob, advance = sizeof(thisNet.neuron[0])*thisNet.numNeurons);
	blob += advance;
	memcpy(thisNet.weight, blob, advance = sizeof(thisNet.weight[0])*thisNet.numWeights);
	blob += advance;
	memcpy(thisNet.source, blob, advance = sizeof(thisNet.source[0])*thisNet.numSource);
	blob += advance;
	return(thisNet);

}
typedef struct {
	int index;
	unsigned int randnum;
}permRecord;
int permCompare(const void *a, const void *b)
{
	return(((permRecord *)a)->randnum - ((permRecord *)b)->randnum);

}

void randperm(int n, int perm[])
{

	permRecord *permArray;

	permArray = (permRecord *)malloc(n * sizeof(permArray[0]));
	for (int i = 0; i < n; i++) {
		permArray[i].index = i;
		permArray[i].randnum = globalRNG::rand();
	}
	qsort(permArray, n, sizeof(permArray[0]), permCompare);
	for (int i = 0; i < n; i++)
		perm[i] = permArray[i].index;

	free(permArray);
}
void destroyNet(netStruct &net)
{
	if (net.neuron)free(net.neuron);
	if (net.source)free(net.source);
	if (net.weight)free(net.weight);
	memset(&net, 0, sizeof(net));
}

MedBP::~MedBP() {
	// free all dynamic resources

	destroyNet(network);
}