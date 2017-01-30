//==============================================================
// micNet2
//==============================================================
// simple deep learning network implementation.
// This variant is the second phase - writing things more ordered and flexible + adding lots of new options
//

#ifndef __micNet__H__
#define __minNet__H__

#include <vector>
#include "MedUtils/MedUtils/MedUtils.h"

using namespace std;

enum NodesTypes {Node_Input=0, Node_LeakyReLU, NodeSoftMax, Node_Normalization, Node_Regression, Node_MaxOut, Node_Chooser};

// next class is used when reading a network structure from the user, 
// and contains the full information about the type of each node , its parameters, and the graph structure
// format to read a node-info: 
// <id=[int]|type=[LR,SoftMax,Norm,Regression,MaxOut,Chooser,Noise]|n_hidden=[int]|sources=[int:int:...]|lr=[float]|drop=[float]|...
// ...keep|A=[float]|B=[float]|lambda=[float]|momentum=[float]|decay=[float]|max_width=[int]|noise_std=[float]|...
// ...wgt_std=[float]|train=[int]>
class NodeInfo {
public:

	// inputs

	int id;						// id 0 is always special and is always the input layer (x)
	int keep_node;				// if 1 will keep current node as is, however one can still change lr,drop,momentum,decay
	string type;				// "Input","LeakyReLU","SoftMax","Normalization","Regression","MaxOut","Chooser"
	int n_hidden;				// number of hidden layers.
	vector<int> sources;		// input nodes, given in order.
	float learn_rate;			// learn_rate
	float momentum;				// momentum for this node
	float rate_decay;			// learn rate decay for this node.
	float drop_in;				// dropout on input
	float A, B;					// variables for LeakyReLU
	float lambda;				// lambda for L2 Regularization whenever relevant
	int max_width;				// width for max out layers
	float noise_std;			// noise level to add to input in training stages
	float wgt_std;				// std for initialization, if == 0 calculated auto as 2/sqrt(n_in)
	int to_train;				// if 1 (default) weights will be learned, if 0 current weights stay

	// calculated from inputs
	vector<int> sinks;						// output nodes for this node
	vector<pair<int, int>> sources_idx;		// indexes of sources into the input matrix
	vector<pair<int, int>> sinks_idx;		// matching indexes in sink matrix
	float drop_out;
	int in_dim, out_dim;					// non bias dimensions for input and output

											// init
	NodeInfo() {
		id=-1; keep_node=0; type=""; n_hidden=0; learn_rate = -1; momentum = -1; rate_decay = -1; drop_in = -1;
		drop_out = -1;  A = -1; B = -1; lambda = -1; max_width = 0; noise_std = -1; wgt_std = -1; to_train = 1;
	}

	void init_from_string(const string &in_str);

};

//
// general base class for a node
// has the basic structures that all nodes share
// each specific node type should expand the missing featues.
// each node must maintain the basic functions such as : init , forward, backward
// and maintaining all that is needed in order that is needed for next nodes to forward (batch_out) and previous nodes to backprop (...)
// 
class micNetNode {
public:
	int type;
	NodeInfo ni;

	int mode;					// either NET_MODE_TRAIN or NET_MODE_TEST

	int n_in, k_out;			// dimensions of input without bias: n x k
	MedMat<float> batch_in;		// matrix of nb x (n + 1) : batch_in(i,n) is always 1 to work with the bias term
	//MedMat<float> *batch_in;	// pointer to batch_in, as sometimes the batch_out of previous node is the batch_in, and this saves copying
	MedMat<float> batch_out;	// matrix of nb x (k + 1) : batch_out(i,k) is always 1 to work with the next bias term
	MedMat<float> wgt;			// size: (n + 1) x (k + 1) : wgt(n,j) is the bias for neuron j, wgt(i,k) is 0, wgt(n,k) is 1 to transfer the 1 bias carrier. 
	MedMat<float> grad_w;		// size: (n + 1) x (k + 1) : gradient for wgt. wgt(i,k) is always 0.
	MedMat<float> F;			// size: nb x (k + 1) :: dL/dZ :: recursion F = delta * wgt^T
	MedMat<float> G;			// size: nb x (k + 1) :: dZ/dS :: the gradient of the actual activation function
	MedMat<float> delta;		// size: nb x (k + 1) :: recursion : delta = F(next) @ G
	MedMat<float> S;			// size: nb x (k + 1) :: linear response before non linear activation :: S = batch_in * wgt

	int get_S_flag; // to be set by each node type , to sign Forward if to get S or not
	int Forward();	// general implementation that will use the forward of the specific node when needed
	int Backward();	// general implementation that will use the backward of the specific node when needed

	virtual int init(NodeInfo _ni) = 0;
	virtual int forward() { return 0; }
	virtual int backward() { return 0; }
};

//
// LeakyRelu Node:
// The LeakyReLU activation function is :
//
//			{ A    ...  if <w,x> >= 0
// f(x) =   {
//          { B    ...  if <w,x> < 0
//
// Typically A=1 B=0 or 0.01 or 0.1
//

class micNodeLeakyReLU : public micNetNode {
	micNodeLeakyReLU() { type = Node_LeakyReLU; get_S_flag = 1; }

	int forward();
};

#endif