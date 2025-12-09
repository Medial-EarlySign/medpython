//==============================================================
// micNet
//==============================================================
// simple deep learning network implementation.
//

#ifndef __micNet__H__
#define __micNet__H__

#include <vector>
#include "MedUtils/MedUtils/MedUtils.h"
#include <MedMat/MedMat/MedMat.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>

using namespace std;

//
// A micNet has several nodes inside it (micNodes)
// some of them (typically one) is the input node, 
// the others have rules to create inputs from previous nodes 
// and then apply the micNode function to get its output.
// There are also output nodes, in which in training phase we also get the y to learn from
// (typicaly we will use softmax for classification and least squares for regression).
// 
// Training:
// (1) We choose an initial state for the network (for example normal random weights and 0 biases)
// (2) We forward the batch through the network, so that each micNode has now batch_in and bactch_out.
// (3) We start from the output nodes back, calculate gradients for the each node weights
//     and propagate the differentials to get a gradient for the weights in each node.
// (4) We step a gradient decent step using the gradients we have and the learning rates.
// (5) We repeat with the next batch until conversion/stop criteria.
// (6) Once an Epoch - we do a full run of the network on all data forward (could again be done in batches to save memory) and get predictions
//     On-Train and On-Test. This monitors the conversion on the network.
//
// In theory we could implement any differntiable function f: R^n -> R^k as a node (!) this gives us much freedom to develop new ideas for functions.
//
// Testing:
// (1) Working in input batches again
// (2) Forwarding the inputs through the network - can be parralelized easily.
//


class micNode;
class micNet;

#define NET_MODE_TRAIN	1
#define NET_MODE_TEST	2

class InputRules {
public:
	int n_input_nodes;
	vector<int> in_node_id;
	vector<string> mode;

	void clear() { n_input_nodes = 0; in_node_id.clear(); mode.clear(); }
	InputRules() { clear(); }

	void push(int node_id, const string &_mode);
};


//
// A micNode implements a function from R^n -> R^k that has M*k optimization weights.
// n is typically the dimension of the previous layer, k will be the dimension for the next
// in other cases we might want to build and spread the input/output dimensions according to other more complex rules (example: convolutional networks)
// Inputs to a micNode can be coordinates from other nodes outputs, or from the input features
//
class micNode : public SerializableObject {

public:

	int id;
	string type;	// "Input","LeakyReLU","SoftMax","Normalization","Regression","MaxOut","Chooser"
	string subtype; // "encoder","decoder"
	string name;	// for prints 
	string loss;	// if "" this is not an output node, if "log" we use logloss (works for softmax), if "least-squares" uses that on input.

	int n_in;		// coming in dimension
	int k_out;		// going out dimension

	MedMat<float> wgt;			// weights. last weight is ALWAYS the bias (in LeakyReLU and Ridge). size (n_in + 1) x (kout + 1), last column is always 0,0,....1 (to transfer 1 forward)
	vector<float> dropout_in;	// if dropout is used, in each batch we randomize a dropout vector, stating which of the inputs are taken, and which of the outputs.
	vector<float> dropout_out;	// this one is for the outputs and must be similar to the next layer in.
	MedMat<int> sparse_bit; // ??

	InputRules ir;
	vector<int> forward_nodes;	// for convience, a list of the node indexes that are forward to this node.
	int	is_terminal;			// easy way to know if this one is a terminal node (output node --> loss node)

	// training related
	MedMat<float> lr_params;	// in LeakyReLU: holds the slopes a,b for each one of the k neurons. size k_out x 2
	MedMat<float> lambda;		// in Ridge and LeakyReLU, this one holds the ridge regularization coefficient for each neuron. size k_out x 1 .
	MedMat<float> learn_rates;	// learning rates for wgts, currently size is k_out x 1 , that is a learn rate for each neuron.

	float momentum;					// momentum to use while learning
	float rate_factor;				// multiplying the learning rate - allowing to manage its decay
	float max_wgt_norm;				// maximal norm for the weights
	float min_wgt_norm;				// minimal norm for the weights
	float dropout_prob_in;			// probability for dropout of the coming in variables
	float dropout_prob_out;			// probability for dropout of the coming out variables 
	float sparse_zero_prob;			// probability for initial constant random sparsness imposed on the weights

	MedMat<float> batch_in;		// current batch inputs. Every batch has a column of "1" as the last: this is the bias column, so the actual size should be batch_size x (n_in + 1)
	MedMat<float> batch_out;	// current batch outputs. size: batch_size x (k_out + 1) (last column always has 1's - to allow bias)

	MedMat<float> grad_w;		// current gradient for w. grad(i,j) = d_f/d_Wij(x) averaged over all x in the batch, size (n_in + 1) x (k_out + 1) , last col always 0
	MedMat<float> grad_s;		// current gradient for x. grad(i,j) = d_f/d_Sij(x) for each sample and neuron ,size (nb) x (k_out + 1) , last col always 0
	MedMat<float> delta;		// current gradient for x. delta(i,j) = d_L/d_Sij(x) for each sample and neuron ,size (nb) x (k_out + 1) , last col always 0

	MedMat<float> prev_grad_w;	// needed for being able to use momentum

	// normalization related
	vector<float> b_mean, b_var;		// keeping last stage mean and variance ... to be updated
	vector<float> alpha, beta;			// actual learnt alpha,beta needed for forward normalization
	float normalization_update_factor;	// in [0,1]
	vector<double> curr_mean, curr_var;		// used for current batch estimators

	// softmax related
	int n_categ;
	int n_per_categ;
	MedMat<float> full_probs; // just a helper array used in calculations , here since needed both in forward and backprop

	// autoencoder related
	int data_node; // node id to decode to
	micNode *data_node_p;


	micNet* my_net; // pointer to the container network

	MedMat<float> y;			// the y values for a batch. This is relevant for output nodes only, and typically will be for a SoftMax
	vector<float> sweights;		// samples weights for this node

	// next is typically 1, however, we may set a micNode (=a layer) not to update weights during a learn cycle
	// This is useful when one needs to keep a layer in the middle with no changes, while letting other layers work.
	// Delta propagations will still flow through the node, to allow lower levels to get correct gradients.
	int update_weights_flag;

	// this one forces ONLY forward passes on the node in learn time. This is relevant when we freeze the first layers, and hence have
	// no need to propagate through them.
	int only_forward_flag;


	micNode() { data_node = -1; data_node_p = NULL; update_weights_flag = 1; only_forward_flag = 0; subtype = ""; }

	// initialize weights random in a uniform segment
	int init_wgts_rand(float min_range, float max_range);

	// initialize weights from a normal distribution
	int init_wgts_rand_normal(float mean, float std);

	int fill_input_node(int *perm, int len, MedMat<float> &x_mat, int last_is_bias_flag);	// copies x into input nodes 
	int fill_output_node(int *perm, int len, MedMat<float> &y_mat, vector<float> &sample_weights);							// copies y into input nodes 

	int get_input_batch(int do_grad_flag);  // sets batch_in for a node that is not an input node
	int forward_batch(int do_grad_flag);

	int forward_batch_leaky_relu(int do_grad_flag);
	int forward_batch_normalization(int do_grad_flag);
	int forward_batch_softmax(int do_grad_flag);
	int forward_batch_regression(int do_grad_flag);

	void forward_batch(const vector<MedMat<float>> &nodes_outputs, MedMat<float> &out) const;
	void get_input_batch(const vector<MedMat<float>> &nodes_out, MedMat<float> &in) const;
	void forward_batch_leaky_relu(const MedMat<float> &in, MedMat<float> &out) const;
	void forward_batch_normalization(const MedMat<float> &in, MedMat<float> &out) const;
	void forward_batch_softmax(const MedMat<float> &in, MedMat<float> &out) const;
	void forward_batch_regression(const MedMat<float> &in, MedMat<float> &out) const;


	int back_propagete_from(micNode *next);
	int get_backprop_delta();

	int weights_gd_step();
	int weights_normalization_step(); // step for a normalization layer


	void print(const string &prefix, int i_state, int i_in);

	std::default_random_engine gen;


	// serializations for a single node (partial... only what's needed by predictions, and not initialized by init_params)

	ADD_CLASS_NAME(micNode)
		ADD_SERIALIZATION_FUNCS(id, wgt, alpha, beta)

};


class NetEval {

public:
	string name;
	int epoch;
	float acc_err;
	float auc_max;
	float auc_exp;
	float corr_max;
	float corr_exp;
	float log_loss;
	float lsq_loss;

	double dt;	// time for epoch (for time performance measurements)

};


//
// micNet - an implementation of a multilayerd NN on top of the micNode class
// 
// currently - fully conected LeakyReLU layers with SoftMax/logloss at the end
//
class micNetParams : public SerializableObject {
public:
	string params_init_string;

	int batch_size;
	int predict_batch_size;
	int n_categ;
	int n_per_categ;
	int nfeat;

	vector<float> samp_ratio; // if size is less than n_categ, just permuting whole data each epoch
							  // if n_categ weights given: chooses by the probability given randomly for each batch 
							  //(with repetitions...but prob for that should be close to 0 on reasonable data)

	float max_wgt_norm;
	float min_wgt_norm;
	float weights_init_std;
	float rate_decay;

	float def_A, def_B; // leaky ReLU defaults
	float def_learning_rate, def_lambda, def_momentum; // params for case of constant params to all layers in net
	vector<float> learning_rates; // learn rate for each layer.

	int n_norm_layers;
	float normalization_factor;
	float sparse_zero_prob;

	vector<int> n_hidden;
	vector<float> dropout_in_probs;

	string net_type; // "fc","autoencoder"
	string loss_type;

	// a few more params needed to comply with MedAlgo
	int min_epochs;
	int max_epochs;
	int n_back;
	float min_improve_n_back; // minimal relative improvemnt on train set looking n_back steps back
	int n_preds_per_sample; // can be 1 or n_categ
	int pred_class;			// the class we will put in preds in case n_preds_per_sample == 1


	// next params are needed for retraining and transfer learning
	// the idea is that we start with a ready network, cut off some layers at the end,
	// and add new layers on top of it.
	// a classic use is train a net with an autoencoder, and then add classification layers on top of it and train them.
	int last_layer_to_keep; // will keep layers 0 (input) up to layers <= last_layer_to_keep (for example should be 1 for the simplest one layer autoencoder)
							// if < 0 then restarting network fresh


//	vector<NodeInfo> node_infos;

	void init_defaults() {
		batch_size = 1024; predict_batch_size = 30000; nfeat = 0; n_categ = 0; n_per_categ = 1; max_wgt_norm = 0; min_wgt_norm = 0;
		weights_init_std = (float)0.01; rate_decay = (float)0.97;
		n_hidden.clear();
		dropout_in_probs.clear();
		loss_type = "log";
		def_A = (float)1.0; def_B = (float)0.01;
		def_learning_rate = (float)0.1; def_lambda = 0; def_momentum = (float)0.9; normalization_factor = (float)0.99;
		n_norm_layers = 0;
		sparse_zero_prob = 0;
		net_type = "fc";
		pred_class = 1;
		n_preds_per_sample = 1;
		min_epochs = 10;
		max_epochs = 100;
		n_back = 10;
		min_improve_n_back = (float)0.001;
		samp_ratio.clear();
		last_layer_to_keep = -1;
		//		node_infos.clear();
	}

	micNetParams() { init_defaults(); }
	int init_from_string(const string &init_str);

	int node_infos_init_finish();
};


class micNet : public SerializableObject {

public:

	int version = 0;
	vector<micNode> nodes;
	micNetParams params;

	vector<micNode> nodes_last_best;

	micNet() { nodes.clear(); }

	void copy_nodes(vector<micNode> &in_nodes);		// needed in order to set up ir pointers correctly

	// adding layers (relying on params to be already initialized)
	int add_input_layer();
	int add_fc_leaky_relu_layer(int in_node, int n_hidden, float dropout_out_p, float sparse_prob, float learn_rate);
	int add_normalization_layer(int in_node);
	int add_softmax_output_layer(int in_node);
	int add_regression_output_layer(int in_node);
	int add_autoencoder_loss(int in_node, int data_node); // in_node is encoder node in this case

	// initialization
	int init_fully_connected(micNetParams &in_params);
	int init_fully_connected(const string &init_str);
	int init_autoencoder(micNetParams &in_params);
	int init_net(const string &init_string); // auto choose initialization (using the net_type param)
	int init_net(micNetParams &in_params);


	// forward and backprop
	int forward_batch(int do_grad_flag);	// assumes Input nodes contain the batch in batch_out
	int back_prop_batch();	// assumes forward batch was run before

	// learn, eval and predict
	int learn(MedMat<float> &x_train, MedMat<float> &y_train, vector<float> &weights, MedMat<float> &x_test, MedMat<float> &y_test, int n_epochs, int eval_freq, int last_is_bias_flag = 0);
	int learn_single_epoch(MedMat<float> &x_train, MedMat<float> &y_train, vector<float> &weights, int last_is_bias_flag = 0);
	int eval(const string &name, MedMat<float> &x, MedMat<float> &y, NetEval &eval, int last_is_bias_flag = 0);
	int predict(MedMat<float> &x, MedMat<float> &preds, int last_is_bias_flag = 0);
	void predict_single(const vector<float> &x, vector<float> &preds) const;

	vector<vector<int>> index_by_categ; // used when choosing a random batch by samp_ratio
	int get_batch_with_samp_ratio(MedMat<float> &y_train, int batch_len, vector<int> &chosen);

	// next is used to estimate the gradient at a specific node numerically and compare to the existing gradient there
	// This is done for debugging
	int test_grad_numerical(int i_node, int i_in, int i_out, float epsilon);


	// API to allow use with MedAlgo
	int init_from_string(string init_str);
	int learn(MedMat<float> &x_train, MedMat<float> &y_train) { vector<float> w; return learn(x_train, y_train, w); }
	int learn(MedMat<float> &x_train, MedMat<float> &y_train, vector<float> &weights);
	int predict(MedMat<float> &x, vector<float> &preds);

	size_t get_size() { return MedSerialize::get_size(version, params.params_init_string, nodes); }
	size_t serialize(unsigned char *blob) { return MedSerialize::serialize(blob, version, params.params_init_string, nodes); }
	size_t deserialize(unsigned char *blob) {
		string init_str;
		size_t size = MedSerialize::deserialize(blob, version, init_str);
		fprintf(stderr, "micNet deserialize init with %s\n", init_str.c_str());
		init_net(init_str);
		size += MedSerialize::deserialize(&blob[size], nodes);
		for (auto &node : nodes) { node.my_net = this; }
		return size;
	}

	int n_preds_per_sample() const { return params.n_preds_per_sample; }

};

//=======================================================
// Joining the MedSerialize Wagon
//=======================================================
MEDSERIALIZE_SUPPORT(micNode);
MEDSERIALIZE_SUPPORT(micNetParams);
MEDSERIALIZE_SUPPORT(micNet);



#endif

