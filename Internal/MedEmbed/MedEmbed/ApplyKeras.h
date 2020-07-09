#pragma once

//
// A simple class to enable:
// (1) reading a (simple) keras model from a text file (prepared in the training scripts in the print_layers function)
// (2) Applying it to a given sparse mat, or a line to get the embedding layer (currently always one before the last)
// (3) serialize/deserialize it via MedSerialize (to enable upper classes put it inside a model)
//
// Currently only supports : 
// (1) dense (fully connected, with bias) , with linear or relu or sigmoid activation
// (2) LeakyReLU activation (typically right AFTER a dense layer)
// (3) Dropout (nothing to do, may need to multiply by a factor)
// (4) batch normalization
// (5) Softmax or Sigmoid
//

#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedSparseMat/MedSparseMat/MedSparseMat.h>
#include <MedMat/MedMat/MedMat.h>
#include <External/Eigen/Core>

//===================================================================================================

enum KerasLayerTypes { K_UNKNOWN = 0, K_DENSE, K_LEAKY, K_DROPOUT, K_BN , K_ACTIVATION};
enum KerasActivations { A_UNKNOWN = 0, A_LINEAR, A_SIGMOID , A_RELU, A_LEAKY, A_SOFTMAX};

//===================================================================================================
class KerasLayer : public SerializableObject {

public:

	int type = K_UNKNOWN;
	string name = "";

	// for dense
	int in_dim = 0;
	int out_dim = 0;
	int n_bias = 0;
	int activation = A_UNKNOWN;

	// bn
	int dim;

	// dropout
	float drop_rate;

	// leaky
	float leaky_alpha;

	// weights (holding also transposed version for ease of use with Eigen)
	MedMat<float> wgts, twgts;

	// biases
	vector<float> bias;

	// helpers
	unordered_map<string, int> name_to_type = { { "dense" , K_DENSE } ,{ "leaky", K_LEAKY } ,{ "dropout" , K_DROPOUT },{ "batch_normalization" , K_BN } , {"activation", K_ACTIVATION} };
	unordered_map<string, int> name_to_activation = { { "linear" , A_LINEAR } ,{ "relu", A_RELU } ,{ "leaky" , A_LEAKY },{ "sigmoid" , A_SIGMOID }, {"softmax", A_SOFTMAX} };


	// appliers for a single sample
	int apply_sparse(vector<pair<int, float>> &sline, vector<float> &output) const;
	int apply_sparse(map<int, float> &sline, vector<float> &output) const;
	int apply(vector<float> &in, vector<float> &out) const;
	int apply_bn(vector<float> &in, vector<float> &out) const;
	int apply_activation(vector<float> &in, vector<float> &out) const; // for in place send same vector for in/out

	// appliers for a batch of samples
	int apply(MedMat<float> &in, MedMat<float>& out) const;
	int apply_bn(MedMat<float> &in, MedMat<float> &out) const;
	int apply_activation(MedMat<float> &in, MedMat<float> &out) const; // for in place send same vector for in/out

	// initialization from string
	int init(map<string, string>& _map);

	ADD_CLASS_NAME(KerasLayer)
	ADD_SERIALIZATION_FUNCS(type, name, in_dim, out_dim, n_bias, activation, dim, drop_rate, leaky_alpha, wgts, twgts, bias)
};


//===================================================================================================
class ApplyKeras : public SerializableObject {

public:
	vector<KerasLayer> layers;

	int init_from_text_file(string layers_file);

	int apply_sparse(vector<pair<int, float>> &sline, vector<float> &output, int to_layer) const;
	int apply_sparse(map<int, float> &sline, vector<float> &output, int to_layer) const;
	int apply(vector<float>& line, vector<float> &output, int to_layer) const;
	int apply(vector<float>& line, vector<float> &output) const { return apply(line, output, (int)(layers.size() - 1)); }
	int apply(MedMat<float> &line, MedMat<float>& output, int to_layer) const;
	int apply(MedMat<float>& line, MedMat<float> &output) const { return apply(line, output, (int)(layers.size() - 1)); }

	int get_all_embeddings(MedSparseMat &smat, int to_layer, MedMat<float> &emat);

	int get_output_dimension() {
		if (layers.size() > 0)
			return layers.back().out_dim;
		return 0;
	}

	ADD_CLASS_NAME(ApplyKeras)
	ADD_SERIALIZATION_FUNCS(layers)
};

MEDSERIALIZE_SUPPORT(KerasLayer)
MEDSERIALIZE_SUPPORT(ApplyKeras)
