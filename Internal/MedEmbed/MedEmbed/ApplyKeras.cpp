#include "ApplyKeras.h"
#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>                                                                                                                                                
#include <thread>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;


//===========================================================================================
// KerasLayer
//===========================================================================================

//-------------------------------------------------------------------------------------------
int KerasLayer::apply_sparse(vector<pair<int, float>> &sline, vector<float> &output) const
{
	if (type == K_DENSE) {

		// first calculate the linear transformations

		// init output with biases
		output.clear();
		output.resize(out_dim, 0);

		// linear parts

#pragma omp parallel for
		for (int d=0; d<out_dim; d++) {
			output[d] = bias[d];
			for (auto &p : sline)
				output[d] += wgts(d,p.first) * p.second;
		}

		if (activation != A_UNKNOWN && activation != A_LINEAR)
			return apply_activation(output, output);


	}
	else if (type == K_DROPOUT) {
		return 0; // seems that keras already multiplied it into the weights
		float factor = drop_rate;
		// in this special case : multiplying in place on the sparse vec
		for (auto &p : sline)
			p.second *= factor;
	}
	else {
		MTHROW_AND_ERR("apply sparse currently only available for Dense or dropout layers\n");
	}

	return 0;
}

//-------------------------------------------------------------------------------------------
int KerasLayer::apply_sparse(map<int, float> &sline, vector<float> &output) const
{
	if (type == K_DENSE) {

		// first calculate the linear transformations

		// init output with biases
		output.clear();
		output.resize(out_dim, 0);

		// linear parts

#pragma omp parallel for
		for (int d = 0; d<out_dim; d++) {
			output[d] = bias[d];
			for (auto &p : sline)
				output[d] += wgts(d,p.first) * p.second;
		}

		if (activation != A_UNKNOWN && activation != A_LINEAR)
			return apply_activation(output, output);


	}
	else if (type == K_DROPOUT) {
		return 0; // seems that keras already multiplied it into the weights
		float factor = drop_rate;
		// in this special case : multiplying in place on the sparse vec
		for (auto &p : sline)
			p.second *= factor;
	}
	else {
		MTHROW_AND_ERR("apply sparse currently only available for Dense or dropout layers\n");
	}

	return 0;
}

//-------------------------------------------------------------------------------------------
int KerasLayer::apply_activation(vector<float> &in, vector<float> &out) const
{
	out.resize(in.size());
	
	// if needed apply sigmoid
	if (activation == A_SIGMOID) {
		for (int d=0; d<in.size(); d++)
			out[d] = 1.0f / (1.0f + exp(-in[d]));
	}

	// if needed apply ReLU
	else if (activation == A_RELU) {
		for (int d=0; d<in.size(); d++)
			if (in[d] < 0)
				out[d] = 0;
			else
				out[d] = in[d];
	}

	else if (activation == A_LEAKY) {
		for (int d=0; d<in.size(); d++) {
			out[d] = in[d];
			if (in[d] < 0)
				out[d] *= leaky_alpha;
		}
	}

	else if (activation == A_SOFTMAX) {
		float sum = 0;
		for (int d = 0; d < in.size(); d++) {
			out[d] = exp(in[d]);
			sum += out[d];
		}
		for (int d = 0; d < in.size(); d++)
			out[d] = out[d] / sum;

	}
	return 0;
}

//-------------------------------------------------------------------------------------------
int KerasLayer::apply_bn(vector<float> &in, vector<float> &out) const
{
	out = in;
#pragma omp parallel for
	for (int d=0; d<dim; d++) {
		out[d] = (out[d] - wgts(2,d))/wgts(3,d);
		out[d] = wgts(0,d) * out[d] + wgts(1,d);
		//MLOG("BN: d=%d : in %f wgts %f %f %f %f out %f\n", d, in[d], wgts(0,d), wgts(1,d), wgts(2,d), wgts(3,d), out[d]);
	}
	return 0;
}

//-------------------------------------------------------------------------------------------
int KerasLayer::apply(vector<float> &in, vector<float> &out) const
{
	if (type == K_DENSE) {
		// first calculate the linear transformations

		// init output with biases
		out = bias;

		// linear parts
#pragma omp parallel for
		for (int d=0; d<out_dim; d++) {
			for (int j=0; j<in_dim; j++)
				out[d] += wgts(d,j) * in[j];
		}

		if (activation != A_UNKNOWN && activation != A_LINEAR)
			return apply_activation(out, out);
	}
	else if (type == K_LEAKY) {
		return apply_activation(in, out);
	}
	else if (type == K_BN) {
		apply_bn(in, out);
	}
	else if (type == K_DROPOUT) {
		out = in; return 0;
		float factor = drop_rate;
		for (int d=0; d<out.size(); d++)
			out[d] *= factor;
	}
	else if (type == K_ACTIVATION) {
		return apply_activation(in, out);
	}

	return 0;
}

//-------------------------------------------------------------------------------------------
// Work on batches
//-------------------------------------------------------------------------------------------
int KerasLayer::apply(MedMat<float> &in, MedMat<float> &out) const
{

	if (type == K_DENSE) {
		out.resize(in.nrows, twgts.ncols);

		Eigen::Map<const Eigen::MatrixXf> eigen_in(in.data_ptr(), in.ncols, in.nrows);
		Eigen::Map<const Eigen::MatrixXf> eigen_w(twgts.data_ptr(), twgts.ncols, twgts.nrows);
		Eigen::Map<const Eigen::VectorXf> eigen_b(bias.data(), bias.size());
		Eigen::Map<Eigen::MatrixXf> eigen_out(out.data_ptr(), out.ncols, out.nrows);

		eigen_out = eigen_w * eigen_in ;
		eigen_out.colwise() += eigen_b;

		if (activation != A_UNKNOWN && activation != A_LINEAR)
			return apply_activation(out, out);
	}
	else if (type == K_LEAKY) {
		return apply_activation(in, out);
	}
	else if (type == K_BN) {
		apply_bn(in, out);
	}
	else if (type == K_DROPOUT) {
		out = in;
	}
	else if (type == K_ACTIVATION) {
		return apply_activation(in, out);
	}

	return 0;
}

//-------------------------------------------------------------------------------------------
int KerasLayer::apply_activation(MedMat<float> &in, MedMat<float> &out) const
{
	out.resize(in.nrows,in.ncols);

	// if needed apply sigmoid
	if (activation == A_SIGMOID) {
#pragma omp parallel for if (in.nrows > 4)
		for (int irow=0; irow<in.nrows; irow++)
			for (int icol=0; icol<in.ncols; icol++)
				out(irow,icol)= 1.0f / (1.0f + exp(-in(irow,icol)));
	}

	// if needed apply softmax
	if (activation == A_SOFTMAX) {
#pragma omp parallel for if (in.nrows > 4)
		for (int irow = 0; irow < in.nrows; irow++) {
			float sum = 0;
			for (int icol = 0; icol < in.ncols; icol++) {
				out(irow, icol) = exp(in(irow, icol));
				sum += out(irow, icol);
			}
			for (int icol = 0; icol < in.ncols; icol++)
				out(irow, icol) = out(irow, icol) / sum;
		}
	}

	// if needed apply ReLU
	else if (activation == A_RELU) {
#pragma omp parallel for if (in.nrows > 4)
		for (int irow = 0; irow<in.nrows; irow++)
			for (int icol = 0; icol < in.ncols; icol++) {
				if (in(irow, icol) < 0)
					out(irow, icol) = 0;
				else
					out(irow, icol) = in(irow, icol);
			}
	}

	else if (activation == A_LEAKY) {
#pragma omp parallel for if (in.nrows > 4)
		for (int irow = 0; irow<in.nrows; irow++)
			for (int icol = 0; icol < in.ncols; icol++) {
				out(irow, icol) = in(irow, icol);
				if (in(irow, icol) < 0)
					out(irow, icol) *= leaky_alpha;
		}
	}

	return 0;
}

//-------------------------------------------------------------------------------------------
int KerasLayer::apply_bn(MedMat<float> &in, MedMat<float> &out) const
{
	out = in;
#pragma omp parallel for
	for (int irow = 0; irow<in.nrows; irow++) {
		for (int d = 0; d < in.ncols; d++) {
			out(irow,d) = (out(irow, d) - wgts(2,d)) / wgts(3,d);
			out(irow, d) = wgts(0,d) * out(irow, d) + wgts(1,d);
		}
	}
	return 0;
}

//-------------------------------------------------------------------------------------------
int KerasLayer::init(map<string, string>& _map)
{
	for (auto entry : _map) {
		//MLOG("es init(): %s -> %s\n", entry.first.c_str(), entry.second.c_str());
		string field = entry.first;
		if (field == "name") name = entry.second;
		else if (field == "in_dim") in_dim = stoi(entry.second);
		else if (field == "out_dim") out_dim = stoi(entry.second);
		else if (field == "n_bias") n_bias = stoi(entry.second);
		else if (field == "dim") dim = stoi(entry.second);
		else if (field == "drop_rate") drop_rate = stof(entry.second);
		else if (field == "leaky_alpha") leaky_alpha = stof(entry.second);
		else if (field == "type") {
			if (name_to_type.find(entry.second) == name_to_type.end())
				MTHROW_AND_ERR("KerasLayer: ERROR: unsupported type : %s\n", entry.second.c_str());
			type = name_to_type[entry.second];
		}
		else if (field == "activation") {
			if (name_to_activation.find(entry.second) == name_to_type.end())
				MTHROW_AND_ERR("KerasLayer: ERROR: unsupported activation : %s\n", entry.second.c_str());
			activation = name_to_activation[entry.second];
		}
	}

	return 0;
}


//===========================================================================================
// ApplyKeras
//===========================================================================================
//-------------------------------------------------------------------------------------------
int ApplyKeras::apply_sparse(vector<pair<int, float>> &sline, vector<float> &output, int to_layer) const
{
	if (to_layer < 0) to_layer = (int)layers.size() - 1 + to_layer;

	int passed_first_dense = 0;
	vector<vector<float>> outs(layers.size(), vector<float>());
	for (int i=0; i<=to_layer; i++) {

		//MLOG("apply_sparse layer %d (passed %d) to_layer=%d (name %s type %d activation %d)\n", i, passed_first_dense, to_layer, layers[i].name.c_str(), layers[i].type, layers[i].activation);
		if (passed_first_dense == 0) {
			layers[i].apply_sparse(sline, outs[i]);
			if (layers[i].type == K_DENSE)
				passed_first_dense = 1;
		}
		else {
			layers[i].apply(outs[i-1], outs[i]);
		}
		//MLOG("out[i] size %d\n", outs[i].size());
	}

	output = outs[to_layer];
	return 0;
}

//-------------------------------------------------------------------------------------------
int ApplyKeras::apply_sparse(map<int, float> &sline, vector<float> &output, int to_layer) const
{
	if (to_layer < 0) to_layer = (int)layers.size() - 1 + to_layer;

	int passed_first_dense = 0;
	vector<vector<float>> outs(layers.size(), vector<float>());
	for (int i = 0; i <= to_layer; i++) {

		//MLOG("apply_sparse layer %d (passed %d) to_layer=%d (name %s type %d activation %d)\n", i, passed_first_dense, to_layer, layers[i].name.c_str(), layers[i].type, layers[i].activation);
		if (passed_first_dense == 0) {
			layers[i].apply_sparse(sline, outs[i]);
			if (layers[i].type == K_DENSE)
				passed_first_dense = 1;
		}
		else {
			layers[i].apply(outs[i - 1], outs[i]);
		}
		//MLOG("out[i] size %d\n", outs[i].size());
	}

	output = outs[to_layer];
	return 0;
}

//-------------------------------------------------------------------------------------------
int ApplyKeras::apply(vector<float>& line, vector<float> &output, int to_layer) const
{
	if (to_layer < 0) to_layer = (int)layers.size() - 1 + to_layer;

	vector<vector<float>> outs(to_layer+1, vector<float>());
	layers[0].apply(line, outs[0]);
	
	for (int i = 1; i <= to_layer; i++) {
		//MLOG("apply layer %d (passed %d) to_layer=%d (name %s type %d activation %d)\n", i, passed_first_dense, to_layer, layers[i].name.c_str(), layers[i].type, layers[i].activation);
		layers[i].apply(outs[i - 1], outs[i]);
		//MLOG("out[i] size %d\n", outs[i].size());
	}

	output = outs[to_layer];
	return 0;

}

//-------------------------------------------------------------------------------------------
int ApplyKeras::apply(MedMat<float>& line, MedMat<float> &output, int to_layer) const
{
	if (to_layer < 0) to_layer = (int)layers.size() - 1 + to_layer;

	vector<MedMat<float>> outs(to_layer + 1, MedMat<float>());
	layers[0].apply(line, outs[0]);

	for (int i = 1; i <= to_layer; i++) {
		//MLOG("apply layer %d (passed %d) to_layer=%d (name %s type %d activation %d)\n", i, passed_first_dense, to_layer, layers[i].name.c_str(), layers[i].type, layers[i].activation);
		layers[i].apply(outs[i - 1], outs[i]);
		//MLOG("out[i] size %d\n", outs[i].size());
	}

	output = outs[to_layer];
	return 0;

}

//-------------------------------------------------------------------------------------------
int ApplyKeras::init_from_text_file(string layers_file)
{
	MLOG("Reading layers file %s\n", layers_file.c_str());
	ifstream inf(layers_file);

	if (!inf) {
		MERR("ApplyKeras: read: Can't open layers file %s\n", layers_file.c_str());
		return -1;
	}

	string curr_line;
	while (getline(inf, curr_line)) {
		if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			if (curr_line[curr_line.size() - 1] == '\r')
				curr_line.erase(curr_line.size() - 1);

			if (curr_line[0] == 'L') {

				MLOG("ApplyKeras: Reading %d : %s\n", (int)layers.size(), curr_line.c_str());
				// parsing a layers line
				vector<string> f;
				boost::split(f, curr_line, boost::is_any_of("\t"));

				if (f[0] != "LAYER") {
					MTHROW_AND_ERR("ApplyKeras: ERROR: explected a LAYER line, got %s\n", curr_line.c_str());
				}

				KerasLayer kl;
				kl.init_from_string(f[1]);

				// read matrices if needed
				if (kl.type == K_DENSE) {
					kl.wgts.resize(kl.out_dim, kl.in_dim);
					kl.twgts.resize(kl.in_dim, kl.out_dim);
					for (int j=0; j<kl.in_dim; j++) {
						getline(inf, curr_line);
						boost::split(f, curr_line, boost::is_any_of(","));
						if (f.size() != kl.out_dim)
							MTHROW_AND_ERR("ApplyKeras: ERROR: non matching sizes for wgts line: %d expected , got %d\n", kl.out_dim, (int)f.size());
						for (int d = 0; d < f.size(); d++) {
							kl.wgts(d, j) = stof(f[d]); // keras prints the matrices transposed.
							kl.twgts(j, d) = stof(f[d]);
						}
					}

					// read bias
					kl.bias.resize(kl.out_dim);
					getline(inf, curr_line);
					boost::split(f, curr_line, boost::is_any_of(","));
					if (f.size() != kl.bias.size())
						MTHROW_AND_ERR("ApplyKeras: ERROR: non matching sizes for bias wgts line: %d expected , got %d\n", (int)kl.bias.size(), (int)f.size());
					for (int j = 0; j < f.size(); j++)
						kl.bias[j] = stof(f[j]);
				}

				if (kl.type == K_BN) {
					kl.wgts.resize(4,kl.dim);
					for (int d=0; d<4; d++) {
						getline(inf, curr_line);
						boost::split(f, curr_line, boost::is_any_of(","));
						if (f.size() != kl.dim)
							MTHROW_AND_ERR("ApplyKeras: ERROR: non matching sizes for batch normalization wgts line: %d expected , got %d\n", kl.dim, (int)f.size());
						for (int j=0; j<f.size(); j++)
							kl.wgts(d,j) = stof(f[j]);
					}
					float epsilon = (float)0.001;
					for (int j=0; j<kl.dim; j++)
						kl.wgts(3, j) = sqrt(kl.wgts(3,j) + epsilon);
					kl.out_dim = kl.dim;
				}

				layers.push_back(kl);
			}
		}
	}

	return 0;
}

// running over all lines in smat and generating matching embedding lines in emat
int ApplyKeras::get_all_embeddings(MedSparseMat &smat, int to_layer, MedMat<float> &emat)
{
	if (to_layer < 0) to_layer = (int)layers.size() - 1 + to_layer;
	int nlines = (int)smat.lines.size();
	int edim = layers[to_layer].out_dim;

	emat.resize(nlines, edim);

#pragma omp parallel
	for (int i = 0; i < nlines; i++) {
		vector<float> i_out;
		apply_sparse(smat.lines[i], i_out, to_layer);
		for (int j = 0; j < i_out.size(); j++)
			emat(i, j) = i_out[j];
	}

	return 0;
}