#define _CRT_SECURE_NO_WARNINGS

#include "MedAlgo.h"
#include "QRF/QRF/QRF.h"


#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <unordered_set>

// Quantized Random Forests

//..............................................................................
void MedQRF::init_defaults()
{
	classifier_type = MODEL_QRF;
	transpose_for_learn = false;
	normalize_for_learn = false;
	normalize_y_for_learn = false;
	transpose_for_predict = false;
	normalize_for_predict = false;

	params.ntrees = MED_QRF_DEF_NTREES;
	params.maxq = MED_QRF_DEF_MAXQ;
	params.min_node = MED_QRF_DEF_MIN_NODE;
	params.type = QRF_CATEGORICAL_ENTROPY_TREE;
	params.n_categ = 2;

	params.learn_nthreads = MED_QRF_DEF_LEARN_NTHREADS;
	params.predict_nthreads = MED_QRF_DEF_PREDICT_NTHREADS;
	params.sampsize = NULL;
	params.ntry = 0;
	params.spread = (float)MED_QRF_DEF_SPREAD;
	params.get_count = PROBS_CATEG_AVG_PROBS;
	params.max_depth = 0;

	params.get_only_this_categ = 1;
	params.max_samp = 0;
	params.samp_factor = 0;

	params.keep_all_values = false;
	params.quantiles.clear();

	params.collect_oob = 0;
}

//..............................................................................
int MedQRF::init(void *_in_params)
{
	init_defaults();

	MedQRFParams *in_params = (MedQRFParams *)_in_params;

	params.ntrees = in_params->ntrees;
	params.maxq = in_params->maxq;
	params.min_node = in_params->min_node;

	params.type = in_params->type;
	params.learn_nthreads = in_params->learn_nthreads;
	params.predict_nthreads = in_params->predict_nthreads;
	params.sampsize = in_params->sampsize;
	params.ntry = in_params->ntry;
	params.spread = in_params->spread;
	params.get_count = in_params->get_count;
	params.n_categ = in_params->n_categ;
	params.max_depth = in_params->max_depth;

	params.get_only_this_categ = in_params->get_only_this_categ;
	params.max_samp = in_params->max_samp;
	params.samp_factor = in_params->samp_factor;

	params.collect_oob = in_params->collect_oob;

	params.keep_all_values = in_params->keep_all_values;
	params.quantiles = in_params->quantiles;

	return 0;
}

//..............................................................................
//int MedQRF::init(const string &init_str)
//{
//	init_defaults();
//
vector<string> fields;
//	split(fields, init_str, boost::is_any_of(",="));
//	
//	for (int i = 0; i < fields.size(); i++) {

//		if (fields[i] == "type") { params.type = (QRF_TreeType)stoi(fields[++i]); }
//		if (fields[i] == "ntrees")	{ params.ntrees = stoi(fields[++i]); }
//		if (fields[i] == "maxq") { params.maxq = stoi(fields[++i]); }
//		if (fields[i] == "min_node") { params.min_node = stoi(fields[++i]); }
//		if (fields[i] == "ntry") { params.ntry = stoi(fields[++i]); }
//		if (fields[i] == "get_count") { params.get_count = stoi(fields[++i]); }
//		if (fields[i] == "get_only_this_categ") { params.get_only_this_categ = stoi(fields[++i]); }
//		if (fields[i] == "max_samp") { params.max_samp = stoi(fields[++i]); }
//		if (fields[i] == "samp_factor") { params.samp_factor = stof(fields[++i]); }
//		if (fields[i] == "n_categ") { params.n_categ = stoi(fields[++i]); }
//		if (fields[i] == "collect_oob") { params.collect_oob = stoi(fields[++i]); }
//		if (fields[i] == "spread") { params.spread = stof(fields[++i]); }
//		if (fields[i] == "sampsize") {
//			vector<string> vals;
//			split(vals, fields[++i], boost::is_any_of(";:"));
//			params.samp_vec.resize(vals.size());
//			for (int j = 0; j < vals.size(); j++)
//				params.samp_vec[j] = stoi(vals[j]);
//			if (vals.size() == 1 && params.samp_vec[0] <= 0) params.sampsize = NULL;
//			else params.sampsize = &params.samp_vec[0];
//		}
//
//	
//	}
//	MLOG("QRF init from string %s :: %d,%d\n", init_str.c_str(), params.type, params.ntrees);
//
//	return 0;
//}


//..............................................................................
int MedQRF::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [MedQRF::init]
		if (field == "ntrees") params.ntrees = stoi(entry.second);
		else if (field == "maxq") params.maxq = stoi(entry.second);
		else if (field == "type") params.type = get_tree_type(entry.second);
		else if (field == "min_node") params.min_node = stoi(entry.second);
		else if (field == "ntry") params.ntry = stoi(entry.second);
		else if (field == "get_count") params.get_count = stoi(entry.second);
		else if (field == "get_only_this_categ") params.get_only_this_categ = stoi(entry.second);
		else if (field == "max_samp") params.max_samp = stoi(entry.second);
		else if (field == "samp_factor") params.samp_factor = stof(entry.second);
		else if (field == "n_categ") params.n_categ = stoi(entry.second);
		else if (field == "spread") params.spread = stof(entry.second);
		else if (field == "learn_nthreads") params.learn_nthreads = stoi(entry.second);
		else if (field == "predict_nthreads") params.predict_nthreads = stoi(entry.second);
		else if (field == "keep_all_values") params.keep_all_values = (bool)(stoi(entry.second) != 0);
		else if (field == "max_depth") params.max_depth = stoi(entry.second);
		else if (field == "quantiles") {
			vector<string> vals;
			split(vals, entry.second, boost::is_any_of(","));
			params.quantiles.resize(vals.size());
			for (int j = 0; j < vals.size(); j++)
				params.quantiles[j] = stof(vals[j]);
		}
		else if (field == "sampsize") {
			vector<string> vals;
			split(vals, entry.second, boost::is_any_of(","));
			params.samp_vec.resize(vals.size());
			for (int j = 0; j < vals.size(); j++)
				params.samp_vec[j] = stoi(vals[j]);

			if (vals.size() == 1 && params.samp_vec[0] <= 0) params.sampsize = NULL;
			else {
				params.sampsize = &params.samp_vec[0];
			}
		}
		else MLOG("Unknown parameter \'%s\' for QRF\n", field.c_str());
		//! [MedQRF::init]

	}

	return 0;
}

//..............................................................................
QRF_TreeType MedQRF::get_tree_type(string name) {

	boost::algorithm::to_lower(name);
	if (name == "binary_tree" || name == "binary")
		return QRF_BINARY_TREE;
	else if (name == "regression_tree" || name == "regression")
		return QRF_REGRESSION_TREE;
	else if (name == "categorial_chi2_tree" || name == "categorical_chi2_tree" || name == "categorial_chi2" || name == "categorical_chi2")
		return QRF_CATEGORICAL_CHI2_TREE;
	else if (name == "categorial_entropy_tree" || name == "categorical_entropy_tree" || name == "categorial_entropy" || name == "categorical_entropy")
		return QRF_CATEGORICAL_ENTROPY_TREE;
	else
		return QRF_LAST;
}

//..............................................................................
void MedQRF::set_sampsize(float *y, int nsamples)
{
	if (params.sampsize != NULL) return;
	if (params.max_samp <= 0 && params.samp_factor <= 0) return;

	if (params.type == QRF_REGRESSION_TREE && params.max_samp <= 0) return;

	params.samp_vec.clear();

	if (params.type == QRF_REGRESSION_TREE) {
		params.samp_vec.push_back(params.max_samp);
		params.sampsize = &params.samp_vec[0];
		return;
	}

	params.samp_vec.resize(params.n_categ, 0);

	for (int i = 0; i < nsamples; i++)
		params.samp_vec[(int)y[i]]++;

	for (int i = 0; i < params.n_categ; i++)
		if (params.samp_vec[i] > params.max_samp) params.samp_vec[i] = params.max_samp;

	params.sampsize = &params.samp_vec[0];

	if (params.samp_factor <= 0) return;

	int max_ind = -1, max_val = 0;
	int max2_ind = -1, max2_val = 0;
	for (int i = 0; i < params.n_categ; i++)
		if (params.samp_vec[i] > max_val) {
			max_ind = i;
			max_val = params.samp_vec[i];
		}
		else if (params.samp_vec[i] > max2_val) {
			max2_ind = i;
			max2_val = params.samp_vec[i];
		}

		if ((float)max_val / (float)(1 + max2_val) > params.samp_factor)
			params.samp_vec[max_ind] = (int)(params.samp_factor*(float)max2_val);

}

//..............................................................................
MedQRF::MedQRF()
{
	init_defaults();
}

//..............................................................................
MedQRF::MedQRF(MedQRFParams& _in_params)
{
	init((void *)&_in_params);
}

//..............................................................................
MedQRF::MedQRF(void *_in_params)
{
	init(_in_params);
}

//..............................................................................
int MedQRF::Learn(float *x, float *y, float *w, int nsamples, int nftrs) {

	if (w != NULL && params.type != QRF_CATEGORICAL_ENTROPY_TREE)
		MWARN("Weights are not implemented for QRF that is not CATEGORICAL_ENTROPY. Ignoring\n");

	// Correct y according to forest type
	map<float, int> y_values;
	int n_categ = 0;
	qf.nthreads = params.learn_nthreads;
	qf.collect_oob = params.collect_oob;
	qf.get_only_this_categ = params.get_only_this_categ;
	qf.get_counts_flag = params.get_count;
	qf.keep_all_values = params.keep_all_values;
	qf.quantiles = params.quantiles;
	qf.max_depth = params.max_depth;

	if (params.type != QRF_REGRESSION_TREE) {
		for (int i = 0; i < nsamples; i++)
			y_values[y[i]] = 1;

		for (auto it = y_values.begin(); it != y_values.end(); it++)
			y_values[it->first] = n_categ++;

		if ((params.type == QRF_BINARY_TREE && n_categ != 2) || ((params.type == QRF_CATEGORICAL_CHI2_TREE || params.type == QRF_CATEGORICAL_ENTROPY_TREE) && n_categ < 2)) {
			MERR("Mismatch between QRF type and number of categories (%d)\n", qf.n_categ);
			return -1;
		}
		if (params.n_categ != n_categ) {
			MERR("Mismatch between requested n_categ: %d and actual n_categ: %d. Discovered following cateogories:\n", params.n_categ, n_categ);
			for (auto i : y_values)
				MERR("%f %d\n", i.first, i.second);
			return -1;
		}
	}

	if (params.type == QRF_CATEGORICAL_CHI2_TREE || params.type == QRF_CATEGORICAL_ENTROPY_TREE) {
		vector<float> qf_y(nsamples);
		for (int i = 0; i < nsamples; i++)
			qf_y[i] = (float)y_values[y[i]];

		if (qf.get_forest_categorical(x, &(qf_y[0]), w, nftrs, nsamples, params.sampsize, params.ntry, params.ntrees, params.maxq, params.min_node, n_categ, params.type) == -1) {
			MERR("Categorial QRF failed\n");
			return -1;
		}
	}
	else if (params.type == QRF_BINARY_TREE) {
		vector<int> qf_y(nsamples);
		for (int i = 0; i < nsamples; i++)
			qf_y[i] = y_values[y[i]];

		if (qf.get_forest(x, &(qf_y[0]), nftrs, nsamples, params.sampsize, params.ntry, params.ntrees, params.maxq) == -1) {
			MERR("Binary QRF failed\n");
			return -1;
		}
	}
	else { // REGRESSION
		int _sampsize = (params.sampsize == NULL) ? nsamples : params.sampsize[0];
		if (qf.get_forest_regression_trees(x, y, nftrs, nsamples, _sampsize, params.ntry, params.ntrees, params.maxq, params.min_node, params.spread) == -1) {
			MERR("Regression QRF failed\n");
			return -1;
		}
	}

	//	if (params.samp_vec.size() > 0) params.sampsize = NULL;

	return 0;
}

//..............................................................................
int MedQRF::Predict(float *x, float *&preds, int nsamples, int nftrs, int _get_count) {
	qf.nthreads = params.predict_nthreads;
	qf.get_only_this_categ = params.get_only_this_categ;
	qf.n_categ = params.n_categ;
	qf.get_only_this_categ = params.get_only_this_categ;
	qf.get_counts_flag = params.get_count;
	qf.quantiles = params.quantiles;
	return qf.score_samples(x, nftrs, nsamples, preds, _get_count);
}

//..............................................................................
int MedQRF::Predict(float *x, float *&preds, int nsamples, int nftrs) {
	return Predict(x, preds, nsamples, nftrs, params.get_count);
}

//..............................................................................
size_t MedQRF::get_size() {
	return qf.get_size()
		+ MedSerialize::get_size(model_features) + MedSerialize::get_size(features_count);
}

//..............................................................................
size_t MedQRF::serialize(unsigned char *blob) {
	size_t size = qf.serialize(blob);
	size += MedSerialize::serialize(blob + size, model_features);
	size += MedSerialize::serialize(blob + size, features_count);
	return size;
}

//..............................................................................
size_t MedQRF::deserialize(unsigned char *blob) {
	size_t s = qf.deserialize(blob);
	if ((int)s < 0)
		return -1;
	params.ntrees = (int)qf.qtrees.size();
	params.n_categ = qf.n_categ;
	params.type = (QRF_TreeType)qf.mode;
	params.get_only_this_categ = qf.get_only_this_categ;
	params.get_count = qf.get_counts_flag;
	params.quantiles = qf.quantiles;

	s += MedSerialize::deserialize(blob + s, model_features);
	s += MedSerialize::deserialize(blob + s, features_count);
	return s;
}

// Printing
void MedQRF::print(FILE *fp, const string& prefix) {
	fprintf(fp, "%s: MedQRF of type %d\n", prefix.c_str(), params.type);
	fprintf(fp, "%s: params = (ntrees=%d, maxq=%d, sampsize=%d, ntry=%d, spread=%f, min_node=%d, n_categ=%d, get_count=%d)\n",
		prefix.c_str(), params.ntrees, params.maxq, params.sampsize == NULL ? -1 : *params.sampsize, params.ntry, params.spread, params.min_node, params.n_categ, params.get_count);

	if (1) {
		for (unsigned int i = 0; i < qf.qtrees.size(); i++) {
			for (unsigned int j = 0; j < qf.qtrees[i].qnodes.size(); j++) {
				QRF_ResNode& node = qf.qtrees[i].qnodes[j];

				fprintf(fp, "Tree %d Node %d ", i, j);
				if (node.is_leaf)
					fprintf(fp, "Prediction %f\n", node.pred);
				else
					fprintf(fp, "Split by %d at %f to %d and %d\n", node.ifeat, node.split_val, node.left, node.right);
			}
		}
	}

}

string printNode(const vector<string> &modelSignalNames, const vector<QRF_ResNode> &nodes, int n, int deepSize, bool leftSide) {
	stringstream out;
	for (size_t i = 0; i < deepSize; ++i)
	{
		out << "\t";
	}
	if (leftSide) {
		out << "<=";
	}
	else {
		out << ">=";
	}
	if (!nodes[n].is_leaf) {
		out << modelSignalNames[nodes[n].ifeat];
	}

	if (nodes[n].counts.size() > 0) { //category, print histogram
		out << "(c0=" << nodes[n].counts[0];
		for (size_t i = 1; i < nodes[n].counts.size(); ++i)
		{
			out << ", c" << i << "=" << nodes[n].counts[i];
		}
	}
	else {
		out << "(pred=" << nodes[n].pred;
	}

	out << ", size=" << nodes[n].size;
	if (!nodes[n].is_leaf) {
		out << ", split=" << nodes[n].split_val << ")";
	}
	out << endl;
	return out.str();
}

void MedQRF::printTrees(const vector<string> &modelSignalNames, const string &outputPath) {
	stringstream treeOut;
	vector<QRF_ResTree> trees = qf.qtrees;
	for (size_t i = 0; i < trees.size(); ++i)
	{
		treeOut << "Tree_" << i << " {" << endl;
		vector<QRF_ResNode> nodes = trees[i].qnodes;
		//print only leaves with full Path
		vector<int> currPath = { 0 };
		//print Root
		treeOut << printNode(modelSignalNames, nodes, 0, 0, true);
		unordered_set<int> completed;

		while (completed.find(0) == completed.end()) { //haven't reached to finish root
			for (int j = (int)currPath.size() - 1; j >= 0; --j) //search in curretn path where to traverse down
			{
				if (!nodes[currPath[j]].is_leaf) {
					int candidateNode = nodes[currPath[j]].left;
					if (completed.find(candidateNode) == completed.end()) {
						currPath.push_back(candidateNode); //iterate left
						string nodeText = printNode(modelSignalNames, nodes, candidateNode, j + 1, true);
						treeOut << nodeText;
						break;
					}
					else {
						candidateNode = nodes[currPath[j]].right;
						if (completed.find(candidateNode) == completed.end()) {
							currPath.push_back(candidateNode); //iterate right
							string nodeText = printNode(modelSignalNames, nodes, candidateNode, j + 1, false);
							treeOut << nodeText;
							break;
						}
						else { //right and left are finished:
							   //mark current as finshed;
							completed.insert(currPath[j]);
							currPath.pop_back();
							break;
						}
					}
				}
				else {
					//mark branch completed, pop to parent
					completed.insert(currPath[j]);
					currPath.pop_back();
					break;
				}

			}
		}
		treeOut << "}";

		treeOut << endl;

		ofstream fw(outputPath);
		if (!fw.good())
			MTHROW_AND_ERR("IO Error: can't read \"%s\"\n", outputPath.c_str());
		fw << treeOut.str();
		fw.close();
	}
}

// Prdictions per sample
int MedQRF::n_preds_per_sample()
{
	if (params.type == QRF_REGRESSION_TREE) {
		if (params.get_count == PREDS_REGRESSION_QUANTILE || params.get_count == PREDS_REGRESSION_WEIGHTED_QUANTILE)
			return (int)params.quantiles.size();
		else
			return 1;
	}
	if (params.get_count == PREDS_CATEG_MAJORITY_AVG || params.get_count == PREDS_CATEG_AVG_PROBS || params.get_count == PREDS_CATEG_AVG_COUNTS)
		return 1;
	if (params.get_only_this_categ >= 0 && params.get_only_this_categ < params.n_categ)
		return 1;
	return (max(1, qf.n_categ));
}

void MedQRF::calc_feature_importance(vector<float> &features_importance_scores,
	const string &general_params) {
	if (qf.qtrees.empty())
		MTHROW_AND_ERR("ERROR:: Requested calc_feature_importance before running learn\n");
	vector<pair<short, double>> res;
	qf.variableImportance(res, model_features.empty() ? features_count : (int)model_features.size());

	features_importance_scores.resize((int)res.size());
	for (size_t i = 0; i < res.size(); ++i)
		features_importance_scores[res[i].first] = (float)res[i].second;
}