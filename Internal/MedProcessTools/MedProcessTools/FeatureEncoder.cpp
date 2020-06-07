#include "FeatureProcess.h"
#include <algorithm>
#include <random>
#include <ctime>

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//-------------------------------------------------------------------------------
// Feature Encoder
//-------------------------------------------------------------------------------

/// check if a set of features is affected by the current processor
//.......................................................................................
bool FeatureEncoder::are_features_affected(unordered_set<string>& out_req_features) {

	// If empty = all features are required
	if (out_req_features.empty())
		return true;
	
	// Otherwise - check in names
	for (string ftr : names) {
		if (out_req_features.find(ftr) != out_req_features.end())
			return true;
	}

	return false;
}

/// update sets of required as input according to set required as output to processor
//.......................................................................................
void FeatureEncoder::update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features) {

	if (are_features_affected(out_req_features))
		// If active, than everything before is required
		in_req_features.clear();
	else
		// If not, do nothing
		in_req_features = out_req_features;
}

//-------------------------------------------------------------------------------
// PCA Learn
//-------------------------------------------------------------------------------
int FeaturePCA::_learn(MedFeatures& features, unordered_set<int>& ids) 
{
	MedMat<float> x;

	features.get_as_matrix(x);

	MedMat<float> pca_base;
	vector<float> varsum;
	MLOG("Running PCA : x: %d x %d\n", x.nrows, x.ncols);
	MedPCA(x, pca_base, varsum); //save linear coef, save varsum for each feature
	MLOG("After PCA: Got base matrix: %d x %d , varsum %d\n", pca_base.nrows, pca_base.ncols, varsum.size());

	for (int i=0; i<params.pca_top; i++) {
		MLOG("PCA base %d : varsum %f\n", i, varsum[i]);
	}
	
	W.resize(params.pca_top, pca_base.ncols);

#pragma omp parallel for
	for (int i=0; i<params.pca_top; i++)
		for (int j=0; j<pca_base.ncols; j++)
			W(i, j) = pca_base(i, j);

	names.clear();
	for (int i=0; i<params.pca_top; i++)
		names.push_back("FTR_ENCODER_PCA_" + to_string(i));

	return 0;
}

#if 0
int FeaturePCA::_learn_old(MedFeatures& features, unordered_set<int>& ids) {
	MedMat<float> init_mat, cov_mat, pca_mat;
	vector<string> feat_names;
	if (!ids.empty()) {
		vector<int> row_ids;
		for (int i = 0; i < features.samples.size(); ++i)
			if (ids.find(features.samples[i].id) != ids.end())
				row_ids.push_back(i);
		features.get_as_matrix(init_mat, feat_names, row_ids);
	}
	else
		features.get_as_matrix(init_mat);
	//subsample base_mat if given:
	MedMat<float> *base_mat = &init_mat;
	MedMat<float> subsample_mat;
	if (params.subsample_count > 0 && params.subsample_count < init_mat.nrows) {
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<> random_index(0, init_mat.nrows - 1);
		vector<bool> seen_index(init_mat.nrows);
		subsample_mat.resize(params.subsample_count, init_mat.ncols);
		for (int i = 0; i < subsample_mat.nrows; ++i)
		{
			int sel = random_index(gen);
			while (seen_index[sel])
				sel = random_index(gen);
			seen_index[sel] = true;
			for (int j = 0; j < init_mat.ncols; ++j)
				subsample_mat(i, j) = init_mat(sel, j);
		}

		base_mat = &subsample_mat;
	}

	vector<float> varsum;
	cov_mat.resize(base_mat->ncols, base_mat->ncols);
	vector<double> mean_vec(base_mat->ncols);
	for (int i = 0; i < base_mat->ncols; ++i)
	{
		double s = 0;
		for (int k = 0; k < base_mat->nrows; ++k)
			s += (*base_mat)(k, i);
		s /= base_mat->nrows;
		mean_vec[i] = s;
	}
	for (int i = 0; i < base_mat->ncols; ++i)
	{
		for (int j = i; j < base_mat->ncols; ++j)
		{
			double s = 0;
			for (int k = 0; k < base_mat->nrows; ++k)
				s += ((*base_mat)(k, i) - mean_vec[i])*((*base_mat)(k, j) - mean_vec[j]);
			s /= base_mat->nrows;

			cov_mat(i, j) = float(s);
			cov_mat(j, i) = cov_mat(i, j);
		}
	}

	//can test for missing value - for speed up assume no missing values
	MedPCA(cov_mat, pca_mat, varsum); //save linear coef, save varsum for each feature
	//filter by thresholds and keep only relevant
	vector<pair<float, int>> eigen_to_index((int)varsum.size());
	for (size_t i = 0; i < varsum.size(); ++i)
	{
		eigen_to_index[i].first = varsum[i];
		eigen_to_index[i].second = (int)i;
	}
	//already sorted! this sort is wrong!:
	/*sort(eigen_to_index.begin(), eigen_to_index.end(),
		[](const pair<float, int> &v1, const pair<float, int> &v2)
	{return (v1.first > v2.first); });*/

	//choose top :
	int final_size = 0;
	vector<bool> selected_index((int)eigen_to_index.size());
	for (size_t i = 0; i < eigen_to_index.size() &&
		(params.pca_top == 0 || i < params.pca_top); ++i)
		if (params.pca_cutoff == 0 || eigen_to_index[i].first >= params.pca_cutoff) {
			++final_size;
			selected_index[i] = true;
			selected_indexes.push_back((int)i);
		}

	W.resize(final_size, (int)eigen_to_index.size());
	//update features:
	int w_i = 0;
	for (int i = 0; i < eigen_to_index.size(); ++i)
		if (selected_index[i]) {
			for (int j = 0; j < eigen_to_index.size(); ++j)
				W(w_i, j) = pca_mat(i, j);
			++w_i;
		}

	names.clear();
	for (size_t i = 0; i < final_size; ++i)
		names.push_back("FTR_ENCODER_PCA_" + to_string(i));

	return 0;
}
#endif

int FeaturePCA::_apply(MedFeatures& features, unordered_set<int>& ids) {
	MedMat<float> x;
	//ids not supported
	features.get_as_matrix(x);
	//apply multiply by W: and add to features

#pragma omp parallel for
	for (int pca_cnt=0; pca_cnt<names.size(); pca_cnt++)
	{
		string name = names[pca_cnt];
#pragma omp critical 
		{
			features.data[name].resize((int)features.samples.size());
			features.attributes[name].imputed = true;
			features.tags[name].insert("pca_encoder");
		}
		//do multiply:
		float *datap = &(features.data[name][0]);
		for (int i=0; i<x.nrows; i++) {
			datap[i] = 0;
			for (int j=0; j<x.ncols; j++)
				datap[i] += x(i, j)*W(pca_cnt, j);
		}
	}

	return 0;
}

int FeaturePCA::init(map<string, string>& mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		//! [FeaturePCA::init]
		if (field == "pca_cutoff") params.pca_cutoff = med_stof(entry.second);
		else if (field == "pca_top") params.pca_top = med_stoi(entry.second);
		else if (field == "subsample_count") params.subsample_count = med_stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeaturePCA\n", field.c_str());
		//! [FeaturePCA::init]
	}

	return 0;
}
