#include "FeatureProcess.h"
#include <algorithm>

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int FeatureEncoder::Learn(MedFeatures& features, unordered_set<int>& ids) {
	// learn encoder...
	if (_learn(features, ids) < 0)
		return -1;
	MLOG("Learn Feature Encoder: Generated %d features\n", (int)names.size());

	return 0;
}

int FeatureEncoder::Apply(MedFeatures& features, unordered_set<int>& ids) {
	if (_apply(features, ids) < 0)
		return -1;
	MLOG("Apply Feature Encoder: Generated %d features\n", (int)names.size());
	return 0;
}

int FeaturePCA::_learn(MedFeatures& features, unordered_set<int>& ids) {
	MedMat<float> base_mat, cov_mat, pca_mat;
	vector<string> feat_names;
	vector<int> row_ids;
	for (int i = 0; i < features.samples.size(); ++i)
		if (ids.find(features.samples[i].id) != ids.end())
			row_ids.push_back(i);
	features.get_as_matrix(base_mat, feat_names, row_ids);
	vector<float> varsum;
	cov_mat.resize(base_mat.ncols, base_mat.ncols);
	vector<double> mean_vec(base_mat.ncols);
	for (int i = 0; i < base_mat.ncols; ++i)
	{
		double s = 0;
		for (int k = 0; k < base_mat.nrows; ++k)
			s += base_mat(k, i);
		s /= base_mat.nrows;
		mean_vec[i] = s;
	}
	for (int i = 0; i < base_mat.ncols; ++i)
	{
		for (int j = i; j < base_mat.ncols; ++j)
		{
			double s = 0;
			for (int k = 0; k < base_mat.nrows; ++k)
				s += base_mat(k, i)*base_mat(k, j);
			s /= base_mat.nrows;
			s -= mean_vec[i] * mean_vec[j];

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
	sort(eigen_to_index.begin(), eigen_to_index.end(),
		[](const pair<float, int> &v1, const pair<float, int> &v2)
	{return (v1.first > v2.first); });

	//choose top :
	int final_size = 0;
	vector<bool> selected_index((int)eigen_to_index.size());
	for (size_t i = 0; i < eigen_to_index.size() &&
		(params.pca_top == 0 || i < params.pca_top); ++i)
		if (params.pca_cutoff == 0 || eigen_to_index[i].first >= params.pca_cutoff) {
			++final_size;
			selected_index[i] = true;
		}

	W.resize(final_size, final_size);
	//update features:
	int w_i = 0, w_j;
	for (int i = 0; i < eigen_to_index.size(); ++i)
		if (selected_index[i]) {
			w_j = 0;
			for (int j = 0; j < eigen_to_index.size(); ++j)
				if (selected_index[j])
				{
					W(w_i, w_j) = pca_mat(i, j);
					++w_j;
				}
			++w_i;
		}

	names.clear();
	for (size_t i = 0; i < final_size; ++i)
		names.push_back("FTR_ENCODER_PCA_" + to_string(i));

	return 0;
}

int FeaturePCA::_apply(MedFeatures& features, unordered_set<int>& ids) {
	MedMat<float> base_mat;
	//ids not supported
	features.get_as_matrix(base_mat);
	//apply multiply by W: and add to features
	int pca_cnt = 0;
	for (string name : names)
	{
		features.data[name].resize((int)features.samples.size());
		features.attributes[name].imputed = true;
		features.tags[name].insert("pca_encoder");
		//do multiply:
		for (int i = 0; i < features.samples.size(); ++i)
			for (int j = 0; j < base_mat.ncols; ++j)
				features.data[name][i] += base_mat(i, j)*W(pca_cnt, j);

		++pca_cnt;
	}

	return 0;
}

int FeaturePCA::init(map<string, string>& mapper) {
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "pca_cutoff") params.pca_cutoff = stof(entry.second);
		else if (field == "pca_top") params.pca_top =stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeaturePCA\n", field.c_str());
	}

	return 0;
}