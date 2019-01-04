#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "FeatureProcess.h"
#include "DoCalcFeatProcessor.h"
#include <MedStat/MedStat/MedPerformance.h>
#include <omp.h>
#include <boost/regex.hpp>
#include "MedStat/MedStat/MedBootstrap.h"

//=======================================================================================
// FeatureSelector
//=======================================================================================
// Learn : Add required to feature selected by inheriting classes
//.......................................................................................
int FeatureSelector::Learn(MedFeatures& features, unordered_set<int>& ids) {
	//MLOG("FeatureSelector::Learn %d features\n", features.data.size());

	// select, possibly ignoring requirments
	if (_learn(features, ids) < 0)
		return -1;

	// Add required signals
	// Collect selected
	vector<string> ftrs;
	features.get_feature_names(ftrs);

	unordered_set<string> selectedFeatures;
	for (string& feature : selected)
		selectedFeatures.insert(feature);

	// Find Missing
	vector<string> missingRequired;
	for (string feature : required) {
		string resolved = resolve_feature_name(features, feature);
		if (selectedFeatures.find(resolved) == selectedFeatures.end()) {
			missingRequired.push_back(resolved);
			MLOG("FeatureSelector::Learn re-inserting removed feature [%s] because it is required!\n", resolved.c_str());
		}
	}

	// Keep maximum numToSelect ...
	if (numToSelect > 0) {
		int nMissing = (int)missingRequired.size();
		int nSelected = (int)selected.size();

		if (nSelected + nMissing < numToSelect)
			selected.resize(nSelected + nMissing, "");
		else
			selected.resize(numToSelect, "");
	}

	// Insert (making sure not to remove required features)
	int insertIndex = (int)selected.size() - 1;
	for (unsigned int i = 0; i < missingRequired.size(); i++) {
		while (required.find(selected[insertIndex]) != required.end()) {
			insertIndex--;
			assert(insertIndex >= 0);
		}
		selected[insertIndex--] = missingRequired[i];
	}

	// Log
	//for (string& feature : selected)
		//MLOG("Feature Selection: Selected %s\n", feature.c_str());

	return 0;
}

// Apply selection : Ignore set of ids
//.......................................................................................
int FeatureSelector::_apply(MedFeatures& features, unordered_set<int>& ids) {

	unordered_set<string> empty_set;
	return _conditional_apply(features, ids, empty_set);
}

//.......................................................................................
int FeatureSelector::_conditional_apply(MedFeatures& features, unordered_set<int>& ids, unordered_set<string>& out_req_features) {

	//MLOG("FeatureSelector::Apply %d features\n", features.data.size());
	unordered_set<string> selectedFeatures;
	for (string& feature : selected) {
		if (out_req_features.empty() || (out_req_features.find(feature) != out_req_features.end()))
			selectedFeatures.insert(feature);
	}

	return features.filter(selectedFeatures);
}

/// check if a set of features is affected by the current processor
//.......................................................................................
bool FeatureSelector::are_features_affected(unordered_set<string>& out_req_features) {

	// Always true . The selected features are all that's left so they must be in the out_req_features
	return true;
}

/// update sets of required as input according to set required as output to processor
//.......................................................................................
void FeatureSelector::update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features) {

	// Only selected features are required ...
	in_req_features.clear();
	for (string ftr : selected) {
		if (out_req_features.empty() || out_req_features.find(ftr) != out_req_features.end())
			in_req_features.insert(ftr);
	}
}

//=======================================================================================
// UnivariateFeatureSelector
//=======================================================================================
// Learn 
//.......................................................................................
int UnivariateFeatureSelector::_learn(MedFeatures& features, unordered_set<int>& ids) {

	// Get Stats
	vector<float> stats;

	// "Correlation" to outcome
	if (params.method == UNIV_SLCT_PRSN)
		getAbsPearsonCorrs(features, ids, stats);
	else if (params.method == UNIV_SLCT_MI) {
		if (getMIs(features, ids, stats) < 0)
			return -1;
	}
	else if (params.method == UNIV_SLCT_DCORR) {
		if (getDistCorrs(features, ids, stats) < 0)
			return -1;
	}
	else {
		MERR("Unknown method %d for univariate feature selection\n", params.method);
		return -1;
	}

	// Select
	vector<pair<string, float >> namedStats(stats.size());
	vector<string> names(stats.size());
	features.get_feature_names(names);
	for (int i = 0; i < names.size(); i++) {
		namedStats[i].first = names[i];
		namedStats[i].second = stats[i];
	}

	sort(namedStats.begin(), namedStats.end(), [](const pair<string, float> &v1, const pair<string, float> &v2) {return (v1.second > v2.second); });

	if (numToSelect == 0) {
		// Select according to minimum value of stat
		for (auto& rec : namedStats) {
			if (rec.second < params.minStat)
				break;
			selected.push_back(rec.first);
		}
	}
	else {
		// Select according to number
		int n = (namedStats.size() > numToSelect) ? numToSelect : (int)namedStats.size();
		selected.resize(n);
		for (int i = 0; i < n; i++)
			selected[i] = namedStats[i].first;
	}

	return 0;
}

// Init
//.......................................................................................
int UnivariateFeatureSelector::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [UnivariateFeatureSelector::init]
		if (field == "missing_value") missing_value = med_stof(entry.second);
		else if (field == "numToSelect") numToSelect = med_stoi(entry.second);
		else if (field == "method") params.method = params.get_method(entry.second);
		else if (field == "minStat") params.minStat = med_stof(entry.second);
		else if (field == "nBins") params.nBins = med_stoi(entry.second);
		else if (field == "binMethod") params.binMethod = params.get_binning_method(entry.second);
		else if (field == "required") boost::split(required, entry.second, boost::is_any_of(","));
		else if (field == "takeSquare") params.takeSquare = med_stoi(entry.second);
		else if (field == "max_samples") params.max_samples = med_stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeatureSelector\n", field.c_str());
		//! [UnivariateFeatureSelector::init]
	}

	return 0;

}

// Utility : Caluclate pearson correlations to a vector
//.......................................................................................
int UnivariateFeatureSelector::getAbsPearsonCorrs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats) {

	// Get outcome
	vector<float> label;
	get_all_outcomes(features, ids, label, params.max_samples);

	int nFeatures = (int)features.data.size();
	vector<string> names(nFeatures);
	features.get_feature_names(names);
	stats.resize(nFeatures);

#pragma omp parallel for 
	for (int i = 0; i < nFeatures; i++) {
		int n;
		vector<float> values;
		get_all_values(features, names[i], ids, values, params.max_samples);

		stats[i] = fabs(medial::performance::pearson_corr(values, label, missing_value, n));
		if (n == 0) stats[i] = 0.0;

		// If required, test also correlation to squared (normalized) values
		if (params.takeSquare) {
			float mean, std;
			int n_cleaned;
			medial::stats::get_mean_and_std(values, missing_value, n, mean, std);
			for (unsigned int j = 0; j < values.size(); j++) {
				if (values[j] != missing_value)
					values[j] = (values[j] - mean)*(values[j] - mean);
			}

			float newStat = fabs(medial::performance::pearson_corr(values, label, missing_value, n_cleaned));
			if (newStat > stats[i])
				stats[i] = newStat;
		}
	}

	return 0;

}

// Utility : Caluclate Mutual Information
//.......................................................................................
int UnivariateFeatureSelector::getMIs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats) {

	// Get outcome
	vector<float> label;
	get_all_outcomes(features, ids, label, params.max_samples);

	vector<int> binnedLabel;
	int nBins;
	if (discretize(label, binnedLabel, nBins, params.nBins, missing_value, params.binMethod) < 0)
		return -1;

	size_t nFeatures = features.data.size();
	vector<string> names(nFeatures);
	features.get_feature_names(names);
	stats.resize(nFeatures);
	vector<vector<int>> binnedValues(nFeatures);

	int RC = 0;
#pragma omp parallel for 
	for (int i = 0; i < names.size(); i++) {
		vector<float> values;
		int nBins;
		get_all_values(features, names[i], ids, values, params.max_samples);
		int rc = discretize(values, binnedValues[i], nBins, params.nBins, missing_value, params.binMethod);
#pragma omp critical
		if (rc < 0)  RC = -1;
	}

#pragma omp parallel for 
	for (int i = 0; i < names.size(); i++) {
		int n;
		stats[i] = medial::performance::mutual_information(binnedValues[i], binnedLabel, n);
		if (stats[i] < 0) stats[i] = 0;
	}

	return 0;

}

// Utility : Caluclate distance correlations
//.......................................................................................
int UnivariateFeatureSelector::getDistCorrs(MedFeatures& features, unordered_set<int>& ids, vector<float>& stats) {

	// Get outcome
	vector<float> label;
	get_all_outcomes(features, ids, label, params.max_samples);

	MedMat<float> labelDistances;
	medial::performance::get_dMatrix(label, labelDistances, missing_value);
	float targetDistVar = medial::performance::get_dVar(labelDistances);
	if (targetDistVar == -1.0) {
		MERR("Cannot calucludate distance Var for target\n");
		return -1;
	}

	size_t nFeatures = features.data.size();
	vector<string> names(nFeatures);
	features.get_feature_names(names);
	stats.resize(nFeatures);

	int RC = 0;
#pragma omp parallel for 
	for (int i = 0; i < names.size(); i++) {

		vector<float> values;
		get_all_values(features, names[i], ids, values, params.max_samples);
		MedMat<float> valueDistances;
		medial::performance::get_dMatrix(values, valueDistances, missing_value);
		float valueDistVar = medial::performance::get_dVar(valueDistances);
		float distCov = medial::performance::get_dCov(labelDistances, valueDistances);
#pragma omp critical
		if (valueDistVar == -1 || distCov == -1) {
			MERR("Cannot calculate distance correlation between label and %s\n", names[i].c_str());
			RC = -1;
		}
		else {
			stats[i] = distCov / sqrt(valueDistVar*targetDistVar);
		}
	}

	return RC;

}

//=======================================================================================
// MRMRFeatureSelector
//=======================================================================================
// Learn 
//.......................................................................................
int MRMRFeatureSelector::_learn(MedFeatures& features, unordered_set<int>& ids) {

	if (numToSelect == 0) {
		MERR("MRMR requires numToSelect>0");
		return -1;
	}

	int nFeatures = (int)features.data.size();
	vector<string> names(nFeatures);
	features.get_feature_names(names);

	// Start filling "Correlation" matrix
	MedMat<float> stats(nFeatures + 1, nFeatures + 1);
	for (int i = 0; i <= nFeatures; i++) {
		for (int j = 0; j <= nFeatures; j++) {
			stats(i, j) = stats(j, i) = -1;
		}
		stats(i, i) = 0;
	}

	if (fillStatsMatrix(features, ids, stats, nFeatures) < 0)
		return -1;

	// Actual selection
	vector <int> selectedIds;
	vector<int> selectFlags(nFeatures, 0);

	for (int iSelect = 0; iSelect < numToSelect; iSelect++) {
		float optScore = missing_value;
		int optFeature = -1;
		for (int i = 0; i < nFeatures; i++) {
			if (selectFlags[i] == 0) {
				float score = stats(i, nFeatures);
				if (iSelect > 0) {
					float penaltyValue = 0.0;
					if (penaltyMethod == MRMR_MAX) {
						for (int j = 0; j < iSelect; j++) {
							if (stats(i, selectedIds[j]) > penaltyValue) {
								penaltyValue = stats(i, selectedIds[j]);
							}
						}
					}
					else if (penaltyMethod = MRMR_MEAN) {
						for (int j = 0; j < iSelect; j++)
							penaltyValue += stats(i, selectedIds[j]);
						penaltyValue /= iSelect;
					}

					score -= penalty * penaltyValue;
				}

				if (optFeature == -1 || score > optScore) {
					optScore = score;
					optFeature = i;
				}
			}
		}
		selectedIds.push_back(optFeature);
		selectFlags[optFeature] = 1;
		fillStatsMatrix(features, ids, stats, optFeature);
	}

	selected.clear();
	for (int id : selectedIds) selected.push_back(names[id]);
	return 0;
}

// Init
//.......................................................................................
int MRMRFeatureSelector::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [MRMRFeatureSelector::init]
		if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "numToSelect") numToSelect = stoi(entry.second);
		else if (field == "method") params.method = params.get_method(entry.second);
		else if (field == "minStat") params.minStat = stof(entry.second);
		else if (field == "nBins") params.nBins = stoi(entry.second);
		else if (field == "binMethod") params.binMethod = params.get_binning_method(entry.second);
		else if (field == "required") boost::split(required, entry.second, boost::is_any_of(","));
		else if (field == "penalty") penalty = stof(entry.second);
		else if (field == "penaltyMethod") penaltyMethod = get_penalty_method(entry.second);
		else if (field == "max_samples") params.max_samples = stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeatureSelector\n", field.c_str());
		//! [MRMRFeatureSelector::init]
	}

	return 0;

}

//.......................................................................................
MRMRPenaltyMethod MRMRFeatureSelector::get_penalty_method(string _method) {

	boost::to_lower(_method);
	if (_method == "max")
		return MRMR_MAX;
	else if (_method == "mean")
		return MRMR_MEAN;
	else
		return MRMR_LAST;

}

//.......................................................................................
void MRMRFeatureSelector::init_defaults() {
	missing_value = MED_MAT_MISSING_VALUE;
	processor_type = FTR_PROCESSOR_MRMR_SELECTOR;
	params.method = UNIV_SLCT_PRSN;
	numToSelect = 50;
	penaltyMethod = MRMR_MAX;
	penalty = 0.5;
}

// Utility : Caluclate  correlations
//.......................................................................................
int MRMRFeatureSelector::fillStatsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index) {

	if (params.method == UNIV_SLCT_PRSN)
		fillAbsPearsonCorrsMatrix(features, ids, stats, index);
	else if (params.method == UNIV_SLCT_MI) {
		if (fillMIsMatrix(features, ids, stats, index) < 0)
			return -1;
	}
	else if (params.method == UNIV_SLCT_DCORR) {
		if (fillDistCorrsMatrix(features, ids, stats, index) < 0)
			return -1;
	}
	else {
		MERR("Unknown method %d for univariate feature selection\n", params.method);
		return -1;
	}

	return 0;
}

// Utility : Caluclate pearson correlations
//.......................................................................................
int MRMRFeatureSelector::fillAbsPearsonCorrsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index) {

	int nFeatures = (int)features.data.size();
	vector<string> names(nFeatures);
	features.get_feature_names(names);
	vector<vector<float>> values(nFeatures);

	// Get outcome
	vector<float> target;
	if (index == nFeatures)
		get_all_outcomes(features, ids, target, params.max_samples);
	else
		get_all_values(features, names[index], ids, target, params.max_samples);

#pragma omp parallel for 
	for (int i = 0; i < nFeatures; i++) {
		if (stats(i, index) == -1)
			get_all_values(features, names[i], ids, values[i], params.max_samples);
	}

#pragma omp parallel for 
	for (int i = 0; i < nFeatures; i++) {
		if (stats(i, index) == -1) {
			int n;
			stats(i, index) = fabs(medial::performance::pearson_corr(values[i], target, missing_value, n));
			if (n == 0) stats(i, index) = 0.0;
			stats(index, i) = stats(i, index);
		}
	}

	return 0;
}

// Utility : Caluclate Mutual Information
//.......................................................................................
int MRMRFeatureSelector::fillMIsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index) {

	size_t nFeatures = features.data.size();
	vector<string> names(nFeatures);
	features.get_feature_names(names);
	vector<vector<int>> binnedValues(nFeatures);

	// Get outcome
	vector<float> target;
	if (index == nFeatures)
		get_all_outcomes(features, ids, target, params.max_samples);
	else
		get_all_values(features, names[index], ids, target, params.max_samples);

	vector<int> binnedTarget;
	int nBins;
	if (discretize(target, binnedTarget, nBins, params.nBins, missing_value, params.binMethod) < 0)
		return -1;
	if (nBins < params.nBins)
		smearBins(binnedTarget, nBins, params.nBins);

	int RC = 0;
#pragma omp parallel for 
	for (int i = 0; i < nFeatures; i++) {
		if (stats(i, index) == -1) {
			vector<float> values;
			int nBins;
			get_all_values(features, names[i], ids, values, params.max_samples);
			int rc = discretize(values, binnedValues[i], nBins, params.nBins, missing_value, params.binMethod);
#pragma omp critical
			if (rc < 0)  RC = -1;

			if (nBins < params.nBins)
				smearBins(binnedValues[i], nBins, params.nBins);
		}
	}

#pragma omp parallel for 
	for (int i = 0; i < nFeatures; i++) {
		if (stats(i, index) == -1) {
			int n;
			stats(i, index) = medial::performance::mutual_information(binnedValues[i], binnedTarget, n);
			if (stats(i, index) < 0) stats(i, index) = 0;
			stats(index, i) = stats(i, index);
		}
	}

	return 0;

}

// Utility : Caluclate distance correlations
//.......................................................................................
int MRMRFeatureSelector::fillDistCorrsMatrix(MedFeatures& features, unordered_set<int>& ids, MedMat<float>& stats, int index) {

	size_t nFeatures = features.data.size();
	vector<string> names(nFeatures);
	features.get_feature_names(names);

	// Get outcome
	vector<float> target;
	if (index == nFeatures)
		get_all_outcomes(features, ids, target, params.max_samples);
	else
		get_all_values(features, names[index], ids, target, params.max_samples);

	MedMat<float> targetDistances;
	medial::performance::get_dMatrix(target, targetDistances, missing_value);
	float targetDistVar = medial::performance::get_dVar(targetDistances);
	if (targetDistVar == -1.0) {
		MERR("Cannot calucludate distance Var for target\n");
		return -1;
	}

	int RC = 0;
#pragma omp parallel for 
	for (int i = 0; i < nFeatures; i++) {
		if (stats(i, index) == -1) {
			vector<float> values;
			get_all_values(features, names[i], ids, values, params.max_samples);
			MedMat<float> valueDistances;
			medial::performance::get_dMatrix(values, valueDistances, missing_value);
			float valueDistVar = medial::performance::get_dVar(valueDistances);
			float distCov = medial::performance::get_dCov(targetDistances, valueDistances);
#pragma omp critical
			if (valueDistVar == -1 || distCov == -1) {
				MERR("Cannot calculate distance correlation between label and %s\n", names[i].c_str());
				RC = -1;
			}
			else {
				stats(index, i) = stats(i, index) = distCov / sqrt(valueDistVar*targetDistVar);
			}
		}
	}

	return RC;

}

//=======================================================================================
// Lasso Feature Selection
//=======================================================================================
// Learn 
//.......................................................................................
int LassoSelector::_learn(MedFeatures& features, unordered_set<int>& ids) {

	vector<string> names;
	features.get_feature_names(names);
	int nFeatures = (int)names.size();

	// Required features
	vector<int> lax_indices;
	unordered_set<string> lax_lasso_features_set(lax_lasso_features.begin(), lax_lasso_features.end());
	for (int i = 0; i < nFeatures; i++) {
		if (lax_lasso_features_set.find(names[i]) != lax_lasso_features_set.end())
			lax_indices.push_back(i);
	}

	if (numToSelect > nFeatures)
		MTHROW_AND_ERR("Cannot select %d out of %d", numToSelect, nFeatures);

	// Labels
	MedMat<float> y((int)features.samples.size(), 1);
	for (int i = 0; i < y.nrows; i++)
		y(i, 0) = features.samples[i].outcome;
	//	y.normalize();

		// Matrix
	MedMat<float> x;
	features.get_as_matrix(x);
	x.missing_value = missing_value;
	vector<float> avg, std;
	x.get_cols_avg_std(avg, std);
	x.normalize(avg, std, 1);

	// Initialize
	int found = 0;
	vector< vector<float> > lambdas(nthreads);
	vector<float> base_lambdas(nthreads);
	float minLambda = 0.0, maxLambda = initMaxLambda;
	//	vector<MedLasso> predictors(nthreads);
	vector<MedGDLM> predictors(nthreads);
	vector<int> nSelected(nthreads);
	vector<float> w(x.nrows, 1.0);

	MLOG("Lasso Feature Selection : From %d to [ %d , %d ] \n", nFeatures, numToSelect - numToSelectDelta, numToSelect + numToSelectDelta);
	selected.clear();
	float lowerBoundLambda = 0.0, upperBoundLambda = -1.0;

	// Prevent being stuck ...
	int nStuck = 0;
	int prevMaxN = -1, prevMinN = -1;

	// Search
	float upperBound = -1;
	while (!found) {
		if (nthreads == 1) {
			base_lambdas[0] = maxLambda;
			lambdas[0].assign(nFeatures, maxLambda);
		}
		else {
			float step = (maxLambda - minLambda) / (nthreads - 1);
			for (int i = 0; i < nthreads; i++) {
				base_lambdas[i] = minLambda + step * i;
				lambdas[i].assign(nFeatures, base_lambdas[i]);
				for (int idx : lax_indices)
					lambdas[i][idx] = base_lambdas[i] * lambdaRatio;
			}
		}

		for (int i = 0; i < nthreads; i++) {
			predictors[i].init_from_string("method = logistic_sgd; last_is_bias = 0; batch_size = 2048; rate_decay = 1; l_ridge = 0; ; err_freq = 10; nthreads = 12;l_lasso = 0" + string(";stop_at_err = ") + to_string(stop_at_err) + string(";rate = ") + to_string(rate) + string(";momentum = ") + to_string(momentum));
			predictors[i].params.ls_lasso = lambdas[i];
		}

#pragma omp parallel for 
		for (int i = 0; i < nthreads; i++) {
			predictors[i].learn(x, y);

			// Identify non-zero parameters
			nSelected[i] = 0;
			for (int j = 0; j < nFeatures; j++) {
				if (predictors[i].b[j] != 0)
					nSelected[i] ++;
			}
		}

		MLOG_V("Lasso Feature Selection: [%f,%f] : nFeatures [%d,%d] nStuck %d\n", base_lambdas[0], base_lambdas[nthreads - 1], nSelected[0], nSelected[nthreads - 1], nStuck);

		if (nthreads == 1) { // Special care

			if ((nSelected[0] >= numToSelect - numToSelectDelta) && (nSelected[0] <= numToSelect + numToSelectDelta)) {
				found = 1;
				for (int j = 0; j < nFeatures; j++) {
					if (predictors[0].b[j] != 0)
						selected.push_back(names[j]);
				}
			}
			else if (nSelected[0] > numToSelect + numToSelectDelta) {
				lowerBoundLambda = maxLambda;
				if (upperBoundLambda != -1.0)
					maxLambda = (upperBoundLambda + maxLambda) / (float)2.0;
				else
					maxLambda = maxLambda * (float)2.0;
			}
			else {
				upperBoundLambda = maxLambda;
				maxLambda = (lowerBoundLambda + maxLambda) / 2.0f;
			}
			minLambda = maxLambda;
		}
		else {
			for (int j = 0; j < nthreads; j++)
				MLOG("N[%.12f] = %d\n", base_lambdas[j], nSelected[j]);

			// float ratio;
			if (nSelected[nthreads - 1] > numToSelect) { // MaxLambda is still too low
				minLambda = maxLambda;
				if (upperBound == -1)
					maxLambda *= 2.0;
				else
					maxLambda = upperBound;
			}
			else { // in between ...
				for (int i = 0; i < nthreads; i++) {
					if (nSelected[i] == numToSelect) {
						found = 1;
						for (int j = 0; j < nFeatures; j++) {
							if (predictors[i].b[j] != 0)
								selected.push_back(names[j]);
						}
						break;
					}
					else if (nSelected[i] < numToSelect) {
						minLambda = base_lambdas[i - 1];
						upperBound = base_lambdas[i];

						// take care of nthreads = 2  
						if (nthreads == 2)
							maxLambda = base_lambdas[i - 1] + (base_lambdas[i] - base_lambdas[i - 1]) / 2.0;
						else
							maxLambda = base_lambdas[i];
						break;
					}
				}
			}
		}

		// Are We stuck ?
		if (nSelected[0] == prevMaxN && nSelected[nthreads - 1] == prevMinN) {
			nStuck++;
			if (nStuck == 3) {

				int minDiff = nFeatures;
				int optimalI = 0;
				for (int i = 0; i < nthreads; i++) {
					if (abs(nSelected[i] - numToSelect) < minDiff) {
						minDiff = abs(nSelected[i] - numToSelect);
						optimalI = i;
					}
				}

				found = 1;
				MLOG("Stuck at same N range for 3 steps. That's enough for now ... Actual NSelected = %d\n", nSelected[optimalI]);
				for (int j = 0; j < nFeatures; j++) {
					if (predictors[optimalI].b[j] != 0)
						selected.push_back(names[j]);
				}
			}
		}
		else {
			nStuck = 0;
			prevMaxN = nSelected[0];
			prevMinN = nSelected[nthreads - 1];
		}


		if (!found)
			MLOG("Lasso Feature Selection: about to try lambdas [%f,%f]\n", minLambda, maxLambda);

	}

	for (auto &f : selected)
		MLOG("Lasso Selected Feature %s \n", f.c_str());


	return 0;
}

// Init
//.......................................................................................
int LassoSelector::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [LassoSelector::init]
		if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "numToSelect") numToSelect = stoi(entry.second);
		else if (field == "numToSelectDelta") numToSelectDelta = stoi(entry.second);
		else if (field == "initMaxLambda") initMaxLambda = stof(entry.second);
		else if (field == "nthreads") nthreads = stoi(entry.second);
		else if (field == "required") boost::split(required, entry.second, boost::is_any_of(","));
		else if (field == "lax_lasso_features") boost::split(lax_lasso_features, entry.second, boost::is_any_of(","));
		else if (field == "lambdaRatio") lambdaRatio = stof(entry.second);
		else if (field == "rate") rate = stof(entry.second);
		else if (field == "momentum") momentum = stof(entry.second);
		else if (field == "stop_at_err") stop_at_err = stof(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeatureSelector\n", field.c_str());
		//! [LassoSelector::init]
	}

	return 0;

}


//=======================================================================================
// Remove Degenerate Features
//=======================================================================================
// Learn 
//.......................................................................................
int DgnrtFeatureRemvoer::_learn(MedFeatures& features, unordered_set<int>& ids) {

	selected.clear();

	for (auto& rec : features.data) {
		string name = rec.first;
		vector<float>& data = rec.second;

		map<float, int> counters;
		for (float& val : data)
			counters[val] ++;

		int maxCount = 0;
		float maxCountValue = missing_value;
		for (auto rec : counters) {
			if (rec.second > maxCount) {
				maxCount = rec.second;
				maxCountValue = rec.first;
			}
		}

		float p = ((float)maxCount) / (float)data.size();
		if (p >= percentage)
			MLOG("DgnrtFeatureRemvoer::_learn removing %s : %f of values is %f\n", name.c_str(), p, maxCountValue);
		else
			selected.push_back(name);
	}

	return 0;
}

// Init
//.......................................................................................
int DgnrtFeatureRemvoer::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [DgnrtFeatureRemvoer::init]
		if (field == "missing_value") missing_value = stof(entry.second);
		if (field == "percentage") percentage = stof(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeatureSelector\n", field.c_str());
		//! [DgnrtFeatureRemvoer::init]
	}

	assert(percentage >= 0 && percentage <= 1.0);
	return 0;

}

int TagFeatureSelector::init(map<string, string>& mapper) {
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		string val = boost::trim_copy(entry.second);
		//! [TagFeatureSelector::init]
		if (field == "missing_value") missing_value = med_stof(entry.second);
		else if (field == "selected_tags") { if (!val.empty()) boost::split(selected_tags, val, boost::is_any_of(",")); }
		else if (field == "removed_tags") { if (!val.empty()) boost::split(removed_tags, val, boost::is_any_of(",")); }
		else if (field == "verbose") verbose = med_stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for TagFeatureSelector\n", field.c_str());
		//! [TagFeatureSelector::init]
	}

	return 0;
}

int TagFeatureSelector::_learn(MedFeatures& features, unordered_set<int>& ids) {
	selected.clear();
	//fix clean names
	unordered_set<string> s(selected_tags.begin(), selected_tags.end());
	unordered_set<string> r(removed_tags.begin(), removed_tags.end());

	if (r.empty())
		MLOG("TagFeatureSelector not removing any features\n");
	else
		for (string sub : r)
			if (verbose)
				MLOG("TagFeatureSelector removing features with tag [%s]\n", sub.c_str());
	if (s.empty())
		MLOG("TagFeatureSelector selecting all features\n");
	else
		for (string sub : s)
			if (verbose)
				MLOG("TagFeatureSelector selecting features with tag [%s]\n", sub.c_str());
	int removed_features = 0, selected_features = 0;
	for (auto it = features.tags.begin(); it != features.tags.end(); ++it) {
		string feature_name = it->first;
		unordered_set<string> feature_tags = it->second;
		bool found_remove;
		if (r.empty())
			found_remove = false; // empty removed_tags mean do not remove any
		else {
			found_remove = false;
			auto start_it = feature_tags.begin();
			while (!found_remove && start_it != feature_tags.end()) {
				for (const string& substring : r) {
					boost::regex regi(substring);
					if (boost::regex_match(*start_it, regi)) {
						found_remove = true;
						if (verbose)
							MLOG("TagFeatureSelector removing [%s] because of tag [%s] that contains [%s]\n",
								feature_name.c_str(), (*start_it).c_str(), substring.c_str());
						break;
					}
				}
				++start_it;
			}
		}
		if (found_remove) {
			removed_features++;
			continue;
		}

		bool found_match;
		if (s.empty())
			found_match = true;// empty selected_tags mean select all
		else {
			found_match = false;
			auto start_it = feature_tags.begin();
			//MLOG("considering feature [%s]\n", feature_name.c_str());
			while (!found_match && start_it != feature_tags.end()) {
				//MLOG("considering tag [%s]\n",(*start_it).c_str());
				for (const string& substring : s) {
					boost::regex regi(substring);
					if (boost::regex_match(*start_it, regi)) {
						found_match = true;
						if (verbose)
							MLOG("TagFeatureSelector selecting [%s] because of tag [%s] that contains [%s]\n",
								feature_name.c_str(), (*start_it).c_str(), substring.c_str());
						break;
					}
				}
				++start_it;
			}
		}
		if (found_match) {
			selected_features++;
			selected.push_back(feature_name);
		}
	}
	MLOG("TagFeatureSelector selected_features %d removed_features %d\n", selected_features, removed_features);
	return 0;
}

//=======================================================================================
// Importance Feature Selection
//=======================================================================================
// Init 
//.......................................................................................
int ImportanceFeatureSelector::init(map<string, string>& mapper) {
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [ImportanceFeatureSelector::init]
		if (field == "missing_value") missing_value = med_stof(entry.second);
		else if (field == "predictor") predictor = entry.second;
		else if (field == "predictor_params") predictor_params = entry.second;
		else if (field == "importance_params") importance_params = entry.second;
		else if (field == "minStat") minStat = med_stof(entry.second);
		else if (field == "verbose") verbose = med_stoi(entry.second) > 0;
		else if (field == "numToSelect") numToSelect = med_stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknown parameter \'%s\' for ImportanceFeatureSelector\n", field.c_str());
		//! [ImportanceFeatureSelector::init]
	}

	return 0;
}

// Learn 
//.......................................................................................
int ImportanceFeatureSelector::_learn(MedFeatures& features, unordered_set<int>& ids) {
	MedPredictor *model = MedPredictor::make_predictor(predictor, predictor_params);
	vector<float> feat_importance;
	model->learn(features);

	model->calc_feature_importance(feat_importance, importance_params);

	vector<pair<string, float>> features_scores((int)feat_importance.size());
	map<string, vector<float>>::iterator it = features.data.begin();
	for (size_t i = 0; i < feat_importance.size(); ++i)
	{
		features_scores[i].first = it->first;
		features_scores[i].second = feat_importance[i];
		++it;
	}
	//sort features by scores:
	sort(features_scores.begin(), features_scores.end(), [](const pair<string, float> &v1, const pair<string, float> &v2)
	{return (v1.second > v2.second); });
	if (verbose) {
		it = features.data.begin();
		for (size_t i = 0; i < features_scores.size(); ++i)
		{
			MLOG("FEATURE %s : %2.3f\n", features_scores[i].first.c_str(), features_scores[i].second);
			++it;
		}
	}

	if (numToSelect == 0) {
		// Select according to minimum value of stat
		for (auto& rec : features_scores) {
			if (rec.second < minStat)
				break;
			selected.push_back(rec.first);
		}
	}
	else {
		// Select according to number
		int n = (features_scores.size() > numToSelect) ? numToSelect : (int)features_scores.size();
		selected.resize(n);
		for (int i = 0; i < n; i++) {
			MLOG("SELECTED FEATURE %s : %2.3f\n", features_scores[i].first.c_str(), features_scores[i].second);
			selected[i] = features_scores[i].first;
		}
	}

	return 0;
}

//=======================================================================================
// Iterative Feature Selection
//=======================================================================================
// Init 
//.......................................................................................
int IterativeFeatureSelector::init(map<string, string>& mapper) {
	init_defaults();

	string folds_s;
	for (auto entry : mapper) {
		string field = entry.first;
		//! [ImportanceFeatureSelector::init]
		if (field == "missing_value") missing_value = med_stof(entry.second);
		else if (field == "predictor") predictor = entry.second;
		else if (field == "predictor_params") predictor_params = entry.second;
		else if (field == "predictor_params_file") predictor_params_file = entry.second;
		else if (field == "nfolds") nfolds = stoi(entry.second);
		else if (field == "folds") folds_s = entry.second; 
		else if (field == "mode") mode = entry.second;
		else if (field == "rates") rates = entry.second;
		else if (field == "cohort_params") cohort_params = entry.second;
		else if (field == "msr_params") msr_params = entry.second;
		else if (field == "bootstrap_params") bootstrap_params = entry.second;
		else if (field == "verbose") verbose = med_stoi(entry.second) > 0;
		else if (field == "work_on_sets") work_on_sets = med_stoi(entry.second) > 0;
		else if (field == "numToSelect") numToSelect = stoi(entry.second);
		else if (field == "required") boost::split(required, entry.second, boost::is_any_of(","));
		else if (field == "ungrouped") boost::split(ungroupd_names, entry.second, boost::is_any_of(","));
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknown parameter \'%s\' for IterativeFeatureSelector\n", field.c_str());
		//! [ImportanceFeatureSelector::init]
	}

	// Parse folds
	vector<string> folds_v;
	boost::split(folds_v, folds_s, boost::is_any_of(","));
	for (string& s : folds_v)
		folds.push_back(stoi(s));

	// Parse rates
	get_rates_vec();

	// Read paramters
	read_params_vec();

	return 0;
}

// Learn 
//.......................................................................................
int IterativeFeatureSelector::_learn(MedFeatures& features, unordered_set<int>& ids) {

	int nSamples = (int)features.samples.size();

	// Divide features into families based on signals
	map<string, vector<string> > featureFamilies;
	get_features_families(features, featureFamilies, required, work_on_sets);

	// Collect original splits
	vector<int> orig_folds(nSamples);
	for (int i = 0; i < nSamples;i++)
		orig_folds[i] = features.samples[i].split;

	// Override splits
	map<int, int> id2fold;
	for (MedSample& sample : features.samples) {
		int id = sample.id;
		if (id2fold.find(id) == id2fold.end())
			id2fold[id] = globalRNG::rand() % nfolds;
		sample.split = id2fold[id];
	}

	// Boostrapping
	MedBootstrapResult bootstrapper;

	bootstrapper.bootstrap_params.loopCnt = 0;
	string bootstrap_init = bootstrap_params;
	boost::replace_all(bootstrap_init, "/", ";");
	boost::replace_all(bootstrap_init, ":", "=");
	bootstrapper.bootstrap_params.init_from_string(bootstrap_init);

	init_bootstrap_cohort(bootstrapper, cohort_params);
	init_bootstrap_params(bootstrapper, msr_params);

	// Optimize
	if (ids.empty()) {
		if (mode == "top2bottom")
			doTop2BottomSelection(features, featureFamilies, bootstrapper);
		else
			doBottom2TopSelection(features, featureFamilies, bootstrapper);
	}
	else {
		MedFeatures filteredFeatures = features;
		vector<int> indices;
		for (unsigned int i = 0; i < filteredFeatures.samples.size(); i++) {
			if (ids.find(filteredFeatures.samples[i].id) != ids.end())
				indices.push_back(i);
		}
		medial::process::filter_row_indexes(filteredFeatures, indices);

		if (mode == "top2bottom")
			doTop2BottomSelection(filteredFeatures, featureFamilies, bootstrapper);
		else
			doBottom2TopSelection(filteredFeatures, featureFamilies, bootstrapper);
	}

	// Reinstall splits
	for (int i = 0; i < nSamples;i++)
		features.samples[i].split = orig_folds[i];

}

// Report to file 
//.......................................................................................
void IterativeFeatureSelector::print_report(string& fileName) {
	
	ofstream of(fileName);
	if (!of)
		MTHROW_AND_ERR("Cannot open %s for writing\n", fileName.c_str());

	for (string& line : report)
		of << line << "\n";
}
