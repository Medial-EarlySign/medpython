#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "FeatureProcess.h"
#include "DoCalcFeatProcessor.h"
#include <omp.h>

//=======================================================================================
// FeatureSelector
//=======================================================================================
// Learn : Add required to feature selected by inheriting classes
//.......................................................................................
int FeatureSelector::Learn(MedFeatures& features, unordered_set<int>& ids) {

	// select, ignoring requirments
	if (_learn(features, ids) < 0)
		return -1;

	// Add required signals
	// Collect selected
	unordered_set<string> selectedFeatures;
	for (string& feature : selected)
		selectedFeatures.insert(feature);

	// Find Missing
	vector<string> missingRequired;
	for (string feature : required) {
		if (selectedFeatures.find(feature) == selectedFeatures.end())
			missingRequired.push_back(feature);
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
	for (string& feature : selected)
		MLOG("Feature Selection: Selected %s\n", feature.c_str());

	return 0;
}

// Apply selection : Ignore set of ids
//.......................................................................................
int FeatureSelector::Apply(MedFeatures& features, unordered_set<int>& ids) {

	unordered_set<string> selectedFeatures;
	for (string& feature : selected)
		selectedFeatures.insert(feature);

	return features.filter(selectedFeatures);
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

		if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "numToSelect") numToSelect = stoi(entry.second);
		else if (field == "method") params.method = params.get_method(entry.second);
		else if (field == "minStat") params.minStat = stof(entry.second);
		else if (field == "nBins") params.nBins = stoi(entry.second);
		else if (field == "binMethod") params.binMethod = params.get_binning_method(entry.second);
		else if (field == "required") boost::split(required, entry.second, boost::is_any_of(","));
		else if (field == "takeSquare") params.takeSquare = stoi(entry.second);
		else if (field == "max_samples") params.max_samples = stoi(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeatureSelector\n", field.c_str());
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
	for (int i = 0; i <nFeatures; i++) {
		int n;
		vector<float> values;
		get_all_values(features, names[i], ids, values, params.max_samples);

		stats[i] = fabs(get_pearson_corr(values, label, n, missing_value));
		if (n == 0) stats[i] = 0.0;

		// If required, test also correlation to squared (normalized) values
		if (params.takeSquare) {
			float mean, std;
			get_mean_and_std(values, mean, std, missing_value);
			for (unsigned int j = 0; j < values.size(); j++) {
				if (values[j] != missing_value)
					values[j] = (values[j] - mean)*(values[j] - mean);
			}

			float newStat = fabs(get_pearson_corr(values, label, n, missing_value));
			if (newStat > stats[i])
				stats[i] = newStat;
		}
	}

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
		get_mutual_information(binnedValues[i], binnedLabel, n, stats[i]);
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
	get_dMatrix(label, labelDistances, missing_value);
	float targetDistVar = get_dVar(labelDistances);
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
		get_dMatrix(values, valueDistances, missing_value);
		float valueDistVar = get_dVar(valueDistances);
		float distCov = get_dCov(labelDistances, valueDistances);
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

					score -= penalty*penaltyValue;
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
			stats(i, index) = fabs(get_pearson_corr(values[i], target, n, missing_value));
			if (n == 0) stats(i, index) = 0.0;
			stats(index, i) = stats(i, index);
		}
	}
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
			get_mutual_information(binnedValues[i], binnedTarget, n, stats(i, index));
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
	get_dMatrix(target, targetDistances, missing_value);
	float targetDistVar = get_dVar(targetDistances);
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
			get_dMatrix(values, valueDistances, missing_value);
			float valueDistVar = get_dVar(valueDistances);
			float distCov = get_dCov(targetDistances, valueDistances);
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
	vector<double> lambdas(nthreads);
	float minLambda = 0.0, maxLambda = initMaxLambda;
//	vector<MedLasso> predictors(nthreads);
	vector<MedGDLM> predictors(nthreads);
	vector<int> nSelected(nthreads);
	vector<float> w(x.nrows, 1.0);

	MLOG("Lasso Feature Selection : From %d to %d\n", nFeatures, numToSelect);
	selected.clear();
	float lowerBoundLambda = 0.0, upperBoundLambda = -1.0;

	// Prevent being stuck ...
	int nStuck = 0; 
	int prevMaxN = -1, prevMinN = -1;

	// Search
	while (!found) {
		if (nthreads == 1)
			lambdas[0] = maxLambda;
		else {
			float step = (maxLambda - minLambda) / (nthreads - 1);
			for (int i = 0; i < nthreads; i++)
				lambdas[i] = minLambda + step *i;
		}
	
		for (int i = 0; i < nthreads; i++) {
			predictors[i].init_from_string("method = logistic_sgd; last_is_bias = 0; stop_at_err = 1e-4; batch_size = 2048; momentum = 0.95; rate = 0.01; rate_decay = 1; l_ridge = 0; l_lasso = " + to_string(lambdas[i]) + "; err_freq = 10; nthreads = 12");
			predictors[i].params.l_lasso = (float)lambdas[i];
//			predictors[i].params.lambda = lambdas[i];
//			predictors[i].initialize_vars(x.data_ptr(), y.data_ptr(), &(w[0]), predictors[i].b, x.nrows, x.ncols);
		}

#pragma omp parallel for 
		for (int i = 0; i < nthreads; i++) {
//			predictors[i].lasso_regression(predictors[i].b, x.nrows, x.ncols, lambdas[i], predictors[i].params.num_iterations);
			predictors[i].learn(x, y);

			// Identify non-zero parameters
			nSelected[i] = 0;
			for (int j = 0; j < nFeatures; j++) {
				if (predictors[i].b[j] != 0)
					nSelected[i] ++;
			}
		}

		MLOG_V("Lasso Feature Selection: [%f,%f] : nFeatures [%d,%d] nStuck %d\n", lambdas[0], lambdas[nthreads - 1], nSelected[0], nSelected[nthreads - 1],nStuck);

		if (nthreads == 1) { // Special care
			if (nSelected[0] == numToSelect) {
				found = 1;
				for (int j = 0; j < nFeatures; j++) {
					if (predictors[0].b[j] != 0)
						selected.push_back(names[j]);
				}
			} 
			else if (nSelected[0] > numToSelect) {
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
				MLOG("N[%.12f] = %d\n", lambdas[j], nSelected[j]);

//			float ratio;
			if (nSelected[nthreads - 1] > numToSelect) { // MaxLambda is still too low
				minLambda = maxLambda;
				maxLambda *= 2.0;
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
						minLambda = (float)lambdas[i - 1];
						maxLambda = (float)lambdas[i];
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
				MLOG("Stuck at same N range for 3 steps. That's enough for now ... Actual NSelected = %d\n",nSelected[optimalI]);
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

	map<float, string> b2Name;
	int i = 0;
	for (auto &feat : features.data) {
		b2Name[predictors[0].b[i++]] = feat.first;
	}

	for (auto &b : b2Name) {
		if (b.first != 0) {
			MLOG("Feature %s :: B %f\n", b.second.c_str(), b.first);
		}
	}


	return 0;
}

// Init
//.......................................................................................
int LassoSelector::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "numToSelect") numToSelect = stoi(entry.second);
		else if (field == "initMaxLambda") initMaxLambda = stof(entry.second);
		else if (field == "nthreads") nthreads = stoi(entry.second);
		else if (field == "required") boost::split(required, entry.second, boost::is_any_of(","));
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeatureSelector\n", field.c_str());
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
			MLOG("Removing %s : %f of values is %f\n", name.c_str(), p, maxCountValue);
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

		if (field == "missing_value") missing_value = stof(entry.second);
		if (field == "percentage") percentage = stof(entry.second);
		else if (field != "names" && field != "fp_type" && field != "tag")
			MLOG("Unknonw parameter \'%s\' for FeatureSelector\n", field.c_str());
	}

	assert(percentage >= 0 && percentage <= 1.0);
	return 0;

}
