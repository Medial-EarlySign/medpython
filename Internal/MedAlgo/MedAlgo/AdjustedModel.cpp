//
// AdjustedModel : "adjusted" , learn a model on (preferably) matched set, and the adjust according to given priors
//
// input example: 
// "predictor=gdlm;predictor_params={method=logistic_sgd;last_is_bias=0;stop_at_err=1e-4;batch_size=2048;momentum=0.95;rate=0.1};priors=priorsFileName;onlyPriorsFeatures=Race,Education;odds=0.02
//
#include <MedAlgo/MedAlgo/AdjustedModel.h>
#include <iostream>
#include <string>

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL LOG_DEF_LEVEL
extern MedLogger global_logger;


//==================================================================
int MedAdjustedModel::set_params(map<string, string>& mapper)
{
	//! [MedAdjustedModel::init]
	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "predictor_params")
			predictor_params = entry.second;
		else if (field == "predictor")
			predictor_name = entry.second;
		else if (field == "priors")
			priorsFile = entry.second;
		else if (field == "odds")
			odds = stof(entry.second);
		else if (field == "onlyPriorsFeatures")
			boost::split(onlyPriorsFeatures, entry.second, boost::is_any_of(","));
		//! [MedAdjustedModel::init]

	}

	// initializing predictor
	predictor = MedPredictor::make_predictor(predictor_name, predictor_params);
	if (predictor == NULL)
		MTHROW_AND_ERR("ERROR: MedAdjustedModel : failed initializing predictor: %s with params: %s\n", predictor_name.c_str(), predictor_params.c_str());

	return 0;
}

//==================================================================
void MedAdjustedModel::readPriors() {

	// Read
	ifstream inf(priorsFile);
	if (!inf)
		MTHROW_AND_ERR("Cannot open %s for reading\n", priorsFile.c_str());

	string line;
	vector<string> fields;
	vector<vector<int>> values;
	vector<float> probs;
	int nValues = 0;
	while (getline(inf, line)) {
		boost::split(fields, line, boost::is_any_of("\t "));
		if (nValues == 0) {
			if (fields.back() != "prob")
				MTHROW_AND_ERR("Expecting last column of header to be prob");
			for (size_t i = 0; i < fields.size() - 1; i++)
				priors.names.push_back(fields[i]);
			nValues = priors.names.size();
		}
		else {
			if (fields.size() != priors.names.size() + 1)
				MTHROW_AND_ERR("Cannot parse prior line %s in %s\n", line.c_str(), priorsFile.c_str());
			vector<int> newValues(fields.size() - 1);
			for (int i = 0; i < nValues; i++)
				newValues[i] = stoi(fields[i]);
			values.push_back(newValues);
			probs.push_back(stof(fields.back()));
		}
	}
	if (probs.size() == 0)
		MTHROW_AND_ERR("No lines in priors file ?\n");

	// Prepare
	priors.min.resize(nValues);
	priors.max.resize(nValues);
	for (int i = 0; i < nValues; i++) {
		priors.min[i] = values[0][i];
		priors.max[i] = values[0][i];
		for (size_t j = 0; j < values.size(); j++) {
			if (values[j][i] > priors.max[i])
				priors.max[i] = values[j][i];
			if (values[j][i] < priors.min[i])
				priors.min[i] = values[j][i];
		}
		MLOG("Value %s range = %d - %d\n", priors.names[i].c_str(), priors.min[i], priors.max[i]);
	}

	priors.factors.assign(nValues, 1);
	int totNum = (priors.max[0] - priors.min[0] + 1);
	for (int i = 1; i < nValues; i++) {
		priors.factors[i] = priors.factors[i - 1] * (priors.max[i - 1] - priors.min[i - 1] + 1);
		totNum *= (priors.max[i] - priors.min[i] + 1);
	}

	// Fill
	priors.probs.assign(totNum, -1.0);
	for (size_t i = 0; i < values.size(); i++) {
		int index = 0;
		for (int j = 0; j < nValues; j++)
			index += (values[i][j] - priors.min[j])*priors.factors[j];
		priors.probs[index] = probs[i];
	}

	// Check
	for (int i = 0; i < totNum; i++) {
		if (priors.probs[i] == -1.0)
			MTHROW_AND_ERR("Missing entry with index=%d\n", i);
	}
}

//==================================================================
void MedAdjustedModel::getPosteriors(vector<float>& preds, MedMat<float>& priorsX) const 
{

#pragma omp parallel for
	for (size_t i = 0; i < preds.size(); i+=n_preds) {
		// Prior
		int index = 0;
		for (size_t j = 0; j < priors.names.size(); j++) {
			int value = (int)priorsX(i, j);
			if (value < priors.min[j] || value > priors.max[j])
				MTHROW_AND_ERR("Value %d of %s is outside priors range [%d,%d]\n", value, priors.names[j].c_str(), priors.min[j], priors.max[j]);
			index += (value - priors.min[j])*priors.factors[j];
		}
		float prior = priors.probs[index];

		float pred = preds[i];
		for (int iPred = 0; iPred < n_preds; iPred++)
			preds[i + iPred] = (pred*prior) / (pred*prior + (1.0 - pred)*(1.0 - prior)*odds);
	}
}

//==================================================================
void MedAdjustedModel::getOdds(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts) {

	double sumProb = 0, sum = 0;
	for (int i = 0; i < x.nrows; i++) {
		int index = 0;
		for (size_t j = 0; j < priors.names.size(); j++) {
			int value = (int)x(i, j);
			if (value < priors.min[j] || value > priors.max[j])
				MTHROW_AND_ERR("Value %d of %s is outside priors range [%d,%d]\n", value, priors.names[j].c_str(), priors.min[j], priors.max[j]);
			index += (value - priors.min[j])*priors.factors[j];
		}
		sumProb += priors.probs[index] * wgts[i];
		sum += wgts[i];
	}

	odds =  sumProb / sum;
}

//==================================================================
void MedAdjustedModel::getSubMatrices(MedMat<float> &x, MedMat<float> &priorsX, MedMat<float> &predictorX) const 
{
	unordered_set<int> priorColiIds;
	for (int id : priors.colIds)
		priorColiIds.insert(id);

	int nPriorCols = priorColiIds.size();
	int nPredictorCols = x.ncols - nPriorCols;
	priorsX.resize(x.nrows, nPriorCols);
	predictorX.resize(x.nrows, nPredictorCols);

#pragma parallel for
	int priorsCol = 0, predictorCol = 0;
	for (int icol = 0; icol < x.ncols; icol++) {
		if (priorColiIds.find(icol) != priorColiIds.end()) {
			for (int irow = 0; irow < x.nrows; irow++)
				priorsX(irow, priorsCol) = x(irow, icol);
			priorsCol++;
		}
		else {
			for (int irow = 0; irow < x.nrows; irow++)
				predictorX(irow, predictorCol) = x(irow, icol);
			predictorCol++;
		}
	}
}
//====================================================================================================
int MedAdjustedModel::learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts)
{

	// Read Priors from file
	readPriors();

	// Resolve feature names
	for (size_t i = 0; i < onlyPriorsFeatures.size(); i++)
		resolvedOnlyPriorsFeatures.insert(x.signals[find_in_feature_names(x.signals, onlyPriorsFeatures[i])]);

	priors.colIds.resize(priors.names.size());
	for (size_t i = 0; i < priors.names.size(); i++)
		priors.colIds[i] = find_in_feature_names(x.signals, priors.names[i]);

	// Get Submatrices
	MedMat<float> priorsX, predictorX;
	getSubMatrices(x, priorsX, predictorX);

	// Learn odds
	if (odds < 0)
		getOdds(priorsX,y, wgts);

	// Learn internal predictor
	predictor->learn(predictorX,y,wgts);

	return 0;
}

//====================================================================================================
int MedAdjustedModel::predict(MedMat<float> &x, vector<float> &preds) const
{
	// Get Submatrices
	MedMat<float> priorsX, predictorX;
	getSubMatrices(x, priorsX, predictorX);

	// Learn internal predictor
	predictor->predict(predictorX,preds);

	// Performa adjustment
	getPosteriors(preds, priorsX);

	return 0;
}