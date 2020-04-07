#include "ProbAdjustPostProcessor.h"
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL LOG_DEF_LEVEL

//==================================================================
int ProbAdjustPostProcessor::init(map<string, string>& mapper)
{
	//! [MedAdjustedModel::init]
	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "priors")
			priorsFile = entry.second;
		else if (field == "json")
			priorsJson = entry.second;
		else if (field == "use_p")
			use_p = stof(entry.second);
		else if (field == "odds") {
			vector<string> odds_s;
			boost::split(odds_s, entry.second, boost::is_any_of(","));
			for (string& _odds : odds_s)
				odds.push_back(stof(_odds));
		}
		//! [MedAdjustedModel::init]

	}

	return 0;
}

//==================================================================
void ProbAdjustPostProcessor::readPriors() {

	// Read
	ifstream inf(priorsFile);
	if (!inf)
		MTHROW_AND_ERR("Cannot open %s for reading\n", priorsFile.c_str());

	string line;
	vector<string> fields;
	vector<vector<int>> values;
	vector<float> _probs;
	int nValues = 0;
	while (getline(inf, line)) {
		boost::split(fields, line, boost::is_any_of("\t "));
		if (nValues == 0) {
			if (fields.back() != "prob")
				MTHROW_AND_ERR("Expecting last column of header to be prob");
			for (size_t i = 0; i < fields.size() - 1; i++)
				names.push_back(fields[i]);
			nValues = names.size();
		}
		else {
			if (fields.size() != names.size() + 1)
				MTHROW_AND_ERR("Cannot parse prior line %s in %s\n", line.c_str(), priorsFile.c_str());
			vector<int> newValues(fields.size() - 1);
			for (int i = 0; i < nValues; i++)
				newValues[i] = stoi(fields[i]);
			values.push_back(newValues);
			_probs.push_back(stof(fields.back()));
		}
	}
	if (_probs.size() == 0)
		MTHROW_AND_ERR("No lines in priors file ?\n");

	// Prepare
	min.resize(nValues);
	max.resize(nValues);
	for (int i = 0; i < nValues; i++) {
		min[i] = values[0][i];
		max[i] = values[0][i];
		for (size_t j = 0; j < values.size(); j++) {
			if (values[j][i] > max[i])
				max[i] = values[j][i];
			if (values[j][i] < min[i])
				min[i] = values[j][i];
		}
		MLOG("Value %s range = %d - %d\n", names[i].c_str(), min[i], max[i]);
	}

	factors.assign(nValues, 1);
	int totNum = (max[0] - min[0] + 1);
	for (int i = 1; i < nValues; i++) {
		factors[i] = factors[i - 1] * (max[i - 1] - min[i - 1] + 1);
		totNum *= (max[i] - min[i] + 1);
	}

	// Fill
	probs.assign(totNum, -1.0);
	for (size_t i = 0; i < values.size(); i++) {
		int index = 0;
		for (int j = 0; j < nValues; j++)
			index += (values[i][j] - min[j])*factors[j];
		probs[index] = _probs[i];
	}

	// Check
	for (int i = 0; i < totNum; i++) {
		if (probs[i] == -1.0)
			MTHROW_AND_ERR("Missing entry with index=%d\n", i);
	}
}

//==================================================================
void ProbAdjustPostProcessor::getOdds(const MedFeatures &matrix) {

	if (matrix.samples.size() == 0)
		MTHROW_AND_ERR("No Matrix given. Cannot learn odds\n");
	if (matrix.samples[0].prediction.size() == 0)
		MTHROW_AND_ERR("No Predictions given. Cannot learn odds\n");

	int nPreds = (int)matrix.samples[0].prediction.size();
	vector<float> sumProbs(nPreds), sums(nPreds);
	for (size_t i = 0; i < matrix.samples.size(); i++) {
		for (int j = 0; j < nPreds; j++) {
			float w = (matrix.weights.empty()) ? 1.0 : matrix.weights[i];
			sums[j] += w;
			sumProbs[j] += w * matrix.samples[i].prediction[j];
		}
	}

	odds.resize(nPreds);
	for (int j = 0; j < nPreds; j++) {
		float meanProb = sumProbs[j] / sums[j];
		odds[j] = meanProb / (1 - meanProb);
	}

}

//====================================================================================================
void ProbAdjustPostProcessor::Learn(const MedFeatures &matrix)
{

	// Create model for generating features
	priorsModel = new MedModel;
	priorsModel->init_from_json_file(priorsJson);
	MedSamples _samples;
	_samples.import_from_sample_vec(matrix.samples);
	priorsModel->learn(*p_rep, _samples, MED_MDL_LEARN_REP_PROCESSORS, MED_MDL_APPLY_FTR_PROCESSORS);

	// Read Priors from file
	readPriors();

	// Resolve feature names
	resolvedNames.resize(names.size());
	for (size_t i = 0; i < names.size(); i++)
		resolvedNames[i] = priorsModel->features.resolve_name(names[i]);

	// Learn odds
	if (odds.empty())
		getOdds(priorsModel->features);

}

//====================================================================================================
void ProbAdjustPostProcessor::Apply(MedFeatures &matrix) const {

	MedSamples _samples;
	_samples.import_from_sample_vec(matrix.samples);
	priorsModel->apply(*p_rep, _samples, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_APPLY_FTR_PROCESSORS);

	vector<const float *> data_p(resolvedNames.size());
	for (size_t j = 0; j < data_p.size(); j++)
		data_p[j] = priorsModel->features.data.at(resolvedNames[j]).data();

#pragma omp parallel for
	for (size_t i = 0; i < priorsModel->features.samples.size(); i++) {
		// Prior
		int index = 0;
		for (size_t j = 0; j < resolvedNames.size(); j++) {
			int value = (int)priorsModel->features.data.at(resolvedNames[j])[i];
			if (value < min[j] || value > max[j])
				MTHROW_AND_ERR("Value %d of %s is outside priors range [%d,%d]\n", value, resolvedNames[j].c_str(), min[j], max[j]);
			index += (value - min[j])*factors[j];
		}
		float prior = probs[index];

		for (size_t j = 0; j < odds.size(); j++)
			matrix.samples[i].prediction[j] = (matrix.samples[i].prediction[j] * prior) / (matrix.samples[i].prediction[j] * prior + (1.0 - matrix.samples[i].prediction[j])*(1.0 - prior)*odds[j]);
	}
}



