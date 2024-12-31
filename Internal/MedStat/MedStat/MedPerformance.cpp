#include "MedPerformance.h"

#include <MedUtils/MedUtils/MedGlobalRNG.h>
#include <MedUtils/MedUtils/MedGenUtils.h>
#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//.........................................................................................................................................
// Class : Measurement
// Constructors

// Default is full AUC
Measurement::Measurement() {
	setParam = "NONE";
	setValue = 1.0;
	queriedParam = "AUC";
}

Measurement::Measurement(const string& qParam, const string& sParam, float sValue) {
	setParam = sParam;
	setValue = sValue;
	queriedParam = qParam;
}

Measurement::Measurement(const string& qParam, float sValue) {
	setParam = "NONE";
	setValue = sValue;
	queriedParam = qParam; ;
}

Measurement::Measurement(const string& qParam) {
	queriedParam = qParam;
	if (queriedParam == "AUC")
		setValue = 1.0;
	else
		setValue = -1.0;
	setParam = "NONE";
}

//.........................................................................................................................................
// Get name of measurement
string Measurement::name() {

	if (queriedParam == "AUC") {
		if (setValue == 1.0)
			return "AUC";
		else
			return "Partial_AUC_" + to_string(setValue);
	}
	else if (setParam == "NONE")
		return queriedParam + "_" + to_string(setValue);
	else
		return queriedParam + "@" + setParam + to_string(setValue);
}

//.........................................................................................................................................
// Init from name
void Measurement::get_from_name(string& name, Measurement& msr) {

	vector<string> fields;
	boost::split(fields, name, boost::is_any_of(" \t"));
	if (fields.size() == 1) {
		string qParam = fields[0];
		Measurement newMsr(qParam);
		msr = newMsr;
	}
	else if (fields.size() == 2) {
		string qParam = fields[0];
		float sValue = stof(fields[1]);
		Measurement newMsr(qParam, sValue);
		msr = newMsr;
	}
	else if (fields.size() == 3) {
		string qParam = fields[0];
		string sParam = fields[1];
		float sValue = stof(fields[2]);
		Measurement newMsr(qParam, sParam, sValue);
		msr = newMsr;
	}
	else
		MTHROW_AND_ERR("Cannot parse measurement entry \'%s\'\n", name.c_str());
}

//.........................................................................................................................................
// Read from file
void Measurement::read_from_file(string& fileName, vector<Measurement>& msrs) {

	ifstream inf(fileName);
	if (!inf)
		MTHROW_AND_ERR("Cannot open %s for reading\n", fileName.c_str());

	string line;
	while (getline(inf, line)) {
		Measurement tmpMsr;
		Measurement::get_from_name(line, tmpMsr);
		msrs.push_back(tmpMsr);
	}

}

//.........................................................................................................................................
// Initialize
void MedClassifierPerformance::init() {

	nneg.clear();
	npos.clear();
	tps.clear();
	fps.clear();
	MeasurementValues.clear();
}

//.........................................................................................................................................
// Constructors
// Default
MedClassifierPerformance::MedClassifierPerformance() {
	preds.clear();
	init();
}

//.........................................................................................................................................
// Loaders

// Load from Smplaes
void MedClassifierPerformance::_load(MedSamples& samples) {

	// First pass to check for splits
	set<int> splits;
	int missingSplit = -1;
	for (auto& idSamples : samples.idSamples) {
		int split = idSamples.split;
		if (split == -1)
			missingSplit = 1;
		else
			splits.insert(split);
	}

	if (missingSplit != -1 && (!(splits.empty())))
		MTHROW_AND_ERR("Cannot have both split and non-split samples\n");

	// Fill
	int nSplits = (int)splits.size();
	map<int, int> split_idx;
	int idx = 0;
	for (int split : splits)
		split_idx[split] = ++idx;
	split_idx[-1] = 0;

	preds.resize(nSplits+1);
	for (auto& idSamples : samples.idSamples) {
		int split = idSamples.split;
		for (auto& sample : idSamples.samples)
			preds[split_idx[split]].push_back({ sample.prediction.back(),sample.outcome });
	}

	// Are there any splits ?
	if (nSplits>0)
		SplitsToComplete();
}

// Load from vectors
void MedClassifierPerformance::_load(vector<pair<float, float> >& in_preds) {
	preds.push_back(in_preds);
}

void MedClassifierPerformance::_load(vector<vector<pair<float, float> > >& in_split_preds) {

	preds.resize(in_split_preds.size() + 1);
	for (unsigned int i = 0; i < in_split_preds.size(); i++)
		preds[i + 1] = in_split_preds[i];

	SplitsToComplete();

}

// Load from predictor data
void MedClassifierPerformance::_load(MedFeatures& ftrs) {

	// First pass to check for splits
	int maxSplits = -1, missingSplit = -1;
	for (auto& sample: ftrs.samples) {
		if (sample.split > maxSplits)
			maxSplits = sample.split;
		if (sample.split == -1)
			missingSplit = 1;
	}
	assert(missingSplit == -1 || maxSplits == -1); // Can't have both splits and non-splits !

	// Fill
	preds.resize(maxSplits + 2);
	for (auto& sample : ftrs.samples) {
		int split = sample.split;
		preds[split + 1].push_back({ sample.prediction.back(),sample.outcome });
	}

	// Are there any splits ?
	if (maxSplits > -1)
		SplitsToComplete();
}

void MedClassifierPerformance::post_load() {
	init();
	ShuffleSort();
	getPerformanceValues();

	PerformancePointers.resize(preds.size());
}

//.........................................................................................................................................
// Helpers

// From splits to complete
void MedClassifierPerformance::SplitsToComplete() {
	preds[0].clear();
	for (int i = 1; i < preds.size(); i++)
		for (int j = 0; j < preds[i].size(); j++)
			preds[0].push_back(preds[i][j]);

}

// Shuffle + Sort
void MedClassifierPerformance::ShuffleSort() {

	for (unsigned int i = 0; i < preds.size(); i++) {
		shuffle(preds[i].begin(), preds[i].end(), globalRNG::get_engine());
		sort(preds[i].begin(), preds[i].end(), _PredsCompare());
	}
}

void MedClassifierPerformance::Count() {

	npos.resize(preds.size(), 0);
	nneg.resize(preds.size(), 0);
	tps.resize(preds.size());
	fps.resize(preds.size());

	for (unsigned int i = 0; i < preds.size(); i++) {
		tps[i].resize(preds[i].size());
		fps[i].resize(preds[i].size());

		for (unsigned j = 0; j < preds[i].size(); j++) {
			if (preds[i][j].second <= 0)
				nneg[i] ++;
			else
				npos[i] ++;

			tps[i][j] = npos[i];
			fps[i][j] = nneg[i];
		}
	}

}

void MedClassifierPerformance::getPerformanceValues() {

	Count();
	PerformanceValues.resize(preds.size());

	for (unsigned int i = 0; i < preds.size(); i++) {
		PerformanceValues[i]["Score"].resize(preds[i].size(), -1);
		PerformanceValues[i]["Sens"].resize(preds[i].size(), -1);
		PerformanceValues[i]["Spec"].resize(preds[i].size(), -1);
		PerformanceValues[i]["PPV"].resize(preds[i].size(), -1);
		PerformanceValues[i]["NPV"].resize(preds[i].size(), -1);
		PerformanceValues[i]["OR"].resize(preds[i].size(), -1);

		for (unsigned j = 0; j < preds[i].size(); j++) {
			// Score
			PerformanceValues[i]["Score"][j] = preds[i][j].first;
			// Sens
			if (npos[i] > 0)
				PerformanceValues[i]["Sens"][j] = ((float)tps[i][j]) / npos[i];
			// Spec
			if (nneg[i] > 0)
				PerformanceValues[i]["Spec"][j] = (float) 1.0 - ((float)fps[i][j]) / nneg[i];
			// PPV
			if (tps[i][j] + fps[i][j] > 0)
				PerformanceValues[i]["PPV"][j] = ((float)tps[i][j]) / (tps[i][j] + fps[i][j]);
			// NPV
			if (nneg[i] - fps[i][j] + npos[i] - tps[i][j] > 0)
				PerformanceValues[i]["NPV"][j] = ((float)nneg[i] - fps[i][j]) / (nneg[i] - fps[i][j] + npos[i] - tps[i][j]);
			// OR
			if (fps[i][j] > 0 && nneg[i] - fps[i][j] > 0 && (npos[i] - tps[i][j]) / ((float)nneg[i] - fps[i][j]) > 0)
				PerformanceValues[i]["OR"][j] = (((float)tps[i][j]) / fps[i][j]) / (((float)npos[i] - tps[i][j]) / (nneg[i] - fps[i][j]));
		}
	}
}

// Queries
// Parameter at point determined by another parameters (e.g. PPV at Specificity = 0.99 is GetPerformanceParam("PPV","Spec",0.99,outPPV). setParams = (Score,Sens,Spec), queriedParams = (Score,Sens,Spec,PPV,NPV,OR)
int MedClassifierPerformance::GetPerformanceParam(const string& setParam, const string& queriedParam, float setValue) {
	
	pair<string, float> set(setParam, setValue);
	Measurement inMeasurement(queriedParam, setParam, setValue);
	if (MeasurementValues.find(inMeasurement) != MeasurementValues.end())
		return 0;

	MeasurementValues[inMeasurement].resize(preds.size(), MED_MAT_MISSING_VALUE);

	for (unsigned int i = 0; i < preds.size(); i++) {
		if (getPerformanceValues(set, queriedParam, i, MeasurementValues[inMeasurement]) < 0)
			return -1;
	}

	return 0;
}

// General performance parameter, with optional value (e.g. AUC = GetPerformanceParam("AUC",outAuc) or GetPerformanceParam("AUC",1.0,outAUC). Partial AUC = GetPerformanceParam("AUC",0.2,partAUC)
int MedClassifierPerformance::GetPerformanceParam(const string& qParam, float sValue) {

	Measurement inMeasurement(qParam, sValue);

	if (MeasurementValues.find(inMeasurement) != MeasurementValues.end())
		return 0;

	MeasurementValues[inMeasurement].resize(preds.size(), MED_MAT_MISSING_VALUE);
	if (qParam == "AUC") {
		getAUC(sValue, MeasurementValues[inMeasurement]);
		return 0;
	}
	else {
		fprintf(stderr, "Unknown required performance parameter %s\n", qParam.c_str());
		return -1;
	}
}

int MedClassifierPerformance::GetPerformanceParam(const string& qParam) {

	if (qParam == "AUC") {
		return GetPerformanceParam(qParam, 1.0);
	}
	else {
		fprintf(stderr, "Unknown required performance parameter %s\n", qParam.c_str());
		return -1;
	}
}

int MedClassifierPerformance::GetPerformanceParam(Measurement& inMeasurement) {

	if (inMeasurement.setParam == "NONE")
		return GetPerformanceParam(inMeasurement.queriedParam, inMeasurement.setValue);
	else
		return GetPerformanceParam(inMeasurement.setParam, inMeasurement.queriedParam, inMeasurement.setValue);
}

vector<float> MedClassifierPerformance::operator() (Measurement& inMeasurement) {
	vector<float> outValues;
	if (GetPerformanceParam(inMeasurement) != -1)
		outValues = MeasurementValues[inMeasurement];

	return outValues;
}

// Helpers for queries
// Get Pointers ...
int MedClassifierPerformance::getPerformancePointer(pair<string, float>& set, int index) {

	if (set.first == "Sens")
		return getPointer(set.first, set.second, index, 1);
	else if (set.first == "Spec")
		return getPointer(set.first, set.second, index, -1);
	else if (set.first == "Score")
		return getPointer(set.first, set.second, index, -1);
	else {
		fprintf(stderr, "Unknown set parameters %s\n", set.first.c_str());
		return -1;
	}
}

int MedClassifierPerformance::getPointer(const string& param, float value, int index, int direction) {

	pair<string, float> pointer(param, value);

	// Find Value >= or <= TargetValue
	int targetIdx = -1;
	for (unsigned int i = 0; i < preds[index].size(); i++) {
		if ((direction == 1 && PerformanceValues[index][param][i] >= value) || (direction == -1 && PerformanceValues[index][param][i] <= value)) {
			targetIdx = i;
			break;
		}
	}

	// Is target outside the range ?
	if (targetIdx == -1 || (targetIdx == 0 && PerformanceValues[index][param][targetIdx] != value))
		return -1;

	if (PerformanceValues[index][param][targetIdx] == value) {
		// Are we exactly at target ?
		PerformancePointers[index][pointer].first = targetIdx;
		while (targetIdx < preds[index].size() && PerformanceValues[index][param][targetIdx] == value)
			targetIdx++;
		PerformancePointers[index][pointer].second = targetIdx - 1;
	}
	else {
		// Have we passed the target ?
		PerformancePointers[index][pointer] = pair<int, int>(targetIdx - 1, targetIdx + 1);
	}

	return 0;
}

// Get Values
int MedClassifierPerformance::getPerformanceValues(pair<string, float>& set, const string &queriedParam, int index, vector<float>& queriedValues) {

	if (PerformanceValues[index].find(queriedParam) == PerformanceValues[index].end()) {
		fprintf(stderr, "Cannot query parameter %s\n", queriedParam.c_str());
		return -1;
	}

	if (PerformancePointers[index].find(set) == PerformancePointers[index].end()) {
		if (getPerformancePointer(set, index) < 0) {
			fprintf(stderr, "Cannot find pointer for %s:%f for index = %d\n", set.first.c_str(), set.second, index);
			return -1;
		}
	}

	int start = PerformancePointers[index][set].first;
	int end = PerformancePointers[index][set].second;

	if (PerformanceValues[index][set.first][start] != set.second) {
		// Are we in-between ?
		float d1 = fabs(PerformanceValues[index][set.first][start] - set.second);
		float v1 = PerformanceValues[index][queriedParam][start];

		float d2 = fabs(PerformanceValues[index][set.first][end] - set.second);
		float v2 = PerformanceValues[index][queriedParam][end];
		queriedValues[index] = (d1*v2 + d2*v1) / (d2 + d1);
	}
	else {
		// Are we exactly at the value ?
		float sum = 0;
		for (int i = start; i <= end; i++)
			sum += PerformanceValues[index][queriedParam][i];
		queriedValues[index] = sum / (end - start + 1);
	}

	return 0;
}

// Get AUC
void MedClassifierPerformance::getAUC(float maxFPR, vector<float>& qValues) {

	qValues.resize(preds.size());
	for (unsigned int i = 0; i < preds.size(); i++)
		qValues[i] = getAUC(maxFPR, i);
}

float MedClassifierPerformance::getAUC(float maxFPR, int idx) {

	double auc = 0;
	int targetFP = (int)(maxFPR * nneg[idx] + 0.5);
	for (unsigned int i = 0; i < preds[idx].size(); i++) {
		if (preds[idx][i].second <= 0) {
			auc += tps[idx][i];

			if (fps[idx][i] == targetFP)
				break;
		}
	}

	return (float)(auc / ((double)npos[idx] * (double)targetFP));
}

// Performance Graph
int MedClassifierPerformance::GetPrformanceGraph(const string& xParam, const string& yParam, vector<vector<float> >& x, vector<vector<float> >& y) {

	x.resize(preds.size());
	y.resize(preds.size());

	for (unsigned int i = 0; i < preds.size(); i++) {
		if (PerformanceValues[i].find(xParam) != PerformanceValues[i].end())
			x[i] = PerformanceValues[i][xParam];
		else if (xParam == "TPR") {
			x[i].resize(preds[i].size());
			for (unsigned j = 0; j < preds[i].size(); j++)
				x[i][j] = (float) 1.0 - PerformanceValues[i]["Spec"][j];
		}
		else
			fprintf(stderr, "Unknown parameters %s\n", xParam.c_str());

		if (PerformanceValues[i].find(yParam) != PerformanceValues[i].end())
			y[i] = PerformanceValues[i][yParam];
		else if (yParam == "TPR") {
			y[i].resize(preds[i].size());
			for (unsigned j = 0; j < preds[i].size(); j++)
				y[i][j] = (float) 1.0 - PerformanceValues[i]["Spec"][j];
		}
		else
			fprintf(stderr, "Unknown parameters %s\n", xParam.c_str());
	}

	return 0;
}

// Comparison (1 - "Current is Better", 0 - "Current is NOT Better", -1 - Problem)
int MedClassifierPerformance::compare(MedClassifierPerformance& other) {

	if (MeasurementValues.find(comparePoint) == MeasurementValues.end())
		GetPerformanceParam(comparePoint);

	vector<float>& values = MeasurementValues[comparePoint];

	if (other.MeasurementValues.find(comparePoint) == other.MeasurementValues.end())
		other.GetPerformanceParam(comparePoint);

	vector<float>& otherValues = other.MeasurementValues[comparePoint];

	if (compareMode == PRF_COMPARE_FULL)
		return (values[0] > otherValues[0]) ? 1 : 0;
	else if (compareMode == PRF_COMPARE_ALL) {
		for (unsigned int i = 0; i < values.size(); i++) {
			if (otherValues[i] >= values[i])
				return 0;
		}
		return 1;
	}
	else if (compareMode == PRF_COMPARE_SPLITS) {
		if (values.size() == 1) {
			fprintf(stderr, "No splits to compare\n");
			return -1;
		}
		for (unsigned int i = 1; i < values.size(); i++) {
			if (otherValues[i] >= values[i])
				return 0;
		}
		return 1;
	}
	else if (compareMode == PRF_COMPARE_FULL_AND_PART_SPLITS) {
		if (values.size() == 1) {
			fprintf(stderr, "No splits to compare\n");
			return -1;
		}
		if (otherValues[0] >= values[0])
			return 0;

		float goodCount = 0;
		for (unsigned int i = 1; i < values.size(); i++) {
			if (values[i] > otherValues[i])
				goodCount++;
		}
		return (goodCount / (values.size() - 1) >= partialCompareRatio) ? 1 : 0;
	}
	else
		return -1;
};


