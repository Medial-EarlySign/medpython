#include "MedPerformance.h"

#include "Logger/Logger/Logger.h"
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
	int maxSplits = -1, missingSplit = -1;
	for (auto& idSamples : samples.idSamples) {
		int split = idSamples.split;
		if (split > maxSplits)
			maxSplits = split;
		if (split == -1)
			missingSplit = 1;
	}
	assert(missingSplit == -1 || maxSplits == -1); // Can't have both splits and non-splits !

	// Fill
	preds.resize(maxSplits + 2);
	for (auto& idSamples : samples.idSamples) {
		int split = idSamples.split;
		for (auto& sample : idSamples.samples)
			preds[split+1].push_back({ sample.prediction.back(),sample.outcome });
	}

	// Are there any splits ?
	if (maxSplits > -1)
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

// Load from predictor data
void MedClassifierPerformance::_load(MedFeaturesData& predictor_data) {

	preds.resize(predictor_data.nsplits + 1);

	vector<int> split_pos(predictor_data.nsplits, predictor_data.n_preds_per_sample - 1); // Note: if there's more than 1 pred_per_sample, we measure performance vs the LAST label ("1")
	for (unsigned int i = 0; i < predictor_data.label.size(); i++) {
		int isplit = predictor_data.splits[i];
		preds[isplit + 1].push_back(pair<float, float>(predictor_data.split_preds[isplit][split_pos[isplit]], predictor_data.label[i]));
		split_pos[isplit] += predictor_data.n_preds_per_sample;
	}

	SplitsToComplete();

}

void MedClassifierPerformance::post_load() {
	init();
	ShuffleSort();
	getPerformanceValues();

	PerformancePointers.resize(preds.size());
}

// Load train predictions
void MedClassifierPerformance::load_preds_on_train(MedFeaturesData& predictor_data) {

	preds.resize(predictor_data.nsplits);

	vector<int> split_pos(predictor_data.nsplits, predictor_data.n_preds_per_sample - 1); // Note: if there's more than 1 pred_per_sample, we measure performance vs the LAST label ("1")
	for (unsigned int i = 0; i < predictor_data.label.size(); i++) {
		int isplit = predictor_data.splits[i];
		for (int tr_split = 0; tr_split < predictor_data.nsplits; tr_split++) {
			if (tr_split != isplit) {
				preds[tr_split].push_back(pair<float, float>(predictor_data.split_preds_on_train[tr_split][split_pos[tr_split]],
					predictor_data.label[i]));
				split_pos[tr_split] += predictor_data.n_preds_per_sample;
			}
		}
	}

	post_load();
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
		random_shuffle(preds[i].begin(), preds[i].end(), rand_N);
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
		if (getPerformancePointer(set, index) < 0)
			return -1;
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

//quantize AUC calculation - calcs auc when the scores of preds are quantized
//.........................................................................................................................................
float get_preds_auc_q(const vector<float> &preds, const vector<float> &y) {
	vector<float> pred_threshold;
	map<float, vector<int>> pred_indexes;
	int tot_true_labels = 0;
	for (size_t i = 0; i < preds.size(); ++i)
	{
		pred_indexes[preds[i]].push_back((int)i);
		tot_true_labels += int(y[i] > 0);
	}
	int tot_false_labels = (int)y.size() - tot_true_labels;
	pred_threshold = vector<float>((int)pred_indexes.size());
	map<float, vector<int>>::iterator it = pred_indexes.begin();
	for (size_t i = 0; i < pred_threshold.size(); ++i)
	{
		pred_threshold[i] = it->first;
		++it;
	}
	sort(pred_threshold.begin(), pred_threshold.end());

	//From up to down sort:
	int t_cnt = 0;
	int f_cnt = 0;
	vector<float> true_rate = vector<float>((int)pred_indexes.size());
	vector<float> false_rate = vector<float>((int)pred_indexes.size());
	int st_size = (int)pred_threshold.size() - 1;
	for (int i = st_size; i >= 0; --i)
	{
		vector<int> indexes = pred_indexes[pred_threshold[i]];
		//calc AUC status for this step:
		for (int ind : indexes)
		{
			bool true_label = y[ind] > 0;
			t_cnt += int(true_label);
			f_cnt += int(!true_label);
		}
		true_rate[st_size - i] = float(t_cnt) / tot_true_labels;
		false_rate[st_size - i] = float(f_cnt) / tot_false_labels;
	}

	float auc = false_rate[0] * true_rate[0] / 2;
	for (size_t i = 1; i < true_rate.size(); ++i)
	{
		auc += (false_rate[i] - false_rate[i - 1]) * (true_rate[i - 1] + true_rate[i]) / 2;
	}
	return auc;
}

// AUC
//.........................................................................................................................................
float get_preds_auc(vector<float> &preds, vector<float> &y) {
	vector<pair<float, float>> preds_y;

	if (preds.size() != y.size())
		return -1;

	preds_y.resize(preds.size());
	for (int i = 0; i < preds.size(); i++) {
		preds_y[i].first = preds[i];
		preds_y[i].second = y[i];
		//preds_y[i].second = (y[i] > 3);
	}

	// Sort from high score to low
	sort(preds_y.begin(), preds_y.end(), [](const pair<float, float>& left, const pair<float, float>& right) { return left.first > right.first; });

	// Identify cutting-point for partial AUC
	int tot_nneg = 0;
	for (unsigned int i = 0; i < preds.size(); i++) {
		if (preds_y[i].second <= 0)
			tot_nneg++;
	}
	int target_nneg = (int)(tot_nneg);

	// Loop - for each Negative, count all positives above it ...
	unsigned long long auc = 0;
	unsigned long long nneg = 0, npos = 0;
	for (unsigned int i = 0; i < preds.size(); i++) {
		if (preds_y[i].second > 0)
			npos++;
		else {
			auc += npos;
			nneg++;
			if (nneg == target_nneg)
				break;
		}
	}

	return (float)(((double)auc) / ((double)npos*(double)target_nneg));
}

// Collect cnts : TP,FP,FN,TN per positive rate (as given by size)
//.........................................................................................................................................
int get_preds_perf_cnts(vector<float> &preds, vector<float> &y, vector<float> &size, int direction, vector<vector<int>> &cnts)
{
	vector<pair<float, float>> preds_y;

	if (preds.size() != y.size())
		return -1;

	preds_y.resize(preds.size());
	for (int i = 0; i < preds.size(); i++) {
		if (direction > 0)
			preds_y[i].first = preds[i];
		else
			preds_y[i].first = -preds[i];

		preds_y[i].second = y[i];
		//preds_y[i].second = (y[i] > 3);
	}

	// Sort from high score to low
	sort(preds_y.begin(), preds_y.end(), [](const pair<float, float>& left, const pair<float, float>& right) { return left.first > right.first; });

	cnts.clear();
	for (auto sz : size) {
		int pos = (int)(sz * (float)preds_y.size());
		int a = 0, b = 0, c = 0, d = 0;
		for (int i = 0; i <= pos; i++) {
			if (preds_y[i].second > 0)
				a++;
			else
				b++;
		}
		for (int i = pos + 1; i < preds_y.size(); i++) {
			if (preds_y[i].second > 0)
				c++;
			else
				d++;
		}
		cnts.push_back({ a,b,c,d });
		//MLOG("sz %f pos %d/%d cnts %d %d %d %d\n", sz, pos, preds_y.size(), a, b, c, d);
	}

	return 0;
}

// Translate counts (TP,FP,FN,TN) to performance measurements (snes,spec,ppv,rr)
//.........................................................................................................................................
int cnts_to_perf(vector<int> &cnt, float &sens, float &spec, float &ppv, float &rr)
{
	float a = (float)cnt[0];
	float b = (float)cnt[1];
	float c = (float)cnt[2];
	float d = (float)cnt[3];

	float epsilon = (float)1e-5;
	sens = a / (a + c + epsilon);
	spec = d / (b + d + epsilon);
	ppv = a / (a + b + epsilon);
	rr = (a / (a + b + epsilon)) / (c / (c + d + epsilon));

	return 0;
}

//.........................................................................................................................................

// given multicategory probs (or probs-like) predictions generates a single prediction of the categ with max prob for each sample
int multicateg_get_max_pred(vector<float> &probs, int nsamples, int ncateg, vector<float> &max_pred)
{
	int i, j;

	max_pred.resize(nsamples);
	for (i = 0; i < nsamples; i++) {
		float max = probs[i*ncateg];
		int max_j = 0;
		for (j = 1; j<ncateg; j++) {
			if (probs[i*ncateg + j] > max) {
				max = probs[i*ncateg + j];
				max_j = j;
			}
		}
		max_pred[i] = (float)max_j;
	}

	return 0;
}

// given multicategory probs (or probs-like) predictions generates a single prediction of a weighted average of categories
int multicateg_get_avg_pred(vector<float> &probs, int nsamples, int ncateg, vector<float> &avg_pred)
{
	int i, j;

	avg_pred.resize(nsamples);
	for (i = 0; i < nsamples; i++) {
		float sum = 0, sum_p = 0;
		for (j = 0; j < ncateg; j++) {
			sum += (float)j*probs[i*ncateg + j];
			sum_p += probs[i*ncateg + j];
		}
		if (sum_p <= 0) sum_p = (float)1e-15;
		avg_pred[i] = sum / sum_p;
	}

	return 0;
}

// given multicategory probs (or probs-like) predictions gets the classification error rate and the rms (also for the avg preds)
int multicateg_get_error_rate(vector<float> &probs, vector<float> &y, int nsamples, int ncateg, float &err_rate, float &rms, float &avg_rms)
{
	vector<float> max_preds, avg_preds;

	multicateg_get_max_pred(probs, nsamples, ncateg, max_preds);
	multicateg_get_avg_pred(probs, nsamples, ncateg, avg_preds);

	err_rate = 0;
	rms = 0;
	avg_rms = 0;

	int i;

	float fact = (float)1 / (float)nsamples;
	for (i = 0; i < nsamples; i++) {

		if (max_preds[i] != y[i]) err_rate++;
		rms += fact*(max_preds[i] - y[i])*(max_preds[i] - y[i]);
		avg_rms += fact*(avg_preds[i] - y[i])*(avg_preds[i] - y[i]);

	}

	err_rate /= (float)nsamples;

	return 0;
}

//.........................................................................................................................................
int get_quantized_breakdown(vector<float> &preds, vector<float> &y, vector<float> &bounds, MedMat<int> &counts)
{
	int i, ip, iy;
	int nb = (int)bounds.size() - 1;

	counts.resize(nb, nb);
	counts.zero();

	for (i = 0; i < preds.size(); i++) {

		ip = 0;
		while (ip<nb && preds[i]>bounds[ip + 1]) ip++;
		iy = 0;
		while (iy<nb && y[i]>bounds[iy + 1]) iy++;

		if (preds[i] <= bounds[ip + 1] && y[i] <= bounds[iy + 1]) {
			counts(ip, iy)++;
		}

	}

	return 0;
}

//.........................................................................................................................................
void print_quantized_breakdown(MedMat<int> &cnt, vector<float> &bounds)
{
	MOUT("Quantized results distribution:\n");
	for (int i = 0; i < cnt.nrows; i++) {
		MOUT("preds %7.2f - %7.2f ::", bounds[i], bounds[i + 1]);
		for (int j = 0; j < cnt.ncols; j++) {
			MOUT(" %6d", cnt(i, j));
		}
		MOUT("\n");
	}
}

// Chi-Square
//.........................................................................................................................................
double get_chi2_n_x_m(vector<int> &cnts, int n, int m, vector<double> &exp)
{
	int i, j;

	if (cnts.size() != n*m)
		return -1;

	vector<double> s_n(n, 0);
	vector<double> s_m(m, 0);
	double sum = 0;

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++) {
			double c = (double)cnts[i*m + j];
			s_n[i] += c;
			s_m[j] += c;
			sum += c;
		}

	if (sum <= 0) return -1;

	exp.resize(n*m);
	double score = 0;
	double epsilon = 1e-5;
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++) {
			double e = (s_n[i] / sum)*s_m[j];
			exp[i*m + j] = e;
			if (e > epsilon) {
				double d = (double)cnts[i*m + j] - e;
				score += (d / e)*d;
			}
		}

	return score;
}

//.........................................................................................................................................
double get_chi2_n_x_m(vector<int> &cnts, int n, int m)
{
	vector<double> exp;
	return(get_chi2_n_x_m(cnts, n, m, exp));
}


// Pearson Correlations
//.........................................................................................................................................
float get_pearson_corr(float *v1, float *v2, int len)
{
	if (len == 0)
		return -2.0;

	double sx, sy, sxy, sxx, syy, n;

	sx = sy = sxy = sxx = syy = 0;

	double fact = 1e-5;
	for (int i = 0; i < len; i++) {
		sx += fact*v1[i];
		sy += fact*v2[i];
		sxx += fact*v1[i] * v1[i];
		syy += fact*v2[i] * v2[i];
		sxy += fact*v1[i] * v2[i];
	}

	n = (double)len;

	sx /= fact*n;
	sy /= fact*n;
	sxx /= fact*n;
	syy /= fact*n;
	sxy /= fact*n;

	double c1 = sxy - sx*sy;
	double c2 = sxx - sx*sx;
	double c3 = syy - sy*sy;

	double epsilon = 1e-8;
	if (c2 < epsilon || c3 < epsilon)
		return 0;
	return (float)(c1 / (sqrt(c2)*sqrt(c3)));
}

float get_pearson_corr(vector<float> &v1, vector<float> &v2) {
	return get_pearson_corr(VEC_DATA(v1), VEC_DATA(v2), (int)v1.size());
}

// Pearson Correaltion after removing missing values. return number of values left in n.
float get_pearson_corr(vector<float> &v1, vector<float> &v2, int &n, float missing_val) {

	vector<float> clean_v1, clean_v2;
	for (unsigned int i = 0; i < v1.size(); i++) {
		if (v1[i] != missing_val && v2[i] != missing_val) {
			clean_v1.push_back(v1[i]);
			clean_v2.push_back(v2[i]);
		}
	}

	n = (int)clean_v1.size(); 
	return get_pearson_corr(VEC_DATA(clean_v1), VEC_DATA(clean_v2), (int)clean_v1.size());
}

// Mutual information of binned-vectors
//.........................................................................................................................................
int get_mutual_information(vector<int>& x, vector<int>& y, int &n, float& mi) {

	// Sanity
	if (y.size() != x.size()) {
		MERR("Size mismatch. Quitting\n");
		return -1;
	}

	// Count bins
	int nXbins = 0, nYbins = 0;
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] + 1 > nXbins)
			nXbins = x[i] + 1;
		if (y[i] + 1 > nYbins)
			nYbins = y[i] + 1;
	}

	// Collect
	vector<int> xCounts(nXbins,0), yCounts(nYbins,0), coCounts(nXbins*nYbins,0);
	n = 0;
	for (unsigned int i = 0; i<x.size(); i++) {
		if (x[i] >= 0 && y[i] >= 0) {
			xCounts[x[i]]++;
			yCounts[y[i]]++;
			coCounts[y[i] * nXbins + x[i]]++;
			n++;
		}
	}

	if (n < 2) {
		MLOG_V("Not enough common non-missing entries for mutual information.\n");
		mi = -1;
	}

	mi = get_mutual_information(xCounts, yCounts, coCounts, n);

	return 0;
}

// Mutual information from counts
float get_mutual_information(vector<int>& xCounts, vector<int>& yCounts, vector<int> coCounts, int n) {

	double mi = 0;
	int nXbins = (int)xCounts.size(); 
	int nYbins = (int)yCounts.size();

	for (int iX = 0; iX < nXbins; iX++) {
		for (int iY = 0; iY < nYbins; iY++) {
			if (coCounts[iY*nXbins + iX] != 0) {
				double p = (coCounts[iY*nXbins + iX] + 0.0) / n;
				double px = (xCounts[iX] + 0.0) / n;
				double py = (yCounts[iY] + 0.0) / n;

				mi += p * log(p / px / py) / log(2.0);
			}
		}
	}

	return (float) mi;
}

// Distance Correlations
//.........................................................................................................................................
// Get Distances matrix
void get_dMatrix(vector<float>& values, MedMat<float>& dMatrix, float missing_value) {

	int n = (int)values.size();
	dMatrix.resize(n,n);

	// Matrix + norms
	vector<double> norm(n, 0);
	vector<int> counts(n, 0);
	double totNorm = 0;
	int totCount = 0;

	for (int i = 1; i < n; i++) {
		if (values[i] == missing_value) {
			for (int j = 0; j < i; j++)
				dMatrix(i, j) = -1;
		}
		else {
			for (int j = 0; j < i; j++) {
				if (values[j] == missing_value)
					dMatrix(i, j) = -1;
				else {
					dMatrix(i, j) = fabs(values[i] - values[j]);
					norm[i] += 2 * dMatrix(i, j); counts[i] += 2;
					norm[j] += 2 * dMatrix(i, j); counts[j] += 2;
					totNorm += 2 * dMatrix(i, j); totCount += 2;
				}
			}
		}
	}

	// Normalize
	for (int i = 0; i < n; i++)
		norm[i] /= counts[i];
	totNorm /= totCount;

	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			if (dMatrix(i, j) != -1)
				dMatrix(i, j) = dMatrix(i, j) - (float)(norm[i] + norm[j] - totNorm);
		}
	}
}

// Get Distance variance
float get_dVar(MedMat<float>& dMatrix) {

	int n = dMatrix.nrows;

	double sum = 0;
	int num = 0;
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			if (dMatrix(i, j) != -1) {
				sum += 2 * dMatrix(i, j)*dMatrix(i, j);
				num += 2;
			}
		}
	}

	if (num)
		return (float)(sum / num);
	else
		return -1;
}

// Get Distance covariance
float get_dCov(MedMat<float>& xDistMat, MedMat<float>& yDistMat) {

	int n = xDistMat.nrows;

	double sum = 0;
	int num = 0;
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			if (xDistMat(i, j) != -1 && yDistMat(i, j) != -1) {
				sum += 2 * xDistMat(i, j)*yDistMat(i, j);
				num += 2;
			}
		}
	}

	if (num)
		return (float)(sum / num);
	else
		return -1;
}


