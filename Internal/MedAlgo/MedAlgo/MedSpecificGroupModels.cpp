#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedAlgo/MedAlgo/MedSpecificGroupModels.h>
#include <MedAlgo/MedAlgo/MedQRF.h>

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL

MedSpecificGroupModels::MedSpecificGroupModels() {
	transpose_for_learn = false;
	transpose_for_predict = false;
	normalize_for_learn = false;
	normalize_y_for_learn = false;
	normalize_for_predict = false;

	classifier_type = MODEL_SPECIFIC_GROUPS_MODELS;
};
MedSpecificGroupModels::~MedSpecificGroupModels() {};

void MedSpecificGroupModels::set_predictors(const vector<MedPredictor *> &predictors)
{
	this->predictors = predictors;
	if (!feat_ths.empty() && predictors.size() != feat_ths.size() + 1)
		MTHROW_AND_ERR("initialize error - predictors size should by feat_ths size + 1");
}

void MedSpecificGroupModels::set_group_selection(int featNum, const vector<float> &feat_ths) {
	this->featNum = featNum;
	this->feat_ths = feat_ths;
	sort(this->feat_ths.begin(), this->feat_ths.end());
	if (!predictors.empty() && predictors.size() != feat_ths.size() + 1)
		MTHROW_AND_ERR("initialize error - predictors size should by feat_ths size + 1");
}

int MedSpecificGroupModels::selectPredictor(const float *x) const {
	float num = x[featNum];
	int i = 0;
	while (i < feat_ths.size() && num > feat_ths[i])
		++i;

	return i;
}


int MedSpecificGroupModels::Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
	if (predictors.empty() || feat_ths.empty())
		MTHROW_AND_ERR("using uninitialize predictor - please set predictors & feat_ths");
	vector<vector<float>> predictors_x(predictors.size());
	vector<vector<float>> predictors_y(predictors.size());
	vector<vector<float>> predictors_w(predictors.size());
	vector<int> predictors_nsamples(predictors.size());

	for (size_t i = 0; i < nsamples; ++i)
	{
		int predictor_index = selectPredictor(&x[i*nftrs]);
		predictors_x[predictor_index].insert(predictors_x[predictor_index].end(), &x[i*nftrs], &x[(i + 1)*nftrs]);
		predictors_y[predictor_index].insert(predictors_y[predictor_index].end(), &y[i], &y[i + 1]);
		if (w != NULL)
			predictors_w[predictor_index].insert(predictors_w[predictor_index].end(), &w[i], &w[i + 1]);
		++predictors_nsamples[predictor_index];
	}

	//now learn each predictor:
	for (size_t i = 0; i < predictors.size(); ++i) {
		predictors[i]->learn(predictors_x[i], predictors_y[i], predictors_w[i], predictors_nsamples[i], nftrs);
	}

	return 0;
}

int MedSpecificGroupModels::Predict(float *x, float *&preds, int nsamples, int nftrs) const {
	if (predictors.empty() || feat_ths.empty())
		MTHROW_AND_ERR("using uninitialize predictor - please set predictors & feat_ths");
	vector<vector<float>> predictors_x(predictors.size());
	vector<int> predictors_nsamples(predictors.size());
	vector<int> predictors_order(nsamples);

	for (size_t i = 0; i < nsamples; ++i)
	{
		int predictor_index = selectPredictor(&x[i*nftrs]);
		predictors_x[predictor_index].insert(predictors_x[predictor_index].end(), &x[i*nftrs], &x[(i + 1)*nftrs]);
		++predictors_nsamples[predictor_index];
		predictors_order[i] = predictor_index;
	}
	vector<vector<float>> predictors_preds(predictors.size());
	for (size_t i = 0; i < predictors.size(); ++i)
	{
		if (predictors_x.size() == 0)
			continue;
		predictors[i]->predict(predictors_x[i], predictors_preds[i], predictors_nsamples[i], nftrs);
	}

	//now collectby order:
	vector<int> predictors_progress(predictors.size());
	for (size_t i = 0; i < predictors_order.size(); ++i)
	{
		int p_index = predictors_order[i];
		preds[i] = predictors_preds[p_index][predictors_progress[p_index]];
		++predictors_progress[p_index];
	}

	return 0;
}

MedPredictor *cloneModel(MedPredictor *model) {
	MedQRFParams pr_qrf;
	map<string, string> empty_m;
	MedPredictor *newM;
	switch (model->classifier_type) {
	case MODEL_QRF:
		pr_qrf = MedQRFParams((*(MedQRF *)model).params);
		newM = new MedQRF(pr_qrf);
		break;
	default:
		throw invalid_argument("Unsupported Type init");
	}

	return newM;
}

MedSpecificGroupModels *MedSpecificGroupModels::clone() const {
	MedSpecificGroupModels *res = new MedSpecificGroupModels;

	res->set_group_selection(featNum, feat_ths);
	vector<MedPredictor *> cp_predictors((int)predictors.size());
	for (size_t i = 0; i < predictors.size(); ++i)
		cp_predictors[i] = cloneModel(predictors[i]);
	res->set_predictors(cp_predictors);

	return res;
}

MedPredictor *MedSpecificGroupModels::get_model(int ind) {
	return predictors[ind];
}

int MedSpecificGroupModels::model_cnt() const {
	return (int)predictors.size();
}