//
// MedAlgo - unified wrappers for prediction algorithms
//

#include "MedAlgo.h"
#include "MedXGB.h"
#include "MedDeepBit.h"
#include "MedUtils/MedUtils/MedIO.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedLightGBM.h"
#include "MedLinearModel.h"

#if NEW_COMPLIER
#include "MedVW.h"
#endif

#include <thread>

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

unordered_map<int, string> predictor_type_to_name = {
	{ MODEL_LINEAR_MODEL , "linear_model"} ,
	{ MODEL_GD_LINEAR , "gdlm" },
	{ MODEL_QRF , "qrf" },
	{ MODEL_KNN , "knn" },
	{ MODEL_MARS , "mars" },
	{ MODEL_GBM , "gbm" },
	{ MODEL_BP , "BP" },
	{ MODEL_MULTI_CLASS , "multi_class" },
	{ MODEL_XGB , "xgb" },
	{ MODEL_LASSO , "lasso" },
	{ MODEL_MIC_NET , "micNet" },
	{ MODEL_BOOSTER , "booster" },
	{ MODEL_DEEP_BIT , "deep_bit" },
	{ MODEL_LIGHTGBM , "lightgbm" },
	{ MODEL_SVM , "svm" },
	{ MODEL_LIGHTGBM , "lightgbm" },
	{ MODEL_LINEAR_SGD , "linear_sgd" },
	{ MODEL_SPECIFIC_GROUPS_MODELS, "multi_models" },
	{ MODEL_VW, "vw" },
	{ MODEL_TQRF, "tqrf"}
};
//=======================================================================================
// MedPredictor
//=======================================================================================
// Model types
MedPredictorTypes predictor_name_to_type(const string& model_name) {
	for (auto it = predictor_type_to_name.begin(); it != predictor_type_to_name.end(); ++it)
		if (it->second == model_name) {
			return MedPredictorTypes(it->first);
		}
	return MODEL_LAST;
}

// Initialization
//.......................................................................................
MedPredictor * MedPredictor::make_predictor(string model_type) {

	return make_predictor(predictor_name_to_type(model_type));
}

//.......................................................................................
MedPredictor * MedPredictor::make_predictor(string model_type, string init_string) {

	return make_predictor(predictor_name_to_type(model_type), init_string);
}

//.......................................................................................
MedPredictor * MedPredictor::make_predictor(MedPredictorTypes model_type) {

	if (model_type == MODEL_LINEAR_MODEL)
		return new MedLM;
	else if (model_type == MODEL_GD_LINEAR)
		return new MedGDLM;
	else if (model_type == MODEL_QRF)
		return new MedQRF;
	else if (model_type == MODEL_KNN)
		return new MedKNN;
	else if (model_type == MODEL_MARS)
		return new MedMars;
	else if (model_type == MODEL_GBM)
		return new MedGBM;
	else if (model_type == MODEL_BP)
		return new MedGBM;
	else if (model_type == MODEL_MULTI_CLASS)
		return new MedMultiClass;
	else if (model_type == MODEL_XGB)
		return new MedXGB;
	else if (model_type == MODEL_LASSO)
		return new MedLasso;
	else if (model_type == MODEL_MIC_NET)
		return new MedMicNet;
	else if (model_type == MODEL_BOOSTER)
		return new MedBooster;
	else if (model_type == MODEL_DEEP_BIT)
		return new MedDeepBit;
	else if (model_type == MODEL_LIGHTGBM)
		return new MedLightGBM;
	else if (model_type == MODEL_SPECIFIC_GROUPS_MODELS)
		return new MedSpecificGroupModels;
	else if (model_type == MODEL_SVM)
		return new MedSvm;
	else if (model_type == MODEL_TQRF)
		return new MedTQRF;
#if NEW_COMPLIER
	else if (model_type == MODEL_VW)
		return new MedVW;
#endif
	else
		return NULL;

}
MedPredictor * MedPredictor::make_predictor(MedPredictorTypes model_type, string init_string) {

	MedPredictor *newPred = make_predictor(model_type);
	newPred->init_from_string(init_string);

	return newPred;
}
//.......................................................................................
int MedPredictor::init_from_string(string text) {

	MLOG("MedPredictor init from string [%s] (classifier type is %d)\n", text.c_str(), classifier_type);

	// parse text of the format "Name = Value ; Name = Value ; ..."

	if (classifier_type == MODEL_MIC_NET) {
		cerr << "But we are going to call mic net version directly\n";
		MedMicNet *mic = (MedMicNet *)this;
		cerr << "before\n";
		int rc = mic->init_from_string(text);
		cerr << "after " << rc << "\n";
		return rc;
	}

	if (classifier_type == MODEL_BOOSTER) {
		cerr << "But we are going to call booster version directly\n";
		MedBooster *med_b = (MedBooster *)this;
		cerr << "before\n";
		int rc = med_b->init_from_string(text);
		cerr << "after " << rc << "\n";
		return rc;
	}

	if (classifier_type == MODEL_LIGHTGBM) {
		MedLightGBM *med_light = (MedLightGBM *)this;
		return med_light->init_from_string(text);
	}
#if NEW_COMPLIER
	if (classifier_type == MODEL_VW) {
		MedVW *vw = (MedVW *)this;
		return vw->init_from_string(text);
	}
#endif
	// remove white spaces
	text.erase(remove_if(text.begin(), text.end(), ::isspace), text.end());

	map<string, string> init_map;
	if (initialization_text_to_map(text, init_map) == -1)
		return -1;

	for (auto rec : init_map)
		MLOG("Initializing predictor with \'%s\' = \'%s\'\n", rec.first.c_str(), rec.second.c_str());

	init(init_map);

	return 0;
}

//.......................................................................................
void MedPredictor::prepare_x_mat(MedMat<float> &x, vector<float> &wgts, int &nsamples, int &nftrs, bool transpose_needed)
{
	if ((transpose_needed && !x.transposed_flag) || (!transpose_needed && x.transposed_flag)) {
		//		MLOG("transposing matrix\n");
		x.transpose();
	}

	if (x.transposed_flag) {
		nsamples = x.ncols;
		nftrs = x.nrows;
	}
	else {
		nsamples = x.nrows;
		nftrs = x.ncols;
	}
}

string norm_feature_name(const string &feat_name) {
	return  feat_name.substr(0, 3) != "FTR" || feat_name.find_first_of('.') == string::npos ? feat_name :
		feat_name.substr(feat_name.find_first_of('.') + 1);
}

//.......................................................................................
int MedPredictor::learn(MedMat<float> &x, MedMat<float> &y, vector<float> &wgts)
{
	if (!x.signals.empty())
		model_features = x.signals;
	else {
		model_features.clear();
		features_count = x.ncols;
	}
	int nsamples, nftrs;

	// patch for micNet
	if (classifier_type == MODEL_MIC_NET) {
		MedMicNet *mic = (MedMicNet *)this;
		cerr << "running micNet learn()\n";
		return mic->learn(x, y);
	}
	// patch for booster
	if (classifier_type == MODEL_BOOSTER) {
		MedBooster *med_b = (MedBooster *)this;
		cerr << "running MedBooster learn()\n";
		return med_b->learn(x, y);
	}

	// ToDo : sanity check of sizes (nonzero, matching x,y dimensions)

	if (normalize_for_learn && !x.normalized_flag) {
		MERR("Learner Requires normalized matrix. Quitting\n");
		return -1;
	}

	prepare_x_mat(x, wgts, nsamples, nftrs, transpose_for_learn);
	if (normalize_y_for_learn && !y.normalized_flag)
		y.normalize();

	return Learn(x.data_ptr(), y.data_ptr(), VEC_DATA(wgts), nsamples, nftrs);
}

//.......................................................................................
int MedPredictor::learn(MedMat<float> &x, vector<float> &y, vector<float> &wgts) {
	if (!x.signals.empty())
		model_features = x.signals;
	else {
		model_features.clear();
		features_count = x.ncols;
	}
	int nsamples, nftrs;

	if (normalize_for_learn && !x.normalized_flag) {
		MERR("Learner Requires normalized matrix. Quitting\n");
		return -1;
	}

	prepare_x_mat(x, wgts, nsamples, nftrs, transpose_for_learn);

	return Learn(x.data_ptr(), y.data(), VEC_DATA(wgts), nsamples, nftrs);
}

//.......................................................................................
int MedPredictor::learn(vector<float> &x, vector<float> &y, vector<float> &w, int n_samples, int n_ftrs)
{
	features_count = n_ftrs;
	return(Learn(VEC_DATA(x), VEC_DATA(y), VEC_DATA(w), n_samples, n_ftrs));
}

//.......................................................................................
int MedPredictor::predict(MedMat<float> &x, vector<float> &preds) {
	if (!model_features.empty()) {//test names of entered matrix:
		if (model_features.size() != x.ncols)
			MTHROW_AND_ERR("(1) Learned Feature model size was %d, request feature size for predict was %d\n",
				(int)model_features.size(), (int)x.ncols);

		if (!x.signals.empty()) //can compare names
			for (int feat_num = 0; feat_num < model_features.size(); ++feat_num)
				if (norm_feature_name(x.signals[feat_num]) != norm_feature_name(model_features[feat_num]))
					MTHROW_AND_ERR("Learned Features are the same. feat_num=%d. in learning was %s, now recieved %s\n",
						(int)model_features.size(), model_features[feat_num].c_str(), x.signals[feat_num].c_str());
	}
	else if (features_count > 0 && features_count != x.ncols)
		MTHROW_AND_ERR("(2) Learned Feature model size was %d, request feature size for predict was %d\n",
			features_count, (int)x.ncols);

	int nsamples, nftrs;
	vector<float> w;

	// patch for micNet
	if (classifier_type == MODEL_MIC_NET) {
		MedMicNet *mic = (MedMicNet *)this;
		cerr << "running micNet predict()\n";
		return mic->predict(x, preds);
	}
	// patch for booster
	if (classifier_type == MODEL_BOOSTER) {
		MedBooster *med_b = (MedBooster *)this;
		cerr << "running MedBooster predict()\n";
		return med_b->predict(x, preds);
	}

	if (classifier_type == MODEL_TQRF) {
		return ((MedTQRF *)this)->Predict(x, preds);
	}

	if (normalize_for_predict && !x.normalized_flag) {
		MERR("Predictor Requires normalized matrix. Quitting\n");
		return -1;
	}


	prepare_x_mat(x, w, nsamples, nftrs, transpose_for_predict);

	preds.resize(nsamples*n_preds_per_sample());
	float *_preds = &(preds[0]);
	//	MLOG("MedMat,vector call: preds size is %d n_preds_per_sample %d nsamples %d\n",preds.size(),n_preds_per_sample(),nsamples);
	return Predict(x.data_ptr(), _preds, nsamples, nftrs);
}

struct pred_thread_info {
	int id;
	int from_sample;
	int to_sample;
	float *preds;
	MedMat<float> *x;
	int nftrs;
	int n_preds_per_sample;
	int rc;
	int state;
};

void MedPredictor::predict_thread(void *p)
//void MedPredictor::predict_thread()
{

	pred_thread_info *tp = (pred_thread_info *)p;

	//MLOG("Start thread %d : from: %d to: %d\n",tp->id,tp->from_sample,tp->to_sample);
	float *x = &(tp->x->m[tp->from_sample * tp->nftrs]);
	float *preds = &(tp->preds[tp->from_sample * n_preds_per_sample()]);
	int nsamples = tp->to_sample - tp->from_sample + 1;
	int nftrs = tp->nftrs;

	tp->rc = Predict(x, preds, nsamples, nftrs);
	//MLOG("End thread %d : from: %d to: %d\n",tp->id,tp->from_sample,tp->to_sample);

	// signing job ended
	tp->state = 1;
}

//.......................................................................................
int MedPredictor::threaded_predict(MedMat<float> &x, vector<float> &preds, int nthreads) {
	if (!model_features.empty()) {//test names of entered matrix:
		if (model_features.size() != x.ncols)
			MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
				(int)model_features.size(), (int)x.ncols);

		if (!x.signals.empty()) //can compare names
			for (int feat_num = 0; feat_num < model_features.size(); ++feat_num)
				if (norm_feature_name(x.signals[feat_num]) != norm_feature_name(model_features[feat_num]))
					MTHROW_AND_ERR("Learned Features are the same. feat_num=%d. in learning was %s, now recieved %s\n",
						(int)model_features.size(), model_features[feat_num].c_str(), x.signals[feat_num].c_str());
	}
	else if (features_count > 0 && features_count != x.ncols)
		MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
			features_count, (int)x.ncols);

	int nsamples, nftrs;
	vector<float> w;

	if (transpose_for_predict) {
		MERR("!!!!!! UNSUPORTED !!!! --> Currently threaded_predict does not support transposed matrices for predictions");
		return -1;
	}

	if (normalize_for_predict && !x.normalized_flag) {
		MERR("Predictor Requires normalized matrix. Quitting\n");
		return -1;
	}

	prepare_x_mat(x, w, nsamples, nftrs, transpose_for_predict);
	preds.resize(nsamples*n_preds_per_sample());

	int th_nsamples = nsamples / nthreads;
	vector<pred_thread_info> tp(nthreads);
	for (int i = 0; i < nthreads; i++) {
		tp[i].id = i;
		tp[i].from_sample = i*th_nsamples;
		tp[i].to_sample = min((i + 1)*th_nsamples - 1, nsamples - 1);
		tp[i].preds = VEC_DATA(preds);
		tp[i].x = &x;
		tp[i].nftrs = nftrs;
		tp[i].rc = 0;
		tp[i].state = 0;
	}


	// sending threads
	vector<thread> th_handle(nthreads);
	for (int i = 0; i < nthreads; i++) {
		//		MLOG("Sending Thread %d\n",i);
		//		th_handle[i] = std::thread(&MedPredictor::predict_thread, (void *)&tp[i]);
		th_handle[i] = std::thread(&MedPredictor::predict_thread, this, (void *)&tp[i]);
	}

	int n_state = 0;
	while (n_state < nthreads) {
		this_thread::sleep_for(chrono::milliseconds(10));
		n_state = 0;
		for (int i = 0; i < nthreads; i++)
			n_state += tp[i].state;
	}
	for (int i = 0; i < nthreads; i++)
		th_handle[i].join();

	for (int i = 0; i < nthreads; i++)
		if (tp[i].rc != 0)
			return -1;
	return 0;

}


//.......................................................................................
int MedPredictor::predict(vector<float> &x, vector<float> &preds, int n_samples, int n_ftrs) {
	if (!model_features.empty()) {
		if (model_features.size() != n_ftrs)
			MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
				(int)model_features.size(), (int)x.size());
	}
	else if (features_count > 0 && features_count != n_ftrs)
		MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
			features_count, n_ftrs);

	preds.resize(n_samples*n_preds_per_sample());
	float *_preds = &(preds[0]);
	return Predict(VEC_DATA(x), _preds, n_samples, n_ftrs);
}


// (De)Serialize
//.......................................................................................
size_t MedPredictor::get_predictor_size() {
	return sizeof(classifier_type) + get_size();
}

//.......................................................................................
size_t MedPredictor::predictor_serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &classifier_type, sizeof(MedPredictorTypes)); ptr += sizeof(MedPredictorTypes);
	ptr += serialize(blob + ptr);

	return ptr;
}

//.......................................................................................

void MedPredictor::print(FILE *fp, const string& prefix) { fprintf(fp, "%s : No Printing method defined\n", prefix.c_str()); }

//===========================================================================================
// MedPredictor MedFeaturesData API
//===========================================================================================
// Cross Validation
int MedPredictor::cross_validate_splits(MedFeaturesData& data) {

	data.split_preds.resize(data.nsplits);
	data.split_preds_on_train.resize(data.nsplits);

	for (int isplit = 0; isplit < data.nsplits; isplit++) {
		if (learn(data, isplit) == -1) {
			fprintf(stderr, "Learning failed\n");
			return -1;
		}

		if (data.predict_on_train) {
			if (predict_on_train(data, isplit) == -1) {
				fprintf(stderr, "Predict on train failed\n");
				return -1;
			}
		}

		if (predict(data, isplit) == -1) {
			fprintf(stderr, "Predict failed\n");
			return -1;
		}
	}

	return 0;
}

void MedPredictor::build_learning_x_mat_for_split(MedFeaturesData& ftrs_data, vector<float>& signal, int isplit, MedMat<float>& x) {
	x.normalized_flag = 1;

	for (int icol = 0; icol < x.ncols; icol++) {
		string& name = ftrs_data.signals[icol];

		int split_irow = 0;
		for (int irow = 0; irow < ftrs_data.splits.size(); irow++) {
			if (ftrs_data.splits[irow] != isplit)
				signal[split_irow++] = ftrs_data.data[name][irow];
		}

		ftrs_data.cleaners[name][isplit].clear(signal);
		ftrs_data.cleaners[name][isplit].normalize(signal);

		for (unsigned irow = 0; irow < signal.size(); irow++)
			x(irow, icol) = signal[irow];
	}
}
// Learning ....
int MedPredictor::learn(MedFeaturesData& ftrs_data, int isplit) {

	MedMat<float> x((int)ftrs_data.get_learning_nrows(isplit), (int)ftrs_data.signals.size());
	MedMat<float> y(x.nrows, 1);
	vector<float> signal(x.nrows);
	//save model features names
	model_features.clear();
	for (auto it = ftrs_data.data.begin(); it != ftrs_data.data.end(); ++it)
		model_features.push_back(it->first);

	// Build X
	build_learning_x_mat_for_split(ftrs_data, signal, isplit, x);

	// Build Y
	int split_irow = 0;
	for (int irow = 0; irow < ftrs_data.splits.size(); irow++) {
		if (ftrs_data.splits[irow] != isplit)
			signal[split_irow++] = ftrs_data.label[irow];
	}
	if (ftrs_data.label_cleaners.size() > 0) {
		ftrs_data.label_cleaners[isplit].clear(signal);
		ftrs_data.label_cleaners[isplit].normalize(signal);
	}

	for (unsigned irow = 0; irow < signal.size(); irow++)
		y(irow, 0) = signal[irow];
	y.normalized_flag = 1;

	// Learn
	if (learn(x, y) == -1)
		return -1;

	ftrs_data.n_preds_per_sample = n_preds_per_sample();

	return 0;
}

int MedPredictor::learn(MedFeaturesData& ftrs_data) {

	// Build X
	MedMat<float> x((int)ftrs_data.label.size(), (int)ftrs_data.signals.size());
	MedMat<float> y(x.nrows, 1);
	x.normalized_flag = 1;
	//save model features names
	model_features.clear();
	for (auto it = ftrs_data.data.begin(); it != ftrs_data.data.end(); ++it)
		model_features.push_back(it->first);

	vector<float> signal(x.nrows);
	for (int icol = 0; icol < x.ncols; icol++) {
		string& name = ftrs_data.signals[icol];

		signal = ftrs_data.data[name];
		ftrs_data.cleaners[name][ftrs_data.nsplits].clear(signal);
		ftrs_data.cleaners[name][ftrs_data.nsplits].normalize(signal);

		for (unsigned irow = 0; irow < signal.size(); irow++)
			x(irow, icol) = signal[irow];
	}

	// Build Y
	signal = ftrs_data.label;
	ftrs_data.label_cleaners[ftrs_data.nsplits].clear(signal);
	ftrs_data.label_cleaners[ftrs_data.nsplits].normalize(signal);

	for (unsigned irow = 0; irow < signal.size(); irow++)
		y(irow, 0) = signal[irow];
	y.normalized_flag = 1;

	// Learn
	if (learn(x, y) == -1)
		return -1;

	return 0;
}

int MedPredictor::learn(MedFeatures& ftrs_data) {

	vector<string> dummy_names;
	return learn(ftrs_data, dummy_names);
}

int MedPredictor::learn(MedFeatures& ftrs_data, vector<string>& names) {
	//save model features names
	model_features.clear();
	for (auto it = ftrs_data.data.begin(); it != ftrs_data.data.end(); ++it)
		model_features.push_back(it->first);


	if (classifier_type == MODEL_TQRF) {
		return (((MedTQRF *)this)->Learn(ftrs_data));
	}

	// Build X
	MedMat<float> x;
	ftrs_data.get_as_matrix(x, names);

	MLOG("MedPredictor::learn() from MedFeatures, got train matrix of %d x %d\n", x.nrows, x.ncols);

	// Labels
	MedMat<float> y(x.nrows, 1);
	for (int i = 0; i < y.nrows; i++)
		y(i, 0) = ftrs_data.samples[i].outcome;

	// Weights
	if (ftrs_data.weights.size())
		return learn(x, y, ftrs_data.weights);
	else
		return learn(x, y);
}

int MedPredictor::learn_prob_calibration(MedMat<float> &x, vector<float> &y,
	vector<float> &min_range, vector<float> &max_range, vector<float> &map_prob,
	int min_bucket_size, float min_score_jump, float min_prob_jump, bool fix_prob_order) {
	// > min and <= max

	//add mapping from model score to probabilty based on big enough bins of score
	//get prediction for X:
	vector<float> preds;
	predict(x, preds);

	unordered_map<float, vector<int>> score_to_indexes;
	vector<float> unique_scores;
	for (size_t i = 0; i < preds.size(); ++i)
	{
		if (score_to_indexes.find(preds[i]) == score_to_indexes.end())
			unique_scores.push_back(preds[i]);
		score_to_indexes[preds[i]].push_back((int)i);
	}
	sort(unique_scores.begin(), unique_scores.end());
	int sz = (int)unique_scores.size();

	float curr_max = (float)INT32_MAX; //unbounded
	float curr_min = curr_max;
	int pred_sum = 0;
	int curr_cnt = 0;
	vector<int> bin_cnts;
	for (int i = sz - 1; i >= 0; --i)
	{
		//update values curr_cnt, pred_avg
		for (int ind : score_to_indexes[unique_scores[i]])
			pred_sum += int(y[ind] > 0);
		curr_cnt += (int)score_to_indexes[unique_scores[i]].size();

		if (curr_cnt > min_bucket_size && curr_max - unique_scores[i] > min_score_jump) {
			//flush buffer
			curr_min = unique_scores[i];
			max_range.push_back(curr_max);
			min_range.push_back(curr_min);
			map_prob.push_back(float(double(pred_sum) / curr_cnt));
			bin_cnts.push_back(curr_cnt);

			//init new buffer:
			curr_cnt = 0;
			curr_max = unique_scores[i];
			pred_sum = 0;
		}
	}
	if (curr_cnt > 0) {
		//flush last buffer
		curr_min = (float)INT32_MIN;
		max_range.push_back(curr_max);
		min_range.push_back(curr_min);
		map_prob.push_back(float(double(pred_sum) / curr_cnt));
		bin_cnts.push_back(curr_cnt);
	}

	//unite similar prob bins:
	vector<int> ind_to_unite;
	for (int i = (int)map_prob.size() - 1; i >= 1; --i)
		if (abs(map_prob[i] - map_prob[i - 1]) < min_prob_jump ||
			(fix_prob_order && map_prob[i] > map_prob[i - 1])) { //unite bins:
			ind_to_unite.push_back(i);
			int new_count = bin_cnts[i] + bin_cnts[i - 1];
			float new_prob = (map_prob[i] * bin_cnts[i] + map_prob[i - 1] * bin_cnts[i - 1]) / new_count;
			float max_th = max_range[i - 1];
			float min_th = min_range[i];
			min_range[i - 1] = min_th;
			max_range[i - 1] = max_th;
			map_prob[i - 1] = new_prob;
			bin_cnts[i - 1] = new_count;
		}

	//unite from end to start:
	for (int i = 0; i < ind_to_unite.size(); ++i)
	{
		int unite_index = ind_to_unite[i];
		//delete old records:
		min_range.erase(min_range.begin() + unite_index);
		max_range.erase(max_range.begin() + unite_index);
		map_prob.erase(map_prob.begin() + unite_index);
		bin_cnts.erase(bin_cnts.begin() + unite_index);
	}

	MLOG("Created %d bins for mapping prediction scores to probabilities\n", map_prob.size());
	for (size_t i = 0; i < map_prob.size(); ++i)
		MLOG_D("Range: [%2.4f, %2.4f] => %2.4f\n", min_range[i], max_range[i], map_prob[i]);

	return 0;
}

int MedPredictor::convert_scores_to_prob(const vector<float> &preds, const vector<float> &min_range,
	const vector<float> &max_range, const vector<float> &map_prob, vector<float> &probs) {
	probs.resize(preds.size());

	for (size_t i = 0; i < probs.size(); ++i)
	{
		//search for right range:
		int pos = 0;
		while (pos < map_prob.size() &&
			!((preds[i] > min_range[pos] || pos == map_prob.size() - 1) && (preds[i] <= max_range[pos] || pos == 0)))
			++pos;
		probs[i] = map_prob[pos];
	}

	return 0;
}

template<class T, class L> int MedPredictor::convert_scores_to_prob(const vector<T> &preds, const vector<double> &params, vector<L> &converted) {
	converted.resize((int)preds.size());
	for (size_t i = 0; i < converted.size(); ++i)
	{
		double val = params[0];
		for (size_t k = 1; k < params.size(); ++k)
			val += params[k] * pow(double(preds[i]), double(k));
		val = 1 / (1 + exp(val));//Platt Scale technique for probabilty calibaration
		converted[i] = (L)val;
	}

	return 0;
}
template int MedPredictor::convert_scores_to_prob<double, double>(const vector<double> &preds, const vector<double> &params, vector<double> &converted);
template int MedPredictor::convert_scores_to_prob<double, float>(const vector<double> &preds, const vector<double> &params, vector<float> &converted);
template int MedPredictor::convert_scores_to_prob<float, double>(const vector<float> &preds, const vector<double> &params, vector<double> &converted);
template int MedPredictor::convert_scores_to_prob<float, float>(const vector<float> &preds, const vector<double> &params, vector<float> &converted);

int MedPredictor::learn_prob_calibration(MedMat<float> &x, vector<float> &y,
	int poly_rank, vector<double> &params, int min_bucket_size, float min_score_jump) {
	vector<float> min_range, max_range, map_prob;
	vector<float> preds;
	predict(x, preds);
	learn_prob_calibration(x, y, min_range, max_range, map_prob, min_bucket_size, min_score_jump);

	vector<float> probs;
	convert_scores_to_prob(preds, min_range, max_range, map_prob, probs);
	//probs is the new Y - lets learn A, B:
	MedLinearModel lm(poly_rank); //B is param[0], A is param[1]

	lm.loss_function = [](const vector<double> &prds, const vector<float> &y) {
		double res = 0;
		//L2 on 1 / (1 + exp(A*score + B)) vs Y. prds[i] = A*score+B: 1 / (1 + exp(prds))
		for (size_t i = 0; i < y.size(); ++i)
		{
			double conv_prob = 1 / (1 + exp(prds[i]));
			res += (conv_prob - y[i]) * (conv_prob - y[i]);
		}
		res /= y.size();
		res = sqrt(res);
		return res;
	};
	lm.loss_function_step = [](const vector<double> &prds, const vector<float> &y, const vector<double> &params) {
		double res = 0;
		double reg_coef = 0;
		//L2 on 1 / (1 + exp(A*score + B)) vs Y. prds[i] = A*score+B: 1 / (1 + exp(prds))
		for (size_t i = 0; i < y.size(); ++i)
		{
			double conv_prob = 1 / (1 + exp(prds[i]));
			res += (conv_prob - y[i]) * (conv_prob - y[i]);
		}
		res /= y.size();
		res = sqrt(res);
		//Reg A,B:
		if (reg_coef > 0)
			res += reg_coef * sqrt((params[0] * params[0] + params[1] * params[1]) / 2);
		return res;
	};
	lm.block_num = float(10 * sqrt(poly_rank + 1));
	lm.sample_count = 1000;
	lm.tot_steps = 500000;
	lm.learning_rate = 3 * 1e-1;

	vector<float> poly_preds_params(preds.size() * poly_rank);
	for (size_t j = 0; j < poly_rank; ++j)
		for (size_t i = 0; i < preds.size(); ++i)
			poly_preds_params[i * poly_rank + j] = (float)pow(preds[i], j + 1);

	lm.learn(poly_preds_params, probs, (int)preds.size(), poly_rank);
	vector<float> factors(poly_rank), mean_shifts(poly_rank);
	lm.get_normalization(mean_shifts, factors);

	//put normalizations inside params:
	params.resize(poly_rank + 1);
	params[0] = lm.model_params[0];
	for (size_t i = 1; i < params.size(); ++i) {
		params[i] = lm.model_params[i] / factors[i - 1];
		params[0] -= lm.model_params[i] * mean_shifts[i - 1] / factors[i - 1];
	}

	vector<double> converted((int)preds.size()), prior_score((int)preds.size());
	convert_scores_to_prob(preds, params, converted);
	int tot_pos = 0;
	for (size_t i = 0; i < y.size(); ++i)
		tot_pos += int(y[i] > 0);
	for (size_t i = 0; i < converted.size(); ++i)
		prior_score[i] = double(tot_pos) / y.size();

	double loss_model = _linear_loss_target_rmse(converted, probs);
	double loss_prior = _linear_loss_target_rmse(prior_score, probs);

	MLOG("Platt Scale prior=%2.5f. loss_model=%2.5f, loss_prior=%2.5f\n",
		double(tot_pos) / y.size(), loss_model, loss_prior);

	return 0;
}

// Predicting
int MedPredictor::predict_on_train(MedFeaturesData& ftrs_data, int isplit) {
	MedMat<float> x((int)ftrs_data.get_learning_nrows(isplit), (int)ftrs_data.signals.size());
	vector<float> signal(x.nrows);
	if (!model_features.empty()) {//test names of entered matrix:
		if (model_features.size() != ftrs_data.data.size())
			MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
				(int)model_features.size(), (int)ftrs_data.data.size());
		int feat_num = 0;
		for (auto it = ftrs_data.data.begin(); it != ftrs_data.data.end(); ++it)
			if (norm_feature_name(it->first) != norm_feature_name(model_features[feat_num++]))
				MTHROW_AND_ERR("Learned Features are the same. feat_num=%d. in learning was %s, now recieved %s\n",
					(int)model_features.size(), model_features[feat_num - 1].c_str(), it->first.c_str());
	}


	// Build X
	build_learning_x_mat_for_split(ftrs_data, signal, isplit, x);

	// Predict
	if (predict(x, ftrs_data.split_preds_on_train[isplit]) == -1)
		return -1;

	if (ftrs_data.label_cleaners.size() > 0) {
		for (unsigned int i = 0; i < ftrs_data.split_preds_on_train[isplit].size(); i++) {
			ftrs_data.split_preds_on_train[isplit][i] += ftrs_data.label_cleaners[isplit].mean;
		}
	}

	ftrs_data.n_preds_per_sample = n_preds_per_sample();

	return 0;
}


// Predicting
int MedPredictor::predict(MedFeaturesData& ftrs_data, int isplit) {
	if (!model_features.empty()) {//test names of entered matrix:
		if (model_features.size() != ftrs_data.data.size())
			MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
				(int)model_features.size(), (int)ftrs_data.data.size());
		int feat_num = 0;
		for (auto it = ftrs_data.data.begin(); it != ftrs_data.data.end(); ++it)
			if (norm_feature_name(it->first) != norm_feature_name(model_features[feat_num++]))
				MTHROW_AND_ERR("Learned Features are the same. feat_num=%d. in learning was %s, now recieved %s\n",
					(int)model_features.size(), model_features[feat_num - 1].c_str(), it->first.c_str());
	}
	else if (features_count > 0 && features_count != ftrs_data.data.size())
		MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
			features_count, (int)ftrs_data.data.size());

	// Build X
	MedMat<float> x((int)ftrs_data.get_testing_nrows(isplit), (int)ftrs_data.signals.size());

	vector<float> signal(x.nrows);
	for (int icol = 0; icol < ftrs_data.signals.size(); icol++) {
		string& name = ftrs_data.signals[icol];

		int split_irow = 0;
		for (int irow = 0; irow < ftrs_data.splits.size(); irow++) {
			if (ftrs_data.splits[irow] == isplit)
				signal[split_irow++] = ftrs_data.data[name][irow];
		}

		ftrs_data.cleaners[name][isplit].clear(signal);
		ftrs_data.cleaners[name][isplit].normalize(signal);

		for (unsigned irow = 0; irow < signal.size(); irow++)
			x(irow, icol) = signal[irow];

	}

	x.normalized_flag = 1;

	// Predict
	if (predict(x, ftrs_data.split_preds[isplit]) == -1)
		return -1;

	if (ftrs_data.label_cleaners.size() > 0) {
		for (unsigned int i = 0; i < ftrs_data.split_preds[isplit].size(); i++) {
			ftrs_data.split_preds[isplit][i] += ftrs_data.label_cleaners[isplit].mean;
		}
	}

	ftrs_data.n_preds_per_sample = n_preds_per_sample();

	return 0;

}

int MedPredictor::predict(MedFeaturesData& ftrs_data) {
	if (!model_features.empty()) {//test names of entered matrix:
		if (model_features.size() != ftrs_data.data.size())
			MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
				(int)model_features.size(), (int)ftrs_data.data.size());
		int feat_num = 0;
		for (auto it = ftrs_data.data.begin(); it != ftrs_data.data.end(); ++it)
			if (norm_feature_name(it->first) != norm_feature_name(model_features[feat_num++]))
				MTHROW_AND_ERR("Learned Features are the same. feat_num=%d. in learning was %s, now recieved %s\n",
					(int)model_features.size(), model_features[feat_num - 1].c_str(), it->first.c_str());
	}
	else if (features_count > 0 && features_count != ftrs_data.data.size())
		MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
			features_count, (int)ftrs_data.data.size());

	// Build X
	MedMat<float> x((int)ftrs_data.label.size(), (int)ftrs_data.signals.size());
	MedMat<float> y(x.nrows, 1);

	vector<float> signal(x.nrows);
	for (int icol = 0; icol < x.ncols; icol++) {
		string& name = ftrs_data.signals[icol];

		signal = ftrs_data.data[name];
		ftrs_data.cleaners[name][ftrs_data.nsplits].clear(signal);

		for (unsigned irow = 0; irow < signal.size(); irow++)
			x(icol, irow) = signal[irow];

	}
	x.normalized_flag = 1;

	// Predict
	ftrs_data.n_preds_per_sample = n_preds_per_sample();
	return predict(x, ftrs_data.preds);
}

int MedPredictor::predict(MedFeatures& ftrs_data) {
	if (!model_features.empty()) {//test names of entered matrix:
		if (model_features.size() != ftrs_data.data.size())
			MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
				(int)model_features.size(), (int)ftrs_data.data.size());
		int feat_num = 0;
		for (auto it = ftrs_data.data.begin(); it != ftrs_data.data.end(); ++it)
			if (norm_feature_name(it->first) != norm_feature_name(model_features[feat_num++]))
				MTHROW_AND_ERR("Learned Features are the same. feat_num=%d. in learning was %s, now recieved %s\n",
					(int)model_features.size(), model_features[feat_num - 1].c_str(), it->first.c_str());
	}
	else if (features_count > 0 && features_count != ftrs_data.data.size())
		MTHROW_AND_ERR("Learned Feature model size was %d, request feature size for predict was %d\n",
			features_count, (int)ftrs_data.data.size());

	// Build X
	MedMat<float> x;
	ftrs_data.get_as_matrix(x);

	// Predict
	vector<float> preds;
	if (predict(x, preds) < 0)
		return -1;

	int n = n_preds_per_sample();
	ftrs_data.samples.resize(preds.size() / n);
	for (int i = 0; i < x.nrows; i++) {
		ftrs_data.samples[i].prediction.resize(n);
		for (int j = 0; j < n; j++)
			ftrs_data.samples[i].prediction[j] = preds[i*n + j];
	}

	return 0;
}

void convertXMat(const vector<vector<float>> x, MedMat<float> &xMat) {
	xMat.resize((int)x[0].size(), (int)x.size());
	for (size_t i = 0; i < x.size(); ++i)
	{
		vector<float> xx = x[i];
		for (size_t k = 0; k < xx.size(); ++k)
		{
			xMat((int)k, (int)i) = xx[k];
		}
	}
	xMat.missing_value = MED_MAT_MISSING_VALUE;
}

string medial::models::getParamsInfraModel(void *model) {
	MedPredictor *m = (MedPredictor *)model;
	MedQRFParams pr_qrf;
	MedLightGBMParams pr_lightGBM;
	MedXGBParams pr_xgb;
	map<string, string> empty_m;
	MedSpecificGroupModels *model_specific;
	MedSvm *svm;
	MedLinearModel *lm;
	MedGDLM *gdlm;
	char buff[2000];
	string l1_str;

	switch (m->classifier_type) {
	case MODEL_QRF:
		pr_qrf = ((MedQRF *)model)->params;
		l1_str = to_string(pr_qrf.type);
		if (pr_qrf.type == QRF_BINARY_TREE)
			l1_str = "binary";
		else if (pr_qrf.type == QRF_REGRESSION_TREE)
			l1_str = "regression";
		else if (pr_qrf.type == QRF_CATEGORICAL_CHI2_TREE)
			l1_str = "categorial_chi2";
		else if (pr_qrf.type == QRF_CATEGORICAL_ENTROPY_TREE)
			l1_str = "categorial_entropy";
		snprintf(buff, 2000, "%s: ntrees=%d; maxq=%d; min_node=%d; ntry=%d; spread=%2.3f; type=%s; max_depth=%d",
			predictor_type_to_name[m->classifier_type].c_str(), pr_qrf.ntrees, pr_qrf.maxq, pr_qrf.min_node, pr_qrf.ntry, pr_qrf.spread, l1_str.c_str(), pr_qrf.max_depth);
		break;
	case MODEL_LIGHTGBM:
		pr_lightGBM = ((MedLightGBM *)model)->params;
		snprintf(buff, 2000, "%s: %s",
			predictor_type_to_name[m->classifier_type].c_str(), pr_lightGBM.user_params.c_str());
		break;
	case MODEL_XGB:
		pr_xgb = ((MedXGB *)model)->params;
		snprintf(buff, 2000, "%s: tree_method=%s; booster=%s; objective=%s; eta=%2.3f; alpha=%2.3f; lambda=%2.3f; gamma=%2.3f; max_depth=%d; colsample_bytree=%2.3f; colsample_bylevel=%2.3f; min_child_weight=%d; num_round=%d; subsample=%2.3f",
			predictor_type_to_name[m->classifier_type].c_str(), pr_xgb.tree_method.c_str(), pr_xgb.booster.c_str(), pr_xgb.objective.c_str(), pr_xgb.eta,
			pr_xgb.alpha, pr_xgb.lambda, pr_xgb.gamma, pr_xgb.max_depth, pr_xgb.colsample_bytree, pr_xgb.colsample_bylevel, pr_xgb.min_child_weight, pr_xgb.num_round, pr_xgb.subsample);
		break;
	case MODEL_SPECIFIC_GROUPS_MODELS:
		model_specific = ((MedSpecificGroupModels *)model);
		snprintf(buff, 2000, "%s: model=%s x %d",
			predictor_type_to_name[m->classifier_type].c_str(), predictor_type_to_name[model_specific->get_model(0)->classifier_type].c_str(),
			model_specific->model_cnt());
		break;
	case MODEL_SVM:
		svm = ((MedSvm *)model);
		snprintf(buff, 2000, "%s: kernal_type=%d; C=%2.3f; coef0=%2.3f; degree=%d; gamma=%2.3f; eps=%2.3f",
			predictor_type_to_name[m->classifier_type].c_str(), svm->params.kernel_type
			, svm->params.C, svm->params.coef0, svm->params.degree, svm->params.gamma, svm->params.eps);
		break;
	case MODEL_LINEAR_SGD:
		lm = ((MedLinearModel *)model);
		l1_str = lm->norm_l1 ? "(L1)" : "(L2)";
		snprintf(buff, 2000, "%s: name=%s; num_params=%d; block_num=%2.3f%s; learning_rate=%2.3f; sample_count=%d; tot_steps=%d",
			predictor_type_to_name[m->classifier_type].c_str(), lm->model_name.c_str(),
			(int)lm->model_params.size(), lm->block_num, l1_str.c_str(),
			lm->learning_rate, lm->sample_count, lm->tot_steps);
		break;
	case MODEL_GD_LINEAR:
		gdlm = (MedGDLM*)model;
		snprintf(buff, 2000, "%s: method=%s; batch-size=%d; l_lasso=%2.3f; l_ridge=%2.3f; rate=%2.3f; rate_decay=%2.3f; momentum=%2.3f; stop_at_err=%2.3f; max_iters=%d",
			predictor_type_to_name[m->classifier_type].c_str(), gdlm->params.method.c_str(), gdlm->params.batch_size,
			gdlm->params.l_lasso, gdlm->params.l_ridge, gdlm->params.rate, gdlm->params.rate_decay,
			gdlm->params.momentum, gdlm->params.stop_at_err, gdlm->params.max_iter);
		break;
	default:
		throw invalid_argument("Unsupported Type init for model " + to_string(m->classifier_type));
	}

	return string(buff);
}

void *medial::models::copyInfraModel(void *model, bool delete_old) {
	MedPredictor *m = (MedPredictor *)model;
	MedQRFParams pr_qrf;
	MedLightGBMParams pr_lightGBM;
	MedXGBParams pr_xgb;
	map<string, string> empty_m;
	MedSpecificGroupModels *model_specific;
	MedSvm *svm;
	MedLinearModel *lm;
	MedGDLM *gdlm;
	void *newM;

	switch (m->classifier_type) {
	case MODEL_QRF:
		pr_qrf = MedQRFParams(((MedQRF *)model)->params);
		newM = new MedQRF(pr_qrf);
		((MedQRF *)newM)->qf = QRF_Forest(); //Erase forest
		if (delete_old)
			delete ((MedQRF *)model);
		break;
	case MODEL_LIGHTGBM:
		pr_lightGBM = MedLightGBMParams(((MedLightGBM *)model)->params);
		if (delete_old)
			delete ((MedLightGBM *)model);
		newM = new MedLightGBM;
		((MedLightGBM *)newM)->params = pr_lightGBM;
		((MedLightGBM *)newM)->init(empty_m);
		break;
	case MODEL_XGB:
		pr_xgb = MedXGBParams(((MedXGB *)model)->params);
		if (delete_old)
			delete ((MedXGB *)model);
		newM = new MedXGB;
		((MedXGB *)newM)->params = pr_xgb;
		break;
	case MODEL_SPECIFIC_GROUPS_MODELS:
		model_specific = ((MedSpecificGroupModels *)model);
		newM = model_specific->clone();
		for (size_t i = 0; i < model_specific->model_cnt(); ++i)
			(*((MedQRF *)model_specific->get_model((int)i))).qf = QRF_Forest(); //to release memory

		break;
	case MODEL_SVM:
		svm = ((MedSvm *)model);
		newM = new MedSvm(svm->params);
		if (delete_old)
			delete (MedSvm *)model;
		break;
	case MODEL_LINEAR_SGD:
		lm = ((MedLinearModel *)model);
		newM = (MedLinearModel *)lm->clone();
		if (delete_old)
			delete (MedLinearModel *)model;
		break;
	case MODEL_GD_LINEAR:
		gdlm = (MedGDLM*)model;
		newM = new MedGDLM(gdlm->params);
		if (delete_old)
			delete (MedGDLM *)model;
		break;
	default:
		throw invalid_argument("Unsupported Type init for model " + to_string(m->classifier_type) + " (copy)");
	}

	return newM;
}

void medial::models::initInfraModel(void *&model) {
	void *newM = copyInfraModel(model, true);

	model = newM;
}

void medial::models::learnInfraModel(void *model, const vector<vector<float>> &xTrain, vector<float> &y, vector<float> &weights) {
	MedMat<float> xTrain_m;
	convertXMat(xTrain, xTrain_m);
	MedPredictor *m = (MedPredictor *)model;
	if (m->normalize_for_learn)
		xTrain_m.normalize();
	m->learn(xTrain_m, y, weights);
}

vector<float> medial::models::predictInfraModel(void *model, const vector<vector<float>> &xTest) {
	MedMat<float> xTest_m;
	convertXMat(xTest, xTest_m);
	MedPredictor *m = (MedPredictor *)model;
	if (m->normalize_for_predict)
		xTest_m.normalize();
	vector<float> preds;
	m->predict(xTest_m, preds);
	return preds;
}