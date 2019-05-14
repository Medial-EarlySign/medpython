#if 1
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE

#include "MedAlgo.h"
#include "MedXGB.h"
#include <boost/lexical_cast.hpp>


#include <dmlc/timer.h>

#include <omp.h>

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

using namespace xgboost;
using namespace std;

MedXGB::~MedXGB() {
	if (my_learner != NULL) {
		XGBoosterFree(my_learner);
		my_learner = NULL;
	}
	for (size_t i = 0; i < learner_per_thread.size(); ++i) {
		if (learner_per_thread[i] != NULL) {
			XGBoosterFree(learner_per_thread[i]);
			learner_per_thread[i] = NULL;
		}
	}
}

int MedXGB::n_preds_per_sample() const {
	if (params.objective == "multi:softprob")
		return params.num_class;
	return 1;
}

/*
#if defined(_MSC_VER) || defined(_WIN32)
__declspec(dllexport) int XGDMatrixCreateFromMat(const float *data,
	xgboost::bst_ulong nrow,
	xgboost::bst_ulong ncol,
	float  missing,
	DMatrixHandle *out) {
	MERR("xgboost not supported in windows as MSVC can not build it");
	throw std::exception();
}

Learner* Learner::Create(const std::vector<DMatrix*>& cache_data) {
	MERR("xgboost not supported in windows as MSVC can not build it");
	throw std::exception();
}
#endif
*/
int MedXGB::Predict(float *x, float *&preds, int nsamples, int nftrs) const {
	DMatrixHandle h_test;
	if (XGDMatrixCreateFromMat(x, nsamples, nftrs, params.missing_value, &h_test) == -1)
		MTHROW_AND_ERR("failed to XGDMatrixCreateFromMat");

	xgboost::bst_ulong out_len;
	const float *out_preds;
	XGBoosterPredict(my_learner, h_test, 0, 0, &out_len, &out_preds);

	int64_t len_res = nsamples * n_preds_per_sample();
	if (preds == NULL) preds = new float[len_res];
	for (int i = 0; i < out_len; i++)
		preds[i] = out_preds[i];

	XGDMatrixFree(h_test);
	return 0;
}


void MedXGB::calc_feature_contribs(MedMat<float> &mat_x, MedMat<float> &mat_contribs) {
	int nsamples, nftrs;
	vector<float> w;
	prepare_x_mat(mat_x, w, nsamples, nftrs, transpose_for_predict);

	mat_contribs.resize(nsamples, nftrs + 1);
	// copy metadata
	mat_contribs.signals.insert(mat_contribs.signals.end(), mat_x.signals.begin(), mat_x.signals.end());
	mat_contribs.signals.push_back("b0");
	mat_contribs.recordsMetadata.insert(mat_contribs.recordsMetadata.end(), mat_x.recordsMetadata.begin(), mat_x.recordsMetadata.end());

	DMatrixHandle h_test;
	if (XGDMatrixCreateFromMat(mat_x.data_ptr(), nsamples, nftrs, params.missing_value, &h_test) == -1)
		MTHROW_AND_ERR("failed to XGDMatrixCreateFromMat");

	xgboost::bst_ulong out_len;
	const float *out_preds;
	const int PRED_CONTRIBS = 4; //, APPROX_CONTRIBS = 8;  , INTERACTION_SHAP = 16;
	int flags = feat_contrib_flags;
	//if (flags == 0) // default value is now not APPROX_CONTRIBS since nan bug was solved.
	//	flags = APPROX_CONTRIBS;

	flags |= PRED_CONTRIBS;

	XGBoosterPredict(my_learner, h_test, flags, 0, &out_len, &out_preds);
	for (int i = 0; i < nsamples; i++) {
		for (int j = 0; j < nftrs; j++) {
			float v = out_preds[i*(nftrs + 1) + j];
			if (isnan(v))
				MTHROW_AND_ERR("got nan in (%d,%d)\n", i, j);
			mat_contribs.set(i, j) = v;
		}
		mat_contribs.set(i, nftrs) = out_preds[i*(nftrs + 1) + nftrs];
	}
	XGDMatrixFree(h_test);
}

int MedXGB::Learn(float *x, float *y, int nsamples, int nftrs) {
	vector<float> w;
	for (int i = 0; i < nsamples; i++)
		w.push_back(1.0);
	return Learn(x, y, &w[0], nsamples, nftrs);
}

int MedXGB::validate_me_while_learning(float *x, float *y, int nsamples, int nftrs) {
	DMatrix *out;
	if (XGDMatrixCreateFromMat(x, nsamples, nftrs, params.missing_value, (DMatrixHandle*)&out) == -1) {
		MERR("failed to XGDMatrixCreateFromMat");
		throw std::exception();
	}

	if (this->dvalidate != NULL)
		delete this->dvalidate;

	dvalidate = out;

	for (int i = 0; i < nsamples; i++)
		dvalidate->Info().labels_.HostVector().push_back(y[i]);

	return 0;
}

void MedXGB::prepare_mat_handle(float *x, float *y, const float *w, int nsamples, int nftrs, DMatrixHandle &matrix_handle)
{
	if (XGDMatrixCreateFromMat(x, nsamples, nftrs, params.missing_value, &matrix_handle) != 0)
		MTHROW_AND_ERR("failed to XGDMatrixCreateFromMat");

	if (XGDMatrixSetFloatInfo(matrix_handle, "label", y, nsamples) != 0)
		MTHROW_AND_ERR("failed XGDMatrixSetFloatInfo label");

	if (w != NULL) {
		if (XGDMatrixSetFloatInfo(matrix_handle, "weight", w, nsamples) != 0)
			MTHROW_AND_ERR("failed XGDMatrixSetFloatInfo weight");
	}
}

int MedXGB::Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
	DMatrixHandle matrix_handles[2];

	if ((params.verbose_eval > 0) & (params.validate_frac > 0))
	{
		//divide to x_train and x_test
		if ((params.validate_frac < 0) || (params.validate_frac > 1))
			MTHROW_AND_ERR("Validation fraction %f is invalid \n", params.validate_frac);

		int nsamples_test = int(params.validate_frac*nsamples);
		int nsamples_train = nsamples - nsamples_test;

		prepare_mat_handle(x, y, w, nsamples_train, nftrs, matrix_handles[0]);
		prepare_mat_handle(x + nsamples_train * nftrs, y + nsamples_train, (w == NULL) ? NULL : w + nsamples_train, nsamples_test, nftrs, matrix_handles[1]);
	}
	else {
		prepare_mat_handle(x, y, w, nsamples, nftrs, matrix_handles[0]);
	}

	BoosterHandle h_booster;
	if (XGBoosterCreate(&matrix_handles[0], 1, &h_booster) != 0)
		MTHROW_AND_ERR("failed XGBoosterCreate weight");

	XGBoosterSetParam(h_booster, "seed", boost::lexical_cast<std::string>(params.seed).c_str());
	XGBoosterSetParam(h_booster, "booster", params.booster.c_str());
	XGBoosterSetParam(h_booster, "objective", params.objective.c_str());
	XGBoosterSetParam(h_booster, "eta", boost::lexical_cast<std::string>(params.eta).c_str());
	XGBoosterSetParam(h_booster, "gamma", boost::lexical_cast<std::string>(params.gamma).c_str());
	XGBoosterSetParam(h_booster, "min_child_weight", boost::lexical_cast<std::string>(params.min_child_weight).c_str());
	XGBoosterSetParam(h_booster, "max_depth", boost::lexical_cast<std::string>(params.max_depth).c_str());
	XGBoosterSetParam(h_booster, "silent", boost::lexical_cast<std::string>(params.silent).c_str());
	XGBoosterSetParam(h_booster, "colsample_bytree", boost::lexical_cast<std::string>(params.colsample_bytree).c_str());
	XGBoosterSetParam(h_booster, "colsample_bylevel", boost::lexical_cast<std::string>(params.colsample_bylevel).c_str());
	XGBoosterSetParam(h_booster, "subsample", boost::lexical_cast<std::string>(params.subsample).c_str());
	XGBoosterSetParam(h_booster, "scale_pos_weight", boost::lexical_cast<std::string>(params.scale_pos_weight).c_str());
	XGBoosterSetParam(h_booster, "lambda", boost::lexical_cast<std::string>(params.lambda).c_str());
	XGBoosterSetParam(h_booster, "alpha", boost::lexical_cast<std::string>(params.alpha).c_str());
	XGBoosterSetParam(h_booster, "num_class", boost::lexical_cast<std::string>(params.num_class).c_str());
	XGBoosterSetParam(h_booster, "tree_method", boost::lexical_cast<std::string>(params.tree_method).c_str());
	XGBoosterSetParam(h_booster, "verbose_eval", boost::lexical_cast<std::string>(params.verbose_eval).c_str());

	for (auto it : params.eval_metric)
	{
		XGBoosterSetParam(h_booster, "eval_metric", it.c_str());
	}

	string split_penalties_s;
	translate_split_penalties(split_penalties_s);
	XGBoosterSetParam(h_booster, "split_penalties_s", boost::lexical_cast<std::string>(split_penalties_s).c_str());

	const double start = dmlc::GetTime();
	const char *evnames[2] = { "train", "test" };
	const char *out_result;

#pragma omp critical
	XGBoosterUpdateOneIter(h_booster, 0, matrix_handles[0]);
	for (int iter = 1; iter < params.num_round; iter++)
	{
		XGBoosterUpdateOneIter(h_booster, iter, matrix_handles[0]);
		if (params.verbose_eval > 0 && iter % params.verbose_eval == 0)
		{
			if (params.validate_frac > 0) { XGBoosterEvalOneIter(h_booster, iter, matrix_handles, evnames, 2, &out_result); }
			else { XGBoosterEvalOneIter(h_booster, iter, matrix_handles, evnames, 1, &out_result); }

			MLOG("Performance: %s \n", out_result);
		}
	}

	double elapsed = dmlc::GetTime() - start;
	if (params.silent == 0)
		MLOG("update end, %d sec overall", elapsed);

	if (this->my_learner != NULL)
		XGBoosterFree(this->my_learner);
	this->my_learner = h_booster;
	_mark_learn_done = true;

	XGDMatrixFree(matrix_handles[0]);
	return 0;
}

void MedXGB::translate_split_penalties(string& split_penalties_s) {

	vector<string> elems;
	boost::split(elems, params.split_penalties, boost::is_any_of(",:"));
	if (elems.size() < 2)
		return;

	vector<string> out_elems;
	for (unsigned int i = 0; i < elems.size(); i += 2) {
		int index = find_in_feature_names(model_features, elems[i]);
		out_elems.push_back(to_string(index) + ":" + elems[i + 1]);
	}

	split_penalties_s = boost::join(out_elems, ",");
}

typedef rabit::utils::MemoryFixSizeBuffer MemoryFixSizeBuffer;
typedef rabit::utils::MemoryBufferStream MemoryBufferStream;

void MedXGB::print(FILE *fp, const string& prefix, int level) const {

	if (level==0)
		fprintf(fp, "%s: MedXGB ()\n", prefix.c_str());
	else {
		xgboost::bst_ulong num_trees;
		const char **out_models;
		XGBoosterDumpModel(my_learner, "", 0, &num_trees, &out_models);
		for (int i = 0; i < num_trees; i++)
			fprintf(fp, "%s xgboost tree %d : %s\n", prefix.c_str(), i, out_models[i]);
	}

}

bool split_token(const string &str, const string &search, bool first
	, string &result) {
	if (str.find(search) == string::npos)
		return false;
	result = first ? str.substr(0, str.find(search)) : str.substr(str.find(search) + search.size());
	return true;
}

/*
Importance type can be defined as:
'weight' - the number of times a feature is used to split the data across all trees.
'gain' - the average gain of the feature when it is used in trees
'cover' - the average coverage of the feature when it is used in trees
'gain_total' - sum of gain the of the feature when it is used in trees (not normalized by number)
'shap' - mean absolute value of shap values
*/
void MedXGB::calc_feature_importance(vector<float> &features_importance_scores,
	const string &general_params, const MedFeatures *features) {
	if (!_mark_learn_done)
		MTHROW_AND_ERR("ERROR:: Requested calc_feature_importance before running learn\n");
	map<string, string> params;

	unordered_set<string> local_types = { "weight", "gain","cover","gain_total" };
	unordered_set<string> legal_types = { "weight", "gain","cover","gain_total", "shap" };

	MedSerialize::initialization_text_to_map(general_params, params);
	string importance_type = "gain";
	for (auto it = params.begin(); it != params.end(); ++it)
		if (it->first == "importance_type")
			importance_type = it->second;
		else
			MTHROW_AND_ERR("Unsupported calc_feature_importance param \"%s\"\n", it->first.c_str());

	if (legal_types.find(importance_type) == legal_types.end())
		MTHROW_AND_ERR("Ilegal importance_type value \"%s\" "
			"- should by one of [weight, gain, cover, gain_total]\n", importance_type.c_str());

	features_importance_scores.resize(model_features.empty() ? features_count : (int)model_features.size());

	if (local_types.count(importance_type) > 0)
		calc_feature_importance_local(features_importance_scores, importance_type);


	if (importance_type == "shap")
		calc_feature_importance_shap(features_importance_scores, importance_type, features);
	
}

void MedXGB::calc_feature_importance_local(vector<float> &features_importance_scores, string &importance_type)
{
	bool do_average = false;
	vector<float> fCnt;
	const char** out_models;
	int with_stats = importance_type != "weight"; //if weight than 0

	xgboost::bst_ulong num_trees;
	XGBoosterDumpModel(my_learner, "", with_stats, &num_trees, &out_models);
	vector<string> arr;
	string mid_token, fids;
	int fid;
	float g;

	if ((importance_type == "gain") || (importance_type == "cover"))
		do_average = true;

	if (importance_type != "weight")
		fCnt.resize((int)features_importance_scores.size());

	string search_str = ((importance_type == "gain_total") ? "gain" : importance_type) + "=";

	for (xgboost::bst_ulong tree_num = 0; tree_num < num_trees; tree_num++)
	{
		vector<string> lines;
		string tree = out_models[tree_num];
		boost::split(lines, tree, boost::is_any_of("\n"));
		for (string &line : lines)
		{
			arr.clear();
			boost::split(arr, line, boost::is_any_of("["));
			if (arr.size() == 1)
				continue;
			if (!split_token(arr[1], "]", true, mid_token))
				MTHROW_AND_ERR("format error in line \"%s\"\n", line.c_str());

			if (!split_token(mid_token, "<", true, fids))
				MTHROW_AND_ERR("format error in line \"%s\"\n", line.c_str());
			fid = stoi(fids.substr(1));

			if (importance_type == "weight")
				features_importance_scores[fid] += 1;
			else {
				fCnt[fid] += 1;
				split_token(arr[1], "]", false, mid_token);
				if (!split_token(mid_token, search_str, false, fids))
					MTHROW_AND_ERR("format error in line \"%s\"\n", line.c_str());
				if (!split_token(fids, ",", true, mid_token))
					MTHROW_AND_ERR("format error in line \"%s\"\n", line.c_str());

				g = stof(mid_token);
				features_importance_scores[fid] += g;
			}
		}
	}

	if (do_average)
		for (size_t i = 0; i < features_importance_scores.size(); ++i)
			if (fCnt[i])
			{
				features_importance_scores[i] /= fCnt[i];
			}
}

//void MedXGB::calc_feature_importance_shap(vector<float> &features_importance_scores, string &importance_type, const MedFeatures *features)
//{
//	MedMat<float> feat_mat,contribs_mat; 
//	if (features == NULL)
//		MTHROW_AND_ERR("SHAP values feature importance requires features \n");
//	
//	features->get_as_matrix(feat_mat);
//	calc_feature_contribs(feat_mat, contribs_mat);
//	for (int j = 0; j < contribs_mat.ncols; ++j)
//	{
//		float col_sum = 0;
//		
//		for (int i = 0; i < contribs_mat.nrows; ++i)
//		{
//			col_sum += abs(contribs_mat.get(i, j));
//		}
//		features_importance_scores[j] = col_sum/(float)contribs_mat.nrows;
//	}
//}

void MedXGB::init_defaults()
{
	classifier_type = MODEL_XGB;
	this->my_learner = NULL;

	transpose_for_learn = false;
	transpose_for_predict = false;
	normalize_for_learn = false;
	normalize_for_predict = false;
	normalize_y_for_learn = false;

	_mark_learn_done = false;
	prepared_single = false;
}

int MedXGB::set_params(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [MedXGB::init]
		if (field == "booster") params.booster = entry.second;
		else if (field == "objective") params.objective = entry.second;
		else if (field == "eval_metric") split(params.eval_metric, entry.second, boost::is_any_of(","));
		else if (field == "eta") params.eta = stof(entry.second);
		else if (field == "gamma") params.gamma = stof(entry.second);
		else if (field == "min_child_weight") params.min_child_weight = stoi(entry.second);
		else if (field == "max_depth") params.max_depth = stoi(entry.second);
		else if (field == "num_round") params.num_round = stoi(entry.second);
		else if (field == "silent") params.silent = stoi(entry.second);
		else if (field == "num_class") params.num_class = stoi(entry.second);
		else if (field == "missing_value") params.missing_value = stof(entry.second);
		else if (field == "colsample_bytree") params.colsample_bytree = stof(entry.second);
		else if (field == "colsample_bylevel") params.colsample_bylevel = stof(entry.second);
		else if (field == "subsample") params.subsample = stof(entry.second);
		else if (field == "scale_pos_weight") params.scale_pos_weight = stof(entry.second);
		else if (field == "tree_method") params.tree_method = entry.second;
		else if (field == "lambda") params.lambda = stof(entry.second);
		else if (field == "alpha") params.alpha = stof(entry.second);
		else if (field == "seed") params.seed = stoi(entry.second);
		else if (field == "verbose_eval") params.verbose_eval = stoi(entry.second);
		else if (field == "split_penalties") params.split_penalties = entry.second;
		else if (field == "validate_frac") params.validate_frac = stof(entry.second);

		else MLOG("Unknonw parameter \'%s\' for XGB\n", field.c_str());
		//! [MedXGB::init]
	}

	return 0;
}


void MedXGB::prepare_predict_single() {
	if (prepared_single)
		return;
	prepared_single = true;
	int nftrs = features_count;
	if (nftrs == 0)
		nftrs = (int)model_features.size();
	XGBBooster *xgb_mdl = static_cast<XGBBooster*>(my_learner);
	xgb_mdl->LazyInit();
	int N_tot_threads = omp_get_max_threads();

	const char* out_dptr;
	xgboost::bst_ulong len;
	if (XGBoosterGetModelRaw(my_learner, &len, &out_dptr) != 0)
		throw runtime_error("failed XGBoosterGetModelRaw\n");

	learner_per_thread.resize(N_tot_threads);
	for (size_t i = 0; i < N_tot_threads; ++i)
	{
		BoosterHandle temp_handler;
		DMatrixHandle h_train_empty[1];
		if (XGBoosterCreate(h_train_empty, 0, &temp_handler) != 0)
			throw runtime_error("failed XGBoosterCreate\n");
		if (XGBoosterLoadModelFromBuffer(temp_handler, out_dptr, len) != 0)
			throw runtime_error("failed XGBoosterLoadModelFromBuffer\n");
		static_cast<XGBBooster*>(temp_handler)->LazyInit();
		learner_per_thread[i] = temp_handler;
		//learner_per_thread[i].LoadModel(dmlc::Stream(out_dptr));
	}
}



void MedXGB::predict_single(const vector<float> &x, vector<float> &preds) const {
	int n_ftrs = (int)x.size();
	std::unique_ptr<data::SimpleCSRSource> p_mat(new data::SimpleCSRSource());
	data::SimpleCSRSource &mat = *p_mat; //copy memory to alter values (const object)
	auto &offset_vec = mat.page_.offset.HostVector();
	vector<xgboost::Entry> &vals = mat.page_.data.HostVector();
	offset_vec.resize(2);
	mat.info.num_row_ = 1;
	mat.info.num_col_ = n_ftrs;

	// count elements for sizing data
	xgboost::bst_ulong nelem = 0;
	for (xgboost::bst_ulong j = 0; j < n_ftrs; ++j)
		if (x[j] != params.missing_value)
			++nelem;
	offset_vec[1] = offset_vec[0] + nelem;
	vals.resize(vals.size() + offset_vec.back());
	xgboost::bst_ulong matj = 0;
	for (xgboost::bst_ulong j = 0; j < n_ftrs; ++j) {
		if (x[j] != params.missing_value) {
			vals[offset_vec[0] + matj] = xgboost::Entry(j, x[j]);
			++matj;
		}
	}
	mat.info.num_nonzero_ = vals.size();

	DMatrix *mat_gen = DMatrix::Create(move(p_mat));

	//int len_res = n_preds_per_sample();
	xgboost::HostDeviceVector<float> wrapper;
	int n_th = omp_get_thread_num();
	//xgboost::Learner *xgb_mdl = static_cast<XGBBooster*>(my_learner)->learner();
	xgboost::Learner *xgb_mdl = static_cast<XGBBooster*>(learner_per_thread[n_th])->learner();

	//for each thread learner
	xgb_mdl->Predict(mat_gen, false, &wrapper, 0, false, false, false, false);

	preds = move(wrapper.HostVector());

	delete mat_gen;
}

#endif
