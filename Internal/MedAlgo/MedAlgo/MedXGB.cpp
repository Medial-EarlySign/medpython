#if 1
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE

#include "MedAlgo.h"
#include "MedXGB.h"
#include <boost/lexical_cast.hpp>
#include <xgboost/learner.h>
#include <dmlc/timer.h>
#include "xgboost/c_api.h"

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

using namespace xgboost;
using namespace std;
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
int MedXGB::Predict(float *x, float *&preds, int nsamples, int nftrs) {
	DMatrix *out;
	if (XGDMatrixCreateFromMat(x, nsamples, nftrs, params.missing_value, (DMatrixHandle*)&out) == -1) {
		MERR("failed to XGDMatrixCreateFromMat");
		throw std::exception();
	}
	std::unique_ptr<DMatrix> dtest(out);

	vector<float> out_preds;
	my_learner->Predict(dtest.get(), false, &out_preds, 0);

	for (int i = 0; i < out_preds.size(); i++)
		preds[i] = out_preds[i];

	return 0;
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
		dvalidate->info().labels.push_back(y[i]);

	return 0;
}
int MedXGB::Learn(float *x, float *y, float *w, int nsamples, int nftrs) {
	std::vector<std::pair<std::string, std::string> > cfg;
	cfg.push_back(std::make_pair("seed", boost::lexical_cast<std::string>(params.seed)));
	cfg.push_back(std::make_pair("booster", params.booster));
	cfg.push_back(std::make_pair("objective", params.objective));
	cfg.push_back(std::make_pair("eta", boost::lexical_cast<std::string>(params.eta)));
	cfg.push_back(std::make_pair("gamma", boost::lexical_cast<std::string>(params.gamma)));
	cfg.push_back(std::make_pair("min_child_weight", boost::lexical_cast<std::string>(params.min_child_weight)));
	cfg.push_back(std::make_pair("max_depth", boost::lexical_cast<std::string>(params.max_depth)));
	cfg.push_back(std::make_pair("silent", boost::lexical_cast<std::string>(params.silent)));
	cfg.push_back(std::make_pair("eval_metric", params.eval_metric));
	cfg.push_back(std::make_pair("colsample_bytree", boost::lexical_cast<std::string>(params.colsample_bytree))); 
	cfg.push_back(std::make_pair("colsample_bylevel", boost::lexical_cast<std::string>(params.colsample_bylevel)));
	cfg.push_back(std::make_pair("subsample", boost::lexical_cast<std::string>(params.subsample))); 
	cfg.push_back(std::make_pair("scale_pos_weight", boost::lexical_cast<std::string>(params.scale_pos_weight)));
	cfg.push_back(std::make_pair("lambda", boost::lexical_cast<std::string>(params.lambda)));
	cfg.push_back(std::make_pair("alpha", boost::lexical_cast<std::string>(params.alpha)));
	cfg.push_back(std::make_pair("tree_method", boost::lexical_cast<std::string>(params.tree_method)));
	//cfg.push_back(std::make_pair("num_class", boost::lexical_cast<std::string>(params.num_class)));

	DMatrix *out;
	if (XGDMatrixCreateFromMat(x, nsamples, nftrs, params.missing_value, (DMatrixHandle*)&out) == -1) {
		MERR("failed to XGDMatrixCreateFromMat");
		throw std::exception();
	}

	std::unique_ptr<DMatrix> dtrain(out);

	for (int i = 0; i < nsamples; i++) {
		dtrain.get()->info().labels.push_back(y[i]);
		dtrain.get()->info().weights.push_back(w[i]);
	}

	std::vector<std::unique_ptr<DMatrix> > deval;
	std::vector<DMatrix*> cache_mats, eval_datasets;
	cache_mats.push_back(dtrain.get());
	if (dvalidate != NULL)
		cache_mats.push_back(dvalidate);

	std::vector<std::string> eval_data_names;
	eval_datasets.push_back(dtrain.get());
	eval_data_names.push_back(std::string("train"));
	if (dvalidate != NULL) {
		eval_datasets.push_back(dvalidate);
		eval_data_names.push_back(std::string("validate"));
	}

	// initialize the learner.
	Learner* learner = Learner::Create(cache_mats);
	int version = 0;

	learner->Configure(cfg);
	learner->InitModel();

	// start training.
	const double start = dmlc::GetTime();
	for (int i = version / 2; i < params.num_round; ++i) {
		double elapsed = dmlc::GetTime() - start;
		if (version % 2 == 0) {
			if (params.silent == 0)
				MLOG("boosting round %d, %d sec elapsed", i, elapsed);
			learner->UpdateOneIter(i, dtrain.get());
			version += 1;
		}
		std::string res = learner->EvalOneIter(i, eval_datasets, eval_data_names);
		if (params.silent == 0) 
			MLOG("%s", res.c_str());
		version += 1;
	}

	double elapsed = dmlc::GetTime() - start;
	if (params.silent == 0)
		MLOG("update end, %d sec in all", elapsed);

	if (this->my_learner != NULL)
		delete this->my_learner;
	this->my_learner = learner;
	return 0;
}

typedef rabit::utils::MemoryFixSizeBuffer MemoryFixSizeBuffer;
typedef rabit::utils::MemoryBufferStream MemoryBufferStream;

size_t MedXGB::get_size() {

	std::string raw_str;
	raw_str.resize(0);

	MemoryBufferStream fo(&raw_str);
	my_learner->Save(&fo);
	return raw_str.size() + sizeof(int);

}

void MedXGB::print(FILE *fp, const string& prefix) {
	fprintf(fp, "%s: MedXGB ()\n", prefix.c_str());
}

size_t MedXGB::serialize(unsigned char *blob) {
	std::string raw_str;
	raw_str.resize(0);

	MemoryBufferStream fo(&raw_str);
	my_learner->Save(&fo);

	*((int*)blob) = (int)raw_str.size();
	memcpy(blob + sizeof(int), raw_str.c_str(), raw_str.size());
	return raw_str.size() + sizeof(int);
}

size_t MedXGB::deserialize(unsigned char *blob) {
	int size = *((int*)blob);
	MemoryFixSizeBuffer fs((blob + sizeof(int)), size);
	my_learner->Load(&fs);
	return size + sizeof(int);
}

void MedXGB::init_defaults()
{
	classifier_type = MODEL_XGB;
	std::vector<DMatrix*> mats;
	this->my_learner = Learner::Create(mats);

	transpose_for_learn = false;
	transpose_for_predict = false;
	normalize_for_learn = false;
	normalize_for_predict = false;
	normalize_y_for_learn = false;

}

int MedXGB::init(map<string, string>& mapper) {
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "booster") params.booster = entry.second;
		else if (field == "objective") params.objective = entry.second;
		else if (field == "eval_metric") params.eval_metric = entry.second;
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
		else MLOG("Unknonw parameter \'%s\' for XGB\n", field.c_str());
	}

	return 0;
}

#endif
