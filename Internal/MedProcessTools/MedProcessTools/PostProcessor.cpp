#include "PostProcessor.h"
#include "Calibration.h"
#include <boost/algorithm/string.hpp>
#include "ExplainWrapper.h"
#include "AggregatePredsPostProcessor.h"
#include "ProbAdjustPostProcessor.h"
#include "FairnessPostProcessor.h"

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL LOG_DEF_LEVEL

PostProcessorTypes post_processor_name_to_type(const string& post_processor) {
	string lower_p = boost::to_lower_copy(post_processor);
	if (lower_p == "multi")
		return FTR_POSTPROCESS_MULTI;
	else if (lower_p == "calibrator")
		return FTR_POSTPROCESS_CALIBRATOR;
	else if (lower_p == "tree_shap")
		return FTR_POSTPROCESS_TREE_SHAP;
	else if (lower_p == "shapley")
		return FTR_POSTPROCESS_SHAPLEY;
	else if (lower_p == "missing_shap")
		return FTR_POSTPROCESS_MISSING_SHAP;
	else if (lower_p == "lime_shap")
		return FTR_POSTPROCESS_LIME_SHAP;
	else if (lower_p == "linear")
		return FTR_POSTPROCESS_LINEAR;
	else if (lower_p == "knn")
		return FTR_POSTPROCESS_KNN_EXPLAIN;
	else if (lower_p == "iterative_set")
		return FTR_POSTPROCESS_ITERATIVE_SET;
	else if (lower_p == "aggregate_preds")
		return FTR_POSTPROCESS_AGGREGATE_PREDS;
	else if (lower_p == "adjust_probs")
		return FTR_POSTPROCESS_ADJUST;
	else if (lower_p == "fairness_adjust")
		return FTR_POSTPROCESS_FAIRNESS;
	else
		MTHROW_AND_ERR("Unsupported PostProcessor %s\n", post_processor.c_str());
}

PostProcessor *PostProcessor::make_processor(const string &processor_name, const string &params) {
	return make_processor(post_processor_name_to_type(processor_name), params);
}

void PostProcessor::dprint(const string &pref) const {
	MLOG("%s :: PP type %d(%s)\n", pref.c_str(), processor_type, my_class_name().c_str());
}

PostProcessor *PostProcessor::make_processor(PostProcessorTypes type, const string &params) {
	PostProcessor *prc;
	if (type == FTR_POSTPROCESS_MULTI)
		prc = new MultiPostProcessor;
	else if (type == FTR_POSTPROCESS_CALIBRATOR)
		prc = new Calibrator;
	else if (type == FTR_POSTPROCESS_TREE_SHAP)
		prc = new TreeExplainer;
	else if (type == FTR_POSTPROCESS_SHAPLEY)
		prc = new ShapleyExplainer;
	else if (type == FTR_POSTPROCESS_MISSING_SHAP)
		prc = new MissingShapExplainer;
	else if (type == FTR_POSTPROCESS_LIME_SHAP)
		prc = new LimeExplainer;
	else if (type == FTR_POSTPROCESS_LINEAR)
		prc = new LinearExplainer;
	else if (type == FTR_POSTPROCESS_KNN_EXPLAIN)
		prc = new KNN_Explainer;
	else if (type == FTR_POSTPROCESS_ITERATIVE_SET)
		prc = new IterativeSetExplainer;
	else if (type == FTR_POSTPROCESS_AGGREGATE_PREDS)
		prc = new AggregatePredsPostProcessor;
	else if (type == FTR_POSTPROCESS_ADJUST)
		prc = new ProbAdjustPostProcessor;
	else if (type == FTR_POSTPROCESS_FAIRNESS)
		prc = new FairnessPostProcessor;
	else
		MTHROW_AND_ERR("Unsupported PostProcessor %d\n", type);

	prc->init_from_string(params);

	return prc;
}

// Create processor from params string (type must be given within string)
//.......................................................................................
PostProcessor *PostProcessor::create_processor(string &params)
{
	string pp_type;
	get_single_val_from_init_string(params, "pp_type", pp_type);
	return (make_processor(pp_type, params));
}


void PostProcessor::Learn(const MedFeatures &matrix) {
	MTHROW_AND_ERR("Learn Not implemented in class %s\n", my_class_name().c_str());
}
void PostProcessor::Apply(MedFeatures &matrix) {
	MTHROW_AND_ERR("Apply Not implemented in class %s\n", my_class_name().c_str());
}

void *PostProcessor::new_polymorphic(string dname)
{
	CONDITIONAL_NEW_CLASS(dname, MultiPostProcessor);
	CONDITIONAL_NEW_CLASS(dname, Calibrator);
	CONDITIONAL_NEW_CLASS(dname, TreeExplainer);
	CONDITIONAL_NEW_CLASS(dname, ShapleyExplainer);
	CONDITIONAL_NEW_CLASS(dname, MissingShapExplainer);
	CONDITIONAL_NEW_CLASS(dname, LimeExplainer);
	CONDITIONAL_NEW_CLASS(dname, LinearExplainer);
	CONDITIONAL_NEW_CLASS(dname, KNN_Explainer);
	CONDITIONAL_NEW_CLASS(dname, IterativeSetExplainer);
	CONDITIONAL_NEW_CLASS(dname, AggregatePredsPostProcessor);
	CONDITIONAL_NEW_CLASS(dname, ProbAdjustPostProcessor);
	CONDITIONAL_NEW_CLASS(dname, FairnessPostProcessor);
	MWARN("Warning in PostProcessor::new_polymorphic - Unsupported class %s\n", dname.c_str());
	return NULL;
}

void MultiPostProcessor::get_input_fields(vector<Effected_Field> &fields) const {
	unordered_set<Effected_Field, Effected_Field::HashFunction> full_ls;
	for (size_t i = 0; i < post_processors.size(); ++i) {
		vector<Effected_Field> f;
		post_processors[i]->get_input_fields(f);
		full_ls.insert(f.begin(), f.end());
	}
	fields.insert(fields.end(), full_ls.begin(), full_ls.end());
}
void MultiPostProcessor::get_output_fields(vector<Effected_Field> &fields) const {
	unordered_set<Effected_Field, Effected_Field::HashFunction> full_ls;
	for (size_t i = 0; i < post_processors.size(); ++i) {
		vector<Effected_Field> f;
		post_processors[i]->get_output_fields(f);
		full_ls.insert(f.begin(), f.end());
	}
	fields.insert(fields.end(), full_ls.begin(), full_ls.end());
}

void MultiPostProcessor::init_post_processor(MedModel& model) {
#pragma omp parallel for if (call_parallel_learn && post_processors.size()>1)
	for (int i = 0; i < post_processors.size(); ++i)
		post_processors[i]->init_post_processor(model);

}

void MultiPostProcessor::Learn(const MedFeatures &matrix) {
	if (call_parallel_learn) {
#pragma omp parallel for
		for (int i = 0; i < post_processors.size(); ++i)
			post_processors[i]->Learn(matrix);
	}
	else
		for (int i = 0; i < post_processors.size(); ++i)
			post_processors[i]->Learn(matrix);
}

void MultiPostProcessor::Apply(MedFeatures &matrix) {
	if (call_parallel_apply) {
#pragma omp parallel for if (post_processors.size() > 1)
		for (int i = 0; i < post_processors.size(); ++i)
			post_processors[i]->Apply(matrix);
	}
	else
		for (int i = 0; i < post_processors.size(); ++i)
			post_processors[i]->Apply(matrix);
}

void MultiPostProcessor::dprint(const string &pref) const {
	MLOG("%s :: %s\n", pref.c_str(), my_class_name().c_str());
	for (size_t i = 0; i < post_processors.size(); ++i)
		post_processors[i]->dprint("\t" + pref + "-in-MULTI[" + to_string(i) + "]");
}

MultiPostProcessor::~MultiPostProcessor() {
	for (size_t i = 0; i < post_processors.size(); ++i)
		if (post_processors[i] != NULL) {
			delete post_processors[i];
			post_processors[i] = NULL;
		}
	post_processors.clear();
}

float MultiPostProcessor::get_use_p() {

	if (post_processors.size() == 0)
		use_p = 0.0;
	else {
		use_p = post_processors[0]->use_p;
		for (size_t i = 1; i < post_processors.size(); i++) {
			if (post_processors[i]->use_p != use_p)
				MTHROW_AND_ERR("MultiPostProcessor: use_p inconsistecny (%f vs %f)\n", use_p, post_processors[i]->use_p);
		}
	}

	return use_p;
}

int MultiPostProcessor::get_use_split() {

	if (post_processors.size() == 0)
		use_split = -1;
	else {
		use_split = post_processors[0]->use_split;
		for (size_t i = 1; i < post_processors.size(); i++) {
			if (post_processors[i]->use_split != use_split)
				MTHROW_AND_ERR("MultiPostProcessor: use_split inconsistecny (%d vs %d)\n", use_split, post_processors[i]->use_split);
		}
	}

	return use_split;
}