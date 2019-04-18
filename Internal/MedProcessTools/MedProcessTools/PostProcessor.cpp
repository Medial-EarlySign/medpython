#include "PostProcessor.h"
#include "Calibration.h"
#include <boost/algorithm/string.hpp>
#include "ExplainWrapper.h"

#define LOCAL_SECTION LOG_MED_MODEL

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
	else
		MTHROW_AND_ERR("Unsupported PostProcessor %s\n", post_processor.c_str());
}

PostProcessor *PostProcessor::make_processor(const string &processor_name, const string &params) {
	return make_processor(post_processor_name_to_type(processor_name), params);
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
	else
		MTHROW_AND_ERR("Unsupported PostProcessor %d\n", type);

	prc->init_from_string(params);

	return prc;
}

void PostProcessor::Learn(MedModel &model, MedPidRepository& rep, const MedFeatures &matrix) {
	MTHROW_AND_ERR("Learn Not implemented in class %s\n", my_class_name().c_str());
}
void PostProcessor::Apply(MedFeatures &matrix) const {
	MTHROW_AND_ERR("Apply Not implemented in class %s\n", my_class_name().c_str());
}

void *PostProcessor::new_polymorphic(string dname)
{
	CONDITIONAL_NEW_CLASS(dname, MultiPostProcessor);
	CONDITIONAL_NEW_CLASS(dname, Calibrator);
	CONDITIONAL_NEW_CLASS(dname, TreeExplainer);
	CONDITIONAL_NEW_CLASS(dname, ShapleyExplainer);
	CONDITIONAL_NEW_CLASS(dname, MissingShapExplainer);
	return NULL;
}

void MultiPostProcessor::Learn(MedModel &model, MedPidRepository& rep, const MedFeatures &matrix) {
	if (call_parallel_learn) {
#pragma omp parallel for
		for (int i = 0; i < post_processors.size(); ++i)
			post_processors[i]->Learn(model, rep, matrix);
	}
	else
		for (int i = 0; i < post_processors.size(); ++i)
			post_processors[i]->Learn(model, rep, matrix);
}

void MultiPostProcessor::Apply(MedFeatures &matrix) const {
#pragma omp parallel for
	for (int i = 0; i < post_processors.size(); ++i)
		post_processors[i]->Apply(matrix);
}