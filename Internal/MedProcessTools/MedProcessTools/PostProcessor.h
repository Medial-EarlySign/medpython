/// @file
#ifndef __POST_PROCESSOR_H__
#define __POST_PROCESSOR_H__

#include <vector>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include "MedSamples.h"
#include "MedModel.h"

/** @enum
* Post Processors types enum
*/
typedef enum {
	FTR_POSTPROCESS_MULTI, ///<"multi_processor" or "multi" to create MultiPostProcessor
	FTR_POSTPROCESS_CALIBRATOR, ///<"calibrator" to create Calibrator
	FTR_POSTPROCESS_TREE_SHAP, ///< "tree_shap" to create TreeExplainer to explain tree mode or mimic generic model with trees model
	FTR_POSTPROCESS_SHAPLEY, ///< "shapley" to create ShapleyExplainer - model agnostic shapley explainer for model. sample masks using gibbs or GAN
	FTR_POSTPROCESS_MISSING_SHAP, ///< "missing_shap" to create MissingShapExplainer - model agnostic shapley algorithm on trained model to handle prediciton with missing values(retrains new model). much faster impl, because gibbs/GAN is not needed
	FTR_POSTPROCESS_LIME_SHAP, ///< "lime_shap" to create LimeExplainer - model agnostic shapley algorithm with lime on shap values sampling
	FTR_POSTPROCESS_KNN_EXPLAIN,///< "knn" Explainer built on knn principles KNN_Explainer
	FTR_POSTPROCESS_LINEAR, ///< "linear" to create LinearExplainer to explain linear model - importance is score change when putting zero in the feature/group of features
	FTR_POSTPROCESS_LAST
} PostProcessorTypes;

using namespace std;

class MedModel;

/**
* An Abstract PostProcessor class
*/
class PostProcessor : public SerializableObject {
public:
	PostProcessorTypes processor_type = PostProcessorTypes::FTR_POSTPROCESS_LAST;

	virtual void Learn(MedModel &model, MedPidRepository& rep, const MedFeatures &matrix);
	virtual void Apply(MedFeatures &matrix) const;

	void *new_polymorphic(string dname);

	static PostProcessor *make_processor(const string &processor_name, const string &params = "");
	static PostProcessor *make_processor(PostProcessorTypes type, const string &params = "");

	static PostProcessor *create_processor(string &params);

	virtual void init_model(MedModel *mdl) {};

	virtual void dprint(const string &pref) const;

	virtual ~PostProcessor() {};

	ADD_CLASS_NAME(PostProcessor)
		ADD_SERIALIZATION_FUNCS(processor_type)
};
PostProcessorTypes post_processor_name_to_type(const string& post_processor);

/**
* A wrapper for parallel call to post_processors group
*/
class MultiPostProcessor : public PostProcessor {
public:
	vector<PostProcessor *> post_processors;
	bool call_parallel_learn = false;
	bool call_parallel_apply = false;

	MultiPostProcessor() { processor_type = PostProcessorTypes::FTR_POSTPROCESS_MULTI; }

	void Learn(MedModel &model, MedPidRepository& rep, const MedFeatures &matrix);
	void Apply(MedFeatures &matrix) const;

	void init_model(MedModel *mdl);

	void dprint(const string &pref) const;

	~MultiPostProcessor();

	ADD_CLASS_NAME(MultiPostProcessor)
		ADD_SERIALIZATION_FUNCS(post_processors)
};

MEDSERIALIZE_SUPPORT(PostProcessor)
MEDSERIALIZE_SUPPORT(MultiPostProcessor)

#endif
