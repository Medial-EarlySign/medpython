/// @file
#ifndef __POST_PROCESSOR_H__
#define __POST_PROCESSOR_H__

#include <vector>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include "MedSamples.h"
#include "MedModel.h"
#include "Effected_Field.h"

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
	FTR_POSTPROCESS_ITERATIVE_SET, ///< "iterative_set" to create IterativeSetExplainer - model agnostic iterative explainer for model. sample masks using gibbs or GAN
	FTR_POSTPROCESS_AGGREGATE_PREDS, ///< "aggregate_preds" to create AggregatePredsPostProcessor - averaging model predictions after resampling
	FTR_POSTPROCESS_ADJUST, ///< "adjust_probs" to adjust model calibrated predictions according to priors. Creates ProbAdjustPostProcessor
	FTR_POSTPROCESS_FAIRNESS, ///< "fairness_adjust" to adjust model calibrated predictions according to priors. Creates ProbAdjustPostProcessor
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

	// The following variables are used for enabling model to put aside a subset of the learning
	// set for post-processor learning. Either use_split (put aside all ids within a given split)
	// or use_p (put aside a proportion 0 < p < 1 of the ids) can be given. but not both
	int use_split = -1;
	float use_p = 0.0;

	/// List of fields that are used by this post_processor
	virtual void get_input_fields(vector<Effected_Field> &fields) const {};
	/// List of fields that are being effected by this post_processor. Options:
	/// "prediction:X", "attr:$NAME", "str_attr:$NAME", "feature:$NAME", "json"
	virtual void get_output_fields(vector<Effected_Field> &fields) const {};

	virtual void init_post_processor(MedModel& mdl) {};
	virtual void Learn(const MedFeatures &matrix);
	virtual void Apply(MedFeatures &matrix);

	void *new_polymorphic(string dname);

	static PostProcessor *make_processor(const string &processor_name, const string &params = "");
	static PostProcessor *make_processor(PostProcessorTypes type, const string &params = "");

	static PostProcessor *create_processor(string &params);

	virtual void dprint(const string &pref) const;

	virtual ~PostProcessor() {};

	virtual float get_use_p() { return use_p; }
	virtual int get_use_split() { return use_split; }

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

	void Learn(const MedFeatures &matrix);
	void Apply(MedFeatures &matrix);

	void get_input_fields(vector<Effected_Field> &fields) const;
	void get_output_fields(vector<Effected_Field> &fields) const;

	void init_post_processor(MedModel& mdl);

	void dprint(const string &pref) const;

	~MultiPostProcessor();

	float get_use_p();
	int get_use_split();

	ADD_CLASS_NAME(MultiPostProcessor)
		ADD_SERIALIZATION_FUNCS(post_processors)
};

MEDSERIALIZE_SUPPORT(PostProcessor)
MEDSERIALIZE_SUPPORT(MultiPostProcessor)

#endif
