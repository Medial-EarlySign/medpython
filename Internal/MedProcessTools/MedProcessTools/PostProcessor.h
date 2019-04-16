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
	bool call_parallel_learn = true;

	MultiPostProcessor() { processor_type = PostProcessorTypes::FTR_POSTPROCESS_MULTI; }

	void Learn(MedModel &model, MedPidRepository& rep, const MedFeatures &matrix);
	void Apply(MedFeatures &matrix) const;

	ADD_CLASS_NAME(MultiPostProcessor)
		ADD_SERIALIZATION_FUNCS(post_processors)
};

MEDSERIALIZE_SUPPORT(PostProcessor)
MEDSERIALIZE_SUPPORT(MultiPostProcessor)

#endif
