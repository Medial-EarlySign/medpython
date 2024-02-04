#ifndef __AVERAGE_PREDS_POST_PROCESSOR_H__
#define __AVERAGE_PREDS_POST_PROCESSOR_H__

#include <vector>
#include <Logger/Logger/Logger.h>
#include "PostProcessor.h"

using namespace std;

/**
* Post processor to resample input and apply predictor multiple time with aggregation.
*/
class AggregatePredsPostProcessor : public PostProcessor {
private: 
	MedPredictor *model_predictor = NULL; ///< predictor we're trying to explain
	FeatureProcessor *feature_processor = NULL; ///< feature processor that manipulate the features - for example imputations for missing values. can be new one without learn
	vector<FeatureProcessor *> before_processors; ///< stores all processors need to be applied after ours
	vector<FeatureProcessor *> after_processors; ///< stores all processors need to be applied after ours

	MedModel *p_model; ///< stores on init_post_processor to learn and apply

	void generate_matrix_till_feature_process(const MedFeatures &input_mat, MedFeatures &res) const;
public:
	string feature_processor_type; ///< the feature processor type. if used from model, reffer with prefix MODEL::%s where %s is feature processpr type - currently support single type if exists
	string feature_processor_args; ///< the feature processor type

	bool force_cancel_imputations; ///< if true will force removal of imputations
	bool use_median; ///< if true will fetch median instead of mean
	int resample_cnt; ///< how much to resample
	int batch_size; ///< how many predictions to process together
	bool print_missing_cnt;
	//Add: batch_size, resample_cnt, store_mean (or median)

	AggregatePredsPostProcessor();

	void get_input_fields(vector<Effected_Field> &fields) const;
	void get_output_fields(vector<Effected_Field> &fields) const;

	 /// Global init for general args in all explainers
	int init(map<string, string> &mapper);

	///Learns from predictor and train_matrix (PostProcessor API)
	void Learn(const MedFeatures &train_mat);
	void Apply(MedFeatures &matrix); 

	void init_post_processor(MedModel& model);

	void dprint(const string &pref) const;

	virtual ~AggregatePredsPostProcessor();

	ADD_CLASS_NAME(AggregatePredsPostProcessor)
		ADD_SERIALIZATION_FUNCS(model_predictor, feature_processor, feature_processor_type, 
			feature_processor_args, use_median, resample_cnt, batch_size, force_cancel_imputations, 
			before_processors, after_processors, print_missing_cnt)
};

MEDSERIALIZE_SUPPORT(AggregatePredsPostProcessor)

#endif