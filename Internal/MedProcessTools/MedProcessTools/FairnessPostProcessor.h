#ifndef __FAIRNESS_POSTPROCESSOR_H__
#define __FAIRNESS_POSTPROCESSOR_H__

#include "PostProcessor.h"
#include <Logger/Logger/Logger.h>

enum Cutoff_Type {
	Score = 0,
	PR = 1,
	Sens = 2
};

extern unordered_map<string, int> map_target_type;
enum Fairness_Target_Type {
	SENS = 0,
	SPEC = 1
};

/// class that hold constrain that translate on reference group into raw score cutoff to measure fainess
class Cutoff_Constraint : public SerializableObject {
public:
	Cutoff_Type type; ///< constraint type
	float value; ///< the value

	void set_type(const string &t);

	ADD_CLASS_NAME(Cutoff_Constraint)
		ADD_SERIALIZATION_FUNCS(type, value)
};

/// A post-processor to adjust probability to fairness between groups.
/// The calibration will fit between constraints thresholds linear transformation A*X+B, to optimize some value
/// The constraint will reduce 1 degree of freedom from the equation.
class FairnessPostProcessor : public PostProcessor {
private:
	unordered_map<float, vector<float>> group_to_score_cutoffs_ranges;
	unordered_map<float, vector<float>> group_to_factors;
	unordered_map<float, vector<float>> group_to_bias;
	bool feature_gen_init = false;
public:
	string feature_name; ///< feautre name to search in matrix created by model_json to generate group for fairness
	string model_json; ///< model json path - important for learn
	float reference_group_val = MED_MAT_MISSING_VALUE; ///< the value for the feature used as refernce group for fairness
	Fairness_Target_Type fairness_target_type; ///< fairness target - SENS of SPEC
	vector<Cutoff_Constraint> constraints; ///< list of constraints cutoffs. Init with comma seperated list for each constraint. The type is prefix with ":". For example Score: PR: SENS:
	double resulotion = 0.1; ///< resulotion for target matching. effect speed/accuracy
	double allow_distance_score = 1.0; ///< max distance allow between score
	double allow_distance_target = 5.0; ///< max distance allow between target
	double allow_distance_cutoff_constraint = 1.0; ///< max distance allow between constraint
	int score_bin_count = 5000; ///< how much bins for score. 0 means no binning
	float score_resulotion = 0; ///< if >0 will apply score resulotion for speedup

	MedPidRepository *p_rep; ///< required for building model for generating model (set by process)
	MedModel group_feature_gen_model; ///< model for generating features for priors (set in learn)
	string resolved_name; ///< resolved feature name (value is set after learn)

	FairnessPostProcessor() { processor_type = PostProcessorTypes::FTR_POSTPROCESS_FAIRNESS; };
	~FairnessPostProcessor() {};
	void init_post_processor(MedModel& mdl) { p_rep = mdl.p_rep; }

	void parse_constrains(const string &s); ///< parses the constrains

	void get_input_fields(vector<Effected_Field> &fields) const;
	void get_output_fields(vector<Effected_Field> &fields) const;

	/// Global init for general args in all explainers
	int init(map<string, string> &mapper);

	///Learns from predictor and train_matrix (PostProcessor API)
	void Learn(const MedFeatures &matrix);
	void Apply(MedFeatures &matrix);

	void dprint(const string &pref) const;

	ADD_CLASS_NAME(FairnessPostProcessor)
		ADD_SERIALIZATION_FUNCS(processor_type, feature_name, resolved_name, reference_group_val,
			group_feature_gen_model, constraints, fairness_target_type, resulotion,
			group_to_score_cutoffs_ranges, group_to_factors, group_to_bias, allow_distance_cutoff_constraint,
			allow_distance_target, allow_distance_score, score_bin_count, score_resulotion)
};

MEDSERIALIZE_SUPPORT(Cutoff_Constraint)
MEDSERIALIZE_SUPPORT(FairnessPostProcessor)

#endif