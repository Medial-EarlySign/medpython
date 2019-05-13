/// @file

#ifndef _MED_MODEL_H_
#define _MED_MODEL_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/RepProcess.h"
#include "MedProcessTools/MedProcessTools/FeatureProcess.h"
#include "MedProcessTools/MedProcessTools/DoCalcFeatProcessor.h"
#include "MedProcessTools/MedProcessTools/FeatureGenerator.h"
#include "MedProcessTools/MedProcessTools/PostProcessor.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include "MedProcessTools/MedProcessTools/MedModelExceptions.h"
#include <boost/property_tree/ptree.hpp>

using namespace boost::property_tree;

/// MedModel learn/apply stages
typedef enum {
	MED_MDL_LEARN_REP_PROCESSORS, ///<Start from learning rep processors
	MED_MDL_LEARN_FTR_GENERATORS, ///<Start from learning feature generators
	MED_MDL_APPLY_FTR_GENERATORS, ///<Start from apply feature generators (already learned)
	MED_MDL_LEARN_FTR_PROCESSORS, ///<Start from learning feature processors
	MED_MDL_APPLY_FTR_PROCESSORS, ///<Start from apply feature processors (already learned)
	MED_MDL_LEARN_PREDICTOR, ///<We have the matrix - learn predcitor
	MED_MDL_APPLY_PREDICTOR, ///<We have trained predcitor, do predict
	MED_MDL_INSERT_PREDS, ///<We have done predict - save results
	MED_MDL_LEARN_POST_PROCESSORS, ///<Start learning of post_processors
	MED_MDL_APPLY_POST_PROCESSORS, ///<start apply of postProcessors
	MED_MDL_END ///<All Done
} MedModelStage;

class PostProcessor;
/// A model = repCleaner + featureGenerator + featureProcessor + MedPredictor
class MedModel final : public SerializableObject {
public:
	/// remember learning set
	int serialize_learning_set = 0;
	int model_json_version = 1; ///< the json version

	/// Repostiroy-level cleaners; to be applied sequentially 
	vector<RepProcessor *> rep_processors;

	/// Feature Generators 
	vector<FeatureGenerator *> generators;

	/// Features-level cleaners; to be applied sequentially 
	vector<FeatureProcessor *> feature_processors;

	///Post Process level - calibrators, explainers...
	vector<PostProcessor *> post_processors;

	/// Predictor
	MedPredictor *predictor = NULL;

	/// Learning samples
	MedSamples *LearningSet = NULL;

	/// Safe Mode for train/test intersection
	int safe_mode = 0;

	/// All required signal names + ids
	unordered_set<string> required_signal_names;
	unordered_set<int> required_signal_ids;

	/// all collected virtual signals (name to type)
	map<string, int> virtual_signals;

	// Constructor/Destructor
	MedModel() { safe_mode = 0; serialize_learning_set = 0; };
	~MedModel() { clear(); };

	void clear();

	MedFeatures features;	///< the created matrix - no need to serialize

	int verbosity = 1; ///< verbosity 0 -> much less printouts in predict

	int generate_masks_for_features = 0;

	// initialize from configuration files
	//int init_rep_processors(const string &fname);
	//int init_feature_generators(const string &fname);

	// staging
	static MedModelStage get_med_model_stage(const string& stage);

	// Add Rep Processorsep
	void add_rep_processor(RepProcessor *processor) { rep_processors.push_back(processor); };
	void add_rep_processors_set(RepProcessorTypes type, vector<string>& signals);
	void add_rep_processors_set(RepProcessorTypes type, vector<string>& signals, string init_string);
	void insert_rep_processor(string init_string, int index);

	// Add Feature Generators
	void add_feature_generator(FeatureGenerator *generator) { generators.push_back(generator); }
	void add_feature_generators(FeatureGeneratorTypes type, vector<string>& signals);
	void add_feature_generators(FeatureGeneratorTypes type, vector<string>& signals, string init_string);
	void add_feature_generator(FeatureGeneratorTypes type, string& signal) { vector<string> signals(1, signal); add_feature_generators(type, signals); }
	void add_feature_generators(FeatureGeneratorTypes type, string& signal, string init_string) { vector<string> signals(1, signal); add_feature_generators(type, signals, init_string); }

	void add_feature_generators(string& name, vector<string>& signals) { add_feature_generators(ftr_generator_name_to_type(name), signals); }
	void add_feature_generators(string& name, vector<string>& signals, string init_string) { add_feature_generators(ftr_generator_name_to_type(name), signals, init_string); }
	void add_feature_generator(string& name, string& signal) { vector<string> signals(1, signal); add_feature_generators(name, signals); }
	void add_feature_generators(string& name, string& signal, string init_string) { vector<string> signals(1, signal); add_feature_generators(name, signals, init_string); }

	void add_age() { generators.push_back(new AgeGenerator); }
	void add_gender() { generators.push_back(new GenderGenerator); }

	void get_all_features_names(vector<string> &feat_names, int before_process_set);

	// Add Feature Processors
	void add_feature_processor(FeatureProcessor *processor) { feature_processors.push_back(processor); };

	void add_feature_processors_set(FeatureProcessorTypes type);
	void add_feature_processors_set(FeatureProcessorTypes type, string init_string);
	void add_feature_processors_set(FeatureProcessorTypes type, vector<string>& features);
	void add_feature_processors_set(FeatureProcessorTypes type, vector<string>& features, string init_string);

	void add_normalizers() { add_feature_processors_set(FTR_PROCESS_NORMALIZER); }
	void add_normalizers(string init_string) { add_feature_processors_set(FTR_PROCESS_NORMALIZER, init_string); }
	void add_normalizers(vector<string>& features) { add_feature_processors_set(FTR_PROCESS_NORMALIZER, features); }
	void add_normalizers(vector<string>& features, string init_string) { add_feature_processors_set(FTR_PROCESS_NORMALIZER, features, init_string); }

	void add_imputers() { add_feature_processors_set(FTR_PROCESS_IMPUTER); }
	void add_imputers(string init_string) { add_feature_processors_set(FTR_PROCESS_IMPUTER, init_string); }
	void add_imputers(vector<string>& features) { add_feature_processors_set(FTR_PROCESS_IMPUTER, features); }
	void add_imputers(vector<string>& features, string init_string) { add_feature_processors_set(FTR_PROCESS_IMPUTER, features, init_string); }

	//post procesors:
	void add_post_processor(PostProcessor *processor) { post_processors.push_back(processor); };

	// general adders for easier handling of config files/lines
	// the idea is to add to a specific set and let the adder create a multi if needed
	int init_from_json_string(string& json_string, const string& fname);
	void init_from_json_file(const string& fname) { vector<string> dummy;  init_from_json_file_with_alterations(fname, dummy); }
	void init_from_json_file_with_alterations_version_1(const string& fname, vector<string>& alterations);
	void init_from_json_file_with_alterations(const string& fname, vector<string>& alterations);
	void add_pre_processors_json_string_to_model(string in_json, string fname);
	void add_rep_processor_to_set(int i_set, const string &init_string);		// rp_type and signal are must have parameters in this case
	void add_feature_generator_to_set(int i_set, const string &init_string);	// fg_type and signal are must have parameters
	void add_feature_processor_to_set(int i_set, int duplicate, const string &init_string);	// fp_type and feature name are must have parameters
	void add_process_to_set(int i_set, int duplicate, const string &init_string); // will auto detect type by which type param is used (rp_type, fg_type OR fp_type)
																				  // and will call the relavant function
	void add_process_to_set(int i_set, const string &init_string) { add_process_to_set(i_set, 0, init_string); }
	void add_post_processor_to_set(int i_set, const string &init_string);

	// Add Predictor
	void set_predictor(MedPredictor *_predictor) { predictor = _predictor; };
	void set_predictor(MedPredictorTypes type) { predictor = MedPredictor::make_predictor(type); }
	void set_predictor(string name) { predictor = MedPredictor::make_predictor(name); }
	void set_predictor(MedPredictorTypes type, string init_string) { predictor = MedPredictor::make_predictor(type, init_string); }
	void set_predictor(string name, string init_string) { predictor = MedPredictor::make_predictor(name, init_string); }

	// signal ids
	void set_required_signal_ids(MedDictionarySections& dict);
	void set_affected_signal_ids(MedDictionarySections& dict);

	// Required signals back-propogation
	void get_required_signal_names(unordered_set<string>& signalNames);
	void get_required_signal_names(vector<string>& signalNames); // same, but get as vector
	// Required signals to generate processed values of target-signals
	void get_required_signal_names_for_processed_values(unordered_set<string>& targetSignalNames, unordered_set<string>& signalNames);
	void get_required_signal_names_for_processed_values(unordered_set<string>& targetSignalNames, vector<string>& signalNames);
	// Get list of Features generated by the model, after everything was applied (RPs, FGs, FPs)
	void get_generated_features_names(vector<string> &feat_names);

	int collect_and_add_virtual_signals(MedRepository &rep);

	/// Initialization : signal ids and tables
	void init_all(MedDictionarySections& dict, MedSignals& sigs);

	// Learn/Apply
	int learn(MedPidRepository& rep, MedSamples* samples) { return learn(rep, samples, MED_MDL_LEARN_REP_PROCESSORS, MED_MDL_END); }
	int learn(MedPidRepository& rep, MedSamples* samples, MedModelStage start_stage, MedModelStage end_stage);
	int apply(MedPidRepository& rep, MedSamples& samples) { return apply(rep, samples, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_END); }
	int apply(MedPidRepository& rep, MedSamples& samples, MedModelStage start_stage, MedModelStage end_stage);

	// Apply on a given Rec : this is needed when someone outside the model runs on Records. No matching method for learn.
	// The process is started with an initialization using init_for_apply_rec
	// Then for each record use : apply_rec , there's a flag for using a copy of the rec rather than the record itself.
	// Calling apply_rec is thread safe, and hence each call returns its own MedFeatures.
	int init_for_apply_rec(MedPidRepository &rep);
	int apply_rec(PidDynamicRec &drec, MedIdSamples idSamples, MedFeatures &_feat, bool copy_rec_flag, int end_stage = MED_MDL_APPLY_FTR_PROCESSORS);

	// De(Serialize)
	virtual void pre_serialization() { if (!serialize_learning_set && LearningSet != NULL) LearningSet = NULL; /*no need to clear(), as this was given by the user*/ }
	ADD_CLASS_NAME(MedModel)
	ADD_SERIALIZATION_FUNCS(rep_processors, generators, feature_processors, predictor, post_processors, generate_masks_for_features, serialize_learning_set, LearningSet)

	int quick_learn_rep_processors(MedPidRepository& rep, MedSamples& samples);
	int learn_rep_processors(MedPidRepository& rep, MedSamples& samples);
	void filter_rep_processors();
	int learn_feature_generators(MedPidRepository &rep, MedSamples *learn_samples);
	int generate_features(MedPidRepository &rep, MedSamples *samples, vector<FeatureGenerator *>& _generators, MedFeatures &features);
	int generate_all_features(MedPidRepository &rep, MedSamples *samples, MedFeatures &features, unordered_set<string>& req_feature_generators);
	int learn_and_apply_feature_processors(MedFeatures &features);
	int learn_feature_processors(MedFeatures &features);
	int apply_feature_processors(MedFeatures &features);
	int apply_feature_processors(MedFeatures &features, vector<unordered_set<string>>& req_features_vec);
	void build_req_features_vec(vector<unordered_set<string>>& req_features_vec);
	void get_applied_generators(unordered_set<string>& req_feature_generators, vector<FeatureGenerator *>& _generators);

	void learn_post_processors(MedPidRepository &rep, MedSamples &post_samples);
	void apply_post_processors(MedFeatures &matrix_after_pred);

	/// following is for debugging, it gets a prefix, and prints it along with information on rep_processors, feature_generators, or feature_processors
	void dprint_process(const string &pref, int rp_flag, int fg_flag, int fp_flag, int predictor_flag, int pp_flag);

	/// following is for debugging : writing the feature to a csv file as a matrix.
	int write_feature_matrix(const string mat_fname);

private:
	void concatAllCombinations(const vector<vector<string> > &allVecs, size_t vecIndex, string strSoFar, vector<string>& result);
	string parse_key_val(string key, string val);
	void fill_list_from_file(const string& fname, vector<string>& list);
	string make_absolute_path(const string& main_file, const string& small_file, bool use_cwd = false);
	void alter_json(string &json_contents, vector<string>& alterations);
	string json_file_to_string(int recursion_level, const string& main_file, vector<string>& alterations, const string& small_file = "", bool add_change_path = false);
	void parse_action(basic_ptree<string, string>& action, vector<vector<string>>& all_action_attrs, int& duplicate, ptree& root, const string& fname);
};

void filter_rep_processors(const vector<string> &current_req_signal_names, vector<RepProcessor *> *rep_processors);
//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(MedModel)

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	*  \brief repository namespace
	*/
	namespace repository {
		/// \brief removes uneeded rep_processors based on needed_sigs and prepares the repository
		/// returns the signal id's neede to read in the repository. MedRepository must be init to read dicts
		vector<int> prepare_repository(MedPidRepository &rep, const vector<int> &needed_sigs,
			vector<int> &phisical_signal_read, vector<RepProcessor *> *rep_processors);
		/// \brief removes uneeded rep_processors based on needed_sigs and prepares the repository
		/// returns the signal id's neede to read in the repository. MedRepository must be init to read dicts
		vector<string> prepare_repository(MedPidRepository &rep, const vector<string> &needed_sigs,
			vector<string> &phisical_signal_read, vector<RepProcessor *> *rep_processors = NULL);
	}
}

#endif
