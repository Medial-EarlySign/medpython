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
#include "Effected_Field.h"

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

class ChangeModelInfo :public SerializableObject {
public:
	string change_name = ""; ///< documentation of change name
	string object_type_name = ""; ///< Object type to change - should match my_class_name() of Serialization Object
	vector<string> json_query_whitelist; ///< array of AND condition by regex search for those patterns in the json - to whitelist
	vector<string> json_query_blacklist; ///< array of AND condition by regex search for those patterns in the json - to blacklist
	string change_command = ""; ///< the command, might be "DELETE" ot delete element, "PRINT" to print element ot otherwise pass into object init command
	int verbose_level = 2; ///< 0 - no output, 1 - only warnings, 2 - also info

	int init(map<string, string>& mapper);

	static void parse_json_string(const string &json_content, vector<ChangeModelInfo> &res);

	ADD_CLASS_NAME(ChangeModelInfo)
		ADD_SERIALIZATION_FUNCS(change_name, object_type_name, json_query_whitelist, json_query_blacklist, change_command,
			verbose_level)
};
/// A model = repCleaner + featureGenerator + featureProcessor + MedPredictor
class MedModel final : public SerializableObject {
public:
	string version_info = ""; ///< a string that represents version info - filled in compile time from git info
	/// remember learning set
	int serialize_learning_set = 0;
	int model_json_version = 1; ///< the json version

	long long max_data_in_mem = 0; ///< maximal data size in memory. when <=0 mean auto mode
	/// hint for feature counts to consider in batching. The user can specify count for 
	/// model after feature selection, otherwise it will count all the features from the 
	/// generators without filtering. The number of features needs to be defined pre ahead
	/// of batching process of the samples.
	int model_feature_count_hint = 0; 

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

	/// when having multiple prediction for same samples - how to aggregate preds - mean or median?
	bool take_mean_pred = true;

	/// All required signal names + ids
	unordered_set<string> required_signal_names;
	unordered_set<int> required_signal_ids;
	vector<unordered_set<string> > required_features_vec;
	unordered_set<string> required_feature_generators;

	/// all collected virtual signals (name to type)
	map<string, int> virtual_signals;
	map<string, string> virtual_signals_generic;

	// Constructor/Destructor
	MedModel() { serialize_learning_set = 0; };
	~MedModel() { clear(); };

	void clear();

	MedFeatures features;	///< the created matrix - no need to serialize

	int verbosity = 1; ///< verbosity 0 -> much less printouts in predict

	int generate_masks_for_features = 0;

	// initialize from configuration files
	//int init_rep_processors(const string &fname);
	//int init_feature_generators(const string &fname);

	/// <summary>
	/// change model object in run time
	/// @param change_request.object_type_name - object type name to search for in rep_processors,generators,feature_processors, etc. For example "FeatureNormalizer"
	/// @param change_request.json_query query on the object json to filter on specific attributes. leave empty to operate on all
	/// @param change_request.change_command - The command to send each matched object. use "DELETE" to remove the object, "PRINT" to print object josn into stdout. otherwise it will pass the argument into init function
	/// </summary>
	void change_model(const ChangeModelInfo &change_request);
	/// <summary>
	/// change model object in run time - multiple requests, one by one
	/// @param change_request.object_type_name - object type name to search for in rep_processors,generators,feature_processors, etc. For example "FeatureNormalizer"
	/// @param change_request.json_query query on the object json to filter on specific attributes. leave empty to operate on all
	/// @param change_request.change_command - The command to send each matched object. use "DELETE" to remove the object, "PRINT" to print object josn into stdout. otherwise it will pass the argument into init function
	/// </summary>
	void change_model(const vector<ChangeModelInfo> &change_request);

	// staging
	static MedModelStage get_med_model_stage(const string& stage);

	/// Add Rep Processorsep
	void add_rep_processor(RepProcessor *processor) { rep_processors.push_back(processor); };
	void add_rep_processors_set(RepProcessorTypes type, vector<string>& signals);
	void add_rep_processors_set(RepProcessorTypes type, vector<string>& signals, string init_string);
	void insert_rep_processor(string init_string, int index);

	/// Add Feature Generators
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

	/// Add Feature Processors
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
	int add_pre_processors_json_string_to_model(string in_json, string fname) { vector<string> dummy; return add_pre_processors_json_string_to_model(in_json, fname, dummy); }
	int add_pre_processors_json_string_to_model(string in_json, string fname, vector<string> &alterations, bool add_rep_first = false);
	int add_post_processors_json_string_to_model(string in_json, string fname) { vector<string> dummy; return add_post_processors_json_string_to_model(in_json, fname, dummy); }
	int add_post_processors_json_string_to_model(string in_json, string fname, vector<string> &alterations);
	void add_rep_processor_to_set(int i_set, const string &init_string);		// rp_type and signal are must have parameters in this case
	void add_feature_generator_to_set(int i_set, const string &init_string);	// fg_type and signal are must have parameters
	void add_feature_processor_to_set(int i_set, int duplicate, const string &init_string);	// fp_type and feature name are must have parameters
	void add_process_to_set(int i_set, int duplicate, const string &init_string); // will auto detect type by which type param is used (rp_type, fg_type OR fp_type)
																				  // and will call the relavant function
	void add_process_to_set(int i_set, const string &init_string) { add_process_to_set(i_set, 0, init_string); }
	void add_post_processor_to_set(int i_set, const string &init_string);

	/// Add Predictor
	void set_predictor(MedPredictor *_predictor) { predictor = _predictor; };
	void set_predictor(MedPredictorTypes type) { predictor = MedPredictor::make_predictor(type); }
	void set_predictor(string name) { predictor = MedPredictor::make_predictor(name); }
	void set_predictor(MedPredictorTypes type, string init_string) { predictor = MedPredictor::make_predictor(type, init_string); }
	void set_predictor(string name, string init_string) { predictor = MedPredictor::make_predictor(name, init_string); }
	void replace_predictor_with_json_predictor(string f_json); // given a loaded model and a json file, replaces the model predictor definition with the one in the json.

	/// signal ids
	void set_required_signal_ids(MedDictionarySections& dict, vector<RepProcessor *> &applied_rep_processors,
		vector<FeatureGenerator *> &applied_generators);
	void set_affected_signal_ids(MedDictionarySections& dict, vector<RepProcessor *> &applied_rep_processors);

	// Required signals back-propogation
	void get_required_signal_names(unordered_set<string>& signalNames) const;
	void get_required_signal_names(vector<string>& signalNames) const; // same, but get as vector
	/// Required signals to generate processed values of target-signals
	void get_required_signal_names_for_processed_values(unordered_set<string>& targetSignalNames, unordered_set<string>& signalNames);
	void get_required_signal_names_for_processed_values(unordered_set<string>& targetSignalNames, vector<string>& signalNames);
	/// Get list of Features generated by the model, after everything was applied (RPs, FGs, FPs)
	void get_generated_features_names(vector<string> &feat_names);
	/// get model features size esstimation - number of columns
	int get_nfeatures();
	///returns the duplicate factor of the model - by how much each sample is duplicated by
	int get_duplicate_factor() const;

	int collect_and_add_virtual_signals(MedRepository &rep);

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

	/// Initialization : signal ids and tables
	void init_all(MedDictionarySections& dict, MedSignals& sigs,vector<RepProcessor *> &applied_rep_processors,
		vector<FeatureGenerator *> &applied_generators);
	void init_all(MedDictionarySections& dict, MedSignals& sigs) { init_all(dict, sigs, rep_processors, generators); }
	// Learn/Apply
	int learn(MedPidRepository& rep, MedSamples* samples) { return learn(rep, samples, MED_MDL_LEARN_REP_PROCESSORS, MED_MDL_END); }
	int learn(MedPidRepository& rep, MedSamples* samples, MedModelStage start_stage, MedModelStage end_stage);
	int learn_skip_matrix_train(MedPidRepository &rep, MedSamples *samples, MedModelStage end_stage);
	int apply(MedPidRepository& rep, MedSamples& samples) { return apply(rep, samples, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_END); }
	int apply(MedPidRepository& rep, MedSamples& samples, MedModelStage start_stage, MedModelStage end_stage);

	void no_init_apply_partial(MedPidRepository& rep, MedSamples& samples, const vector<Effected_Field> &requested_outputs);

	// follows are apply methods separating the initialization of the model from the actual apply
	int init_model_for_apply(MedPidRepository &rep, MedModelStage start_stage, MedModelStage end_stage);
	int no_init_apply(MedPidRepository& rep, MedSamples& samples, MedModelStage start_stage, MedModelStage end_stage);

	// Learn with a vector of samples - one for the actual learning, and additional one for each post-processor.
	// PostProcessors that do not require samples, can be assigned empty samples.
	int learn(MedPidRepository& rep, MedSamples& model_learning_set, vector<MedSamples>& post_processors_learning_sets) {
		return learn(rep, model_learning_set, post_processors_learning_sets, MED_MDL_LEARN_REP_PROCESSORS, MED_MDL_END);
	}
	int learn(MedPidRepository& rep, MedSamples& model_learning_set, vector<MedSamples>& post_processors_learning_sets, MedModelStage start_stage, MedModelStage end_stage);

	// Envelopes for normal calling 
	int learn(MedPidRepository& rep, MedSamples& samples) { return learn(rep, &samples); }
	int learn(MedPidRepository& rep, MedSamples& samples, MedModelStage start_stage, MedModelStage end_stage) { return learn(rep, &samples, start_stage, end_stage); }

	// Apply on a given Rec : this is needed when someone outside the model runs on Records. No matching method for learn.
	// The process is started with an initialization using init_for_apply_rec
	// Then for each record use : apply_rec , there's a flag for using a copy of the rec rather than the record itself.
	// Calling apply_rec is thread safe, and hence each call returns its own MedFeatures.
	int init_for_apply_rec(MedPidRepository &rep);
	int apply_rec(PidDynamicRec &drec, MedIdSamples idSamples, MedFeatures &_feat, bool copy_rec_flag, int end_stage = MED_MDL_APPLY_FTR_PROCESSORS);

	// De(Serialize)
	virtual void pre_serialization() { if (!serialize_learning_set && LearningSet != NULL) LearningSet = NULL; /*no need to clear(), as this was given by the user*/ }
	ADD_CLASS_NAME(MedModel)
		ADD_SERIALIZATION_FUNCS(rep_processors, generators, feature_processors, predictor, post_processors, generate_masks_for_features, serialize_learning_set, LearningSet, take_mean_pred, version_info)

		int quick_learn_rep_processors(MedPidRepository& rep, MedSamples& samples);
	int learn_rep_processors(MedPidRepository& rep, MedSamples& samples);
	int learn_all_rep_processors(MedPidRepository& rep, MedSamples& samples);
	void filter_rep_processors();
	int learn_feature_generators(MedPidRepository &rep, MedSamples *learn_samples);
	int generate_features(MedPidRepository &rep, MedSamples *samples, vector<FeatureGenerator *>& _generators, MedFeatures &features);
	int generate_all_features(MedPidRepository &rep, MedSamples *samples, MedFeatures &features, unordered_set<string>& req_feature_generators);
	int learn_and_apply_feature_processors(MedFeatures &features);
	int learn_feature_processors(MedFeatures &features);
	int apply_feature_processors(MedFeatures &features, bool learning);
	int apply_feature_processors(MedFeatures &features, vector<unordered_set<string>>& req_features_vec, bool learning);
	void build_req_features_vec(vector<unordered_set<string>>& req_features_vec) const;
	void get_applied_generators(unordered_set<string>& req_feature_generators, vector<FeatureGenerator *>& _generators) const;

	/// following is for debugging, it gets a prefix, and prints it along with information on rep_processors, feature_generators, or feature_processors
	void dprint_process(const string &pref, int rp_flag, int fg_flag, int fp_flag, int predictor_flag, int pp_flag);

	/// following is for debugging : writing the feature to a csv file as a matrix.
	int write_feature_matrix(const string mat_fname, bool write_attributes = false, bool append = false);

	/// loading a repository (optionally allowing for adjustment to model according to available signals)
	void load_repository(const string& configFile, MedPidRepository& rep, bool allow_adjustment = false);
	void load_repository(const string& configFile, vector<int> ids, MedPidRepository& rep, bool allow_adjustment = false);
	void fit_for_repository(MedPidRepository& rep);
	/// Read binary model from file + json changes req for run-time (empty string for no changes)
	void read_from_file_with_changes(const string &model_binary_path, const string &path_to_json_changes);

	///clones this object into out
	void clone_model(MedModel &out);
	/// copy in modle into this object
	void copy_from_model(MedModel &in);

	/// returns how many samples are done in a single batche (in apply)
	int get_apply_batch_count();

	MedPidRepository *p_rep = NULL; ///< not serialized. stores pointer to rep used in Learn or Apply after call.
private:
	void concatAllCombinations(const vector<vector<string> > &allVecs, size_t vecIndex, string strSoFar, vector<string>& result);
	string parse_key_val(string key, string val);
	void fill_list_from_file(const string& fname, vector<string>& list);
	string make_absolute_path(const string& main_file, const string& small_file, bool use_cwd = false);
	void alter_json(string &json_contents, vector<string>& alterations);
	void insert_environment_params_to_json(string& json_content);
	string json_file_to_string(int recursion_level, const string& main_file, vector<string>& alterations, const string& small_file = "", bool add_change_path = false);
	void parse_action(basic_ptree<string, string>& action, vector<vector<string>>& all_action_attrs, int& duplicate, ptree& root, const string& fname);
	int apply_predictor(MedSamples &samples);

	// Handle learning sets for model/post-processors
	void split_learning_set(MedSamples& inSamples, vector<MedSamples>& post_processors_learning_sets, MedSamples& model_learning_set);
	void split_learning_set(MedSamples& inSamples_orig, vector<MedSamples>& post_processors_learning_sets_orig, vector<MedSamples>& post_processors_learning_sets, MedSamples& model_learning_set);

	void clean_model();
	void find_object(RepProcessor *c, vector<RepProcessor *> &res, vector<RepProcessor **> &res_pointer);
	void find_object(FeatureGenerator *c, vector<FeatureGenerator *> &res, vector<FeatureGenerator **> &res_pointer);
	void find_object(FeatureProcessor *c, vector<FeatureProcessor *> &res, vector<FeatureProcessor **> &res_pointer);
	void find_object(PostProcessor *c, vector<PostProcessor *> &res, vector<PostProcessor **> &res_pointer);
	void find_object(MedPredictor *c, vector<MedPredictor *> &res, vector<MedPredictor **> &res_pointer);


	template <class T> void apply_change(const ChangeModelInfo &change_request, void *obj);

	vector<FeatureGenerator *> applied_generators_to_use;
	vector<RepProcessor *> applied_rep_processors_to_use;
	void get_applied_pipeline(vector<unordered_set<string>> &req_features_vec, unordered_set<string> &required_feature_generators, vector<RepProcessor *> &applied_rep_processors,
		vector<FeatureGenerator *> &applied_generators) const;
	void get_applied_all(vector<unordered_set<string>> &req_features_vec, unordered_set<string> &required_feature_generators, vector<RepProcessor *> &applied_rep_processors,
		vector<FeatureGenerator *> &applied_generators, unordered_set<string>& signalNames) const;
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

		/// \brief removes uneeded rep_processors based on model
		void prepare_repository(const vector<int> &pids, const string &RepositoryPath,
			MedModel &mod, MedPidRepository &rep);

		/// \brief removes uneeded rep_processors based on model
		void prepare_repository(const MedSamples &samples, const string &RepositoryPath,
			MedModel &mod, MedPidRepository &rep);
	}

	namespace medmodel {

		/// \brief given a medmodel object, a rep and samples, do the apply , throws upon a problem
		void apply(MedModel &model, MedSamples &samples, string rep_fname, MedModelStage to_stage = MED_MDL_INSERT_PREDS); // returns just the model : model.features is updated. no need to read samples/already read
		void apply(MedModel &model, string rep_fname, string f_samples, MedSamples &samples, MedModelStage to_stage = MED_MDL_INSERT_PREDS); // returns also a MedSamples object
		void apply(MedModel &model, string rep_fname, string f_samples, MedModelStage to_stage = MED_MDL_INSERT_PREDS); // returns just the model : model.features is updated

	}

	namespace print {
		void medmodel_logging(bool turn_on);
	}
}

#endif
