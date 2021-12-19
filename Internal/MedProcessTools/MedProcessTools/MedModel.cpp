#include "MedModel.h"
#include "MedProcessUtils.h"
#include <omp.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/find.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <cmath>
#include <string>
#include "StripComments.h"
#include "DuplicateProcessor.h"
#include "medial_utilities/medial_utilities/globalRNG.h"

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#define CHECK_CRC 0

using namespace boost::property_tree;

//=======================================================================================
// MedModelStage
//=======================================================================================

map<string, MedModelStage> med_mdl_stage_name_to_stage = {
	{ "learn_rep_processors",MED_MDL_LEARN_REP_PROCESSORS },{ "learn_ftr_generators",MED_MDL_LEARN_FTR_GENERATORS },
	{ "apply_ftr_generators",MED_MDL_APPLY_FTR_GENERATORS },{ "learn_ftr_processors",MED_MDL_LEARN_FTR_PROCESSORS },
	{ "apply_ftr_processors",MED_MDL_APPLY_FTR_PROCESSORS },{ "learn_predictor",MED_MDL_LEARN_PREDICTOR },
	{ "apply_predictor",MED_MDL_APPLY_PREDICTOR },{ "insert_preds",MED_MDL_INSERT_PREDS },
	{ "learn_post_processors",MED_MDL_LEARN_POST_PROCESSORS },{ "apply_post_processors",MED_MDL_APPLY_POST_PROCESSORS },
	{ "end",MED_MDL_END } };

MedModelStage MedModel::get_med_model_stage(const string& stage) {

	string _stage = stage;
	boost::to_lower(_stage);

	if (med_mdl_stage_name_to_stage.find(_stage) == med_mdl_stage_name_to_stage.end()) {
		unordered_map<int, string> op;
		for (const auto &it : med_mdl_stage_name_to_stage)
			op[it.second] = it.first;
		string all_opts = medial::io::get_list_op(op, "\n");
		MTHROW_AND_ERR("unknown stage %s\nOptions are:\n%s", stage.c_str(), all_opts.c_str())
	}
	else
		return med_mdl_stage_name_to_stage[_stage];
}

//=======================================================================================
// MedModel
//=======================================================================================
// Learn with a single MedSamples
//.......................................................................................
int MedModel::learn(MedPidRepository& rep, MedSamples* _samples, MedModelStage start_stage, MedModelStage end_stage) {
	p_rep = &rep;
	MLOG("MedModel() : starting learn process on %d samples, stages %d - %d \n", _samples->nSamples(), start_stage, end_stage);
	// Stage Sanity
	if (start_stage > end_stage) {
		MERR("MedModel learn() : Illegal start and end\n");
		return -1;
	}

	// Set aside parts of the learning set required for post-processors training
	vector<MedSamples> post_processors_learning_sets;
	//MedSamples model_learning_set;
	//split_learning_set(*_samples, post_processors_learning_sets, model_learning_set);

	return learn(rep, *_samples, post_processors_learning_sets, start_stage, end_stage);
}


// special common case of wanting to keep the model matrix training info and retrain predictor
// note that there will be some mess at the moment regarding post_processors and their exact training list - this is TBD !!!
int MedModel::learn_skip_matrix_train(MedPidRepository &rep, MedSamples *samples, MedModelStage end_stage)
{
	p_rep = &rep;
	MLOG("Starting a learn process but skipping the train of the matrix !!!!\n");


	// first applying and generating the matrix
	MLOG("Step 1: Generating matrix with current model without training (using direct apply)\n");
	apply(rep, (*samples), MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_APPLY_FTR_PROCESSORS);

	// second: going for a training session for the predictor and on
	vector<MedSamples> dummy;
	return(learn(rep, (*samples), dummy, MED_MDL_LEARN_PREDICTOR, end_stage));
}

void aggregate_samples(MedFeatures &features, const vector<int> &sample_ids, bool take_mean, bool get_attr) {
	if (features.samples.size() != sample_ids.size())
		MTHROW_AND_ERR("Error aggregate_samples - mismatch size (%zu, %zu)\n", features.samples.size(), sample_ids.size());
	MLOG("INFO: model is aggregating samples\n");
	int max_id = *max_element(sample_ids.begin(), sample_ids.end());
	int min_id = 0;
	vector<bool> seen_id(max_id - min_id + 1);
	MedFeatures final_f;
	final_f.attributes = features.attributes;
	final_f.medf_missing_value = features.medf_missing_value;
	final_f.tags = features.tags;
	final_f.time_unit = features.time_unit;
	vector<vector<float>> col_preds(max_id - min_id + 1);
	vector<int> ids_ord;
	for (size_t i = 0; i < sample_ids.size(); ++i)
	{
		if (!seen_id[sample_ids[i]]) {
			seen_id[sample_ids[i]] = true;
			//take features and create samples. update data,samples,weights,masks:
			final_f.samples.push_back(features.samples[i]);
			if (!features.weights.empty())
				final_f.weights.push_back(features.weights[i]);
			for (const auto &it : features.data)
				final_f.data[it.first].push_back(it.second[i]);
			for (const auto &it : features.masks)
				final_f.masks[it.first].push_back(it.second[i]);
			ids_ord.push_back(sample_ids[i]);
		}
		//collect prediction of sample to aggregated:
		col_preds[sample_ids[i]].push_back(features.samples[i].prediction[0]);
	}
	final_f.init_pid_pos_len();
	//aggregate preds:
	for (size_t i = 0; i < ids_ord.size(); ++i)
	{
		int id = ids_ord[i];
		vector<float> &collected = col_preds[id]; //at least 1 element
		float res;
		if (take_mean)
			res = medial::stats::mean_without_cleaning(collected);
		else
			res = medial::stats::median_without_cleaning(collected, true);
		//override pred
		final_f.samples[i].prediction[0] = res;
		if (get_attr) {
			vector<float> prctile_list = { (float)0.05, (float)0.5, (float)0.95 };
			float std, mean;
			vector<float> res_p;
			if (take_mean) {
				std = medial::stats::std_without_cleaning(collected, res);
				mean = res;
				medial::stats::get_percentiles(collected, prctile_list, res_p);
			}
			else {
				medial::stats::get_mean_and_std_without_cleaning(collected, mean, std);
				medial::stats::get_percentiles(collected, prctile_list, res_p);
			}
			map<string, float> &s = final_f.samples[i].attributes;
			s["pred.std"] = std;
			s["pred.ci_lower"] = res_p[0];
			s["pred.ci_upper"] = res_p[2];
			s["pred.mean"] = mean;
			s["pred.median"] = res_p[1];
		}
	}

	features = move(final_f);
}

void aggregate_samples(MedFeatures &features, bool take_mean, bool get_attr = false) {
	//if sample_ids is empty use pid+time:
	vector<int> sample_ids(features.samples.size());
	int id = 0;
	unordered_map<int, unordered_set<int>> pid_to_times;
	for (size_t i = 0; i < sample_ids.size(); ++i)
	{
		if (pid_to_times.find(features.samples[i].id) == pid_to_times.end()
			|| pid_to_times[features.samples[i].id].find(features.samples[i].time) == pid_to_times[features.samples[i].id].end()) {
			pid_to_times[features.samples[i].id].insert(features.samples[i].time);
			++id;
		}
		sample_ids[i] = id;
	}

	aggregate_samples(features, sample_ids, take_mean, get_attr);
}

// Learn with multiple MedSamples
//.......................................................................................
int MedModel::learn(MedPidRepository& rep, MedSamples& model_learning_set_orig, vector<MedSamples>& post_processors_learning_sets_orig, MedModelStage start_stage, MedModelStage end_stage) {
	p_rep = &rep;
	MedTimer timer;

	// preparing learning sets for model and for post processors (mainly making sure we do the use_p correctly)
	vector<MedSamples> post_processors_learning_sets;
	MedSamples model_learning_set;
	if (post_processors_learning_sets_orig.size() == 0 && post_processors.size() > 0) {
		// if no post processors learning sets were given AND there are post processors, we will use the use_p and use_split rules to create 
		// new sets to train post_processors with and adjust the learning set as well
		split_learning_set(model_learning_set_orig, post_processors_learning_sets_orig, post_processors_learning_sets, model_learning_set);
	}
	else {
		// in this case, either there are no post processors, OR post_processors learning sets were given by the user
		model_learning_set = model_learning_set_orig;
		post_processors_learning_sets = post_processors_learning_sets_orig;
	}


	LearningSet = &model_learning_set;

	// Sanity
	if (post_processors_learning_sets.size() != post_processors.size())
		MTHROW_AND_ERR("MedModel::Learn - Not enough samples given for post-processors learning");

	//init to check we have removed all we can (or if need to create virtual signals?):
	fit_for_repository(rep);
	// init virtual signals
	if (collect_and_add_virtual_signals(rep) < 0) {
		MERR("FAILED collect_and_add_virtual_signals\n");
		return -1;
	}

	// Filter un-needed repository processors
	filter_rep_processors();

	// Set of signals
	if (start_stage <= MED_MDL_APPLY_FTR_GENERATORS) {
		init_all(rep.dict, rep.sigs);

		// Required signals
		required_signal_names.clear();
		required_signal_ids.clear();

		get_required_signal_names(required_signal_names);
		for (string signal : required_signal_names)
			required_signal_ids.insert(rep.dict.id(signal));

		for (int signalId : required_signal_ids) {
			if ((!rep.in_mem_mode_active()) && rep.index.index_table[signalId].is_loaded != 1)
				MLOG("MedModel::learn WARNING signal [%d] = [%s] is required by model but not loaded in rep\n",
					signalId, rep.dict.name(signalId).c_str());;
		}
	}
	//dprint_process("==> In Learn (1) <==", 2, 0, 0);

	// Learn RepProcessors
	if (start_stage <= MED_MDL_LEARN_REP_PROCESSORS) {
		MLOG("MedModel() : starting learn rep processors on %d samples\n", model_learning_set.nSamples(), start_stage, end_stage);
		timer.start();
		if (learn_rep_processors(rep, model_learning_set) < 0) { //??? why are rep processors initialized for ALL time points in an id??
			MERR("MedModel learn() : ERROR: Failed learn_rep_processors()\n");
			return -1;
		}
		timer.take_curr_time();
		MLOG("MedModel::learn() : learn rep processors time %g ms\n", timer.diff_milisec());
	}
	if (end_stage <= MED_MDL_LEARN_REP_PROCESSORS)
		return 0;

	// Learn Feature Generators
	if (start_stage <= MED_MDL_LEARN_FTR_GENERATORS) {
		timer.start();
		if (learn_feature_generators(rep, &model_learning_set) < 0) {
			MERR("MedModel learn() : ERROR: Failed learn_feature_generators\n");
			return -1;
		}
		timer.take_curr_time();
		MLOG("MedModel::learn() : learn feature generators %g ms\n", timer.diff_milisec());
	}

	//dprint_process("==> In Learn (2) <==", 0, 2, 0);
	if (end_stage <= MED_MDL_LEARN_FTR_GENERATORS)
		return 0;

	// Generate features
	unordered_set<string> empty_set;
	if (start_stage <= MED_MDL_APPLY_FTR_GENERATORS) {
		features.clear();
		features.set_time_unit(model_learning_set.time_unit);

		timer.start();
		if (generate_all_features(rep, &model_learning_set, features, empty_set) < 0) {
			MERR("MedModel learn() : ERROR: Failed generate_all_features()\n");
			return -1;
		}
		timer.take_curr_time();
		if (CHECK_CRC)
			MLOG("MedModel::learn() : generating learn matrix time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());
		else
			MLOG("MedModel::learn() : generating learn matrix time %g ms\n", timer.diff_milisec());
	}
	if (end_stage <= MED_MDL_APPLY_FTR_GENERATORS)
		return 0;

	// Learn and/or apply Feature processors and apply
	if (start_stage <= MED_MDL_LEARN_FTR_PROCESSORS && end_stage >= MED_MDL_APPLY_FTR_PROCESSORS) {
		// Learn + Apply
		timer.start();
		if (learn_and_apply_feature_processors(features) < 0) {
			MERR("MedModel::learn() : ERROR: Failed learn_and_apply_feature_processors()\n");
			return -1;
		}
		timer.take_curr_time();
		if (CHECK_CRC)
			MLOG("MedModel::learn() : feature processing learn and apply time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());
		else
			MLOG("MedModel::learn() : feature processing learn and apply time %g ms\n", timer.diff_milisec());
	}
	else if (start_stage <= MED_MDL_LEARN_FTR_GENERATORS && end_stage < MED_MDL_APPLY_FTR_PROCESSORS) {
		// Just learn feature processors
		timer.start();
		if (learn_feature_processors(features) < 0) {
			MERR("MedModel::learn() : ERROR: Failed learn_feature_processors()\n");
			return -1;
		}
		timer.take_curr_time();
		if (CHECK_CRC)
			MLOG("MedModel::learn() : feature processing learn time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());
		else
			MLOG("MedModel::learn() : feature processing learn time %g ms\n", timer.diff_milisec());
	}
	else if (start_stage <= MED_MDL_APPLY_FTR_PROCESSORS) {
		// Just apply feature processors
		timer.start();
		if (generate_masks_for_features) features.mark_imputed_in_masks();
		if (apply_feature_processors(features, true) < 0) {
			MERR("MedModel::apply() : ERROR: Failed apply_feature_processors()\n");
			return -1;
		}
		timer.take_curr_time();
		if (CHECK_CRC)
			MLOG("MedModel::learn() : feature processing time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());
		else
			MLOG("MedModel::learn() : feature processing time %g ms\n", timer.diff_milisec());
	}
	if (end_stage <= MED_MDL_APPLY_FTR_PROCESSORS)
		return 0;

	// Learn predictor
	if (start_stage <= MED_MDL_LEARN_PREDICTOR && predictor != NULL) {
		timer.start();
		//MLOG("features: %d : \n", features.data.size());
		//for (auto &e : features.data) { MLOG("%s\n", e.first.c_str()); };
		int rc = predictor->learn(features);
		timer.take_curr_time();
		MLOG("MedModel::learn() : model train time: %g ms\n", timer.diff_milisec());
		if (rc != 0)
			return rc;
	}
	if (end_stage < MED_MDL_LEARN_POST_PROCESSORS)
		return 0;

	// Learn post-processesors. Possibly on different subset of samples
	// A VERY IMPORTANT NOTE: 
	// Currently, post-processors are assumed to be independent of each other. Thus,
	// we do not apply a post-processor after learning it and before learning the next
	// one, and all post processor used the pre-processed scores (and matrix). This
	// saves us a lot of cross-validation mess ...
	if (start_stage <= MED_MDL_LEARN_POST_PROCESSORS) {
		MLOG("MedModel::learn() : learn post_processors\n");
		//aggregate preds from learn process - features might be resampled.
		if (!post_processors.empty()) {
			bool need_agg = false;
			vector<string> m_tags = { "" };
			for (const FeatureProcessor * fp : feature_processors)
			{
				need_agg = fp->select_learn_matrix(m_tags) != "";
				if (need_agg)
					break;
			}
			if (need_agg) //to save time - check is need to aggregate - has some FP that generates new matrix
				aggregate_samples(features, take_mean_pred);
		}

		for (size_t i = 0; i < post_processors.size(); ++i) {
			post_processors[i]->init_post_processor(*this);

			// Prepare learning matrix for post-processor and learn
			if (post_processors_learning_sets[i].idSamples.empty())
				post_processors[i]->Learn(features);
			else {
				MedFeatures origFeatures = move(features);
				apply(rep, post_processors_learning_sets[i], MedModelStage::MED_MDL_APPLY_FTR_GENERATORS, MedModelStage::MED_MDL_APPLY_PREDICTOR);
				cerr << "#####################################HERE###############################################\n";
				post_processors[i]->Learn(features);
				features = move(origFeatures);
			}
		}
	}

	return 0;
}

//.......................................................................................
// Apply
int MedModel::apply(MedPidRepository& rep, MedSamples& samples, MedModelStage start_stage, MedModelStage end_stage) {

	//maximal number of samples to apply together in a batch. takes into account duplicate factor of samples, # of features
	// the goal is to have a matrix with less than MAX_INT elements. can be changed later to other number.
	int max_smp_batch;
	if (start_stage < MED_MDL_APPLY_FTR_PROCESSORS)
		max_smp_batch = get_apply_batch_count();

	init_model_for_apply(rep, start_stage, end_stage);

	if (start_stage >= MED_MDL_APPLY_FTR_PROCESSORS || samples.nSamples() <= max_smp_batch)
		return no_init_apply(rep, samples, start_stage, end_stage);
	else {
		//Do in batches:
		vector<MedSample> flat_samples;
		samples.export_to_sample_vec(flat_samples);
		int original_size = (int)flat_samples.size();
		int num_of_batches = int(ceil(double(original_size) / max_smp_batch));
		MLOG("INFO:: MedModel::apply - split into %d batches, each of size %d\n",
			num_of_batches, max_smp_batch);

		MedProgress progress("MedModel::Apply_IN_BATCHES", num_of_batches, 30, 1);
		int pos = 0;
		while (pos < original_size) {
			int batch_size = original_size - pos;
			if (batch_size > max_smp_batch)
				batch_size = max_smp_batch;
			vector<MedSample> flat_batch(flat_samples.data() + pos, flat_samples.data() + pos + batch_size);
			MedSamples batch;
			//1. store in batch
			batch.import_from_sample_vec(flat_batch);
			//2. Apply on batch
			no_init_apply(rep, batch, start_stage, end_stage);
			//3. fetch results from batch into samples
			int inner_p = pos;
			for (size_t i = 0; i < batch.idSamples.size(); ++i)
				for (size_t j = 0; j < batch.idSamples[i].samples.size(); ++j) {
					flat_samples[inner_p] = move(batch.idSamples[i].samples[j]);
					++inner_p;
				}
			//4. advance into next batch:
			pos += batch_size;
			progress.update();
		}
		//convert from flat_samples back to samples:
		samples.import_from_sample_vec(flat_samples);

		//Clean MedFeature samples and data:
		features.samples = move(flat_samples);
		for (auto &it : features.data) //Earse all - it's only partial data, clear memory
			it.second.clear();
	}
	return 0;
}

int MedModel::get_apply_batch_count() {
	//maximal number of samples to apply together in a batch. takes into account duplicate factor of samples, # of features
	// the goal is to have a matrix with less than MAX_INT elements. can be changed later to other number.
	long long max_sz = max_data_in_mem;
	if (max_sz <= 0) {
		//TODO: change to use size to suit free memory in the machine
		max_sz = INT_MAX;
	}
	long long num_of_features = model_feature_count_hint;
	int max_model_feature_count = get_nfeatures();
	if (num_of_features <= 0 || max_model_feature_count < model_feature_count_hint)
		num_of_features = max_model_feature_count;

	long long max_smp_batch = (long long)(((max_sz) / (get_duplicate_factor()*num_of_features)) * 0.95) - 1;
	if (max_smp_batch > INT_MAX)
		max_smp_batch = INT_MAX;
	return (int)max_smp_batch;
}

//-------------------------------------------------------------------------------------------------------------------------------
int MedModel::init_model_for_apply(MedPidRepository &rep, MedModelStage start_stage, MedModelStage end_stage)
{
	p_rep = &rep;
	// Stage Sanity
	if (end_stage < MED_MDL_APPLY_FTR_GENERATORS) {
		MERR("MedModel apply() : Illegal end stage %d\n", end_stage);
		return -1;
	}

	//only perform when needed
	if (start_stage < MED_MDL_APPLY_PREDICTOR) {
		//init to check we have remove all we can (or if need to create virtual signals?):
		fit_for_repository(rep);
		// init virtual signals
		if (collect_and_add_virtual_signals(rep) < 0) {
			MERR("FAILED collect_and_add_virtual_signals\n");
			return -1;
		}

	}
	//dprint_process("==> In Apply (1) <==", 2, 0, 0);

	// Build sets of required features at each stage of processing
	// The last entry tells us which features to generate
	required_features_vec.clear();
	build_req_features_vec(required_features_vec);
	required_feature_generators = required_features_vec[feature_processors.size()];

	if (start_stage <= MED_MDL_APPLY_FTR_GENERATORS) {

		// Initialize
		init_all(rep.dict, rep.sigs);

		// Required signals
		required_signal_names.clear();
		required_signal_ids.clear();

		get_required_signal_names(required_signal_names);
		for (string signal : required_signal_names)
			required_signal_ids.insert(rep.dict.id(signal));

		for (int signalId : required_signal_ids) {
			if ((!rep.in_mem_mode_active()) && rep.index.index_table[signalId].is_loaded != 1)
				MLOG("MedModel::apply WARNING signal [%d] = [%s] is required by model but not loaded in rep\n",
					signalId, rep.dict.name(signalId).c_str());;
		}

		//dprint_process("==> In Apply (2) <==", 2, 0, 0);
	}

	if (end_stage <= MED_MDL_APPLY_FTR_PROCESSORS)
		return 0;

	if (start_stage <= MED_MDL_APPLY_PREDICTOR && predictor != NULL) {
		predictor->prepare_predict_single();
	}

	if (start_stage <= MED_MDL_APPLY_POST_PROCESSORS) {

		if (verbosity > 0) MLOG("Initializing %d postprocessors\n", (int)post_processors.size());
		for (size_t i = 0; i < post_processors.size(); ++i)
			post_processors[i]->init_post_processor(*this);

	}

	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------------
int MedModel::no_init_apply(MedPidRepository& rep, MedSamples& samples, MedModelStage start_stage, MedModelStage end_stage)
{
	p_rep = &rep;
	// Stage Sanity
	if (end_stage < MED_MDL_APPLY_FTR_GENERATORS) {
		MERR("MedModel apply() : Illegal end stage %d\n", end_stage);
		return -1;
	}

	if (start_stage <= MED_MDL_APPLY_FTR_GENERATORS) {

		// Generate features
		features.clear();
		features.set_time_unit(samples.time_unit);
		if (verbosity > 0) MLOG("MedModel apply() : before generate_all_features() samples of %d ids\n", samples.idSamples.size());
		if (generate_all_features(rep, &samples, features, required_feature_generators) < 0) {
			MERR("MedModel apply() : ERROR: Failed generate_all_features()\n");
			return -1;
		}
		if (verbosity > 0) MLOG("MedModel apply() : after generate_all_features() samples of %d ids\n", samples.idSamples.size());
		if (samples.copy_attributes(features.samples) != 0) {
			MERR("Insertion of Feature Generators attributes to samples failed\n");
			return -1;
		}
	}

	if (end_stage <= MED_MDL_APPLY_FTR_GENERATORS) {
		if (samples.insert_preds(features) != 0) {
			MERR("Insertion of predictions to samples failed\n");
			return -1;
		}
		return 0;
	}

	// Process Features
	if (start_stage <= MED_MDL_APPLY_FTR_PROCESSORS) {
		if (verbosity > 0) MLOG("MedModel apply() on %d samples : before applying feature processors : generate_masks = %d\n", samples.idSamples.size(), generate_masks_for_features);
		if (generate_masks_for_features) features.mark_imputed_in_masks();
		if (apply_feature_processors(features, required_features_vec, false) < 0) {
			MERR("MedModel::apply() : ERROR: Failed apply_feature_cleaners()\n");
			return -1;
		}
		if (verbosity > 0) MLOG("MedModel apply() : after applying feature processors\n", samples.idSamples.size());
	}

	if (end_stage <= MED_MDL_APPLY_FTR_PROCESSORS || predictor == NULL)
		return 0;

	// Apply predictor
	if (start_stage <= MED_MDL_APPLY_PREDICTOR) {
		if (verbosity > 0) MLOG("before predict: for MedFeatures of: %d x %d\n", features.data.size(), features.samples.size());
		if (features.samples.size() == 1 && !predictor->predict_single_not_implemented()) {
			vector<float> pred_res, features_vec(features.data.size());
			int i_feat = 0;
			for (const auto &it : features.data)
			{
				features_vec[i_feat] = it.second[0];
				++i_feat;
			}
			predictor->predict_single(features_vec, pred_res);
			features.samples[0].prediction = move(pred_res);
		}
		else {
			if (predictor->predict(features) < 0) {
				MERR("Predictor failed\n");
				return -1;
			}
			//MLOG("samples %d features.samples %d n_preds %d\n", samples.nSamples(), features.samples.size(), predictor->n_preds_per_sample());
			bool need_agg = samples.nSamples() != features.samples.size();
			if (need_agg) //to save time - check is need to aggregate - has some FP that generates new matrix
				aggregate_samples(features, take_mean_pred, true);
		}
	}

	if (end_stage <= MED_MDL_APPLY_PREDICTOR)
		return 0;

	if (start_stage <= MED_MDL_INSERT_PREDS && end_stage < MED_MDL_APPLY_POST_PROCESSORS) { //insert preds now only if has no post_processors
		if (verbosity > 0) MLOG("Inserting predictions\n");
		if (samples.insert_preds(features) != 0) {
			MERR("Insertion of predictions to samples failed\n");
			return -1;
		}
		return 0;
	}

	if (start_stage <= MED_MDL_APPLY_POST_PROCESSORS) {

		if (verbosity > 0) MLOG("Applying %d postprocessors\n", (int)post_processors.size());
		MedTimer pp_timer("post_processors"); pp_timer.start();
		for (size_t i = 0; i < post_processors.size(); ++i)
			post_processors[i]->Apply(features);
		pp_timer.take_curr_time();
		if (verbosity > 0) MLOG("Finished postprocessors within %2.1f seconds\n", pp_timer.diff_sec());

		if (samples.insert_preds(features) != 0) {
			MERR("Insertion of predictions to samples failed\n");
			return -1;
		}

		if (samples.copy_attributes(features.samples) != 0) {
			MERR("Insertion of post_process to samples failed\n");
			return -1;
		}
	}

	return 0;
}

//.......................................................................................
// Learn rep-cleaning
int MedModel::quick_learn_rep_processors(MedPidRepository& rep, MedSamples& samples) {

	vector<int> rc(rep_processors.size(), 0);

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < rep_processors.size(); i++)
		rc[i] = rep_processors[i]->learn(rep, samples);

	for (auto RC : rc) if (RC < 0)	return -1;
	return 0;
}

//.......................................................................................
int MedModel::learn_feature_generators(MedPidRepository &rep, MedSamples *learn_samples)
{

	vector<int> rc(generators.size(), 0);
	//omp_set_nested(true);

	MedProgress progress("MedModel::learn_feature_generators", (int)generators.size(), 60, 1);
	//#pragma omp parallel for //num_threads(4) //schedule(dynamic)
	for (int i = 0; i < generators.size(); i++) {
		rc[i] = generators[i]->learn(rep, *learn_samples, rep_processors);
		progress.update();
	}

	for (auto RC : rc) if (RC < 0)	return -1;
	return 0;
}

//.......................................................................................
void MedModel::get_applied_generators(unordered_set<string>& req_feature_generators, vector<FeatureGenerator *>& _generators) {

	for (auto& generator : generators) {
		if (req_feature_generators.empty() || generator->filter_features(req_feature_generators) != 0)
			_generators.push_back(generator);
	}
}

//.......................................................................................
int MedModel::generate_all_features(MedPidRepository &rep, MedSamples *samples, MedFeatures &features, unordered_set<string>& req_feature_generators) {

	vector<FeatureGenerator *> _generators;
	get_applied_generators(req_feature_generators, _generators);

	int res = generate_features(rep, samples, _generators, features);
	//print rep_processors_summary:
	for (unsigned int i = 0; i < rep_processors.size(); i++)
		rep_processors[i]->make_summary();

	return res;
}

//.......................................................................................
int MedModel::generate_features(MedPidRepository &rep, MedSamples *samples, vector<FeatureGenerator *>& _generators, MedFeatures &features)
{
	vector<int> req_signals;
	for (int signalId : required_signal_ids)
		req_signals.push_back(signalId);
	for (auto &vsig : virtual_signals) {
		//		MLOG("GENERATE: vsig %s %d\n", vsig.first.c_str(), vsig.second);
		req_signals.push_back(rep.sigs.sid(vsig.first));
	}

	for (auto &vsig : virtual_signals_generic) {
		//		MLOG("GENERATE: vsig %s %d\n", vsig.first.c_str(), vsig.second);
		req_signals.push_back(rep.sigs.sid(vsig.first));
	}

	// prepare for generation
	for (auto& generator : _generators)
		generator->prepare(features, rep, *samples);

	// preparing records and features for threading
	int N_tot_threads = omp_get_max_threads();
	//	MLOG("MedModel::learn/apply() : feature generation with %d threads\n", N_tot_threads);
	vector<PidDynamicRec> idRec(N_tot_threads);
	features.init_all_samples(samples->idSamples);

	// if attr_train_weight exists (testing on first sample), we enter weights to features
	if (features.samples.size() > 0) {
		if (features.samples[0].attributes.find("train_weight") != features.samples[0].attributes.end()) {
			features.weights.clear();
			for (auto &s : features.samples)
				features.weights.push_back(s.attributes["train_weight"]);
		}
	}

	// Resize data vectors and collect pointers
	int samples_size = (int)features.samples.size();
	for (auto &generator : _generators) {
		generator->p_data.clear();
		if (generator->iGenerateWeights) {
			features.weights.resize(samples_size, 0);
			if (generator->names.size() != 1)
				MTHROW_AND_ERR("Cannot generate weights using a multi-feature generator (type %d generates %d features)\n", generator->generator_type, (int)generator->names.size());
			generator->p_data.push_back(&(features.weights[0]));
		}
		else {
			for (string& name : generator->names)
				features.data[name].resize(samples_size, 0);
			generator->get_p_data(features);
		}
	}

	// Get Required signals per processor (set)
	vector<unordered_set<int> > current_req_signal_ids(rep_processors.size());
	for (unsigned int i = 0; i < rep_processors.size(); i++)
		get_all_required_signal_ids(current_req_signal_ids[i], rep_processors, i, generators);

	// Loop on ids
	int RC = 0;
	int thrown = 0;

	try {
#pragma omp parallel for schedule(dynamic) if (samples->idSamples.size() > 1)
		for (int j = 0; j < samples->idSamples.size(); j++) {
			try {
				MedIdSamples& pid_samples = samples->idSamples[j];
				int n_th = omp_get_thread_num();
				int rc = 0;

				// Generate DynamicRec with all relevant signals
				if (idRec[n_th].init_from_rep(std::addressof(rep), pid_samples.id, req_signals, (int)pid_samples.samples.size()) < 0) rc = -1;

				// Apply rep-processing
				for (unsigned int i = 0; i < rep_processors.size(); i++) {
					rep_processors[i]->dprint(to_string(i), 0);
					if (rep_processors[i]->conditional_apply(idRec[n_th], pid_samples, current_req_signal_ids[i]) < 0) rc = -1;
				}

				// Generate Features
				for (auto& generator : _generators)
					if (generator->generate(idRec[n_th], features) < 0)	rc = -1;


#pragma omp critical 
				if (rc < 0) RC = -1;
			}
			catch (...) {
				// have to catch each thread
				thrown = -1;
				MERR("!!! Got thrown set to %d\n", thrown);
				//throw std::runtime_error("thrown from openmp");
			}
		}
	}
	catch (...) {
		MERR("Caught An error in generate_feature()\n");
		throw;
	}
	// call summary for generation
	for (auto& generator : _generators)
		generator->make_summary();

	if (thrown < 0) MTHROW_AND_ERR("Thrown in generate_features()\n");
	//throw thrown; // throwing if needed
	return RC;

}

//.......................................................................................
int MedModel::learn_and_apply_feature_processors(MedFeatures &features)
{
	vector<string> matrix_versions = { "" };
	vector<MedFeatures> allocated_mem_mat;

	unordered_map<string, MedFeatures *> ver_set;
	ver_set[""] = &features;
	for (auto& processor : feature_processors) {
		string ver = processor->select_learn_matrix(matrix_versions);
		if (ver_set.find(ver) == ver_set.end()) {
			//create new version:
			matrix_versions.push_back(ver);
			MedFeatures new_mat = features;
			allocated_mem_mat.push_back(move(new_mat));
			ver_set[ver] = &allocated_mem_mat.back();
		}

		if (processor->learn(*ver_set.at(ver)) < 0) return -1;
		//Apply on all versions:
		for (const auto &it : ver_set)
			if (processor->apply(*it.second, true) < 0) return -1;
	}
	//update feature with last version for predcitor - if needed and has more than one version
	if (matrix_versions.size() > 1) {
		features = move(*ver_set.at(matrix_versions.back()));
		MLOG("INFO:: has more than one version of matrix in learn_and_apply\n");
	}
	return 0;
}

//.......................................................................................
int MedModel::learn_feature_processors(MedFeatures &features)
{
	vector<string> matrix_versions = { "" };
	vector<MedFeatures> allocated_mem_mat;

	unordered_map<string, MedFeatures *> ver_set;
	ver_set[""] = &features;
	for (auto& processor : feature_processors) {
		string ver = processor->select_learn_matrix(matrix_versions);
		if (ver_set.find(ver) == ver_set.end()) {
			//create new version:
			matrix_versions.push_back(ver);
			MedFeatures new_mat = features;
			allocated_mem_mat.push_back(move(new_mat));
			ver_set[ver] = &allocated_mem_mat.back();
		}

		if (processor->learn(*ver_set.at(ver)) < 0) return -1;
	}
	if (matrix_versions.size() > 1) {
		features = move(*ver_set.at(matrix_versions.back()));
		MLOG("INFO:: has more than one version of matrix in learn\n");
	}

	return 0;
}

//.......................................................................................
int MedModel::apply_feature_processors(MedFeatures &features, bool learning)
{

	vector<unordered_set<string> > req_features_vec(feature_processors.size());
	for (auto& _set : req_features_vec)
		_set.clear();
	return apply_feature_processors(features, req_features_vec, learning);

}
//.......................................................................................
int MedModel::apply_feature_processors(MedFeatures &features, vector<unordered_set<string>>& req_features_vec, bool learning)
{
	int n = (int)feature_processors.size();
	for (int i = 0; i < n; i++) {
		if (feature_processors[i]->apply(features, req_features_vec[n - i - 1], learning) < 0)
			return -1;
	}

	return 0;
}

//.......................................................................................
void MedModel::build_req_features_vec(vector<unordered_set<string>>& req_features_vec) {

	req_features_vec.resize(feature_processors.size() + 1);
	req_features_vec[0] = {};
	for (int i = 1; i <= feature_processors.size(); i++) {
		size_t idx = feature_processors.size() - i;
		feature_processors[idx]->update_req_features_vec(req_features_vec[i - 1], req_features_vec[i]);
	}

	return;
}

// Learn rep-processors iteratively, must be serial...
//.......................................................................................
int MedModel::learn_rep_processors(MedPidRepository& rep, MedSamples& samples) {

	vector<RepProcessor *> temp_processors;
	for (unsigned int i = 0; i < rep_processors.size(); i++) {
		unordered_set<int> current_req_signal_ids;
		get_all_required_signal_ids(current_req_signal_ids, rep_processors, i, generators);
		if (rep_processors[i]->conditional_learn(rep, samples, temp_processors, current_req_signal_ids) < 0) return -1;
		temp_processors.push_back(rep_processors[i]);
	}

	return 0;
}

// Learn rep-processors iteratively, must be serial...
//.......................................................................................
int MedModel::learn_all_rep_processors(MedPidRepository& rep, MedSamples& samples) {

	vector<RepProcessor *> temp_processors;
	for (unsigned int i = 0; i < rep_processors.size(); i++) {
		if (rep_processors[i]->learn(rep, samples, temp_processors) < 0) return -1;
		temp_processors.push_back(rep_processors[i]);
	}

	return 0;
}


// Filter rep-processors that are not used, iteratively
//.......................................................................................
void MedModel::filter_rep_processors() {

	vector<RepProcessor *> filtered_processors;
	bool did_something = false;
	for (unsigned int i = 0; i < rep_processors.size(); i++) {
		unordered_set<string> current_req_signal_names;
		get_all_required_signal_names(current_req_signal_names, rep_processors, i, generators);
		if (!rep_processors[i]->filter(current_req_signal_names))
			filtered_processors.push_back(rep_processors[i]);
		else {//cleaning uneeded rep_processors!:
			delete rep_processors[i];
			rep_processors[i] = NULL;
			did_something = true;
		}
	}
	if (did_something)
		MLOG("Filtering unneeded rep_processors. keeping %zu rep_proccessors out of %zu\n",
			filtered_processors.size(), rep_processors.size());

	rep_processors.swap(filtered_processors);
}

// Set ids of required signals
//.......................................................................................
void MedModel::set_required_signal_ids(MedDictionarySections& dict) {

	required_signal_ids.clear();

	for (RepProcessor *processor : rep_processors)
		processor->set_required_signal_ids(dict);

	for (FeatureGenerator *generator : generators)
		generator->set_required_signal_ids(dict);

}

void MedModel::concatAllCombinations(const vector<vector<string> > &allVecs, size_t vecIndex, string strSoFar, vector<string>& result)
{
	if (vecIndex >= allVecs.size())
	{
		result.push_back(strSoFar.substr(0, strSoFar.length() - 1));
		return;
	}
	for (size_t i = 0; i < allVecs[vecIndex].size(); i++)
		concatAllCombinations(allVecs, vecIndex + 1, strSoFar + allVecs[vecIndex][i] + ";", result);
}
string MedModel::parse_key_val(string key, string val) {
	if (val.length() > 0 && val[0] != '{' && val.find('=') != string::npos) {
		if (val.find(';') == string::npos)
			MLOG("found as-is literal string [%s]\n", val.c_str());
		return val;
	}
	else return key + "=" + val;
}
void MedModel::fill_list_from_file(const string& fname, vector<string>& list) {
	ifstream inf(fname);
	if (!inf)
		MTHROW_AND_ERR("can't open file %s for read\n", fname.c_str());

	string curr_line;
	int lines = 0;
	while (getline(inf, curr_line)) {
		if (curr_line[curr_line.size() - 1] == '\r')
			curr_line.erase(curr_line.size() - 1);
		size_t npos = curr_line.find("//");
		if (npos == 0)
			continue;
		lines++;
		if (npos < string::npos)
			list.push_back(curr_line.substr(0, npos));
		else list.push_back(curr_line);
	}
	MLOG("read %d lines from: %s\n", lines, fname.c_str());
	inf.close();

}
string MedModel::make_absolute_path(const string& main_file, const string& small_file, bool use_cwd) {
	boost::filesystem::path p(main_file);
	string main_file_path = p.parent_path().string();
	if (use_cwd)
		main_file_path = run_current_path;

	if (
		(small_file.size() > 2 && (small_file[0] == '/' || small_file[1] == ':')) ||
		(main_file_path.size() == 0)
		)
		return small_file;

	string abs;
	if (main_file_path.empty())
		abs = small_file;
	else
		abs = main_file_path + path_sep() + small_file;
	if (use_cwd)
		MLOG_D("resolved relative path using cwd [%s] to [%s]\n", small_file.c_str(), abs.c_str());
	else
		MLOG_D("resolved relative path [%s] to [%s]\n", small_file.c_str(), abs.c_str());
	return abs;
}
void MedModel::alter_json(string &json_contents, vector<string>& alterations) {

	if (alterations.size() == 0) return;

	// Alterations strings are of the format from::to
	vector<string> fields;
	MLOG_D("Json : replacing ");
	for (string& alt : alterations) {
		boost::algorithm::split_regex(fields, alt, boost::regex("::"));
		if (fields.size() != 2)
			MTHROW_AND_ERR("Cannot parse alteration string [%s] \n", alt.c_str());
		vector<string> res;
		boost::find_all(res, json_contents, fields[0]);
		if (res.size() > 0)
			MLOG_D("[%s]*%d -> [%s] ", fields[0].c_str(), res.size(), fields[1].c_str());
		else
			MLOG_D("[%s]*%d ", fields[0].c_str(), res.size());
		boost::replace_all(json_contents, fields[0], fields[1]);
	}
	MLOG_D("\n");
}

void MedModel::insert_environment_params_to_json(string& json_content) {

	boost::regex expr{ "ENV\\{(\\S+)\\}" };
	boost::sregex_iterator it(json_content.begin(), json_content.end(), expr);
	boost::sregex_iterator end;

	map<string, string> mapping;
	for (; it != end; ++it) {
		string name = it->str(1);
		char *val = getenv(name.c_str());
		mapping["ENV{" + name + "}"] = string(val);
	}

	for (auto& rec : mapping) {
		vector<string> res;
		boost::find_all(res, json_content, rec.first);
		if (res.size() > 0)
			MLOG_D("[%s]*%d -> [%s] ", rec.first.c_str(), res.size(), rec.second.c_str());
		else
			MLOG_D("[%s]*%d ", rec.first.c_str(), res.size());
		boost::replace_all(json_content, rec.first, rec.second);
	}
	MLOG_D("\n");

}

string MedModel::json_file_to_string(int recursion_level, const string& main_file, vector<string>& alterations,
	const string& small_file, bool add_change_path) {
	if (recursion_level > 3)
		MTHROW_AND_ERR("main file [%s] referenced file [%s], recusion_level 3 reached", main_file.c_str(), small_file.c_str());
	string fname;
	if (small_file == "")
		fname = main_file;
	else
		fname = make_absolute_path(main_file, small_file, false);

	ifstream inf(fname);
	if (!inf)
		MTHROW_AND_ERR("can't open json file [%s] for read\n", fname.c_str());
	stringstream sstr;
	sstr << inf.rdbuf();
	inf.close();
	string orig = stripComments(sstr.str());
	alter_json(orig, alterations);
	insert_environment_params_to_json(orig);

	const char* pattern = "\\\"[[:blank:]]*json\\:(.+?)[[:blank:]]*?\\\"";
	boost::regex ip_regex(pattern);

	boost::sregex_iterator it(orig.begin(), orig.end(), ip_regex);
	boost::sregex_iterator end;
	int last_char = 0;
	string out_string = "";
	string add_path;
	char buff[5000];
	for (; it != end; ++it) {
		string json_ref = it->str(1);
		if (!small_file.empty())
			MLOG_D("Json : found %s, parent %s\n", json_ref.c_str(), small_file.c_str());
		else
			MLOG_D("Json : found %s\n", json_ref.c_str());
		vector<string> tokens;
		boost::split(tokens, json_ref, boost::is_any_of(";"));
		if (tokens.empty())
			MTHROW_AND_ERR("could not parse [%s]", it->str(0).c_str());
		string small_file_inc = tokens[0];
		vector<string> my_alterations;
		for (int i = 1; i < tokens.size(); i++)
			my_alterations.push_back(tokens[i]);
		for (string alt : alterations) {
			vector<string> fields;
			boost::algorithm::split_regex(fields, alt, boost::regex("::"));
			if (fields.size() != 2)
				MTHROW_AND_ERR("Cannot parse alteration string [%s] \n", alt.c_str());
			bool overriden = false;
			for (string existing_alt : my_alterations) {
				vector<string> existing_fields;
				boost::algorithm::split_regex(existing_fields, existing_alt, boost::regex("::"));
				if (fields[0] == existing_fields[0]) {
					MLOG_D("alteration [%s] overriden in the context of [%s] to [%s]\n",
						fields[0].c_str(), small_file_inc.c_str(), existing_fields[1].c_str());
					overriden = true;
				}
			}
			if (!overriden)
				my_alterations.push_back(alt);
		}
		out_string += orig.substr(last_char, it->position() - last_char);
		if (add_change_path) {
			boost::filesystem::path json_p(small_file_inc);
			boost::filesystem::path json_par(fname);
			string pth = boost::filesystem::absolute(json_p.parent_path(), json_par.parent_path()).string();
			snprintf(buff, sizeof(buff), "{\"action_type\":\"change_path:%s\"},\n", pth.c_str());
			add_path = string(buff);
			out_string += add_path;
		}
		out_string += json_file_to_string(recursion_level + 1, main_file, my_alterations, small_file_inc, add_change_path);
		if (add_change_path) {
			snprintf(buff, sizeof(buff), "\n,{\"action_type\":\"change_path:%s\"}\n", run_current_path.c_str());
			add_path = string(buff);
			out_string += add_path;
		}
		last_char = (int)it->position() + (int)it->str(0).size();
	}
	out_string += orig.substr(last_char);
	return out_string;
}

void MedModel::init_from_json_file_with_alterations_version_1(const string &fname, vector<string>& alterations) {
	MWARN("USING DEPRECATED MODEL JSON VERSION 1, PLEASE UPGRADE TO model_json_version: 2\n");
	string json_contents = json_file_to_string(0, fname, alterations);
	istringstream no_comments_stream(json_contents);

	MLOG("MedModel:: init model from json file [%s]\n", fname.c_str());

	ptree pt;
	read_json(no_comments_stream, pt);
	string ser = pt.get<string>("serialize_learning_set", to_string(this->serialize_learning_set).c_str());
	this->serialize_learning_set = stoi(ser);

	for (ptree::value_type &p : pt.get_child("processes"))
	{
		int process_set = -1;
		int duplicate = 1;
		vector<vector<string>> all_attr_values;
		for (ptree::value_type &attr : p.second) {
			string attr_name = attr.first;
			string single_attr_value = attr.second.data();
			if (attr_name == "process_set")
				process_set = stoi(single_attr_value);
			else if (attr_name == "duplicate") {
				boost::algorithm::to_lower(single_attr_value);
				if (single_attr_value == "no" || single_attr_value == "n" || single_attr_value == "0")
					duplicate = 0;
				else if (single_attr_value != "yes" && single_attr_value != "y" && single_attr_value != "1")
					MWARN("NOTE: cannot parse duplicate information \'%s\'. Ignoring\n", single_attr_value.c_str());
			}
			else {
				vector<string> current_attr_values;
				if (single_attr_value.length() > 0) {
					//MLOG("attr %s %s\n", attr_name.c_str(), single_attr_value.c_str());
					if (boost::starts_with(single_attr_value, "file:")) {
						//e.g. "signal": "file:my_list.txt" - file can be relative
						vector<string> my_list;
						string small_file = single_attr_value.substr(5);
						fill_list_from_file(make_absolute_path(fname, small_file), my_list);
						for (string s : my_list)
							current_attr_values.push_back(parse_key_val(attr_name, s));
					}
					else if (boost::starts_with(single_attr_value, "ref:")) {
						auto my_ref = pt.get_child(single_attr_value.substr(4));
						for (auto &r : my_ref)
							//e.g. "signal": "ref:signals"
							current_attr_values.push_back(parse_key_val(attr_name, r.second.data()));
					}
					else
						// e.g. "fg_type": "gender"
						current_attr_values.push_back(parse_key_val(attr_name, single_attr_value));
				}
				else
					//e.g. "type": ["last", "slope"]
					for (ptree::value_type &attr_value : attr.second)
						current_attr_values.push_back(parse_key_val(attr_name, attr_value.second.data()));
				all_attr_values.push_back(current_attr_values);
			}
		}

		vector<string> all_combinations;
		concatAllCombinations(all_attr_values, 0, "", all_combinations);
		for (string c : all_combinations) {
			//MLOG("MedModel::init [%s]\n", c.c_str());
			add_process_to_set(process_set, duplicate, c);
		}
	}
	if (pt.count("predictor") > 0) {
		auto my_pred = pt.get_child("predictor");
		auto my_pred_params = pt.get_child("predictor_params");
		set_predictor(my_pred.data(), my_pred_params.data());
	}
	else MWARN("NOTE: no [predictor] node found in file\n");

}

// inserting a new rep-processor from string at a given position
//.......................................................................................
void MedModel::insert_rep_processor(string init_string, int idx)
{
	RepProcessor *rep_proc = RepProcessor::create_processor(init_string);
	rep_processors.insert(rep_processors.begin() + idx, rep_proc);
}

// generalized adder

// generalized adder
// type and signal are must have parameters in this case
//.......................................................................................
void MedModel::add_rep_processor_to_set(int i_set, const string &init_string)
{
	// check if i_set already initialized, and if not a multiprocessor change it into one
	if (i_set < rep_processors.size()) {
		// exists 
		if (rep_processors[i_set] == NULL) {
			// NULL ... in that case init an empty MultiProcessor in i_set
			MLOG_D("Adding new rep_processor set [%d]\n", i_set);
			RepMultiProcessor *processor = new RepMultiProcessor;
			rep_processors[i_set] = processor;
		}
		else if (rep_processors[i_set]->processor_type != REP_PROCESS_MULTI) {
			// the processor was not multi, and hence we create one switch it , and push the current into it
			RepProcessor *curr_p = rep_processors[i_set];
			RepMultiProcessor *mprocessor = new RepMultiProcessor;
			rep_processors[i_set] = mprocessor;
			mprocessor->processors.push_back(curr_p);
		}
	}
	else {
		// resize rep_processors
		rep_processors.resize(i_set + 1, NULL);
		for (int i = 0; i < i_set + 1; i++)
			// put a new empty multi in i_set
			if (rep_processors[i] == NULL) {
				MLOG_D("Adding new rep_processor set [%d]\n", i);
				RepMultiProcessor *processor = new RepMultiProcessor;
				rep_processors[i] = processor;
			}
	}

	// Now we are at a state in which we have a multi processor at i_set and need to create a new processor and push it in
	string in = init_string;
	RepProcessor *rep_proc = RepProcessor::create_processor(in);

	// push it in
	((RepMultiProcessor *)rep_processors[i_set])->processors.push_back(rep_proc);

}

//.......................................................................................
// fp_type and feature name are must have parameters
void MedModel::add_feature_processor_to_set(int i_set, int duplicate, const string &init_string)
{
	// if init_string does not have a names list (names parameter empty) it means a feature processor
	// will be added to each of the currently initialized features.
	// This means that this is order dependent.
	// One has also to be careful not to enter the same feature twice.

	// check if i_set already initialized, and if not a multiprocessor change it into one
	if (i_set < feature_processors.size()) {
		// exists 
		if (feature_processors[i_set] == NULL) {
			// NULL ... in that case init an empty MultiProcessor in i_set
			MLOG_D("Adding new feature_processor set [%d]\n", i_set);
			MultiFeatureProcessor *mfprocessor = new MultiFeatureProcessor;
			if (mfprocessor->init_from_string(init_string) < 0)
				MTHROW_AND_ERR("Cannot init MultiFeatureProcessor  with init string \'%s\'\n", init_string.c_str());
			feature_processors[i_set] = mfprocessor;
		}
		else if (feature_processors[i_set]->processor_type != FTR_PROCESS_MULTI) {
			// the processor was not multi, and hence we create one switch it , and push the current into it
			FeatureProcessor *curr_fp = feature_processors[i_set];
			MultiFeatureProcessor *mfprocessor = new MultiFeatureProcessor;
			if (mfprocessor->init_from_string(init_string) < 0)
				MTHROW_AND_ERR("Cannot init MultiFeatureProcessor  with init string \'%s\'\n", init_string.c_str());
			feature_processors[i_set] = mfprocessor;
			mfprocessor->processors.push_back(curr_fp);

		}
	}
	else {
		// resize feature_processors
		feature_processors.resize(i_set + 1, NULL);

		for (int i = 0; i < i_set + 1; i++)
			// put a new empty multi in i_set
			if (feature_processors[i] == NULL) {
				MLOG_D("Adding new feature_processor set [%d]\n", i);
				MultiFeatureProcessor *mfprocessor = new MultiFeatureProcessor;
				if (mfprocessor->init_from_string(init_string) < 0)
					MTHROW_AND_ERR("Cannot init MultiFeatureProcessor  with init string \'%s\'\n", init_string.c_str());
				feature_processors[i] = mfprocessor;
			}
	}

	// Now we are at a state in which we have a multi feature processor at i_set and need to create a new processor or processors and push it in

	// get all relevant features names
	string feat_names;
	get_single_val_from_init_string(init_string, "names", feat_names);

	// get type of feature processor
	string fp_type;
	get_single_val_from_init_string(init_string, "fp_type", fp_type);
	FeatureProcessorTypes type = feature_processor_name_to_type(fp_type);


	if (feat_names != "" && feat_names != "All") { // Are features given ?
		vector<string> features;
		boost::split(features, feat_names, boost::is_any_of(","));
		MLOG("fp_type [%s] acting on [%d] features\n", fp_type.c_str(), int(features.size()));
		((MultiFeatureProcessor *)feature_processors[i_set])->add_processors_set(type, features, init_string);
	}
	else if (feat_names == "All" || (feat_names == "" && duplicate)) { // Work on all features. Will be created at Learn and will adhere to "tag" filtering
		((MultiFeatureProcessor *)feature_processors[i_set])->init_string = init_string;
		((MultiFeatureProcessor *)feature_processors[i_set])->members_type = type;
		((MultiFeatureProcessor *)feature_processors[i_set])->duplicate = 1;
	}
	else { // No duplicating and no feature name given (e.g. selector)
		FeatureProcessor *processor = FeatureProcessor::make_processor(type, init_string);
		((MultiFeatureProcessor *)feature_processors[i_set])->processors.push_back(processor);
	}

}

//.......................................................................................
void MedModel::add_feature_generator_to_set(int i_set, const string &init_string)
{
	// currently there's NO multi feature generator .... (TBD)
	// hence currently we simply ignore i_set, and pile up generators into generators

	string in = init_string;
	FeatureGenerator *feat_gen = FeatureGenerator::create_generator(in);

	// push it in
	generators.push_back(feat_gen);
}

void MedModel::add_post_processor_to_set(int i_set, const string &init_string)
{
	// check if i_set already initialized, and if not a multiprocessor change it into one
	if (i_set < post_processors.size()) {
		// exists 
		if (post_processors[i_set] == NULL) {
			// NULL ... in that case init an empty MultiProcessor in i_set
			MLOG_D("Adding new rep_processor set [%d]\n", i_set);
			MultiPostProcessor *processor = new MultiPostProcessor;
			post_processors[i_set] = processor;
		}
		else if (post_processors[i_set]->processor_type != FTR_POSTPROCESS_MULTI) {
			// the processor was not multi, and hence we create one switch it , and push the current into it
			PostProcessor *curr_p = post_processors[i_set];
			MultiPostProcessor *mprocessor = new MultiPostProcessor;
			post_processors[i_set] = mprocessor;
			mprocessor->post_processors.push_back(curr_p);
		}
	}
	else {
		// resize rep_processors
		post_processors.resize(i_set + 1, NULL);
		for (int i = 0; i < i_set + 1; i++)
			// put a new empty multi in i_set
			if (post_processors[i] == NULL) {
				MLOG_D("Adding new post_processor set [%d]\n", i);
				MultiPostProcessor *processor = new MultiPostProcessor;
				post_processors[i] = processor;
			}
	}

	// Now we are at a state in which we have a multi processor at i_set and need to create a new processor and push it in
	string in = init_string;
	string pp_type;
	get_single_val_from_init_string(in, "pp_type", pp_type);
	PostProcessor *post_proc = PostProcessor::make_processor(pp_type, init_string);

	// push it in
	((MultiPostProcessor *)post_processors[i_set])->post_processors.push_back(post_proc);
}

//.......................................................................................
void MedModel::add_process_to_set(int i_set, int duplicate, const string &init_string)
{
	if (init_string.find("rp_type") != string::npos) return add_rep_processor_to_set(i_set, init_string);
	if (init_string.find("fg_type") != string::npos) return add_feature_generator_to_set(i_set, init_string);
	if (init_string.find("fp_type") != string::npos) return add_feature_processor_to_set(i_set, duplicate, init_string);
	if (init_string.find("pp_type") != string::npos) return add_post_processor_to_set(i_set, init_string);

	MTHROW_AND_ERR("add_process_to_set():: Can't process line %s\n", init_string.c_str());
}


// Add multi processors
//.......................................................................................
void MedModel::add_rep_processors_set(RepProcessorTypes type, vector<string>& signals) {

	RepMultiProcessor *processor = new RepMultiProcessor;
	processor->add_processors_set(type, signals);
	add_rep_processor(processor);

}

// Affected Signals
//.......................................................................................
void MedModel::set_affected_signal_ids(MedDictionarySections& dict) {

	for (RepProcessor *processor : rep_processors)
		processor->set_affected_signal_ids(dict);
}

// initialization :  find signal ids, init tables
//.......................................................................................
void MedModel::init_all(MedDictionarySections& dict, MedSignals& sigs) {

	// signal ids
	set_affected_signal_ids(dict);
	set_required_signal_ids(dict);

	for (RepProcessor *processor : rep_processors)
		processor->set_signal_ids(sigs);

	for (FeatureGenerator *generator : generators)
		generator->set_signal_ids(sigs);

	// tables
	for (RepProcessor *processor : rep_processors)
		processor->init_tables(dict, sigs);

	for (FeatureGenerator *generator : generators)
		generator->init_tables(dict);

	// attributes
	for (RepProcessor *processor : rep_processors)
		processor->init_attributes();
}

// Create a required signal names set by back propograting : First find what's required by
// the feature generators, and then find add signals required by the rep_porcessors that
// are required ....
//.......................................................................................
void MedModel::get_required_signal_names(unordered_set<string>& signalNames) {

	// Identify required generators
	vector<unordered_set<string> > req_features_vec;
	build_req_features_vec(req_features_vec);
	unordered_set<string>& req_feature_generators = req_features_vec[feature_processors.size()];

	vector<FeatureGenerator *> applied_generators;
	get_applied_generators(req_feature_generators, applied_generators);

	// Identify required signals
	get_all_required_signal_names(signalNames, rep_processors, -1, applied_generators);

	// collect virtuals
	for (RepProcessor *processor : rep_processors) {
		if (verbosity) MLOG_D("MedModel::get_required_signal_names adding virtual signals from rep type %d\n", processor->processor_type);
		processor->add_virtual_signals(virtual_signals, virtual_signals_generic);
	}

	if (verbosity) MLOG_D("MedModel::get_required_signal_names %d signalNames %d virtual_signals\n", signalNames.size(), virtual_signals.size());


	// Erasing virtual signals !
	for (auto &vsig : virtual_signals) {
		if (verbosity) MLOG_D("check virtual %s\n", vsig.first.c_str());
		if (signalNames.find(vsig.first) != signalNames.end())
			signalNames.erase(vsig.first);
	}

	// Erasing virtual signals !
	for (auto &vsig : virtual_signals_generic) {
		if (verbosity) MLOG_D("check virtual %s\n", vsig.first.c_str());
		if (signalNames.find(vsig.first) != signalNames.end())
			signalNames.erase(vsig.first);
	}

	if (verbosity) MLOG_D("MedModel::get_required_signal_names %d signalNames %d virtual_signals after erasing\n", signalNames.size(), virtual_signals.size() + virtual_signals_generic.size());

}

// Get required names as a vector
//.......................................................................................
void MedModel::get_required_signal_names(vector<string>& signalNames) {
	unordered_set<string> sigs;
	get_required_signal_names(sigs);
	signalNames.clear();
	for (auto &s : sigs)
		signalNames.push_back(s);
}

// Collect required signal names to generate processed values of signals in target-set
//.......................................................................................
void MedModel::get_required_signal_names_for_processed_values(unordered_set<string>& targetSignalNames, unordered_set<string>& signalNames) {

	signalNames = targetSignalNames;

	// Collect from processors itertively
	for (int i = (int)rep_processors.size() - 1; i >= 0; i--) {
		rep_processors[i]->get_required_signal_names(signalNames, signalNames);
	}

	// collect virtuals
	for (RepProcessor *processor : rep_processors) {
		if (verbosity) MLOG_D("MedModel::get_required_signal_names adding virtual signals from rep type %d\n", processor->processor_type);
		processor->add_virtual_signals(virtual_signals, virtual_signals_generic);
	}

	if (verbosity) MLOG_D("MedModel::get_required_signal_names %d signalNames %d virtual_signals\n", signalNames.size(), virtual_signals.size());

	// Erasing virtual signals !
	for (auto &vsig : virtual_signals) {
		if (verbosity) MLOG("check virtual %s\n", vsig.first.c_str());
		if (signalNames.find(vsig.first) != signalNames.end())
			signalNames.erase(vsig.first);
	}

	for (auto &vsig : virtual_signals_generic) {
		if (verbosity) MLOG("check virtual %s\n", vsig.first.c_str());
		if (signalNames.find(vsig.first) != signalNames.end())
			signalNames.erase(vsig.first);
	}

	if (verbosity) MLOG_D("MedModel::get_required_signal_names %d signalNames %d virtual_signals after erasing\n", signalNames.size(), virtual_signals.size() + virtual_signals_generic.size());

}

// Get required names as a vector
//.......................................................................................
void MedModel::get_required_signal_names_for_processed_values(unordered_set<string>& targetSignalNames, vector<string>& signalNames) {

	unordered_set<string> sigs;
	get_required_signal_names_for_processed_values(targetSignalNames, sigs);

	signalNames.clear();
	for (auto &s : sigs)
		signalNames.push_back(s);
}

//.......................................................................................
void collect_and_add_virtual_signals_static(MedRepository &rep, vector<RepProcessor *> &rep_processors,
	map<string, int> *virtual_signals = NULL, map<string, string> *virtual_signals_generic = NULL, bool verbosity = true)
{
	map<string, int> temp_m;
	map<string, string> temp_mg;
	if (virtual_signals == NULL)
		virtual_signals = &temp_m;
	if (virtual_signals_generic == NULL)
		virtual_signals_generic = &temp_mg;
	// collecting
	for (RepProcessor *processor : rep_processors)
		processor->add_virtual_signals(*virtual_signals, *virtual_signals_generic);
	//register section_name to section_id if needed in each rep_processor virtual signal
	for (RepProcessor *processor : rep_processors)
		processor->register_virtual_section_name_id(rep.dict);

	vector<string> all_virtual_signals_names;

	// adding to rep
	for (auto &vsig : *virtual_signals) {
		//MLOG("Attempting to add virtual signal %s type %d (%d)\n", vsig.first.c_str(), vsig.second, rep.sigs.sid(vsig.first));
		if (rep.sigs.sid(vsig.first) < 0) {
			int new_id = rep.sigs.insert_virtual_signal(vsig.first, vsig.second);
			all_virtual_signals_names.push_back(vsig.first);
			if (verbosity > 0)
				MLOG_D("Added Virtual Signal %s type %d : got id %d\n", vsig.first.c_str(), vsig.second, new_id);
		}
		else {
			if (rep.sigs.sid(vsig.first) < 100)
				MTHROW_AND_ERR("Failed defining virtual signal %s (type %d)...(curr sid for it is: %d)\n", vsig.first.c_str(), vsig.second, rep.sigs.sid(vsig.first));
		}
	}

	// adding generic virual signals to rep
	for (auto &vsig : *virtual_signals_generic) {
		if (rep.sigs.sid(vsig.first) < 0) {
			int new_id = rep.sigs.insert_virtual_signal(vsig.first, vsig.second);
			all_virtual_signals_names.push_back(vsig.first);
			if (verbosity > 0)
				MLOG_D("Added Generic Virtual Signal %s spec='%s' : got id %d\n", vsig.first.c_str(), vsig.second.c_str(), new_id);
		}
		else {
			if (rep.sigs.sid(vsig.first) < 100)
				MTHROW_AND_ERR("Failed defining virtual signal %s (spec=%s)...(curr sid for it is: %d)\n", vsig.first.c_str(), vsig.second.c_str(), rep.sigs.sid(vsig.first));
		}
	}

	// update dict
	for (const auto &vsig_name : all_virtual_signals_names) {
		int vsig_id = rep.sigs.Name2Sid[vsig_name];
		int add_section = rep.dict.section_id(vsig_name);
		rep.dict.dicts[add_section].Name2Id[vsig_name] = vsig_id;
		rep.dict.dicts[0].Name2Id[vsig_name] = vsig_id;
		rep.dict.dicts[add_section].Id2Name[vsig_id] = vsig_name;
		rep.dict.dicts[add_section].Id2Names[vsig_id] = { vsig_name };
		rep.sigs.Sid2Info[vsig_id].time_unit = rep.sigs.my_repo->time_unit;
		//rep.dict.SectionName2Id[vsig_name] = 0;
		MLOG_D("updated dict %d : %d\n", add_section, rep.dict.dicts[add_section].id(vsig_name));
	}

}

int MedModel::collect_and_add_virtual_signals(MedRepository &rep)
{
	collect_and_add_virtual_signals_static(rep, rep_processors, &virtual_signals, &virtual_signals_generic, verbosity);
	return 0;
}

//.......................................................................................
void  MedModel::add_rep_processors_set(RepProcessorTypes type, vector<string>& signals, string init_string) {

	RepMultiProcessor *processor = new RepMultiProcessor;
	processor->add_processors_set(type, signals, init_string);
	add_rep_processor(processor);

}

//.......................................................................................
void MedModel::add_feature_processors_set(FeatureProcessorTypes type) {

	vector<string> features;
	get_all_features_names(features, int(feature_processors.size()));

	//	for (auto &name : features) 
	//		MLOG("Adding %s to processors of type %d\n", name.c_str(),type);

	add_feature_processors_set(type, features);

}

//.......................................................................................
void MedModel::add_feature_processors_set(FeatureProcessorTypes type, string init_string) {

	vector<string> features;
	get_all_features_names(features, int(feature_processors.size()));

	add_feature_processors_set(type, features, init_string);
}
//.......................................................................................
void MedModel::add_feature_processors_set(FeatureProcessorTypes type, vector<string>& features) {

	MultiFeatureProcessor *fProcessor = new MultiFeatureProcessor;
	fProcessor->add_processors_set(type, features);
	add_feature_processor(fProcessor);

}

//.......................................................................................
void MedModel::add_feature_processors_set(FeatureProcessorTypes type, vector<string>& features, string init_string) {

	MultiFeatureProcessor *fProcessor = new MultiFeatureProcessor;
	fProcessor->add_processors_set(type, features, init_string);
	add_feature_processor(fProcessor);

}

// Add sets of generators
//.......................................................................................
void MedModel::add_feature_generators(FeatureGeneratorTypes type, vector<string>& signals) {

	for (string& signal : signals) {
		FeatureGenerator *generator = FeatureGenerator::make_generator(type, "signalName=" + signal);
		add_feature_generator(generator);
	}
}

//.......................................................................................
void MedModel::add_feature_generators(FeatureGeneratorTypes type, vector<string>& signals, string init_string) {

	for (string& signal : signals) {
		FeatureGenerator *generator;
		if (signal != "")
			generator = FeatureGenerator::make_generator(type, init_string + ";signalName=" + signal);
		else
			generator = FeatureGenerator::make_generator(type, init_string);
		add_feature_generator(generator);
	}

	if (signals.size() == 0) {
		FeatureGenerator *generator = FeatureGenerator::make_generator(type, init_string);
		add_feature_generator(generator);
	}
}

//.......................................................................................
void MedModel::get_all_features_names(vector<string> &feat_names, int before_process_set)
{
	unordered_set<string> uniq_feat_names;
	for (unsigned int i = 0; i < generators.size(); i++) {
		for (string& name : generators[i]->names)
			uniq_feat_names.insert(name);
	}
	assert((before_process_set <= feature_processors.size()) || (before_process_set == 0 && feature_processors.size() == 0));
	for (int i = 0; i < before_process_set; i++) {
		//if (feature_processors[i]->is_selector()) {
		//	FeatureSelector *fs = (FeatureSelector *)feature_processors[i];
		//	uniq_feat_names = fs->selected;
		//}
		vector<string> my_names;
		feature_processors[i]->get_feature_names(my_names);
		for (string& name : my_names)
			uniq_feat_names.insert(name);
	}
	feat_names.clear();
	for (auto f : uniq_feat_names)
		feat_names.push_back(f);
}


void MedModel::clear()
{
	if (rep_processors.size() > 0) {
		for (auto prep : rep_processors)
			if (prep != NULL) {
				delete prep;
				prep = NULL;
			}
		rep_processors.clear();
	}

	if (generators.size() > 0) {
		for (auto pgen : generators)
			if (pgen != NULL) {
				delete pgen;
				pgen = NULL;
			}
		generators.clear();
	}

	if (feature_processors.size() > 0) {
		for (auto pfeat : feature_processors)
			if (pfeat != NULL) {
				delete pfeat;
				pfeat = NULL;
			}
		feature_processors.clear();
	}

	if (post_processors.size() > 0) {
		for (auto postprocc : post_processors)
			if (postprocc != NULL) {
				delete postprocc;
				postprocc = NULL;
			}
		post_processors.clear();
	}

	if (predictor != NULL) {
		delete predictor;
		predictor = NULL;
	}

	//if (LearningSet != NULL) {
	//	delete LearningSet;
	//	LearningSet = NULL;
	//}
}

//.......................................................................................
void MedModel::dprint_process(const string &pref, int rp_flag, int fg_flag, int fp_flag, int predictor_flag, int pp_flag)
{
	unordered_set<string> sigs;

	get_required_signal_names(sigs);
	MLOG("%s : MedModel with rp(%d) fg(%d) fp(%d)\n", pref.c_str(), rep_processors.size(), generators.size(), feature_processors.size());
	MLOG("%s : MedModel with required_signal_names(%d) : ", pref.c_str(), sigs.size());
	for (auto& s : sigs) MLOG("%s,", s.c_str());
	MLOG("\n");
	if (rp_flag > 0) for (auto& rp : rep_processors) rp->dprint(pref, rp_flag);
	if (fg_flag > 0) for (auto& fg : generators) fg->dprint(pref, fg_flag);
	if (fp_flag > 0) for (auto& fp : feature_processors) fp->dprint(pref, fp_flag);
	if (predictor_flag > 0 && predictor != NULL) predictor->print(stderr, pref, predictor_flag);
	if (pp_flag > 0) for (auto& pp : post_processors) pp->dprint(pref);
}

void MedModel::get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const {
	for (auto& proc : rep_processors) {
		unordered_map<string, vector<string>> local_use;
		proc->get_required_signal_categories(local_use);
		//do merge with signal_categories_in_use:
		for (const auto &it : local_use)
		{
			if (signal_categories_in_use.find(it.first) == signal_categories_in_use.end())
				signal_categories_in_use[it.first] = move(it.second);
			else {
				//merge with existing:
				unordered_set<string> existing_sets(signal_categories_in_use.at(it.first).begin(),
					signal_categories_in_use.at(it.first).end());
				existing_sets.insert(it.second.begin(), it.second.end());;
				vector<string> uniq_vec(existing_sets.begin(), existing_sets.end());
				signal_categories_in_use[it.first] = move(uniq_vec);
			}
		}
	}

	for (const auto& proc : generators) {
		unordered_map<string, vector<string>> local_use;
		proc->get_required_signal_categories(local_use);
		//do merge with signal_categories_in_use:
		for (auto &it : local_use)
		{
			if (signal_categories_in_use.find(it.first) == signal_categories_in_use.end())
				signal_categories_in_use[it.first] = move(it.second);
			else {
				//merge with existing:
				unordered_set<string> existing_sets(signal_categories_in_use.at(it.first).begin(),
					signal_categories_in_use.at(it.first).end());
				existing_sets.insert(it.second.begin(), it.second.end());;
				vector<string> uniq_vec(existing_sets.begin(), existing_sets.end());
				signal_categories_in_use[it.first] = move(uniq_vec);
			}
		}
	}
}

//.......................................................................................
void filter_rep_processors(const vector<string> &current_req_signal_names, vector<RepProcessor *> *rep_processors) {
	unordered_set<string> req_signal_names(current_req_signal_names.begin(), current_req_signal_names.end());
	vector<RepProcessor *> filtered_processors;
	bool did_something = false;

	for (unsigned int i = 0; i < rep_processors->size(); i++) {
		unordered_set<string> current_req_signal_names(req_signal_names.begin(), req_signal_names.end());
		for (int k = (int)rep_processors->size() - 1; k > i; k--)
			rep_processors->at(k)->get_required_signal_names(current_req_signal_names, current_req_signal_names);

		if (!(*rep_processors)[i]->filter(current_req_signal_names))
			filtered_processors.push_back((*rep_processors)[i]);
		else {//cleaning uneeded rep_processors!:
			delete (*rep_processors)[i];
			(*rep_processors)[i] = NULL;
			did_something = true;
		}
	}
	if (did_something)
		MLOG("Filtering uneeded rep_processors. keeping %zu rep_proccessors out of %zu\n",
			filtered_processors.size(), rep_processors->size());

	rep_processors->swap(filtered_processors);
}

void medial::repository::prepare_repository(const vector<int> &pids, const string &RepositoryPath,
	MedModel &mod, MedPidRepository &rep) {
	MLOG("Reading repo file [%s]\n", RepositoryPath.c_str());
	unordered_set<string> req_names;
	if (rep.init(RepositoryPath) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", RepositoryPath.c_str());
	mod.fit_for_repository(rep);
	mod.filter_rep_processors();

	mod.get_required_signal_names(req_names);

	vector<string> sigs = { "BDATE", "GENDER", "TRAIN" };
	for (string s : req_names)
		sigs.push_back(s);
	sort(sigs.begin(), sigs.end());
	auto it = unique(sigs.begin(), sigs.end());
	sigs.resize(std::distance(sigs.begin(), it));

	if (rep.read_all(RepositoryPath, pids, sigs) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", RepositoryPath.c_str());
}

void medial::repository::prepare_repository(const MedSamples &samples, const string &RepositoryPath,
	MedModel &mod, MedPidRepository &rep) {
	vector<int> pids;
	samples.get_ids(pids);
	prepare_repository(pids, RepositoryPath, mod, rep);
}

vector<string> medial::repository::prepare_repository(MedPidRepository &rep, const vector<string> &needed_sigs,
	vector<string> &phisical_signal_read, vector<RepProcessor *> *rep_processors) {

	vector<unordered_set<string>> current_req_signal_names;
	if (rep_processors != NULL && !rep_processors->empty()) {
		for (RepProcessor *processor : *rep_processors)
			processor->fit_for_repository(rep);
		collect_and_add_virtual_signals_static(rep, *rep_processors);
		//init to check if need to remove (may seem it can remove after init)
		filter_rep_processors(needed_sigs, rep_processors);
		for (RepProcessor *processor : *rep_processors) {
			processor->set_affected_signal_ids(rep.dict);
			processor->set_signal_ids(rep.sigs);
			processor->set_required_signal_ids(rep.dict);
			processor->init_tables(rep.dict, rep.sigs);
			processor->init_attributes();
		}

		//vector<RepProcessor *> temp_processors;
		for (int i = 0; i < rep_processors->size(); i++) {
			//unordered_set<int> current_req_signal_ids;
			//for (int k = (int)rep_processors->size() - 1; k > i; --k)
			//	(*rep_processors)[i]->get_required_signal_ids(current_req_signal_ids, current_req_signal_ids);
			if ((*rep_processors)[i]->learn(rep) < 0)
				MTHROW_AND_ERR("Unable to learn rep_processor\n");
			//temp_processors.push_back((*rep_processors)[i]);
		}

		current_req_signal_names.resize(rep_processors->size());
		for (unsigned int i = 0; i < rep_processors->size(); i++)
			(*rep_processors)[i]->get_required_signal_names(current_req_signal_names[i]);
	}
	unordered_set<string> all_rep_sigs(needed_sigs.begin(), needed_sigs.end());
	for (size_t i = 0; i < current_req_signal_names.size(); ++i)
		all_rep_sigs.insert(current_req_signal_names[i].begin(), current_req_signal_names[i].end());
	vector<string> final_use_sigs(all_rep_sigs.begin(), all_rep_sigs.end());
	phisical_signal_read.clear();
	phisical_signal_read.reserve(final_use_sigs.size());
	for (size_t i = 0; i < final_use_sigs.size(); ++i)
	{
		int sid = rep.sigs.sid(final_use_sigs[i]);
		if (sid < 0)
			MWARN("Warning: signal %s not found in prepared repository.\n", final_use_sigs[i].c_str());
		if (sid >= 0 && !rep.sigs.Sid2Info[sid].virtual_sig)
			phisical_signal_read.push_back(final_use_sigs[i]);
	}

	return final_use_sigs;
}

vector<int> medial::repository::prepare_repository(MedPidRepository &rep, const vector<int> &needed_sigs,
	vector<int> &phisical_signal_read, vector<RepProcessor *> *rep_processors) {
	vector<string> conv_in(needed_sigs.size());
	for (size_t i = 0; i < needed_sigs.size(); ++i)
		conv_in[i] = rep.sigs.name(needed_sigs[i]);
	vector<string> phisical_sigs;
	vector<string> res = prepare_repository(rep, conv_in, phisical_sigs, rep_processors);
	vector<int> conv_res(res.size());
	for (size_t i = 0; i < conv_res.size(); ++i)
		conv_res[i] = rep.sigs.sid(res[i]);
	phisical_signal_read.resize(phisical_sigs.size());
	for (size_t i = 0; i < phisical_sigs.size(); ++i)
		phisical_signal_read[i] = rep.sigs.sid(phisical_sigs[i]);

	return conv_res;
}

//.......................................................................................
int MedModel::write_feature_matrix(const string mat_fname)
{
	MedMat<float> x;

	features.get_as_matrix(x);

	return x.write_to_csv_file(mat_fname);
}

#if 1
//...................................................................................................
int MedModel::init_for_apply_rec(MedPidRepository &rep)
{
	// init virtual signals
	if (collect_and_add_virtual_signals(rep) < 0) {
		MERR("FAILED collect_and_add_virtual_signals\n");
		return -1;
	}

	// Initialize
	init_all(rep.dict, rep.sigs);

	// Required signals
	required_signal_names.clear();
	required_signal_ids.clear();

	get_required_signal_names(required_signal_names);
	for (string signal : required_signal_names)
		required_signal_ids.insert(rep.dict.id(signal));

	for (int signalId : required_signal_ids) {
		if ((!rep.in_mem_mode_active()) && rep.index.index_table[signalId].is_loaded != 1)
			MLOG("MedModel::apply WARNING signal [%d] = [%s] is required by model but not loaded in rep\n",
				signalId, rep.dict.name(signalId).c_str());;
	}

	// prepare for generation
	MedSamples samples;
	for (auto& generator : generators) {
		if (generator->generator_type == FTR_GEN_MODEL)
			MTHROW_AND_ERR("MedModel::init_for_apply_rec() does not support FTR_GEN_MODEL type in this version\n");
		generator->prepare(features, rep, samples);
	}

	return 0;
}

//...................................................................................................
// important:
// currently this implementation DOES NOT support weights and runs ALL rep_processors without
// the mechanisms for conditional_apply. Also Does not support attributes at the moment.
// Both should be added in the future.
// Main usage for now : allow most well defined models to be run rec by rec, with a natural way of 
// combining them in other models.
int MedModel::apply_rec(PidDynamicRec &drec, MedIdSamples idSamples, MedFeatures &_feat, bool copy_rec_flag, int end_stage)
{
	int start_stage = MED_MDL_APPLY_FTR_GENERATORS;
	//int end_stage = MED_MDL_APPLY_FTR_PROCESSORS;// MED_MDL_END;

	// Stage Sanity
	if (end_stage < MED_MDL_APPLY_FTR_GENERATORS) {
		MERR("MedModel apply() : Illegal end stage %d\n", end_stage);
		return -1;
	}


	// drec copy if needed
	PidDynamicRec &rec = drec;
	PidDynamicRec copy_rec;
	if (copy_rec_flag) {
		copy_rec = drec;
		copy_rec.set_data_to_buffer();
		rec = copy_rec;
	}


	if (start_stage <= MED_MDL_APPLY_FTR_GENERATORS) {


		// Currently this implementation ignores the savings that can be made using conditional apply, 
		// Mainly for clarity and simplicity.
		// It is to be done later

		// _feat initializations
		// not supporting weights attributes in this code... again for simplicity
		_feat = features; // copy template features that was already initialized
		_feat.set_time_unit(rec.my_base_rep->time_unit);
		_feat.append_samples(idSamples);
		_feat.init_pid_pos_len();

		vector<int> times;
		idSamples.get_times(times);

		// Resize data vectors and collect pointers
		int samples_size = (int)idSamples.samples.size();
		vector<vector<float *>> _pdata(generators.size());
		float jj = 0;
		for (auto &generator : generators) {
			for (string& name : generator->names)
				_feat.data[name].resize(samples_size, jj++);
		}

		int k = 0;
		for (auto &generator : generators)
			generator->get_p_data(_feat, _pdata[k++]);

		// Apply rep-processing
		for (unsigned int i = 0; i < rep_processors.size(); i++)
			if (rep_processors[i]->_apply_simple(rec, times) < 0)
				return -1;

		// Generate Features
		k = 0;
		for (auto& generator : generators) {
			if (generator->_generate(rec, _feat, 0, samples_size, _pdata[k++]) < 0)
				return -1;
		}
	}

	if (end_stage <= MED_MDL_APPLY_FTR_GENERATORS)
		return 0;

	// Process Features
	for (auto& processor : feature_processors) {
		if (processor->apply(_feat, false) < 0)
			return -1;
	}

	if (end_stage <= MED_MDL_APPLY_FTR_PROCESSORS)
		return 0;

	// Apply predictor
	if (predictor->predict(_feat) < 0) {
		MERR("Predictor failed\n");
		return -1;
	}

	return 0;
}

//-----------------------------------------------------------------------------------------------------------------
int MedModel::get_nfeatures()
{
	int res = 0;
	// in this case we collect all feature generator names
	vector<unordered_set<string> > req_features_vec;
	build_req_features_vec(req_features_vec); //from all feature processor

	unordered_set<string> ftr_names = req_features_vec[feature_processors.size()];
	if (ftr_names.empty()) {
		// in this case we collect all feature generator names
		for (FeatureGenerator *generator : generators)
			generator->get_generated_features(ftr_names);
	}

	unordered_set<string> names;
	for (FeatureGenerator *generator : generators) {
		if (generator->names.empty()) //if the generation is dynamic names will be empty - fetch count from nfeatures
			res += generator->nfeatures();
	}
	res += (int)ftr_names.size();

	// next is done in order to return feature list sorted in the order it will be in the final matrix
	return res;
}
//-----------------------------------------------------------------------------------------------------------------
int MedModel::get_duplicate_factor() const {
	int factor = 1;
	//search for DuplicateProcessor (AggregatePredsPostProcessor has a built in batch processing, skip him)
	for (FeatureProcessor *fp : feature_processors)
	{
		if (fp->processor_type == FTR_PROCESS_MULTI) {
			for (FeatureProcessor *multi_fp : static_cast<MultiFeatureProcessor *>(fp)->processors)
				if (multi_fp->processor_type == FTR_PROCESS_DUPLICATE)
					factor *= static_cast<DuplicateProcessor *>(multi_fp)->resample_cnt;
		}
		else if (fp->processor_type == FTR_PROCESS_DUPLICATE)
			factor *= static_cast<DuplicateProcessor *>(fp)->resample_cnt;
	}

	return factor;
}
//-----------------------------------------------------------------------------------------------------------------

void MedModel::get_generated_features_names(vector<string> &feat_names)
{
	vector<unordered_set<string> > req_features_vec;
	build_req_features_vec(req_features_vec);

	unordered_set<string> ftr_names = req_features_vec[feature_processors.size()];
	if (ftr_names.empty()) {
		// in this case we collect all feature generator names
		for (auto &generator : generators)
			generator->get_generated_features(ftr_names);
	}

	// next is done in order to return feature list sorted in the order it will be in the final matrix
	map<string, int> sort_me;
	for (auto &e : ftr_names) sort_me[e] = 1;
	feat_names.clear();
	for (auto &e : sort_me) feat_names.push_back(e.first);
}

// Handle learning sets for model/post-processors
//-----------------------------------------------------------------------------------------------------------------
void MedModel::split_learning_set(MedSamples& inSamples, vector<MedSamples>& post_processors_learning_sets, MedSamples& model_learning_set) {

	int nIds = (int)inSamples.idSamples.size();
	vector<int> assignments(nIds, 0);
	post_processors_learning_sets.resize(post_processors.size());

	// Assign to post-processors
	int idx = 1;
	int nFreeIds = nIds;
	for (PostProcessor *processor : post_processors) {
		float use_p = processor->get_use_p();
		int use_split = processor->get_use_split();
		if (use_p > 0 && use_split >= 0)
			MTHROW_AND_ERR("Split_Learning_Set: At most one of use_p (%f) & use_split (%d) allowed for post-processor\n", use_p, use_split);
		if (use_p > 0) {
			// Adjust use_p according to free ids
			float eff_use_p = (use_p * nIds) / nFreeIds;
			if (eff_use_p > 1.0)
				MTHROW_AND_ERR("Split_Learning_Set: Inconsistency at selection of subset for post-process learning : Not enough ids left for post-processor #%d\n", idx);

			// Assign
			int nAssigned = 0;
			for (int iId = 0; iId < nIds; iId++) {
				if (assignments[iId] == 0 && (globalRNG::rand() / (globalRNG::max() + 1.0)) < eff_use_p) {
					assignments[iId] = idx;
					nAssigned++;
				}
			}

			if (nAssigned == 0)
				MTHROW_AND_ERR("Split_Learning_Set:: Failed to assign any ids to post-processor #%d - use-p = %f , total ids = %d, before assignment, left with %d ids\n",
					idx, processor->use_p, nIds, nFreeIds);

			MLOG("Split_Learning_Set: Assigned %d ids out of %d to post-processor #%d with use-p = %f\n", nAssigned, nIds, idx, processor->use_p);
			nFreeIds -= nAssigned;
		}
		else if (use_split >= 0) {
			int nAssigned = 0;
			for (int iId = 0; iId < nIds; iId++) {
				if (assignments[iId] == 0 && inSamples.idSamples[iId].split == use_split) {
					assignments[iId] = idx;
					nAssigned++;
				}
			}

			if (nAssigned == 0)
				MTHROW_AND_ERR("Split_Learning_Set:: Failed to assign any ids to post-processor #%d - use-split = %d\n", idx, use_split);

			MLOG("Split_Learning_Set: Assigned %d ids out of %d to post-processor #%d with use-split = %d\n", nAssigned, nIds, idx, use_split);
			nFreeIds -= nAssigned;
		}

		if (nFreeIds == 0)
			MTHROW_AND_ERR("Split_Learning_Set: Left with no ids after processor #%d\n", idx);
		idx++;
	}

	MLOG("Split_Learning_Set: Assigned %d ids out of %d to model learning set\n", nFreeIds, nIds);

	// Create MedSamples
	model_learning_set.time_unit = inSamples.time_unit;
	for (auto& learning_set : post_processors_learning_sets)
		learning_set.time_unit = inSamples.time_unit;

	for (int i = 0; i < nIds; i++) {
		if (assignments[i] == 0)
			model_learning_set.idSamples.push_back(inSamples.idSamples[i]);
		else
			post_processors_learning_sets[assignments[i] - 1].idSamples.push_back(inSamples.idSamples[i]);
	}
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void MedModel::split_learning_set(MedSamples& inSamples, vector<MedSamples>& post_processors_learning_sets_orig, vector<MedSamples>& post_processors_learning_sets, MedSamples& model_learning_set)
{
	post_processors_learning_sets.resize(post_processors.size());

	// deciding which ids should be randomized/splut
	// we use all ids from inSamples, and add also all ids from the post_processors learning sets that have use_p > 0 or use_split>=0
	unordered_set<int> all_ids;
	for (auto &ids : inSamples.idSamples) all_ids.insert(ids.id);
	for (int j = 0; j < post_processors.size(); j++)
		if (post_processors[j]->use_p > 0 || post_processors[j]->use_split >= 0) {
			for (auto &s : post_processors_learning_sets_orig[j].idSamples)
				all_ids.insert(s.id);
		}

	vector<int> v_ids(all_ids.begin(), all_ids.end());
	int nIds = (int)v_ids.size();
	vector<int> assignments(nIds, 0);
	unordered_map<int, int> id2assignment;

	// Assign to post-processors
	int idx = 1;
	int nFreeIds = nIds;
	for (int iId = 0; iId < nIds; iId++)
		id2assignment[v_ids[iId]] = 0;
	for (PostProcessor *processor : post_processors) {
		float use_p = processor->get_use_p();
		int use_split = processor->get_use_split();
		if (use_p > 0 && use_split >= 0)
			MTHROW_AND_ERR("Split_Learning_Set: At most one of use_p (%f) & use_split (%d) allowed for post-processor\n", use_p, use_split);
		if (use_p > 0) {
			// Adjust use_p according to free ids
			float eff_use_p = (use_p * nIds) / nFreeIds;
			if (eff_use_p > 1.0)
				MTHROW_AND_ERR("Split_Learning_Set: Inconsistency at selection of subset for post-process learning : Not enough ids left for post-processor #%d\n", idx);

			// Assign
			int nAssigned = 0;
			for (int iId = 0; iId < nIds; iId++) {
				if (assignments[iId] == 0 && (globalRNG::rand() / (globalRNG::max() + 1.0)) < eff_use_p) {
					assignments[iId] = idx;
					id2assignment[v_ids[iId]] = idx;
					nAssigned++;
				}
			}

			if (nAssigned == 0)
				MTHROW_AND_ERR("Split_Learning_Set:: Failed to assign any ids to post-processor #%d - use-p = %f , total ids = %d, before assignment, left with %d ids\n",
					idx, processor->use_p, nIds, nFreeIds);

			MLOG("Split_Learning_Set: Assigned %d ids out of %d to post-processor #%d with use-p = %f\n", nAssigned, nIds, idx, processor->use_p);
			nFreeIds -= nAssigned;
		}
		else if (use_split >= 0) {
			int nAssigned = 0;
			for (int iId = 0; iId < nIds; iId++) {
				if (assignments[iId] == 0 && inSamples.idSamples[iId].split == use_split) {
					assignments[iId] = idx;
					id2assignment[v_ids[iId]] = idx;
					nAssigned++;
				}
			}

			if (nAssigned == 0)
				MTHROW_AND_ERR("Split_Learning_Set:: Failed to assign any ids to post-processor #%d - use-split = %d\n", idx, use_split);

			MLOG("Split_Learning_Set: Assigned %d ids out of %d to post-processor #%d with use-split = %d\n", nAssigned, nIds, idx, use_split);
			nFreeIds -= nAssigned;
		}

		if (nFreeIds == 0)
			MTHROW_AND_ERR("Split_Learning_Set: Left with no ids after processor #%d\n", idx);
		idx++;
	}

	MLOG("Split_Learning_Set: Assigned %d ids out of %d to model learning set\n", nFreeIds, nIds);

	// Create MedSamples
	model_learning_set.time_unit = inSamples.time_unit;
	for (auto& learning_set : post_processors_learning_sets)
		learning_set.time_unit = inSamples.time_unit;

	for (auto &s : inSamples.idSamples)
		if (id2assignment[s.id] == 0)
			model_learning_set.idSamples.push_back(s);


	if (post_processors.size() > 0) {
		if (post_processors_learning_sets_orig.size() == post_processors.size()) {
			MLOG("Split_Learning_Set: Building lists from given post processors lists\n");
			// case user gave lists to work with
			// In this case we have several options:
			// (1) No use_p and use_split are given : in this case we take the orig list as is (!!)
			// (2) use_p given : we make sure not to use ids from orig that were not selected for this case in idSamples (but leave the others).
			// (3) use_split given : choose by split
			for (int j = 0; j < post_processors.size(); j++) {


				if (post_processors[j]->use_p > 0 || post_processors[j]->use_split >= 0) {

					for (int i = 0; i < post_processors_learning_sets_orig[j].idSamples.size(); i++) {
						int id = post_processors_learning_sets_orig[j].idSamples[i].id;
						if (id2assignment.find(id) == id2assignment.end() || id2assignment[id] == j + 1)
							post_processors_learning_sets[j].idSamples.push_back(post_processors_learning_sets_orig[j].idSamples[i]);
					}

				}
				else {
					post_processors_learning_sets[j] = post_processors_learning_sets_orig[j];
				}

				MLOG("Split_Learning_Set: post processor %d : orig %d ids : selected %d ids ( use_p is %f , use_split is %d )\n",
					j, post_processors_learning_sets[j].idSamples.size(), post_processors_learning_sets_orig[j].idSamples.size(), post_processors[j]->use_p, post_processors[j]->use_split);
			}
		}
		else
		{
			// case we need to build lists from inSamples

			MLOG("Split_Learning_Set: Building lists from selected learning set members\n");
			for (int i = 0; i < nIds; i++) {
				if (assignments[i] != 0)
					post_processors_learning_sets[assignments[i] - 1].idSamples.push_back(inSamples.idSamples[i]);
			}

		}
	}



}


// Adjust model according to signals available in repository
//--------------------------------------------------------------------------------------------------------
void MedModel::fit_for_repository(MedPidRepository& rep) {

	// Currently - only RepProcessors are adjustable
	for (RepProcessor *processor : rep_processors)
		processor->fit_for_repository(rep);
}

// loading a repository (optionally allowing for adjustment to model according to available signals)
//--------------------------------------------------------------------------------------------------------
void MedModel::load_repository(const string& configFile, MedPidRepository& rep, bool allow_adjustment) {

	vector<int> empty_ids_list;
	load_repository(configFile, empty_ids_list, rep, allow_adjustment);
}

void MedModel::load_repository(const string& configFile, vector<int> ids, MedPidRepository& rep, bool allow_adjustment) {

	// Adjust Model
	if (allow_adjustment) {
		if (rep.init(configFile) < 0)
			MTHROW_AND_ERR("Cannot initialize repository from %s\n", configFile.c_str());
		fit_for_repository(rep);
	}

	// Get Required signals
	vector<string> req_signals;
	get_required_signal_names(req_signals);

	// Read Repository
	MLOG("Reading Repository from %s\n", configFile.c_str());
	if (rep.read_all(configFile, ids, req_signals) != 0)
		MTHROW_AND_ERR("Read repository from %s failed\n", configFile.c_str());
}


//========================================================================================================
// medial::medmodel:: functions
//========================================================================================================

//--------------------------------------------------------------------------------------------------------
void medial::medmodel::apply(MedModel &model, MedSamples &samples, string rep_fname, MedModelStage to_stage)
{
	unordered_set<string> req_sigs;
	vector<string> rsigs;

	MedPidRepository rep;
	if (rep.init(rep_fname) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", rep_fname.c_str());
	model.fit_for_repository(rep);
	model.get_required_signal_names(req_sigs);
	for (auto &s : req_sigs) rsigs.push_back(s);

	vector<int> pids;

	samples.get_ids(pids);

	if (rep.read_all(rep_fname, pids, rsigs) < 0)
		MTHROW_AND_ERR("medial::medmodel::apply() ERROR :: could not read repository %s\n", rep_fname.c_str());

	if (model.apply(rep, samples, MED_MDL_APPLY_FTR_GENERATORS, to_stage) < 0)
		MTHROW_AND_ERR("medial::medmodel::apply() ERROR :: could not apply model\n");
}
//--------------------------------------------------------------------------------------------------------
void medial::medmodel::apply(MedModel &model, string rep_fname, string f_samples, MedSamples &samples, MedModelStage to_stage)
{
	if (samples.read_from_file(f_samples) < 0)
		MTHROW_AND_ERR("medial::medmodel::apply() ERROR :: could not read samples file %s\n", f_samples.c_str());

	medial::medmodel::apply(model, samples, rep_fname, to_stage);
}
//--------------------------------------------------------------------------------------------------------
void medial::medmodel::apply(MedModel &model, string rep_fname, string f_samples, MedModelStage to_stage)
{
	// returns just the model : model.features is updated
	MedSamples samples;
	medial::medmodel::apply(model, rep_fname, f_samples, samples, to_stage);
}

void MedModel::change_model(const vector<ChangeModelInfo> &change_request) {
	for (const ChangeModelInfo &req : change_request)
		change_model(req);
}

void MedModel::clean_model() {
	vector<RepProcessor *> final_res;
	for (size_t i = 0; i < rep_processors.size(); ++i)
	{
		if (rep_processors[i] != NULL) {
			if (rep_processors[i]->processor_type == REP_PROCESS_MULTI) {
				RepMultiProcessor *multi = static_cast<RepMultiProcessor *>(rep_processors[i]);
				vector<RepProcessor *> final_res_internal;
				for (size_t j = 0; j < multi->processors.size(); ++j)
					if (multi->processors[j] != NULL)
						final_res_internal.push_back(multi->processors[j]);
				multi->processors = move(final_res_internal);
				if (multi->processors.empty()) {
					delete multi;
					rep_processors[i] = NULL;
				}
			}
			if (rep_processors[i] != NULL)
				final_res.push_back(rep_processors[i]);
		}
	}
	rep_processors = move(final_res);

	vector<FeatureGenerator *> final_res_g;
	for (size_t i = 0; i < generators.size(); ++i)
		if (generators[i] != NULL)
			final_res_g.push_back(generators[i]);
	generators = move(final_res_g);

	vector<FeatureProcessor *> final_res_fp;
	for (size_t i = 0; i < feature_processors.size(); ++i)
	{
		if (feature_processors[i] != NULL) {
			if (feature_processors[i]->processor_type == FTR_PROCESS_MULTI) {
				MultiFeatureProcessor *multi = static_cast<MultiFeatureProcessor *>(feature_processors[i]);
				vector<FeatureProcessor *> final_res_internal;
				for (size_t j = 0; j < multi->processors.size(); ++j)
					if (multi->processors[j] != NULL)
						final_res_internal.push_back(multi->processors[j]);
				multi->processors = move(final_res_internal);
				if (multi->processors.empty()) {
					delete multi;
					feature_processors[i] = NULL;
				}
			}
			if (feature_processors[i] != NULL)
				final_res_fp.push_back(feature_processors[i]);
		}
	}
	feature_processors = move(final_res_fp);

	vector<PostProcessor *> final_res_pp;
	for (size_t i = 0; i < post_processors.size(); ++i)
	{
		if (post_processors[i] != NULL) {
			if (post_processors[i]->processor_type == FTR_POSTPROCESS_MULTI) {
				MultiPostProcessor *multi = static_cast<MultiPostProcessor *>(post_processors[i]);
				vector<PostProcessor *> final_res_internal;
				for (size_t j = 0; j < multi->post_processors.size(); ++j)
					if (multi->post_processors[j] != NULL)
						final_res_internal.push_back(multi->post_processors[j]);
				multi->post_processors = move(final_res_internal);
				if (multi->post_processors.empty()) {
					delete multi;
					post_processors[i] = NULL;
				}
			}
			if (post_processors[i] != NULL)
				final_res_pp.push_back(post_processors[i]);
		}
	}
	post_processors = move(final_res_pp);
}


void MedModel::find_object(RepProcessor *c, vector<RepProcessor *> &res, vector<RepProcessor **> &res_pointer) {
	for (size_t i = 0; i < rep_processors.size(); ++i)
	{
		if (typeid(*rep_processors[i]) == typeid(*c)) {
			res.push_back(rep_processors[i]);
			res_pointer.push_back(&rep_processors[i]);
		}
		else if (rep_processors[i]->processor_type == REP_PROCESS_MULTI) {
			RepMultiProcessor *multi = static_cast<RepMultiProcessor *>(rep_processors[i]);
			for (size_t j = 0; j < multi->processors.size(); ++j) {
				if (typeid(*multi->processors[j]) == typeid(*c)) {
					res.push_back(multi->processors[j]);
					res_pointer.push_back(&multi->processors[j]);
				}
			}
		}
	}
}

void MedModel::find_object(FeatureGenerator *c, vector<FeatureGenerator *> &res, vector<FeatureGenerator **> &res_pointer) {
	for (size_t i = 0; i < generators.size(); ++i)
		if (typeid(*generators[i]) == typeid(*c)) {
			res.push_back(generators[i]);
			res_pointer.push_back(&generators[i]);
		}
}

void MedModel::find_object(FeatureProcessor *c, vector<FeatureProcessor *> &res, vector<FeatureProcessor **> &res_pointer) {
	for (size_t i = 0; i < feature_processors.size(); ++i)
	{
		if (typeid(*feature_processors[i]) == typeid(*c)) {
			res.push_back(feature_processors[i]);
			res_pointer.push_back(&feature_processors[i]);
		}
		else if (feature_processors[i]->processor_type == FTR_PROCESS_MULTI) {
			MultiFeatureProcessor *multi = static_cast<MultiFeatureProcessor *>(feature_processors[i]);
			for (size_t j = 0; j < multi->processors.size(); ++j) {
				if (typeid(*multi->processors[j]) == typeid(*c)) {
					res.push_back(multi->processors[j]);
					res_pointer.push_back(&multi->processors[j]);
				}
			}
		}
	}
}

void MedModel::find_object(PostProcessor *c, vector<PostProcessor *> &res, vector<PostProcessor **> &res_pointer) {
	for (size_t i = 0; i < post_processors.size(); ++i)
	{
		if (typeid(*post_processors[i]) == typeid(*c)) {
			res.push_back(post_processors[i]);
			res_pointer.push_back(&post_processors[i]);
		}
		else if (post_processors[i]->processor_type == FTR_POSTPROCESS_MULTI) {
			MultiPostProcessor *multi = static_cast<MultiPostProcessor *>(post_processors[i]);
			for (size_t j = 0; j < multi->post_processors.size(); ++j) {
				if (typeid(*multi->post_processors[j]) == typeid(*c)) {
					res.push_back(multi->post_processors[j]);
					res_pointer.push_back(&multi->post_processors[j]);
				}
			}
		}
	}
}

void MedModel::find_object(MedPredictor *c, vector<MedPredictor *> &res, vector<MedPredictor **> &res_pointer) {
	if (typeid(*predictor) == typeid(*c)) {
		res.push_back(predictor);
		res_pointer.push_back(&predictor);
	}
}

bool filter_by_json_query(SerializableObject &obj, const vector<string> &json_query_whitelist,
	const vector<string> &json_query_blacklist) {
	if (json_query_whitelist.empty() && json_query_blacklist.empty())
		return true;

	string obj_json = obj.object_json();
	bool status = true;
	//check whitelist:
	for (const string &s : json_query_whitelist)
	{
		boost::regex reg_pat(s);
		status = boost::regex_search(obj_json, reg_pat);
		if (!status)
			break;
	}
	//check blacklist if passed whitelist:
	if (status) {
		for (const string &s : json_query_blacklist)
		{
			boost::regex reg_pat(s);
			status = !boost::regex_search(obj_json, reg_pat);
			if (!status)
				break;
		}
	}

	return status;
}

template <class T> void MedModel::apply_change(const ChangeModelInfo &change_request, void *obj) {
	T *c = static_cast<T *>(obj);
	T *cl_name = new T;
	string class_name = cl_name->my_class_name();
	delete cl_name;
	if (change_request.verbose_level >= 2)
		MLOG("Change Request \"%s\" :: Identified object \"%s\" as %s\n",
			change_request.change_name.c_str(), change_request.object_type_name.c_str(), class_name.c_str());
	vector<T *> found_res;
	vector<T **> found_p_res;
	find_object(c, found_res, found_p_res);

	int did_cnt = 0;
	int succ_cnt = 0;
	for (size_t i = 0; i < found_res.size(); ++i)
	{
		if (!filter_by_json_query(*found_res[i], change_request.json_query_whitelist, change_request.json_query_blacklist))
			continue;
		++did_cnt;
		if (change_request.change_command == "DELETE") {
			delete found_res[i];
			//need to erase from original pos:
			*found_p_res[i] = NULL;
			++succ_cnt;
		}
		else {
			if (change_request.change_command == "PRINT") {
				string str = found_res[i]->object_json();
				cout << "PRINT:\n" << str << endl;
				++succ_cnt;
			}
			else {
				try {
					found_res[i]->init_from_string(change_request.change_command);
					++succ_cnt;
				}
				catch (exception &exp) {
					MWARN("Error in sending init, got %s\n", exp.what());
				}
			}
		}
	}
	if (did_cnt > 0 && change_request.change_command == "DELETE")
		clean_model();
	if (change_request.verbose_level >= 2)
		MLOG("Change Request \"%s\" :: Found object as %s - matched to %d objects - succeed in %d\n",
			change_request.change_name.c_str(), class_name.c_str(), did_cnt, succ_cnt);
	else if (change_request.verbose_level >= 1 && (did_cnt == 0 || succ_cnt != did_cnt))
		MLOG("Change Request \"%s\" :: Found object as %s - matched to %d objects - succeed in %d\n",
			change_request.change_name.c_str(), class_name.c_str(), did_cnt, succ_cnt);
	delete c;
}

void MedModel::change_model(const ChangeModelInfo &change_request) {
	//non interactive:
	RepProcessor test_rep;
	FeatureGenerator test_gen;
	FeatureProcessor ftr_processor; //needs try catch
	PostProcessor pp_processor;
	MedPredictor predictor;

	if (change_request.object_type_name == "MedModel") {
		//Global model changes:
		map<string, string> mapper;
		if (MedSerialize::init_map_from_string(change_request.change_command, mapper) < 0)
			MTHROW_AND_ERR("Error Init from String %s\n", change_request.change_command.c_str());
		for (const auto &it : mapper)
		{
			if (it.first == "max_data_in_mem")
				max_data_in_mem = med_stoi(it.second);
			else if (it.first == "serialize_learning_set")
				serialize_learning_set = med_stoi(it.second) > 0;
			else if (it.first == "take_mean_pred")
				take_mean_pred = med_stoi(it.second) > 0;
			else if (it.first == "generate_masks_for_features")
				generate_masks_for_features = med_stoi(it.second) > 0;
			else if (it.first == "verbosity" || it.first == "verbose")
				verbosity = med_stoi(it.second) > 0;
			else
				MTHROW_AND_ERR("Error - unkown model param \"%s\"", it.first.c_str());
		}
		return;
	}

	global_logger.levels[LOG_REPCLEANER] = MAX_LOG_LEVEL + 1;
	void *obj = test_rep.new_polymorphic(change_request.object_type_name);
	global_logger.levels[LOG_REPCLEANER] = LOG_DEF_LEVEL;
	if (obj != NULL) {
		apply_change<RepProcessor>(change_request, obj);
		return;
	}

	global_logger.levels[LOG_FTRGNRTR] = MAX_LOG_LEVEL + 1;
	void *obj_gen = test_gen.new_polymorphic(change_request.object_type_name);
	global_logger.levels[LOG_FTRGNRTR] = LOG_DEF_LEVEL;
	if (obj_gen != NULL) {
		apply_change<FeatureGenerator>(change_request, obj_gen);
		return;
	}

	void *obj_ftr_processor = NULL;
	global_logger.levels[LOG_FEATCLEANER] = MAX_LOG_LEVEL + 1;
	try {
		obj_ftr_processor = ftr_processor.new_polymorphic(change_request.object_type_name);
	}
	catch (...) {
		obj_ftr_processor = NULL;
	}
	global_logger.levels[LOG_FEATCLEANER] = LOG_DEF_LEVEL;
	if (obj_ftr_processor != NULL) {
		apply_change<FeatureProcessor>(change_request, obj_ftr_processor);
		return;
	}

	global_logger.levels[LOG_MED_MODEL] = MAX_LOG_LEVEL + 1;
	void *obj_post_processor = pp_processor.new_polymorphic(change_request.object_type_name);
	global_logger.levels[LOG_MED_MODEL] = LOG_DEF_LEVEL;
	if (obj_post_processor != NULL) {
		apply_change<PostProcessor>(change_request, obj_post_processor);
		return;
	}

	global_logger.levels[LOG_MEDALGO] = MAX_LOG_LEVEL + 1;
	void *obj_predictor = predictor.new_polymorphic(change_request.object_type_name);
	global_logger.levels[LOG_MEDALGO] = LOG_DEF_LEVEL;
	if (obj_predictor != NULL) {
		apply_change<MedPredictor>(change_request, obj_predictor);
		return;
	}

	if (change_request.verbose_level >= 1)
		MLOG("Warn in change request \"%s\" - can't find object %s - no change\n",
			change_request.change_name.c_str(), change_request.object_type_name.c_str());
}

int ChangeModelInfo::init(map<string, string>& mapper) {
	vector<string> field_names;
	serialized_fields_name(field_names);
	string option_fields = medial::io::get_list(field_names);

	for (auto &it : mapper)
	{
		//! [ChangeModelInfo::init]
		if (it.first == "change_command")
			change_command = it.second;
		else if (it.first == "change_name")
			change_name = it.second;
		else if (it.first == "object_type_name")
			object_type_name = it.second;
		else if (it.first == "verbose_level")
			verbose_level = med_stoi(it.second);
		else if (it.first == "json_query_whitelist")
			boost::split(json_query_whitelist, it.second, boost::is_any_of("~"));
		else if (it.first == "json_query_blacklist")
			boost::split(json_query_blacklist, it.second, boost::is_any_of("~"));
		else if (it.first == "from_json_file") {
			vector<ChangeModelInfo> res;
			ifstream inf(it.second);
			if (!inf)
				MTHROW_AND_ERR("can't open json file [%s] for read\n", it.second.c_str());
			stringstream sstr;
			sstr << inf.rdbuf();
			inf.close();
			string json_content = stripComments(sstr.str());
			ChangeModelInfo::parse_json_string(json_content, res);
			if (res.size() != 1)
				MTHROW_AND_ERR("Can read only 1 ChangeModelInfo from ChangeModelInfo::init\n");
			*this = move(res[0]);
		}
		//! [ChangeModelInfo::init]
		else
			MTHROW_AND_ERR("Error ChangeModelInfo::init - unsupported param \"%s\", options: [%s]\n",
				it.first.c_str(), option_fields.c_str());
	}
	return 0;
}

void MedModel::read_from_file_with_changes(const string &model_binary_path, const string &path_to_json_changes) {
	read_from_file(model_binary_path);

	if (!path_to_json_changes.empty()) {
		vector<string> alts;
		string json_content = json_file_to_string(0, path_to_json_changes, alts);
		vector<ChangeModelInfo> change_reqs;
		ChangeModelInfo::parse_json_string(json_content, change_reqs);
		change_model(change_reqs);
	}
}

void MedModel::clone_model(MedModel &out) {
	out.clear();
	vector<unsigned char> blob;
	serialize_vec(blob);
	out.deserialize_vec(blob);
}

void MedModel::copy_from_model(MedModel &in) {
	clear();
	vector<unsigned char> blob;
	in.serialize_vec(blob);
	deserialize_vec(blob);
}

#endif
