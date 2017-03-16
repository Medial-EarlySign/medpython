#include "MedModel.h"
#include "MedProcessUtils.h"
#include <omp.h>

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// MedModel
//=======================================================================================
// Learn
//.......................................................................................
int MedModel::learn(MedPidRepository& rep, MedSamples* samples, MedModelStage start_stage, MedModelStage end_stage) {

	DummyFeatsSelector selector;
	return learn(rep, samples, selector, start_stage, end_stage);

}

//.......................................................................................
int MedModel::learn(MedPidRepository& rep, MedSamples* _samples, FeatureSelector& selector, MedModelStage start_stage, MedModelStage end_stage) {

	MedTimer timer;

	// Stage Sanity
	if (start_stage > end_stage) {
		MERR("MedModel learn() : Illegal start and end\n");
		return -1;
	}

	// Set of signals
	init_signal_ids(rep.dict);

	LearningSet = _samples;

	// Set of ids
	vector<int> ids;
	LearningSet->get_ids(ids);

	// Learn RepCleaners
	if (start_stage <= MED_MDL_REP_PROCESSORS) {
		timer.start();
		if (learn_rep_processors(rep, ids) < 0) { //??? why are rep processors initialized for ALL time points in an id??
			MERR("MedModel learn() : ERROR: Failed learn_rep_processors()\n");
			return -1;
		}
		timer.take_curr_time();
		MLOG("MedModel::learn() : learn rep processors time %g ms\n", timer.diff_milisec());
	}
	if (end_stage <= MED_MDL_REP_PROCESSORS)
		return 0;

	// Learn Feature Generators
	if (start_stage <= MED_MDL_FTR_GENERATORS) {
		timer.start();
		if (learn_feature_generators(rep, LearningSet) < 0) {
			MERR("MedModel learn() : ERROR: Failed learn_feature_generators\n");
			return -1;
		}
		timer.take_curr_time();
		MLOG("MedModel::learn() : learn feature generators %g ms\n", timer.diff_milisec());
	}
	if (end_stage <= MED_MDL_FTR_GENERATORS)
		return 0;

	features.clear();
	features.set_time_unit(LearningSet->time_unit);

	// Generate features
	timer.start();
	if (generate_all_features(rep, LearningSet, features) < 0) {
		MERR("MedModel learn() : ERROR: Failed generate_all_features()\n");
		return -1;
	}
	timer.take_curr_time();
	MLOG("MedModel::learn() : generating learn matrix time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());

	// Learn Feature processors and apply
	if (start_stage <= MED_MDL_FTR_PROCESSORS) {
		timer.start();
		if (learn_and_apply_feature_processors(features) < 0) {
			MERR("MedModel::learn() : ERROR: Failed learn_and_apply_feature_cleaners()\n");
			return -1;
		}
		timer.take_curr_time();
		MLOG("MedModel::learn() : feature processing learn and apply time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());
	}
	else {
		// Just apply feature processors
		timer.start();
		if (apply_feature_processors(features) < 0) {
			MERR("MedModel::apply() : ERROR: Failed apply_feature_cleaners()\n");
			return -1;
		}
		timer.take_curr_time();
		MLOG("MedModel::learn() : feature processing time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());
	}
	if (end_stage <= MED_MDL_FTR_PROCESSORS)
		return 0;

	// Select Features
	selector.select(generators);

	// Learn predictor
	if (start_stage <= MED_MDL_PREDICTOR) {
		timer.start();
		int rc = predictor->learn(features);
		timer.take_curr_time();
		MLOG("MedModel::learn() : model train time: %g ms\n", timer.diff_milisec());
		if (rc != 0)
			return rc;
	}
	if (end_stage <= MED_MDL_PREDICTOR)
		return 0;


	return 0;
}

//.......................................................................................
// Apply
int MedModel::apply(MedPidRepository& rep, MedSamples& samples, MedModelStage end_stage) {

	// Stage Sanity
	if (end_stage == MED_MDL_REP_PROCESSORS) {
		MERR("MedModel apply() : Illegal end stage %d\n",end_stage);
		return -1;
	}

	// Set of signals
	init_signal_ids(rep.dict);

	vector<int> req_signals;
	for (int signalId : required_signals)
		req_signals.push_back(signalId);

	// Generate features
	//MedFeatures features(samples.time_unit);
	features.clear();
	features.set_time_unit(samples.time_unit);

	if (generate_all_features(rep, &samples, features) < 0) {
		MERR("MedModel apply() : ERROR: Failed generate_all_features()\n");
		return -1;
	}
	if (end_stage <= MED_MDL_FTR_GENERATORS)
		return 0;

	// Process Features
	if (apply_feature_processors(features) < 0) {
		MERR("MedModel::apply() : ERROR: Failed apply_feature_cleaners()\n");
		return -1;
	}
	if (end_stage <= MED_MDL_FTR_PROCESSORS)
		return 0;

	// Apply predictor
	if (predictor->predict(features) < 0) {
		MERR("Predictor failed\n");
		return -1;
	}

	samples.insert_preds(features);
	if (end_stage <= MED_MDL_PREDICTOR)
		return 0;
	return 0;
}

//.......................................................................................
// Learn rep-cleaning
int MedModel::quick_learn_rep_processors(MedPidRepository& rep, vector<int>& ids) {

	vector<int> rc(rep_processors.size(), 0);

#pragma omp parallel for schedule(dynamic)
	for (int i=0; i<rep_processors.size(); i++)
		rc[i] = rep_processors[i]->learn(rep,ids);

	for (auto RC : rc) if (RC < 0)	return -1;
	return 0;
}

//.......................................................................................
int MedModel::learn_feature_generators(MedPidRepository &rep, MedSamples *learn_samples)
{
	vector<int> ids;
	learn_samples->get_ids(ids);

	vector<int> rc(generators.size(), 0);
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<generators.size(); i++)
		rc[i] = generators[i]->learn(rep, ids,rep_processors); //??? why is this done for ALL time points in an id???

	for (auto RC : rc) if (RC < 0)	return -1;
	return 0;
}

//.......................................................................................
int MedModel::generate_features(MedPidRepository &rep, MedSamples *samples, vector<FeatureGenerator *>& _generators, MedFeatures &features)
{
	vector<int> req_signals;
	for (int signalId : required_signals)
		req_signals.push_back(signalId);

	// init features attributes
	for (auto& generator : _generators)
		generator->init(features);

	// preparing records and features for threading
	int N_tot_threads = omp_get_max_threads();
	MLOG("MedModel::learn/apply() : feature generation with %d threads\n", N_tot_threads);
	vector<PidDynamicRec> idRec(N_tot_threads);
	features.init_all_samples(samples->idSamples);
	int samples_size = (int)features.samples.size();
	for (auto &generator : _generators) {
		for (string& name : generator->names)
			features.data[name].resize(samples_size, 0);
	}

	// Loop on ids
	int RC = 0;
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j<samples->idSamples.size(); j++) {
		MedIdSamples& pid_samples = samples->idSamples[j];
		int n_th = omp_get_thread_num();
		int rc = 0;

		// Generate DynamicRec with all relevant signals
		if (idRec[n_th].init_from_rep(std::addressof(rep), pid_samples.id, req_signals, (int)pid_samples.samples.size()) < 0) rc = -1;

		// Apply rep-cleaning
		for (auto& processor : rep_processors)
			if (processor->apply(idRec[n_th], pid_samples) < 0) rc = -1;

		// Generate Features
		for (auto& generator : _generators)
			if (generator->generate(idRec[n_th], features) < 0) rc = -1;
#pragma omp critical 
		if (rc < 0) RC = -1;
	}

	return RC;
}

//.......................................................................................
int MedModel::learn_and_apply_feature_processors(MedFeatures &features)
{
	for (auto& processor : feature_processors) {
		if (processor->learn(features) < 0) return -1;
		if (processor->apply(features) < 0) return -1;
	}
	return 0;
}

//.......................................................................................
int MedModel::apply_feature_processors(MedFeatures &features)
{
	for (auto& processor : feature_processors) {
		if (processor->apply(features) < 0) return -1;
	}
	return 0;
}

//.......................................................................................
// Learn rep-cleaning iteratively, must be serial...
int MedModel::learn_rep_processors(MedPidRepository& rep, vector<int>& ids) {

	vector<RepProcessor *> temp_processors;
	for (RepProcessor *processor : rep_processors) {
		if (processor->learn(rep, ids, temp_processors) < 0) return -1;
		temp_processors.push_back(processor);
	}

	return 0;
}

// Required Signals
//.......................................................................................
void MedModel::get_required_signals(MedDictionarySections& dict) {

	required_signals.clear();

	for (RepProcessor *processor : rep_processors)
		processor->get_required_signal_ids(required_signals,dict);

	int ii = 0;
	for (FeatureGenerator *generator : generators) 
		generator->get_required_signal_ids(required_signals, dict);

}

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
		rep_processors.resize(i_set+1, NULL);

		// put a new empty multi in i_set
		RepMultiProcessor *processor = new RepMultiProcessor;
		rep_processors[i_set] = processor;
	}

	// Now we are at a state in which we have a multi processor at i_set and need to create a new processor and push it in
	string in = init_string;
	RepProcessor *rep_proc = RepProcessor::create_processor(in);

	// push it in
	((RepMultiProcessor *)rep_processors[i_set])->processors.push_back(rep_proc);

}

//.......................................................................................
// fp_type and feature name are must have parameters
void MedModel::add_feature_processor_to_set(int i_set, const string &init_string)
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
			MultiFeatureProcessor *mfprocessor = new MultiFeatureProcessor;
			feature_processors[i_set] = mfprocessor;
		}
		else if (feature_processors[i_set]->processor_type != FTR_PROCESS_MULTI) {
			// the processor was not multi, and hence we create one switch it , and push the current into it
			FeatureProcessor *curr_fp = feature_processors[i_set];
			MultiFeatureProcessor *mfprocessor = new MultiFeatureProcessor;
			feature_processors[i_set] = mfprocessor;
			mfprocessor->processors.push_back(curr_fp);

		}
	}
	else {
		// resize feature_processors
		feature_processors.resize(i_set+1, NULL);

		// put a new empty multi in i_set
		MultiFeatureProcessor *mfprocessor = new MultiFeatureProcessor;
		feature_processors[i_set] = mfprocessor;
	}

	// Now we are at a state in which we have a multi feature processor at i_set and need to create a new processor or processors and push it in

	// get all relevant features names
	string feat_names;
	get_single_val_from_init_string(init_string, "names", feat_names);

	vector<string> features;
	if (feat_names == "" || feat_names == "All")
		get_all_features_names(features);
	else
		boost::split(features, feat_names, boost::is_any_of(","));

	// get type of feature processor
	string fp_type;
	get_single_val_from_init_string(init_string, "fp_type", fp_type);
	FeatureProcessorTypes type = feature_processor_name_to_type(fp_type);
	//MLOG("fp_type = %s %d\n", fp_type.c_str(), (int)type);

	// actual adding to relevant MultiFeatureProcessor
	((MultiFeatureProcessor *)feature_processors[i_set])->add_processors_set(type, features, init_string);
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

//.......................................................................................
void MedModel::add_process_to_set(int i_set, const string &init_string)
{
	if (init_string.find("rp_type") != string::npos) return add_rep_processor_to_set(i_set, init_string);
	if (init_string.find("fg_type") != string::npos) return add_feature_generator_to_set(i_set, init_string);
	if (init_string.find("fp_type") != string::npos) return add_feature_processor_to_set(i_set, init_string);

	MERR("add_process_to_set():: Can't process line %s\n", init_string.c_str());
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
void MedModel::get_affected_signals(MedDictionarySections& dict) {

	for (RepProcessor *processor : rep_processors)
		processor->get_affected_signal_ids(dict);

}

// All signal ids
//.......................................................................................
void MedModel::init_signal_ids(MedDictionarySections& dict) {

	get_affected_signals(dict);
	get_required_signals(dict);

	for (RepProcessor *processor : rep_processors)
		processor->get_signal_ids(dict);

	for (FeatureGenerator *generator : generators)
		generator->get_signal_ids(dict);
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
	get_all_features_names(features);

//	for (auto &name : features) 
//		MLOG("Adding %s to processors of type %d\n", name.c_str(),type);

	add_feature_processors_set(type, features);

}

//.......................................................................................
void MedModel::add_feature_processors_set(FeatureProcessorTypes type, string init_string) {

	vector<string> features;
	get_all_features_names(features);

	add_feature_processors_set(type, features,init_string);
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
		FeatureGenerator *generator = FeatureGenerator::make_generator(type, init_string + ";signalName="+signal);
		add_feature_generator(generator);
	}
}

//.......................................................................................
void MedModel::get_all_features_names(vector<string> &feat_names)
{
	feat_names.clear();
	for (unsigned int i = 0; i < generators.size(); i++) {
		for (string& name : generators[i]->names)
			feat_names.push_back(name);
	}

}

//.......................................................................................
// De(Serialize)
size_t MedModel::get_size() {

	size_t size = 0;

	// Rep-Cleaners
	size += sizeof(size_t); 
	for (auto& processor : rep_processors)
		size += processor->get_processor_size();

	// Feature Generators 
	size += sizeof(size_t);
	for (auto& generator : generators)
		size += generator->get_generator_size();

	// Features-level cleaners;
	size += sizeof(size_t);
	for (auto& processor : feature_processors) 
		size += processor->get_processor_size();

	// Predictor
	size += predictor->get_predictor_size();

	// Learning samples
	size += LearningSet->get_size();
	  
	return size;
}

//.......................................................................................
size_t MedModel::serialize(unsigned char *blob) {

	size_t ptr = 0;

	// Rep-Cleaners
	size_t n = rep_processors.size();
	memcpy(blob + ptr, &n, sizeof(size_t)); ptr += sizeof(size_t);

	for (auto& processor : rep_processors)
		ptr += processor->processor_serialize(blob + ptr);

	// Feature Generators
	n = generators.size();
	memcpy(blob + ptr, &n, sizeof(size_t)); ptr += sizeof(size_t);

	for (auto& generator : generators)
		ptr += generator->generator_serialize(blob + ptr);

	// Features-level processors;
	n = feature_processors.size();
	memcpy(blob + ptr, &n, sizeof(size_t)); ptr += sizeof(size_t);

	for (auto& processor : feature_processors) 
		ptr += processor->processor_serialize(blob + ptr);

	// Predictor
	ptr += predictor->predictor_serialize(blob+ptr);

	// Learning samples
	ptr += LearningSet->serialize(blob + ptr);

	return ptr;
}

//.......................................................................................
size_t MedModel::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	// Rep-Processors
	size_t n; 
	memcpy(&n, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);

	rep_processors.resize(n);
	for (size_t i = 0; i < n; i++) {
		RepProcessorTypes type;
		memcpy(&type, blob + ptr, sizeof(RepProcessorTypes)); ptr += sizeof(RepProcessorTypes);
		rep_processors[i] = RepProcessor::make_processor(type);
		ptr += rep_processors[i]->deserialize(blob + ptr);

	}

	// Feature Generators
	memcpy(&n, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);

	generators.resize(n);
	for (size_t i = 0; i < n; i++) {
		FeatureGeneratorTypes type;
		memcpy(&type, blob + ptr, sizeof(FeatureGeneratorTypes)); ptr += sizeof(FeatureGeneratorTypes);
		generators[i] = FeatureGenerator::make_generator(type);
		ptr += generators[i]->deserialize(blob + ptr);

	}

	// Features-level processors;
	memcpy(&n, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t); 

	feature_processors.resize(n);
	for (size_t i = 0; i < n; i++) {
		FeatureProcessorTypes type;
		memcpy(&type, blob + ptr, sizeof(FeatureProcessorTypes)); ptr += sizeof(FeatureProcessorTypes);
		feature_processors[i] = FeatureProcessor::make_processor(type);
		ptr += feature_processors[i]->deserialize(blob + ptr);
	}

	// Predictor
	MedPredictorTypes type;
	memcpy(&type, blob + ptr, sizeof(MedPredictorTypes)); ptr += sizeof(MedPredictorTypes);
	predictor = MedPredictor::make_predictor(type);
	ptr += predictor->deserialize(blob + ptr);

	// Learning samples
	LearningSet = new MedSamples;
	ptr += LearningSet->deserialize(blob + ptr);

	return ptr;
}