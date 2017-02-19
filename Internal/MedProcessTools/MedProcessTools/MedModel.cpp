#include "MedModel.h"
#include <omp.h>

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// MedModel
//=======================================================================================
// Learn
int MedModel::learn(MedPidRepository& rep, MedSamples* samples) {
	
	DummyFeatsSelector selector;
	return learn(rep, samples, selector);

}

//.......................................................................................
int MedModel::learn(MedPidRepository& rep, MedSamples* _samples, FeatureSelector& selector) {

	MedTimer timer;

	LearningSet = _samples;

	// Set of ids
	vector<int> ids;
	LearningSet->get_ids(ids);

	// Learn RepCleaners
	timer.start();
	if (learn_rep_processors(rep, ids) < 0) {
		MERR("MedModel learn() : ERROR: Failed learn_rep_processors()\n");
		return -1;
	}
	timer.take_curr_time();
	MLOG("MedModel::learn() : learn rep processors time %g ms\n", timer.diff_milisec());

	// Learn Feature Generators
	timer.start();
	if (learn_feature_generators(rep, LearningSet) < 0) {
		MERR("MedModel learn() : ERROR: Failed learn_feature_generators\n");
		return -1;
	}
	timer.take_curr_time();
	MLOG("MedModel::learn() : learn feature generators %g ms\n", timer.diff_milisec());

	MedFeatures features;

	// Generate features
	timer.start();
	if (generate_all_features(rep, LearningSet, features) < 0) {
		MERR("MedModel learn() : ERROR: Failed generate_all_features()\n");
		return -1;
	}
	timer.take_curr_time();
	MLOG("MedModel::learn() : generating learn matrix time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());

	// Learn Feature cleaners and clean
	timer.start();
	if (learn_and_apply_feature_processors(features) < 0) {
		MERR("MedModel::learn() : ERROR: Failed learn_and_apply_feature_cleaners()\n");
		return -1;
	}
	timer.take_curr_time();
	MLOG("MedModel::learn() : feature processing learn and apply time %g ms :: features crc %08x\n", timer.diff_milisec(), features.get_crc());

	// Select Features
	selector.select(generators);

	// Learn predictor
	timer.start();
	int rc = predictor->learn(features);
	timer.take_curr_time();
	MLOG("MedModel::learn() : model train time: %g ms\n", timer.diff_milisec());

	return rc;
}

//.......................................................................................
// Apply
int MedModel::apply(MedPidRepository& rep, MedSamples& samples) {

	// Set of signals
	get_required_signals(rep.dict);

	vector<int> req_signals;
	for (int signalId : required_signals)
		req_signals.push_back(signalId);

	// Generate features
	MedFeatures features;
	if (generate_all_features(rep, &samples, features) < 0) {
		MERR("MedModel apply() : ERROR: Failed generate_all_features()\n");
		return -1;
	}

	// Clean Features
	if (apply_feature_processors(features) < 0) {
		MERR("MedModel::apply() : ERROR: Failed apply_feature_cleaners()\n");
		return -1;
	}

	// Apply predictor
	if (predictor->predict(features) < 0) {
		MERR("Predictor failed\n");
		return -1;
	}

	samples.insert_preds(features);
	return 0;
}

//.......................................................................................
// Learn rep-cleaning
int MedModel::quick_learn_rep_processors(MedPidRepository& rep, vector<int>& ids) {

	vector<int> rc(rep_processors.size(), 0);

#pragma omp parallel for
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

#pragma omp parallel for
	for (int i = 0; i<generators.size(); i++)
		rc[i] = generators[i]->learn(rep, ids,rep_processors);

	for (auto RC : rc) if (RC < 0)	return -1;
	return 0;
}

//.......................................................................................
int MedModel::generate_all_features(MedPidRepository &rep, MedSamples *samples, MedFeatures &features)
{
	get_required_signals(rep.dict);
	vector<int> req_signals;
	for (int signalId : required_signals)
		req_signals.push_back(signalId);

	// init features attributes
	for (auto& generator : generators)
		generator->init(features);

	// preparing records and features for threading
	int N_tot_threads = omp_get_max_threads();
	MLOG("MedModel::apply() : feature generation with %d threads\n", N_tot_threads);
	vector<PidDynamicRec> idRec(N_tot_threads);
	features.init_all_samples(samples->idSamples);
	int samples_size = (int)features.samples.size();
	for (auto &generator : generators) {
		for (string& name : generator->names)
			features.data[name].resize(samples_size, 0);
	}

	// Loop on ids
	int RC = 0;
#pragma omp parallel for
	for (int j=0; j<samples->idSamples.size(); j++) {
		MedIdSamples& pid_samples = samples->idSamples[j];
		int n_th = omp_get_thread_num();
		int rc = 0;

		// Generate DynamicRec with all relevant signals
		if (idRec[n_th].init_from_rep(std::addressof(rep), pid_samples.id, req_signals, (int)pid_samples.samples.size()) < 0) rc = -1;

		// Apply rep-cleaning
		for (auto& processor : rep_processors)
			if (processor->apply(idRec[n_th], pid_samples) < 0) rc = -1;

		// Generate Features
		for (auto& generator : generators)
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

//.......................................................................................
// Required Signals
void MedModel::get_required_signals(MedDictionarySections& dict) {

	required_signals.clear();

	for (RepProcessor *processor : rep_processors)
		processor->get_required_signal_ids(required_signals,dict);

	int ii = 0;
	for (FeatureGenerator *generator : generators) 
		generator->get_required_signal_ids(required_signals, dict);

}

// Add multi processors
//.......................................................................................
void MedModel::add_rep_processors_set(RepProcessorTypes type, vector<string>& signals) {
	
	RepMultiProcessor *processor = new RepMultiProcessor;
	processor->add_processors_set(type, signals);
	add_rep_processor(processor);

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
	for (unsigned int i = 0; i < generators.size(); i++) {
		for (string& name : generators[i]->names) {
			fprintf(stderr, "Adding %s to processors of type %d\n", name.c_str(),type);
			features.push_back(name);
		}
	}

	add_feature_processors_set(type, features);

}

//.......................................................................................
void MedModel::add_feature_processors_set(FeatureProcessorTypes type, string init_string) {

	vector<string> features;
	for (unsigned int i = 0; i < generators.size(); i++) {
		for (string& name : generators[i]->names)
			features.push_back(name);
	}

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