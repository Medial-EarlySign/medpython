// A model = repCleaner + featureGenerator + featureCleaner + MedPredictor

#ifndef _MED_MODEL_H_
#define _MED_MODEL_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/RepProcess.h"
#include "MedProcessTools/MedProcessTools/FeatureProcess.h"
#include "MedProcessTools/MedProcessTools/FeatureGenerator.h"
#include "MedProcessTools/MedProcessTools/FeatureSelector.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"

//.......................................................................................
//.......................................................................................
// A model = repCleaner + featureGenerator + featureCleaner + MedPredictor
// Feature selection can be applied in the learning stage but only the set of selected features should be maintained (?)
//.......................................................................................
//.......................................................................................
class MedModel : public SerializableObject {
public:

	// Repostiroy-level cleaners; to be applied sequentially 
	vector<RepProcessor *> rep_processors;

	// Feature Generators 
	vector<FeatureGenerator *> generators;

	// Features-level cleaners; to be applied sequentially 
	vector<FeatureProcessor *> feature_processors;

	// Predictor
	MedPredictor *predictor;

	// Learning samples
	MedSamples *LearningSet;

	// Safe Mode for train/test intersection
	int safe_mode;

	// All required signals
	unordered_set<int> required_signals;

	// Constructor/Destructor
	MedModel() { safe_mode = 0;  };
	~MedModel() {};

	// Add Rep Processors
	void add_rep_processor(RepProcessor *processor) { rep_processors.push_back(processor); };
	void add_rep_processors_set(RepProcessorTypes type, vector<string>& signals);
	void add_rep_processors_set(RepProcessorTypes type, vector<string>& signals, string init_string);

	// Add Feature Generators
	void add_feature_generator(FeatureGenerator *generator) { generators.push_back(generator); }
	void add_feature_generators(FeatureGeneratorTypes type, vector<string>& signals);
	void add_feature_generators(FeatureGeneratorTypes type, vector<string>& signals, string init_string);
	void add_feature_generator(FeatureGeneratorTypes type, string& signal) { vector<string> signals(1,signal) ; add_feature_generators(type, signals); }
	void add_feature_generators(FeatureGeneratorTypes type, string& signal, string init_string) { vector<string> signals(1,signal) ; add_feature_generators(type, signals, init_string); }

	void add_feature_generators(string& name, vector<string>& signals) { add_feature_generators(ftr_generator_name_to_type(name), signals); }
	void add_feature_generators(string& name, vector<string>& signals, string init_string) { add_feature_generators(ftr_generator_name_to_type(name), signals, init_string); }
	void add_feature_generator(string& name, string& signal) { vector<string> signals(1,signal) ; add_feature_generators(name, signals); }
	void add_feature_generators(string& name, string& signal, string init_string) { vector<string> signals(1,signal) ; add_feature_generators(name, signals, init_string); }

	void add_age() {generators.push_back(new AgeGenerator);}
	void add_gender() { generators.push_back(new GenderGenerator); }


	// Add Feature Processors
	void add_feature_processor(FeatureProcessor *processor) { feature_processors.push_back(processor); };

	void add_feature_processors_set(FeatureProcessorTypes type);
	void add_feature_processors_set(FeatureProcessorTypes type, string init_string);
	void add_feature_processors_set(FeatureProcessorTypes type, vector<string>& features);
	void add_feature_processors_set(FeatureProcessorTypes type, vector<string>& features, string init_string);

	void add_normalizers() {add_feature_processors_set(FTR_PROCESS_NORMALIZER);}
	void add_normalizers(string init_string) { add_feature_processors_set(FTR_PROCESS_NORMALIZER, init_string); }
	void add_normalizers(vector<string>& features) { add_feature_processors_set(FTR_PROCESS_NORMALIZER, features); }
	void add_normalizers(vector<string>& features, string init_string) { add_feature_processors_set(FTR_PROCESS_NORMALIZER, features, init_string); }

	void add_imputers() { add_feature_processors_set(FTR_PROCESS_IMPUTER); }
	void add_imputers(string init_string) { add_feature_processors_set(FTR_PROCESS_IMPUTER, init_string); }
	void add_imputers(vector<string>& features) { add_feature_processors_set(FTR_PROCESS_IMPUTER, features); }
	void add_imputers(vector<string>& features, string init_string) { add_feature_processors_set(FTR_PROCESS_IMPUTER, features, init_string); }

	// Add Predictor
	void set_predcitor(MedPredictor *_predictor) { predictor = _predictor; };
	void set_predictor(MedPredictorTypes type) { predictor = MedPredictor::make_predictor(type); }
	void set_predictor(string name) { predictor = MedPredictor::make_predictor(name); }
	void set_predictor(MedPredictorTypes type, string init_string) { predictor = MedPredictor::make_predictor(type,init_string); }
	void set_predictor(string name, string init_string) { predictor = MedPredictor::make_predictor(name,init_string); }

	// signal ids
	void get_required_signals(MedDictionarySections& dict);
	void get_affected_signals(MedDictionarySections& dict);
	void init_signal_ids(MedDictionarySections& dict);

	// Apply
	int learn(MedPidRepository& rep, MedSamples* samples);
	int learn(MedPidRepository& rep, MedSamples* samples, FeatureSelector& selector);
	int apply(MedPidRepository& rep, MedSamples& samples) ;

	// De(Serialize)
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
		
	int quick_learn_rep_processors(MedPidRepository& rep, vector<int>& ids);
	int learn_rep_processors(MedPidRepository& rep, vector<int>& ids);
	int learn_feature_generators(MedPidRepository &rep, MedSamples *learn_samples);
	int generate_all_features(MedPidRepository &rep, MedSamples *learn_samples, MedFeatures &features);
	int learn_and_apply_feature_processors(MedFeatures &features);
	int apply_feature_processors(MedFeatures &features);

};

#endif
