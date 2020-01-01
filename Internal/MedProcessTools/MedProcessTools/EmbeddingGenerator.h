#pragma once
#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>
#include <MedEmbed/MedEmbed/MedEmbed.h>
#include <MedEmbed/MedEmbed/ApplyKeras.h>
#include <MedSparseMat/MedSparseMat/MedSparseMat.h>

//
// Allow the usage of a trained embedding to generate features
// needed are the .scheme and .layers file created using the Embedding app and Keras script.
// This version only APPLIES and does not LEARN the embedding.
//

class EmbeddingGenerator : public FeatureGenerator {
public:

	// Feature Descrption
	string signalName;
	int signalId;

	// parameters
	string f_scheme = "";
	string f_layers = "";

	// serialized parameters
	string name_prefix = "Embedding"; // used when naming the features, names will be FTR_<i>.<name_prefix>.col_<j>
	int to_layer = 9; // default embedding layer in our keras generated models

	// basic embedding structures, serialized (init from scheme, layers files upon creation from init() )
	EmbedMatCreator emc;
	ApplyKeras embedder;

	// helper variables
	int e_dim = 0; // we get it from embedder, but serialize it for simpler code


	// Signal to determine allowed time-range (e.g. current stay/admission for inpatients)

														// Constructor/Destructor
	EmbeddingGenerator() : FeatureGenerator() {};
	~EmbeddingGenerator() {};

	/// The parsed fields from init command.
	/// @snippet DrugIntakeGenerator.cpp DrugIntakeGenerator::init
	int init(map<string, string>& mapper);
	//void init_defaults();

	// Naming
	void set_names();

	// preparing a batch result of model results if an internal MedModel generator is used
	void prepare(MedFeatures & features, MedPidRepository& rep, MedSamples& samples);

	// Learn a generator
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors);

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	int generate_by_rec(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);

	// Signal Ids
	void set_signal_ids(MedSignals& sigs);

	// Init required tables
	void init_tables(MedDictionarySections& dict);

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

	// in case of selection
	int filter_features(unordered_set<string>& validFeatures) { return 1; }; // TODO: improve to push only selected columns

	// Serialization
	ADD_CLASS_NAME(EmbeddingGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, tags, names, req_signals, signalName, iGenerateWeights, name_prefix, to_layer, e_dim, emc, embedder)
};

MEDSERIALIZE_SUPPORT(EmbeddingGenerator);
