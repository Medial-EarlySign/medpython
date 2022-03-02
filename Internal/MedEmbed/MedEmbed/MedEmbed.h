#pragma once
//
// MedEmbed :
// Contains a set of structures and tools to allow for preparing a training of an embedding,
// and using it to create features by generating the sparse categorial line and running it
// through a pre trained model (currently we use keras models)
// 

#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedTime/MedTime/MedTime.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedSparseMat/MedSparseMat/MedSparseMat.h>

#include <unordered_map>
#include <unordered_set>

enum EmbeddedCodeType { ECTYPE_CATEGORIAL = 0, ECTYPE_CONTINUOUS, ECTYPE_AGE, ECTYPE_DUMMY , ECTYPE_MODEL, ECTYPE_UNDEFINED};

//======================================================================================
// needed info on each signal participating
// containing all the information needed to create the codes
//======================================================================================
class EmbeddingSig : public SerializableObject {

public:
	string sig;
	EmbeddedCodeType type = ECTYPE_CATEGORIAL;
	int add_hierarchy = 1;
	int do_shrink = 1; // keeping the bit of wheather to shrink this one when shrinking is done.
	int do_counts = 1; // collecting data as counts if =1 , or just as 0/1 if =0

	// next is used for categorial signals. It holds the string names of the categories we want to embed.
	// initialization is with the categories= parameter , it is a comma (,) separeted list of names (best initialize with the "list:" option)
	vector<string> categories_to_embed;
	string regex_filter = "";

	// for categorial: ranges to use (empty is all). The sequence is:
	// a value and all its hierarchy (if asked for) are created. Then only those within the asked for ranges are actually added and used.

	// for non categorial: defines value ranges r[i][0] <= val < r[i][1] means it will be counted into original category i
	// This allows creating categorial features out of continous signals
	vector<vector<float>> ranges;

	int time_chan = 0;
	int val_chan = 0;
	int win_from = 0;
	int win_to = 365;
	int sig_time_unit = global_default_time_unit;
	int win_time_unit = global_default_windows_time_unit;

	// next are for the model type
	// a model can be initiated with a model file (you have to pretrain it) .
	// It is highly recommended to use models creating matrices that are imputed AND normalized.
	// the features generated will be copied into the sparse matrix
	// features generated like this are never shrunk
	MedModel *model = NULL; // only object needed for serializations
	string model_file = "";
	vector<string> model_req_sigs;
	vector<string> model_features_names;
	vector<float *> feat_ptrs; // after matrix was created, holds pointers to each model feature for faster access
	map<pair<int, int>, int> pidtime2idx; // since order in batch prep is not necessarily the same as the one in the generation, we need this mapping


	// for categorials only : all sets of a value : we need to build this table before each apply using the given dictionary (and the known categories)
	unordered_map<int, vector<int>> sig_members2sets;

	// for categorials: after limiting to sets in range only ( = orig values)
	unordered_map<int, vector<int>> sig_members2sets_in_range;

	// orig to code and name
	map<string, int> Name2Id; // Only relevant for categorial cases : keeping the sub Name2Id table that was used to build the Orig2X tables.
							  // when using a new different repository we need to translate the ids to these ones.
							  // this is done in the categs_convert table
	vector<int> categ_convert; // relevant also for the shrunk case
	map<int, int> Orig2Code;
	map<int, string> Orig2Name;

	// orig to shrunk code
	map<int, int> Orig2ShrunkCode;

	void clear_tables() { sig_members2sets.clear(); sig_members2sets_in_range.clear(); Name2Id.clear(); Orig2Code.clear(); Orig2Name.clear(); Orig2ShrunkCode.clear(); }

	// simple API's

	// appends the orig values to the given codes vector , returns number of elements added
	int get_categ_orig(int val, vector<int> &codes) const;
	
	// appends the codes to the given codes vector , returns number of elements added
	int get_categ_codes(int val, vector<int> &codes, int use_shrink = 1) const;

	// appends the shrunk codes to the given codes vector , returns number of elements added
	int get_categ_shrunk_codes(int val, vector<int> &codes) const;

	// appends the orig values to the given codes vector , returns number of elements added
	int get_continuous_orig(float val) const;

	// appends the codes to the given codes vector , returns number of elements added
	int get_continuous_codes(float val, int use_shrink = 1) const;

	// appends the shrunk codes to the given codes vector , returns number of elements added
	int get_continuous_shrunk_codes(float val) const;

	// helper and not needed to serialize params
	int sid = -1;

	// initialization from string
	int init(map<string, string>& _map);


	// initializing a categorial case : need to get a dictionary and init  : Orig2Code, Orig2Name (see also init_categotial_tables)
	// this is needed in order to make the embedding independent of the actual values given in the directory and rely on names only.
	// This has the potential of allowing to transfer embeddings between different data sets, as long as they use the same signal names with the same category names in the dictionary.
	int init_categorial(MedDictionarySections &dict, int &curr_code);

	// the next is special for the categorial case:
	// We need to initialize the Name2Id table (only if it is not empty !! , as it may be full from the original mapping that was used to build the Orig tables)
	// once we have that table, we need to initialize the following tables:
	// sig_members2sets, sig_members2sets_in_range, and also categ_convert
	int init_categorial_tables(MedDictionarySections &dict);

	// initialize a continous or age case : preparing the Orig2X tables based on the given ranges.
	int init_continous(int &curr_code);

	// initialize a dummy case : simple constant variable always added to make sure we have at least one entry per sample (helps in some cases)
	int init_dummy();


	// actually collecting matrix lines
	int add_sig_to_lines(UniversalSigVec &usv, int pid, int time, int use_shrink, map<int, map<int, float>> &out_lines) const;
	int get_codes(UniversalSigVec &usv, int pid, int time, int use_shrink, vector<int> &codes) const;
	int add_codes_to_line(vector<int> &codes, map<int, float> &out_line) const;
	int add_to_line(UniversalSigVec &usv, int pid, int time, int use_shrink, map<int, float> &out_line) const;

	// preparing a batch of model results (will also initialize the feat_ptr vector, and the pidtime2idx map)
	int prep_model_batch(MedPidRepository &rep, MedSamples &samples);


	EmbeddedCodeType type_name_to_code(string name);

	string print_to_string(int verbosity);

	// next can be used after shrinking was done
	// it keeps the minimal structures needed in order to allow matrix creation and lower scheme file size.
	int minimize();

	ADD_CLASS_NAME(EmbeddingSig)
	ADD_SERIALIZATION_FUNCS(sig, type, add_hierarchy, do_shrink, ranges, time_chan, val_chan, win_from, win_to, categories_to_embed, Name2Id, Orig2Code, Orig2Name, Orig2ShrunkCode, model)
};


//============================================================================================================================
// EmbedMatsCreator : major class for creating sparse embedding matrices for a given setup + list of times, window_lens, etc
//============================================================================================================================
class EmbedMatCreator : public SerializableObject {

public:
	vector<string> sigs_to_load;

	int rep_time_unit = MedTime::Date;
	int win_time_unit = MedTime::Days;
	int byear_time_unit = MedTime::Years;

	vector<EmbeddingSig> embed_sigs;			// containing all the information on each the sigs to embed

	// general high level operations

	// prepare needs to be run before creating a matrix for the first time:
	// (1) initializes the sigs_to_load vector
	// (2) initializes the embed_sigs objects up to the pre shrinking stage
	// When starting with a serialized object there's no need to call this one.
	int prepare(MedPidRepository &rep);
	
	// adding all the needed lines for a pid. Written for a dynamic record, to allow easy connection to MedProcessTools
	int add_pid_lines(PidDynamicRec &pdr, MedSparseMat &smat, vector<int> &times, int use_shrink);
	int get_pid_out_line(PidDynamicRec &pdr, int ver, int time, int use_shrink, map<int, float> &out_line);

	// another api to generate a matrix given a list of pids and times, that promises the SAME order as in the given input
	// the input pair vector has pids on first, and times on second
	// works directly through the rep (not the PidDynamicRec path)
	int get_sparse_mat(MedPidRepository &rep, vector<pair<int, int>> &pids_times, int use_shrink, MedSparseMat &smat);

	// sometimes easier to use, BUT the ORDER of lines in the matrix is the order of normalized samples,
	// this makes it a problem when needing to produce a matrix with a different order for lines.
	int get_sparse_mat(MedPidRepository &rep, MedSamples &samples , int use_outcome_time, int use_shrink, MedSparseMat &smat);

	// helper for es preparation 
	void prep_memebers_to_sets(MedPidRepository &rep, EmbeddingSig &es);

	// shrinking calculation
	// gets an smat that had been produced with the non shrinked dictionary,
	// then selects the columns that will stay (es with do_shrink = 0, or those with at least min_p-max_p rows containing it.
	int get_shrinked_dictionary(MedSparseMat &smat, float min_p, float max_p);

	// apply shrinking to a given matrix
	// (other better option is to build it with the use_shrink=1 flag)
	int shrink_mat(MedSparseMat &smat, MedSparseMat &shrunk_smat);

	// initialization from string
	int init(map<string, string>& _map);

	// needed before we start using the class on a specific rep, but AFTER params and embed_sigs were initialized.
	void init_sids(MedSignals &sigs);

	// next must be called after coming from serialization, at the moment we get hold of dict.
	void init_tables(MedDictionarySections &dict) { for (auto &es : embed_sigs) es.init_categorial_tables(dict); }


	// API to write the dictionary to a file, to have a readable interpretation of the codes.
	int write_dict_to_file(string fname, int only_shrink);


	// printing object to string
	string print_to_string(int verbosity);

	// minimizing size of shrunk categorials for smaller scheme files
	// if this is run before serialization one will only be able to create the shrunk version (which is what is needed...)
	int minimize() { for (auto &es : embed_sigs) es.minimize(); return 0; };

	// next is needed in order to allow for batch preparations of model es
	void prep_models_batches(MedPidRepository &rep, MedSamples &samples) { for (auto &es : embed_sigs) es.prep_model_batch(rep, samples); }


	ADD_CLASS_NAME(EmbedMatCreator)
	ADD_SERIALIZATION_FUNCS(sigs_to_load, rep_time_unit, win_time_unit, byear_time_unit, embed_sigs)

private:
	// helpers
	int curr_code = 1; // needed in codes allocation process

};


//============================================================================================================================
struct EmbedXYRecord {
	int pid;
	int x_time;
	int y_time;
};


//============================================================================================================================
// train matrices creation class
//============================================================================================================================
class EmbedTrainCreator : public SerializableObject
{

public:

	// params
	string x_params;
	string y_params;
	int use_same_dictionaries = 1; // if on : x,y must have the SAME es order, the same Orig2Code, Orig2Name in each, and we will copy 
								   // the x shrinking dictionary to y.								

	// next params are to generate an xy list
	int min_time = 20060101;
	int max_time = 20160101;
	int min_age = 30;
	int max_age = 100;
	int npoints_per_pid = 1;
	float min_p = (float)0.001;
	float max_p = (float)0.95;
	vector<int> time_dist_range;
	vector<int> time_dist_points ={ -365, 0, 365 };

	// general technical params needed for production
	int rep_time_unit = MedTime::Date;
	int win_time_unit = MedTime::Days;
	int byear_time_unit = MedTime::Years;
	string prefix = "smat";

	// matrices params
	float p_train = (float)0.8;

	EmbedMatCreator x_creator, y_creator;

	// generate x,y matrices for a given xy-file, and write them to files (including dictionaries)
	int generate_from_xy_file(string xy_fname, string rep_fname, string out_prefix);

	// generate an xy list and write it to file, input is a list of pids and a repository
	int generate_xy_list(string xy_fname, string pids_fname, string rep_fname);

	// helpers: read/write a file of <pid> <xtime> <ytime> records
	int read_xy_records(string xy_fname, vector<EmbedXYRecord> &xy);
	int write_xy_records(string xy_fname, vector<EmbedXYRecord> &xy);

	// init
	int init(map<string, string>& _map);
};




//=================================================================
// Joining the MedSerialize Wagon
//=================================================================

MEDSERIALIZE_SUPPORT(EmbeddingSig)
MEDSERIALIZE_SUPPORT(EmbedMatCreator)


