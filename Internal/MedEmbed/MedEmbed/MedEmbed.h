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
	// initialization is with the categories= parameter , it is a comma (,) separeted list of names, if a name starts with file: then
	// it means that what's after the : is a file name containing a list of the categories to use.
	// When creating features the order is : start with the categorial signal value, and add all its hierarchy if asked for, then check if any of those are
	// in the categories_to_embed list, and if so, create a feature for it.
	unordered_set<string> categories_to_embed;

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
	MedModel *model = NULL;
	string model_file = "";
	vector<int> model_sids;
	MedFeatures feat;
	map<pair<int, int>, int> pid_time2idx;

	int get_feat_for_model(MedPidRepository &rep, vector<pair<int, int>> &pids_times);

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
	int get_categ_orig(int val, vector<int> &codes);
	
	// appends the codes to the given codes vector , returns number of elements added
	int get_categ_codes(int val, vector<int> &codes, int use_shrink = 1);

	// appends the shrunk codes to the given codes vector , returns number of elements added
	int get_categ_shrunk_codes(int val, vector<int> &codes);

	// appends the orig values to the given codes vector , returns number of elements added
	int get_continuous_orig(float val);

	// appends the codes to the given codes vector , returns number of elements added
	int get_continuous_codes(float val, int use_shrink = 1);

	// appends the shrunk codes to the given codes vector , returns number of elements added
	int get_continuous_shrunk_codes(float val);

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
	int add_sig_to_lines(UniversalSigVec &usv, int pid, int time, int use_shrink, map<int, map<int, float>> &out_lines);

	EmbeddedCodeType type_name_to_code(string name);

	string print_to_string(int verbosity);

	ADD_CLASS_NAME(EmbeddingSig)
	ADD_SERIALIZATION_FUNCS(sig, type, add_hierarchy, do_shrink, ranges, time_chan, val_chan, win_from, win_to, Name2Id, Orig2Code, Orig2Name, Orig2ShrunkCode, model)
};


//============================================================================================================================
// EmbedMatsCreator : major class for creating sparse embedding matrices for a given setup + list of times, window_lens, etc
//============================================================================================================================
class EmbedMatCreator : public SerializableObject {

public:
	vector<string> sigs_to_load;
	int start_sid = -1, end_sid = -1, death_sid = -1;

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
	
	// major helper func:
	// gets a usv for a signal, a time and a win_len, and adds the relevant data to output lines
	// later and elsewhere these lines will be added atomically and by order to a sparse mat (to allow for easy threading code)
	int add_model_feats_to_lines(EmbeddingSig &es, PidDynamicRec &pdr, vector<int> &times, int use_shrink, map<int, map<int, float>> &out_lines);
	int add_model_feats_to_lines(EmbeddingSig &es, int pid, vector<int> &times, int use_shrink, map<int, map<int, float>> &out_lines);
	
	// adding all the needed lines for a pid. Written for a dynamic record, to allow easy connection to MedProcessTools
	int add_pid_lines(PidDynamicRec &pdr, MedSparseMat &smat, vector<int> &times, int use_shrink);

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
	void init_sids(MedPidRepository &rep);

	// next must be called after coming from serialization, at the moment we get hold of dict.
	void init_tables(MedDictionarySections &dict) { for (auto &es : embed_sigs) es.init_categorial_tables(dict); }


	// API to write the dictionary to a file, to have a readable interpretation of the codes.
	int write_dict_to_file(string fname, int only_shrink);


	// printing object to string
	string print_to_string(int verbosity);

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


