#include "MedProcessTools/MedProcessTools/FeatureGenerator.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"

class DrugIntakeGenerator : public FeatureGenerator {
public:

	// Feature Descrption
	string signalName;
	int signalId;

	// parameters (should be serialized)
	int win_from = 0, win_to = 360000;			// time window for feature: date-win_to <= t < date-win_from. both are AFTER 
	int time_unit_win = MedTime::Undefined;			// the time unit in which the windows are given. Default: Undefined
	vector<string> sets;						// for FTR_CATEGORY_SET_* , the list of sets 
	int time_unit_sig = MedTime::Undefined;		// the time init in which the signal is given. (set correctly from Repository in learn and Generate)
	string in_set_name = "";					// set name (if not given - take list of members)

	// helpers
	vector<unsigned char> lut;							// to be used when generating FTR_CATEGORY_SET_*

	// Constructor/Destructor
	DrugIntakeGenerator() : FeatureGenerator() { init_defaults(); };
	~DrugIntakeGenerator() {};

	/// The parsed fields from init command.
	/// @snippet DrugIntakeGenerator.cpp DrugIntakeGenerator::init
	int init(map<string, string>& mapper);
	void init_defaults();

	// Naming
	void set_names();

	// Copy
	virtual void copy(FeatureGenerator *generator) { *this = *(dynamic_cast<DrugIntakeGenerator *>(generator)); }

	// Learn a generator
	int _learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) { time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(signalName)].time_unit; return 0; }

	// generate a new feature
	int Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num);
	float get_value(PidDynamicRec &rec, UniversalSigVec &usv, int time, int sig_outcomeTime);

	// Signal Ids
	void set_signal_ids(MedDictionarySections& dict) { signalId = dict.id(signalName); }

	// Init required tables
	void init_tables(MedDictionarySections& dict);

	// Serialization
	ADD_SERIALIZATION_FUNCS(generator_type, tags, serial_id, win_from, win_to,time_unit_win, signalName, sets, names, req_signals, in_set_name, iGenerateWeights)
};

MEDSERIALIZE_SUPPORT(DrugIntakeGenerator);