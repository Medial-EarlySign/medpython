#pragma once
#include <MedProcessTools/MedProcessTools/FeatureGenerator.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>

/**
* calculate drug coverage of prescription time in defined the time window. a value between 0 to 1.
*/
class DrugIntakeGenerator : public FeatureGenerator {
public:

	// Feature Descrption
	string signalName;
	int signalId;

	// parameters (should be serialized)
	int win_from = 0; ///< time window for feature: from is the minimal time before prediciton time
	int win_to = 360000;			///< time window for feature: from is the maximal time before prediciton time
	int time_unit_win = MedTime::Undefined;			///< the time unit in which the windows are given. Default: Undefined
	vector<string> sets;						///< for FTR_CATEGORY_SET_* , the list of sets 
	int time_unit_sig = MedTime::Undefined;		///< the time init in which the signal is given. (set correctly from Repository in learn and Generate)
	string in_set_name = "";					///< set name (if not given - take list of members)
	bool bound_outcomeTime = false;


	// Signal to determine allowed time-range (e.g. current stay/admission for inpatients)
	string timeRangeSignalName = "";
	int timeRangeSignalId;
	TimeRangeTypes timeRangeType = TIME_RANGE_CURRENT;
	int time_unit_range_sig = MedTime::Undefined;		///< the time unit in which the range signal is given. (set correctly from Repository in learn and _generate)

	// helpers
	vector<unsigned char> lut;							///< to be used when generating FTR_CATEGORY_SET_*

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
	int _learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) { time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(signalName)].time_unit; return 0; }

	// generate a new feature
	int _generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data);
	float get_value(PidDynamicRec &rec, int idx, int time, int sig_outcomeTime);

	// Signal Ids
	void set_signal_ids(MedSignals& sigs) { signalId = sigs.sid(signalName); }

	// Init required tables
	void init_tables(MedDictionarySections& dict);

	void get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const;

	// Serialization
	ADD_CLASS_NAME(DrugIntakeGenerator)
	ADD_SERIALIZATION_FUNCS(generator_type, tags, serial_id, win_from, win_to,time_unit_win, signalName, sets, names, req_signals, in_set_name, iGenerateWeights, timeRangeSignalName, timeRangeType , bound_outcomeTime)
};

MEDSERIALIZE_SUPPORT(DrugIntakeGenerator);