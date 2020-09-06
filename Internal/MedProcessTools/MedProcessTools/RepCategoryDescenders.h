#ifndef _REP_DESCENDERS_H_
#define _REP_DESCENDERS_H_

#include "RepProcess.h"
using namespace std;

/**
* RepCategoryDescenders - creates all descenders values for each values in the original signal value
*/
class RepCategoryDescenders : public RepProcessor {
private:
	int signalId;	///< id of signal 
	unordered_map<int, vector<int>> _set2Members; ///< for hierarchy
public:
	string signalName = ""; ///< name of signal
	int val_channel = 0; ///< value channel


	RepCategoryDescenders() { processor_type = REP_PROCESS_CATEGORY_DESCENDERS; }

	/// @snippet RepCategoryDescenders.cpp RepCategoryDescenders::init
	int init(map<string, string>& mapper);
	void init_tables(MedDictionarySections& dict, MedSignals& sigs);

	void set_signal(const string& _signalName) { signalId = -1; signalName = _signalName; init_lists(); }

	void init_lists();
	void set_signal_ids(MedSignals& sigs);
	// Applying
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points : Should be 
	/// implemented for all inheriting classes 
	/// </summary>
	int _apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat);


	ADD_CLASS_NAME(RepCategoryDescenders)
		ADD_SERIALIZATION_FUNCS(signalName, val_channel, req_signals, aff_signals)
};

MEDSERIALIZE_SUPPORT(RepCategoryDescenders)

#endif

