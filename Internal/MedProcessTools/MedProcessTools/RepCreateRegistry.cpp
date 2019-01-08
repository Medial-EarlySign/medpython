#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "RepProcess.h"

//=======================================================================================
// RepCreateRegistry for creating repositories as signals
//=======================================================================================

/// Init from map
//=======================================================================================
int RepCreateRegistry::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepCreateRegistry::init]
		if (field == "registry") registry_name = entry.second;
		else if (field == "names") boost::split(names, entry.second, boost::is_any_of(","));
		else if (field == "signals") boost::split(signals, entry.second, boost::is_any_of(","));
		else if (field == "signals_time_unit") signals_time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "rp_type") {}
		else MTHROW_AND_ERR("Error in RepCreateRegistry::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepCreateRegistry::init]
	}

	registry = name2type.at(registry_name);

	if (signals_time_unit == -1 || signals_time_unit == MedTime::Undefined) {
		MWARN("Warning in RepCreateRegistry::init - using signals_time_unit = Days as defualt time unit\n");
		signals_time_unit = MedTime::Days;
	}

	virtual_signals = type2Virtuals.at(registry);
	if (!names.empty()) {
		if (names.size() != type2Virtuals.at(registry).size())
			MTHROW_AND_ERR("Wrong number of names supplied for RepCreateRegistry::%s - supplied %zd, required %zd\n", registry_name.c_str(), names.size(),
				type2Virtuals.at(registry).size());
		for (size_t i = 0; i < names.size(); i++)
			virtual_signals[i].first = names[i];
	}

	// required/affected signals
	init_lists();
}

/// Required/Affected signals
//=======================================================================================
void RepCreateRegistry::init_lists() {

	req_signals.clear();
	for (string signalName : type2reqSigs.at(registry))
		req_signals.insert(signalName);

	aff_signals.clear();
	for (auto& rec : virtual_signals)
		aff_signals.insert(rec.first);
}


// making sure V_ids and sigs_ids are initialized
//=======================================================================================
void RepCreateRegistry::init_tables(MedDictionarySections& dict, MedSignals& sigs) {

	virtual_ids.clear();
	for (auto& rec : virtual_signals)
		virtual_ids.push_back(sigs.sid(rec.first));

	req_signal_ids.clear();
	for (auto &rsig : signals) {
		int sid = sigs.sid(rsig);
		sig_ids_s.insert(sid);
		sig_ids.push_back(sid);
	}

	aff_signal_ids.clear();
	aff_signal_ids.insert(virtual_ids.begin(), virtual_ids.end());

	// Registry specific tables
	if (registry == REP_REGISTRY_HT)
		init_ht_registry_tables(dict, sigs);
	else if (registry == REP_REGISTRY_DM)
		init_dm_registry_tables(dict, sigs);

	// Allocate
	usvs.resize(sig_ids.size());
	all_v_times.resize(virtual_ids.size());
	all_v_vals.resize(virtual_ids.size());
	final_sizes.reserve(virtual_ids.size());
}

// Applying
/// <summary> apply processing on a single PidDynamicRec at a set of time-points : Should be implemented for all inheriting classes </summary>
int RepCreateRegistry::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {

	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int len;
	allVersionsIterator vit(rec, sig_ids_s);

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		for (size_t isig = 0; isig < sig_ids.size(); isig++)
			rec.uget(sig_ids[isig], iver, usvs[isig]);

		if (registry == REP_REGISTRY_HT)
			ht_registry_apply(rec, time_points, iver);
		else if (registry == REP_REGISTRY_DM)
			dm_registry_apply(rec, time_points, iver);
		

		// pushing virtual data into rec
		for (size_t ivir = 0; ivir<virtual_ids.size(); ivir++)
			rec.set_version_universal_data(virtual_ids[ivir], iver, &(all_v_times[ivir][0]), &(all_v_vals[ivir][0]), final_sizes[ivir]);
	}



}

// Registry-Specific functions : HyperTension
//=======================================================================================
void RepCreateRegistry::init_ht_registry_tables(MedDictionarySections& dict, MedSignals& sigs) {

}

void RepCreateRegistry::ht_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver) {

}


// Registry-Specific functions : Diabetes
//=======================================================================================
void RepCreateRegistry::init_dm_registry_tables(MedDictionarySections& dict, MedSignals& sigs) {

}

void RepCreateRegistry::dm_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver) {

}