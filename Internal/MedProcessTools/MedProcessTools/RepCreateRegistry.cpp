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
		else if (field == "dm_drug_sig") dm_drug_sig = entry.second;
		else if (field == "dm_drug_sets") boost::split(dm_drug_sets, entry.second, boost::is_any_of(","));
		else if (field == "dm_diagnoses_sig") dm_diagnoses_sig = entry.second;
		else if (field == "dm_diagnoses_sets") boost::split(dm_diagnoses_sets, entry.second, boost::is_any_of(","));

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

	allVersionsIterator vit(rec, sig_ids_s);
	vector<UniversalSigVec> usvs;

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		for (size_t isig = 0; isig < sig_ids.size(); isig++)
			rec.uget(sig_ids[isig], iver, usvs[isig]);

		if (registry == REP_REGISTRY_HT)
			ht_registry_apply(rec, time_points, iver);
		else if (registry == REP_REGISTRY_DM)
			dm_registry_apply(rec, time_points, iver, usvs);
		

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
void RepCreateRegistry::init_dm_registry_tables(MedDictionarySections& dict, MedSignals& sigs)
{

	int i = 0;
	for (auto &rsig : signals) {
		if (rsig == dm_drug_sig) dm_drug_idx = i;
		if (rsig == dm_diagnoses_sig) dm_diagnoses_idx = i;
		if (rsig == dm_glucose_sig) dm_glucose_idx = i;
		if (rsig == dm_hba1c_sig) dm_hba1c_idx = i;
		i++;
	}

	// lookup tables
	if (dm_drug_sig != "")	dict.dicts[dict.section_id(dm_drug_sig)].prep_sets_lookup_table(dm_drug_sets, dm_drug_lut);
	if (dm_diagnoses_sig != "") dict.dicts[dict.section_id(dm_diagnoses_sig)].prep_sets_lookup_table(dm_diagnoses_sets, dm_diagnoses_lut);

}

void RepCreateRegistry::dm_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec> &usvs)
{
	vector<RegistryEvent> evs;

	int time = -1;
	if (time_points.size() > 0) time = time_points[iver];

	// step 1 collect events

	// glucose events
	if (dm_glucose_idx >= 0) {

		UniversalSigVec &glu_usv = usvs[dm_glucose_idx];

		for (int i = 0; i < glu_usv.len; i++) {
			int i_time = glu_usv.Time(i);
			if (time > 0 && time > i_time) break;
			float i_val = glu_usv.Val(i);
			int severity = 0;
			if (i_val > 99.99f) severity = 1;	// 2 tests enough to define pre diabetic
			if (i_val > 109.99f) severity = 2;	// single test enough to define pre diabetic
			if (i_val > 125.99f) severity = 3;	// 2 tests enough to define diabetic
			if (i_val > 199.99f) severity = 4;  // 1 test enough to define diabetic

			evs.push_back(RegistryEvent(i_time, REG_EVENT_DM_GLUCOSE, i_val, severity));
		}
	}

	// HbA1C events
	if (dm_hba1c_idx >= 0) {
		UniversalSigVec &hba1c_usv = usvs[dm_hba1c_idx];

		for (int i = 0; i < hba1c_usv.len; i++) {
			int i_time = hba1c_usv.Time(i);
			if (time > 0 && time > i_time) break;
			float i_val = hba1c_usv.Val(i);
			int severity = 0;
			//if (i_val > 99.99f) severity = 1;	// 2 tests enough to define pre diabetic
			if (i_val > 5.69f) severity = 2;	// single test enough to define pre diabetic
			if (i_val > 6.49f) severity = 3;	// 2 tests enough to define diabetic
			if (i_val > 7.99f) severity = 4;  // 1 test enough to define diabetic

			evs.push_back(RegistryEvent(i_time, REG_EVENT_DM_HBA1C, i_val, severity));
		}
	}

	// Drug Events
	if (dm_drug_idx >= 0) {
		UniversalSigVec &drug_usv = usvs[dm_drug_idx];

		for (int i = 0; i < drug_usv.len; i++) {
			int i_time = drug_usv.Time(i);
			if (time > 0 && time > i_time) break;
			float i_val = drug_usv.Val(i);
			if (dm_drug_lut[(int)i_val]) {
				int severity = 4; // currently the first diabetic drug usage makes you diabetic for life.... this is extreme, but given this, we only need the first.
				evs.push_back(RegistryEvent(i_time, REG_EVENT_DM_DRUG, 1, severity)); 
				break;
			}

		}
	}


	// Diagnoses Events
	if (dm_diagnoses_idx >= 0) {
		UniversalSigVec &diag_usv = usvs[dm_diagnoses_idx];

		for (int i = 0; i < diag_usv.len; i++) {
			int i_time = diag_usv.Time(i);
			if (time > 0 && time > i_time) break;
			float i_val = diag_usv.Val(i);
			if (dm_diagnoses_lut[(int)i_val]) {
				int severity = dm_diagnoses_severity; 
				evs.push_back(RegistryEvent(i_time, REG_EVENT_DM_DIAGNOSES, 1, severity));
				if (dm_diagnoses_severity >= 4) break;
			}

		}
	}

	// collection of events done

	// sorting events
	sort(evs.begin(), evs.end(), [](const RegistryEvent &v1, const RegistryEvent &v2) { return v1.time < v2.time; });

	// applying rules
	vector<pair<int, int>> ranges(3, pair<int,int>(-1,-1)); // 0: for healthy, 1: for prediabetic , 2: for diabetic

	for (int j = 0; j < evs.size(); j++) {

		auto &ev = evs[j];

		// rules:
		// (1) to be Diabetic: (a) a single severity 4 (b) adjacent or within 2 years: 2 severity 3 (real mode: the second time, biological mode: the first time)
		// (2) to be PreDiabetic: (a) a single 2,3,4 (b) adjacent or within 2 years: 2 severity 1 (real mode: the second, bio_mode: the first)
		// (3) to be Healthy: severity 0 , or severity 1 after 2 years of not developing into diabetic or pre.

		if (ranges[2].first > 0) {
			// person is diabetic, we have nothing else to do
			if (time > 0) {
				ranges[2].second = time; // Diabetic up to current time point.
				break;
			}
			else {
				ranges[2].second = evs.back().time;
				break;
			}
		}

		if (ev.event_severity == 4) {
			ranges[2].first = ev.time;
			ranges[2].second = time;
			continue;
		}

		if (ev.event_severity == 3) {

			// need to check for severity 3 2 years back
			int back_time = med_time_converter.add_subtract_time(ev.time, signals_time_unit, -730, MedTime::Days);
			int found = 0;
			for (int k = j - 1; k >= 0; k--) {
				if (evs[k].time < back_time) break;
				if (evs[k].event_severity == 3) {
					found = 1;
					break;
				}
			}
			if (found) {
				ranges[2].first = ev.time;
				ranges[2].second = time;
				continue;
			}

		}

		if (ranges[1].first > 0) {
			// the person is pre diabetic and the current severity is not yetleading to diabetic
			// therefore we simply elongate the prediabetic period
			ranges[1].second = ev.time;
			continue;
		}

		if (ev.event_severity >= 2) { // >= as may be severity 3 that wasn't yet enough for becoming diabetic

			ranges[1].first = ev.time;
			ranges[1].second = ev.time;
			continue;
		}

		if (ev.event_severity == 1) {
			int back_time = med_time_converter.add_subtract_time(ev.time, signals_time_unit, -730, MedTime::Days);
			int found = 0;
			for (int k = j - 1; k >= 0; k--) {
				if (evs[k].time < back_time) break;
				if (evs[k].event_severity == 1) {
					found = 1;
					break;
				}
			}
			if (found) {
				ranges[1].first = ev.time;
				ranges[1].second = ev.time;
				continue;
			}
		}

		if (ev.event_severity == 0) {
			int back_time = med_time_converter.add_subtract_time(ev.time, signals_time_unit, -730, MedTime::Days);
			int found = 0;
			for (int k = j - 1; k >= 0; k--) {
				if (evs[k].time < back_time) break;
				if (evs[k].event_severity >= 1) {
					found = 1;
					break;
				}
			}
			if (found) continue;

			// we are for certain in a point of health (current severity 0 and no severity 1 in the last 2 years. We can not get here with severity 2 and above !
			if (ranges[0].first < 0) ranges[0].first = ev.time;
			ranges[0].second = ev.time;
		}
	}

}
