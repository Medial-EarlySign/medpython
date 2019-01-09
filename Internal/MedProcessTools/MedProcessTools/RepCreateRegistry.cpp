#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "RepProcess.h"
#include <MedUtils/MedUtils/MedUtils.h>
#include <queue>

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
		else if (field == "time_unit") time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "registry_values") boost::split(registry_values, entry.second, boost::is_any_of(","));

		// Hypertension
		else if (field == "ht_identifiers") { boost::split(ht_identifiers, entry.second, boost::is_any_of(",")); ht_identifiers_given = true; }
		else if (field == "chf_identifiers") { boost::split(chf_identifiers, entry.second, boost::is_any_of(",")); chf_identifiers_given = true; }
		else if (field == "mi_identifiers") { boost::split(mi_identifiers, entry.second, boost::is_any_of(",")); mi_identifiers_given = true; }
		else if (field == "af_identifiers") { boost::split(af_identifiers, entry.second, boost::is_any_of(",")); af_identifiers_given = true; }
		else if (field == "ht_drugs") boost::split(ht_drugs, entry.second, boost::is_any_of(","));
		else if (field == "ht_chf_drugs") boost::split(ht_chf_drugs, entry.second, boost::is_any_of(","));
		else if (field == "ht_dm_drugs") boost::split(ht_dm_drugs, entry.second, boost::is_any_of(","));
		else if (field == "ht_extra_drugs") boost::split(ht_extra_drugs, entry.second, boost::is_any_of(","));
		else if (field == "ht_drugs_gap") ht_drugs_gap = stoi(entry.second);

		else if (field == "rp_type") {}
		else MTHROW_AND_ERR("Error in RepCreateRegistry::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepCreateRegistry::init]
	}

	registry = name2type.at(registry_name);

	// Time unit
	if (time_unit == -1)
		time_unit = global_default_time_unit;

	// Input/Output signal names
	virtual_signals = type2Virtuals.at(registry);
	if (!names.empty()) {
		if (names.size() != type2Virtuals.at(registry).size())
			MTHROW_AND_ERR("Wrong number of names supplied for RepCreateRegistry::%s - supplied %zd, required %zd\n", registry_name.c_str(), names.size(),
				type2Virtuals.at(registry).size());
		for (size_t i = 0; i < names.size(); i++)
			virtual_signals[i].first = names[i];
	}

	if (!signals.empty()) {
		if (signals.size() != type2reqSigs.size())
			MTHROW_AND_ERR("Wrong number of signals supplied for RepCreateRegistry::%s - supplied %zd, required %zd\n", registry_name.c_str(), signals.size(),
				type2reqSigs.at(registry).size());
	}
	else
		signals = type2reqSigs.at(registry);

	// required/affected signals
	init_lists();

	// Default initializations
	if (registry == REP_REGISTRY_HT)
		ht_init_defaults();

	return 0;
}

/// Required/Affected signals
//=======================================================================================
void RepCreateRegistry::init_lists() {

	req_signals.clear();
	for (string signalName : signals) 
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
		req_signal_ids.insert(sid);
	}

	aff_signal_ids.clear();
	aff_signal_ids.insert(virtual_ids.begin(), virtual_ids.end());

	// Dictionary
	if (!registry_values.empty()) {
		dict.add_section(virtual_signals[0].first);
		int newSectionId = dict.section_id(virtual_signals[0].first);
		for (size_t i = 1; i < virtual_signals.size(); i++)
			dict.connect_to_section(virtual_signals[i].first, newSectionId);

		for (size_t i = 0; i < registry_values.size(); i++)
			dict.dicts[newSectionId].push_new_def(registry_values[i], i);
	}

	// Time units
	signal_time_units.resize(sig_ids.size());
	for (int i = 0; i < sig_ids.size(); i++)
		signal_time_units[i] = sigs.Sid2Info[sig_ids[i]].time_unit;

	// Registry specific tables
	if (registry == REP_REGISTRY_HT)
		init_ht_registry_tables(dict, sigs);
	else if (registry == REP_REGISTRY_DM)
		init_dm_registry_tables(dict, sigs);
}

// Applying
/// <summary> apply processing on a single PidDynamicRec at a set of time-points : Should be implemented for all inheriting classes </summary>
int RepCreateRegistry::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {

	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	// Utility vector
	vector<UniversalSigVec> usvs(sig_ids.size());
	vector<vector<float>> all_v_vals(virtual_ids.size());
	vector<vector<int>> all_v_times(virtual_ids.size());
	vector<int> final_sizes(virtual_ids.size());

	int len;
	allVersionsIterator vit(rec, sig_ids_s);

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		for (size_t isig = 0; isig < sig_ids.size(); isig++)
			rec.uget(sig_ids[isig], iver, usvs[isig]);

		if (registry == REP_REGISTRY_HT)
			ht_registry_apply(rec, time_points, iver, usvs, all_v_vals, all_v_times, final_sizes);
		else if (registry == REP_REGISTRY_DM)
			dm_registry_apply(rec, time_points, iver, usvs, all_v_vals, all_v_times, final_sizes);
		

		// pushing virtual data into rec
		for (size_t ivir = 0; ivir < virtual_ids.size(); ivir++) 
			rec.set_version_universal_data(virtual_ids[ivir], iver, &(all_v_times[ivir][0]), &(all_v_vals[ivir][0]), final_sizes[ivir]);
	}

	return 0;

}

// Registry-Specific functions : HyperTension
//=======================================================================================
void RepCreateRegistry::init_ht_registry_tables(MedDictionarySections& dict, MedSignals& sigs) {

	// Look up tables for HT/CHF/MI/AF
	int sigId = sig_ids[rc_idx];
	int sectionId = dict.section_id(signals[rc_idx]);

	dict.prep_sets_lookup_table(sectionId, ht_identifiers, htLut);
	dict.prep_sets_lookup_table(sectionId, chf_identifiers, chfLut);
	dict.prep_sets_lookup_table(sectionId, mi_identifiers, miLut);
	dict.prep_sets_lookup_table(sectionId, af_identifiers, afLut);

	// build drug-specific look : 0 = not relevant to HT. 1 = inidcative of HT. 2 = indicative of HT unless CHF. 3 = indicative of HT unless diabetes. 4 = indicative of HT unless CHF/MI/AF.
	sectionId = dict.section_id(signals[drug_idx]);
	buildLookupTableForHTDrugs(dict.dicts[sectionId], htDrugLut);

}

// Build a look up table for HT drugs
void RepCreateRegistry::fillLookupTableForHTDrugs(MedDictionary& dict, vector<char>& lut, vector<string>& sets, char val) {

	// convert names to ids
	vector<int> sig_ids;
	for (auto &name : sets) {
		int myid = dict.id(name);
		if (myid > 0)
			sig_ids.push_back(myid);
		else
			fprintf(stderr, "prep_sets_lookup_table() : Found bad name %s :: not found in dictionary()\n", name.c_str());
	}

	for (int j = 0; j<sig_ids.size(); j++) {
		queue<int> q;
		q.push(sig_ids[j]);

		while (q.size() > 0) {
			int s = q.front();
			q.pop();
			lut[s] = val;
			for (auto elem : dict.Set2Members[s])
				if (lut[elem] == 0)
					q.push(elem);

		}

	}
	return;
}

void RepCreateRegistry::buildLookupTableForHTDrugs(MedDictionary& dict, vector<char>& lut) {

	int maxId = dict.Id2Name.rbegin()->first;
	lut.assign(maxId + 1, 0);

	fillLookupTableForHTDrugs(dict, lut, ht_drugs, (char)1);
	fillLookupTableForHTDrugs(dict, lut, ht_chf_drugs, (char)2);
	fillLookupTableForHTDrugs(dict, lut, ht_dm_drugs, (char)3);
	fillLookupTableForHTDrugs(dict, lut, ht_extra_drugs, (char)4);

}

void RepCreateRegistry::ht_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, 
	vector<int>& final_sizes)
{

	int byear = usvs[byear_idx].Val(0);
	vector<pair<int, int> > data; // 0 = Normal BP ; 1 = High BP ; 20 + X = HT Drug ; 3 = HT Read Code (4/5/6 = CHF/MI/AF Read Codes ; 7 = DM)

	// Blood Pressure
	for (int i = 0; i < usvs[bp_idx].len; i++) {
		int time = usvs[bp_idx].Time(i);
		if (time_points.size() != 0 && time > time_points[iver])
			break;

		int age = 1900 + med_time_converter.convert_times(signal_time_units[byear_idx], MedTime::Years, time) - byear;
		int bpFlag = ((age >= 60 && usvs[bp_idx].Val(i,1) > 150) || (age < 60 && usvs[bp_idx].Val(i, 1)  > 140) || usvs[bp_idx].Val(i, 0)  > 90) ? 1 : 0;
		data.push_back({ med_time_converter.convert_times(signal_time_units[byear_idx], MedTime::Days, time) , bpFlag });
	}

	// Drugs
	for (int i = 0; i < usvs[drug_idx].len; i++) {
		int time = usvs[drug_idx].Time(i);
		if (time_points.size() != 0 && time > time_points[iver])
			break;

		if (htDrugLut[(int)usvs[drug_idx].Val(i,0)])
			data.push_back({ med_time_converter.convert_times(signal_time_units[drug_idx], MedTime::Days, time) , 20 + htDrugLut[usvs[drug_idx].Val(i,0)] });
	}

	// Identifiers (ReadCodes)
	for (int i = 0; i < usvs[rc_idx].len; i++) {
		int time = usvs[rc_idx].Time(i);
		if (time_points.size() != 0 && time > time_points[iver])
			break;

		int days = med_time_converter.convert_times(signal_time_units[rc_idx], MedTime::Days, time);
		if (htLut[(int)usvs[rc_idx].Val(i, 0)])
			data.push_back({ days, 3 });

		if (chfLut[(int)usvs[rc_idx].Val(i, 0)])
			data.push_back({ days, 4 });


		if (miLut[(int)usvs[rc_idx].Val(i, 0)])
			data.push_back({ days, 5 });

		if (afLut[(int)usvs[rc_idx].Val(i, 0)])
			data.push_back({ days, 6 });
	}

	// Diabetes
	for (int i = 0; i < usvs[dm_registry_idx].len; i++) {
		int time = usvs[dm_registry_idx].Time(i, 0);
		if (time_points.size() != 0 && time > time_points[iver])
			break;
		if (usvs[dm_registry_idx].Val(i) == 2)
			data.push_back({ med_time_converter.convert_times(signal_time_units[dm_registry_idx], MedTime::Days, time), 7 });
	}

	// Sort and analyze
	stable_sort(data.begin(), data.end(), [](const pair<int, int> &v1, const pair<int, int> &v2) {return (v1.first < v2.first); });

	int bpStatus = -1;
	vector<int> bpStatusVec;
	int lastBP = -1;
	int lastDrugDays = -1;
	int chfStatus = 0, miStatus = 0, afStatus = 0, dmStatus = 0;

	for (auto& irec : data) {
		int days = irec.first;
		int info = irec.second;

		int bpStatusToPush = -1;
		
		// Background : CHF/MI/AF/DM
		if (info == 4)
			chfStatus = 1;
		else if (info == 5)
			miStatus = 1;
		else if (info == 6)
			afStatus = 1;
		else if (info == 7)
			dmStatus = 1;
		else { // HT Indications 
			if (bpStatus <= 0) { // Non HyperTensinve (or no info yet)
				if (info == 0) { // Normal BP , still non hypertensive
					lastBP = 0;
					bpStatus = 0;
				}
				else if (info == 1) { // High BP, move to HT given previous indication
					if (lastBP == 1 || (lastDrugDays != -1 && days - lastDrugDays < ht_drugs_gap))
						bpStatus = 2;
					lastBP = 1;
				}
				else if (info > 20) { // HT Drug, move to unclear (depending on background)
									  // build drug-specific look : 1 = inidcative of HT. 2 = indicative of HT unless CHF. 3 = indicative of HT unless diabetes. 
								      //4 = indicative of HT unless CHF/MI/AF.
					if (info == 21 || (info == 22 && !chfStatus) || (info == 23 && !dmStatus) || (info == 24 && !chfStatus && !miStatus && !afStatus)) {
						bpStatus = 1;
						lastDrugDays = days;
					}
				}
				else if (info == 3) { // Read Code. if last bp is normal, mark as unclear, otherwise, mark as HT
					if (lastBP == 0)
						bpStatus = 1;
					else
						bpStatus = 2;
				}
			}
			else if (bpStatus == 1) { // Unclear.
				if (info == 0)  // Normal BP, still unclear
					lastBP = 0;
				else if (info == 1) { // High BP. move to HT if previous BP was also high
					if (lastBP == 1)
						bpStatus = 2;
					lastBP = 1;
				}
				else if (info > 20) { // HT Drug. move to HT if last BP was high or HT Drug was taken within the last 6 months
					if (info == 21 || (info == 22 && !chfStatus) || (info == 23 && !dmStatus) || (info == 24 && !chfStatus && !miStatus && !afStatus)) {
						if (lastBP == 1 || (lastDrugDays != -1 && days - lastDrugDays < ht_drugs_gap && days - lastDrugDays > 0))
							bpStatus = 2;
						lastDrugDays = days;
					}
				}
				else if (info == 3) // ReadCode. Move to HT
					bpStatus = 2;
			}

			bpStatusToPush = bpStatus;
		}
		bpStatusVec.push_back(bpStatusToPush);
//		MLOG("id %d ver %d (%d) . Date %d . Info = %d => %d\n", rec.pid, iver, time_points[iver], med_time_converter.convert_days(MedTime::Date, days), info, bpStatusToPush);
	}

	// Collect
	int firstNorm = -1, lastNorm = -1, firstHT = -1, lastHT = -1;

	int status = 0;
	for (unsigned int i = 0; i < bpStatusVec.size(); i++) {
		if (bpStatusVec[i] == 2) {
			lastHT = data[i].first;
			if (firstHT == -1)
				firstHT = data[i].first;

			// If end-time is given, we can stop now.
			if (time_points.size() != 0) {
				lastHT = med_time_converter.convert_times(global_default_time_unit,MedTime::Days,time_points[iver]);
				break;
			}
		}
		else if (bpStatusVec[i] == 0) {
			lastNorm = data[i].first;
			if (firstNorm == -1)
				firstNorm = lastNorm;
		}
	}

	if (lastNorm > 0) {
		final_sizes[0]++;
		all_v_vals[0].push_back(0);
		all_v_times[0].push_back(med_time_converter.convert_times(MedTime::Days,time_unit,firstNorm));
		all_v_times[0].push_back(med_time_converter.convert_times(MedTime::Days, time_unit, lastNorm));
	}

	if (lastHT > 0) {
		final_sizes[0]++;
		all_v_vals[0].push_back(1);
		all_v_times[0].push_back(med_time_converter.convert_times(MedTime::Days, time_unit, firstHT));
		all_v_times[0].push_back(med_time_converter.convert_times(MedTime::Days, time_unit, lastHT));
	}

}

// Get lists of identifiers from default files, if not given
void read_identifiers_list(char *pPath, string fileName, vector<string>& list) {

	if (pPath == NULL)
		MTHROW_AND_ERR("Cannot find root path for reading default file \'%s\'\n", fileName.c_str());

	string fullPathFileName = string(pPath) + "/Tools/Registries/Lists/" + fileName;
	medial::io::read_codes_file(fullPathFileName, list);
}

void RepCreateRegistry::ht_init_defaults() {

	// Read files
	char* pPath;
	pPath = getenv("MR_ROOT");

	if (!ht_identifiers_given)
		read_identifiers_list(pPath,"hyper_tension.desc", ht_identifiers);

	if (!chf_identifiers_given)
		read_identifiers_list(pPath, "heart_failure_events.desc", chf_identifiers);

	if (!mi_identifiers_given)
		read_identifiers_list(pPath, "mi.desc", mi_identifiers);

	if (!af_identifiers_given)
		read_identifiers_list(pPath, "AtrialFibrilatioReadCodes.desc", af_identifiers);

	// Registry values
	if (registry_values.empty())
		registry_values = ht_def_values;

}

// Registry-Specific functions : Diabetes
//=======================================================================================
void RepCreateRegistry::init_dm_registry_tables(MedDictionarySections& dict, MedSignals& sigs) {

}

void RepCreateRegistry::dm_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times,
	vector<int>& final_sizes) {

}