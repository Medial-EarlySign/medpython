#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "RepProcess.h"
#include "RepCreateRegistry.h"
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
		else if (field == "verbose") verbose = stoi(entry.second);

		// Diabetes
		else if (field == "dm_drug_sig") dm_drug_sig = entry.second;
		else if (field == "dm_drug_sets") {
			//MLOG("drug_sets: %s\n", entry.second.c_str());
			boost::split(dm_drug_sets, entry.second, boost::is_any_of(","));
		}
		else if (field == "dm_diagnoses_sig") dm_diagnoses_sig = entry.second;
		else if (field == "dm_diagnoses_sets") boost::split(dm_diagnoses_sets, entry.second, boost::is_any_of(","));
		else if (field == "dm_glucose_sig") dm_glucose_sig = entry.second;
		else if (field == "dm_hba1c_sig") dm_hba1c_sig = entry.second;
		else if (field == "dm_diagnoses_severity") dm_diagnoses_severity = stoi(entry.second);
		else if (field == "dm_bio_mode") dm_bio_mode = stoi(entry.second);

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
		else if (field == "ht_systolic_first") ht_systolic_first = stoi(entry.second) > 0;

		// Prteinuria
		else if (field == "urine_tests_categories") boost::split(urine_tests_categories, entry.second, boost::is_any_of("/"));

		// ckd
		else if (field == "ckd_egfr_sig") ckd_egfr_sig = entry.second;
		else if (field == "ckd_proteinuria_sig") ckd_proteinuria_sig = entry.second;

		//custom
		else if (field == "registry_custom_type") registry_custom_type = entry.second;
		else if (field == "registry_custom_args") registry_custom_args = entry.second;

		else if (field == "rp_type") {}
		else MTHROW_AND_ERR("Error in RepCreateRegistry::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepCreateRegistry::init]
	}

	registry = name2type.at(registry_name);

	// Time unit
	if (time_unit == -1)
		time_unit = global_default_time_unit;

	// Input/Output signal names]
	virtual_signals.clear();
	virtual_signals_generic.clear();
	virtual_signals_generic = type2Virtuals.at(registry);
	if (!names.empty()) {
		if (names.size() != type2Virtuals.at(registry).size())
			MTHROW_AND_ERR("Wrong number of names supplied for RepCreateRegistry::%s - supplied %zd, required %zd\n", registry_name.c_str(), names.size(),
				type2Virtuals.at(registry).size());
		for (size_t i = 0; i < names.size(); i++)
			virtual_signals_generic[i].first = names[i];
	}

	MLOG_D("virtual 0 : %s\n", virtual_signals_generic[0].first.c_str());

	if (registry == REP_REGISTRY_CUSTOM) {
		//update signals myself:
		//only supports keep alive
		custom_registry = MedRegistry::make_registry(registry_custom_type, registry_custom_args);
		custom_registry->get_registry_creation_codes(signals);

		registry_values.clear(); //no output categories
	}
	else {
		if (!signals.empty()) {
			if (signals.size() != type2reqSigs.at(registry).size())
				MTHROW_AND_ERR("Wrong number of signals supplied for RepCreateRegistry::%s - supplied %zd, required %zd\n", registry_name.c_str(), signals.size(),
					type2reqSigs.at(registry).size());
		}
		else {
			if (type2reqSigs.find(registry) == type2reqSigs.end())
				MTHROW_AND_ERR("Error - RepCreateRegistry::init - must provide signals for this registry type\n");
			signals = type2reqSigs.at(registry);
		}
	}

	// required/affected signals
	init_lists();

	// Default initializations
	if (registry == REP_REGISTRY_HT)
		ht_init_defaults();

	if (registry == REP_REGISTRY_DM) {
		if (registry_values.empty())
			registry_values = dm_reg_values;
	}

	if (registry == REP_REGISTRY_PROTEINURIA) {
		if (registry_values.empty())
			registry_values = proteinuria_reg_values;
	}

	if (registry == REP_REGISTRY_CKD) {
		if (registry_values.empty())
			registry_values = ckd_reg_values;
	}

	return 0;
}

/// Required/Affected signals
//=======================================================================================
void RepCreateRegistry::init_lists() {

	req_signals.clear();
	for (string signalName : signals)
		req_signals.insert(signalName);

	aff_signals.clear();
	for (auto& rec : virtual_signals_generic)
		aff_signals.insert(rec.first);
}

// making sure virtual_signals_generic exists
//=======================================================================================
void RepCreateRegistry::post_deserialization() {
	if (virtual_signals_generic.empty())
		MTHROW_AND_ERR("Error virtual_signals_generic not serialized - please call ModelConvertor to update model.  /server/Work/FrozenTools/ModelConvertor/ModelConvertor --input_model <input_model>  --output_model <output_model> \n");
}


// making sure V_ids and sigs_ids are initialized
//=======================================================================================
void RepCreateRegistry::init_tables(MedDictionarySections& dict, MedSignals& sigs) {

	virtual_ids.clear();
	for (auto& rec : virtual_signals_generic)
		virtual_ids.push_back(sigs.sid(rec.first));

	req_signal_ids.clear();
	sig_ids_s.clear();
	sig_ids.clear();
	for (auto &rsig : signals) {
		int sid = sigs.sid(rsig);
		if (sid < 0)
			MTHROW_AND_ERR("Error in RepCreateRegistry::init_tables - unrecognized signal %s\n", rsig.c_str());
		sig_ids_s.insert(sid);
		sig_ids.push_back(sid);
		req_signal_ids.insert(sid);
	}

	aff_signal_ids.clear();
	aff_signal_ids.insert(virtual_ids.begin(), virtual_ids.end());

	// Dictionary
	if (!registry_values.empty()) {
		dict.add_section(virtual_signals_generic[0].first);
		int newSectionId = dict.section_id(virtual_signals_generic[0].first);
		for (size_t i = 1; i < virtual_signals_generic.size(); i++)
			dict.connect_to_section(virtual_signals_generic[i].first, newSectionId);

		for (size_t i = 0; i < registry_values.size(); i++)
			dict.dicts[newSectionId].push_new_def(registry_values[i], (int)i);
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
	else if (registry == REP_REGISTRY_PROTEINURIA)
		init_proteinuria_registry_tables(dict, sigs);
	else if (registry == REP_REGISTRY_CKD)
		init_ckd_registry_tables(dict, sigs);
	else if (registry == REP_REGISTRY_CUSTOM) {
		bDateCode = sigs.sid("BDATE");
		vector<string> create_custom_input_sigs_names;
		custom_registry->get_registry_use_codes(create_custom_input_sigs_names);
		create_custom_input_sigs.resize(create_custom_input_sigs_names.size());
		for (size_t i = 0; i < create_custom_input_sigs_names.size(); ++i)
		{
			create_custom_input_sigs[i] = sigs.sid(create_custom_input_sigs_names[i]);
			if (create_custom_input_sigs[i] < 0)
				MTHROW_AND_ERR("Error in RepCreateRegistry::init_tables - signal %s was not found\n",
					create_custom_input_sigs_names[i].c_str());
		}
	}
}

void RepCreateRegistry::get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const {
	if (registry == REP_REGISTRY_HT) {
		string signal_name_diag = signals[rc_idx];
		unordered_set<string> uniq_set_codes(ht_identifiers.begin(), ht_identifiers.end());
		uniq_set_codes.insert(chf_identifiers.begin(), chf_identifiers.end());
		uniq_set_codes.insert(mi_identifiers.begin(), mi_identifiers.end());
		uniq_set_codes.insert(af_identifiers.begin(), af_identifiers.end());
		vector<string> uniq_set_ls(uniq_set_codes.begin(), uniq_set_codes.end());
		signal_categories_in_use[signal_name_diag] = move(uniq_set_ls);

		string signal_name_drugs = signals[drug_idx];
		unordered_set<string> uniq_set_drugs(ht_drugs.begin(), ht_drugs.end());
		uniq_set_drugs.insert(ht_chf_drugs.begin(), ht_chf_drugs.end());
		uniq_set_drugs.insert(ht_dm_drugs.begin(), ht_dm_drugs.end());
		uniq_set_drugs.insert(ht_extra_drugs.begin(), ht_extra_drugs.end());
		vector<string> uniq_set_ls_drugs(uniq_set_drugs.begin(), uniq_set_drugs.end());
		signal_categories_in_use[signal_name_drugs] = move(uniq_set_ls_drugs);
	}
	else if (registry == REP_REGISTRY_DM) {
		if (dm_drug_sig != "")
			signal_categories_in_use[dm_drug_sig] = dm_drug_sets;
		if (dm_diagnoses_sig != "")
			signal_categories_in_use[dm_diagnoses_sig] = dm_diagnoses_sets;
	}
	else if (registry == REP_REGISTRY_PROTEINURIA) {
		for (auto &c : urine_tests_categories) {
			vector<string> f;
			boost::split(f, c, boost::is_any_of(":"));
			RegistryDecisionRanges rdr;
			rdr.sig_name = f[0];
			rdr.is_numeric = stoi(f[1]);
			for (int j = 2; j < f.size(); j++) {
				vector<string> f2;
				boost::split(f2, f[j], boost::is_any_of(","));
				if (rdr.is_numeric) {
					rdr.ranges.push_back(pair<float, float>(stof(f2[0]), stof(f2[1])));
				}
				else {
					rdr.categories.push_back(f2);
				}
			}
			if (!rdr.is_numeric) {
				unordered_set<string> flat_set;
				for (vector<string> &e : rdr.categories)
					flat_set.insert(e.begin(), e.end());
				vector<string> flat_ls(flat_set.begin(), flat_set.end());
				signal_categories_in_use[rdr.sig_name] = move(flat_ls);
			}
		}
	}
	else if (registry == REP_REGISTRY_CKD) {
		//empty - no categories
	}
	else if (registry == REP_REGISTRY_CUSTOM) {
		//empty - no categories
	}
	else
		MTHROW_AND_ERR("Unsupported registry type %d\n", registry);
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

	allVersionsIterator vit(rec, sig_ids_s);


	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		for (size_t isig = 0; isig < sig_ids.size(); isig++)
			rec.uget(sig_ids[isig], iver, usvs[isig]);

		if (registry == REP_REGISTRY_HT)
			ht_registry_apply(rec, time_points, iver, usvs, all_v_vals, all_v_times, final_sizes);
		else if (registry == REP_REGISTRY_DM)
			dm_registry_apply(rec, time_points, iver, usvs, all_v_vals, all_v_times, final_sizes);
		else if (registry == REP_REGISTRY_PROTEINURIA)
			proteinuria_registry_apply(rec, time_points, iver, usvs, all_v_vals, all_v_times, final_sizes);
		else if (registry == REP_REGISTRY_CKD)
			ckd_registry_apply(rec, time_points, iver, usvs, all_v_vals, all_v_times, final_sizes);
		else if (registry == REP_REGISTRY_CUSTOM)
			custom_registry_apply(rec, time_points, iver, usvs, all_v_vals, all_v_times, final_sizes);

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

	for (int j = 0; j < sig_ids.size(); j++) {
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
	int sys_ch = 1;
	if (ht_systolic_first)
		sys_ch = 0;

	int bdate = usvs[bdate_idx].Val(0);
	int byear = bdate;
	if (signals[bdate_idx] == "BDATE")
		byear = int(bdate / 10000);
	vector<pair<int, int> > data; // 0 = Normal BP ; 1 = High BP ; 20 + X = HT Drug ; 3 = HT Read Code (4/5/6 = CHF/MI/AF Read Codes ; 7 = DM)

	// Blood Pressure
	int _num_ht_measurements = 0, _found_op_ht = 0;
	for (int i = 0; i < usvs[bp_idx].len; i++) {
		int time = usvs[bp_idx].Time(i);
		if (time_points.size() != 0 && time > time_points[iver])
			break;
		if (usvs[bp_idx].Val(i, sys_ch) < usvs[bp_idx].Val(i, 1 - sys_ch))
			++_found_op_ht;
		++_num_ht_measurements;
		int age = 1900 + med_time_converter.convert_times(signal_time_units[bdate_idx], MedTime::Years, time) - byear;
		int bpFlag = ((age >= 60 && usvs[bp_idx].Val(i, sys_ch) > 150) || (age < 60 && usvs[bp_idx].Val(i, sys_ch)  > 140) || usvs[bp_idx].Val(i, 1 - sys_ch) > 90) ? 1 : 0;
		data.push_back({ med_time_converter.convert_times(signal_time_units[bdate_idx], MedTime::Days, time) , bpFlag });
		// cout << "BP" << i << " time " << time << " " << med_time_converter.convert_times(signal_time_units[bdate_idx], MedTime::Days, time) << " " << usvs[bp_idx].Val(i, 0) << " " << usvs[bp_idx].Val(i, 1) << " " << bpFlag << "\n";
	}
#pragma omp critical 
	{
	num_ht_measurements += _num_ht_measurements;
	found_op_ht += _found_op_ht;
	}

	// Drugs
	for (int i = 0; i < usvs[drug_idx].len; i++) {
		int time = usvs[drug_idx].Time(i);
		if (time_points.size() != 0 && time > time_points[iver])
			break;

		if (htDrugLut[(int)usvs[drug_idx].Val(i, 0)])
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

		// Ignore illegal date
		if (days < 0)
			continue;

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
					if (lastBP == 1)
						bpStatus = 2;
					lastBP = 1;
				}
				else if (info > 20) { // HT Drug, move to unclear (depending on background)
									  // using drug-specific info : 21 = inidcative of HT. 22 = indicative of HT unless CHF. 23 = indicative of HT unless diabetes. 
									  // 24 = indicative of HT unless CHF/MI/AF.
					if (info == 21 || (info == 22 && !chfStatus) || (info == 23 && !dmStatus) || (info == 24 && !chfStatus && !miStatus && !afStatus)) {
						bpStatus = 1;
						lastDrugDays = days;
					}
				}
				else if (info == 3) { // Read Code. if last bp is normal, mark as unclear, otherwise, mark as HT
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
		if (verbose)
			MLOG("id %d ver %d (%d) . Date %d . Info = %d => %d\n", rec.pid, iver, time_points[iver], med_time_converter.convert_days(MedTime::Date, days), info, bpStatusToPush);
	}

	// Collect
	int firstNorm = -1, lastNorm = -1, firstHT = -1, lastHT = -1;

	for (unsigned int i = 0; i < bpStatusVec.size(); i++) {
		if (bpStatusVec[i] == 2) {
			lastHT = data[i].first;
			if (firstHT == -1)
				firstHT = data[i].first;

			// If end-time is given, we can stop now.
			if (time_points.size() != 0) {
				lastHT = med_time_converter.convert_times(global_default_time_unit, MedTime::Days, time_points[iver]);
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
		all_v_times[0].push_back(med_time_converter.convert_times(MedTime::Days, time_unit, firstNorm));
		all_v_times[0].push_back(med_time_converter.convert_times(MedTime::Days, time_unit, lastNorm));
	}

	if (lastHT > 0) {
		final_sizes[0]++;
		all_v_vals[0].push_back(1);
		all_v_times[0].push_back(med_time_converter.convert_times(MedTime::Days, time_unit, firstHT));
		all_v_times[0].push_back(med_time_converter.convert_times(MedTime::Days, time_unit, lastHT));
	}
}

void RepCreateRegistry::custom_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, const vector<UniversalSigVec>& usvs,
	vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes) {
	vector<MedRegistryRecord> vals;
	custom_registry->get_pid_records(rec, bDateCode, create_custom_input_sigs, vals);

	//update outputs: all_v_vals, all_v_times, final_sizes: - 1 output signal:
	final_sizes = { (int)vals.size() };
	all_v_vals.resize(1);
	all_v_times.resize(1);
	vector<float> &out_vals = all_v_vals[0];
	vector<int> &out_times = all_v_times[0];
	//update out_vals, out_times:
	out_vals.resize(vals.size());
	out_times.resize(2 * vals.size()); //2 time channels
	for (size_t i = 0; i < vals.size(); ++i)
	{
		out_vals[i] = vals[i].registry_value;
		out_times[2 * i + 0] = vals[i].start_date;
		out_times[2 * i + 1] = vals[i].end_date;
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
		read_identifiers_list(pPath, "hyper_tension.desc", ht_identifiers);

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
void RepCreateRegistry::init_dm_registry_tables(MedDictionarySections& dict, MedSignals& sigs)
{

	int i = 0;
	for (auto &rsig : signals) {
		if (rsig == dm_drug_sig) dm_drug_idx = i;
		if (rsig == dm_diagnoses_sig) dm_diagnoses_idx = i;
		if (rsig == dm_glucose_sig) dm_glucose_idx = i;
		if (rsig == dm_hba1c_sig) dm_hba1c_idx = i;
		//MLOG("dm_reg: indexes : rsig %s : i %d : drug %d diag %d glu %d hb %d\n", rsig.c_str(), i, dm_drug_idx, dm_diagnoses_idx, dm_glucose_idx, dm_hba1c_idx);
		i++;
	}

	// lookup tables
	if (dm_drug_sig != "")	dict.dicts[dict.section_id(dm_drug_sig)].prep_sets_lookup_table(dm_drug_sets, dm_drug_lut);
	if (dm_diagnoses_sig != "") dict.dicts[dict.section_id(dm_diagnoses_sig)].prep_sets_lookup_table(dm_diagnoses_sets, dm_diagnoses_lut);

}

//==================================================================================================================================================
void RepCreateRegistry::dm_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes)
{
	vector<RegistryEvent> evs;

	int time = -1;
	if (time_points.size() > 0) time = time_points[iver];
	int time_unit = signal_time_units[0]; // taking the assumption that all our signals use THE SAME time unit

	// step 1 collect events

	// glucose events
	if (dm_glucose_idx >= 0) {

		UniversalSigVec &glu_usv = usvs[dm_glucose_idx];
		//MLOG("Glucose for pid %d -> %d len\n", rec.pid, glu_usv.len);
		for (int i = 0; i < glu_usv.len; i++) {
			//MLOG("Glucose : pid %d %d : i %d : %d , %f\n", rec.pid, time, i, glu_usv.Time(i), glu_usv.Val(i));
			int i_time = glu_usv.Time(i);
			if (time > 0 && time < i_time) break;
			float i_val = glu_usv.Val(i);
			int severity = 0;
			if (i_val > 99.99f) severity = 1;	// 2 tests enough to define pre diabetic
			if (i_val > 109.99f) severity = 2;	// single test enough to define pre diabetic
			if (i_val > 125.99f) severity = 3;	// 2 tests enough to define diabetic
			if (i_val > 199.99f) severity = 4;  // 1 test enough to define diabetic

			evs.push_back(RegistryEvent(i_time, REG_EVENT_DM_GLUCOSE, i_val, severity));
			// cout << "Glucose " << i_time << " " << i_val << " " << severity << "\n";
		}
	}

	// HbA1C events
	if (dm_hba1c_idx >= 0) {
		UniversalSigVec &hba1c_usv = usvs[dm_hba1c_idx];

		for (int i = 0; i < hba1c_usv.len; i++) {
			int i_time = hba1c_usv.Time(i);
			if (time > 0 && time < i_time) break;
			float i_val = hba1c_usv.Val(i);
			int severity = 0;
			//if (i_val > 99.99f) severity = 1;	// 2 tests enough to define pre diabetic
			if (i_val > 5.69f) severity = 2;	// single test enough to define pre diabetic
			if (i_val > 6.49f) severity = 3;	// 2 tests enough to define diabetic
			if (i_val > 7.99f) severity = 4;  // 1 test enough to define diabetic

			evs.push_back(RegistryEvent(i_time, REG_EVENT_DM_HBA1C, i_val, severity));
			// cout << "hba1c " << i_time << " " << i_val << " " << severity << "\n";
		}
	}

	// Drug Events
	if (dm_drug_idx >= 0) {
		UniversalSigVec &drug_usv = usvs[dm_drug_idx];

		for (int i = 0; i < drug_usv.len; i++) {
			int i_time = drug_usv.Time(i);
			if (time > 0 && time < i_time) break;
			int i_val = (int)drug_usv.Val(i);
			if (i_val < 0 || i_val > dm_drug_lut.size())
				MTHROW_AND_ERR("ERROR in dm Registry drug_idx : got i_val %d while lut size is %d\n", i_val, (int)dm_drug_lut.size());
			if (dm_drug_lut.size() > 0 && dm_drug_lut[i_val]) {
				int severity = 4; // currently the first diabetic drug usage makes you diabetic for life.... this is extreme, but given this, we only need the first.
				evs.push_back(RegistryEvent(i_time, REG_EVENT_DM_DRUG, 1, severity));
				// cout << "Drug " << i_time << " " << severity << "\n";
				break;
			}

		}
	}


	// Diagnoses Events
	if (dm_diagnoses_idx >= 0) {
		UniversalSigVec &diag_usv = usvs[dm_diagnoses_idx];

		for (int i = 0; i < diag_usv.len; i++) {
			int i_time = diag_usv.Time(i);
			if (time > 0 && time < i_time) break;
			int i_val = (int)diag_usv.Val(i);
			if (i_val < 0 || i_val > dm_diagnoses_lut.size())
				MTHROW_AND_ERR("ERROR in dm Registry diagnoses_idx : got i_val %d while lut size is %d\n", i_val, (int)dm_diagnoses_lut.size());
			if (dm_diagnoses_lut.size() > 0 && dm_diagnoses_lut[i_val]) {
				int severity = dm_diagnoses_severity;
				evs.push_back(RegistryEvent(i_time, REG_EVENT_DM_DIAGNOSES, 1, severity));
				// cout << "Diagnosis " << i_time << " " << severity << "\n";
				if (dm_diagnoses_severity >= 4) break;
			}

		}
	}

	// collection of events done

	// sorting events
	sort(evs.begin(), evs.end(), [](const RegistryEvent &v1, const RegistryEvent &v2) { return v1.time < v2.time; });

	// applying rules
	vector<pair<int, int>> ranges(3, pair<int, int>(-1, -1)); // 0: for healthy, 1: for prediabetic , 2: for diabetic

	for (int j = 0; j < evs.size(); j++) {

		auto &ev = evs[j];

		//MLOG("diabetes reg : j %d ev: time %d type %d val %f severity %d\n", j, ev.time, ev.event_type, ev.event_val, ev.event_severity);

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
			int back_time = med_time_converter.add_subtract_time(ev.time, time_unit, -730, MedTime::Days);

			int found = 0;
			int first_index = 0;
			for (int k = j - 1; k >= 0; k--) {
				if (evs[k].time < back_time) break;
				if (evs[k].event_severity == 3) {
					found = 1;
					first_index = k;
					break;
				}
			}
			if (found) {
				// found a diabetic, several cases now :
				// (1) dm_bio_mode = 1 : we take the time of the first indication
				// (2) the type of first is REG_EVENT_DM_DIAGNOSES : we take the time of the first 
				// (3) other cases : we take the second
				if (dm_bio_mode || evs[first_index].event_type == REG_EVENT_DM_DIAGNOSES) {
					ranges[2].first = evs[first_index].time;
					ranges[2].second = time;
				}
				else {
					ranges[2].first = ev.time;
					ranges[2].second = time;
				}
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
			int back_time = med_time_converter.add_subtract_time(ev.time, time_unit, -730, MedTime::Days);

			int found = 0;
			int first_index = 0;
			for (int k = j - 1; k >= 0; k--) {
				if (evs[k].time < back_time) break;
				if (evs[k].event_severity == 1) {
					found = 1;
					first_index = k;
					break;
				}
			}
			if (found) {
				if (dm_bio_mode) {
					ranges[1].first = evs[first_index].time;
					ranges[1].second = ev.time;
				}
				else {
					ranges[1].first = ev.time;
					ranges[1].second = ev.time;
				}
				continue;
			}
		}

		if (ev.event_severity == 0) {
			int back_time = med_time_converter.add_subtract_time(ev.time, time_unit, -730, MedTime::Days);
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

	// now preparing for this line :
	// 			rec.set_version_universal_data(virtual_ids[ivir], iver, &(all_v_times[ivir][0]), &(all_v_vals[ivir][0]), final_sizes[ivir]);
	// first dimension is initialized already

	all_v_times[0].clear();
	all_v_vals[0].clear();
	final_sizes[0] = 0;

	for (int j = 0; j < 3; j++)
		if (ranges[j].first > 0) {
			// push Healthy, Pre, or DM
			all_v_vals[0].push_back((float)j);
			all_v_times[0].push_back(ranges[j].first);
			all_v_times[0].push_back(ranges[j].second);
			final_sizes[0]++;
		}
#if 0
	// debug print
	int c = 0;
	for (auto &ev : evs) {
		MLOG("pid %d %d : ev %d : time %d type %d val %f severity %d\n", rec.pid, time, c++, ev.time, ev.event_type, ev.event_val, ev.event_severity);
	}
	MLOG("DM_registry calculation: pid %d %d : Healthy %d %d : Pre %d %d : Diabetic %d %d\n", rec.pid, time, ranges[0].first, ranges[0].second, ranges[1].first, ranges[1].second, ranges[2].first, ranges[2].second);
#endif
}

//===============================================================================================================================
// Proteinuria (3 levels) code
//===============================================================================================================================
// proteinuria tables
void RepCreateRegistry::init_proteinuria_registry_tables(MedDictionarySections& dict, MedSignals& sigs)
{
	proteinuria_ranges.clear();
	for (auto &c : urine_tests_categories) {
		//MLOG("Parsing %s\n", c.c_str());
		vector<string> f;
		boost::split(f, c, boost::is_any_of(":"));
		RegistryDecisionRanges rdr;
		rdr.sig_name = f[0];
		rdr.is_numeric = stoi(f[1]);
		for (int j = 2; j < f.size(); j++) {
			vector<string> f2;
			boost::split(f2, f[j], boost::is_any_of(","));
			if (rdr.is_numeric) {
				rdr.ranges.push_back(pair<float, float>(stof(f2[0]), stof(f2[1])));
			}
			else {
				rdr.categories.push_back(f2);
			}
		}
		rdr.usv_idx = -1;
		//MLOG("rdr %s is_n %d ", rdr.sig_name.c_str(), rdr.is_numeric);
		//for (auto &e : rdr.ranges) MLOG(" %f-%f ", e.first, e.second);
		//for (auto &e : rdr.categories) for (auto &s : e) MLOG(" %s ", s.c_str());
		//MLOG("\n");
		proteinuria_ranges.push_back(rdr);

	}

	for (auto &r : proteinuria_ranges) {
		r.sig_id = sigs.sid(r.sig_name);
		if (!r.is_numeric) {
			r.categories_i.resize(r.categories.size());
			for (int j = 0; j < r.categories.size(); j++) {
				r.categories_i[j].clear();
				int section_id = dict.section_id(r.sig_name);
				unordered_set<int> vals_for_j;
				for (auto &c : r.categories[j]) {
					if (dict.dicts[section_id].Name2Id.find(c) != dict.dicts[section_id].Name2Id.end()) {
						vector<int> m;
						dict.dicts[section_id].get_set_members(c, m);
						vals_for_j.insert(m.begin(), m.end());
					}
				}
				r.categories_i[j].insert(r.categories_i[j].begin(), vals_for_j.begin(), vals_for_j.end());
			}
		}

		int j = 0;
		for (auto &rsig : signals) {
			if (r.sig_name == rsig) {
				r.usv_idx = j;
				break;
			}
			j++;
		}
	}

}

// proteinuria apply
void RepCreateRegistry::proteinuria_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes)
{
	int time = -1;
	if (time_points.size() > 0) time = time_points[iver];

	vector<pair<int, int>> proteinuria_ev;

	// collecting events
	for (auto &r : proteinuria_ranges) {
		UniversalSigVec &rusv = usvs[r.usv_idx];
		//MLOG("Proteinuria: pid %d,%d : sig %s , %d : idx %d : len %d\n", rec.pid, time, r.sig_name.c_str(), r.sig_id, r.usv_idx, usvs[r.usv_idx].len);
		for (int j = 0; j < rusv.len; j++) {

			int i_time = rusv.Time(j);
			if (time > 0 && i_time > time) break;

			int found = -1;
			if (r.is_numeric) {
				float f_val = rusv.Val(j);
				//MLOG("pid %d,%d : sig %s : time %d val %f\n", rec.pid, time, r.sig_name.c_str(), i_time, f_val);
				for (int k = 0; found < 0 && k < r.ranges.size(); k++) {
					if (f_val >= r.ranges[k].first && f_val < r.ranges[k].second)
						found = k;
				}
			}
			else {
				int i_val = (int)rusv.Val(j);
				//MLOG("pid %d,%d : sig %s : time %d val %d\n", rec.pid, time, r.sig_name.c_str(), i_time, i_val);
				for (int k = 0; found < 0 && k < r.categories_i.size(); k++) {
					for (auto &ci : r.categories_i[k]) {
						//MLOG("testing ci %d i_val %d\n", ci, i_val);
						if (i_val == ci) {
							found = k;
							break;
						}
					}
				}

			}
			if (found >= 0) {
				//MLOG("pid %d,%d : sig %s : =====> time %d found %d\n", rec.pid, time, r.sig_name.c_str(), i_time, found);
				proteinuria_ev.push_back(pair<int, int>(i_time, found));
			}
		}
	}

	// get the max value for each day, sort and unique all in one using map
	map<int, int> time2val;
	for (auto &p : proteinuria_ev) {
		if (time2val.find(p.first) == time2val.end())
			time2val[p.first] = p.second;
		else
			if (p.second > time2val[p.first])
				time2val[p.first] = p.second;
	}

	// loading into all_v_times, all_v_vals, final_sizes
	all_v_times[0].clear();
	all_v_vals[0].clear();
	final_sizes[0] = 0;

	for (auto &e : time2val) {
		// push Healthy, Pre, or DM
		all_v_vals[0].push_back((float)e.second);
		all_v_times[0].push_back(e.first);
		final_sizes[0]++;
	}

#if 0
	// debug
	MLOG("Proteinuria State : pid %d %d : ", rec.pid, time);
	for (auto &e : time2val)
		MLOG(" %d,%d :", e.first, e.second);
	MLOG("\n");
#endif
}

//===============================================================================================================================
// CKD (5 levels) code
//===============================================================================================================================
// ckd tables
void RepCreateRegistry::init_ckd_registry_tables(MedDictionarySections& dict, MedSignals& sigs)
{
	int i = 0;
	for (auto &rsig : signals) {
		if (rsig == ckd_egfr_sig) ckd_egfr_idx = i;
		if (rsig == ckd_proteinuria_sig) ckd_proteinuria_idx = i;
		//MLOG("CKD: rsig %s ( %s , %s ) i %d idx: e %d p %d\n", rsig.c_str(), ckd_egfr_sig.c_str(), ckd_proteinuria_sig.c_str(), i, ckd_egfr_idx, ckd_proteinuria_idx);
		i++;
	}
}

// ckd apply
void RepCreateRegistry::ckd_registry_apply(PidDynamicRec& rec, vector<int>& time_points, int iver, vector<UniversalSigVec>& usvs, vector<vector<float>>& all_v_vals, vector<vector<int>>& all_v_times, vector<int>& final_sizes)
{
	int time = -1;
	if (time_points.size() > 0) time = time_points[iver];

	map<int, int> ckd_ev;

	UniversalSigVec &p_usv = usvs[ckd_proteinuria_idx];
	UniversalSigVec &e_usv = usvs[ckd_egfr_idx];

	int ip = 0, ie = 0;
	float last_p = -1, last_e = -1;
	int e_time = -1, p_time = -1;
	float e_val = -1, p_val = -1;
	int i_time = -1;
	while (ip < p_usv.len || ie < e_usv.len) {

		if (ie < e_usv.len) {
			e_time = e_usv.Time(ie);
			e_val = e_usv.Val(ie);
		}
		else
			e_time = -1;

		if (ip < p_usv.len) {
			p_time = p_usv.Time(ip);
			p_val = p_usv.Val(ip);
		}
		else
			p_time = -1;

		i_time = -1;

		if (e_time > 0 && (p_time < 0 || e_time <= p_time)) {
			i_time = e_time;
			last_e = e_val;
			ie++;
		}

		if (p_time > 0 && (e_time < 0 || p_time <= e_time)) {
			i_time = p_time;
			last_p = p_val;
			ip++;
		}

		if (time > 0 && i_time > time) break;

		//MLOG("CKD: pid %d %d : ip %d/%d ie %d/%d i_time %d last_e %f last_p %f\n", rec.pid, time, ip, p_usv.len, ie, e_usv.len, i_time, last_e, last_p);

		// we now have last_e, last_p and i_time , and can insert a new event
		if (i_time > 0) {
			pair<int, int> ev(i_time, -1);

			if (last_e <= 15 && last_e >= 0) ev.second = 4;

			else if ((last_e < 0 && last_p <= 0) || (last_e > 60 && last_p <= 0))
				ev.second = 0;

			else if ((last_e > 45 && last_e <= 60 && last_p <= 0) ||
				(last_e > 60 && last_p == 1) ||
				(last_e < 0 && last_p == 1))
				ev.second = 1;

			else if ((last_e > 30 && last_e <= 45 && last_p <= 0) ||
				(last_e > 45 && last_e <= 60 && last_p == 1) ||
				(last_e > 60 && last_p == 2) ||
				(last_e < 0 && last_p == 2))
				ev.second = 2;

			else if ((last_e > 15 && last_e <= 30) ||
				(last_e > 30 && last_e <= 45 && last_p >= 1) ||
				(last_e > 45 && last_e <= 60 && last_p == 2))
				ev.second = 3;

			if (ev.second >= 0) {
				if (ckd_ev.find(i_time) == ckd_ev.end())
					ckd_ev[i_time] = ev.second;
				else
					if (ev.second > ckd_ev[i_time])
						ckd_ev[i_time] = ev.second;
			}
		}
	}

	// loading into all_v_times, all_v_vals, final_sizes
	all_v_times[0].clear();
	all_v_vals[0].clear();
	final_sizes[0] = 0;

	for (auto &e : ckd_ev) {
		all_v_vals[0].push_back((float)e.second);
		all_v_times[0].push_back(e.first);
		final_sizes[0]++;
	}

#if 0
	// debug
	MLOG("CKD State : pid %d %d : ", rec.pid, time);
	for (auto &e : ckd_ev)
		MLOG(" %d,%d :", e.first, e.second);
	MLOG("\n");
#endif


}

void RepCreateRegistry::make_summary() {
	if (num_ht_measurements > 0 && float(found_op_ht) / num_ht_measurements > 0.1)
		MTHROW_AND_ERR("Error in RepCreateRegistry - hypertenstion seems to be opposite (%2.2f) percentage\n",
			100 * float(found_op_ht) / num_ht_measurements);
}