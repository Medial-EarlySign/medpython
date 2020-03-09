#include "DiabetesFinderGenerator.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int DiabetesFinderGenerator::_resolve(vector<DiabetesEvent>& dm_events, int calc_time, json& json_out) {
	int ret = 0;

	unordered_set<string> past_evidence_text, recent_evidence_text;
	int ret_recent = 0;
	int ret_past = 0;

	//mark events by the 2 events rull
	int _latest_noteable_event_time = -1;
	for (auto& de : dm_events) {
		if (de.time > calc_time)
			continue;
		if ((de.de_type == DFG_DIABETES_EVENT_GLUCOSE && de.val >= dm_by_second_glucose)
			|| (de.de_type == DFG_DIABETES_EVENT_HBA1C && de.val >= dm_by_second_hba1c)) {
			if (_latest_noteable_event_time != -1 && de.time - _latest_noteable_event_time < dm_by_second_time_delta) {
				de.is_second = true;
			}
			_latest_noteable_event_time = de.time;
		}
	}

	for (const auto& de : dm_events) {
		//MLOG("DE: (%d,%d,%f)\n", de.de_type, de.time, de.val);
		//MLOG("DE:    calc_time = %d\n", calc_time);
		if (de.time > calc_time)
			continue;
		bool is_recent = de.time >= med_time_converter.add_subtruct_days(calc_time, -1* dm_past_event_days);
		//MLOG("DE:    add_subtruct_days(calc_time, -1 * dm_past_event_days) = %d\n", med_time_converter.add_subtruct_days(calc_time, -1 * dm_past_event_days));
		//MLOG("DE:    is_recent = %d\n", (int)is_recent);
		//MLOG("DE:    dm_by_single_glucose = %d\n", (int)dm_by_single_glucose);
		//MLOG("DE:    (de.de_type == DFG_DIABETES_EVENT_GLUCOSE) = %d\n", (int)(de.de_type == DFG_DIABETES_EVENT_GLUCOSE));
		string reason_str = "";
		if (de.de_type == DFG_DIABETES_EVENT_GLUCOSE && de.val >= dm_by_single_glucose) {
			reason_str = string("Single measurement of glucose >=") + to_string((int)dm_by_single_glucose) + " mg/dL";
			if (is_recent)
				ret_recent |= REASON_RECENT_LABS;
			else ret_past |= REASON_PAST_LABS;
		}
		if (de.is_second) {
			reason_str = string("Second measurement of glucose >= ")
				+ to_string((int)dm_by_second_glucose)
				+ " mg/dL or HbA1C >= "
				+ to_string(dm_by_second_hba1c)
				+ "% whithin "
				+ to_string(dm_by_second_time_delta_days)
				+ " days";
			if (is_recent)
				ret_recent |= REASON_RECENT_LABS;
			else ret_past |= REASON_PAST_LABS;
		}
		if (de.de_type == DFG_DIABETES_EVENT_DRUG) {
			reason_str = string("Drugs evidence in ")+to_string(de.time);
			if (is_recent)
				ret_recent |= REASON_RECENT_DRUGS;
			else ret_past |= REASON_PAST_DRUGS;
		}
		if (de.de_type == DFG_DIABETES_EVENT_DIAGNOSIS) {
			reason_str = string("Diagnosis evidence in "+to_string(de.time));
			if (is_recent)
				ret_recent |= REASON_RECENT_DRUGS;
			else ret_past |= REASON_PAST_DRUGS;
		}
		//MLOG("reason_str = %s\n", reason_str.c_str());
		if (reason_str.length() != 0)
		{
			if (is_recent)
				recent_evidence_text.insert(reason_str);
			else
				past_evidence_text.insert(reason_str);
		}
	}
	
	json json_evidence = json::array();

	if (ret_recent != 0) {
		ret = ret_recent;
		for (const auto &v : recent_evidence_text)
			json_evidence.push_back(v);
	}
	else if (ret_past != 0) {
		ret = ret_past;
		for (const auto &v : past_evidence_text)
			json_evidence.push_back(v);
	}
	else ret = 1000;
	
	json_out = json::object();
	json_out["Evidence"] = json_evidence;

	//string json_str = json_out.dump();
	//MLOG("json_out = %s\nret = %d\n", json_str.c_str(), ret);


	return ret;
}


//=======================================================================================
// DiabetesFinderGenerator : Implementation for the Diabetes Finder AM 
//=======================================================================================


// Generate
//.......................................................................................
int DiabetesFinderGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {

	vector<DiabetesEvent> dm_events;
	UniversalSigVec usv;

	rec.uget(dm_glucose_sig, 0, usv);
	for (int i = 0; i < usv.len; ++i) {
		dm_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_GLUCOSE, usv.Time(i), usv.Val(i)));
	}

	rec.uget(dm_hba1c_sig, 0, usv);
	for (int i = 0; i < usv.len; ++i) {
		dm_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_HBA1C, usv.Time(i), usv.Val(i)));
	}

	rec.uget(dm_drug_sig, 0, usv);
	for (int i = 0; i < usv.len; ++i) {
		int i_val = usv.Val<int>(i);
		if (i_val < 0 || i_val > dm_drug_lut.size())
			MTHROW_AND_ERR("ERROR in DiabetesFinderGenerator: got i_val=%d while DRUG lut size is %d\n", i_val, (int)dm_drug_lut.size());
		if (dm_drug_lut.size() > 0 && dm_drug_lut[i_val]) {
			//TODO: Add Exception for metformin??
			dm_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_DRUG, usv.Time(i), usv.Val(i)));
		}
	}
	
	rec.uget(dm_diagnosis_sig, 0, usv);
	for (int i = 0; i < usv.len; ++i) {
		int i_val = usv.Val<int>(i);
		if (i_val < 0 || i_val > dm_diagnosis_lut.size())
			MTHROW_AND_ERR("ERROR in DiabetesFinderGenerator: got i_val=%d while RC lut size is %d\n", i_val, (int)dm_diagnosis_lut.size());
		if (dm_diagnosis_lut.size() > 0 && dm_diagnosis_lut[i_val]) {
			dm_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_DIAGNOSIS, usv.Time(i), usv.Val(i)));
		}
	}

	// sorting events
	sort(dm_events.begin(), dm_events.end(), [](const DiabetesEvent &v1, const DiabetesEvent &v2) { return v1.time < v2.time; });

	float *p_feat = _p_data[0] + index;

	for (int i = 0; i < num; i++)
		p_feat[i] = 0;

	for (int i = 0; i < num; i++) {
		int s_time = features.samples[index + i].time;
		json json_out;
		features.samples[index + i].prediction.push_back(_resolve(dm_events, s_time, json_out));
		features.samples[index + i].str_attributes["Explanations"] = json_out.dump();
	}

/*
	float *p_feat = _p_data[0] + index;

	for (int i = 0; i < num; i++)
		p_feat[i] = 0;

	for (int i = 0; i < num; i++) {
		int s_time = features.samples[index + i].time;
		int last_glucose_idx = -1;
		for (int j = usv.len-1; j >= 0 ; --j) {
			if (usv.Time(j) <= s_time)
			{
				last_glucose_idx = j;
				break;
			}
		}
		if (last_glucose_idx == -1)
			continue;
		float glucose = usv.Val(last_glucose_idx);

		if (glucose < 100.0f)
		{
			p_feat[i] = 1.0f;
			features.samples[index + i].str_attributes["DiabetesStatus"] = "Non-diabetic";
			features.samples[index + i].str_attributes["Urgency"] = "None";
		}
		else if (glucose < 110) {
			p_feat[i] = 2.0f;
			features.samples[index + i].str_attributes["DiabetesStatus"] = "Pre-diabetic-level1";
			features.samples[index + i].str_attributes["Urgency"] = "Low";
		}
		else if (glucose < 120) {
			p_feat[i] = 3.0f;
			features.samples[index + i].str_attributes["DiabetesStatus"] = "Pre-diabetic-level2";
			features.samples[index + i].str_attributes["Urgency"] = "Intermediate";
		}
		else { 
			p_feat[i] = 4.0f; 
			features.samples[index + i].str_attributes["DiabetesStatus"] = "Diabetic";
			features.samples[index + i].str_attributes["Urgency"] = "High";
		}

	}
*/

#if 0
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

#endif // 0



	return 0;
}

// Init
//.......................................................................................
int DiabetesFinderGenerator::init(map<string, string>& mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "dm_diagnosis_sets") boost::split(dm_diagnosis_sets, entry.second, boost::is_any_of(","));
		else if (field == "dm_coded_sets") boost::split(dm_coded_sets, entry.second, boost::is_any_of(","));
		else if (field == "dm_drug_sets") boost::split(dm_drug_sets, entry.second, boost::is_any_of(","));
		else if (field == "dm_diagnosis_sig") dm_diagnosis_sig = entry.second; // "RC";
		else if (field == "dm_coded_sig") dm_coded_sig = entry.second; // "RC";
		else if (field == "dm_glucose_sig") dm_glucose_sig = entry.second; // "Glucose";
		else if (field == "dm_hba1c_sig") dm_hba1c_sig = entry.second; // "HbA1C";
		else if (field == "dm_drug_sig") dm_drug_sig = entry.second; // "Drug";
		else if (field == "dm_past_event_days") dm_past_event_days = med_stoi(entry.second); //(365) * 3;
		else if (field == "dm_by_single_glucose") dm_by_single_glucose = med_stof(entry.second); //200.0f;
		else if (field == "dm_by_second_glucose") dm_by_second_glucose = med_stof(entry.second); //126.0f;
		else if (field == "dm_by_second_hba1c") dm_by_second_hba1c = med_stof(entry.second); //6.5f;
		else if (field == "dm_by_second_time_delta_days") dm_by_second_time_delta_days = med_stoi(entry.second); //(365) * 2;
		else if (field == "dm_by_second_time_delta") dm_by_second_time_delta = med_stoi(entry.second); //-1;
		
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter [%s] for DiabetesFinderGenerator\n", field.c_str());
	}
	set_names();

	req_signals.clear();
	req_signals.push_back(dm_glucose_sig);
	req_signals.push_back(dm_hba1c_sig);
	req_signals.push_back(dm_diagnosis_sig);
	req_signals.push_back(dm_drug_sig);
	if(dm_coded_sig != dm_diagnosis_sig)
		req_signals.push_back(dm_coded_sig);

	return 0;
}

void DiabetesFinderGenerator::init_tables(MedDictionarySections& dict) {
	if (dm_drug_lut.size() == 0) {
		dict.prep_sets_indexed_lookup_table(dict.section_id(dm_drug_sig), dm_drug_sets, dm_drug_lut);
	}
	if (dm_diagnosis_lut.size() == 0) {
		dict.prep_sets_indexed_lookup_table(dict.section_id(dm_diagnosis_sig), dm_diagnosis_sets, dm_diagnosis_lut);
	}

	return;
}

void DiabetesFinderGenerator::init_defaults() {

}

