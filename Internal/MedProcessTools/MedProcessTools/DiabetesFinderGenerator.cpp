#include "DiabetesFinderGenerator.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int DiabetesFinderGenerator::_resolve(vector<DiabetesEvent>& df_events, int calc_time, json& json_out) {
	int ret = 0;

	unordered_set<string> past_evidence_text, recent_evidence_text;
	int ret_recent = 0;
	int ret_past = 0;

	//mark events by the 2 events rull
	int _latest_noteable_event_time = -1;
	for (auto& de : df_events) {
		if (de.time > calc_time)
			continue;
		if ((de.de_type == DFG_DIABETES_EVENT_GLUCOSE && de.val >= df_by_second_glucose)
			|| (de.de_type == DFG_DIABETES_EVENT_HBA1C && de.val >= df_by_second_hba1c)) {
			if (_latest_noteable_event_time != -1 && de.time - _latest_noteable_event_time < df_by_second_time_delta) {
				de.is_second = true;
			}
			_latest_noteable_event_time = de.time;
		}
	}

	for (const auto& de : df_events) {
		//MLOG("DE: (%d,%d,%f)\n", de.de_type, de.time, de.val);
		//MLOG("DE:    calc_time = %d\n", calc_time);
		if (de.time > calc_time)
			continue;
		bool is_recent = de.time >= med_time_converter.add_subtruct_days(calc_time, -1* df_past_event_days);
		//MLOG("DE:    add_subtruct_days(calc_time, -1 * df_past_event_days) = %d\n", med_time_converter.add_subtruct_days(calc_time, -1 * df_past_event_days));
		//MLOG("DE:    is_recent = %d\n", (int)is_recent);
		//MLOG("DE:    df_by_single_glucose = %d\n", (int)df_by_single_glucose);
		//MLOG("DE:    (de.de_type == DFG_DIABETES_EVENT_GLUCOSE) = %d\n", (int)(de.de_type == DFG_DIABETES_EVENT_GLUCOSE));
		string reason_str = "";
		if (de.de_type == DFG_DIABETES_EVENT_GLUCOSE && de.val >= df_by_single_glucose) {
			reason_str = string("Single measurement of glucose >=") + to_string((int)df_by_single_glucose) + " mg/dL";
			if (is_recent)
				ret_recent |= REASON_RECENT_LABS;
			else ret_past |= REASON_PAST_LABS;
		}
		if (de.is_second) {
			reason_str = string("Second measurement of glucose >= ")
				+ to_string((int)df_by_second_glucose)
				+ " mg/dL or HbA1C >= "
				+ to_string(df_by_second_hba1c)
				+ "% whithin "
				+ to_string(df_by_second_time_delta_days)
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
	
	json_out = json::object();
	json_out["Evidence"] = json_evidence;

	//string json_str = json_out.dump();
	//MLOG("json_out = %s\nret = %d\n", json_str.c_str(), ret);

	if (df_score_is_flag) {
		if (ret_recent != 0)
			return 1;
		if (ret_past != 0)
			return 2;
		return 0;
	}

	return ret;
}


//=======================================================================================
// DiabetesFinderGenerator : Implementation for the Diabetes Finder AM 
//=======================================================================================


// Generate
//.......................................................................................
int DiabetesFinderGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {

	vector<DiabetesEvent> df_events;
	UniversalSigVec usv;

	rec.uget(df_glucose_sig, 0, usv);
	for (int i = 0; i < usv.len; ++i) {
		df_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_GLUCOSE, usv.Time(i), usv.Val(i)));
	}

	rec.uget(df_hba1c_sig, 0, usv);
	for (int i = 0; i < usv.len; ++i) {
		df_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_HBA1C, usv.Time(i), usv.Val(i)));
	}

	rec.uget(df_drug_sig, 0, usv);
	for (int i = 0; i < usv.len; ++i) {
		int i_val = usv.Val<int>(i);
		if (i_val < 0 || i_val > df_drug_lut.size())
			MTHROW_AND_ERR("ERROR in DiabetesFinderGenerator: got i_val=%d while DRUG lut size is %d\n", i_val, (int)df_drug_lut.size());
		if (df_drug_lut.size() > 0 && df_drug_lut[i_val]) {
			//TODO: Add Exception for metformin??
			df_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_DRUG, usv.Time(i), usv.Val(i)));
		}
	}
	
	rec.uget(df_diagnosis_sig, 0, usv);
	for (int i = 0; i < usv.len; ++i) {
		int i_val = usv.Val<int>(i);
		if (i_val < 0 || i_val > df_diagnosis_lut.size())
			MTHROW_AND_ERR("ERROR in DiabetesFinderGenerator: got i_val=%d while RC lut size is %d\n", i_val, (int)df_diagnosis_lut.size());
		if (df_diagnosis_lut.size() > 0 && df_diagnosis_lut[i_val]) {
			df_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_DIAGNOSIS, usv.Time(i), usv.Val(i)));
		}
	}

	// sorting events
	sort(df_events.begin(), df_events.end(), [](const DiabetesEvent &v1, const DiabetesEvent &v2) { return v1.time < v2.time; });

	float *p_feat = _p_data[0] + index;

	for (int i = 0; i < num; i++)
		p_feat[i] = 0;

	for (int i = 0; i < num; i++) {
		int s_time = features.samples[index + i].time;
		json json_out;
		features.samples[index + i].prediction.push_back(_resolve(df_events, s_time, json_out));
		features.samples[index + i].str_attributes["Explanations"] = json_out.dump();
	}

	return 0;
}

// Init
//.......................................................................................
int DiabetesFinderGenerator::init(map<string, string>& mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "df_score_mode") {
			df_score_is_flag = (boost::to_upper_copy(entry.second) == "FLAG");
			df_score_is_bitmask = (boost::to_upper_copy(entry.second) == "BITMASK");
		}
		else if (field == "df_diagnosis_sets") boost::split(df_diagnosis_sets, entry.second, boost::is_any_of(","));
		else if (field == "df_coded_sets") boost::split(df_coded_sets, entry.second, boost::is_any_of(","));
		else if (field == "df_drug_sets") boost::split(df_drug_sets, entry.second, boost::is_any_of(","));
		else if (field == "df_diagnosis_sig") df_diagnosis_sig = entry.second; // "RC";
		else if (field == "df_coded_sig") df_coded_sig = entry.second; // "RC";
		else if (field == "df_glucose_sig") df_glucose_sig = entry.second; // "Glucose";
		else if (field == "df_hba1c_sig") df_hba1c_sig = entry.second; // "HbA1C";
		else if (field == "df_drug_sig") df_drug_sig = entry.second; // "Drug";
		else if (field == "df_past_event_days") df_past_event_days = med_stoi(entry.second); //(365) * 3;
		else if (field == "df_by_single_glucose") df_by_single_glucose = med_stof(entry.second); //200.0f;
		else if (field == "df_by_second_glucose") df_by_second_glucose = med_stof(entry.second); //126.0f;
		else if (field == "df_by_second_hba1c") df_by_second_hba1c = med_stof(entry.second); //6.5f;
		else if (field == "df_by_second_time_delta_days") df_by_second_time_delta_days = med_stoi(entry.second); //(365) * 2;
		else if (field == "df_by_second_time_delta") df_by_second_time_delta = med_stoi(entry.second); //-1;
		
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter [%s] for DiabetesFinderGenerator\n", field.c_str());
	}
	set_names();

	req_signals.clear();
	req_signals.push_back(df_glucose_sig);
	req_signals.push_back(df_hba1c_sig);
	req_signals.push_back(df_diagnosis_sig);
	req_signals.push_back(df_drug_sig);
	if(df_coded_sig != df_diagnosis_sig)
		req_signals.push_back(df_coded_sig);

	return 0;
}

void DiabetesFinderGenerator::init_tables(MedDictionarySections& dict) {
	if (df_drug_lut.size() == 0) {
		dict.prep_sets_indexed_lookup_table(dict.section_id(df_drug_sig), df_drug_sets, df_drug_lut);
	}
	if (df_diagnosis_lut.size() == 0) {
		dict.prep_sets_indexed_lookup_table(dict.section_id(df_diagnosis_sig), df_diagnosis_sets, df_diagnosis_lut);
	}

	return;
}

void DiabetesFinderGenerator::init_defaults() {

}

