#include "DiabetesFinderGenerator.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//-------------------------------------------------------------------------------------------------------------
int DiabetesFinderGenerator::_resolve(PidDynamicRec& rec, vector<DiabetesEvent>& df_events, int coded_date, int coded_val, int calc_time, json& json_out) {
	
	
	//
	// we return:
	// 0: if no coding <= calc_time found or if already coded (coded_date <= calc_time)
	// 1: if new coding that lasts into the range [calc_time - df_past_event_days] happened.
	// 2: if a new coding happened with no indication in the last df_past_event_days..
	//

	json j_explain = json({});
	if (coded_date <= calc_time) {
		json_out += {"code", 0};
		j_explain += {"sig", df_coded_sig};
		j_explain += {"date", coded_date};
		string name = rec.my_base_rep->dict.dicts[df_coded_section_id].Id2LongestName[coded_val];
		j_explain += {"val", name};
		j_explain += {"reason", "already coded"};

		json_out += {"explanation", j_explain};
		return 0;
	}

	int last_indication_date = -1;
	int first_coding_date = -1;
	if (df_output_verbosity >= 3)
		json_out += { "indications", json::array() };

	int n_good_indications = 0;
	int start_good_date = -1;
	int pre_d_date = -1; // last day of pre diabetes and up
	for (auto de : df_events) {

		if (de.time > calc_time)
			break;

		if ((de.de_type == DFG_DIABETES_EVENT_GLUCOSE && de.val >= df_pre_d_glucose)
			|| (de.de_type == DFG_DIABETES_EVENT_HBA1C && de.val >= df_pre_d_hba1c)) {
			pre_d_date = de.time;
		}

		if (first_coding_date < 0) {

			if (de.de_type == DFG_DIABETES_EVENT_GLUCOSE && (de.is_first || de.is_second)) {
				j_explain += {"sig", df_glucose_sig};
				j_explain += { "date", de.time };
				j_explain += { "value", (int)de.val };
				j_explain += { "reason", de.reason };
				first_coding_date = de.time;
			}

			if (de.de_type == DFG_DIABETES_EVENT_HBA1C && (de.is_first || de.is_second)) {
				j_explain += {"sig", df_hba1c_sig};
				j_explain += { "date", de.time };
				j_explain += { "value", de.val };
				j_explain += { "reason", de.reason };
				first_coding_date = de.time;
			}

			if (de.de_type == DFG_DIABETES_EVENT_DRUG) {
				string name = rec.my_base_rep->dict.dicts[df_drug_section_id].Id2LongestName[de.val];
				j_explain += {"sig", df_drug_sig};
				j_explain += { "date", de.time };
				j_explain += { "value", name };
				j_explain += { "reason", de.reason };
				first_coding_date = de.time;
			}

			if (de.de_type == DFG_DIABETES_EVENT_DIAGNOSIS) {
				string name = rec.my_base_rep->dict.dicts[df_drug_section_id].Id2LongestName[de.val];
				j_explain += {"sig", df_diagnosis_sig};
				j_explain += { "date", de.time };
				j_explain += { "value", name };
				j_explain += { "reason", de.reason };
				first_coding_date = de.time;
			}

		}

		if (de.de_type == DFG_DIABETES_EVENT_DRUG || de.de_type == DFG_DIABETES_EVENT_DIAGNOSIS || de.is_first || de.is_second) {
			last_indication_date = de.time;
			n_good_indications = 0;
		}
		else {
			if (n_good_indications == 0) start_good_date = de.time;
			n_good_indications++;
		}

		if (!de.is_non_dm && df_output_verbosity>=3) {
			json j_ind;
			switch (de.de_type) {
			case DFG_DIABETES_EVENT_GLUCOSE:
				j_ind = { {"sig", df_glucose_sig}, { "date", de.time }, { "value" , (int)de.val }};
				break;
			case DFG_DIABETES_EVENT_HBA1C:
				j_ind = { {"sig", df_hba1c_sig}, { "date", de.time }, { "value" , de.val }};
				break;
			case DFG_DIABETES_EVENT_DRUG:
				j_ind = { {"sig", df_drug_sig}, { "date", de.time }, { "value" , rec.my_base_rep->dict.dicts[df_drug_section_id].Id2LongestName[de.val] }};
				break;
			case DFG_DIABETES_EVENT_DIAGNOSIS:
				j_ind = { {"sig", df_diagnosis_sig}, { "date", de.time }, { "value" , rec.my_base_rep->dict.dicts[df_drug_section_id].Id2LongestName[de.val] }};
				break;
			}
			
			json_out["indications"].push_back(j_ind);
			
		}
	}
	
	int code = 0;
	if (last_indication_date > 0) {
		int past_date = med_time_converter.add_subtruct_days(calc_time, -df_past_event_days);
		if (last_indication_date >= past_date)
			code = 1;
		else
			code = 2;
	}

	json_out += {"code", code};

	if (code == 0) {
		if (pre_d_date > 0) json_out += {"pre diabetic", pre_d_date};
		j_explain += {"reason", "no indication"};
	}
	json_out += {"explanation", j_explain};
	if (code > 0)
		json_out += {"dm_date", first_coding_date};
	if (last_indication_date > 0)
		json_out += {"last_indication_date", last_indication_date};
	if (n_good_indications > 0 && df_output_non_dm_period >= 2) {
		json_out += {"n_last_non_dm_indications", n_good_indications};
		json_out += {"start_non_dm_date", start_good_date};
	}

	return code;

	/*
	int ret = 0;

	unordered_set<string> past_evidence_text, recent_evidence_text;
	int ret_recent = 0;
	int ret_past = 0;





	for (const auto& de : df_events) {
		//MLOG("DE: (%d,%d,%f)\n", de.de_type, de.time, de.val);
		//MLOG("DE:    calc_time = %d\n", calc_time);
		if (de.time > calc_time)
			continue;

		bool is_recent = de.time >= med_time_converter.add_subtruct_days(calc_time, -df_past_event_days);
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
		if (de.de_type == DFG_DIABETES_EVENT_HBA1C && de.val >= df_by_single_hba1c) {
			reason_str = string("Single measurement of HbA1C >=") + to_string((int)df_by_single_hba1c) + " %";
			if (is_recent)
				ret_recent |= REASON_RECENT_LABS;
			else ret_past |= REASON_PAST_LABS;
		}
		if (de.is_second) {
			reason_str = string("Second measurement of glucose >= ")
				+ to_string((int)df_by_second_glucose)
				+ " mg/dL or HbA1C >= "
				+ to_string(df_by_second_hba1c)
				+ "% within "
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
	
	//json_out = json::object();
	json_out = json_evidence;

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
	*/
}


//=======================================================================================
// DiabetesFinderGenerator : Implementation for the Diabetes Finder AM 
//=======================================================================================


// Generate
//.......................................................................................
int DiabetesFinderGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {

	vector<DiabetesEvent> df_events;
	UniversalSigVec usv;

	//MLOG("In here with pid %d\n", rec.pid);

	// we do here a pass over all data, we will later do the time calculations for each prediction point.
	if (df_glucose_sid > 0) {
		rec.uget(df_glucose_sid, 0, usv);
		for (int i = 0; i < usv.len; ++i) {
			df_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_GLUCOSE, usv.Time(i), usv.Val(i)));
		}
	}

	if (df_hba1c_sid > 0) {
		rec.uget(df_hba1c_sid, 0, usv);
		for (int i = 0; i < usv.len; ++i) {
			df_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_HBA1C, usv.Time(i), usv.Val(i)));
		}
	}

	if (df_drug_sid > 0) {
		rec.uget(df_drug_sid, 0, usv);
		for (int i = 0; i < usv.len; ++i) {
			int i_val = usv.Val<int>(i);
			if (i_val < 0 || (i_val > (int)df_drug_lut.size()))
				MTHROW_AND_ERR("ERROR in DiabetesFinderGenerator: drug got i_val=%d while DRUG lut size is %d\n", i_val, (int)df_drug_lut.size());
			if (df_drug_lut.size() > 0 && df_drug_lut[i_val]) {
				//TODO: Add Exception for metformin??
				df_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_DRUG, usv.Time(i), usv.Val(i)));
				//break; // only first usage is interesting
			}
		}
	}
	

	if (df_diagnosis_sid > 0) {
		rec.uget(df_diagnosis_sid, 0, usv);
		for (int i = 0; i < usv.len; ++i) {
			int i_val = usv.Val<int>(i);
			if (i_val < 0 || (i_val > (int)df_diagnosis_lut.size()))
				MTHROW_AND_ERR("ERROR in DiabetesFinderGenerator: diagnosis got i_val=%d while RC lut size is %d\n", i_val, (int)df_diagnosis_lut.size());
			if (df_diagnosis_lut.size() > 0 && df_diagnosis_lut[i_val]) {
				df_events.push_back(DiabetesEvent(DFG_DIABETES_EVENT_DIAGNOSIS, usv.Time(i), usv.Val(i)));
				//break; // only first usage is interesting
			}
		}
	}

	// sorting events
	sort(df_events.begin(), df_events.end(), [](const DiabetesEvent &v1, const DiabetesEvent &v2) { return v1.time < v2.time; });

	// find first coded
	int coded_date = 29999999; // assuming this code will not be used in the year 3000.
	int coded_val = -1;
	if (df_coded_sid > 0) {
		rec.uget(df_coded_sid, 0, usv);
		for (int i = 0; i < usv.len; ++i) {
			int i_val = usv.Val<int>(i);
			if (i_val < 0 || (i_val > (int)df_coded_lut.size()))
				MTHROW_AND_ERR("ERROR in DiabetesFinderGenerator: coded got i_val=%d while RC lut size is %d\n", i_val, (int)df_diagnosis_lut.size());
			if (df_coded_lut.size() > 0 && df_diagnosis_lut[i_val]) {
				coded_date = usv.Time(i);
				coded_val = i_val;
				break; // only first usage is interesting
			}
		}
	}


	// mark the 2 events rule cases (if there is one)
	int _latest_noteable_event_time = -1;
	for (auto& de : df_events) {

		if ((de.de_type == DFG_DIABETES_EVENT_GLUCOSE && de.val >= df_by_single_glucose)
			|| (de.de_type == DFG_DIABETES_EVENT_HBA1C && de.val >= df_by_single_hba1c)) {
			de.is_first = true;
			de.reason = "single bad test";
		}

		if ((de.de_type == DFG_DIABETES_EVENT_GLUCOSE && de.val >= df_by_second_glucose)
			|| (de.de_type == DFG_DIABETES_EVENT_HBA1C && de.val >= df_by_second_hba1c)) {
			if (_latest_noteable_event_time != -1 &&
				med_time_converter.add_subtruct_days(de.time, -df_by_second_time_delta_days) <= _latest_noteable_event_time)
			{
				de.is_second = true;
			}
			//MLOG("type: %d, time: %d, val: %f , is_second: %d\n  _latest_noteable_event_time=%d\n  de.time - _latest_noteable_event_time = %d\n",de.de_type, de.time, de.val, (int)de.is_second, _latest_noteable_event_time, de.time - _latest_noteable_event_time);
			_latest_noteable_event_time = de.time;
			if (de.reason == "")
				de.reason = "second bad test";
		}

		if ((de.de_type == DFG_DIABETES_EVENT_GLUCOSE && de.val < df_by_second_glucose)
			|| (de.de_type == DFG_DIABETES_EVENT_HBA1C && de.val < df_by_second_hba1c)) {
			de.is_non_dm = true;
		}

		if (de.de_type == DFG_DIABETES_EVENT_DIAGNOSIS) de.reason = "indicative diagnosis";
		if (de.de_type == DFG_DIABETES_EVENT_DRUG) de.reason = "indicative drug";
	}

	float *p_feat = _p_data[0] + index;

	for (int i = 0; i < num; i++)
		p_feat[i] = 0;

	for (int i = 0; i < num; i++) {
		int s_time = features.samples[index + i].time;
		json json_out = json({});
		p_feat[i] = _resolve(rec, df_events, coded_date, coded_val, s_time, json_out);
		features.samples[index + i].prediction.push_back(p_feat[i]);
		features.samples[index + i].jrec += {"DiabetesCoder", json_out};
		features.samples[index + i].str_attributes["DiabetesCoder"] = json_out.dump();
	}

	return 0;
}

// Init
//.......................................................................................
int DiabetesFinderGenerator::init(map<string, string>& mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
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
		else if (field == "df_by_single_hba1c") df_by_single_hba1c = med_stof(entry.second); //7.0f;
		else if (field == "df_by_second_glucose") df_by_second_glucose = med_stof(entry.second); //126.0f;
		else if (field == "df_by_second_hba1c") df_by_second_hba1c = med_stof(entry.second); //6.5f;
		else if (field == "df_pre_d_hba1c") df_pre_d_hba1c = med_stof(entry.second); //5.8f;
		else if (field == "df_pre_d_glucose") df_pre_d_glucose = med_stof(entry.second); //101f;
		else if (field == "df_by_second_time_delta_days") df_by_second_time_delta_days = med_stoi(entry.second); //(365) * 2;
		else if (field == "df_output_verbosity") df_output_verbosity = med_stoi(entry.second); // 2;
		else if (field == "df_output_non_dm_period") df_output_non_dm_period = med_stoi(entry.second); // 0;

		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter [%s] for DiabetesFinderGenerator\n", field.c_str());
	}
	set_names();


	if (df_drug_sig == "NONE") df_drug_sig = "";
	if (df_diagnosis_sig == "NONE") df_diagnosis_sig = "";
	if (df_coded_sig == "NONE") df_coded_sig = "";


	req_signals.clear();
	req_signals.push_back(df_glucose_sig);
	req_signals.push_back(df_hba1c_sig);
	if (df_diagnosis_sig != "") req_signals.push_back(df_diagnosis_sig);
	if (df_drug_sig != "") req_signals.push_back(df_drug_sig);
	if(df_coded_sig != df_diagnosis_sig)
		if (df_coded_sig != "") req_signals.push_back(df_coded_sig);

	return 0;
}

//-------------------------------------------------------------------------------------------------------------
void DiabetesFinderGenerator::init_tables(MedDictionarySections& dict) {
	if (df_drug_lut.size() == 0 && df_drug_sig != "") {
		dict.prep_sets_indexed_lookup_table(dict.section_id(df_drug_sig), df_drug_sets, df_drug_lut);
	}

	if (df_diagnosis_lut.size() == 0 && df_diagnosis_sig != "") {
		dict.prep_sets_indexed_lookup_table(dict.section_id(df_diagnosis_sig), df_diagnosis_sets, df_diagnosis_lut);
	}

	if (df_coded_lut.size() == 0 && df_coded_sig != "") {
		dict.prep_sets_indexed_lookup_table(dict.section_id(df_coded_sig), df_coded_sets, df_coded_lut);
	}

	df_drug_section_id = dict.section_id(df_drug_sig);
	df_diagnosis_section_id = dict.section_id(df_diagnosis_sig);
	df_coded_section_id = dict.section_id(df_coded_sig);


	return;
}


//-------------------------------------------------------------------------------------------------------------
void DiabetesFinderGenerator::set_signal_ids(MedSignals& sigs)
{
	df_drug_sid = sigs.sid(df_drug_sig);
	df_glucose_sid = sigs.sid(df_glucose_sig);
	df_hba1c_sid = sigs.sid(df_hba1c_sig);
	df_diagnosis_sid = sigs.sid(df_diagnosis_sig);
	df_coded_sid = sigs.sid(df_coded_sig);
}

//-------------------------------------------------------------------------------------------------------------
void DiabetesFinderGenerator::init_defaults() {

}

