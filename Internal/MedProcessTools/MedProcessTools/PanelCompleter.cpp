#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "RepProcess.h"

//=======================================================================================
// RepPanelCompleter fills-in calculatable signal values. Enriching existing signals
//=======================================================================================
//.......................................................................................
// Init from map
int RepPanelCompleter::init(map<string, string>& mapper)
{

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepPanelCompleter::init]
		if (panel2type.find(field) != panel2type.end()) {
			if (update_signal_names(field, entry.second) != 0) return -1;
		}
		else if (field == "missing") missing_val = stof(entry.second);
		else if (field == "config" || field == "metadata") metadata_file = entry.second;
		else if (field == "panels") {
			if (update_panels(entry.second) != 0) return -1;
		}
		//! [RepPanelCompleter::init]
	}

	init_lists();
	return 0;
}

//.......................................................................................
// Update signal names
int RepPanelCompleter::update_signal_names(string panel, string& names) {

	vector<string> signals;
	boost::split(signals, names, boost::is_any_of(","));

	if (signals.size() != 0 && signals.size() != panel2signals[panel].size()) {
		MERR("Wrong number of signals given for panel %s", panel.c_str());
		return -1;
	}
	else {
		panel_signal_names[panel2type[panel]] = signals;
		return 0;
	}
}

//.......................................................................................
// Update panels to work on.
int RepPanelCompleter::update_panels(string& panels) {

	vector<string> list;
	boost::split(list, panels, boost::is_any_of(","));

	unordered_set<string> selected;
	for (string& panel : list) {
		if (panel2signals.find(panel) == panel2signals.end()) {
			MERR("Required unknown panel %s for completion\n", panel.c_str());
			return -1;
		}
		else
			selected.insert(panel);
	}

	for (auto& panel : panel2signals) {
		if (selected.find(panel.first) == selected.end())
			panel_signal_names[panel2type[panel.first]].clear();
	}
	return 0;
}

//.......................................................................................
// Set signal names to default
void RepPanelCompleter::init_defaults() {

	panel_signal_names.resize(panel2type.size());
	for (auto& panel : panel2type)
		panel_signal_names[panel.second] = panel2signals[panel.first];

	ageDirectlyGiven = med_rep_type.ageDirectlyGiven;
	genderSignalName = med_rep_type.genderSignalName;
}

//.......................................................................................
// initialize signal ids
void RepPanelCompleter::set_signal_ids(MedDictionarySections& dict) {

	panel_signal_ids.resize(panel_signal_names.size());
	for (int iPanel = 0; iPanel < panel_signal_ids.size(); iPanel++) {
		panel_signal_ids[iPanel].resize(panel_signal_names[iPanel].size());
		for (int iSig = 0; iSig < panel_signal_ids[iPanel].size(); iSig++) {
			panel_signal_ids[iPanel][iSig] = dict.id(panel_signal_names[iPanel][iSig]);
			if (panel_signal_ids[iPanel][iSig] == -1)
				MTHROW_AND_ERR("Cannot find signal-id for %s\n", panel_signal_names[iPanel][iSig].c_str());
		}
	}

	// EGFR Requires age and gender
	if (panel_signal_names[REP_CMPLT_EGFR_PANEL].size()) {
		genderId = dict.id(genderSignalName);

		if (ageDirectlyGiven)
			ageId = dict.id("Age");
		else
			byearId = dict.id("BYEAR");

	}
}

//.......................................................................................
// Fill required and affected signals
void RepPanelCompleter::init_lists() {

	req_signals.clear();
	aff_signals.clear();

	for (int iPanel = 0; iPanel < panel_signal_names.size(); iPanel++) {
		for (int iSig = 0; iSig < panel_signal_names[iPanel].size(); iSig++) {
			req_signals.insert(panel_signal_names[iPanel][iSig]);
			aff_signals.insert(panel_signal_names[iPanel][iSig]);
		}
	}

	if (panel_signal_names[REP_CMPLT_EGFR_PANEL].size()) {
		req_signals.insert(genderSignalName);

		if (ageDirectlyGiven)
			req_signals.insert("Age");
		else
			req_signals.insert("BYEAR");

	}
}

// Apply completions (no relevant attributes)
//.......................................................................................
int RepPanelCompleter::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {

	// Check that we have the correct number of dynamic-versions : one per time-point (if given)
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int rc = 0;
	if (panel_signal_ids[REP_CMPLT_RED_LINE_PANEL].size())
		rc = apply_red_line_completer(rec, time_points);
	if (rc < 0)
		return -1;

	if (panel_signal_ids[REP_CMPLT_WHITE_LINE_PANEL].size())
		rc = apply_white_line_completer(rec, time_points);
	if (rc < 0)
		return -1;

	if (panel_signal_ids[REP_CMPLT_PLATELETS_PANEL].size())
		rc = apply_platelets_completer(rec, time_points);
	if (rc < 0)
		return -1;

	if (panel_signal_ids[REP_CMPLT_LIPIDS_PANEL].size())
		rc = apply_lipids_completer(rec, time_points);
	if (rc < 0)
		return -1;

	if (panel_signal_ids[REP_CMPLT_EGFR_PANEL].size())
		rc = apply_eGFR_completer(rec, time_points);
	if (rc < 0)
		return -1;

	if (panel_signal_ids[REP_CMPLT_BMI_PANEL].size())
		rc = apply_BMI_completer(rec, time_points);
	if (rc < 0)
		return -1;

	return 0;
}

// Completion of red blood line panels
//.......................................................................................
int RepPanelCompleter::apply_red_line_completer(PidDynamicRec& rec, vector<int>& time_points) {

	vector<int>& sigs_ids = panel_signal_ids[REP_CMPLT_RED_LINE_PANEL];
	int n_sigs = (int)sigs_ids.size();

	vector<float>& orig_res = original_sig_res[REP_CMPLT_RED_LINE_PANEL];
	vector<float>& final_res = final_sig_res[REP_CMPLT_RED_LINE_PANEL];
	vector<float>& conv = sig_conversion_factors[REP_CMPLT_RED_LINE_PANEL];

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end());

	allVersionsIterator vit(rec, iteratorSignalIds);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		int time_limit = (time_points.size()) ? time_points[iver] : -1;

		// Get Signals
		rec.usvs.resize(n_sigs);
		for (size_t i = 0; i < n_sigs; ++i)
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);

		// Get Panels
		vector<int> panel_times;
		vector<vector<float> > panels;
		get_panels(rec.usvs, panel_times, panels, time_limit, RED_PNL_LAST);

		// Complete
		vector<int> changed_signals(n_sigs, 0);
		for (size_t i = 0; i < panels.size(); i++) {
			int complete = 1;
			while (complete) {
				complete = 0;

				// MCV = 10*HCT/RBC
				complete += triplet_complete(panels[i], 10, RED_PNL_MCV, RED_PNL_HCT, RED_PNL_RBC, orig_res, final_res, conv, changed_signals);

				// MCH = 10 * HGB / RBC
				complete += triplet_complete(panels[i], 10, RED_PNL_MCH, RED_PNL_HGB, RED_PNL_RBC, orig_res, final_res, conv, changed_signals);

				// MCHC = 100 * HGB / HCT
				complete += triplet_complete(panels[i], 100, RED_PNL_MCHC, RED_PNL_HGB, RED_PNL_HCT, orig_res, final_res, conv, changed_signals);
			}
		}

		// Update changed signals
		if (update_signals(rec, iver, panels, panel_times, sigs_ids, changed_signals) < 0)
			return -1;
	}

	return 0;
}

// Completion of white blood line panels
//.......................................................................................
int RepPanelCompleter::apply_white_line_completer(PidDynamicRec& rec, vector<int>& time_points) {

	vector<int>& sigs_ids = panel_signal_ids[REP_CMPLT_WHITE_LINE_PANEL];
	int n_sigs = (int)sigs_ids.size();

	vector<float>& orig_res = original_sig_res[REP_CMPLT_WHITE_LINE_PANEL];
	vector<float>& final_res = final_sig_res[REP_CMPLT_WHITE_LINE_PANEL];
	vector<float>& conv = sig_conversion_factors[REP_CMPLT_WHITE_LINE_PANEL];

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end());

	allVersionsIterator vit(rec, iteratorSignalIds);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		int time_limit = (time_points.size()) ? time_points[iver] : -1;

		// Get Signals
		rec.usvs.resize(n_sigs);
		for (size_t i = 0; i < n_sigs; ++i)
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);

		// Get Panels
		vector<int> panel_times;
		vector<vector<float> > panels;
		get_panels(rec.usvs, panel_times, panels, time_limit, WHITE_PNL_LAST);

		// Complete
		vector<int> changed_signals(n_sigs, 0);
		for (size_t i = 0; i < panels.size(); i++) {
			int complete = 1;
			while (complete) {
				complete = 0;

				// WBC = SUM(#s)
				complete += sum_complete(panels[i], WHITE_PNL_WBC, white_panel_nums, orig_res, final_res, conv, changed_signals);

				// White subtypes - 
				for (int j = 0; j<white_panel_nums.size(); j++) {
					int num_idx = white_panel_nums[j];
					int perc_idx = white_panel_precs[j];

					// Perc = 100 * Num/WBC ;
					complete += triplet_complete(panels[i], 100, perc_idx, num_idx, WHITE_PNL_WBC, orig_res, final_res, conv, changed_signals);
				}
			}
		}

		// Update changed signals
		if (update_signals(rec, iver, panels, panel_times, sigs_ids, changed_signals) < 0)
			return -1;

	}

	return 0;
}

// Completion of platelets panels
//.......................................................................................
int RepPanelCompleter::apply_platelets_completer(PidDynamicRec& rec, vector<int>& time_points) {
	vector<int>& sigs_ids = panel_signal_ids[REP_CMPLT_PLATELETS_PANEL];
	int n_sigs = (int)sigs_ids.size();

	vector<float>& orig_res = original_sig_res[REP_CMPLT_PLATELETS_PANEL];
	vector<float>& final_res = final_sig_res[REP_CMPLT_PLATELETS_PANEL];
	vector<float>& conv = sig_conversion_factors[REP_CMPLT_PLATELETS_PANEL];

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end());

	allVersionsIterator vit(rec, iteratorSignalIds);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		int time_limit = (time_points.size()) ? time_points[iver] : -1;

		// Get Signals
		rec.usvs.resize(n_sigs);
		for (size_t i = 0; i < n_sigs; ++i)
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);

		// Get Panels
		vector<int> panel_times;
		vector<vector<float> > panels;
		get_panels(rec.usvs, panel_times, panels, time_limit, PLT_PNL_LAST);

		// Complete
		vector<int> changed_signals(n_sigs, 0);
		for (size_t i = 0; i < panels.size(); i++)
			triplet_complete(panels[i], 100, PLT_PNL_MPV, PLT_PNL_PLT_HCT, PLT_PNL_PLTS, orig_res, final_res, conv, changed_signals);


		// Update changed signals
		if (update_signals(rec, iver, panels, panel_times, sigs_ids, changed_signals) < 0)
			return -1;
	}

	return 0;
}

// Completion of lipids panels
//.......................................................................................
int RepPanelCompleter::apply_lipids_completer(PidDynamicRec& rec, vector<int>& time_points) {

	vector<int>& sigs_ids = panel_signal_ids[REP_CMPLT_LIPIDS_PANEL];
	int n_sigs = (int)sigs_ids.size();

	vector<float>& orig_res = original_sig_res[REP_CMPLT_LIPIDS_PANEL];
	vector<float>& final_res = final_sig_res[REP_CMPLT_LIPIDS_PANEL];
	vector<float>& conv = sig_conversion_factors[REP_CMPLT_LIPIDS_PANEL];

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end());

	allVersionsIterator vit(rec, iteratorSignalIds);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		int time_limit = (time_points.size()) ? time_points[iver] : -1;

		// Get Signals
		rec.usvs.resize(n_sigs);
		for (size_t i = 0; i < n_sigs; ++i)
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);

		// Get Panels
		vector<int> panel_times;
		vector<vector<float> > panels;
		get_panels(rec.usvs, panel_times, panels, time_limit, LIPIDS_PNL_LAST);

		// Tryglicerids -> VLDL
		for (int iPanel = 0; iPanel < panels.size(); iPanel++) {
			if (panels[iPanel][LIPIDS_PNL_TRGS] != missing_val)
				panels[iPanel][LIPIDS_PNL_VLDL] = panels[iPanel][LIPIDS_PNL_TRGS] / 5.0F;
		}

		// Complete
		vector<int> changed_signals(LIPIDS_PNL_LAST, 0);
		for (size_t i = 0; i < panels.size(); i++) {
			int complete = 1;
			while (complete) {
				complete = 0;

				complete += triplet_complete(panels[i], 1, LIPIDS_PNL_HDL_OVER_CHOL, LIPIDS_PNL_HDL, LIPIDS_PNL_CHOL, orig_res, final_res, conv, changed_signals);
				complete += triplet_complete(panels[i], 1, LIPIDS_PNL_CHOL_OVER_HDL, LIPIDS_PNL_CHOL, LIPIDS_PNL_HDL, orig_res, final_res, conv, changed_signals);
				complete += triplet_complete(panels[i], 1, LIPIDS_PNL_HDL_OVER_LDL, LIPIDS_PNL_HDL, LIPIDS_PNL_LDL, orig_res, final_res, conv, changed_signals);
				complete += triplet_complete(panels[i], 1, LIPIDS_PNL_LDL_OVER_HDL, LIPIDS_PNL_LDL, LIPIDS_PNL_HDL, orig_res, final_res, conv, changed_signals);
				complete += triplet_complete(panels[i], 1, LIPIDS_PNL_HDL_OVER_NON_HDL, LIPIDS_PNL_HDL, LIPIDS_PNL_NON_HDL_CHOL, orig_res, final_res, conv, changed_signals);

				complete += reciprocal_complete(panels[i], 1, LIPIDS_PNL_HDL_OVER_LDL, LIPIDS_PNL_LDL_OVER_HDL, orig_res, final_res, conv, changed_signals);
				complete += reciprocal_complete(panels[i], 1, LIPIDS_PNL_HDL_OVER_CHOL, LIPIDS_PNL_CHOL_OVER_HDL, orig_res, final_res, conv, changed_signals);

				complete += sum_complete(panels[i], LIPIDS_PNL_CHOL, chol_types1, orig_res, final_res, conv, changed_signals);
				complete += sum_complete(panels[i], LIPIDS_PNL_CHOL, chol_types2, orig_res, final_res, conv, changed_signals);
			}
		}

		// VLDL -> Tryglicerids
		if (changed_signals[LIPIDS_PNL_VLDL]) {
			changed_signals[LIPIDS_PNL_TRGS] = 1;
			for (int iPanel = 0; iPanel < panels.size(); iPanel++) {
				if (panels[iPanel][LIPIDS_PNL_TRGS] == missing_val && panels[iPanel][LIPIDS_PNL_VLDL] != missing_val)
					panels[iPanel][LIPIDS_PNL_TRGS] = 5.0F * panels[iPanel][LIPIDS_PNL_VLDL];
			}
		}

		// Update changed signals
		if (update_signals(rec, iver, panels, panel_times, sigs_ids, changed_signals) < 0)
			return -1;

	}

	return 0;
}

// Completion of eGFR panels
//.......................................................................................
int RepPanelCompleter::apply_eGFR_completer(PidDynamicRec& rec, vector<int>& time_points) {

	vector<int>& sigs_ids = panel_signal_ids[REP_CMPLT_EGFR_PANEL];
	int n_sigs = (int)sigs_ids.size();

	vector<float>& orig_res = original_sig_res[REP_CMPLT_EGFR_PANEL];
	vector<float>& final_res = final_sig_res[REP_CMPLT_EGFR_PANEL];
	vector<float>& conv = sig_conversion_factors[REP_CMPLT_EGFR_PANEL];

	//  Age & Gender
	int age, bYear, gender;
	if (perpare_for_age_and_gender(rec, age, bYear, gender) < 0)
		return -1;

	int sig_time_unit = rec.my_rep->sigs.Sid2Info[sigs_ids[EGFR_PNL_CRT]].time_unit;

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end());

	allVersionsIterator vit(rec, iteratorSignalIds);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		int time_limit = (time_points.size()) ? time_points[iver] : -1;

		// Get Signals
		rec.usvs.resize(n_sigs);
		for (size_t i = 0; i < n_sigs; ++i)
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);

		// Get Panels
		vector<int> panel_times;
		vector<vector<float> > panels;
		get_panels(rec.usvs, panel_times, panels, time_limit, EGFR_PNL_LAST);

		// Complete
		vector<int> changed_signals(n_sigs, 0);
		float current_age;
		for (size_t i = 0; i < panels.size(); i++) {
			// Age
			if (ageDirectlyGiven)
				current_age = (float)age;
			else
				current_age = (float)(1900 + med_time_converter.convert_date(MedTime::Years, panel_times[i]) - bYear);

			egfr_complete(panels[i], current_age, gender, orig_res, final_res, conv, changed_signals);
		}


		// Update changed signals
		if (update_signals(rec, iver, panels, panel_times, sigs_ids, changed_signals) < 0)
			return -1;
	}

	return 0;
}

// Age/Gender util
//.......................................................................................
int RepPanelCompleter::perpare_for_age_and_gender(PidDynamicRec& rec, int& age, int& bYear, int& gender) {

	rec.uget(genderId, 0);
	if (rec.usv.len == 0) {
		MERR("No Gender given for %d\n", rec.pid);
		return -1;
	}
	gender = (int)(rec.usv.Val(0));

	if (ageDirectlyGiven) {
		rec.uget(ageId, 0);
		if (rec.usv.len == 0) {
			MERR("No AGE given for %d\n", rec.pid);
			return -1;
		}
		age = (int)(rec.usv.Val(0));
	}
	else {
		rec.uget(byearId, 0);
		if (rec.usv.len == 0) {
			MERR("No BYEAR given for %d\n", rec.pid);
			return -1;
		}
		bYear = (int)(rec.usv.Val(0));
	}

	return 0;
}

// Completion of BMI panels
//.......................................................................................
int RepPanelCompleter::apply_BMI_completer(PidDynamicRec& rec, vector<int>& time_points) {

	vector<int>& sigs_ids = panel_signal_ids[REP_CMPLT_BMI_PANEL];
	int n_sigs = (int)sigs_ids.size();

	vector<float>& orig_res = original_sig_res[REP_CMPLT_BMI_PANEL];
	vector<float>& final_res = final_sig_res[REP_CMPLT_BMI_PANEL];
	vector<float>& conv = sig_conversion_factors[REP_CMPLT_BMI_PANEL];

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end());

	allVersionsIterator vit(rec, iteratorSignalIds);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		int time_limit = (time_points.size()) ? time_points[iver] : -1;

		// Get Signals
		rec.usvs.resize(n_sigs);
		for (size_t i = 0; i < n_sigs; ++i)
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);

		// Get Panels
		vector<int> panel_times;
		vector<vector<float> > panels;
		get_panels(rec.usvs, panel_times, panels, time_limit, BMI_PNL_LAST);

		// Get Square of height (in meters)
		for (int iPanel = 0; iPanel < panels.size(); iPanel++) {
			if (panels[iPanel][BMI_PNL_HGT] != missing_val)
				panels[iPanel][BMI_PNL_HGT_SQR] = (panels[iPanel][BMI_PNL_HGT] / 100.0F)*(panels[iPanel][BMI_PNL_HGT] / 100.0F);
		}

		// Complete
		vector<int> changed_signals(BMI_PNL_LAST, 0);
		for (size_t i = 0; i < panels.size(); i++)
			triplet_complete(panels[i], 1, BMI_PNL_BMI, BMI_PNL_WGT, BMI_PNL_HGT_SQR, orig_res, final_res, conv, changed_signals);


		// Get height in cm
		if (changed_signals[BMI_PNL_HGT_SQR]) {
			changed_signals[BMI_PNL_HGT] = 1;
			for (int iPanel = 0; iPanel < panels.size(); iPanel++) {
				if (panels[iPanel][BMI_PNL_HGT] == missing_val && panels[iPanel][BMI_PNL_HGT_SQR] != missing_val)
					panels[iPanel][BMI_PNL_HGT] = completer_round(sqrt(panels[iPanel][BMI_PNL_HGT_SQR]) * 100.0F, original_sig_res[REP_CMPLT_BMI_PANEL][BMI_PNL_HGT],
						final_sig_res[REP_CMPLT_BMI_PANEL][BMI_PNL_HGT], sig_conversion_factors[REP_CMPLT_BMI_PANEL][BMI_PNL_HGT]);
			}
		}

		// Update changed signals
		if (update_signals(rec, iver, panels, panel_times, sigs_ids, changed_signals) < 0)
			return -1;
	}

	return 0;
}


// Utilities
// Generate panels from signals
//.......................................................................................

struct candidate_compare {
	bool operator() (const pair<int, int>& lhs, const pair<int, int>& rhs) const {
		return (lhs.first < rhs.first) || (lhs.first == rhs.first && lhs.second < rhs.second);
	}
};

void RepPanelCompleter::get_panels(vector<UniversalSigVec>& usvs, vector<int>& panel_times, vector<vector<float>> &panels, int time_limit, int panel_size) {

	// Prepare
	panels.clear();
	panel_times.clear();

	set<pair<int, int>, candidate_compare> candidates; // Which signal should we insert next ;

	vector<int> idx(usvs.size(), 0);
	for (int i = 0; i < usvs.size(); i++) {
		if (usvs[i].len > 0)
			candidates.insert({ usvs[i].Time(0),i });
	}

	int currentTime = -1;
	while (!candidates.empty()) {
		// Get the next element
		std::set<pair<int, int>, candidate_compare>::iterator it = candidates.begin();

		// New Time Point
		if (it->first != currentTime) {
			currentTime = it->first;
			panel_times.push_back(it->first);
			panels.push_back(vector<float>(panel_size, missing_val));
		}
		int sig_idx = it->second;
		candidates.erase(it);

		// Update panel
		panels.back()[sig_idx] = usvs[sig_idx].Val(idx[sig_idx]++);

		// New candidate
		if (usvs[sig_idx].len > idx[sig_idx]) {
			int time = usvs[sig_idx].Time(idx[sig_idx]);
			if (time_limit == -1 || time <= time_limit)
				candidates.insert({ time ,sig_idx });
		}
	}

	/* Print Panels 
	for (int i = 0; i < panels.size(); i++) {
		cerr << panel_times[i];
		for (int j = 0; j < panels[i].size(); j++)
			cerr << " " << panels[i][j];
			cerr << "\n";
	}
	*/

}

// Complete X = factor*Y/Z ; Y = X*Z/factor ; Z = factor*Y/X
//.......................................................................................
int RepPanelCompleter::triplet_complete(vector<float>& panel, float factor, int x_idx, int y_idx, int z_idx, vector<float>& orig_res, vector<float>& final_res, vector<float>& conv, vector<int>& changed) {

	// Try completing ...
	if (panel[x_idx] == missing_val && panel[y_idx] != missing_val && panel[z_idx] != missing_val && panel[z_idx] != 0.0) {
		panel[x_idx] = factor*panel[y_idx] / panel[z_idx];
		if (x_idx < orig_res.size())
			panel[x_idx] = completer_round(panel[x_idx], orig_res[x_idx], final_res[x_idx], conv[x_idx]);
		changed[x_idx] = 1;
		return 1;
	}
	else if (panel[y_idx] == missing_val && panel[x_idx] != missing_val && panel[z_idx] != missing_val) {
		panel[y_idx] = panel[x_idx] * panel[z_idx] / factor;
		if (y_idx < orig_res.size())
			panel[y_idx] = completer_round(panel[y_idx], orig_res[y_idx], final_res[y_idx], conv[y_idx]);
		changed[y_idx] = 1;
		return 1;
	}
	else if (panel[z_idx] == missing_val &&panel[y_idx] != missing_val && panel[x_idx] != missing_val && panel[x_idx] != 0) {
		panel[z_idx] = factor*panel[y_idx] / panel[x_idx];
		if (z_idx < orig_res.size())
			panel[z_idx] = completer_round(panel[z_idx], orig_res[z_idx], final_res[z_idx], conv[z_idx]);
		changed[z_idx] = 1;
		return 1;
	}
	else
		return 0;
}

// Complete sum = Sum of summands
//.......................................................................................
int RepPanelCompleter::sum_complete(vector<float>& panel, int sum, vector<int>& summands, vector<float>& orig_res, vector<float>& final_res, vector<float>& conv, vector<int>& changed) {

	int npresent = 0;
	float sumVal = 0.0;
	int missing = -1;
	for (int i = 0; i < summands.size(); i++) {
		if (panel[summands[i]] != missing_val) {
			npresent++;
			sumVal += panel[summands[i]];
		}
		else
			missing = summands[i];
	}

	// Can we complete ?
	if (npresent == summands.size() && panel[sum] == missing_val) {
		panel[sum] = sumVal;
		if (sum < orig_res.size())
			panel[sum] = completer_round(panel[sum], orig_res[sum], final_res[sum], conv[sum]);
		changed[sum] = 1;
		return 1;
	}
	else if (npresent == summands.size() - 1 && panel[sum] != missing_val) {
		float val = panel[sum] - sumVal;
		if (val >= 0) {
			panel[missing] = val;
			if (missing < orig_res.size())
				panel[missing] = completer_round(panel[missing], orig_res[missing], final_res[missing], conv[missing]);
			changed[missing] = 1;
			return 1;
		}
	}

	return 0;
}

// Complete x = factor/y
//.......................................................................................
int RepPanelCompleter::reciprocal_complete(vector<float>& panel, float factor, int x_idx, int y_idx, vector<float>& orig_res, vector<float>& final_res, vector<float>& conv, vector<int>& changed) {

	// Can We complete ?
	if (panel[x_idx] == missing_val && panel[y_idx] != missing_val && panel[y_idx] != 0.0) {
		panel[x_idx] = factor / panel[y_idx];
		if (x_idx < orig_res.size())
			panel[x_idx] = completer_round(panel[x_idx], orig_res[x_idx], final_res[x_idx], conv[x_idx]);
		changed[x_idx] = 1;
		return 1;
	}
	else if (panel[y_idx] == missing_val && panel[x_idx] != missing_val && panel[x_idx] != 0.0) {
		panel[y_idx] = factor / panel[x_idx];
		if (y_idx < orig_res.size())
			panel[y_idx] = completer_round(panel[y_idx], orig_res[y_idx], final_res[y_idx], conv[y_idx]);
		changed[y_idx] = 1;
		return 1;
	}

	return 0;
}

// Complete eGFRs
//.......................................................................................
int RepPanelCompleter::egfr_complete(vector<float>& panel, float age, int gender, vector<float>& orig_res, vector<float>& final_res, vector<float>& conv, vector<int>& changed) {

	int complete = 0;
	float egfr;
	if (panel[EGFR_PNL_CRT] != missing_val && panel[EGFR_PNL_CRT] != 0.0) {
		egfr = completer_round(get_eGFR_CKD_EPI(age, panel[EGFR_PNL_CRT], gender), orig_res[EGFR_PNL_CKD_EPI], final_res[EGFR_PNL_CKD_EPI], conv[EGFR_PNL_CKD_EPI]);
		if (isfinite(egfr)) {
			panel[EGFR_PNL_CKD_EPI] = egfr;
			changed[EGFR_PNL_CKD_EPI] = 1;
			complete = 1;
		}

		egfr = completer_round(get_eGFR_MDRD(age, panel[EGFR_PNL_CRT], gender), orig_res[EGFR_PNL_MDRD], final_res[EGFR_PNL_MDRD], conv[EGFR_PNL_MDRD]);
		if (isfinite(egfr)) {
			panel[EGFR_PNL_MDRD] = egfr;
			changed[EGFR_PNL_MDRD] = 1;
			complete = 1;
		}
	}

	return complete;
}

// Updating signals in dynamic-rec
//.......................................................................................
int RepPanelCompleter::update_signals(PidDynamicRec& rec, int iver, vector<vector<float>>& panels, vector<int>& panel_times, vector<int>& sigs_ids, vector<int>& changed) {

	for (int iSig = 0; iSig < sigs_ids.size(); iSig++) {
		if (changed[iSig]) {
			cerr << "Changing " << iSig << endl;
			vector<float> values(panels.size());
			vector<int> times(panels.size());

			int trueSize = 0;
			for (int iPanel = 0; iPanel < panels.size(); iPanel++) {
				if (panels[iPanel][iSig] != missing_val) {
					values[trueSize] = panels[iPanel][iSig];
					times[trueSize] = panel_times[iPanel];
					trueSize++;
				}
			}


			if (rec.set_version_universal_data(sigs_ids[iSig], iver, &(times[0]), &(values[0]), trueSize) < 0)
				return -1;
		}
	}

	return 0;
}

// Read Signals metadata - extract resolutions and conversions from a csv
//.......................................................................................
void RepPanelCompleter::read_metadata() {

	if (metadata_file.empty())
		MTHROW_AND_ERR("No metadata file given\n");

	// Open
	ifstream infile;
	infile.open(metadata_file.c_str(), ifstream::in);
	if (!infile.is_open())
		MTHROW_AND_ERR("Cannot open %s for reading\n", metadata_file.c_str());

	// Read
	int header = 1;
	string thisLine;
	map<string, int> columns;
	map<string, float> all_original_res, all_final_res, all_conversion_factors;
	vector<string> required = { "Name","FinalFactor","OrigResolution","FinalResolution" };

	while (!infile.eof()) {
		getline(infile, thisLine);
		if (thisLine.empty() || thisLine.substr(0, 1) == "#")
			continue;

		vector<string> fields;
		boost::split(fields, thisLine, boost::is_any_of(","));

		if (header == 1) {
			for (int iCol = 0; iCol < fields.size(); iCol++)
				columns[fields[iCol]] = iCol;

			for (string& req : required) {
				if (columns.find(req) == columns.end())
					MTHROW_AND_ERR("Cannot find %s in meta-data file \'%s\'", req.c_str(), metadata_file.c_str());
			}

			header = 0;
		}
		else {
			string sigName = fields[columns["Name"]];
			all_original_res[sigName] = stof(fields[columns["OrigResolution"]]);
			all_final_res[sigName] = stof(fields[columns["FinalResolution"]]);
			all_conversion_factors[sigName] = stof(fields[columns["FinalFactor"]]);
		}
	}

	// Fill in
	original_sig_res.resize(panel_signal_names.size());
	final_sig_res.resize(panel_signal_names.size());
	sig_conversion_factors.resize(panel_signal_names.size());

	for (int iPanel = 0; iPanel < panel_signal_names.size(); iPanel++) {
		original_sig_res[iPanel].resize(panel_signal_names[iPanel].size());
		final_sig_res[iPanel].resize(panel_signal_names[iPanel].size());
		sig_conversion_factors[iPanel].resize(panel_signal_names[iPanel].size());

		for (int iSig = 0; iSig < panel_signal_names[iPanel].size(); iSig++) {

			if (all_original_res.find(panel_signal_names[iPanel][iSig]) == all_original_res.end())
				MTHROW_AND_ERR("Cannot find metadata for signal %s\n", panel_signal_names[iPanel][iSig].c_str());

			original_sig_res[iPanel][iSig] = all_original_res[panel_signal_names[iPanel][iSig]];
			final_sig_res[iPanel][iSig] = all_final_res[panel_signal_names[iPanel][iSig]];
			sig_conversion_factors[iPanel][iSig] = all_conversion_factors[panel_signal_names[iPanel][iSig]];
		}
	}


	infile.close();
}
