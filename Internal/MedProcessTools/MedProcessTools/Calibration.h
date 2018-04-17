#pragma once
#include <vector>
#include "MedFeat/MedFeat/MedOutcome.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "Logger/Logger/Logger.h"

using namespace std;

namespace medial {
	/*!
	*  \brief process namespace
	*/
	namespace calibration {
		struct calibration_entry {
			int bin;
			float min_pred, max_pred;
			int cnt_cases, cnt_controls;
			float mean_pred, mean_outcome;
			float cumul_pct;
			vector<int> controls_per_time_slot;
			vector<int> cases_per_time_slot;
			float kaplan_meier;
		};
		double calc_kaplan_meier(vector<int> controls_per_time_slot, vector<int> cases_per_time_slot);
		void do_calibration(string binning_method, vector<MedSample>& all_samples, int bins_num,
			int pos_sample_min_days_before_case, int pos_sample_max_days_before_case, int km_time_resolution_in_days,
			vector<calibration_entry>& res);
		void smooth_calibration_entries(const vector<calibration_entry>& cals, int min_cases_for_calibration_pct, vector<calibration_entry>& smooth_cals);
		void write_calibration_table(vector<calibration_entry>& cals, const string & calibration_table_file);
		void read_calibration_table(const string& fname, vector<calibration_entry>& cals);
		calibration_entry calibrate_pred(vector<calibration_entry>& cals, float pred);
	}
}