#pragma once
#include <vector>
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedFeat/MedFeat/MedOutcome.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "Logger/Logger/Logger.h"

using namespace std;

class PostProcessor : public SerializableObject {
};

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

class Calibrator : public PostProcessor {
public:
	string estimator_type = "kaplan_meier";
	string binning_method = "equal_num_of_samples_per_bin";	
	int bins_num = 1000;
	int pos_sample_min_days_before_case = 0;
	int pos_sample_max_days_before_case = 360;
	int km_time_resolution_in_days = 1;
	int min_cases_for_calibration_smoothing_pct = 10;
	int do_calibration_smoothing = 1;

	vector<calibration_entry> cals;

	virtual int init(map<string, string>& mapper);
	virtual int Learn(const MedSamples& samples);
	virtual int Learn(const vector <MedSample>& samples);
	virtual int Apply(MedSamples& samples);

	calibration_entry calibrate_pred(float pred);

	void write_calibration_table(const string & calibration_table_file);
	void read_calibration_table(const string& fname);

	ADD_SERIALIZATION_FUNCS(estimator_type, binning_method, bins_num, pos_sample_min_days_before_case, pos_sample_max_days_before_case,
		km_time_resolution_in_days, min_cases_for_calibration_smoothing_pct, do_calibration_smoothing, cals)

protected:
	double calc_kaplan_meier(vector<int> controls_per_time_slot, vector<int> cases_per_time_slot);
	void smooth_calibration_entries(const vector<calibration_entry>& cals, vector<calibration_entry>& smooth_cals);
};
