#pragma once
#include <vector>
#include <Logger/Logger/Logger.h>
#include <MedProcessTools/MedProcessTools/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedFeat/MedFeat/MedOutcome.h>

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

enum CaliberationTypes {
	probabilty_time_window = 0, ///< "time_window" bining, but also doing time window\kaplan meir
	probabilty_binning = 1, ///< "binning" - binning
	probabilty_platt_scale = 2 ///< "platt_scale" - platt scale method - sigmoid on pred score and optimize factors
};

extern unordered_map<int, string> caliberation_method_to_name;
static CaliberationTypes cliberation_name_to_type(const string& caliberation_name);

class Calibrator : public PostProcessor {
public:
	CaliberationTypes caliberation_type = probabilty_time_window;

	string estimator_type = "kaplan_meier";
	string binning_method = "equal_num_of_samples_per_bin";
	int bins_num = 1000;
	int pos_sample_min_days_before_case = 0;
	int pos_sample_max_days_before_case = 360;
	int km_time_resolution_in_days = 1;
	int min_cases_for_calibration_smoothing_pct = 10;
	int do_calibration_smoothing = 1;

	int min_preds_in_bin = 100; ///< minimal number of obseravtion to create bin
	float min_score_res = 0; ///< score resulotion value to round to and merge similar
	float min_prob_res = 0; ///< final probality resulotion value to round to and merge similar
	bool fix_pred_order = false; ///< If true will not allow higher scores to have lower probabilites
	int poly_rank = 1; ///< Only in platt_scale - the polynon rank for optimizing sigmoid of prob

	vector<calibration_entry> cals; ///< for "time_window"
	vector<float> min_range, max_range, map_prob; ///< for "binning"
	vector<double> platt_params; ///< for "platt_scale"

	virtual int init(map<string, string>& mapper);
	virtual int Learn(const MedSamples& samples);
	virtual int Learn(const vector <MedSample>& samples);
	virtual int Apply(MedSamples& samples);
	virtual int Apply(vector <MedSample>& samples);

	calibration_entry calibrate_pred(float pred);

	void write_calibration_table(const string & calibration_table_file);
	void read_calibration_table(const string& fname);

	ADD_SERIALIZATION_FUNCS(caliberation_type, estimator_type, binning_method, bins_num, pos_sample_min_days_before_case, pos_sample_max_days_before_case,
		km_time_resolution_in_days, min_cases_for_calibration_smoothing_pct, do_calibration_smoothing,
		min_preds_in_bin, min_score_res, min_prob_res, fix_pred_order, poly_rank,
		cals, min_range, max_range, map_prob, platt_params)

protected:
	double calc_kaplan_meier(vector<int> controls_per_time_slot, vector<int> cases_per_time_slot);
	void smooth_calibration_entries(const vector<calibration_entry>& cals, vector<calibration_entry>& smooth_cals);

private:
	int learn_time_window(const vector<MedSample>& orig_samples);
	int apply_time_window(MedSamples& samples);
	int apply_time_window(vector<MedSample>& samples);
	void write_calibration_time_window(const string & calibration_table_file);
	void read_calibration_time_window(const string& fname);
};
