#ifndef __CALIBRATION_H__
#define __CALIBRATION_H__

#include <vector>
#include <Logger/Logger/Logger.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <MedFeat/MedFeat/MedOutcome.h>
#include "PostProcessor.h"

using namespace std;

/**
* Calibrator are post processocrs using for recalibration of a model
*/

// Calibration entry is used for transoforming a score into probabibilty
class calibration_entry : public SerializableObject {
public:
	int bin;
	float min_pred, max_pred;
	double cnt_cases, cnt_controls;
	double cnt_cases_no_w, cnt_controls_no_w;
	float mean_pred, mean_outcome;
	float cumul_pct;
	vector<double> controls_per_time_slot;
	vector<double> cases_per_time_slot;
	float kaplan_meier;

	string str() const;

	ADD_CLASS_NAME(calibration_entry)
		ADD_SERIALIZATION_FUNCS(bin, min_pred, max_pred, cnt_cases, cnt_controls, cnt_cases_no_w, cnt_controls_no_w, mean_pred, mean_outcome, controls_per_time_slot, cases_per_time_slot, kaplan_meier)
};
MEDSERIALIZE_SUPPORT(calibration_entry)

enum CalibrationTypes {
	probability_time_window = 0, ///< "time_window" bining, but also doing time window\kaplan meir
	probability_binning = 1, ///< "binning" - binning
	probability_platt_scale = 2, ///< "platt_scale" - platt scale method - sigmoid on pred score and optimize factors
	probability_isotonic = 3 ///< "isotonic_regression" - piece-wise constance isotonic function from score to probability
};

extern unordered_map<int, string> calibration_method_to_name;
static CalibrationTypes clibration_name_to_type(const string& calibration_name);

class Calibrator : public PostProcessor {
public:
	Calibrator() { processor_type = PostProcessorTypes::FTR_POSTPROCESS_CALIBRATOR; }

	CalibrationTypes calibration_type = probability_time_window;

	int time_unit = MedTime::Days;

	string estimator_type = "kaplan_meier";
	string binning_method = "equal_num_of_samples_per_bin";
	int bins_num = 1000;
	int pos_sample_min_time_before_case = 0;
	int pos_sample_max_time_before_case = 360;
	int km_time_resolution = 1;
	int min_cases_for_calibration_smoothing_pct = 10;
	int do_calibration_smoothing = 1;
	int censor_controls = 0; ///< censor controls without long-enough followup even in mean-outcome mode
	string weights_attr_name = "weight"; //weight attr to look for in samples attributes
	int min_control_bins = -1;
	bool use_isotonic = false; ///< If true will use isotonic on time_window

	int min_preds_in_bin = 100; ///< minimal number of obseravtion to create bin
	float min_score_res = 0; ///< score resulotion value to round to and merge similar
	float min_prob_res = 0; ///< final probality resulotion value to round to and merge similar
	bool fix_pred_order = false; ///< If true will not allow higher scores to have lower probabilites
	int poly_rank = 1; ///< Only in platt_scale - the polynon rank for optimizing sigmoid of prob
	double control_weight_down_sample = 1; ///< factor weight for controls when downsampling controls by this factor
	bool verbose = true; ///< If true will print verbose information for calibration

	int n_top_controls = 0; ///< number of controls to add with maximal-score for regularization of isotonic regression
	int n_bottom_cases = 0; ///< number of cases to add with minimal-score for regularization of isotonic regression

	vector<calibration_entry> cals; ///< for "time_window"
	vector<float> min_range, max_range, map_prob; ///< for "binning/isotonic-regression"
	vector<double> platt_params; ///< for "platt_scale"

	/// @snippet Calibration.cpp Calibrator::init
	virtual int init(map<string, string>& mapper);
	virtual int Learn(const MedSamples& samples);
	virtual int Learn(const vector<MedSample>& samples) { return Learn(samples, global_default_time_unit); }
	virtual int Learn(const vector <MedSample>& samples, const int samples_time_unit);
	virtual int Apply(MedSamples& samples);
	virtual int Apply(vector <MedSample>& samples);
	void Apply(const vector<float> &preds, vector<float> &probs) const;
	float Apply(float pred) const;

	void get_input_fields(vector<Effected_Field> &fields) const;
	void get_output_fields(vector<Effected_Field> &fields) const;

	//PostProcessor functions:
	void Learn(const MedFeatures &matrix) {Learn(matrix.samples); }
	void Apply(MedFeatures &matrix);

	calibration_entry calibrate_pred(float pred);
	float calibrate_pred(float pred, int type) const;

	void write_calibration_table(const string & calibration_table_file);
	void read_calibration_table(const string& fname);

	void dprint(const string &pref) const;
	void learn_isotonic_regression(const vector<float> &x, const vector<float> &y, const vector<float> &weights, vector<float> &min_range, vector<float> &max_range, vector<float> &map_prob, int n_top_controls, int n_bottom_cases,
		bool verbose);

	ADD_CLASS_NAME(Calibrator)
	ADD_SERIALIZATION_FUNCS(calibration_type, estimator_type, binning_method, bins_num, time_unit, pos_sample_min_time_before_case, pos_sample_max_time_before_case,
		km_time_resolution, min_cases_for_calibration_smoothing_pct, do_calibration_smoothing, censor_controls,
		min_preds_in_bin, min_score_res, min_prob_res, fix_pred_order, poly_rank, control_weight_down_sample,
		cals, min_range, max_range, map_prob, platt_params, use_isotonic, n_top_controls, n_bottom_cases)

protected:
	double calc_kaplan_meier(vector<double> controls_per_time_slot, vector<double> cases_per_time_slot, double controls_factor);
	void smooth_calibration_entries(const vector<calibration_entry>& cals, vector<calibration_entry>& smooth_cals, double controls_factor);

private:
	int learn_time_window(const vector<MedSample>& orig_samples, const int samples_time_unit);
	int apply_time_window(MedSamples& samples) const;
	int apply_time_window(vector<MedSample>& samples) const;
	void write_calibration_time_window(const string & calibration_table_file);
	void read_calibration_time_window(const string& fname);
};

MEDSERIALIZE_SUPPORT(Calibrator)

#endif