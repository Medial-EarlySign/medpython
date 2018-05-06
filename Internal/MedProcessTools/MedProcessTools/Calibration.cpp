#include "Calibration.h"

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL

int Calibrator::init(map<string, string>& mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		//! [Calibrator::init]
		if (field == "estimator_type") estimator_type = entry.second;
		else if (field == "binning_method") binning_method = entry.second;
		else if (field == "bins_num") bins_num = stoi(entry.second);
		else if (field == "pos_sample_min_days_before_case") pos_sample_min_days_before_case = stoi(entry.second);
		else if (field == "pos_sample_max_days_before_case") pos_sample_max_days_before_case = stoi(entry.second);
		else if (field == "km_time_resolution_in_days") km_time_resolution_in_days = stoi(entry.second);
		else if (field == "do_calibration_smoothing") do_calibration_smoothing = stoi(entry.second);
		else if (field == "min_cases_for_calibration_smoothing_pct") min_cases_for_calibration_smoothing_pct = stoi(entry.second);	
		else MTHROW_AND_ERR("unknown init option [%s] for Calibrator\n", field.c_str());
		//! [Calibrator::init]
	}

	return 0;
}


double Calibrator::calc_kaplan_meier(vector<int> controls_per_time_slot, vector<int> cases_per_time_slot) {
	double prob = 1.0;
	int total_controls_all = 0;
	for (int i = 0; i < controls_per_time_slot.size(); i++)
		total_controls_all += controls_per_time_slot[i];
	//MLOG("size %d total_controls_all %d\n", controls_per_time_slot.size(), total_controls_all);
	for (int i = 0; i < controls_per_time_slot.size(); i++) {
		if (total_controls_all <= 0)
			MTHROW_AND_ERR("calc_kaplan_meier left without controls in time slot [%d]\n", i);
		prob *= (double)total_controls_all / (cases_per_time_slot[i] + total_controls_all);
		total_controls_all -= controls_per_time_slot[i];
	}
	return 1.0 - prob;
}

// expand to neighbor calibration entries, until finding enough cases
void Calibrator::smooth_calibration_entries(const vector<calibration_entry>& cals, vector<calibration_entry>& smooth_cals) {
	smooth_cals.clear();
	int cases = 0;
	for (auto& c : cals)
		cases += c.cnt_cases;
	int min_cases_for_calibration = min_cases_for_calibration_smoothing_pct * cases / 100;
	MLOG("smooth_calibration_entries requiring min_cases_for_calibration_smoothing = [%d * %d / 100 = %d]\n", min_cases_for_calibration_smoothing_pct, cases,
		min_cases_for_calibration);

	for (int s = 0; s < cals.size(); s++) {
		int end = s, start = s;
		int controls = cals[start].cnt_controls;
		int cases = cals[start].cnt_cases;
		vector<int> controls_per_time_slot;
		vector<int> cases_per_time_slot;
		for (size_t j = 0; j < cals[start].controls_per_time_slot.size(); j++)
		{
			controls_per_time_slot.push_back(cals[start].controls_per_time_slot[j]);
			cases_per_time_slot.push_back(cals[start].cases_per_time_slot[j]);
		}

		while (cases < min_cases_for_calibration) {
			if (end == (cals.size() - 1) && start == 0)
				// the entire calibration table holds less than min_cases_for_calibration cases
				break;
			if (end < cals.size() - 1) {
				end++;
				controls += cals[end].cnt_controls;
				cases += cals[end].cnt_cases;
				for (size_t j = 0; j < cals[end].controls_per_time_slot.size(); j++)
				{
					controls_per_time_slot[j] += cals[end].controls_per_time_slot[j];
					cases_per_time_slot[j] += cals[end].cases_per_time_slot[j];
				}
			}
			if (start > 0) {
				start--;
				controls += cals[start].cnt_controls;
				cases += cals[start].cnt_cases;
				for (size_t j = 0; j < cals[start].controls_per_time_slot.size(); j++)
				{
					controls_per_time_slot[j] += cals[start].controls_per_time_slot[j];
					cases_per_time_slot[j] += cals[start].cases_per_time_slot[j];
				}
			}
		}
		calibration_entry res = cals[s];
		res.cnt_controls = controls;
		res.cnt_cases = cases;
		res.mean_outcome = 1.0f * cases / (cases + controls);
		res.cumul_pct = cals[start].cumul_pct;
		res.kaplan_meier = (float)calc_kaplan_meier(controls_per_time_slot, cases_per_time_slot);
		smooth_cals.push_back(res);
	}
}

int Calibrator::Learn(const MedSamples& orig_samples) {
	vector<MedSample> samples;
	for (const auto& pat : orig_samples.idSamples)
		for (const auto& s : pat.samples)
			samples.push_back(s);
	Learn(samples);
	return 0;
}
int Calibrator::Learn(const vector<MedSample>& orig_samples ) {
	cals.clear();
	float min_pred = 100000.0, max_pred = -100000.0;
	int cases = 0;
	set<float> unique_preds;
	vector<MedSample> samples;
	for (auto e : orig_samples) {
		if (unique_preds.find(e.prediction[0]) == unique_preds.end())
			unique_preds.insert(e.prediction[0]);
		if (e.prediction[0] < min_pred)
			min_pred = e.prediction[0];
		if (e.prediction[0] > max_pred)
			max_pred = e.prediction[0];
		int gap = get_day_approximate(e.outcomeTime) - get_day_approximate(e.time);
		if (gap < pos_sample_min_days_before_case)
			// too close to outcome date or censor date (chance for peeking, or even beyond the outcome date)
			continue;
		if (e.outcome >= 1 && gap > pos_sample_max_days_before_case)
			// too far case is considered as control
			e.outcome = 0;
		if (e.outcome >= 1)
			cases++;
		samples.push_back(e);
	}
	std::sort(samples.begin(), samples.end(), comp_sample_pred);
	MLOG("eligible samples [%d] cases [%d]\n", int(samples.size()), cases);
	int max_samples_per_bin = 0;
	int max_cases_per_bin = 0;
	float max_delta_in_bin = 0.0;
	if (binning_method == "unique_score_per_bin") {
		bins_num = (int)unique_preds.size();
		MLOG("unique_score_per_bin, bins_num [%d] \n", bins_num);
	}
	else if (binning_method == "equal_num_of_samples_per_bin") {
		max_samples_per_bin = max((int)samples.size()/ bins_num, 10);
		MLOG("equal_num_of_samples_per_bin bins_num: %d max_samples_per_bin: %d \n",
			bins_num, max_samples_per_bin);
	}
	else if (binning_method == "equal_num_of_cases_per_bin") {
		max_cases_per_bin = max(cases / bins_num, 10);
		MLOG("equal_num_of_cases_per_bin bins_num: %d max_cases_per_bin: %d \n",
			bins_num, max_cases_per_bin);
	}
	else if (binning_method == "equal_score_delta_per_bin") {
		max_delta_in_bin = (max_pred - min_pred) / bins_num;
		MLOG("equal_score_delta_per_bin min_pred: %f max_pred: %f max_delta_in_bin: %f \n",
			min_pred, max_pred, (int)samples.size(), max_delta_in_bin);
	}
	else MTHROW_AND_ERR("unknown binning method [%s]\n", binning_method.c_str());

	vector<int> cnt_cases;
	vector<int> cnt_controls;
	vector<float> bin_max_preds;
	vector<float> bin_min_preds;
	vector<float> bin_sum_preds;
	vector<vector<int>> bin_controls_per_time_slot;
	vector<vector<int>> bin_cases_per_time_slot;
	int km_time_slots = (pos_sample_max_days_before_case - pos_sample_min_days_before_case) / km_time_resolution_in_days;
	if (km_time_slots <= 0)
		km_time_slots = 1;
	MLOG("km_time_slots [%d] \n", km_time_slots);
	for (size_t i = 0; i < (bins_num + 1) * 2; i++)
	{
		cnt_cases.push_back(0);
		cnt_controls.push_back(0);
		bin_sum_preds.push_back(0.0);
		bin_min_preds.push_back(min_pred);
		bin_max_preds.push_back(0.0);
		vector<int> controls_per_time_slot;
		vector<int> cases_per_time_slot;
		for (size_t j = 0; j <= km_time_slots; j++) {
			controls_per_time_slot.push_back(0);
			cases_per_time_slot.push_back(0);
		}
		bin_controls_per_time_slot.push_back(controls_per_time_slot);
		bin_cases_per_time_slot.push_back(cases_per_time_slot);
	}
	int cnt = 0, bin = 1;
	float prev_pred = min_pred;
	for (auto &o : samples) {
		// TODO: use medtime
		int gap = get_day_approximate(o.outcomeTime) - get_day_approximate(o.time);
		if (
			(bin < bins_num )
			&&
			(
			(binning_method == "unique_score_per_bin" && prev_pred != o.prediction[0]) ||
			(binning_method == "equal_num_of_samples_per_bin" && cnt_controls[bin] + cnt_cases[bin] >= max_samples_per_bin) ||
			(binning_method == "equal_score_delta_per_bin" && (o.prediction[0] - bin_min_preds[bin]) > max_delta_in_bin) ||
			(binning_method == "equal_num_of_cases_per_bin" && cnt_cases[bin] >= max_cases_per_bin))
			)
		{
			bin++;
			bin_min_preds[bin] = o.prediction[0];
		}
		int time_slot;
		if (gap > pos_sample_max_days_before_case)
			time_slot = km_time_slots;
		else
			time_slot = (gap - pos_sample_min_days_before_case) / km_time_resolution_in_days;

		if (o.outcome >= 1) {
			cnt_cases[bin] ++;
			bin_cases_per_time_slot[bin][time_slot]++;
		}
		else {
			cnt_controls[bin]++;
			bin_controls_per_time_slot[bin][time_slot]++;
		}

		bin_sum_preds[bin] += o.prediction[0];
		bin_max_preds[bin] = o.prediction[0];
		prev_pred = o.prediction[0];
	}
	int cumul_cnt = 0;
	for (int i = 1; i <= bin; i++)
	{
		calibration_entry ce;
		ce.bin = i;
		ce.min_pred = bin_min_preds[i];
		ce.max_pred = bin_max_preds[i];
		ce.cnt_controls = cnt_controls[i]; ce.cnt_cases = cnt_cases[i];
		ce.mean_pred = 1.0f * bin_sum_preds[i] / (cnt_controls[i] + cnt_cases[i]);
		ce.mean_outcome = 1.0f * cnt_cases[i] / (cnt_controls[i] + cnt_cases[i]);
		ce.cumul_pct = 1.0f * (cumul_cnt + ((cnt_controls[i] + cnt_cases[i]) / 2)) / (float)samples.size();
		ce.controls_per_time_slot = bin_controls_per_time_slot[i];
		ce.cases_per_time_slot = bin_cases_per_time_slot[i];
		ce.kaplan_meier = (float)calc_kaplan_meier(bin_controls_per_time_slot[i], bin_cases_per_time_slot[i]);
		cumul_cnt += ce.cnt_controls + ce.cnt_cases;
		cals.push_back(ce);
	}
	if (do_calibration_smoothing) {
		vector<calibration_entry> smooth_cals;
		smooth_calibration_entries(cals, smooth_cals);
		cals = smooth_cals;
	}
}

void Calibrator::write_calibration_table(const string & calibration_table_file) {
	ofstream of;
	of.open(calibration_table_file, ios::out);
	if (!of) {
		MLOG("can't open file %s for write\n", calibration_table_file.c_str());
		throw exception();
	}
	float total_hosmer = 0.0;
	of << "bin,min_pred,max_pred,cnt_controls,cnt_cases,mean_pred,mean_outcome,cumul_pct,kaplan_meier\n";
	for (calibration_entry& ce : cals)
	{
		of << ce.bin << ","
			<< ce.min_pred << ","
			<< ce.max_pred << ","
			<< ce.cnt_controls << ","
			<< ce.cnt_cases << ","
			<< ce.mean_pred << ","
			<< ce.mean_outcome << ","
			<< ce.cumul_pct << ","
			<< ce.kaplan_meier
			<< endl;
	}

	of.close();
	MLOG("wrote [%d] bins into [%s]\n", cals.size(), calibration_table_file.c_str());
}

void Calibrator::read_calibration_table(const string& fname) {
	ifstream inf(fname);
	if (!inf) {
		MLOG("can't open file %s for read\n", fname.c_str());
		throw exception();
	}
	MLOG("reading from: [%s]\n", fname.c_str());

	string curr_line;
	int cnt = 0;
	while (getline(inf, curr_line)) {
		if (curr_line[curr_line.size() - 1] == '\r')
			curr_line.erase(curr_line.size() - 1);

		vector<string> fields;
		split(fields, curr_line, boost::is_any_of(","));
		if (fields[0] == "bin")
			continue; // header

		calibration_entry ce;
		istringstream ls(curr_line);
		char delim;
		ls >> ce.bin >> delim
			>> ce.min_pred >> delim
			>> ce.max_pred >> delim
			>> ce.cnt_controls >> delim
			>> ce.cnt_cases >> delim
			>> ce.mean_pred >> delim
			>> ce.mean_outcome >> delim
			>> ce.cumul_pct >> delim
			>> ce.kaplan_meier;
		cals.push_back(ce);
	}
	MLOG("read %d entries from [%s]\n", cals.size(), fname.c_str());
	inf.close();
}

calibration_entry Calibrator::calibrate_pred(float pred) {
	int start = 0;
	for (int i = 0; i < cals.size(); i++) {
		if (pred >= cals[i].min_pred)
			start = i;
	}
	return cals[start];
}

int Calibrator::Apply(MedSamples& samples) {
	int do_km;
	if (estimator_type == "kaplan_meier") {
		do_km = 1;
		MLOG("calibrating [%d] samples using kaplan_meier estimator");
	}
	else if (estimator_type == "mean_cases") {
		do_km = 0;
		MLOG("calibrating [%d] samples using mean_cases estimator");
	}
	else MTHROW_AND_ERR("unknown estimator type [%s]", estimator_type.c_str());

	for (auto &pat : samples.idSamples) {
		for (auto& s : pat.samples) {
			calibration_entry best = calibrate_pred(s.prediction[0]);
			if (do_km)
				s.prediction[0] = best.kaplan_meier;
			else
				s.prediction[0] = best.mean_outcome;
		}
	}
}