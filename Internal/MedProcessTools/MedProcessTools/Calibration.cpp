#include "Calibration.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedAlgo/MedAlgo/MedLinearModel.h"

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL

unordered_map<int, string> caliberation_method_to_name = {
	{probabilty_time_window, "time_window"},
	{ probabilty_binning, "binning" },
	{ probabilty_platt_scale, "platt_scale" }
};

CaliberationTypes cliberation_name_to_type(const string& caliberation_name) {
	for (auto it = caliberation_method_to_name.begin(); it != caliberation_method_to_name.end(); ++it)
		if (it->second == caliberation_name)
			return CaliberationTypes(it->first);
	string opts = medial::io::get_list_op(caliberation_method_to_name);
	MTHROW_AND_ERR("unknown caliberation_name \"%s\"\nOptions Are:%s\n",
		caliberation_name.c_str(), opts.c_str());
}

int Calibrator::init(map<string, string>& mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		//! [Calibrator::init]
		if (field == "caliberation_type") caliberation_type = cliberation_name_to_type(entry.second);
		else if (field == "estimator_type") estimator_type = entry.second;
		else if (field == "binning_method") binning_method = entry.second;
		else if (field == "bins_num") bins_num = stoi(entry.second);
		else if (field == "pos_sample_min_days_before_case") pos_sample_min_days_before_case = stoi(entry.second);
		else if (field == "pos_sample_max_days_before_case") pos_sample_max_days_before_case = stoi(entry.second);
		else if (field == "km_time_resolution_in_days") km_time_resolution_in_days = stoi(entry.second);
		else if (field == "do_calibration_smoothing") do_calibration_smoothing = stoi(entry.second);
		else if (field == "min_cases_for_calibration_smoothing_pct") min_cases_for_calibration_smoothing_pct = stoi(entry.second);
		else if (field == "min_preds_in_bin") min_preds_in_bin = stoi(entry.second);
		else if (field == "min_score_res") min_score_res = stof(entry.second);
		else if (field == "min_prob_res") min_prob_res = stof(entry.second);
		else if (field == "fix_pred_order") fix_pred_order = stoi(entry.second) > 0;
		else if (field == "poly_rank") poly_rank = stoi(entry.second);
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

void collect_preds_labels(const MedSamples &samples,
	vector<float> &preds, vector<float> &labels) {
	for (size_t i = 0; i < samples.idSamples.size(); ++i)
	{
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j) {
			if (samples.idSamples[i].samples[j].prediction.empty())
				MTHROW_AND_ERR("no prediction for samples %d, %d\n", (int)i, (int)j);
			preds[i] = samples.idSamples[i].samples[j].prediction[0];
			labels[i] = samples.idSamples[i].samples[j].outcome;
		}
	}
}

void collect_preds_labels(const vector<MedSample>& orig_samples,
	vector<float> &preds, vector<float> &labels) {
	preds.resize(orig_samples.size());
	labels.resize(orig_samples.size());
	for (size_t i = 0; i < orig_samples.size(); ++i)
	{
		if (orig_samples[i].prediction.empty())
			MTHROW_AND_ERR("no prediction for samples %d\n", (int)i);
		preds[i] = orig_samples[i].prediction[0];
		labels[i] = orig_samples[i].outcome;
	}
}

void apply_binned_prob(const vector<float> &preds, const vector<float> &min_range,
	const vector<float> &max_range, const vector<float> &map_prob, vector<float> &probs) {
	probs.resize(preds.size());

	for (size_t i = 0; i < probs.size(); ++i)
	{
		//search for right range:
		int pos = 0;
		while (pos < map_prob.size() &&
			!((preds[i] > min_range[pos] || pos == map_prob.size() - 1) && (preds[i] <= max_range[pos] || pos == 0)))
			++pos;
		probs[i] = map_prob[pos];
	}
}

template<class T, class L> void apply_platt_scale(const vector<T> &preds, const vector<double> &params, vector<L> &converted) {
	converted.resize((int)preds.size());
	for (size_t i = 0; i < converted.size(); ++i)
	{
		double val = params[0];
		for (size_t k = 1; k < params.size(); ++k)
			val += params[k] * pow(double(preds[i]), double(k));
		val = 1 / (1 + exp(val));//Platt Scale technique for probabilty calibaration
		converted[i] = (L)val;
	}
}
template void apply_platt_scale<double, double>(const vector<double> &preds, const vector<double> &params, vector<double> &converted);
template void apply_platt_scale<double, float>(const vector<double> &preds, const vector<double> &params, vector<float> &converted);
template void apply_platt_scale<float, double>(const vector<float> &preds, const vector<double> &params, vector<double> &converted);
template void apply_platt_scale<float, float>(const vector<float> &preds, const vector<double> &params, vector<float> &converted);

int Calibrator::apply_time_window(MedSamples& samples) {
	int do_km;
	if (estimator_type == "kaplan_meier") {
		do_km = 1;
		MLOG("calibrating [%d] samples using kaplan_meier estimator\n",samples.nSamples());
	}
	else if (estimator_type == "mean_cases") {
		do_km = 0;
		MLOG("calibrating [%d] samples using mean_cases estimator\n", samples.nSamples());
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
	return 0;
}

int Calibrator::apply_time_window(vector<MedSample>& samples) {
	int do_km;
	if (estimator_type == "kaplan_meier") {
		do_km = 1;
		MLOG("calibrating [%d] samples using kaplan_meier estimator\n",samples.size());
	}
	else if (estimator_type == "mean_cases") {
		do_km = 0;
		MLOG("calibrating [%d] samples using mean_cases estimator\n", samples.size());
	}
	else MTHROW_AND_ERR("unknown estimator type [%s]", estimator_type.c_str());

	for (auto& s : samples) {
		calibration_entry best = calibrate_pred(s.prediction[0]);
		if (do_km)
			s.prediction[0] = best.kaplan_meier;
		else
			s.prediction[0] = best.mean_outcome;
	}

	return 0;
}


void write_to_predicition(MedSamples& samples, vector<float> &probs) {
	int idx = 0;
	for (auto &pat : samples.idSamples) {
		for (auto& s : pat.samples) {
			s.prediction.resize(1);
			s.prediction[0] = probs[idx];
			++idx;
		}
	}
}

void write_to_predicition(vector<MedSample>& samples, vector<float> &probs) {
	int idx = 0;
	for (auto& s : samples) {
		s.prediction.resize(1);
		s.prediction[0] = probs[idx];
		++idx;
	}
}


int Calibrator::Apply(MedSamples& samples) {
	vector<float> preds, labels, probs;
	switch (caliberation_type)
	{
	case CaliberationTypes::probabilty_time_window:
		return apply_time_window(samples);
		break;
	case CaliberationTypes::probabilty_binning:
		collect_preds_labels(samples, preds, labels);
		apply_binned_prob(preds, min_range, max_range, map_prob, probs);
		write_to_predicition(samples, probs);

		break;
	case CaliberationTypes::probabilty_platt_scale:
		collect_preds_labels(samples, preds, labels);
		apply_platt_scale(preds, platt_params, probs);
		write_to_predicition(samples, probs);
		break;
	default:
		MTHROW_AND_ERR("Unsupported implementation for learning caliberation method %s\n",
			caliberation_method_to_name[caliberation_type].c_str());
	}
	return 0;
}

int Calibrator::Apply(vector <MedSample>& samples) {
	vector<float> preds, labels, probs;
	switch (caliberation_type)
	{
	case CaliberationTypes::probabilty_time_window:
		return apply_time_window(samples);
		break;
	case CaliberationTypes::probabilty_binning:
		collect_preds_labels(samples, preds, labels);
		apply_binned_prob(preds, min_range, max_range, map_prob, probs);
		write_to_predicition(samples, probs);

		break;
	case CaliberationTypes::probabilty_platt_scale:
		collect_preds_labels(samples, preds, labels);
		apply_platt_scale(preds, platt_params, probs);
		write_to_predicition(samples, probs);
		break;
	default:
		MTHROW_AND_ERR("Unsupported implementation for learning caliberation method %s\n",
			caliberation_method_to_name[caliberation_type].c_str());
	}
	return 0;
}

int Calibrator::Learn(const MedSamples& orig_samples) {
	vector<MedSample> samples;
	for (const auto& pat : orig_samples.idSamples)
		for (const auto& s : pat.samples)
			samples.push_back(s);
	Learn(samples);
	return 0;
}

int Calibrator::learn_time_window(const vector<MedSample>& orig_samples) {

	int do_km;
	if (estimator_type == "kaplan_meier")
		do_km = 1;
	else if (estimator_type == "mean_cases")
		do_km = 0;
	else MTHROW_AND_ERR("unknown estimator type [%s]", estimator_type.c_str());

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
		max_samples_per_bin = max((int)samples.size() / bins_num, 10);
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
			(bin < bins_num)
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
		ce.cumul_pct = 1.0f * (cumul_cnt + ((cnt_controls[i] + cnt_cases[i]) / 2)) / (float)samples.size();
		ce.controls_per_time_slot = bin_controls_per_time_slot[i];
		ce.cases_per_time_slot = bin_cases_per_time_slot[i];
		if (do_km) {
			ce.kaplan_meier = (float)calc_kaplan_meier(bin_controls_per_time_slot[i], bin_cases_per_time_slot[i]);
			ce.mean_outcome = 0.0;
		}
		else {
			ce.kaplan_meier = 0.0;
			ce.mean_outcome = 1.0F * cnt_cases[i] / (cnt_controls[i] + cnt_cases[i]);
		}
		cumul_cnt += ce.cnt_controls + ce.cnt_cases;
		cals.push_back(ce);
	}
	if (do_calibration_smoothing) {
		vector<calibration_entry> smooth_cals;
		smooth_calibration_entries(cals, smooth_cals);
		cals = smooth_cals;
	}
	return 0;
}

void learn_binned_probs(vector<float> &x, const vector<float> &y,
	int min_bucket_size, float min_score_jump, float min_prob_jump, bool fix_prob_order,
	vector<float> &min_range, vector<float> &max_range, vector<float> &map_prob) {
	unordered_map<float, vector<int>> score_to_indexes;
	vector<float> unique_scores;
	for (size_t i = 0; i < x.size(); ++i)
	{
		if (score_to_indexes.find(x[i]) == score_to_indexes.end())
			unique_scores.push_back(x[i]);
		score_to_indexes[x[i]].push_back((int)i);
	}
	sort(unique_scores.begin(), unique_scores.end());
	int sz = (int)unique_scores.size();

	float curr_max = (float)INT32_MAX; //unbounded
	float curr_min = curr_max;
	int pred_sum = 0;
	int curr_cnt = 0;
	vector<int> bin_cnts;
	for (int i = sz - 1; i >= 0; --i)
	{
		//update values curr_cnt, pred_avg
		for (int ind : score_to_indexes[unique_scores[i]])
			pred_sum += int(y[ind] > 0);
		curr_cnt += (int)score_to_indexes[unique_scores[i]].size();

		if (curr_cnt > min_bucket_size && curr_max - unique_scores[i] > min_score_jump) {
			//flush buffer
			curr_min = unique_scores[i];
			max_range.push_back(curr_max);
			min_range.push_back(curr_min);
			map_prob.push_back(float(double(pred_sum) / curr_cnt));
			bin_cnts.push_back(curr_cnt);

			//init new buffer:
			curr_cnt = 0;
			curr_max = unique_scores[i];
			pred_sum = 0;
		}
	}
	if (curr_cnt > 0) {
		//flush last buffer
		curr_min = (float)INT32_MIN;
		max_range.push_back(curr_max);
		min_range.push_back(curr_min);
		map_prob.push_back(float(double(pred_sum) / curr_cnt));
		bin_cnts.push_back(curr_cnt);
	}

	//unite similar prob bins:
	vector<int> ind_to_unite;
	for (int i = (int)map_prob.size() - 1; i >= 1; --i)
		if (abs(map_prob[i] - map_prob[i - 1]) < min_prob_jump ||
			(fix_prob_order && map_prob[i] > map_prob[i - 1])) { //unite bins:
			ind_to_unite.push_back(i);
			int new_count = bin_cnts[i] + bin_cnts[i - 1];
			float new_prob = (map_prob[i] * bin_cnts[i] + map_prob[i - 1] * bin_cnts[i - 1]) / new_count;
			float max_th = max_range[i - 1];
			float min_th = min_range[i];
			min_range[i - 1] = min_th;
			max_range[i - 1] = max_th;
			map_prob[i - 1] = new_prob;
			bin_cnts[i - 1] = new_count;
		}

	//unite from end to start:
	for (int i = 0; i < ind_to_unite.size(); ++i)
	{
		int unite_index = ind_to_unite[i];
		//delete old records:
		min_range.erase(min_range.begin() + unite_index);
		max_range.erase(max_range.begin() + unite_index);
		map_prob.erase(map_prob.begin() + unite_index);
		bin_cnts.erase(bin_cnts.begin() + unite_index);
	}

	MLOG("Created %d bins for mapping prediction scores to probabilities\n", map_prob.size());
	for (size_t i = 0; i < map_prob.size(); ++i)
		MLOG_D("Range: [%2.4f, %2.4f] => %2.4f | %1.2f%%(%d / %d)\n", 
			min_range[i], max_range[i], map_prob[i],
			100 * double(bin_cnts[i]) / y.size(), bin_cnts[i], (int)y.size());
}

void learn_platt_scale(vector<float> x, vector<float> &y,
	int poly_rank, vector<double> &params, int min_bucket_size, float min_score_jump
	, float min_prob_jump, bool fix_pred_order) {
	vector<float> min_range, max_range, map_prob;

	learn_binned_probs(x, y, min_bucket_size, min_score_jump, min_prob_jump, fix_pred_order,
		min_range, max_range, map_prob);

	vector<float> probs;
	apply_binned_prob(x, min_range, max_range, map_prob, probs);
	//probs is the new Y - lets learn A, B:
	MedLinearModel lm(poly_rank); //B is param[0], A is param[1]

	lm.loss_function = [](const vector<double> &prds, const vector<float> &y) {
		double res = 0;
		//L2 on 1 / (1 + exp(A*score + B)) vs Y. prds[i] = A*score+B: 1 / (1 + exp(prds))
		for (size_t i = 0; i < y.size(); ++i)
		{
			double conv_prob = 1 / (1 + exp(prds[i]));
			res += (conv_prob - y[i]) * (conv_prob - y[i]);
		}
		res /= y.size();
		res = sqrt(res);
		return res;
	};
	lm.loss_function_step = [](const vector<double> &prds, const vector<float> &y, const vector<double> &params) {
		double res = 0;
		double reg_coef = 0;
		//L2 on 1 / (1 + exp(A*score + B)) vs Y. prds[i] = A*score+B: 1 / (1 + exp(prds))
		for (size_t i = 0; i < y.size(); ++i)
		{
			double conv_prob = 1 / (1 + exp(prds[i]));
			res += (conv_prob - y[i]) * (conv_prob - y[i]);
		}
		res /= y.size();
		res = sqrt(res);
		//Reg A,B:
		if (reg_coef > 0)
			res += reg_coef * sqrt((params[0] * params[0] + params[1] * params[1]) / 2);
		return res;
	};
	lm.block_num = float(10 * sqrt(poly_rank + 1));
	lm.sample_count = 1000;
	lm.tot_steps = 500000;
	lm.learning_rate = 3 * 1e-1;

	vector<float> poly_preds_params(x.size() * poly_rank);
	for (size_t j = 0; j < poly_rank; ++j)
		for (size_t i = 0; i < x.size(); ++i)
			poly_preds_params[i * poly_rank + j] = (float)pow(x[i], j + 1);

	lm.learn(poly_preds_params, probs, (int)x.size(), poly_rank);
	vector<float> factors(poly_rank), mean_shifts(poly_rank);
	lm.get_normalization(mean_shifts, factors);

	//put normalizations inside params:
	params.resize(poly_rank + 1);
	params[0] = lm.model_params[0];
	for (size_t i = 1; i < params.size(); ++i) {
		params[i] = lm.model_params[i] / factors[i - 1];
		params[0] -= lm.model_params[i] * mean_shifts[i - 1] / factors[i - 1];
	}

	vector<double> converted((int)x.size()), prior_score((int)x.size());
	apply_platt_scale(x, params, converted);
	int tot_pos = 0;
	for (size_t i = 0; i < y.size(); ++i)
		tot_pos += int(y[i] > 0);
	for (size_t i = 0; i < converted.size(); ++i)
		prior_score[i] = double(tot_pos) / y.size();

	double loss_model = _linear_loss_target_rmse(converted, probs);
	double loss_prior = _linear_loss_target_rmse(prior_score, probs);

	MLOG("Platt Scale prior=%2.5f. loss_model=%2.5f, loss_prior=%2.5f\n",
		double(tot_pos) / y.size(), loss_model, loss_prior);
}

int Calibrator::Learn(const vector<MedSample>& orig_samples) {
	vector<float> preds, labels;
	switch (caliberation_type)
	{
	case CaliberationTypes::probabilty_time_window:
		learn_time_window(orig_samples);
		break;
	case CaliberationTypes::probabilty_binning:
		collect_preds_labels(orig_samples, preds, labels);
		learn_binned_probs(preds, labels, min_preds_in_bin,
			min_score_res, min_prob_res, fix_pred_order, min_range, max_range, map_prob);
		break;
	case CaliberationTypes::probabilty_platt_scale:
		collect_preds_labels(orig_samples, preds, labels);
		learn_platt_scale(preds, labels, poly_rank, platt_params, min_preds_in_bin,
			min_score_res, min_prob_res, fix_pred_order);
		break;
	default:
		MTHROW_AND_ERR("Unsupported implementation for learning caliberation method %s\n",
			caliberation_method_to_name[caliberation_type].c_str());
	}
	return 0;
}

void Calibrator::write_calibration_time_window(const string & calibration_table_file) {
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

void Calibrator::write_calibration_table(const string & calibration_table_file) {
	ofstream of;
	switch (caliberation_type)
	{
	case probabilty_time_window:
		write_calibration_time_window(calibration_table_file);
		break;
	case probabilty_binning:
		of.open(calibration_table_file, ios::out);
		if (!of)
			MTHROW_AND_ERR("can't open file %s for write\n", calibration_table_file.c_str());
		of << "bin,min_pred,max_pred,mean_outcome\n";
		for (size_t i = 0; i < min_range.size(); ++i)
			of << i << "," << min_range[i] << "," << max_range[i] << "," << map_prob[i] << "\n";
		of.close();
		MLOG("wrote [%d] bins into [%s]\n", (int)min_range.size(), calibration_table_file.c_str());
		break;
	case probabilty_platt_scale:
		of.open(calibration_table_file, ios::out);
		if (!of)
			MTHROW_AND_ERR("can't open file %s for write\n", calibration_table_file.c_str());
		of << "bin,coeff\n";
		for (size_t i = 0; i < platt_params.size(); ++i)
			of << i << "," << platt_params[i] << "\n";
		of.close();
		MLOG("wrote [%d] bins into [%s]\n", (int)platt_params.size(), calibration_table_file.c_str());
		break;
	default:
		MTHROW_AND_ERR("Unsupported implementation for writing caliberation method %s\n",
			caliberation_method_to_name[caliberation_type].c_str());
	}
}

void Calibrator::read_calibration_time_window(const string& fname) {
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

void Calibrator::read_calibration_table(const string& fname) {
	ifstream f;
	string curr_line;
	vector<string> tokens;
	switch (caliberation_type)
	{
	case probabilty_time_window:
		read_calibration_time_window(fname);
		break;
	case probabilty_binning:
		f.open(fname, ios::in);
		if (!f)
			MTHROW_AND_ERR("can't open file %s for write\n", fname.c_str());
		while (getline(f, curr_line)) {
			boost::split(tokens, curr_line, boost::is_any_of(","));
			if (tokens.size() != 4)
				MTHROW_AND_ERR("Bad format in line:\n%s\nexpected 4 tokens in probabilty bining\n",
					curr_line.c_str());
			if (tokens[0] == "bin")
				continue; //skip header
			min_range.push_back(stof(tokens[1]));
			max_range.push_back(stof(tokens[2]));
			map_prob.push_back(stof(tokens[3]));
		}
		f.close();
		MLOG("read [%d] bins into [%s]\n", (int)min_range.size(), fname.c_str());
		break;
	case probabilty_platt_scale:
		f.open(fname, ios::in);
		if (!f)
			MTHROW_AND_ERR("can't open file %s for write\n", fname.c_str());
		while (getline(f, curr_line)) {
			boost::split(tokens, curr_line, boost::is_any_of(","));
			if (tokens.size() != 2)
				MTHROW_AND_ERR("Bad format in line:\n%s\nexpected 2 tokens in probabilty platt scale\n",
					curr_line.c_str());
			if (tokens[0] == "bin")
				continue; //skip header
			platt_params.push_back(stod(tokens[1]));
		}
		f.close();
		MLOG("read [%d] bins into [%s]\n", (int)platt_params.size(), fname.c_str());
		break;
	default:
		MTHROW_AND_ERR("Unsupported implementation for reading caliberation method %s\n",
			caliberation_method_to_name[caliberation_type].c_str());
	}
}

calibration_entry Calibrator::calibrate_pred(float pred) {
	if (caliberation_type == probabilty_time_window) {
		int start = 0;
		for (int i = 0; i < cals.size(); i++) {
			if (pred >= cals[i].min_pred)
				start = i;
		}
		return cals[start];
	}
	else {
		MedSample smp;
		smp.prediction = { pred };
		smp.outcome = 0; //doesn't matter
		smp.id = 1; //doesn't matter
		MedIdSamples smp_id(1);
		smp_id.samples.push_back(smp);
		MedSamples samples;
		samples.idSamples = { smp_id };
		Apply(samples);
		calibration_entry res;
		res.mean_outcome = samples.idSamples[0].samples[0].prediction[0];
		res.min_pred = pred;
		res.max_pred = pred;
		res.mean_pred = pred;
		res.bin = 0;
		res.cnt_cases = 0;
		res.cnt_controls = 0;
		return res;
	}
}


