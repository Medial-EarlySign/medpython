#include "bootstrap.h"
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <random>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <cmath>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
//#define WARN_SKIP_WP

#pragma region Helper Functions
float meanArr(const vector<float> &arr) {
	double res = 0;
	int cnt = 0;
	for (size_t i = 0; i < arr.size(); ++i)
		if (arr[i] != MED_MAT_MISSING_VALUE) {
			res += arr[i];
			++cnt;
		}
	if (cnt > 0)
		res /= cnt;
	else
		res = MED_MAT_MISSING_VALUE;
	return (float)res;
}
float stdArr(const vector<float> &arr, float meanVal) {
	double res = 0;
	int cnt = 0;
	for (size_t i = 0; i < arr.size(); ++i)
		if (arr[i] != MED_MAT_MISSING_VALUE) {
			res += (arr[i] - meanVal) * (arr[i] - meanVal);
			++cnt;
		}
	if (cnt > 0) {
		res /= cnt;
		res = sqrt(res);
	}
	else
		res = MED_MAT_MISSING_VALUE;
	return (float)res;
}
template<class T> string print_obj(T obj, string format) {
	//return to_string((round(num * 1000) / 1000));
	char res[50];
	snprintf(res, sizeof(res), format.c_str(), obj);
	return string(res);
}
string printVec(const vector<float> &v, int from, int to) {
	string res;
	if (from < 0)
		from = 0;
	for (size_t i = from; i < to && i < v.size(); ++i)
	{
		res += to_string(v[i]) + ", ";
	}
	return res;
}
#pragma endregion

map<string, float> booststrap_analyze_cohort(const vector<float> &preds, const vector<float> &y,
	const vector<int> &pids, float sample_ratio, int sample_per_pid, int loopCnt,
	const vector<MeasurementFunctions> &meas_functions, const vector<void *> &function_params,
	ProcessMeasurementParamFunc process_measurments_params,
	const map<string, vector<float>> &additional_info, const vector<float> &y_full,
	const vector<int> &pids_full, FilterCohortFunc cohort_def, void *cohort_params) {
	//for each pid - randomize x sample from all it's tests. do loop_times
	float ci_bound = (float)0.95;
	unordered_map<int, vector<int>> pid_to_inds;
	for (size_t i = 0; i < pids.size(); ++i)
		pid_to_inds[pids[i]].push_back(int(i));

	//initialize measurement params per cohort:
	for (size_t i = 0; i < function_params.size(); ++i)
		if (process_measurments_params != NULL && function_params[i] != NULL)
			process_measurments_params(additional_info, y_full, pids_full, cohort_def, cohort_params
				, function_params[i]);

	random_device rd;
	mt19937 rd_gen(rd());

	//this function called after filter cohort

	map<string, vector<float>> all_measures;
	//save results for all cohort:
	for (size_t k = 0; k < meas_functions.size(); ++k)
	{
		map<string, float> batch_measures;
		batch_measures = meas_functions[k](preds, y, function_params[k]);
		for (auto jt = batch_measures.begin(); jt != batch_measures.end(); ++jt)
			all_measures[jt->first + "_Obs"].push_back(jt->second);
	}

	int cohort_size = int(sample_ratio * pid_to_inds.size());
	int cnt_i = 0;
	vector<int> ind_to_pid((int)pid_to_inds.size());

	int min_pid_start = INT_MAX;
	int max_pid_start = 0;
	for (auto it = pid_to_inds.begin(); it != pid_to_inds.end(); ++it)
	{
		ind_to_pid[cnt_i] = it->first;
		++cnt_i;
		if (it->first < min_pid_start)
			min_pid_start = it->first;
		if (it->first > max_pid_start)
			max_pid_start = it->first;
	}
	int rep_pid_size = max_pid_start - min_pid_start;

	//choose pids:
	uniform_int_distribution<> rand_pids(0, (int)pid_to_inds.size() - 1);
	if (sample_per_pid > 0) {
#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < loopCnt; ++i)
		{
			int curr_ind = 0;
			vector<int> selected_inds(cohort_size * sample_per_pid);
			vector<bool> selected_ind_pid((int)pid_to_inds.size());

			for (size_t k = 0; k < cohort_size; ++k)
			{
				int ind_pid;
				if (sample_ratio < 1) {
					ind_pid = rand_pids(rd_gen);
					while (selected_ind_pid[ind_pid - min_pid_start])
						ind_pid = rand_pids(rd_gen);
					selected_ind_pid[ind_pid - min_pid_start] = true;
				}
				else
					ind_pid = (int)k;

				vector<int> inds = pid_to_inds[ind_to_pid[ind_pid]];
				uniform_int_distribution<> random_num(0, (int)inds.size() - 1);
				for (size_t kk = 0; kk < sample_per_pid; ++kk)
				{
					selected_inds[curr_ind] = inds[random_num(rd_gen)];
					++curr_ind;
				}
			}

			//now calculate measures on cohort of selected indexes for preds, y:
			vector<float> selected_preds(cohort_size * sample_per_pid), selected_y(cohort_size * sample_per_pid);
			for (size_t k = 0; k < selected_inds.size(); ++k)
			{
				selected_preds[k] = preds[selected_inds[k]];
				selected_y[k] = y[selected_inds[k]];
			}

			for (size_t k = 0; k < meas_functions.size(); ++k)
			{
				map<string, float> batch_measures;
				batch_measures = meas_functions[k](selected_preds, selected_y, function_params[k]);
#pragma omp critical
				for (auto jt = batch_measures.begin(); jt != batch_measures.end(); ++jt)
					all_measures[jt->first].push_back(jt->second);
			}
		}
	}
	else { //other sampling - sample pids and take all thier data:
		   //now sample cohort 


#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < loopCnt; ++i)
		{
			vector<bool> selected_ind_pid;
			vector<int> selected_pids(cohort_size);
			for (size_t k = 0; k < cohort_size; ++k)
			{
				int ind_pid = rand_pids(rd_gen);
				while (selected_ind_pid[ind_pid - min_pid_start])
					ind_pid = rand_pids(rd_gen);
				selected_ind_pid[ind_pid - min_pid_start] = true;
				selected_pids[k] = ind_to_pid[ind_pid];
			}

			//create preds, y for all seleceted pids:
			vector<float> selected_preds, selected_y;
			for (size_t k = 0; k < selected_pids.size(); ++k)
			{
				int pid = selected_pids[k];
				vector<int> ind_vec = pid_to_inds[pid];
				for (int ind : ind_vec)
				{
					selected_preds.push_back(preds[ind]);
					selected_y.push_back(y[ind]);
				}
			}

			//calc measures for sample:
			for (size_t k = 0; k < meas_functions.size(); ++k)
			{
				map<string, float> batch_measures;
				batch_measures = meas_functions[k](selected_preds, selected_y, function_params[k]);
#pragma omp critical
				for (auto jt = batch_measures.begin(); jt != batch_measures.end(); ++jt)
					all_measures[jt->first].push_back(jt->second);
			}
		}
	}

	//now calc - mean, std , CI0.95_lower, CI0.95_upper for each measurement in all exp
	map<string, float> all_final_measures;
	for (auto it = all_measures.begin(); it != all_measures.end(); ++it)
	{
		vector<float> measures = it->second;
		sort(measures.begin(), measures.end());
		float meanVal = meanArr(measures);
		float stdVal = stdArr(measures, meanVal);
		float lower_ci = measures[(int)round(((1 - ci_bound) / 2) * measures.size())];
		int max_pos = (int)round((ci_bound + (1 - ci_bound) / 2) * measures.size());
		if (max_pos > measures.size())
			max_pos = (int)measures.size() - 1;
		float upper_ci = measures[max_pos];

		if (it->first.size() > 4 && it->first.substr(it->first.size() - 4) == "_Obs")
			all_final_measures[it->first] = meanVal;
		else {
			all_final_measures[it->first + "_Mean"] = meanVal;
			all_final_measures[it->first + "_Std"] = stdVal;
			all_final_measures[it->first + "_CI.Lower.95"] = lower_ci;
			all_final_measures[it->first + "_CI.Upper.95"] = upper_ci;
		}
	}

	return all_final_measures;
}

map<string, map<string, float>> booststrap_analyze(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
	const map<string, vector<float>> &additional_info, const map<string, FilterCohortFunc> &filter_cohort
	, const vector<MeasurementFunctions> &meas_functions, const map<string, void *> *cohort_params,
	const vector<void *> *function_params, ProcessMeasurementParamFunc process_measurments_params,
	float sample_ratio, int sample_per_pid,
	int loopCnt, bool binary_outcome) {
	//for each pid - randomize x sample from all it's tests. do loop_times
	if (preds.size() != y.size() || preds.size() != pids.size()) {
		cerr << "bootstrap sizes aren't equal preds=" << preds.size() << " y=" << y.size() << endl;
		throw invalid_argument("bootstrap sizes aren't equal");
	}
	vector<void *> params((int)meas_functions.size());
	if (function_params == NULL)
		for (size_t i = 0; i < params.size(); ++i)
			params[i] = NULL;
	MLOG_D("Started Bootstarp Analysis on %d samples with %d cohorts\n", preds.size(), filter_cohort.size());
	time_t start = time(NULL);

	map<string, map<string, float>> all_cohorts_measurments;
	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it)
	{
		void *c_params = NULL;
		if (cohort_params != NULL && (*cohort_params).find(it->first) != (*cohort_params).end())
			c_params = (*cohort_params).at(it->first);

		vector<float> preds_c, y_c;
		vector<int> pids_c;
		vector<int> class_sz(2);
		for (size_t j = 0; j < y.size(); ++j)
			if (it->second(additional_info, (int)j, c_params)) {
				pids_c.push_back(pids[j]);
				y_c.push_back(y[j]);
				preds_c.push_back(preds[j]);
				++class_sz[y[j] > 0];
			}
		//now we have cohort: run analysis:
		string cohort_name = it->first;
		if (y_c.size() < 100) {
			MWARN("WARN: Cohort %s is too small - has %d samples. Skipping\n", cohort_name.c_str(), y_c.size());
			continue;
		}
		if (binary_outcome && (class_sz[0] < 10 || class_sz[1] < 10)) {
			MWARN("WARN: Cohort %s is too small - has %d samples with labels = [%d, %d]. Skipping\n",
				cohort_name.c_str(), y_c.size(), class_sz[0], class_sz[1]);
			continue;
		}
		map<string, float> cohort_measurments = booststrap_analyze_cohort(preds_c, y_c, pids_c,
			sample_ratio, sample_per_pid, loopCnt, meas_functions,
			function_params != NULL ? *function_params : params,
			process_measurments_params, additional_info, y, pids, it->second, c_params);

		all_cohorts_measurments[cohort_name] = cohort_measurments;
	}
	MLOG_D("Finished Bootstarp Analysis. took %2.1f seconds\n", difftime(time(NULL), start));
	return all_cohorts_measurments;
}

void write_bootstrap_results(const string &file_name, const map<string, map<string, float>> &all_cohorts_measurments) {
	string delimeter = ",";
	if (all_cohorts_measurments.empty())
		throw invalid_argument("all_cohorts_measurments can't be empty");
	unordered_set<string> all_columns_uniq;
	for (auto jt = all_cohorts_measurments.begin(); jt != all_cohorts_measurments.end(); ++jt)
		for (auto it = jt->second.begin(); it != jt->second.end(); ++it)
			all_columns_uniq.insert(it->first);
	vector<string> all_columns(all_columns_uniq.begin(), all_columns_uniq.end());
	sort(all_columns.begin(), all_columns.end());
	ofstream fw(file_name);

	fw << "Cohort_Description";
	for (size_t i = 0; i < all_columns.size(); ++i)
		fw << delimeter << all_columns[i];
	fw << endl;

	for (auto it = all_cohorts_measurments.begin(); it != all_cohorts_measurments.end(); ++it)
	{
		string cohort_name = it->first;
		map<string, float> cohort_values = it->second;
		fw << cohort_name;
		for (size_t i = 0; i < all_columns.size(); ++i)
			fw << delimeter <<
			(cohort_values.find(all_columns[i]) != cohort_values.end() ? cohort_values.at(all_columns[i]) : MED_MAT_MISSING_VALUE);
		fw << endl;
	}

	fw.close();
}
void read_bootstrap_results(const string &file_name, map<string, map<string, float>> &all_cohorts_measurments) {
	string delimeter = ",";
	ifstream of(file_name);
	string line, header;
	getline(of, header); //read header
	vector<string> column_names;
	boost::split(column_names, header, boost::is_any_of(delimeter));
	int cohort_name_ind = (int)distance(column_names.begin(), find(column_names.begin(), column_names.end(), "Cohort_Description"));
	if (cohort_name_ind > column_names.size())
		MTHROW_AND_ERR("Couldn't find \"Cohort_Description\" in bootstrap header\n");

	while (getline(of, line)) {
		vector<string> tokens;
		boost::split(tokens, line, boost::is_any_of(delimeter));
		if (tokens.size() != column_names.size())
			MTHROW_AND_ERR("Bad bootstrap format! header has %d columns. got line with %d fields. line=\"%s\"\n",
				(int)column_names.size(), (int)tokens.size(), line.c_str());
		string name = tokens[cohort_name_ind];
		map<string, float> cohort_values;
		for (size_t i = 0; i < tokens.size(); ++i)
		{
			if (i == cohort_name_ind)
				continue;
			cohort_values[column_names[i]] = stof(tokens[i]);
		}
		all_cohorts_measurments[name] = cohort_values;
	}
	of.close();
}


#pragma region Measurements Fucntions

map<string, float> calc_npos_nneg(const vector<float> &preds, const vector<float> &y, void *function_params) {
	map<string, float> res;

	map<float, int> cnts;
	for (size_t i = 0; i < y.size(); ++i)
		cnts[y[i]] += 1;

	res["NPOS"] = (float)cnts[(float)1.0];
	res["NNEG"] = (float)cnts[(float)0];

	return res;
}

map<string, float> calc_only_auc(const vector<float> &preds, const vector<float> &y, void *function_params) {
	map<string, float> res;

	vector<float> pred_threshold;
	map<float, vector<int>> pred_indexes;
	int tot_true_labels = 0;
	for (size_t i = 0; i < preds.size(); ++i)
	{
		pred_indexes[preds[i]].push_back((int)i);
		tot_true_labels += int(y[i] > 0);
	}
	int tot_false_labels = (int)y.size() - tot_true_labels;
	if (tot_true_labels == 0 || tot_false_labels == 0)
		throw invalid_argument("only falses or positives exists in cohort");
	pred_threshold = vector<float>((int)pred_indexes.size());
	map<float, vector<int>>::iterator it = pred_indexes.begin();
	for (size_t i = 0; i < pred_threshold.size(); ++i)
	{
		pred_threshold[i] = it->first;
		++it;
	}
	sort(pred_threshold.begin(), pred_threshold.end());
	//From up to down sort:
	int t_cnt = 0;
	int f_cnt = 0;
	vector<float> true_rate = vector<float>((int)pred_indexes.size());
	vector<float> false_rate = vector<float>((int)pred_indexes.size());
	int st_size = (int)pred_threshold.size() - 1;
	for (int i = st_size; i >= 0; --i)
	{
		vector<int> indexes = pred_indexes[pred_threshold[i]];
		//calc AUC status for this step:
		for (int ind : indexes)
		{
			bool true_label = y[ind] > 0;
			t_cnt += int(true_label);
			f_cnt += int(!true_label);
		}
		true_rate[st_size - i] = float(t_cnt) / tot_true_labels;
		false_rate[st_size - i] = float(f_cnt) / tot_false_labels;
	}

	float auc = false_rate[0] * true_rate[0] / 2;
	for (size_t i = 1; i < true_rate.size(); ++i)
		auc += (false_rate[i] - false_rate[i - 1]) * (true_rate[i - 1] + true_rate[i]) / 2;

	res["AUC"] = auc;

	return res;
}

map<string, float> calc_roc_measures(const vector<float> &preds, const vector<float> &y, void *function_params) {
	map<string, float> res;
	int max_qunt_vals = 10; //below it treat as "binary" bootstrap and choose those working points
	ROC_Params params;
	if (function_params != NULL)
		params = *(ROC_Params *)function_params;
	float max_diff_in_wp = params.max_diff_working_point;

	vector<float> fpr_points(params.working_point_FPR); //Working FR points:
	sort(fpr_points.begin(), fpr_points.end());
	for (size_t i = 0; i < fpr_points.size(); ++i)
		fpr_points[i] /= 100.0;
	vector<float> sens_points(params.working_point_SENS); //Working FR points:
	sort(sens_points.begin(), sens_points.end());
	for (size_t i = 0; i < sens_points.size(); ++i)
		sens_points[i] /= 100.0;

	//AUC, SPEC, SENS, Score

	vector<float> pred_threshold;
	map<float, vector<int>> pred_indexes;
	int tot_true_labels = 0;
	for (size_t i = 0; i < preds.size(); ++i)
	{
		pred_indexes[preds[i]].push_back((int)i);
		tot_true_labels += int(y[i] > 0);
	}
	int tot_false_labels = (int)y.size() - tot_true_labels;
	if (tot_true_labels == 0 || tot_false_labels == 0)
		throw invalid_argument("only falses or positives exists in cohort");
	pred_threshold = vector<float>((int)pred_indexes.size());
	map<float, vector<int>>::iterator it = pred_indexes.begin();
	for (size_t i = 0; i < pred_threshold.size(); ++i)
	{
		pred_threshold[i] = it->first;
		++it;
	}
	sort(pred_threshold.begin(), pred_threshold.end());
	bool use_wp = pred_threshold.size() > max_qunt_vals || params.use_score_working_points; //change all working points

														 //From up to down sort:
	int t_cnt = 0;
	int f_cnt = 0;
	vector<float> true_rate = vector<float>((int)pred_indexes.size());
	vector<float> false_rate = vector<float>((int)pred_indexes.size());
	int st_size = (int)pred_threshold.size() - 1;
	for (int i = st_size; i >= 0; --i)
	{
		vector<int> indexes = pred_indexes[pred_threshold[i]];
		//calc AUC status for this step:
		for (int ind : indexes)
		{
			bool true_label = y[ind] > 0;
			t_cnt += int(true_label);
			f_cnt += int(!true_label);
		}
		true_rate[st_size - i] = float(t_cnt) / tot_true_labels;
		false_rate[st_size - i] = float(f_cnt) / tot_false_labels;
	}

	float auc = false_rate[0] * true_rate[0] / 2;
	for (size_t i = 1; i < true_rate.size(); ++i)
		auc += (false_rate[i] - false_rate[i - 1]) * (true_rate[i - 1] + true_rate[i]) / 2;

	int curr_wp_fpr_ind = 0, curr_wp_sens_ind = 0;
	int i = 0;
	vector<float> wp_fpr_spec, wp_fpr_sens, wp_fpr_score;
	vector<float> wp_sens_spec, wp_sens_score;

	if (use_wp) {
		wp_fpr_spec.resize((int)fpr_points.size());
		wp_fpr_sens.resize((int)fpr_points.size());
		wp_fpr_score.resize((int)fpr_points.size());
		while (curr_wp_fpr_ind < fpr_points.size() && false_rate[i] > fpr_points[curr_wp_fpr_ind])
		{
			wp_fpr_score[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_spec[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_sens[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			++curr_wp_fpr_ind;
		}

		wp_sens_spec.resize((int)sens_points.size());
		wp_sens_score.resize((int)sens_points.size());
		while (curr_wp_sens_ind < sens_points.size() &&
			true_rate[i] > sens_points[curr_wp_sens_ind])
		{
			wp_sens_score[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_spec[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			++curr_wp_sens_ind;
		}

		//fpr points:
		while (i < true_rate.size() && curr_wp_fpr_ind < fpr_points.size())
		{
			if (curr_wp_fpr_ind < fpr_points.size() &&
				false_rate[i] >= fpr_points[curr_wp_fpr_ind]) { //passed work_point - take 2 last points for measure - by distance from wp

				float prev_diff = fpr_points[curr_wp_fpr_ind] - false_rate[i - 1];
				float curr_diff = false_rate[i] - fpr_points[curr_wp_fpr_ind];
				float tot_diff = prev_diff + curr_diff;
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					wp_fpr_score[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_sens[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_spec[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT FPR=%f, prev_FPR=%f, next_FPR=%f, prev_score=%f, next_score=%f\n",
						fpr_points[curr_wp_fpr_ind], false_rate[i - 1], false_rate[i],
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif
					++curr_wp_fpr_ind;
					continue; //skip working point - diff is too big
				}
				wp_fpr_score[curr_wp_fpr_ind] = pred_threshold[st_size - i] * (prev_diff / tot_diff) +
					pred_threshold[st_size - (i - 1)] * (curr_diff / tot_diff);
				wp_fpr_sens[curr_wp_fpr_ind] = true_rate[i] * (prev_diff / tot_diff) +
					true_rate[i - 1] * (curr_diff / tot_diff);
				wp_fpr_spec[curr_wp_fpr_ind] = (1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff);

				++curr_wp_fpr_ind;
				continue;
	}
			++i;
}

		//sens points:
		i = 0;
		while (i < true_rate.size() && curr_wp_sens_ind < sens_points.size())
		{
			if (curr_wp_sens_ind < sens_points.size() &&
				true_rate[i] >= sens_points[curr_wp_sens_ind]) { //passed work_point - take 2 last points for measure - by distance from wp

				float prev_diff = sens_points[curr_wp_sens_ind] - true_rate[i - 1];
				float curr_diff = true_rate[i] - sens_points[curr_wp_sens_ind];
				float tot_diff = prev_diff + curr_diff;
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					wp_sens_score[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_spec[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT SENS=%f, prev_SENS=%f, next_SENS=%f, prev_score=%f, next_score=%f\n",
						sens_points[curr_wp_sens_ind], true_rate[i - 1], true_rate[i],
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif
					++curr_wp_sens_ind;
					continue; //skip working point - diff is too big
				}
				wp_sens_score[curr_wp_sens_ind] = pred_threshold[st_size - i] * (prev_diff / tot_diff) +
					pred_threshold[st_size - (i - 1)] * (curr_diff / tot_diff);
				wp_sens_spec[curr_wp_sens_ind] = (1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff);

				++curr_wp_sens_ind;
				continue;
		}
			++i;
	}
}
	else {
		wp_fpr_spec.resize((int)true_rate.size());
		wp_fpr_sens.resize((int)true_rate.size());
		wp_fpr_score.resize((int)true_rate.size());
		for (i = 0; i < true_rate.size(); ++i)
		{
			wp_fpr_score[i] = pred_threshold[st_size - i];
			wp_fpr_sens[i] = true_rate[i];
			wp_fpr_spec[i] = (1 - false_rate[i]);
		}
	}


	res["AUC"] = auc;
	if (use_wp) {
		for (size_t k = 0; k < fpr_points.size(); ++k)
		{
			res["SPEC@FPR_" + print_obj(fpr_points[k] * 100, "%05.2f")] = wp_fpr_spec[k];
			res["SENS@FPR_" + print_obj(fpr_points[k] * 100, "%05.2f")] = wp_fpr_sens[k];
			res["SCORE@FPR_" + print_obj(fpr_points[k] * 100, "%05.2f")] = wp_fpr_score[k];
		}
		for (size_t k = 0; k < sens_points.size(); ++k)
		{
			res["SPEC@SENS_" + print_obj(sens_points[k] * 100, "%05.2f")] = wp_sens_spec[k];
			res["SCORE@SENS_" + print_obj(sens_points[k] * 100, "%05.2f")] = wp_sens_score[k];
		}
	}
	else
		for (size_t k = 0; k < pred_threshold.size(); ++k)
		{
			res["SPEC@SCORE_" + print_obj(wp_fpr_score[k], "%05.3f")] = wp_fpr_spec[k];
			res["SENS@SCORE_" + print_obj(wp_fpr_score[k], "%05.3f")] = wp_fpr_sens[k];
		}

	return res;
}

map<string, float> calc_roc_measures_full(const vector<float> &preds, const vector<float> &y_prob, void *function_params) {
	map<string, float> res;
	int max_qunt_vals = 10; //below it treat as "binary" bootstrap and choose those working points
	bool censor_removed = true; //wheter or not to count remove data from positive to negative

	ROC_Params params;
	if (function_params != NULL)
		params = *(ROC_Params *)function_params;
	float max_diff_in_wp = params.max_diff_working_point;

	vector<float> fpr_points = params.working_point_FPR;
	sort(fpr_points.begin(), fpr_points.end());
	for (size_t i = 0; i < fpr_points.size(); ++i)
		fpr_points[i] /= 100.0;
	vector<float> sens_points = params.working_point_SENS; //Working FR points:
	sort(sens_points.begin(), sens_points.end());
	for (size_t i = 0; i < sens_points.size(); ++i)
		sens_points[i] /= 100.0;
	vector<float> pr_points = params.working_point_PR; //Working FR points:
	sort(pr_points.begin(), pr_points.end());
	for (size_t i = 0; i < pr_points.size(); ++i)
		pr_points[i] /= 100.0;

	unordered_map<float, vector<int>> thresholds_indexes;
	vector<float> unique_scores;
	for (size_t i = 0; i < preds.size(); ++i) {
		if (thresholds_indexes.find(preds[i]) == thresholds_indexes.end())
			unique_scores.push_back(preds[i]);
		thresholds_indexes[preds[i]].push_back((int)i);
	}
	sort(unique_scores.begin(), unique_scores.end());

	//calc measures on each bucket of scores as possible threshold:
	double t_sum = 0, f_sum = 0;
	int f_cnt = 0;
	int t_cnt = 0;
	vector<float> true_rate((int)unique_scores.size());
	vector<float> false_rate((int)unique_scores.size());
	int st_size = (int)unique_scores.size() - 1;
	for (int i = st_size; i >= 0; --i)
	{
		vector<int> indexes = thresholds_indexes[unique_scores[i]];
		for (int ind : indexes)
		{
			float true_label = y_prob[ind];
			t_sum += true_label;
			if (!censor_removed)
				f_sum += (1 - true_label);
			else
				f_sum += int(true_label <= 0);
			f_cnt += int(true_label <= 0);
			t_cnt += int(true_label > 0);
		}
		true_rate[st_size - i] = float(t_sum);
		false_rate[st_size - i] = float(f_sum);
	}

	if (f_cnt == 0 || t_sum <= 0)
		throw invalid_argument("no falses or no positives exists in cohort");
	for (size_t i = 0; i < true_rate.size(); ++i) {
		true_rate[i] /= float(t_sum);
		false_rate[i] /= float(f_sum);
	}
	//calc maesures based on true_rate and false_rate
	double auc = false_rate[0] * true_rate[0] / 2; //"auc" on expectitions:
	for (size_t i = 1; i < true_rate.size(); ++i)
		auc += (false_rate[i] - false_rate[i - 1]) * (true_rate[i - 1] + true_rate[i]) / 2;

	bool use_wp = unique_scores.size() > max_qunt_vals || params.use_score_working_points; //change all working points
	int curr_wp_fpr_ind = 0, curr_wp_sens_ind = 0, curr_wp_pr_ind = 0;
	int i = 0;
	vector<float> wp_fpr_spec, wp_fpr_sens, wp_fpr_score, wp_fpr_ppv, wp_fpr_pr;
	vector<float> wp_sens_spec, wp_sens_score, wp_sens_ppv, wp_sens_pr, wp_sens_fpr;
	vector<float> wp_pr_spec, wp_pr_sens, wp_pr_score, wp_pr_ppv, wp_pr_fpr;

	if (use_wp) {
		wp_fpr_spec.resize((int)fpr_points.size());
		wp_fpr_sens.resize((int)fpr_points.size());
		wp_fpr_score.resize((int)fpr_points.size());
		wp_fpr_ppv.resize((int)fpr_points.size());
		wp_fpr_pr.resize((int)fpr_points.size());
		while (curr_wp_fpr_ind < fpr_points.size() && false_rate[i] > fpr_points[curr_wp_fpr_ind])
		{
			wp_fpr_score[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_spec[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_sens[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_ppv[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_pr[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			++curr_wp_fpr_ind;
		}

		wp_sens_score.resize((int)sens_points.size());
		wp_sens_spec.resize((int)sens_points.size());
		wp_sens_fpr.resize((int)sens_points.size());
		wp_sens_ppv.resize((int)sens_points.size());
		wp_sens_pr.resize((int)sens_points.size());
		while (curr_wp_sens_ind < sens_points.size() && true_rate[i] > sens_points[curr_wp_sens_ind])
		{
			wp_sens_score[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_spec[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_fpr[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_ppv[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_pr[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			++curr_wp_sens_ind;
		}

		wp_pr_score.resize((int)pr_points.size());
		wp_pr_spec.resize((int)pr_points.size());
		wp_pr_fpr.resize((int)pr_points.size());
		wp_pr_ppv.resize((int)pr_points.size());
		wp_pr_sens.resize((int)pr_points.size());
		while (curr_wp_pr_ind < pr_points.size() && true_rate[i] > pr_points[curr_wp_pr_ind])
		{
			wp_pr_score[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			wp_pr_spec[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			wp_pr_fpr[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			wp_pr_ppv[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			wp_pr_sens[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			++curr_wp_pr_ind;
		}

		//fpr points:
		i = 1;
		while (i < true_rate.size() && curr_wp_fpr_ind < fpr_points.size())
		{
			if (curr_wp_fpr_ind < fpr_points.size() &&
				false_rate[i] >= fpr_points[curr_wp_fpr_ind]) { //passed work_point - take 2 last points for measure - by distance from wp

				float prev_diff = fpr_points[curr_wp_fpr_ind] - false_rate[i - 1];
				float curr_diff = false_rate[i] - fpr_points[curr_wp_fpr_ind];
				float tot_diff = prev_diff + curr_diff;
				if (tot_diff <= 0) {
					curr_diff = 1;
					tot_diff = 1; //take prev - first apeareance
			}
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					wp_fpr_score[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_sens[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_spec[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_pr[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_ppv[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT FPR=%f, prev_FPR=%f, next_FPR=%f, prev_score=%f, next_score=%f\n",
						fpr_points[curr_wp_fpr_ind], false_rate[i - 1], false_rate[i],
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif
					++curr_wp_fpr_ind;
					continue; //skip working point - diff is too big
				}
				wp_fpr_score[curr_wp_fpr_ind] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				wp_fpr_sens[curr_wp_fpr_ind] = true_rate[i] * (prev_diff / tot_diff) +
					true_rate[i - 1] * (curr_diff / tot_diff);
				wp_fpr_spec[curr_wp_fpr_ind] = (1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff);
				float ppv_c = float((true_rate[i] * t_sum) /
					((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
				float ppv_prev = float((true_rate[i - 1] * t_sum) /
					((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				wp_fpr_ppv[curr_wp_fpr_ind] = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				float pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
					(t_sum + f_sum));
				float pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
					(t_sum + f_sum));
				wp_fpr_pr[curr_wp_fpr_ind] = pr_c* (prev_diff / tot_diff) + pr_prev * (curr_diff / tot_diff);

				++curr_wp_fpr_ind;
				continue;
		}
			++i;
	}

		//handle sens points:
		i = 1; //first point is always before
		while (i < true_rate.size() && curr_wp_sens_ind < sens_points.size())
		{
			if (curr_wp_sens_ind < sens_points.size() &&
				true_rate[i] >= sens_points[curr_wp_sens_ind]) { //passed work_point - take 2 last points for measure - by distance from wp

				float prev_diff = sens_points[curr_wp_sens_ind] - true_rate[i - 1];
				float curr_diff = true_rate[i] - sens_points[curr_wp_sens_ind];
				float tot_diff = prev_diff + curr_diff;
				if (tot_diff <= 0) {
					curr_diff = 1;
					tot_diff = 1; //take prev - first apeareance
			}
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					wp_sens_score[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_spec[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_fpr[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_ppv[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_pr[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT SENS=%f, prev_SENS=%f, next_SENS=%f, prev_score=%f, next_score=%f\n",
						sens_points[curr_wp_sens_ind], true_rate[i - 1], true_rate[i],
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif
					++curr_wp_sens_ind;
					continue; //skip working point - diff is too big
				}
				wp_sens_score[curr_wp_sens_ind] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				wp_sens_fpr[curr_wp_sens_ind] = false_rate[i] * (prev_diff / tot_diff) +
					false_rate[i - 1] * (curr_diff / tot_diff);
				wp_sens_spec[curr_wp_sens_ind] = (1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff);
				float ppv_c = float((true_rate[i] * t_sum) /
					((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
				float ppv_prev = float((true_rate[i - 1] * t_sum) /
					((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				wp_sens_ppv[curr_wp_sens_ind] = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				float pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
					(t_sum + f_sum));
				float pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
					(t_sum + f_sum));
				wp_sens_pr[curr_wp_sens_ind] = pr_c* (prev_diff / tot_diff) + pr_prev * (curr_diff / tot_diff);

				++curr_wp_sens_ind;
				continue;
		}
			++i;
}

		//handle pr points:
		i = 1; //first point is always before
		while (i < true_rate.size() && curr_wp_pr_ind < pr_points.size())
		{
			float pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
				(t_sum + f_sum));
			if (curr_wp_pr_ind < pr_points.size() && pr_c >= pr_points[curr_wp_pr_ind]) { //passed work_point - take 2 last points for measure - by distance from wp
				float pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
					(t_sum + f_sum));

				float prev_diff = pr_points[curr_wp_pr_ind] - pr_prev;
				float curr_diff = pr_c - pr_points[curr_wp_pr_ind];
				float tot_diff = prev_diff + curr_diff;
				if (tot_diff <= 0) {
					curr_diff = 1;
					tot_diff = 1; //take prev - first apeareance
			}
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					wp_pr_score[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
					wp_pr_fpr[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
					wp_pr_spec[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
					wp_pr_ppv[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
					wp_pr_sens[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT PR=%f, prev_PR=%f, next_PR=%f, prev_score=%f, next_score=%f\n",
						pr_points[curr_wp_pr_ind], pr_prev, pr_c,
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif //  WARN_SKIP_WP
					++curr_wp_pr_ind;
					continue; //skip working point - diff is too big
				}
				wp_pr_score[curr_wp_pr_ind] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				wp_pr_fpr[curr_wp_pr_ind] = false_rate[i] * (prev_diff / tot_diff) +
					false_rate[i - 1] * (curr_diff / tot_diff);
				wp_pr_spec[curr_wp_pr_ind] = (1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff);
				float ppv_c = float((true_rate[i] * t_sum) /
					((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
				float ppv_prev = float((true_rate[i - 1] * t_sum) /
					((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				wp_pr_ppv[curr_wp_pr_ind] = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				wp_pr_sens[curr_wp_pr_ind] = true_rate[i] * (prev_diff / tot_diff) + true_rate[i - 1] * (curr_diff / tot_diff);

				++curr_wp_pr_ind;
				continue;
		}
			++i;
		}

	}
	else {
		wp_fpr_spec.resize((int)true_rate.size());
		wp_fpr_sens.resize((int)true_rate.size());
		wp_fpr_score.resize((int)true_rate.size());
		wp_fpr_ppv.resize((int)true_rate.size());
		wp_fpr_pr.resize((int)true_rate.size());
		for (i = 0; i < true_rate.size(); ++i)
		{
			wp_fpr_score[i] = unique_scores[st_size - i];
			wp_fpr_sens[i] = true_rate[i];
			wp_fpr_spec[i] = (1 - false_rate[i]);
			wp_fpr_ppv[i] = float((true_rate[i] * t_sum) /
				((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
			wp_fpr_pr[i] = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
				(t_sum + f_sum));
		}
	}


	res["AUC"] = float(auc);
	if (use_wp) {
		for (size_t k = 0; k < fpr_points.size(); ++k)
		{
			res["SPEC@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_spec[k];
			res["SENS@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_sens[k];
			res["SCORE@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_score[k];
			res["PPV@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_ppv[k];
			res["PR@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_pr[k];
		}
		for (size_t k = 0; k < sens_points.size(); ++k)
		{
			res["SPEC@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_spec[k];
			res["FPR@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_fpr[k];
			res["SCORE@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_score[k];
			res["PPV@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_ppv[k];
			res["PR@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_pr[k];
		}
		for (size_t k = 0; k < pr_points.size(); ++k)
		{
			res["SPEC@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_spec[k];
			res["SENS@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_sens[k];
			res["SCORE@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_score[k];
			res["PPV@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_ppv[k];
			res["FPR@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_fpr[k];
		}
	}
	else
		for (size_t k = 0; k < unique_scores.size(); ++k)
		{
			res["SPEC@SCORE_" + print_obj(wp_fpr_score[k], "%06.3f")] = wp_fpr_spec[k];
			res["SENS@SCORE_" + print_obj(wp_fpr_score[k], "%06.3f")] = wp_fpr_sens[k];
			res["PPV@SCORE_" + print_obj(wp_fpr_score[k], "%05.3f")] = wp_fpr_ppv[k];
			res["PR@SCORE_" + print_obj(wp_fpr_score[k], "%05.3f")] = wp_fpr_pr[k];
		}

	res["NEG_SUM"] = float(f_sum);
	res["POS_SUM"] = float(t_sum);
	res["POS_CNT"] = float(t_cnt);
	res["NEG_CNT"] = float(f_cnt);

	return res;
}

class compare_pair
{
	bool reverse;
public:
	compare_pair(const bool& revparam = true)
	{
		reverse = revparam;
	}
	bool operator() (const pair<int, int>& lhs, const pair<int, int>&rhs) const
	{
		if (reverse) return (lhs.first > rhs.first);
		else return (lhs.first < rhs.first);
	}
};

void merge_down(vector<int> &ind_to_size, vector<vector<pair<int, int>>> &size_to_ind, set<int> &sizes,
	const pair<int, int> *index_to_merge, int min_size) {
	pair<int, int> *merge_into = NULL;
	int to_merge_size = ind_to_size[index_to_merge->first - 1];
	//remove index_to_merge.first - 1:
	for (size_t j = 0; j < size_to_ind[to_merge_size].size(); ++j)
		if (size_to_ind[to_merge_size][j].second <= index_to_merge->first - 1 &&
			size_to_ind[to_merge_size][j].first >= index_to_merge->first - 1) {
			merge_into = &size_to_ind[to_merge_size][j];
			size_to_ind[to_merge_size].erase(size_to_ind[to_merge_size].begin() + j);
			break;
		}
	if (merge_into == NULL)
		MTHROW_AND_ERR("Bug couldn't found merge_into\n");

	int new_size = min_size + to_merge_size;
	sizes.insert(new_size);
	//insert new union
	size_to_ind[new_size].push_back(pair<int, int>(merge_into->first, index_to_merge->second));
	//update in min,max:
	ind_to_size[merge_into->first] = new_size;
	ind_to_size[index_to_merge->second] = new_size;
}
void merge_up(vector<int> &ind_to_size, vector<vector<pair<int, int>>> &size_to_ind, set<int> &sizes,
	const pair<int, int> *index_to_merge, int min_size) {
	//merge with +1
	pair<int, int> *merge_into = NULL;
	int to_merge_size = ind_to_size[index_to_merge->second + 1];
	//remove index_to_merge.second + 1:
	for (size_t j = 0; j < size_to_ind[to_merge_size].size(); ++j)
		if (size_to_ind[to_merge_size][j].second <= index_to_merge->second + 1 &&
			size_to_ind[to_merge_size][j].first >= index_to_merge->second + 1) {
			merge_into = &size_to_ind[to_merge_size][j];
			size_to_ind[to_merge_size].erase(size_to_ind[to_merge_size].begin() + j);
			break;
		}
	if (merge_into == NULL)
		MTHROW_AND_ERR("Bug couldn't found merge_into\n");

	int new_size = min_size + to_merge_size;
	sizes.insert(new_size);
	//insert new union set:
	size_to_ind[new_size].push_back(pair<int, int>(index_to_merge->first, merge_into->second));
	//update in min,max:
	ind_to_size[index_to_merge->first] = new_size;
	ind_to_size[merge_into->second] = new_size;
}

map<string, float> calc_roc_measures_with_inc(const vector<float> &preds, const vector<float> &y, void *function_params) {
	map<string, float> res;
	int max_qunt_vals = 10;
	bool censor_removed = true;

	ROC_Params params;
	if (function_params != NULL)
		params = *(ROC_Params *)function_params;
	float max_diff_in_wp = params.max_diff_working_point;
	int scores_bin = params.score_bins;

	vector<float> fpr_points = params.working_point_FPR;
	sort(fpr_points.begin(), fpr_points.end());
	for (size_t i = 0; i < fpr_points.size(); ++i)
		fpr_points[i] /= 100.0;
	vector<float> sens_points = params.working_point_SENS; //Working FR points:
	sort(sens_points.begin(), sens_points.end());
	for (size_t i = 0; i < sens_points.size(); ++i)
		sens_points[i] /= 100.0;
	vector<float> pr_points = params.working_point_PR; //Working FR points:
	sort(pr_points.begin(), pr_points.end());
	for (size_t i = 0; i < pr_points.size(); ++i)
		pr_points[i] /= 100.0;

	unordered_map<float, vector<int>> thresholds_indexes;
	vector<float> unique_scores;
	for (size_t i = 0; i < preds.size(); ++i)
		thresholds_indexes[preds[i]].push_back((int)i);
	unique_scores.resize((int)thresholds_indexes.size());
	int ind_p = 0;
	for (auto it = thresholds_indexes.begin(); it != thresholds_indexes.end(); ++it)
	{
		unique_scores[ind_p] = it->first;
		++ind_p;
	}
	sort(unique_scores.begin(), unique_scores.end());

	bool clear_mem = false;
	if (scores_bin > 0 && unique_scores.size() > scores_bin) {
		clear_mem = true;
		int bin_size_last = (int)unique_scores.size();
		int c = 0;
		vector<vector<pair<int, int>>> size_to_ind(preds.size()); //size, group, index_min_max
		vector<int> ind_to_size(bin_size_last);
		set<int> sizes;
		for (auto it = thresholds_indexes.begin(); it != thresholds_indexes.end(); ++it)
		{
			size_to_ind[(int)it->second.size()].push_back(pair<int, int>(c, c));
			ind_to_size[c] = (int)it->second.size();
			++c;
			sizes.insert((int)it->second.size());
		}

		while (bin_size_last > scores_bin) {
			int min_size = *sizes.begin();
			pair<int, int> *index_to_merge = &size_to_ind[min_size].back();
			size_to_ind[min_size].pop_back();
			if (size_to_ind[min_size].empty()) //erase if left empty after merge
				sizes.erase(sizes.begin());
			pair<int, int> *merge_into = NULL;
			//merge index_to_merge with index_to_merge+-1. and update size_to_ind, ind_to_size, sizes
			if (index_to_merge->second == unique_scores.size() - 1)
				merge_down(ind_to_size, size_to_ind, sizes, index_to_merge, min_size);
			else if (index_to_merge->first == 0)
				merge_up(ind_to_size, size_to_ind, sizes, index_to_merge, min_size);
			else {
				if (ind_to_size[index_to_merge->second + 1] < ind_to_size[index_to_merge->first - 1])
					merge_up(ind_to_size, size_to_ind, sizes, index_to_merge, min_size);
				else
					merge_down(ind_to_size, size_to_ind, sizes, index_to_merge, min_size);
			}

			--bin_size_last;
		}

		//update thresholds_indexes based on: size_to_ind groups -
		//merge all indexes in each group to first index in thresholds_indexes. "mean" other scores to unique_scores

		unordered_map<float, vector<int>> *new_thresholds = new unordered_map<float, vector<int>>;
		for (auto it = sizes.begin(); it != sizes.end(); ++it)
		{
			for (size_t k = 0; k < size_to_ind[*it].size(); ++k)
			{ //merge from first => second
				pair<int, int> *merge = &size_to_ind[*it][k];
				double mean_score = 0, tot_cnt = 0;
				vector<int> *merged_inds = new vector<int>;
				for (int ii = merge->first; ii <= merge->second; ++ii) {
					mean_score += unique_scores[ii] * thresholds_indexes[unique_scores[ii]].size();
					tot_cnt += thresholds_indexes[unique_scores[ii]].size();
					merged_inds->insert(merged_inds->end(),
						thresholds_indexes[unique_scores[ii]].begin(), thresholds_indexes[unique_scores[ii]].end());
				}
				mean_score /= tot_cnt;
				(*new_thresholds)[(float)mean_score] = move(*merged_inds);
			}
		}

		thresholds_indexes = move(*new_thresholds);
		unique_scores.resize((int)thresholds_indexes.size());
		int ind_pp = 0;
		for (auto it = thresholds_indexes.begin(); it != thresholds_indexes.end(); ++it) {
			unique_scores[ind_pp] = it->first;
			++ind_pp;
		}
		sort(unique_scores.begin(), unique_scores.end());
	}

	//calc measures on each bucket of scores as possible threshold:
	double t_sum = 0, f_sum = 0;
	int f_cnt = 0;
	int t_cnt = 0;
	vector<float> true_rate((int)unique_scores.size());
	vector<float> false_rate((int)unique_scores.size());
	int st_size = (int)unique_scores.size() - 1;
	for (int i = st_size; i >= 0; --i)
	{
		vector<int> indexes = thresholds_indexes[unique_scores[i]];
		for (int ind : indexes)
		{
			float true_label = y[ind];
			t_sum += true_label;
			if (!censor_removed)
				f_sum += (1 - true_label);
			else
				f_sum += int(true_label <= 0);
			f_cnt += int(true_label <= 0);
			t_cnt += int(true_label > 0);
		}
		true_rate[st_size - i] = float(t_sum);
		false_rate[st_size - i] = float(f_sum);
	}

	if (f_cnt == 0 || t_sum <= 0)
		throw invalid_argument("no falses or no positives exists in cohort");
	for (size_t i = 0; i < true_rate.size(); ++i) {
		true_rate[i] /= float(t_sum);
		false_rate[i] /= float(f_sum);
	}
	//calc maesures based on true_rate and false_rate
	double auc = false_rate[0] * true_rate[0] / 2; //"auc" on expectitions:
	for (size_t i = 1; i < true_rate.size(); ++i)
		auc += (false_rate[i] - false_rate[i - 1]) * (true_rate[i - 1] + true_rate[i]) / 2;

	bool use_wp = unique_scores.size() > max_qunt_vals || params.use_score_working_points; //change all working points
	int curr_wp_fpr_ind = 0, curr_wp_sens_ind = 0, curr_wp_pr_ind = 0;
	int i = 0;
	vector<float> wp_fpr_spec, wp_fpr_sens, wp_fpr_score, wp_fpr_ppv, wp_fpr_pr;
	vector<float> wp_sens_spec, wp_sens_score, wp_sens_ppv, wp_sens_pr, wp_sens_fpr;
	vector<float> wp_pr_spec, wp_pr_sens, wp_pr_score, wp_pr_ppv, wp_pr_fpr;

	if (use_wp) {
		wp_fpr_spec.resize((int)fpr_points.size());
		wp_fpr_sens.resize((int)fpr_points.size());
		wp_fpr_score.resize((int)fpr_points.size());
		wp_fpr_ppv.resize((int)fpr_points.size());
		wp_fpr_pr.resize((int)fpr_points.size());
		while (curr_wp_fpr_ind < fpr_points.size() && false_rate[i] > fpr_points[curr_wp_fpr_ind])
		{
			wp_fpr_score[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_spec[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_sens[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_ppv[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			wp_fpr_pr[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
			++curr_wp_fpr_ind;
		}

		wp_sens_score.resize((int)sens_points.size());
		wp_sens_spec.resize((int)sens_points.size());
		wp_sens_fpr.resize((int)sens_points.size());
		wp_sens_ppv.resize((int)sens_points.size());
		wp_sens_pr.resize((int)sens_points.size());
		while (curr_wp_sens_ind < sens_points.size() && true_rate[i] > sens_points[curr_wp_sens_ind])
		{
			wp_sens_score[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_spec[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_fpr[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_ppv[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			wp_sens_pr[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
			++curr_wp_sens_ind;
		}

		wp_pr_score.resize((int)pr_points.size());
		wp_pr_spec.resize((int)pr_points.size());
		wp_pr_fpr.resize((int)pr_points.size());
		wp_pr_ppv.resize((int)pr_points.size());
		wp_pr_sens.resize((int)pr_points.size());
		while (curr_wp_pr_ind < pr_points.size() && true_rate[i] > pr_points[curr_wp_pr_ind])
		{
			wp_pr_score[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			wp_pr_spec[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			wp_pr_fpr[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			wp_pr_ppv[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			wp_pr_sens[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
			++curr_wp_pr_ind;
		}

		//fpr points:
		i = 1;
		while (i < true_rate.size() && curr_wp_fpr_ind < fpr_points.size())
		{
			if (curr_wp_fpr_ind < fpr_points.size() &&
				false_rate[i] >= fpr_points[curr_wp_fpr_ind]) { //passed work_point - take 2 last points for measure - by distance from wp

				float prev_diff = fpr_points[curr_wp_fpr_ind] - false_rate[i - 1];
				float curr_diff = false_rate[i] - fpr_points[curr_wp_fpr_ind];
				float tot_diff = prev_diff + curr_diff;
				if (tot_diff <= 0) {
					curr_diff = 1;
					tot_diff = 1; //take prev - first apeareance
			}
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					wp_fpr_score[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_sens[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_spec[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_pr[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
					wp_fpr_ppv[curr_wp_fpr_ind] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT FPR=%f, prev_FPR=%f, next_FPR=%f, prev_score=%f, next_score=%f\n",
						fpr_points[curr_wp_fpr_ind], false_rate[i - 1], false_rate[i],
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif
					++curr_wp_fpr_ind;
					continue; //skip working point - diff is too big
				}
				wp_fpr_score[curr_wp_fpr_ind] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				wp_fpr_sens[curr_wp_fpr_ind] = true_rate[i] * (prev_diff / tot_diff) +
					true_rate[i - 1] * (curr_diff / tot_diff);
				wp_fpr_spec[curr_wp_fpr_ind] = (1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff);
				float ppv_c = float((true_rate[i] * t_sum) /
					((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
				float ppv_prev = float((true_rate[i - 1] * t_sum) /
					((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				if (params.incidence_fix > 0) {
					ppv_c = float(params.incidence_fix*(true_rate[i] * t_sum) /
						(params.incidence_fix*(true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
					ppv_prev = float(params.incidence_fix*(true_rate[i - 1] * t_sum) /
						(params.incidence_fix*(true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				}
				wp_fpr_ppv[curr_wp_fpr_ind] = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				float pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
					(t_sum + f_sum));
				float pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
					(t_sum + f_sum));
				if (params.incidence_fix > 0) {
					pr_c = float((params.incidence_fix*(true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
						(params.incidence_fix*t_sum + f_sum));
					pr_prev = float((params.incidence_fix*(true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
						(params.incidence_fix*t_sum + f_sum));
				}
				wp_fpr_pr[curr_wp_fpr_ind] = pr_c* (prev_diff / tot_diff) + pr_prev * (curr_diff / tot_diff);

				++curr_wp_fpr_ind;
				continue;
		}
			++i;
	}

		//handle sens points:
		i = 1; //first point is always before
		while (i < true_rate.size() && curr_wp_sens_ind < sens_points.size())
		{
			if (curr_wp_sens_ind < sens_points.size() &&
				true_rate[i] >= sens_points[curr_wp_sens_ind]) { //passed work_point - take 2 last points for measure - by distance from wp

				float prev_diff = sens_points[curr_wp_sens_ind] - true_rate[i - 1];
				float curr_diff = true_rate[i] - sens_points[curr_wp_sens_ind];
				float tot_diff = prev_diff + curr_diff;
				if (tot_diff <= 0) {
					curr_diff = 1;
					tot_diff = 1; //take prev - first apeareance
			}
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					wp_sens_score[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_spec[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_fpr[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_ppv[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
					wp_sens_pr[curr_wp_sens_ind] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT SENS=%f, prev_SENS=%f, next_SENS=%f, prev_score=%f, next_score=%f\n",
						sens_points[curr_wp_sens_ind], true_rate[i - 1], true_rate[i],
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif
					++curr_wp_sens_ind;
					continue; //skip working point - diff is too big
				}
				wp_sens_score[curr_wp_sens_ind] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				wp_sens_fpr[curr_wp_sens_ind] = false_rate[i] * (prev_diff / tot_diff) +
					false_rate[i - 1] * (curr_diff / tot_diff);
				wp_sens_spec[curr_wp_sens_ind] = (1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff);
				float ppv_c = float((true_rate[i] * t_sum) /
					((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
				float ppv_prev = float((true_rate[i - 1] * t_sum) /
					((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				if (params.incidence_fix > 0) {
					ppv_c = float(params.incidence_fix*(true_rate[i] * t_sum) /
						(params.incidence_fix*(true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
					ppv_prev = float(params.incidence_fix*(true_rate[i - 1] * t_sum) /
						(params.incidence_fix*(true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				}
				wp_sens_ppv[curr_wp_sens_ind] = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				float pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
					(t_sum + f_sum));
				float pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
					(t_sum + f_sum));
				if (params.incidence_fix > 0) {
					pr_c = float((params.incidence_fix*(true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
						(params.incidence_fix*t_sum + f_sum));
					pr_prev = float((params.incidence_fix*(true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
						(params.incidence_fix*t_sum + f_sum));
				}
				wp_sens_pr[curr_wp_sens_ind] = pr_c* (prev_diff / tot_diff) + pr_prev * (curr_diff / tot_diff);

				++curr_wp_sens_ind;
				continue;
		}
			++i;
}

		//handle pr points:
		i = 1; //first point is always before
		while (i < true_rate.size() && curr_wp_pr_ind < pr_points.size())
		{
			float pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
				(t_sum + f_sum));
			if (curr_wp_pr_ind < pr_points.size() && pr_c >= pr_points[curr_wp_pr_ind]) { //passed work_point - take 2 last points for measure - by distance from wp
				float pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
					(t_sum + f_sum));

				float prev_diff = pr_points[curr_wp_pr_ind] - pr_prev;
				float curr_diff = pr_c - pr_points[curr_wp_pr_ind];
				float tot_diff = prev_diff + curr_diff;
				if (tot_diff <= 0) {
					curr_diff = 1;
					tot_diff = 1; //take prev - first apeareance
			}
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					wp_pr_score[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
					wp_pr_fpr[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
					wp_pr_spec[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
					wp_pr_ppv[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
					wp_pr_sens[curr_wp_pr_ind] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT PR=%f, prev_PR=%f, next_PR=%f, prev_score=%f, next_score=%f\n",
						pr_points[curr_wp_pr_ind], pr_prev, pr_c,
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif //  WARN_SKIP_WP
					++curr_wp_pr_ind;
					continue; //skip working point - diff is too big
				}
				wp_pr_score[curr_wp_pr_ind] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				wp_pr_fpr[curr_wp_pr_ind] = false_rate[i] * (prev_diff / tot_diff) +
					false_rate[i - 1] * (curr_diff / tot_diff);
				wp_pr_spec[curr_wp_pr_ind] = (1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff);
				float ppv_c = float((true_rate[i] * t_sum) /
					((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
				float ppv_prev = float((true_rate[i - 1] * t_sum) /
					((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				if (params.incidence_fix > 0) {
					ppv_c = float(params.incidence_fix*(true_rate[i] * t_sum) /
						(params.incidence_fix*(true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
					ppv_prev = float(params.incidence_fix*(true_rate[i - 1] * t_sum) /
						(params.incidence_fix*(true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
				}
				wp_pr_ppv[curr_wp_pr_ind] = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				wp_pr_sens[curr_wp_pr_ind] = true_rate[i] * (prev_diff / tot_diff) + true_rate[i - 1] * (curr_diff / tot_diff);

				++curr_wp_pr_ind;
				continue;
		}
			++i;
		}

	}
	else {
		wp_fpr_spec.resize((int)true_rate.size());
		wp_fpr_sens.resize((int)true_rate.size());
		wp_fpr_score.resize((int)true_rate.size());
		wp_fpr_ppv.resize((int)true_rate.size());
		wp_fpr_pr.resize((int)true_rate.size());
		for (i = 0; i < true_rate.size(); ++i)
		{
			wp_fpr_score[i] = unique_scores[st_size - i];
			wp_fpr_sens[i] = true_rate[i];
			wp_fpr_spec[i] = (1 - false_rate[i]);
			wp_fpr_ppv[i] = float((true_rate[i] * t_sum) /
				((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
			if (params.incidence_fix > 0) {
				wp_fpr_ppv[i] = float(params.incidence_fix*(true_rate[i] * t_sum) /
					(params.incidence_fix*(true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
			}
			wp_fpr_pr[i] = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
				(t_sum + f_sum));
			if (params.incidence_fix > 0) {
				wp_fpr_pr[i] = float((params.incidence_fix*(true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
					(params.incidence_fix*t_sum + f_sum));
			}
		}
	}


	res["AUC"] = float(auc);
	if (use_wp) {
		for (size_t k = 0; k < fpr_points.size(); ++k)
		{
			res["SPEC@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_spec[k];
			res["SENS@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_sens[k];
			res["SCORE@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_score[k];
			res["PPV@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_ppv[k];
			res["PR@FPR_" + print_obj(fpr_points[k] * 100, "%06.3f")] = wp_fpr_pr[k];
		}
		for (size_t k = 0; k < sens_points.size(); ++k)
		{
			res["SPEC@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_spec[k];
			res["FPR@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_fpr[k];
			res["SCORE@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_score[k];
			res["PPV@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_ppv[k];
			res["PR@SENS_" + print_obj(sens_points[k] * 100, "%06.3f")] = wp_sens_pr[k];
		}
		for (size_t k = 0; k < pr_points.size(); ++k)
		{
			res["SPEC@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_spec[k];
			res["SENS@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_sens[k];
			res["SCORE@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_score[k];
			res["PPV@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_ppv[k];
			res["FPR@PR_" + print_obj(pr_points[k] * 100, "%06.3f")] = wp_pr_fpr[k];
		}
	}
	else
		for (size_t k = 0; k < unique_scores.size(); ++k)
		{
			res["SPEC@SCORE_" + print_obj(wp_fpr_score[k], "%06.3f")] = wp_fpr_spec[k];
			res["SENS@SCORE_" + print_obj(wp_fpr_score[k], "%06.3f")] = wp_fpr_sens[k];
			res["PPV@SCORE_" + print_obj(wp_fpr_score[k], "%06.3f")] = wp_fpr_ppv[k];
			res["PR@SCORE_" + print_obj(wp_fpr_score[k], "%06.3f")] = wp_fpr_pr[k];
		}

	res["NEG_SUM"] = float(f_sum);
	res["POS_SUM"] = float(t_sum);
	res["POS_CNT"] = float(t_cnt);
	res["NEG_CNT"] = float(f_cnt);

	if (clear_mem) {
		for (auto it = thresholds_indexes.begin(); it != thresholds_indexes.end(); ++it)
			delete &it->second;
		delete &thresholds_indexes;
	}

	return res;
}

#pragma endregion

#pragma region Cohort Fucntions
bool time_range_filter(float outcome, int min_time, int max_time, int time, int outcome_time) {
	if (med_time.YearsMonths2Days.empty())
		med_time.init_time_tables();
	int diff_days = (med_time.convert_date(MedTime::Days, outcome_time) -
		med_time.convert_date(MedTime::Days, time));
	return ((outcome > 0 && diff_days >= min_time && diff_days <= max_time) ||
		(outcome <= 0 && diff_days > max_time));
}
bool time_range_filter(float outcome, float min_time, float max_time, float diff_days) {
	return ((outcome > 0 && diff_days >= min_time && diff_days <= max_time) ||
		(outcome <= 0 && diff_days > max_time));
}

bool filter_range_param(const map<string, vector<float>> &record_info, int index, void *cohort_params) {
	Filter_Param *param = (Filter_Param *)cohort_params; //can't be null
	if (param->param_name != "Time-Window")
		return record_info.at(param->param_name)[index] >= param->min_range &&
		record_info.at(param->param_name)[index] <= param->max_range;
	else
		return time_range_filter(record_info.at(param->param_name)[index] > 0, param->min_range,
			param->max_range, abs(record_info.at(param->param_name)[index]));
}

bool filter_range_params(const map<string, vector<float>> &record_info, int index, void *cohort_params) {
	vector<Filter_Param> *param = (vector<Filter_Param> *)cohort_params; //can't be null
	bool res = true;
	int i = 0;
	while (res && i < (*param).size()) {
		if ((*param)[i].param_name != "Time-Window")
			res = record_info.at((*param)[i].param_name)[index] >= (*param)[i].min_range &&
			record_info.at((*param)[i].param_name)[index] <= (*param)[i].max_range;
		else
			res = time_range_filter(record_info.at((*param)[i].param_name)[index] > 0, (*param)[i].min_range,
				(*param)[i].max_range, abs(record_info.at((*param)[i].param_name)[index]));
		++i;
	}
	return res;
}
#pragma endregion

#pragma region Process Measurement Param Functions
void fix_cohort_sample_incidence(const map<string, vector<float>> &additional_info,
	const vector<float> &y, const vector<int> &pids, FilterCohortFunc cohort_def,
	void *cohort_params, void *function_params) {
	ROC_Params *params = (ROC_Params *)function_params;
	//calculating the "fixed" incidence in the cohort giving the true inc. in the general population
	// select cohort - and multiply in the given original incidence
	if (params->inc_stats.sorted_outcome_labels.size() != 2)
		MTHROW_AND_ERR("Category outcome aren't supported for now\n");
	if (additional_info.find("Age") == additional_info.end() || additional_info.find("Gender") == additional_info.end())
		MTHROW_AND_ERR("Age or Gender Signals are missings\n");

	int bin_counts = (int)floor((params->inc_stats.max_age - params->inc_stats.min_age) / params->inc_stats.age_bin_years);
	if (bin_counts * params->inc_stats.age_bin_years >
		(params->inc_stats.max_age - params->inc_stats.min_age) + 0.5)
		++bin_counts; //has at least 0.5 years for last bin to create it
	if (params->inc_stats.male_labels_count_per_age.size() != bin_counts)
		MTHROW_AND_ERR("Male vector has %d members. and need to have %d members\n",
			(int)params->inc_stats.male_labels_count_per_age.size(), bin_counts);

	MLOG_D("Running fix_cohort_sample_incidence...\n");
	vector<vector<double>> male_counts(2), female_counts(2);
	vector<vector<double>> filtered_male_counts(2), filtered_female_counts(2);
	for (size_t i = 0; i < male_counts.size(); ++i)
		male_counts[i].resize(bin_counts);
	for (size_t i = 0; i < female_counts.size(); ++i)
		female_counts[i].resize(bin_counts);
	for (size_t i = 0; i < filtered_male_counts.size(); ++i)
		filtered_male_counts[i].resize(bin_counts);
	for (size_t i = 0; i < filtered_female_counts.size(); ++i)
		filtered_female_counts[i].resize(bin_counts);

	for (size_t i = 0; i < y.size(); ++i)
	{
		if (additional_info.at("Age")[i] < params->inc_stats.min_age ||
			additional_info.at("Age")[i] > params->inc_stats.max_age)
			continue; //skip out of range
		int age_index = (int)floor((additional_info.at("Age")[i] - params->inc_stats.min_age) /
			params->inc_stats.age_bin_years);
		if (age_index >= bin_counts)
			age_index = bin_counts - 1;

		if (additional_info.at("Gender")[i] == 1) { //Male
			++male_counts[y[i] > 0][age_index];
			if (cohort_def(additional_info, (int)i, cohort_params))
				++filtered_male_counts[y[i] > 0][age_index];
		}
		else {//Female
			++female_counts[y[i] > 0][age_index];
			if (cohort_def(additional_info, (int)i, cohort_params))
				++filtered_female_counts[y[i] > 0][age_index];
		}
	}

	//params->cohort_inc_stats.age_bin_years = params->inc_stats.age_bin_years;
	//params->cohort_inc_stats.max_age = params->inc_stats.max_age;
	//params->cohort_inc_stats.min_age = params->inc_stats.min_age;
	//params->cohort_inc_stats.sorted_outcome_labels = params->inc_stats.sorted_outcome_labels;
	Incident_Stats cohort_inc_stats;
	cohort_inc_stats.male_labels_count_per_age.resize(bin_counts);
	cohort_inc_stats.female_labels_count_per_age.resize(bin_counts);
	//recalc new ratio of 1/0 and fix stats
	for (size_t i = 0; i < bin_counts; ++i)
	{
		if (i == 0) {
			cohort_inc_stats.male_labels_count_per_age[i].resize(2);
			cohort_inc_stats.female_labels_count_per_age[i].resize(2);
		}
		if (male_counts[i][0] > 0)
			cohort_inc_stats.male_labels_count_per_age[i][1] =
			(male_counts[i][1] / male_counts[i][0]) * params->inc_stats.male_labels_count_per_age[i][1];
		else {
			MWARN("Warning:: has too little controls\n");
			cohort_inc_stats.male_labels_count_per_age[i][1] = params->inc_stats.male_labels_count_per_age[i][1];
		}

		if (female_counts[i][0] > 0)
			cohort_inc_stats.female_labels_count_per_age[i][1] =
			(female_counts[i][1] / female_counts[i][0]) * params->inc_stats.female_labels_count_per_age[i][1];
		else {
			MWARN("Warning:: has too litle controls\n");
			cohort_inc_stats.female_labels_count_per_age[i][1] = params->inc_stats.female_labels_count_per_age[i][1];
		}
	}

	//calc mean fix:
	params->incidence_fix = 0;
	double controls_count = 0;
	for (size_t i = 0; i < bin_counts; ++i)
	{
		params->incidence_fix += cohort_inc_stats.male_labels_count_per_age[i][1];
		params->incidence_fix += cohort_inc_stats.female_labels_count_per_age[i][1];

		controls_count += cohort_inc_stats.male_labels_count_per_age[i][0];
		controls_count += cohort_inc_stats.female_labels_count_per_age[i][0];
	}

	params->incidence_fix /= controls_count;
}
#pragma endregion

#pragma region Parameter Functions
Filter_Param::Filter_Param(const string &init_string) {
	if (init_string.find(':') == string::npos)
		MTHROW_AND_ERR("Wrong format given \"%s\". expected format is \"PARAM_NAME:min_range,max_range\"\n",
			init_string.c_str());
	param_name = init_string.substr(0, init_string.find(':'));
	string rest = init_string.substr(init_string.find(':') + 1);
	if (rest.find(',') == string::npos)
		MTHROW_AND_ERR("Wrong format given \"%s\". expected format is \"PARAM_NAME:min_range,max_range\"\n",
			init_string.c_str());
	min_range = stof(rest.substr(0, rest.find(',')));
	max_range = stof(rest.substr(rest.find(',') + 1));
}
void Incident_Stats::write_to_text_file(const string &text_file) {
	ofstream fw(text_file);
	string delim = "\t";
	fw << "AGE_BIN" << delim << age_bin_years << endl;
	fw << "AGE_MIN" << delim << min_age << endl;
	fw << "AGE_MAX" << delim << max_age << endl;
	for (size_t i = 0; i < sorted_outcome_labels.size(); ++i)
		fw << "OUTCOME_VALUE" << delim << sorted_outcome_labels[i] << "\n";
	fw.flush();
	for (size_t i = 0; i < male_labels_count_per_age.size(); ++i)
		for (size_t j = 0; j < male_labels_count_per_age[i].size(); ++j)
			fw << "STATS_ROW" << delim << "MALE" << delim << min_age + i*age_bin_years
			<< delim << sorted_outcome_labels[j] << delim << male_labels_count_per_age[i][j] << "\n";
	for (size_t i = 0; i < female_labels_count_per_age.size(); ++i)
		for (size_t j = 0; j < female_labels_count_per_age[i].size(); ++j)
			fw << "STATS_ROW" << delim << "FEMALE" << delim << min_age + i*age_bin_years
			<< delim << sorted_outcome_labels[j] << delim << female_labels_count_per_age[i][j] << "\n";
	fw.flush();
	fw.close();
}
void Incident_Stats::read_from_text_file(const string &text_file) {
	ifstream of(text_file);
	string line;
	while (getline(of, line)) {
		if (line.empty() || boost::starts_with(line, "#"))
			continue;
		vector<string> tokens;
		boost::split(tokens, line, boost::is_any_of("\t"));
		if (tokens.size() < 2)
			MTHROW_AND_ERR("Format Error: got line: \"%s\"\n", line.c_str());
		string command = tokens[0];
		if (command == "AGE_BIN")
			age_bin_years = stoi(tokens[1]);
		else if (command == "AGE_MIN")
			min_age = stof(tokens[1]);
		else if (command == "AGE_MAX")
			max_age = stof(tokens[1]);
		else if (command == "OUTCOME_VALUE")
			sorted_outcome_labels.push_back(stof(tokens[1]));
		else if (command == "STATS_ROW") {
			if (tokens.size() != 5)
				MTHROW_AND_ERR("Unknown lines format \"%s\"\n", line.c_str());
			float age = stof(tokens[2]);
			if (age < min_age || age> max_age) {
				MWARN("Warning:: skip age because out of range in line \"%s\"", line.c_str());
				continue;
			}
			int age_bin = (int)floor((age - min_age) / age_bin_years);
			int max_bins = (int)floor((max_age - min_age) / age_bin_years);
			if (max_bins * age_bin_years > (max_age - min_age) + 0.5)
				++max_bins;
			if (age_bin >= max_bins)
				age_bin = max_bins - 1;
			float outcome_val = stof(tokens[3]);
			int outcome_ind = (int)distance(sorted_outcome_labels.begin(),
				find(sorted_outcome_labels.begin(), sorted_outcome_labels.end(), outcome_val));
			if (outcome_ind > sorted_outcome_labels.size())
				MTHROW_AND_ERR("Couldn't find outcome_value=%2.3f\n", outcome_val);
			if (tokens[1] == "MALE")
				male_labels_count_per_age[age_bin][outcome_ind] = stof(tokens[4]);
			else if (tokens[1] == "FEMALE")
				female_labels_count_per_age[age_bin][outcome_ind] = stof(tokens[4]);
			else
				MTHROW_AND_ERR("Unknown gender \"%s\"\n", tokens[1].c_str());
		}
		else
			MTHROW_AND_ERR("Unknown command \"%s\"\n", command.c_str());
	}
	sort(sorted_outcome_labels.begin(), sorted_outcome_labels.end());
	of.close();
}
void parse_vector(const string &value, vector<float> &output_vec) {
	vector<string> vec;
	boost::split(vec, value, boost::is_any_of(","));
	output_vec.resize((int)vec.size());
	for (size_t i = 0; i < vec.size(); ++i)
		output_vec[i] = stof(vec[i]);
}
ROC_Params::ROC_Params(const string &init_string) {
	max_diff_working_point = (float)0.05;
	use_score_working_points = false;
	working_point_FPR = { (float)0.1, 1, 5, 10,20,30,40,50,55,60,65,70,75,80,85,90,95 };
	score_bins = 0;
	incidence_fix = 0;

	//override default with given string:
	vector<string> tokens;
	boost::split(tokens, init_string, boost::is_any_of(";"));
	for (string token : tokens)
	{
		if (token.find('=') == string::npos)
			MTHROW_AND_ERR("Wrong token. has no value \"%s\"\n", token.c_str());
		string param_name = token.substr(0, token.find('='));
		string param_value = token.substr(token.find('=') + 1);
		boost::to_lower(param_name);

		if (param_name == "max_diff_working_point")
			max_diff_working_point = stof(param_value);
		else if (param_name == "use_score_working_points")
			use_score_working_points = stoi(param_value) > 0;
		else if (param_name == "score_bins")
			score_bins = stoi(param_value);
		else if (param_name == "inc_stats_text")
			inc_stats.read_from_text_file(param_value);
		else if (param_name == "inc_stats_bin")
			inc_stats.read_from_file(param_value);
		else if (param_name == "working_point_fpr")
			parse_vector(param_value, working_point_FPR);
		else if (param_name == "working_point_pr")
			parse_vector(param_value, working_point_PR);
		else if (param_name == "working_point_sens")
			parse_vector(param_value, working_point_SENS);
		else
			MTHROW_AND_ERR("Unknown paramter \"%s\" for ROC_Params\n", param_name.c_str());
	}
}
#pragma endregion
