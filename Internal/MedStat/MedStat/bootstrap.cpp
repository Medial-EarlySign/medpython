#include "bootstrap.h"
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <random>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fenv.h>
#ifndef  __unix__
#pragma float_control( except, on )
#endif

#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
//#define WARN_SKIP_WP
//#define USE_MIN_THREADS

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

random_device Lazy_Iterator::rd;

Lazy_Iterator::Lazy_Iterator(const vector<int> *p_pids, const vector<float> *p_preds,
	const vector<float> *p_y, float p_sample_ratio, int p_sample_per_pid, int max_loops, int seed) {
	sample_per_pid = p_sample_per_pid;
	//sample_ratio = p_sample_ratio;
	sample_ratio = 1.0; //no support for smaller size for now - need to fix Std for smaller sizes
	sample_all_no_sampling = false;
	pids = p_pids;
	y = p_y->data();
	preds = p_preds->data();
	maxThreadCount = max_loops;
	vec_size.resize(maxThreadCount);
	vec_y.resize(maxThreadCount);
	vec_preds.resize(maxThreadCount);
	vec_size.back() = (int)p_pids->size();
	vec_y.back() = y;
	vec_preds.back() = preds;

	unordered_map<int, vector<int>> pid_to_inds;
	for (size_t i = 0; i < pids->size(); ++i)
		pid_to_inds[(*pids)[i]].push_back(int(i));

	pid_index_to_indexes.resize((int)pid_to_inds.size());
	ind_to_pid.resize((int)pid_index_to_indexes.size());
	min_pid_start = INT_MAX;
	int max_pid_start = 0;
	int cnt_i = 0;
	int max_samples = 0;
	for (auto it = pid_to_inds.begin(); it != pid_to_inds.end(); ++it)
	{
		ind_to_pid[cnt_i] = it->first;
		if (max_samples < it->second.size())
			max_samples = (int)it->second.size();
		pid_index_to_indexes[cnt_i].swap(it->second);
		++cnt_i;
		if (it->first < min_pid_start)
			min_pid_start = it->first;
		if (it->first > max_pid_start)
			max_pid_start = it->first;
	}
	int rep_pid_size = max_pid_start - min_pid_start;
	cohort_size = int(sample_ratio * pid_index_to_indexes.size());
	//init:
	rd_gen.resize(maxThreadCount);
	for (size_t i = 0; i < maxThreadCount; ++i)
		if (seed == 0)
			rd_gen[i] = mt19937(rd());
		else
			rd_gen[i] = mt19937(seed);
	rand_pids = uniform_int_distribution<>(0, (int)pid_index_to_indexes.size() - 1);
	internal_random.resize(max_samples + 1);
	for (int i = 1; i <= max_samples; ++i)
		internal_random[i] = uniform_int_distribution<>(0, i - 1);
	//MLOG_D("created %d random gens\n", (int)internal_random.size());

	current_pos.resize(maxThreadCount, 0);
	inner_pos.resize(maxThreadCount, 0);
	sel_pid_index.resize(maxThreadCount, -1);
}

void Lazy_Iterator::set_static(const vector<float> *p_y, const vector<float> *p_preds, int thread_num) {
#pragma omp critical 
{
	vec_size[thread_num] = (int)p_y->size();
	vec_y[thread_num] = p_y->data();
	vec_preds[thread_num] = p_preds->data();
}
}

bool Lazy_Iterator::fetch_next(int thread, float &ret_y, float &ret_pred) {
	if (sample_per_pid > 0) {
		//choose pid:
		int selected_pid_index = int(current_pos[thread] / sample_per_pid);
		if (!sample_all_no_sampling)
			selected_pid_index = rand_pids(rd_gen[thread]);

		vector<int> *inds = &pid_index_to_indexes[selected_pid_index];
		uniform_int_distribution<> *rnd_num = &internal_random[inds->size()];

		int selected_index = (*inds)[(*rnd_num)(rd_gen[thread])];
		ret_y = y[selected_index];
		ret_pred = preds[selected_index];
#pragma omp atomic
		++current_pos[thread];
		return current_pos[thread] < sample_per_pid * cohort_size;
	}
	else { //taking all samples for pid when selected, sample_ratio is less than 1
		if (sample_all_no_sampling) {
			//iterate on all!:
			ret_y = vec_y[thread][current_pos[thread]];
			ret_pred = vec_preds[thread][current_pos[thread]];
#pragma omp atomic
			++current_pos[thread];
			return current_pos[thread] < vec_size[thread];
		}
		if (sel_pid_index[thread] < 0)
		{
			int selected_pid_index = rand_pids(rd_gen[thread]);
#pragma omp critical 
			{
				sel_pid_index[thread] = selected_pid_index;
				inner_pos[thread] = 0;
			}
		}
		vector<int> *inds = &pid_index_to_indexes[sel_pid_index[thread]];
		int final_index = (*inds)[inner_pos[thread]];
		ret_y = y[final_index];
		ret_pred = preds[final_index];
		//take all inds:
#pragma omp atomic
		++inner_pos[thread];
		if (inner_pos[thread] >= inds->size()) {
#pragma omp critical
		{
			sel_pid_index[thread] = -1;
			++current_pos[thread]; //mark pid as done
		}
		}
		return current_pos[thread] < cohort_size;
	}
}

bool Lazy_Iterator::fetch_next_external(int thread, float &ret_y, float &ret_pred) {
	return fetch_next(thread, ret_y, ret_pred);
}

void Lazy_Iterator::restart_iterator(int thread) {

	if (sample_ratio < 1) {
#pragma omp critical 
	{
		current_pos[thread] = 0;
		inner_pos[thread] = 0;
		sel_pid_index[thread] = -1;
	}
	}
	else {
#pragma omp critical 
		current_pos[thread] = 0;
	}

}

Lazy_Iterator::~Lazy_Iterator() {} //do nothing. nothing to clear

inline string format_working_point(const string &init_str, float wp, bool perc = true) {
	char res[100];
	if (perc)
		wp *= 100;
	snprintf(res, sizeof(res), "%s_%06.3f", init_str.c_str(), wp);
	return string(res);
}

template<typename T> inline int binary_search_position(const T *begin, const T *end, T val, bool reversed = false) {
	int maxSize = (int)(end - begin) + 1;
	int mid = int((maxSize - 1) / 2);
	if (maxSize <= 2) {
		if (!reversed) {
			if (val <= *begin) {
				return 0;
			}
			else if (val <= *end) {
				return 1;
			}
			else {
				return maxSize;
			}
		}
		else {
			if (val >= *begin) {
				return 0;
			}
			else if (val >= *end) {
				return 1;
			}
			else {
				return maxSize;
			}
		}
	}

	if (!reversed) {
		if (val <= begin[mid]) {
			return binary_search_position(begin, begin + mid, val, reversed);
		}
		else {
			return mid + binary_search_position(begin + mid, end, val, reversed);
		}
	}
	else {
		if (val >= begin[mid]) {
			return binary_search_position(begin, begin + mid, val, reversed);
		}
		else {
			return mid + binary_search_position(begin + mid, end, val, reversed);
		}
	}
}

template<typename T> inline int binary_search_position_last(const T *begin, const T *end, T val, bool reversed = false) {
	int maxSize = (int)(end - begin) + 1;
	int mid = int((maxSize - 1) / 2);
	if (maxSize <= 2) {
		if (!reversed) {
			if (val < *begin) {
				return 0;
			}
			else if (val < *end) {
				return 1;
			}
			else {
				return maxSize;
			}
		}
		else {
			if (val > *begin) {
				return 0;
			}
			else if (val > *end) {
				return 1;
			}
			else {
				return maxSize;
			}
		}
	}

	if (!reversed) {
		if (val < begin[mid]) {
			return binary_search_position_last(begin, begin + mid, val, reversed);
		}
		else {
			return mid + binary_search_position_last(begin + mid, end, val, reversed);
		}
	}
	else {
		if (val > begin[mid]) {
			return binary_search_position_last(begin, begin + mid, val, reversed);
		}
		else {
			return mid + binary_search_position_last(begin + mid, end, val, reversed);
		}
	}
}


#pragma endregion

int get_checksum(const vector<int> &pids) {
	int checksum = 0;
	for (int pid : pids)
		checksum = (checksum + pid) & 0xFFFF;
	return checksum;
}

map<string, float> booststrap_analyze_cohort(const vector<float> &preds, const vector<float> &y,
	const vector<int> &pids, float sample_ratio, int sample_per_pid, int loopCnt,
	const vector<MeasurementFunctions> &meas_functions, const vector<void *> &function_params,
	ProcessMeasurementParamFunc process_measurments_params,
	const map<string, vector<float>> &additional_info, const vector<float> &y_full,
	const vector<int> &pids_full, const vector<int> &filter_indexes, FilterCohortFunc cohort_def, void *cohort_params, int seed = 0) {
	//this function called after filter cohort
	//for each pid - randomize x sample from all it's tests. do loop_times
	float ci_bound = (float)0.95;

	//initialize measurement params per cohort:
	time_t st = time(NULL);
	for (size_t i = 0; i < function_params.size(); ++i)
		if (process_measurments_params != NULL && function_params[i] != NULL) {
			ROC_And_Filter_Params prm;
			prm.roc_params = (ROC_Params *)function_params[i];
			prm.filter = (vector<Filter_Param> *)cohort_params;
			process_measurments_params(additional_info, y, pids, &prm,
				filter_indexes, y_full, pids_full);
		}
	//MLOG_D("took %2.1f sec to process_measurments_params\n", (float)difftime(time(NULL), st));

#ifdef USE_MIN_THREADS
	Lazy_Iterator iterator(&pids, &preds, &y, sample_ratio, sample_per_pid, omp_get_max_threads(), seed); //for Obs
#else
	Lazy_Iterator iterator(&pids, &preds, &y, sample_ratio, sample_per_pid, loopCnt + 1, seed); //for Obs
#endif
	//MLOG_D("took %2.1f sec till allocate mem\n", (float)difftime(time(NULL), st));

	map<string, vector<float>> all_measures;
	iterator.sample_all_no_sampling = true;
	//iterator.sample_per_pid = 0; //take all samples in Obs
	//iterator.sample_ratio = 1; //take all pids
#ifdef USE_MIN_THREADS
	int main_thread = 0;
#else
	int main_thread = loopCnt;
#endif
	for (size_t k = 0; k < meas_functions.size(); ++k)
	{
		if (k > 0)
			iterator.restart_iterator(main_thread);
		map<string, float> batch_measures = meas_functions[k](&iterator, main_thread, function_params[k]);
		for (auto jt = batch_measures.begin(); jt != batch_measures.end(); ++jt)
			all_measures[jt->first + "_Obs"].push_back(jt->second);
	}
#ifdef USE_MIN_THREADS
	iterator.restart_iterator(0);
#endif

	if (sample_per_pid > 0) {
		//save results for all cohort:
		iterator.sample_all_no_sampling = false;
		//iterator.sample_per_pid = sample_per_pid;
		//iterator.sample_ratio = sample_ratio;

		Lazy_Iterator *iter_for_omp = &iterator;
#pragma omp parallel for schedule(static)
		for (int i = 0; i < loopCnt; ++i)
		{
#ifdef USE_MIN_THREADS
			int th_num = omp_get_thread_num();
#else
			int th_num = i;
#endif
			for (size_t k = 0; k < meas_functions.size(); ++k)
			{
#ifdef USE_MIN_THREADS
				iterator.restart_iterator(th_num);
#else
				if (k > 0)
					iterator.restart_iterator(th_num);
#endif
				map<string, float> batch_measures = meas_functions[k](iter_for_omp, th_num, function_params[k]);
#pragma omp critical
				for (auto jt = batch_measures.begin(); jt != batch_measures.end(); ++jt)
					all_measures[jt->first].push_back(jt->second);
			}
		}
	}
	else {
		//old implementition with memory:
		iterator.sample_all_no_sampling = true;

		mt19937 rd_gen;
		if (seed > 0)
			rd_gen = mt19937(seed);
		else {
			random_device rd;
			rd_gen = mt19937(rd());
		}
		unordered_map<int, vector<int>> pid_to_inds;
		for (size_t i = 0; i < pids.size(); ++i)
			pid_to_inds[pids[i]].push_back(int(i));

		int cohort_size = int(sample_ratio * pid_to_inds.size());
		int cnt_i = 0;
		vector<int> ind_to_pid((int)pid_to_inds.size());
		for (auto it = pid_to_inds.begin(); it != pid_to_inds.end(); ++it)
		{
			ind_to_pid[cnt_i] = it->first;
			++cnt_i;
		}
		//choose pids:
		uniform_int_distribution<> rand_pids(0, (int)pid_to_inds.size() - 1);

		//other sampling - sample pids and take all thier data:
		//now sample cohort 

#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < loopCnt; ++i)
		{
			vector<int> selected_pids(cohort_size);
			for (size_t k = 0; k < cohort_size; ++k)
			{
				int ind_pid = rand_pids(rd_gen);
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

			iterator.set_static(&selected_y, &selected_preds, i);
#ifdef USE_MIN_THREADS
			int th_num = omp_get_thread_num();
#else
			int th_num = i;
#endif
			//calc measures for sample:
			for (size_t k = 0; k < meas_functions.size(); ++k)
			{
				map<string, float> batch_measures;
#ifdef USE_MIN_THREADS
				iterator.restart_iterator(th_num);
#else
				if (k > 0)
					iterator.restart_iterator(th_num);
#endif
				batch_measures = meas_functions[k](&iterator, i, function_params[k]);
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
	all_final_measures["Checksum"] = (float)get_checksum(pids);
	//MLOG_D("took %2.1f sec to cohort\n", (float)difftime(time(NULL), st));

	return all_final_measures;
}

map<string, map<string, float>> booststrap_analyze(const vector<float> &preds, const vector<float> &y, const vector<int> &pids,
	const map<string, vector<float>> &additional_info, const map<string, FilterCohortFunc> &filter_cohort
	, const vector<MeasurementFunctions> &meas_functions, const map<string, void *> *cohort_params,
	const vector<void *> *function_params, ProcessMeasurementParamFunc process_measurments_params,
	PreprocessScoresFunc preprocess_scores, void *preprocess_scores_params, float sample_ratio, int sample_per_pid,
	int loopCnt, int seed, bool binary_outcome) {
#if defined(__unix__)
	//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
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
	const vector<float> *final_preds = &preds;
	vector<float> copy_preds;
	if (preprocess_scores != NULL) {
		copy_preds = vector<float>(preds);
		preprocess_scores(copy_preds, preprocess_scores_params);
		final_preds = &copy_preds;
	}

	map<string, map<string, float>> all_cohorts_measurments;
	vector<float> preds_c, y_c;
	vector<int> pids_c, filtered_indexes;
	vector<int> class_sz;
	pids_c.reserve((int)y.size());
	preds_c.reserve((int)y.size());
	filtered_indexes.reserve((int)y.size());
	y_c.reserve((int)y.size());
	for (auto it = filter_cohort.begin(); it != filter_cohort.end(); ++it)
	{
		void *c_params = NULL;
		if (cohort_params != NULL && (*cohort_params).find(it->first) != (*cohort_params).end())
			c_params = (*cohort_params).at(it->first);

		class_sz.resize(2, 0);
		class_sz[0] = 0, class_sz[1] = 0;
		pids_c.clear();
		preds_c.clear();
		y_c.clear();
		filtered_indexes.clear();
		for (size_t j = 0; j < y.size(); ++j)
			if (it->second(additional_info, (int)j, c_params)) {
				pids_c.push_back(pids[j]);
				y_c.push_back(y[j]);
				preds_c.push_back((*final_preds)[j]);
				filtered_indexes.push_back((int)j);
				++class_sz[y[j] > 0];
			}
		//now we have cohort: run analysis:
		string cohort_name = it->first;

		if (y_c.size() < 10) {
			MWARN("WARN: Cohort [%s] is too small - has %d samples. Skipping\n", cohort_name.c_str(), int(y_c.size()));
			continue;
		}
		if (binary_outcome) {
			if ((class_sz[0] < 1 || class_sz[1] < 1)) {
				MWARN("WARN: Cohort [%s] is too small - has %d samples with labels = [%d, %d]. Skipping\n",
					cohort_name.c_str(), int(y_c.size()), class_sz[0], class_sz[1]);
				continue;
			}
			else MLOG("Cohort [%s] - has %d samples with labels = [%d, %d]\n",
				cohort_name.c_str(), int(y_c.size()), class_sz[0], class_sz[1]);
		}
		map<string, float> cohort_measurments = booststrap_analyze_cohort(preds_c, y_c, pids_c,
			sample_ratio, sample_per_pid, loopCnt, meas_functions,
			function_params != NULL ? *function_params : params,
			process_measurments_params, additional_info, y, pids, filtered_indexes, it->second, c_params, seed);

		all_cohorts_measurments[cohort_name] = cohort_measurments;
	}
	MLOG_D("Finished Bootstarp Analysis. took %2.1f seconds\n", difftime(time(NULL), start));
	return all_cohorts_measurments;
}

void write_bootstrap_results(const string &file_name, const map<string, map<string, float>> &all_cohorts_measurments, const string& run_id) {
	string delimeter = "\t";
	if (all_cohorts_measurments.empty())
		throw invalid_argument("all_cohorts_measurments can't be empty");
	unordered_set<string> all_columns_uniq;
	for (auto jt = all_cohorts_measurments.begin(); jt != all_cohorts_measurments.end(); ++jt)
		for (auto it = jt->second.begin(); it != jt->second.end(); ++it)
			all_columns_uniq.insert(it->first);
	vector<string> all_columns(all_columns_uniq.begin(), all_columns_uniq.end());
	sort(all_columns.begin(), all_columns.end());
	ofstream fw(file_name);
	if (!fw.good())
		MTHROW_AND_ERR("IO Error: can't write \"%s\"\n", file_name.c_str());

	fw << "Cohort_Description";
	for (size_t i = 0; i < all_columns.size(); ++i)
		fw << delimeter << all_columns[i];
	if (!run_id.empty())
		fw << delimeter << "run_id";
	fw << endl;

	for (auto it = all_cohorts_measurments.begin(); it != all_cohorts_measurments.end(); ++it)
	{
		string cohort_name = it->first;
		map<string, float> cohort_values = it->second;
		fw << cohort_name;
		for (size_t i = 0; i < all_columns.size(); ++i)
			fw << delimeter <<
			(cohort_values.find(all_columns[i]) != cohort_values.end() ? cohort_values.at(all_columns[i]) : MED_MAT_MISSING_VALUE);
		if (!run_id.empty())
			fw << delimeter << run_id;
		fw << endl;
	}

	fw.close();
}
void read_bootstrap_results(const string &file_name, map<string, map<string, float>> &all_cohorts_measurments) {
	string delimeter = "\t";
	ifstream of(file_name);
	if (!of.good())
		MTHROW_AND_ERR("IO Error: can't read \"%s\"\n", file_name.c_str());
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

void write_pivot_bootstrap_results(const string &file_name, const map<string, map<string, float>> &all_cohorts_measurments, const string& run_id) {
	string delimeter = "\t";
	if (all_cohorts_measurments.empty())
		throw invalid_argument("all_cohorts_measurments can't be empty");
	map<string, float> flat_map;
	for (auto jt = all_cohorts_measurments.begin(); jt != all_cohorts_measurments.end(); ++jt) {
		char buff[1000];
		for (auto it = jt->second.begin(); it != jt->second.end(); ++it) {
			snprintf(buff, sizeof(buff), "%s%s%s", jt->first.c_str(), delimeter.c_str(), it->first.c_str());
			flat_map[string(buff)] = it->second;
		}
	}

	ofstream fw(file_name);
	if (!fw.good())
		MTHROW_AND_ERR("IO Error: can't write \"%s\"\n", file_name.c_str());

	fw << "Cohort" << delimeter << "Measurement" << delimeter << "Value" << endl;
	for (auto it = flat_map.begin(); it != flat_map.end(); ++it)
	{
		string cohort_measure_name = it->first;
		float value = it->second;
		fw << cohort_measure_name << delimeter << value << "\n";
	}
	if (!run_id.empty())
		for (auto jt = all_cohorts_measurments.begin(); jt != all_cohorts_measurments.end(); ++jt)
			fw << jt->first << delimeter << "run_id" << delimeter << run_id << "\n";

	fw.flush();
	fw.close();
}
void read_pivot_bootstrap_results(const string &file_name, map<string, map<string, float>> &all_cohorts_measurments) {
	string delimeter = "\t";
	if (all_cohorts_measurments.empty())
		throw invalid_argument("all_cohorts_measurments can't be empty");
	map<string, float> flat_map;

	ifstream fr(file_name);
	if (!fr.good())
		MTHROW_AND_ERR("IO Error: can't read \"%s\"\n", file_name.c_str());
	string line;
	getline(fr, line); //skip header
	while (getline(fr, line)) {
		boost::trim(line);
		if (line.empty())
			continue;
		vector<string> tokens;
		boost::split(tokens, line, boost::is_any_of(delimeter));
		if (tokens.size() != 3)
			MTHROW_AND_ERR("format error in line \"%s\"\n", line.c_str());
		string &cohort_name = tokens[0];
		string &measure_name = tokens[1];
		if (measure_name == "run_id")
			continue;
		float value = stof(tokens[2]);
		all_cohorts_measurments[cohort_name][measure_name] = value;
	}

	fr.close();
}

#pragma region Measurements Fucntions

map<string, float> calc_npos_nneg(Lazy_Iterator *iterator, int thread_num, void *function_params) {
	map<string, float> res;

	map<float, int> cnts;
	float y, pred;
	while (iterator->fetch_next(thread_num, y, pred))
		cnts[y] += 1;
	cnts[y] += 1; //last one

	res["NPOS"] = (float)cnts[(float)1.0];
	res["NNEG"] = (float)cnts[(float)0];

	return res;
}

map<string, float> calc_only_auc(Lazy_Iterator *iterator, int thread_num, void *function_params) {
	map<string, float> res;

	vector<float> pred_threshold;
	unordered_map<float, vector<float>> pred_to_labels;
	int tot_true_labels = 0;
	float y, pred;
	int tot_cnt = 0;
	while (iterator->fetch_next(thread_num, y, pred)) {
		pred_to_labels[pred].push_back(y);
		tot_true_labels += int(y > 0);
		++tot_cnt;
	}
	//last one
	pred_to_labels[pred].push_back(y);
	tot_true_labels += int(y > 0);
	++tot_cnt;

	int tot_false_labels = tot_cnt - tot_true_labels;
	if (tot_true_labels == 0 || tot_false_labels == 0)
		throw invalid_argument("only falses or positives exists in cohort");
	pred_threshold = vector<float>((int)pred_to_labels.size());
	unordered_map<float, vector<float>>::iterator it = pred_to_labels.begin();
	for (size_t i = 0; i < pred_threshold.size(); ++i)
	{
		pred_threshold[i] = it->first;
		++it;
	}
	sort(pred_threshold.begin(), pred_threshold.end());
	//From up to down sort:
	int t_cnt = 0;
	int f_cnt = 0;
	vector<float> true_rate = vector<float>((int)pred_to_labels.size());
	vector<float> false_rate = vector<float>((int)pred_to_labels.size());
	int st_size = (int)pred_threshold.size() - 1;
	for (int i = st_size; i >= 0; --i)
	{
		vector<float> *indexes = &pred_to_labels[pred_threshold[i]];
		//calc AUC status for this step:
		for (float y : *indexes)
		{
			bool true_label = y > 0;
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

map<string, float> calc_roc_measures_with_inc(Lazy_Iterator *iterator, int thread_num, void *function_params) {
	map<string, float> res;
	int max_qunt_vals = 10;
	bool censor_removed = true;
	bool trunc_max = false;

	ROC_Params *params = (ROC_Params *)function_params;
	float max_diff_in_wp = params->max_diff_working_point;
	int scores_bin = params->score_bins;

	vector<float> fpr_points = params->working_point_FPR; //Working FPR points:
	sort(fpr_points.begin(), fpr_points.end());
	for (size_t i = 0; i < fpr_points.size(); ++i)
		fpr_points[i] /= 100.0;
	vector<float> sens_points = params->working_point_SENS; //Working SENS points:
	sort(sens_points.begin(), sens_points.end());
	for (size_t i = 0; i < sens_points.size(); ++i)
		sens_points[i] /= 100.0;
	vector<float> pr_points = params->working_point_PR; //Working PR points:
	sort(pr_points.begin(), pr_points.end());
	for (size_t i = 0; i < pr_points.size(); ++i)
		pr_points[i] /= 100.0;

	unordered_map<float, vector<float>> thresholds_labels;
	vector<float> unique_scores;
	float y, pred;
	while (iterator->fetch_next(thread_num, y, pred))
		thresholds_labels[pred].push_back(y);
	thresholds_labels[pred].push_back(y); //last one

	unique_scores.resize((int)thresholds_labels.size());
	int ind_p = 0;
	for (auto it = thresholds_labels.begin(); it != thresholds_labels.end(); ++it)
	{
		unique_scores[ind_p] = it->first;
		++ind_p;
	}
	sort(unique_scores.begin(), unique_scores.end());

	//calc measures on each bucket of scores as possible threshold:
	double t_sum = 0, f_sum = 0, tt_cnt = 0;
	int f_cnt = 0;
	int t_cnt = 0;
	vector<float> true_rate((int)unique_scores.size());
	vector<float> false_rate((int)unique_scores.size());
	int st_size = (int)unique_scores.size() - 1;
	for (int i = st_size; i >= 0; --i)
	{
		vector<float> *labels = &thresholds_labels[unique_scores[i]];
		for (float y : *labels)
		{
			float true_label = params->fix_label_to_binary ? y > 0 : y;
			t_sum += true_label;
			tt_cnt += true_label > 0 ? true_label : 0;
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

	if (f_cnt == 0 || t_sum <= 0) {
		MWARN("no falses or no positives exists in cohort\n");
		return res;
	}
	for (size_t i = 0; i < true_rate.size(); ++i) {
		true_rate[i] /= float(!trunc_max ? t_sum : tt_cnt);
		false_rate[i] /= float(f_sum);
	}
	//calc maesures based on true_rate and false_rate
	double auc = false_rate[0] * true_rate[0] / 2; //"auc" on expectitions:
	for (size_t i = 1; i < true_rate.size(); ++i)
		auc += (false_rate[i] - false_rate[i - 1]) * (true_rate[i - 1] + true_rate[i]) / 2;

	bool use_wp = unique_scores.size() > max_qunt_vals && !params->use_score_working_points; //change all working points
	int curr_wp_fpr_ind = 0, curr_wp_sens_ind = 0, curr_wp_pr_ind = 0;
	int i = 0;

	float ppv_c, pr_prev, ppv_prev, pr_c, npv_c, npv_prev, or_prev, or_c, rr_prev, rr_c;
	if (use_wp) {
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
					res[format_working_point("SCORE@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("SENS@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("PR@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("PPV@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("NPV@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("OR@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("LIFT@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("RR@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT FPR=%f, prev_FPR=%f, next_FPR=%f, prev_score=%f, next_score=%f\n",
						fpr_points[curr_wp_fpr_ind], false_rate[i - 1], false_rate[i],
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif
					++curr_wp_fpr_ind;
					continue; //skip working point - diff is too big
				}
				res[format_working_point("SCORE@FPR", fpr_points[curr_wp_fpr_ind])] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				res[format_working_point("SENS@FPR", fpr_points[curr_wp_fpr_ind])] = 100 * (true_rate[i] * (prev_diff / tot_diff) +
					true_rate[i - 1] * (curr_diff / tot_diff));
				if (params->incidence_fix > 0) {
					ppv_c = float(params->incidence_fix*true_rate[i] / (params->incidence_fix*true_rate[i] + (1 - params->incidence_fix)*false_rate[i]));
					if (true_rate[i - 1] > 0 || false_rate[i - 1] > 0)
						ppv_prev = float(params->incidence_fix*true_rate[i - 1] / (params->incidence_fix*true_rate[i - 1] + (1 - params->incidence_fix)*false_rate[i - 1]));
					else
						ppv_prev = ppv_c;
				}
				else {
					ppv_c = float((true_rate[i] * t_sum) /
						((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
					if (true_rate[i - 1] > 0 || false_rate[i - 1] > 0)
						ppv_prev = float((true_rate[i - 1] * t_sum) /
							((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
					else
						ppv_prev = ppv_c;
				}
				float ppv = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				res[format_working_point("PPV@FPR", fpr_points[curr_wp_fpr_ind])] = 100 * ppv;
				if (params->incidence_fix > 0) {
					pr_c = float(params->incidence_fix*true_rate[i] + (1 - params->incidence_fix)*false_rate[i]);
					pr_prev = float(params->incidence_fix*true_rate[i - 1] + (1 - params->incidence_fix)*false_rate[i - 1]);
				}
				else {
					pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
						(t_sum + f_sum));
					pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
						(t_sum + f_sum));
				}
				res[format_working_point("PR@FPR", fpr_points[curr_wp_fpr_ind])] = 100 * (pr_c* (prev_diff / tot_diff) + pr_prev * (curr_diff / tot_diff));
				if (params->incidence_fix > 0) {
					npv_prev = float(((1 - false_rate[i - 1]) *  (1 - params->incidence_fix)) /
						(((1 - true_rate[i - 1]) *  params->incidence_fix) + ((1 - false_rate[i - 1]) *  (1 - params->incidence_fix))));
					if (true_rate[i] < 1 || false_rate[i] < 1)
						npv_c = float(((1 - false_rate[i]) *  (1 - params->incidence_fix)) /
							(((1 - true_rate[i]) *  params->incidence_fix) + ((1 - false_rate[i]) *  (1 - params->incidence_fix))));
					else
						npv_c = npv_prev;
				}
				else {
					npv_prev = float(((1 - false_rate[i - 1]) * f_sum) /
						(((1 - true_rate[i - 1]) * t_sum) + ((1 - false_rate[i - 1]) * f_sum)));
					if (true_rate[i] < 1 || false_rate[i] < 1)
						npv_c = float(((1 - false_rate[i]) * f_sum) /
							(((1 - true_rate[i]) * t_sum) + ((1 - false_rate[i]) * f_sum)));
					else
						npv_c = npv_prev;
				}
				res[format_working_point("NPV@FPR", fpr_points[curr_wp_fpr_ind])] = 100 * (npv_c * (prev_diff / tot_diff) + npv_prev*(curr_diff / tot_diff));
				if (params->incidence_fix > 0)
					res[format_working_point("LIFT@FPR", fpr_points[curr_wp_fpr_ind])] = float(ppv / params->incidence_fix);
				else
					res[format_working_point("LIFT@FPR", fpr_points[curr_wp_fpr_ind])] = float(ppv /
						(t_sum / (t_sum + f_sum))); //lift of prevalance when there is no inc

				if (false_rate[i] > 0 && false_rate[i] < 1 && true_rate[i] < 1)
					or_c = float(
						(true_rate[i] / false_rate[i]) / ((1 - true_rate[i]) / (1 - false_rate[i])));
				else
					or_c = MED_MAT_MISSING_VALUE;
				if (false_rate[i - 1] > 0 && false_rate[i - 1] < 1 && true_rate[i - 1] < 1)
					or_prev = float(
						(true_rate[i - 1] / false_rate[i - 1]) / ((1 - true_rate[i - 1]) / (1 - false_rate[i - 1])));
				else
					or_prev = MED_MAT_MISSING_VALUE;
				if (or_c != MED_MAT_MISSING_VALUE && or_prev != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@FPR", fpr_points[curr_wp_fpr_ind])] = (or_c * (prev_diff / tot_diff) +
						or_prev * (curr_diff / tot_diff));
				else if (or_c != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@FPR", fpr_points[curr_wp_fpr_ind])] = or_c;
				else if (or_prev != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@FPR", fpr_points[curr_wp_fpr_ind])] = or_prev;
				else
					res[format_working_point("OR@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;

				if (params->incidence_fix > 0) {
					if (true_rate[i - 1] < 1)
						rr_prev = float(ppv_prev + ppv_prev * (1 - params->incidence_fix)* (1 - false_rate[i - 1]) /
							(params->incidence_fix * (1 - true_rate[i - 1])));
					else
						rr_prev = MED_MAT_MISSING_VALUE;

					if (true_rate[i] < 1)
						rr_c = float(ppv_c + ppv_c * (1 - params->incidence_fix)* (1 - false_rate[i]) /
							(params->incidence_fix * (1 - true_rate[i])));
					else
						rr_c = MED_MAT_MISSING_VALUE;
					if (rr_c != MED_MAT_MISSING_VALUE && rr_prev != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@FPR", fpr_points[curr_wp_fpr_ind])] = (rr_c * (prev_diff / tot_diff) +
							rr_prev * (curr_diff / tot_diff));
					else if (rr_c != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@FPR", fpr_points[curr_wp_fpr_ind])] = rr_c;
					else if (rr_prev != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@FPR", fpr_points[curr_wp_fpr_ind])] = rr_prev;
					else
						res[format_working_point("RR@FPR", fpr_points[curr_wp_fpr_ind])] = MED_MAT_MISSING_VALUE;
				}

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
					res[format_working_point("SCORE@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("FPR@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("SPEC@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("PR@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("PPV@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("NPV@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("OR@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("RR@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("LIFT@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT SENS=%f, prev_SENS=%f, next_SENS=%f, prev_score=%f, next_score=%f\n",
						sens_points[curr_wp_sens_ind], true_rate[i - 1], true_rate[i],
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif
					++curr_wp_sens_ind;
					continue; //skip working point - diff is too big
				}
				res[format_working_point("SCORE@SENS", sens_points[curr_wp_sens_ind])] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				res[format_working_point("FPR@SENS", sens_points[curr_wp_sens_ind])] = 100 * (false_rate[i] * (prev_diff / tot_diff) +
					false_rate[i - 1] * (curr_diff / tot_diff));
				res[format_working_point("SPEC@SENS", sens_points[curr_wp_sens_ind])] = 100 * ((1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff));
				if (params->incidence_fix > 0) {
					ppv_c = float(params->incidence_fix*true_rate[i] / (params->incidence_fix*true_rate[i] + (1 - params->incidence_fix)*false_rate[i]));
					if (true_rate[i - 1] > 0 || false_rate[i - 1] > 0)
						ppv_prev = float(params->incidence_fix*true_rate[i - 1] / (params->incidence_fix*true_rate[i - 1] + (1 - params->incidence_fix)*false_rate[i - 1]));
					else
						ppv_prev = ppv_c;
				}
				else {
					ppv_c = float((true_rate[i] * t_sum) /
						((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
					if (true_rate[i - 1] > 0 || false_rate[i - 1] > 0)
						ppv_prev = float((true_rate[i - 1] * t_sum) /
							((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
					else
						ppv_prev = ppv_c;
				}
				float ppv = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				res[format_working_point("PPV@SENS", sens_points[curr_wp_sens_ind])] = 100 * ppv;
				if (params->incidence_fix > 0) {
					pr_c = float(params->incidence_fix*true_rate[i] + (1 - params->incidence_fix)*false_rate[i]);
					pr_prev = float(params->incidence_fix*true_rate[i - 1] + (1 - params->incidence_fix)*false_rate[i - 1]);
				}
				else {
					pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
						(t_sum + f_sum));
					pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
						(t_sum + f_sum));
				}
				res[format_working_point("PR@SENS", sens_points[curr_wp_sens_ind])] = 100 * (pr_c* (prev_diff / tot_diff) + pr_prev * (curr_diff / tot_diff));
				if (params->incidence_fix > 0) {
					npv_prev = float(((1 - false_rate[i - 1]) *  (1 - params->incidence_fix)) /
						(((1 - true_rate[i - 1]) *  params->incidence_fix) + ((1 - false_rate[i - 1]) *  (1 - params->incidence_fix))));
					if (true_rate[i] < 1 || false_rate[i] < 1)
						npv_c = float(((1 - false_rate[i]) *  (1 - params->incidence_fix)) /
							(((1 - true_rate[i]) *  params->incidence_fix) + ((1 - false_rate[i]) *  (1 - params->incidence_fix))));
					else
						npv_c = npv_prev;
				}
				else {
					npv_prev = float(((1 - false_rate[i - 1]) * f_sum) /
						(((1 - true_rate[i - 1]) * t_sum) + ((1 - false_rate[i - 1]) * f_sum)));
					if (true_rate[i] < 1 || false_rate[i] < 1)
						npv_c = float(((1 - false_rate[i]) * f_sum) /
							(((1 - true_rate[i]) * t_sum) + ((1 - false_rate[i]) * f_sum)));
					else
						npv_c = npv_prev;
				}
				res[format_working_point("NPV@SENS", sens_points[curr_wp_sens_ind])] = 100 * (npv_c * (prev_diff / tot_diff) + npv_prev*(curr_diff / tot_diff));
				if (params->incidence_fix > 0)
					res[format_working_point("LIFT@SENS", sens_points[curr_wp_sens_ind])] = float(ppv / params->incidence_fix);
				else
					res[format_working_point("LIFT@SENS", sens_points[curr_wp_sens_ind])] = float(ppv /
						(t_sum / (t_sum + f_sum))); //lift of prevalance when there is no inc

				if (false_rate[i] > 0 && false_rate[i] < 1 && true_rate[i] < 1)
					or_c = float(
						(true_rate[i] / false_rate[i]) / ((1 - true_rate[i]) / (1 - false_rate[i])));
				else
					or_c = MED_MAT_MISSING_VALUE;
				if (false_rate[i - 1] > 0 && false_rate[i - 1] < 1 && true_rate[i - 1] < 1)
					or_prev = float(
						(true_rate[i - 1] / false_rate[i - 1]) / ((1 - true_rate[i - 1]) / (1 - false_rate[i - 1])));
				else
					or_prev = MED_MAT_MISSING_VALUE;
				if (or_c != MED_MAT_MISSING_VALUE && or_prev != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@SENS", sens_points[curr_wp_sens_ind])] = (or_c * (prev_diff / tot_diff) +
						or_prev * (curr_diff / tot_diff));
				else if (or_c != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@SENS", sens_points[curr_wp_sens_ind])] = or_c;
				else if (or_prev != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@SENS", sens_points[curr_wp_sens_ind])] = or_prev;
				else
					res[format_working_point("OR@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;

				if (params->incidence_fix > 0) {
					if (true_rate[i - 1] < 1)
						rr_prev = float(ppv_prev + ppv_prev * (1 - params->incidence_fix)* (1 - false_rate[i - 1]) /
							(params->incidence_fix * (1 - true_rate[i - 1])));
					else
						rr_prev = MED_MAT_MISSING_VALUE;

					if (true_rate[i] < 1)
						rr_c = float(ppv_c + ppv_c * (1 - params->incidence_fix)* (1 - false_rate[i]) /
							(params->incidence_fix * (1 - true_rate[i])));
					else
						rr_c = MED_MAT_MISSING_VALUE;
					if (rr_c != MED_MAT_MISSING_VALUE && rr_prev != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@SENS", sens_points[curr_wp_sens_ind])] = (rr_c * (prev_diff / tot_diff) +
							rr_prev * (curr_diff / tot_diff));
					else if (rr_c != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@SENS", sens_points[curr_wp_sens_ind])] = rr_c;
					else if (rr_prev != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@SENS", sens_points[curr_wp_sens_ind])] = rr_prev;
					else
						res[format_working_point("RR@SENS", sens_points[curr_wp_sens_ind])] = MED_MAT_MISSING_VALUE;
				}

				++curr_wp_sens_ind;
				continue;
			}
			++i;
		}

		//handle pr points:
		i = 1; //first point is always before
		while (i < true_rate.size() && curr_wp_pr_ind < pr_points.size())
		{
			if (params->incidence_fix > 0)
				pr_c = float(params->incidence_fix*true_rate[i] + (1 - params->incidence_fix)*false_rate[i]);
			else
				pr_c = float(((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
					(t_sum + f_sum));

			if (curr_wp_pr_ind < pr_points.size() && pr_c >= pr_points[curr_wp_pr_ind]) { //passed work_point - take 2 last points for measure - by distance from wp
				if (params->incidence_fix > 0)
					pr_prev = float(params->incidence_fix*true_rate[i - 1] + (1 - params->incidence_fix)*false_rate[i - 1]);
				else
					pr_prev = float(((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)) /
						(t_sum + f_sum));

				float prev_diff = pr_points[curr_wp_pr_ind] - pr_prev;
				float curr_diff = pr_c - pr_points[curr_wp_pr_ind];
				float tot_diff = prev_diff + curr_diff;
				if (tot_diff <= 0) {
					curr_diff = 1;
					tot_diff = 1; //take prev - first apeareance
				}
				if (prev_diff > max_diff_in_wp || curr_diff > max_diff_in_wp) {
					res[format_working_point("SCORE@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("FPR@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("SPEC@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("SENS@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("PPV@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("NPV@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("OR@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("RR@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
					res[format_working_point("LIFT@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
#ifdef  WARN_SKIP_WP
					MWARN("SKIP WORKING POINT PR=%f, prev_PR=%f, next_PR=%f, prev_score=%f, next_score=%f\n",
						pr_points[curr_wp_pr_ind], pr_prev, pr_c,
						pred_threshold[st_size - (i - 1)], pred_threshold[st_size - i]);
#endif //  WARN_SKIP_WP
					++curr_wp_pr_ind;
					continue; //skip working point - diff is too big
				}
				res[format_working_point("SCORE@PR", pr_points[curr_wp_pr_ind])] = unique_scores[st_size - i] * (prev_diff / tot_diff) +
					unique_scores[st_size - (i - 1)] * (curr_diff / tot_diff);
				res[format_working_point("FPR@PR", pr_points[curr_wp_pr_ind])] = 100 * (false_rate[i] * (prev_diff / tot_diff) +
					false_rate[i - 1] * (curr_diff / tot_diff));
				res[format_working_point("SPEC@PR", pr_points[curr_wp_pr_ind])] = 100 * ((1 - false_rate[i]) * (prev_diff / tot_diff) +
					(1 - false_rate[i - 1]) * (curr_diff / tot_diff));
				if (params->incidence_fix > 0) {
					ppv_c = float(params->incidence_fix*true_rate[i] / (params->incidence_fix*true_rate[i] + (1 - params->incidence_fix)*false_rate[i]));
					if (false_rate[i - 1] > 0 || true_rate[i - 1] > 0)
						ppv_prev = float(params->incidence_fix*true_rate[i - 1] / (params->incidence_fix*true_rate[i - 1] + (1 - params->incidence_fix)*false_rate[i - 1]));
					else
						ppv_prev = ppv_c;
				}
				else {
					ppv_c = float((true_rate[i] * t_sum) /
						((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
					if (false_rate[i - 1] > 0 || true_rate[i - 1] > 0)
						ppv_prev = float((true_rate[i - 1] * t_sum) /
							((true_rate[i - 1] * t_sum) + (false_rate[i - 1] * f_sum)));
					else
						ppv_prev = ppv_c;
				}
				float ppv = ppv_c * (prev_diff / tot_diff) + ppv_prev*(curr_diff / tot_diff);
				res[format_working_point("PPV@PR", pr_points[curr_wp_pr_ind])] = 100 * ppv;
				res[format_working_point("SENS@PR", pr_points[curr_wp_pr_ind])] = 100 * (true_rate[i] * (prev_diff / tot_diff) + true_rate[i - 1] * (curr_diff / tot_diff));
				if (params->incidence_fix > 0) {
					npv_prev = float(((1 - false_rate[i - 1]) *  (1 - params->incidence_fix)) /
						(((1 - true_rate[i - 1]) *  params->incidence_fix) + ((1 - false_rate[i - 1]) *  (1 - params->incidence_fix))));
					if (true_rate[i] < 1 || false_rate[i] < 1)
						npv_c = float(((1 - false_rate[i]) *  (1 - params->incidence_fix)) /
							(((1 - true_rate[i]) *  params->incidence_fix) + ((1 - false_rate[i]) *  (1 - params->incidence_fix))));
					else
						npv_c = npv_prev;
				}
				else {
					npv_prev = float(((1 - false_rate[i - 1]) * f_sum) /
						(((1 - true_rate[i - 1]) * t_sum) + ((1 - false_rate[i - 1]) * f_sum)));
					if (true_rate[i] < 1 || false_rate[i] < 1)
						npv_c = float(((1 - false_rate[i]) * f_sum) /
							(((1 - true_rate[i]) * t_sum) + ((1 - false_rate[i]) * f_sum)));
					else
						npv_c = npv_prev;
				}
				res[format_working_point("NPV@PR", pr_points[curr_wp_pr_ind])] = 100 * (npv_c * (prev_diff / tot_diff) + npv_prev*(curr_diff / tot_diff));
				if (params->incidence_fix > 0)
					res[format_working_point("LIFT@PR", pr_points[curr_wp_pr_ind])] = float(ppv / params->incidence_fix);
				else
					res[format_working_point("LIFT@PR", pr_points[curr_wp_pr_ind])] = float(ppv /
						(t_sum / (t_sum + f_sum))); //lift of prevalance when there is no inc
				if (false_rate[i] > 0 && false_rate[i] < 1 && true_rate[i] < 1)
					or_c = float(
						(true_rate[i] / false_rate[i]) / ((1 - true_rate[i]) / (1 - false_rate[i])));
				else
					or_c = MED_MAT_MISSING_VALUE;
				if (false_rate[i - 1] > 0 && false_rate[i - 1] < 1 && true_rate[i - 1] < 1)
					or_prev = float(
						(true_rate[i - 1] / false_rate[i - 1]) / ((1 - true_rate[i - 1]) / (1 - false_rate[i - 1])));
				else
					or_prev = MED_MAT_MISSING_VALUE;
				if (or_c != MED_MAT_MISSING_VALUE && or_prev != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@PR", pr_points[curr_wp_pr_ind])] = (or_c * (prev_diff / tot_diff) +
						or_prev * (curr_diff / tot_diff));
				else if (or_c != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@PR", pr_points[curr_wp_pr_ind])] = or_c;
				else if (or_prev != MED_MAT_MISSING_VALUE)
					res[format_working_point("OR@PR", pr_points[curr_wp_pr_ind])] = or_prev;
				else
					res[format_working_point("OR@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;

				if (params->incidence_fix > 0) {
					if (true_rate[i - 1] < 1)
						rr_prev = float(ppv_prev + ppv_prev * (1 - params->incidence_fix)* (1 - false_rate[i - 1]) /
							(params->incidence_fix * (1 - true_rate[i - 1])));
					else
						rr_prev = MED_MAT_MISSING_VALUE;

					if (true_rate[i] < 1)
						rr_c = float(ppv_c + ppv_c * (1 - params->incidence_fix)* (1 - false_rate[i]) /
							(params->incidence_fix * (1 - true_rate[i])));
					else
						rr_c = MED_MAT_MISSING_VALUE;
					if (rr_c != MED_MAT_MISSING_VALUE && rr_prev != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@PR", pr_points[curr_wp_pr_ind])] = (rr_c * (prev_diff / tot_diff) +
							rr_prev * (curr_diff / tot_diff));
					else if (rr_c != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@PR", pr_points[curr_wp_pr_ind])] = rr_c;
					else if (rr_prev != MED_MAT_MISSING_VALUE)
						res[format_working_point("RR@PR", pr_points[curr_wp_pr_ind])] = rr_prev;
					else
						res[format_working_point("RR@PR", pr_points[curr_wp_pr_ind])] = MED_MAT_MISSING_VALUE;
				}

				++curr_wp_pr_ind;
				continue;
			}
			++i;
		}

	}
	else {
		float score_working_point;
		for (i = 0; i < true_rate.size(); ++i)
		{
			score_working_point = unique_scores[st_size - i];
			res[format_working_point("SENS@SCORE", score_working_point, false)] = 100 * true_rate[i];
			res[format_working_point("SPEC@SCORE", score_working_point, false)] = 100 * (1 - false_rate[i]);
			float ppv = MED_MAT_MISSING_VALUE;
			if (true_rate[i] > 0 || false_rate[i] > 0) {
				if (params->incidence_fix > 0)
					ppv = float((true_rate[i] * params->incidence_fix) /
						(params->incidence_fix*(true_rate[i]) +
							(false_rate[i] * (1 - params->incidence_fix))));
				else
					ppv = float((true_rate[i] * t_sum) /
						((true_rate[i] * t_sum) + (false_rate[i] * f_sum)));
				res[format_working_point("PPV@SCORE", score_working_point, false)] = 100 * ppv;
			}
			else
				res[format_working_point("PPV@SCORE", score_working_point, false)] = MED_MAT_MISSING_VALUE;

			if (params->incidence_fix > 0) {
				res[format_working_point("PR@SCORE", score_working_point, false)] = float(100 * ((true_rate[i] * params->incidence_fix) + (false_rate[i] * (1 - params->incidence_fix))));
			}
			else {
				res[format_working_point("PR@SCORE", score_working_point, false)] = float(100 * ((true_rate[i] * t_sum) + (false_rate[i] * f_sum)) /
					(t_sum + f_sum));
			}
			if (true_rate[i] < 1 || false_rate[i] < 1) {
				if (params->incidence_fix > 0) {
					res[format_working_point("NPV@SCORE", score_working_point, false)] = float(100 * ((1 - false_rate[i]) * (1 - params->incidence_fix)) /
						(((1 - true_rate[i]) * params->incidence_fix) + ((1 - false_rate[i]) *  (1 - params->incidence_fix))));
				}
				else {
					res[format_working_point("NPV@SCORE", score_working_point, false)] = float(100 * ((1 - false_rate[i]) * f_sum) /
						(((1 - true_rate[i]) * t_sum) + ((1 - false_rate[i]) * f_sum)));
				}
			}
			else
				res[format_working_point("NPV@SCORE", score_working_point, false)] = MED_MAT_MISSING_VALUE;
			if (params->incidence_fix > 0) {
				if (true_rate[i] > 0 || false_rate[i] > 0 || ppv == MED_MAT_MISSING_VALUE)
					res[format_working_point("LIFT@SCORE", score_working_point, false)] = float(ppv / params->incidence_fix);
				else
					res[format_working_point("LIFT@SCORE", score_working_point, false)] = MED_MAT_MISSING_VALUE;
			}
			else {
				if (true_rate[i] > 0 || false_rate[i] > 0 || ppv == MED_MAT_MISSING_VALUE)
					res[format_working_point("LIFT@SCORE", score_working_point, false)] = float(ppv /
						(t_sum / (t_sum + f_sum)));
				else
					res[format_working_point("LIFT@SCORE", score_working_point, false)] = MED_MAT_MISSING_VALUE;
			}
			if (false_rate[i] > 0 && false_rate[i] < 1 && true_rate[i] < 1)
				res[format_working_point("OR@SCORE", score_working_point, false)] = float(
					(true_rate[i] / false_rate[i]) / ((1 - true_rate[i]) / (1 - false_rate[i])));
			else
				res[format_working_point("OR@SCORE", score_working_point, false)] = MED_MAT_MISSING_VALUE;

			if (params->incidence_fix > 0) {
				if (true_rate[i] < 1 || ppv == MED_MAT_MISSING_VALUE)
					res[format_working_point("RR@SCORE", score_working_point, false)] = float((ppv + ppv * (1 - params->incidence_fix)* (1 - false_rate[i]) /
						(params->incidence_fix * (1 - true_rate[i]))));
				else
					res[format_working_point("RR@SCORE", score_working_point, false)] = MED_MAT_MISSING_VALUE;
			}
			else {
				if (true_rate[i] < 1 || ppv == MED_MAT_MISSING_VALUE) {
					float FOR = float(((1.0 - true_rate[i]) * t_sum) /
						((1 - true_rate[i]) * t_sum + (1 - false_rate[i]) * f_sum));
					res[format_working_point("RR@SCORE", score_working_point, false)] =
						float(ppv / FOR);
				}
				else
					res[format_working_point("RR@SCORE", score_working_point, false)] = MED_MAT_MISSING_VALUE;
			}
		}
	}


	res["AUC"] = float(auc);

	if (abs(t_cnt - t_sum) > 0.01) {
		res["NEG_SUM"] = float(f_sum);
		res["POS_SUM"] = float(t_sum);
		res["POS_CNT"] = float(t_cnt);
		res["NEG_CNT"] = float(f_cnt);
	}
	else {
		res["NNEG"] = float(f_sum);
		res["NPOS"] = float(t_sum);
	}

	return res;
}

map<string, float> calc_kandel_tau(Lazy_Iterator *iterator, int thread_num, void *function_params) {
	map<string, float> res;

	double tau = 0, cnt = 0;
	float y, pred;
	//vector<float> scores, labels;
	unordered_map<float, vector<float>> label_to_scores;
	while (iterator->fetch_next(thread_num, y, pred))
		label_to_scores[y].push_back(pred);
	label_to_scores[y].push_back(pred);// last one
	for (auto it = label_to_scores.begin(); it != label_to_scores.end(); ++it)
		sort(it->second.begin(), it->second.end());

	for (auto it = label_to_scores.begin(); it != label_to_scores.end(); ++it)
	{
		auto bg = it;
		++bg;
		vector<float> *preds = &it->second;
		int pred_i_bigger;
		double pred_i_smaller;
		for (auto jt = bg; jt != label_to_scores.end(); ++jt)
		{
			vector<float> *preds_comp = &jt->second;
			double p_size = (double)preds_comp->size();
			for (float pred : *preds)
			{
				pred_i_bigger = binary_search_position(preds_comp->data(), preds_comp->data() + preds_comp->size() - 1, pred);
				pred_i_smaller = p_size - binary_search_position_last(preds_comp->data(), preds_comp->data() + preds_comp->size() - 1, pred);
				if (it->first > jt->first)
					//tau += pred_i_bigger;
					tau += pred_i_bigger - pred_i_smaller;
				else
					//tau += pred_i_smaller;
					tau += pred_i_smaller - pred_i_bigger;
			}
			cnt += p_size * preds->size();
		}
	}

	if (cnt > 1) {
		tau /= cnt;
		res["Kendall-Tau"] = (float)tau;
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
		(outcome <= 0 && diff_days >= max_time));
}

bool filter_range_param(const map<string, vector<float>> &record_info, int index, void *cohort_params) {
	Filter_Param *param = (Filter_Param *)cohort_params; //can't be null
	if (param->param_name != "Time-Window")
		return record_info.at(param->param_name)[index] >= param->min_range &&
		record_info.at(param->param_name)[index] <= param->max_range;
	else
		return time_range_filter(record_info.at("Label")[index] > 0, param->min_range,
			param->max_range, record_info.at(param->param_name)[index]);
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
			res = time_range_filter(record_info.at("Label")[index] > 0, (*param)[i].min_range,
				(*param)[i].max_range, record_info.at((*param)[i].param_name)[index]);
		++i;
	}
	return res;
}
#pragma endregion

#pragma region Process Measurement Param Functions
void count_stats(int bin_counts, const vector<float> &y, const map<string, vector<float>> &additional_info
	, const vector<int> &filtered_indexes, const ROC_Params *params,
	vector<vector<double>> &male_counts, vector<vector<double>> &female_counts) {

	male_counts.resize(bin_counts);
	female_counts.resize(bin_counts);
	for (size_t i = 0; i < male_counts.size(); ++i)
		male_counts[i].resize(2);
	for (size_t i = 0; i < female_counts.size(); ++i)
		female_counts[i].resize(2);
	int min_age = (int)params->inc_stats.min_age;
	int max_age = (int)params->inc_stats.max_age;
	//if filtered_indexes is empty pass on all y. otherwise traverse over indexes:
	if (filtered_indexes.empty()) {
		for (size_t i = 0; i < y.size(); ++i)
		{
			if (additional_info.at("Age")[i] < min_age ||
				additional_info.at("Age")[i] >= max_age + params->inc_stats.age_bin_years)
				continue; //skip out of range or already case
			int age_index = (int)floor((additional_info.at("Age")[i] - min_age) /
				params->inc_stats.age_bin_years);
			if (age_index >= bin_counts)
				age_index = bin_counts - 1;

			if (additional_info.at("Gender")[i] == GENDER_MALE)  //Male
				++male_counts[age_index][y[i] > 0];
			else //Female
				++female_counts[age_index][y[i] > 0];
		}
	}
	else {
		for (size_t ii = 0; ii < filtered_indexes.size(); ++ii)
		{
			int i = filtered_indexes[ii];
			if (additional_info.at("Age")[i] < min_age ||
				additional_info.at("Age")[i] >= max_age + params->inc_stats.age_bin_years)
				continue; //skip out of range or already case
			int age_index = (int)floor((additional_info.at("Age")[i] - min_age) /
				params->inc_stats.age_bin_years);
			if (age_index >= bin_counts)
				age_index = bin_counts - 1;

			if (additional_info.at("Gender")[i] == GENDER_MALE)  //Male
				++male_counts[age_index][y[i] > 0];
			else //Female
				++female_counts[age_index][y[i] > 0];
		}
	}
}

void fix_cohort_sample_incidence(const map<string, vector<float>> &additional_info,
	const vector<float> &y, const vector<int> &pids, void *function_params,
	const vector<int> &filtered_indexes, const vector<float> &y_full, const vector<int> &pids_full) {
	ROC_And_Filter_Params *pr_full = (ROC_And_Filter_Params *)function_params;
	ROC_Params *params = pr_full->roc_params;
	vector<Filter_Param> *cohort_filt = pr_full->filter;

	if (params->inc_stats.sorted_outcome_labels.empty())
		return; //no inc file
	//calculating the "fixed" incidence in the cohort giving the true inc. in the general population
	// select cohort - and multiply in the given original incidence
	if (params->inc_stats.sorted_outcome_labels.size() != 2)
		MTHROW_AND_ERR("Category outcome aren't supported for now\n");
	if (additional_info.find("Age") == additional_info.end() || additional_info.find("Gender") == additional_info.end())
		MTHROW_AND_ERR("Age or Gender Signals are missings\n");

	int min_age = (int)params->inc_stats.min_age;
	int max_age = (int)params->inc_stats.max_age;
	int bin_counts = (int)floor((max_age - min_age) / params->inc_stats.age_bin_years);
	if (bin_counts * params->inc_stats.age_bin_years >
		(max_age - min_age) + 0.5)
		++bin_counts; //has at least 0.5 years for last bin to create it
	if (params->inc_stats.male_labels_count_per_age.size() != bin_counts)
		MTHROW_AND_ERR("Male vector has %d members. and need to have %d members\n",
			(int)params->inc_stats.male_labels_count_per_age.size(), bin_counts);

	vector<vector<double>> filtered_male_counts, filtered_female_counts;
	vector<vector<double>> all_male_counts, all_female_counts;
	count_stats(bin_counts, y_full, additional_info, filtered_indexes, params,
		filtered_male_counts, filtered_female_counts);
	//always filter Time-Window:
	vector<int> baseline_all;
	Filter_Param *time_win_cond = NULL;
	for (auto i = 0; i < cohort_filt->size() && time_win_cond == NULL; ++i)
		if ((*cohort_filt)[i].param_name == "Time-Window")
			time_win_cond = &(*cohort_filt)[i];
	if (time_win_cond != NULL)
		for (size_t i = 0; i < y_full.size(); ++i)
			if (filter_range_param(additional_info, (int)i, time_win_cond))
				baseline_all.push_back((int)i);
	count_stats(bin_counts, y_full, additional_info, baseline_all, params,
		all_male_counts, all_female_counts);

	//Lets calc the ratio for the incidence in the filter:
	params->incidence_fix = 0;
	double tot_controls = 0;
	//recalc new ratio of #1/(#1+#0) and fix stats
	for (size_t i = 0; i < bin_counts; ++i)
	{
		double tot_cn_filtered = filtered_male_counts[i][0] + filtered_male_counts[i][1];
		if (filtered_male_counts[i][0] > 0 && all_male_counts[i][1] > 0) {
			double general_inc = params->inc_stats.male_labels_count_per_age[i][1] /
				(params->inc_stats.male_labels_count_per_age[i][1] +
					params->inc_stats.male_labels_count_per_age[i][0]);
			tot_controls += filtered_male_counts[i][0];
			double filtered_inc = filtered_male_counts[i][1] / tot_cn_filtered;
			double all_inc = all_male_counts[i][1] / (all_male_counts[i][1] + all_male_counts[i][0]);
			double inc_ratio = filtered_inc / all_inc;

			params->incidence_fix += filtered_male_counts[i][0] * general_inc * inc_ratio;
		}
		tot_cn_filtered = filtered_female_counts[i][0] + filtered_female_counts[i][1];
		if (filtered_female_counts[i][0] > 0 && all_female_counts[i][1]) {
			double general_inc = params->inc_stats.female_labels_count_per_age[i][1] /
				(params->inc_stats.female_labels_count_per_age[i][1] +
					params->inc_stats.female_labels_count_per_age[i][0]);
			tot_controls += filtered_female_counts[i][0];
			double filtered_inc = filtered_female_counts[i][1] / tot_cn_filtered;
			double all_inc = all_female_counts[i][1] / (all_female_counts[i][1] + all_female_counts[i][0]);
			double inc_ratio = filtered_inc / all_inc;

			params->incidence_fix += filtered_female_counts[i][0] * general_inc * inc_ratio;
		}
	}

	if (tot_controls > 0)
		params->incidence_fix /= tot_controls;

	MLOG_D("Running fix_cohort_sample_incidence and got %2.4f%% mean incidence\n",
		100 * params->incidence_fix);
}

void fix_cohort_sample_incidence_old(const map<string, vector<float>> &additional_info,
	const vector<float> &y, const vector<int> &pids, void *function_params,
	const vector<int> &filtered_indexes, const vector<float> &y_full, const vector<int> &pids_full) {
	ROC_And_Filter_Params *pr_full = (ROC_And_Filter_Params *)function_params;
	ROC_Params *params = pr_full->roc_params;
	if (params->inc_stats.sorted_outcome_labels.empty())
		return; //no inc file
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

	vector<vector<double>> filtered_male_counts, filtered_female_counts;
	count_stats(bin_counts, y_full, additional_info, filtered_indexes, params,
		filtered_male_counts, filtered_female_counts);

	params->incidence_fix = 0;
	double tot_controls = 0;
	//recalc new ratio of #1/(#1+#0) and fix stats
	for (size_t i = 0; i < bin_counts; ++i)
	{
		//Males:
		if (filtered_male_counts[i][0] > 0) {
			double general_inc = params->inc_stats.male_labels_count_per_age[i][1] /
				(params->inc_stats.male_labels_count_per_age[i][1] +
					params->inc_stats.male_labels_count_per_age[i][0]);
			tot_controls += filtered_male_counts[i][0];
			params->incidence_fix += filtered_male_counts[i][0] * general_inc;
		}
		//Females:
		if (filtered_female_counts[i][0] > 0) {
			double general_inc = params->inc_stats.female_labels_count_per_age[i][1] /
				(params->inc_stats.female_labels_count_per_age[i][1] +
					params->inc_stats.female_labels_count_per_age[i][0]);
			tot_controls += filtered_female_counts[i][0];
			params->incidence_fix += filtered_female_counts[i][0] * general_inc;
		}
	}

	if (tot_controls > 0)
		params->incidence_fix /= tot_controls;

	MLOG_D("Running fix_cohort_sample_incidence and got %2.4f%% mean incidence\n",
		100 * params->incidence_fix);
}
#pragma endregion

#pragma region Process Scores Functions
void _simple_find(const vector<pair<int, int>> &vec, int &found_pos, int search_pos) {
	found_pos = -1;
	for (int j = (int)vec.size() - 1; j >= 0 && found_pos == -1; --j)
		if (vec[j].second >= search_pos && vec[j].first <= search_pos)
			found_pos = j;
}

void merge_down(vector<int> &ind_to_size, vector<vector<pair<int, int>>> &size_to_ind, set<int> &sizes,
	const pair<int, int> *index_to_merge) {
	pair<int, int> *merge_into = NULL;
	int to_merge_size = ind_to_size[index_to_merge->first - 1];
	int erase_index = -1;
	//remove index_to_merge.first - 1:
	_simple_find(size_to_ind[to_merge_size], erase_index, index_to_merge->first - 1);
	if (erase_index == -1)
		MTHROW_AND_ERR("down: Bug couldn't found merge_into\n");
	merge_into = &size_to_ind[to_merge_size][erase_index];

	int new_size = *sizes.begin() + to_merge_size;
	sizes.insert(new_size);
	//update in min,max:
	ind_to_size[merge_into->first] = new_size;
	ind_to_size[index_to_merge->second] = new_size;
	ind_to_size[merge_into->second] = new_size;
	ind_to_size[index_to_merge->first] = new_size;
	//erase old one
	int first_pos = merge_into->first;
	int second_pos = index_to_merge->second;
	//already popd merged_into element
	size_to_ind[to_merge_size].erase(size_to_ind[to_merge_size].begin() + erase_index);

	//insert new union
	size_to_ind[new_size].push_back(pair<int, int>(first_pos, second_pos));
}
void merge_up(vector<int> &ind_to_size, vector<vector<pair<int, int>>> &size_to_ind, set<int> &sizes,
	const pair<int, int> *index_to_merge) {
	//merge with +1
	pair<int, int> *merge_into = NULL;
	int to_merge_size = ind_to_size[index_to_merge->second + 1];
	int erase_index = -1;
	//remove index_to_merge.second + 1:
	_simple_find(size_to_ind[to_merge_size], erase_index, index_to_merge->second + 1);
	if (erase_index == -1)
		MTHROW_AND_ERR("up: Bug couldn't found merge_into\n");
	merge_into = &size_to_ind[to_merge_size][erase_index];

	int new_size = *sizes.begin() + to_merge_size;
	sizes.insert(new_size);
	//update in min,max:
	ind_to_size[index_to_merge->first] = new_size;
	ind_to_size[merge_into->second] = new_size;
	ind_to_size[index_to_merge->second] = new_size;
	ind_to_size[merge_into->first] = new_size;
	//erase old one:
	int first_pos = index_to_merge->first;
	int second_pos = merge_into->second;
	//already popd merged_into element
	size_to_ind[to_merge_size].erase(size_to_ind[to_merge_size].begin() + erase_index);

	//insert new union set:
	size_to_ind[new_size].push_back(pair<int, int>(first_pos, second_pos));
}

void preprocess_bin_scores(vector<float> &preds, void *function_params) {
	ROC_Params params;
	if (function_params != NULL)
		params = *(ROC_Params *)function_params;
	else
		return;

	if (params.use_score_working_points)
		return;

	if (params.score_resolution != 0)
		for (size_t i = 0; i < preds.size(); ++i)
			preds[i] = (float)round((double)preds[i] / params.score_resolution) *
			params.score_resolution;

	unordered_map<float, vector<int>> thresholds_indexes;
	vector<float> unique_scores;
	for (size_t i = 0; i < preds.size(); ++i)
		thresholds_indexes[preds[i]].push_back((int)i);
	unique_scores.resize((int)thresholds_indexes.size());
	int ind_p = 0, min_size = -1;
	for (auto it = thresholds_indexes.begin(); it != thresholds_indexes.end(); ++it)
	{
		unique_scores[ind_p] = it->first;
		++ind_p;
		if (min_size == -1 || min_size > it->second.size())
			min_size = (int)it->second.size();
	}
	sort(unique_scores.begin(), unique_scores.end());
	int bin_size_last = (int)thresholds_indexes.size();
	if (params.score_bins > 0 && bin_size_last < 10)
		if (params.score_resolution != 0)
			MWARN("Warnning Bootstrap:: requested specific working points, but score vector"
				" is highly quantitize(%d). try canceling preprocess_score by "
				"score_resolution, score_bins. Will use score working points\n",
				bin_size_last);
		else
			MWARN("Warnning Bootstrap:: requested specific working points, but score vector"
				" is highly quantitize(%d). Will use score working points\n",
				bin_size_last);

	if ((params.score_bins > 0 && bin_size_last > params.score_bins) ||
		(params.score_min_samples > 0 && min_size < params.score_min_samples)) {
		int c = 0;
		vector<vector<pair<int, int>>> size_to_ind(preds.size()); //size, group, index_min_max
		vector<int> ind_to_size(bin_size_last);
		set<int> sizes;
		for (auto it = unique_scores.begin(); it != unique_scores.end(); ++it)
		{
			size_to_ind[(int)thresholds_indexes[*it].size()].push_back(pair<int, int>(c, c));
			ind_to_size[c] = (int)thresholds_indexes[*it].size();
			++c;
			sizes.insert((int)thresholds_indexes[*it].size());
		}

		while ((params.score_bins > 0 && bin_size_last > params.score_bins)
			|| (params.score_min_samples > 0 && *sizes.begin() < params.score_min_samples)) {
			min_size = *sizes.begin();
			if (size_to_ind[min_size].empty())
				MTHROW_AND_ERR("Bug couldn't found min_size=%d\n", min_size);

			pair<int, int> index_to_merge = size_to_ind[min_size].back();
			size_to_ind[min_size].pop_back(); //now popback
			pair<int, int> *merge_into = NULL;
			//merge index_to_merge with index_to_merge+-1. and update size_to_ind, ind_to_size, sizes
			if (index_to_merge.second == unique_scores.size() - 1)
				merge_down(ind_to_size, size_to_ind, sizes, &index_to_merge);
			else if (index_to_merge.first == 0)
				merge_up(ind_to_size, size_to_ind, sizes, &index_to_merge);
			else {
				//MLOG("DEBUG: %d,%d\n", index_to_merge.first, index_to_merge.second);
				if (ind_to_size[index_to_merge.second + 1] < ind_to_size[index_to_merge.first - 1])
					merge_up(ind_to_size, size_to_ind, sizes, &index_to_merge);
				else
					merge_down(ind_to_size, size_to_ind, sizes, &index_to_merge);
			}

			while (size_to_ind[min_size].empty()) {//erase if left empty after merge
				sizes.erase(sizes.begin());
				min_size = *sizes.begin();
			}
			--bin_size_last;
		}

		//update thresholds_indexes based on: size_to_ind groups -
		//merge all indexes in each group to first index in thresholds_indexes. "mean" other scores to unique_scores
		unordered_set<float> u_scores;
		for (auto it = sizes.begin(); it != sizes.end(); ++it)
		{
			for (size_t k = 0; k < size_to_ind[*it].size(); ++k)
			{ //merge from first => second
				pair<int, int> *merge = &size_to_ind[*it][k];
				double mean_score = 0, tot_cnt = 0;
				vector<int> merged_inds;
				for (int ii = merge->first; ii <= merge->second; ++ii) {
					mean_score += unique_scores[ii] * thresholds_indexes[unique_scores[ii]].size();
					tot_cnt += thresholds_indexes[unique_scores[ii]].size();
					merged_inds.insert(merged_inds.end(),
						thresholds_indexes[unique_scores[ii]].begin(), thresholds_indexes[unique_scores[ii]].end());
				}
				mean_score /= tot_cnt;
				//update all preds to mean_score in merged_inds:
				for (int ind : merged_inds)
					preds[ind] = (float)mean_score;
				u_scores.insert((float)mean_score);
			}
		}
		if (u_scores.size() < 10) {
			MWARN("Warnning Bootstrap:: requested specific working points, but score vector"
				" is highly quantitize(%d). try canceling preprocess_score by "
				"score_resolution, score_bins. Will use score working points\n",
				(int)u_scores.size());
		}
	}

	MLOG_D("Preprocess_bin_scores Done - left with %d bins!\n", bin_size_last);
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
	if (!fw.good())
		MTHROW_AND_ERR("IO Error: can't write \"%s\"\n", text_file.c_str());
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
	MLOG("Loading Incidence file %s\n", text_file.c_str());
	ifstream of(text_file);
	if (!of.good())
		MTHROW_AND_ERR("IO Error: can't read \"%s\"\n", text_file.c_str());
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
		else if (command == "AGE_MIN") {
			min_age = stof(tokens[1]);
		}
		else if (command == "AGE_MAX")
			max_age = stof(tokens[1]);
		else if (command == "OUTCOME_VALUE")
			sorted_outcome_labels.push_back(stof(tokens[1]));
		else if (command == "STATS_ROW") {
			if (tokens.size() != 5)
				MTHROW_AND_ERR("Unknown lines format \"%s\"\n", line.c_str());
			float age = stof(tokens[2]);
			min_age = float(int(min_age / age_bin_years) * age_bin_years);
			max_age = float(int(ceil(max_age / age_bin_years)) * age_bin_years);

			if (age < min_age || age> max_age) {
				MWARN("Warning:: skip age because out of range in line \"%s\"\n", line.c_str());
				continue;
			}
			int age_bin = (int)floor((age - min_age) / age_bin_years);
			int max_bins = (int)floor((max_age - min_age) / age_bin_years);
			if (max_bins * age_bin_years > (max_age - min_age) + 0.5)
				++max_bins;
			if (age_bin >= max_bins)
				age_bin = max_bins - 1;
			if (male_labels_count_per_age.empty()) {
				male_labels_count_per_age.resize(max_bins);
				for (size_t i = 0; i < male_labels_count_per_age.size(); ++i)
					male_labels_count_per_age[i].resize((int)sorted_outcome_labels.size());
			}
			if (female_labels_count_per_age.empty()) {
				female_labels_count_per_age.resize(max_bins);
				for (size_t i = 0; i < female_labels_count_per_age.size(); ++i)
					female_labels_count_per_age[i].resize((int)sorted_outcome_labels.size());
			}
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
	score_resolution = 0;
	incidence_fix = 0;
	score_min_samples = 0;
	fix_label_to_binary = true;

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
		else if (param_name == "fix_label_to_binary")
			fix_label_to_binary = stoi(param_value) > 0;
		else if (param_name == "score_bins")
			score_bins = stoi(param_value);
		else if (param_name == "score_min_samples")
			score_min_samples = stoi(param_value);
		else if (param_name == "score_resolution")
			score_resolution = stof(param_value);
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
