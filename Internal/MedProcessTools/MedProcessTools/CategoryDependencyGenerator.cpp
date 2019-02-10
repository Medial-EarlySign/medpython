#include "FeatureGenerator.h"
#include <cmath>
#include <regex>
#include <MedUtils/MedUtils/MedRegistry.h>
#include <omp.h>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void pid_data_vec::add_data_point(int time, int val) {
	times.push_back(time);
	values.push_back(val);
}

pid_data_vec::pid_data_vec() {
	pid = -1; byear = -1; gender = -1;
}

pid_data_vec::pid_data_vec(int p, int b, int g) {
	pid = p; byear = b; gender = g;
}

int pid_data_vec::get_time(int idx) const {
	return times[idx];
}
int pid_data_vec::get_val(int idx) const {
	return values[idx];
}
int pid_data_vec::get_age(int idx) const {
	return med_time_converter.convert_times(global_default_time_unit, MedTime::Date, times[idx]) / 10000 - byear;
}

void CategoryDependencyGenerator::init_defaults() {
	generator_type = FTR_GEN_BASIC;
	signalName = "";
	signalId = -1;
	time_channel = 0;
	val_channel = 0;
	win_from = 0;
	win_to = 360000;
	time_unit_win = MedTime::Days;
	regex_filter = "";
	min_age = 0;
	max_age = 120;
	age_bin = 5;
	min_code_cnt = 0;
	fdr = (float)0.01;
	take_top = 50;
	lift_below = (float)0.9;
	lift_above = (float)1.1;
	chi_square_at_least = 0;
	minimal_chi_cnt = 5;
	stat_metric = category_stat_test::mcnemar;
	max_depth = 0;
	max_parents = 1;
	verbose = false;

	req_signals = { "BYEAR", "GENDER" };
}

static unordered_map<string, int> conv_map_stats = {
	{ "chi-square", (int)category_stat_test::chi_square },
	{ "mcnemar", (int)category_stat_test::mcnemar }
};
int CategoryDependencyGenerator::init(map<string, string>& mapper) {

	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [CategoryDependencyGenerator::init]
		if (it->first == "signal")
			signalName = it->second;
		else if (it->first == "val_channel")
			val_channel = med_stoi(it->second);
		else if (it->first == "time_channel")
			time_channel = med_stoi(it->second);
		else if (it->first == "win_from")
			win_from = med_stoi(it->second);
		else if (it->first == "win_to")
			win_to = med_stoi(it->second);
		else if (it->first == "time_unit_win")
			time_unit_win = med_time_converter.string_to_type(it->second);
		else if (it->first == "regex_filter")
			regex_filter = it->second;
		else if (it->first == "min_age")
			min_age = med_stoi(it->second);
		else if (it->first == "max_age")
			max_age = med_stoi(it->second);
		else if (it->first == "age_bin")
			age_bin = med_stoi(it->second);
		else if (it->first == "min_code_cnt")
			min_code_cnt = med_stoi(it->second);
		else if (it->first == "fdr")
			fdr = med_stof(it->second);
		else if (it->first == "take_top")
			take_top = med_stoi(it->second);
		else if (it->first == "lift_below")
			lift_below = med_stof(it->second);
		else if (it->first == "lift_above")
			lift_above = med_stof(it->second);
		else if (it->first == "chi_square_at_least")
			chi_square_at_least = med_stof(it->second);
		else if (it->first == "minimal_chi_cnt")
			minimal_chi_cnt = med_stoi(it->second);
		else if (it->first == "verbose")
			verbose = med_stoi(it->second) > 0;
		else if (it->first == "stat_metric") {
			if (conv_map_stats.find(it->second) != conv_map_stats.end())
				stat_metric = category_stat_test(conv_map_stats.at(it->second));
			else
				MTHROW_AND_ERR("Unknown stat_test \"%s\". options are: %s\n",
					it->second.c_str(), medial::io::get_list(conv_map_stats).c_str());
		}
		else if (it->first == "max_depth")
			max_depth = med_stoi(it->second);
		else if (it->first == "max_parents")
			max_parents = med_stoi(it->second);
		else if (it->first == "fg_type") {}
		else
			MTHROW_AND_ERR("Unknown parameter \'%s\' for CategoryDependencyGenerator\n", it->first.c_str())
			//! [CategoryDependencyGenerator::init]
	}

	if (signalName.empty())
		MTHROW_AND_ERR("Erorr in CategoryDependencyGenerator - Must provide signal\n");

	req_signals.push_back(signalName);

	return 0;
}

void CategoryDependencyGenerator::set_signal_ids(MedDictionarySections& dict) {
	signalId = dict.id(signalName);
	byear_sid = dict.id("BYEAR");
	gender_sid = dict.id("GENDER");
}

void CategoryDependencyGenerator::init_tables(MedDictionarySections& dict) {
	int section_id = dict.section_id(signalName);
	categoryId_to_name = dict.dict(section_id)->Id2Names;
	_member2Sets = dict.dict(section_id)->Member2Sets;

	luts.resize(top_codes.size());
	for (size_t i = 0; i < top_codes.size(); ++i) {
		vector<string> s_names = { top_codes[i] };
		dict.prep_sets_lookup_table(section_id, s_names, luts[i]);
	}
}

int CategoryDependencyGenerator::_learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *> processors) {
	if (signalId == -1 || byear_sid == -1 || gender_sid == -1)
		MTHROW_AND_ERR("Uninitialized signalId,byear_sid or gender_sid - or not loaded\n");

	// Required signals
	vector<int> all_req_signal_ids_v;
	vector<unordered_set<int> > current_required_signal_ids(processors.size());
	vector<FeatureGenerator *> generators = { this };
	unordered_set<int> extra_req_signal_ids;
	handle_required_signals(processors, generators, extra_req_signal_ids, all_req_signal_ids_v, current_required_signal_ids);

	// Collect Data
	PidDynamicRec rec;
	UniversalSigVec usv;
	for (unsigned int i = 0; i < ids.size(); i++) {
		int pid = ids[i];
		int gender = medial::repository::get_value(rep, pid, gender_sid);
		int byear = medial::repository::get_value(rep, pid, byear_sid);
		pid_data_vec p_data(pid, byear, gender);

		if (rep.sigs.Sid2Info[signalId].virtual_sig == 0)
			rep.uget(pid, signalId, usv);
		else {
			vector<int> time_points;
			rec.init_from_rep(std::addressof(rep), pid, all_req_signal_ids_v, 1);

			// Apply Processors
			vector<vector<float>> dummy_attributes_mat;
			for (unsigned int i = 0; i < processors.size(); i++)
				processors[i]->conditional_apply(rec, time_points, current_required_signal_ids[i], dummy_attributes_mat);

			// Collect values and ages
			rec.uget(signalId, 0, usv);
		}
		for (int k = 0; k < usv.len; ++k)
			p_data.add_data_point(usv.Time(k, time_channel), usv.Val(k, val_channel));
		pids_data[pid] = p_data;
	}


	return 0;
}

template<class T> void apply_filter(vector<int> &indexes, const vector<T> &vecCnts
	, double min_val, double max_val) {
	vector<int> filtered_indexes;
	filtered_indexes.reserve(indexes.size());
	for (size_t i = 0; i < indexes.size(); ++i)
		if (vecCnts[indexes[i]] >= min_val && vecCnts[indexes[i]] <= max_val)
			filtered_indexes.push_back(indexes[i]);
	indexes.swap(filtered_indexes);
}

int _count_legal_rows(const  vector<vector<int>> &m, int minimal_balls) {
	int res = 0;
	for (auto i = 0; i < m.size(); ++i)
	{
		int ind = 0;
		bool all_good = true;
		while (all_good && ind < m[i].size()) {
			all_good = m[i][ind] >= minimal_balls;
			++ind;
		}
		res += int(all_good);
	}
	return res;
}

void CategoryDependencyGenerator::get_parents(int codeGroup, vector<int> &parents) const {
	vector<int> last_parents = { codeGroup };
	if (last_parents.front() < 0)
		return; //no parents
	parents = { codeGroup };

	for (size_t k = 0; k < max_depth; ++k) {
		vector<int> new_layer;
		for (int par : last_parents)
			if (_member2Sets.find(par) != _member2Sets.end()) {
				new_layer.insert(new_layer.end(), _member2Sets.at(par).begin(), _member2Sets.at(par).end());
				parents.insert(parents.end(), _member2Sets.at(par).begin(), _member2Sets.at(par).end()); //aggregate all parents
			}
		if (parents.size() >= max_parents)
			break;
		new_layer.swap(last_parents);
		if (last_parents.empty())
			break; //no more parents to loop up
	}
}

void CategoryDependencyGenerator::get_stats(const unordered_map<int, vector<vector<vector<int>>>> &categoryVal_to_stats,
	vector<int> &all_signal_values, vector<int> &signal_indexes, vector<double> &valCnts,
	vector<double> &posCnts, vector<double> &lift, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio) {

	unordered_set<int> all_vals;
	for (auto i = categoryVal_to_stats.begin(); i != categoryVal_to_stats.end(); ++i)
		all_vals.insert(i->first);
	all_signal_values.insert(all_signal_values.end(), all_vals.begin(), all_vals.end());
	signal_indexes.resize(all_signal_values.size());
	for (int i = 0; i < signal_indexes.size(); ++i)
		signal_indexes[i] = i;
	posCnts.resize(all_signal_values.size());
	valCnts.resize(all_signal_values.size());
	lift.resize(all_signal_values.size());
	scores.resize(all_signal_values.size());
	p_values.resize(all_signal_values.size());
	pos_ratio.resize(all_signal_values.size());

	for (int index : signal_indexes)
	{
		int signalVal = all_signal_values[index];
		//check chi-square for this value:
		double totCnt = 0;
		double sig_sum = 0;
		double sum_noSig_reg = 0;
		double sum_noSig_tot = 0;

		const vector<vector<vector<int>>> &code_stats = categoryVal_to_stats.at(signalVal);
		for (size_t i = 0; i < code_stats.size(); ++i) //iterate genders
			for (size_t j = 0; j < code_stats[i].size(); ++j) //iterate age
				if (!code_stats[i][j].empty()) {
					totCnt += code_stats[i][j][2] + code_stats[i][j][3];
					posCnts[index] += code_stats[i][j][1 + 2];
					sig_sum += code_stats[i][j][0 + 2];
					sum_noSig_reg += code_stats[i][j][1 + 0];
					sum_noSig_tot += code_stats[i][j][1 + 0] + code_stats[i][j][0 + 0];
				}

		if (totCnt == 0)
			continue;
		valCnts[index] = totCnt; //for signal apeareance
		sig_sum += posCnts[index];
		if (sig_sum > 0 && sum_noSig_reg > 0)
			lift[index] = (posCnts[index] / sig_sum) / (sum_noSig_reg / sum_noSig_tot);
		if (sig_sum > 0 && sum_noSig_reg <= 0)
			lift[index] = 2 * posCnts[index]; //maximum lift
		pos_ratio[index] = posCnts[index] / totCnt;

		double regScore = 0;
		int dof = -1;
		for (size_t i = 0; i < code_stats.size(); ++i) {//iterate genders
			map<float, vector<int>> all_grps;
			for (size_t j = 0; j < code_stats[i].size(); ++j)
				if (!code_stats[i][j].empty())
					all_grps[float(j)] = code_stats[i][j];
			if (stat_metric == category_stat_test::mcnemar)
				regScore += medial::contingency_tables::calc_mcnemar_square_dist(all_grps);
			else if (stat_metric == category_stat_test::chi_square)
				regScore += medial::contingency_tables::calc_chi_square_dist(all_grps, 0, chi_square_at_least, minimal_chi_cnt);
			dof += _count_legal_rows(code_stats[i], stat_metric == category_stat_test::chi_square ? minimal_chi_cnt : 0);
		}
		scores[index] = (float)regScore;

		double pv = medial::contingency_tables::chisqr(dof, regScore);
		p_values[index] = pv;
	}
}

void CategoryDependencyGenerator::post_learn_from_samples(MedPidRepository& rep, const MedSamples& samples) {
	unordered_map<int, vector<vector<vector<int>>>> categoryVal_to_stats; //stats is gender,age, 4 ints counts:
	int age_bin_cnt = int(ceil((max_age - min_age + 1) / float(age_bin)));
	regex reg_pat;
	if (!regex_filter.empty())
		reg_pat = regex(regex_filter);

	vector<vector<vector<int>>> total_stats(2); //gender, age, label
	for (size_t i = 0; i < 2; ++i) {
		total_stats[i].resize(age_bin_cnt);
		for (size_t j = 0; j < age_bin_cnt; ++j)
			total_stats[i][j].resize(2);
	}

	MedTimer tm;
	tm.start();
	chrono::high_resolution_clock::time_point tm_prog = chrono::high_resolution_clock::now();
	int progress = 0;
	//unordered_map<int, vector<vector<bool>>> pid_label_age_bin;// stores for each pid if saw label,age_bin
	unordered_map<int, unordered_map<int, vector<vector<bool>>>> code_pid_label_age_bin;// stores for each code => pid if saw label,age_bin
	bool nested_state = omp_get_nested();
	omp_set_nested(true);
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < samples.idSamples.size(); ++i)
	{
		vector<vector<bool>> p_lbl_age;
		p_lbl_age.resize(2);
		for (size_t k = 0; k < 2; ++k)
			p_lbl_age[k].resize(age_bin_cnt);
		for (size_t j = 0; j < samples.idSamples[i].samples.size(); ++j)
		{
			if (pids_data.find(samples.idSamples[i].samples[j].id) == pids_data.end())
				MTHROW_AND_ERR("Error CategoryDependencyGenerator::post_learn_from_samples - couldn't find pid %d in post_learn\n",
					samples.idSamples[i].samples[j].id);
			/*vector<vector<bool>> &p_lbl_age = pid_label_age_bin[samples.idSamples[i].samples[j].id];
			if (p_lbl_age.empty()) {
				p_lbl_age.resize(2);
				for (size_t k = 0; k < 2; ++k)
					p_lbl_age[k].resize(age_bin_cnt);
			}*/
			const pid_data_vec &p_d = pids_data.at(samples.idSamples[i].samples[j].id);
			int byear = p_d.byear;
			int age = med_time_converter.convert_times(samples.time_unit, MedTime::Date, samples.idSamples[i].samples[j].time) / 10000 - byear;
			if (age > max_age || age < min_age)
				continue;
			int gend_idx = pids_data.at(samples.idSamples[i].samples[j].id).gender - 1;
			int outcome_idx = samples.idSamples[i].samples[j].outcome > 0;
			int age_idx = (age - min_age) / age_bin;
			if (!p_lbl_age[outcome_idx][age_idx]) {
				p_lbl_age[outcome_idx][age_idx] = true;
#pragma omp atomic
				++total_stats[gend_idx][age_idx][outcome_idx];
			}

			//allocate stats per category value:

			int start_time_win = med_time_converter.convert_times(time_unit_win, samples.time_unit, med_time_converter.convert_times(samples.time_unit, time_unit_win, samples.idSamples[i].samples[j].time) - win_to);
			int end_time_win = med_time_converter.convert_times(time_unit_win, samples.time_unit, med_time_converter.convert_times(samples.time_unit, time_unit_win, samples.idSamples[i].samples[j].time) - win_from);

			for (int k = 0; k < p_d.times.size(); ++k)
				if (p_d.get_time(k) >= start_time_win && p_d.get_time(k) <= end_time_win) { //get values in time window:
					//get filter regex codes:
					int base_code = p_d.get_val(k);
					vector<int> all_codes;
					get_parents(base_code, all_codes);
					for (int code : all_codes)
					{
						bool pass_regex_filter = regex_filter.empty();
						if (!pass_regex_filter) {
							if (categoryId_to_name.find(code) == categoryId_to_name.end())
								MTHROW_AND_ERR("CategoryDependencyGenerator::post_learn_from_samples - code %d wasn't found in dict\n", code);
							const vector<string> &names = categoryId_to_name.at(code);
							int nm_idx = 0;
							while (!pass_regex_filter && nm_idx < names.size())
							{
								pass_regex_filter = regex_match(names[nm_idx], reg_pat);
								++nm_idx;
							}
						}
						if (pass_regex_filter) {
							//process request:
#pragma omp critical 
							{
								vector<vector<vector<int>>> &code_stats = categoryVal_to_stats[code]; //gender,age, 4 counts per state
								if (code_stats.empty()) {
									code_stats.resize(2);
									for (size_t kk = 0; kk < 2; ++kk)
										code_stats[kk].resize(age_bin_cnt);
								}
								if (code_stats[gend_idx][age_idx].empty())
									code_stats[gend_idx][age_idx].resize(4);

								vector<vector<bool>> &p_code_lbl_age = code_pid_label_age_bin[code][samples.idSamples[i].samples[j].id];
								if (p_code_lbl_age.empty()) {
									p_code_lbl_age.resize(2);
									for (size_t kk = 0; kk < 2; ++kk)
										p_code_lbl_age[kk].resize(age_bin_cnt);
								}
								if (!p_code_lbl_age[gend_idx][age_idx]) {
									int stat_idx = 2 + outcome_idx;
									p_code_lbl_age[gend_idx][age_idx] = true;
									++code_stats[gend_idx][age_idx][stat_idx];
								}
							}
						}
					}
				}
		}

#pragma omp atomic
		++progress;
		double duration = (unsigned long long)(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
			- tm_prog).count()) / 1000000.0;
		if (progress % 100 == 0 && duration > 15) {
#pragma omp critical
			tm_prog = chrono::high_resolution_clock::now();
			double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
				- tm.t[0]).count()) / 1000000.0;
			double estimate_time = int(double(samples.idSamples.size() - progress) / double(progress) * double(time_elapsed));
			char buffer[1000];
			snprintf(buffer, sizeof(buffer), "Processed %d out of %d(%2.2f%%) time elapsed: %2.1f Minutes, "
				"estimate time to finish %2.1f Minutes",
				progress, (int)samples.idSamples.size(), 100.0*(progress / float(samples.idSamples.size())), time_elapsed / 60,
				estimate_time / 60.0);
			MLOG("%s\n", buffer);
		}
	}
	tm.take_curr_time();
	MLOG("Took %2.2f seconds to complete\n", tm.diff_sec());
	omp_set_nested(nested_state);

	//complete stats in rows:
	for (auto it = categoryVal_to_stats.begin(); it != categoryVal_to_stats.end(); ++it)
		for (size_t i = 0; i < 2; ++i)
			for (size_t j = 0; j < age_bin_cnt; ++j)
				if (!it->second[i][j].empty()) {
					it->second[i][j][0] = total_stats[i][j][0] - it->second[i][j][3];
					it->second[i][j][1] = total_stats[i][j][1] - it->second[i][j][4];
				}

	//filter and take top:
	vector<int> code_list, indexes;
	vector<double> codeCnts, posCnts, lift, scores, pvalues, pos_ratio;
	get_stats(categoryVal_to_stats, code_list, indexes, codeCnts, posCnts, lift, scores, pvalues, pos_ratio);

	apply_filter(indexes, codeCnts, min_code_cnt, INT_MAX);
	apply_filter(indexes, pvalues, 0, fdr);
	vector<int> top_idx(indexes), bottom_idx(indexes);
	apply_filter(top_idx, lift, lift_above, INT_MAX);
	apply_filter(bottom_idx, lift, -1, lift_below);
	indexes.swap(top_idx);
	indexes.insert(indexes.end(), bottom_idx.begin(), bottom_idx.end());
	//join both results from up and down filters on the lift:

	if (indexes.size() > take_top)
		indexes.resize(take_top);
	//Save all codes in indexes:
	top_codes.resize(indexes.size());
	luts.resize(indexes.size());
	int section_id = rep.dict.section_id(signalName);
	for (size_t i = 0; i < indexes.size(); ++i) {
		top_codes[i] = categoryId_to_name.at(code_list[indexes[i]]).front();
		vector<string> s_names = { top_codes[i] };
		rep.dict.prep_sets_lookup_table(section_id, s_names, luts[i]);
	}
	if (verbose) {
		MLOG("CategoryDependencyGenerator on %s - created %d features\n", signalName.c_str(), top_codes.size());
		for (size_t i = 0; i < indexes.size(); ++i)
		{
			MLOG("#NUM %zu:\t%s", i + 1, categoryId_to_name.at(code_list[indexes[i]]).front().c_str());
			for (size_t j = 1; j < 4 && j < categoryId_to_name.at(code_list[indexes[i]]).size(); ++j)
				MLOG("|%s", categoryId_to_name.at(code_list[indexes[i]])[j].c_str());
			MLOG("\tTOT_CNT:%d\tP_VAL=%f\tLift=%1.3f\n",
				(int)codeCnts[indexes[i]], pvalues[indexes[i]], lift[indexes[i]]);
		}
	}
}

void CategoryDependencyGenerator::set_names() {
	names.resize(top_codes.size());
	for (size_t i = 0; i < top_codes.size(); ++i) {
		char buff[5000];
		snprintf(buff, sizeof(buff), "%s.category_dep_set_%s.win_%d_%d",
			signalName.c_str(), top_codes[i].c_str(), win_from, win_to);
		names[i] = string(buff);
	}
}

int CategoryDependencyGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {
	int time_unit_sig = rec.my_base_rep->sigs.Sid2Info[signalId].time_unit;

	MedSample *p_samples = &(features.samples[index]);
	for (size_t k = 0; k < top_codes.size(); ++k)
	{
		//init lut for signal:
		vector<char> &lut = luts[k];
		float *p_feat = _p_data[k] + index;

		for (int i = 0; i < num; i++) {
			rec.uget(signalId, i);
			int min_time, max_time;
			int time = med_time_converter.convert_times(features.time_unit, time_unit_win, p_samples[i].time);
			get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

			bool val = 0;
			for (int i = 0; i < rec.usv.len; i++) {
				int itime = rec.usv.Time(i, time_channel);
				if (itime > max_time) break;
				if (itime >= min_time && lut[(int)rec.usv.Val(i, val_channel)]) {
					val = 1;
					break;
				}
			}

			p_feat[i] = val;
		}
	}

	return 0;
}