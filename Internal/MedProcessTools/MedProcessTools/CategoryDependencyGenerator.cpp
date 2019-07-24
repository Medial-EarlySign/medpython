#include "FeatureGenerator.h"
#include <cmath>
#include <regex>
#include <MedUtils/MedUtils/MedRegistry.h>
#include <omp.h>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void CategoryDependencyGenerator::init_defaults() {
	generator_type = FTR_GEN_CATEGORY_DEPEND;
	signalName = "";
	signalId = -1;
	time_channel = 0;
	val_channel = 0;
	win_from = 0;
	win_to = 360000;
	time_unit_win = MedTime::Days;
	regex_filter = "";
	remove_regex_filter = "";
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
	max_depth = 100;
	max_parents = 5000;
	filter_child_count_ratio = (float)0.05;
	filter_child_pval_diff = (float)1e-10;
	filter_child_lift_ratio = (float)0.05;
	filter_child_removed_ratio = 1;
	filter_hierarchy = true;
	verbose = false;
	use_fixed_lift = false;
	verbose_full = false;

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
		else if (it->first == "remove_regex_filter")
			remove_regex_filter = it->second;
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
		else if (it->first == "filter_hierarchy")
			filter_hierarchy = med_stoi(it->second) > 0;
		else if (it->first == "lift_below")
			lift_below = med_stof(it->second);
		else if (it->first == "lift_above")
			lift_above = med_stof(it->second);
		else if (it->first == "filter_child_count_ratio")
			filter_child_count_ratio = med_stof(it->second);
		else if (it->first == "filter_child_lift_ratio")
			filter_child_lift_ratio = med_stof(it->second);
		else if (it->first == "filter_child_pval_diff")
			filter_child_pval_diff = med_stof(it->second);
		else if (it->first == "filter_child_removed_ratio")
			filter_child_removed_ratio = med_stof(it->second);
		else if (it->first == "chi_square_at_least")
			chi_square_at_least = med_stof(it->second);
		else if (it->first == "minimal_chi_cnt")
			minimal_chi_cnt = med_stoi(it->second);
		else if (it->first == "use_fixed_lift")
			use_fixed_lift = med_stoi(it->second) > 0;
		else if (it->first == "verbose")
			verbose = med_stoi(it->second) > 0;
		else if (it->first == "verbose_full")
			verbose_full = med_stoi(it->second) > 0;
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

void CategoryDependencyGenerator::set_signal_ids(MedSignals& sigs) {
	signalId = sigs.sid(signalName);
	byear_sid = sigs.sid("BYEAR");
	gender_sid = sigs.sid("GENDER");
}

void CategoryDependencyGenerator::init_tables(MedDictionarySections& dict) {
	int section_id = dict.section_id(signalName);
	categoryId_to_name = dict.dict(section_id)->Id2Names;
	_member2Sets = dict.dict(section_id)->Member2Sets;
	_set2Members = dict.dict(section_id)->Set2Members;

	luts.resize(top_codes.size());
	for (size_t i = 0; i < top_codes.size(); ++i) {
		vector<string> s_names = { top_codes[i] };
		dict.prep_sets_lookup_table(section_id, s_names, luts[i]);
	}
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

void CategoryDependencyGenerator::get_parents(int codeGroup, vector<int> &parents, const regex &reg_pat, const regex &remove_reg_pat) {

	bool cached = false;
#pragma omp critical
	if (_member2Sets_flat_cache.find(codeGroup) != _member2Sets_flat_cache.end()) {
		parents = _member2Sets_flat_cache.at(codeGroup);
		cached = true;
	}
	if (cached)
		return;

	vector<int> last_parents = { codeGroup };
	if (last_parents.front() < 0)
		return; //no parents
	parents = {};

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

	if (!regex_filter.empty() || !remove_regex_filter.empty()) {
		vector<int> filtered_p;
		filtered_p.reserve(parents.size());
		for (int code : parents)
		{
			if (categoryId_to_name.find(code) == categoryId_to_name.end())
				MTHROW_AND_ERR("CategoryDependencyGenerator::post_learn_from_samples - code %d wasn't found in dict\n", code);
			const vector<string> &names = categoryId_to_name.at(code);
			int nm_idx = 0;
			bool pass_regex_filter = false;  
			bool pass_remove_regex_filter = false;
			while (!(pass_regex_filter && pass_remove_regex_filter) && nm_idx < names.size())
			{
				if (!regex_filter.empty())
					pass_regex_filter = regex_match(names[nm_idx], reg_pat);
				else
					pass_regex_filter = true;
				if (!remove_regex_filter.empty())
					pass_remove_regex_filter = regex_match(names[nm_idx], remove_reg_pat);
				++nm_idx;
			}
			if (pass_regex_filter && !pass_remove_regex_filter)
				filtered_p.push_back(code);
		}
		parents.swap(filtered_p);
	}

#pragma omp critical
	_member2Sets_flat_cache[codeGroup] = parents;
}

void CategoryDependencyGenerator::get_stats(const unordered_map<int, vector<vector<vector<int>>>> &categoryVal_to_stats,
	vector<int> &all_signal_values, vector<int> &signal_indexes, vector<double> &valCnts,
	vector<double> &posCnts, vector<double> &lift, vector<double> &scores, vector<double> &p_values, vector<double> &pos_ratio,
	vector<int> &dof, const vector<vector<double>> &prior_per_bin) const {

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
	dof.resize(all_signal_values.size());

	for (int index : signal_indexes)
	{
		int signalVal = all_signal_values[index];
		//check chi-square for this value:
		double totCnt = 0, lift_summed = 0;
		const vector<vector<vector<int>>> &code_stats = categoryVal_to_stats.at(signalVal);
		for (size_t i = 0; i < code_stats.size(); ++i) //iterate genders
			for (size_t j = 0; j < code_stats[i].size(); ++j) //iterate age
				if (!code_stats[i][j].empty()) {
					totCnt += code_stats[i][j][2] + code_stats[i][j][3];
					posCnts[index] += code_stats[i][j][1 + 2];
					if (code_stats[i][j][2] + code_stats[i][j][3] > 0 && prior_per_bin[i][j] > 0)
						lift_summed += (code_stats[i][j][1 + 2] / prior_per_bin[i][j]);
				}
		if (totCnt == 0)
			continue;
		valCnts[index] = totCnt; //for signal apeareance
		lift[index] = lift_summed / totCnt;

		pos_ratio[index] = posCnts[index] / totCnt;

		double regScore = 0;
		int dof_val = -1;
		for (size_t i = 0; i < code_stats.size(); ++i) {//iterate genders
			map<float, vector<int>> all_grps;
			for (size_t j = 0; j < code_stats[i].size(); ++j)
				if (!code_stats[i][j].empty())
					all_grps[float(j)] = code_stats[i][j];
			if (stat_metric == category_stat_test::mcnemar)
				regScore += medial::contingency_tables::calc_mcnemar_square_dist(all_grps);
			else if (stat_metric == category_stat_test::chi_square)
				regScore += medial::contingency_tables::calc_chi_square_dist(all_grps, 0, chi_square_at_least, minimal_chi_cnt);
			dof_val += _count_legal_rows(code_stats[i], stat_metric == category_stat_test::chi_square ? minimal_chi_cnt : 0);
		}
		scores[index] = (float)regScore;
		dof[index] = dof_val;

		double pv = medial::contingency_tables::chisqr(dof_val, regScore);
		p_values[index] = pv;
	}
}

int CategoryDependencyGenerator::_learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) {

	if (signalId == -1 || byear_sid == -1 || gender_sid == -1)
		MTHROW_AND_ERR("Uninitialized signalId,byear_sid or gender_sid - or not loaded\n");
	names.clear();

	// Required signals
	vector<int> all_req_signal_ids_v;
	vector<unordered_set<int> > current_required_signal_ids(processors.size());
	vector<FeatureGenerator *> generators = { this };
	unordered_set<int> extra_req_signal_ids;
	handle_required_signals(processors, generators, extra_req_signal_ids, all_req_signal_ids_v, current_required_signal_ids);
	regex reg_pat;
	regex remove_reg_pat;
	if (!regex_filter.empty())
		reg_pat = regex(regex_filter);
	if (!remove_regex_filter.empty())
		remove_reg_pat = regex(remove_regex_filter);

	// Preparations
	unordered_map<int, vector<vector<vector<int>>>> categoryVal_to_stats; //stats is gender,age, 4 ints counts:
	int age_bin_cnt = int(ceil((max_age - min_age + 1) / float(age_bin)));

	vector<vector<vector<int>>> total_stats(2); //gender, age, label
	for (size_t i = 0; i < 2; ++i) {
		total_stats[i].resize(age_bin_cnt);
		for (size_t j = 0; j < age_bin_cnt; ++j)
			total_stats[i][j].resize(2);
	}

	MedProgress progress("CategoryDependencyGenerator:" + signalName, (int)samples.idSamples.size(), 15,100);
	//unordered_map<int, unordered_map<int, vector<vector<bool>>>> code_pid_label_age_bin;// stores for each code => pid if saw label,age_bin
	//bool nested_state = omp_get_nested();
	//omp_set_nested(true);

	int N_tot_threads = omp_get_max_threads();
	vector<PidDynamicRec> idRec(N_tot_threads);

	// Collect data
#pragma omp parallel for schedule(dynamic,128)
	for (int i = 0; i < samples.idSamples.size(); ++i)
	{
		int nSamples = (int)samples.idSamples[i].samples.size();
		unordered_map<int, vector<vector<vector<bool>>>> pid_categoryVal_to_stats;

		UniversalSigVec usv;
		int n_th = omp_get_thread_num();

		int pid = samples.idSamples[i].id;
		int gend_idx = medial::repository::get_value(rep, pid, gender_sid) - 1;
		int byear = medial::repository::get_value(rep, pid, byear_sid);

		idRec[n_th].init_from_rep(std::addressof(rep), pid, all_req_signal_ids_v, nSamples);

		// Apply Processors
		for (unsigned int j = 0; j < processors.size(); j++)
			processors[j]->conditional_apply_without_attributes(idRec[n_th], samples.idSamples[i], current_required_signal_ids[j]);

		vector<vector<bool>> p_lbl_age;
		p_lbl_age.resize(2);
		for (size_t k = 0; k < 2; ++k)
			p_lbl_age[k].resize(age_bin_cnt);

		for (size_t j = 0; j < nSamples; ++j)
		{
			// Collect Data
			idRec[n_th].uget(signalId, (int)j, usv);

			int age = med_time_converter.convert_times(samples.time_unit, MedTime::Date, samples.idSamples[i].samples[j].time) / 10000 - byear;
			if (age > max_age || age < min_age)
				continue;

			int outcome_idx = samples.idSamples[i].samples[j].outcome > 0;
			int age_idx = (age - min_age) / age_bin;
			if (!p_lbl_age[outcome_idx][age_idx]) {
				p_lbl_age[outcome_idx][age_idx] = true;
#pragma omp critical(sigdep_part1)
				++total_stats[gend_idx][age_idx][outcome_idx];
			}

			//allocate stats per category value:

			int start_time_win = med_time_converter.convert_times(time_unit_win, samples.time_unit, med_time_converter.convert_times(samples.time_unit, time_unit_win, samples.idSamples[i].samples[j].time) - win_to);
			int end_time_win = med_time_converter.convert_times(time_unit_win, samples.time_unit, med_time_converter.convert_times(samples.time_unit, time_unit_win, samples.idSamples[i].samples[j].time) - win_from);

			for (int k = 0; k < usv.len; ++k)
				if (usv.Time(k, time_channel) >= start_time_win && usv.Time(k, time_channel) <= end_time_win) { //get values in time window:
					//get filter regex codes:
					int base_code = (int)usv.Val(k, val_channel);

					//process request:
					vector<vector<vector<bool>>> &code_stats = pid_categoryVal_to_stats[base_code]; //gender,age, 4 counts per state
					if (code_stats.empty()) {
						code_stats.resize(2);
						for (size_t kk = 0; kk < 2; ++kk)
							code_stats[kk].resize(age_bin_cnt);
					}
					if (code_stats[gend_idx][age_idx].empty())
						code_stats[gend_idx][age_idx].resize(2);

					if (!code_stats[gend_idx][age_idx][outcome_idx]) {
						/*MLOG_D("DEBUG: code=%d,gend_idx=%d,age_idx=%d,outcome_idx=%d pid=%d, timepoint_idx=%d, time=%d\n",
							code,gend_idx, age_idx, outcome_idx, pid, k, usv.Time(k, time_channel));*/
						code_stats[gend_idx][age_idx][outcome_idx] = true;
					}

				}
		}

		//process hirarchy:
		vector<int> update_ls;
		update_ls.reserve(pid_categoryVal_to_stats.size());
		for (auto it = pid_categoryVal_to_stats.begin(); it != pid_categoryVal_to_stats.end(); ++it)
			update_ls.push_back(it->first);
		for (int base_code : update_ls) {
			vector<int> all_parents;
			get_parents(base_code, all_parents, reg_pat, remove_reg_pat);

			const vector<vector<vector<bool>>> &base_code_stats = pid_categoryVal_to_stats.at(base_code);
			for (int code : all_parents)
			{
				//process request for code - aggregate stats from [2+0] [2+1] to code:
				vector<vector<vector<bool>>> &code_stats = pid_categoryVal_to_stats[code]; //gender,age, 4 counts per state
				if (code_stats.empty()) {
					code_stats.resize(2);
					for (size_t kk = 0; kk < 2; ++kk)
						code_stats[kk].resize(age_bin_cnt);
				}
				for (size_t ii = 0; ii < base_code_stats.size(); ++ii) {
					for (size_t jj = 0; jj < base_code_stats[ii].size(); ++jj)
						if (!base_code_stats[ii][jj].empty()) {
							if (code_stats[ii][jj].empty())
								code_stats[ii][jj].resize(2);
							for (size_t k = 0; k < code_stats[ii][jj].size(); ++k)
								code_stats[ii][jj][k] = code_stats[ii][jj][k] || base_code_stats[ii][jj][k];
						}
				}
			}

		}
		//update original:
#pragma omp critical(sigdep_part2)
		{
			for (auto it = pid_categoryVal_to_stats.begin(); it != pid_categoryVal_to_stats.end(); ++it)
			{
				int code = it->first;
				const vector<vector<vector<bool>>> &pid_code_stats = it->second; //gender,age, 4 counts per state
				vector<vector<vector<int>>> &code_stats = categoryVal_to_stats[code]; //gender,age, 4 counts per state
				if (code_stats.empty()) {
					code_stats.resize(2);
					for (size_t kk = 0; kk < 2; ++kk)
						code_stats[kk].resize(age_bin_cnt);
				}
				for (size_t ii = 0; ii < code_stats.size(); ++ii)
					for (size_t jj = 0; jj < code_stats[ii].size(); ++jj) {
						if (!pid_code_stats[ii][jj].empty()) {
							if (code_stats[ii][jj].empty())
								code_stats[ii][jj].resize(4);
							for (size_t k = 0; k < pid_code_stats[ii][jj].size(); ++k)
								code_stats[ii][jj][2 + k] += int(pid_code_stats[ii][jj][k]);
						}
					}
			}

			// Some timing printing
			//#pragma omp atomic
			progress.update();
			
		}
	}
	//omp_set_nested(nested_state);

	//complete stats in rows:
	for (auto it = categoryVal_to_stats.begin(); it != categoryVal_to_stats.end(); ++it)
		for (size_t i = 0; i < 2; ++i)
			for (size_t j = 0; j < age_bin_cnt; ++j)
				if (!it->second[i][j].empty()) {
					it->second[i][j][0] = total_stats[i][j][0] - it->second[i][j][2];
					it->second[i][j][1] = total_stats[i][j][1] - it->second[i][j][3];
					if (it->second[i][j][0] < 0 || it->second[i][j][1] < 0)
						MTHROW_AND_ERR("Bug in calc - negative count in stat bin\n");
				}

	//filter regex if given:
	if (!regex_filter.empty() || !remove_regex_filter.empty())
		for (auto it = categoryVal_to_stats.begin(); it != categoryVal_to_stats.end();) {
			int base_code = it->first;
			bool found_match = false;
			bool found_remove_match = false;
			const vector<string> &names = categoryId_to_name.at(base_code);
			int pos_i = 0;
			while (pos_i < names.size() && !found_match) {
				if (!regex_filter.empty())
					found_match = regex_match(names[pos_i], reg_pat);
				else
					found_match = true;
				found_remove_match = regex_match(names[pos_i], remove_reg_pat);
				++pos_i;
			}

			if (found_match && !found_remove_match)
				++it;
			else
				it = categoryVal_to_stats.erase(it);
		}

	if (verbose_full) {
		for (auto it = categoryVal_to_stats.begin(); it != categoryVal_to_stats.end(); ++it) {
			MLOG("code=%d(%s):\n", it->first, categoryId_to_name.at(it->first).back().c_str());
			for (size_t ii = 0; ii < 2; ++ii)
				for (size_t jj = 0; jj < age_bin_cnt; ++jj)
					if (!it->second[ii][jj].empty())
						MLOG("Gender_idx=%zu,Age_idx=%zu\tctrls:%d\tcases:%d\ttot_ctrls:%d\ttot_cases:%d\n",
							ii, jj, it->second[ii][jj][2], it->second[ii][jj][3],
							total_stats[ii][jj][0], total_stats[ii][jj][1]);
		}
	}


	vector<vector<double>> prior_per_bin(2);
	for (size_t i = 0; i < prior_per_bin.size(); ++i)
	{
		prior_per_bin[i].resize(age_bin_cnt);
		for (size_t j = 0; j < prior_per_bin[i].size(); ++j)
			if (total_stats[i][j][1] + total_stats[i][j][0] > 0)
				prior_per_bin[i][j] = double(total_stats[i][j][1]) / (total_stats[i][j][1] + total_stats[i][j][0]);
	}
	//filter and take top:
	vector<int> code_list, indexes, dof;
	vector<double> codeCnts, posCnts, lift, scores, pvalues, pos_ratio;
	get_stats(categoryVal_to_stats, code_list, indexes, codeCnts, posCnts, lift, scores, pvalues, pos_ratio, dof, prior_per_bin);
	if (verbose_full) {
		for (unsigned int i = 0; i < code_list.size(); i++)
			MLOG("Value=%d : Cnts = %f PosCnts = %f Lift = %f Score = %f P-value = %f\n", code_list[i], codeCnts[i], posCnts[i], lift[i], scores[i], pvalues[i]);
	}
	unordered_map<int, double> code_cnts;
	for (size_t i = 0; i < code_list.size(); ++i)
		code_cnts[code_list[i]] = codeCnts[i];

	int before_cnt = (int)indexes.size();
	apply_filter(indexes, codeCnts, min_code_cnt, INT_MAX);
	if (verbose)
		MLOG("CategoryDependencyGenerator on %s - count_filter left %zu(out of %d)\n",
			signalName.c_str(), indexes.size(), before_cnt);
	before_cnt = (int)indexes.size();

	vector<int> top_idx(indexes), bottom_idx(indexes);
	apply_filter(top_idx, lift, lift_above, INT_MAX);
	apply_filter(bottom_idx, lift, -1, lift_below);
	indexes.swap(top_idx);
	indexes.insert(indexes.end(), bottom_idx.begin(), bottom_idx.end());
	if (verbose)
		MLOG("CategoryDependencyGenerator on %s - lift_filter left %zu(out of %d)\n",
			signalName.c_str(), indexes.size(), before_cnt);
	before_cnt = (int)indexes.size();
	//filter hierarchy:
	if (filter_hierarchy) {
		medial::contingency_tables::filterHirarchy(_member2Sets, _set2Members, indexes, code_list, pvalues, codeCnts, lift,
			code_cnts, filter_child_pval_diff, filter_child_lift_ratio, filter_child_count_ratio, filter_child_removed_ratio);
		if (verbose)
			MLOG("CategoryDependencyGenerator on %s - Hirarchy_filter left %zu(out of %d)\n",
				signalName.c_str(), indexes.size(), before_cnt);
		before_cnt = (int)indexes.size();
	}
	//Real FDR filtering:
	vector<double> fixed_lift(lift); //for sorting
	for (int index : indexes)
		if (fixed_lift[index] < 1 && fixed_lift[index] > 0)
			fixed_lift[index] = 1 / fixed_lift[index];
		else if (fixed_lift[index] == 0)
			fixed_lift[index] = 999999;

	medial::contingency_tables::FilterFDR(indexes, scores, pvalues, fixed_lift, fdr);
	if (verbose)
		MLOG("CategoryDependencyGenerator on %s - fdr_filter left %zu(out of %d)\n",
			signalName.c_str(), indexes.size(), before_cnt);

	before_cnt = (int)indexes.size();
	//join both results from up and down filters on the lift:
	//sort before taking top:
	//sort by p_value, lift, score:

	vector<int> indexes_order(indexes.size());
	vector<pair<int, vector<double>>> sort_pars(indexes.size());
	for (size_t i = 0; i < indexes.size(); ++i)
	{
		sort_pars[i].first = (int)indexes[i];
		sort_pars[i].second.resize(3);
		sort_pars[i].second[0] = pvalues[indexes[i]];
		if (use_fixed_lift)
			sort_pars[i].second[1] = -fixed_lift[indexes[i]];
		else
			sort_pars[i].second[1] = -lift[indexes[i]];
		sort_pars[i].second[2] = -scores[indexes[i]];
	}
	sort(sort_pars.begin(), sort_pars.end(), [](pair<int, vector<double>> a, pair<int, vector<double>> b) {
		int pos = 0;
		while (pos < a.second.size() &&
			a.second[pos] == b.second[pos])
			++pos;
		if (pos >= a.second.size())
			return false;
		return b.second[pos] > a.second[pos];
	});
	for (size_t i = 0; i < sort_pars.size(); ++i)
		indexes_order[i] = sort_pars[i].first;
	if (indexes_order.size() > take_top)
		indexes_order.resize(take_top);
	//Save all codes in indexes:
	top_codes.resize(indexes_order.size());
	luts.resize(indexes_order.size());
	int section_id = rep.dict.section_id(signalName);
	for (size_t i = 0; i < indexes_order.size(); ++i) {
		top_codes[i] = categoryId_to_name.at(code_list[indexes_order[i]]).front();
		vector<string> s_names = { top_codes[i] };
		rep.dict.prep_sets_lookup_table(section_id, s_names, luts[i]);
	}
	if (verbose) {
		MLOG("CategoryDependencyGenerator on %s - created %d(out of %zu passed filters) features\n",
			signalName.c_str(), top_codes.size(), sort_pars.size());
		for (size_t i = 0; i < indexes_order.size(); ++i)
		{
			MLOG("#NUM %zu:\t%s", i + 1, categoryId_to_name.at(code_list[indexes_order[i]]).front().c_str());
			for (size_t j = 1; j < 4 && j < categoryId_to_name.at(code_list[indexes_order[i]]).size(); ++j)
				MLOG("|%s", categoryId_to_name.at(code_list[indexes_order[i]])[j].c_str());
			MLOG("\tTOT_CNT:%d\tP_VAL=%.12g\tScore=%.3f\tLift=%1.3f\tDOF=%d\n",
				(int)codeCnts[indexes_order[i]], pvalues[indexes_order[i]], scores[indexes_order[i]], lift[indexes_order[i]],
				dof[indexes_order[i]]);
		}
	}

	return 0;
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

int CategoryDependencyGenerator::nfeatures() {
	if (!names.empty())
		return (int)names.size();
	if (take_top > 0)
		return take_top;
	return 0;
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

//.......................................................................................
// Filter generated features according to a set. return number of valid features (does not affect single-feature genertors, just returns 1/0 if feature name in set)
int CategoryDependencyGenerator::filter_features(unordered_set<string>& validFeatures) {

	vector<int> selected;
	for (int i = 0; i < names.size(); i++) {
		if (validFeatures.find(names[i]) != validFeatures.end())
			selected.push_back(i);
	}

	for (int i = 0; i < selected.size(); i++)
		top_codes[i] = top_codes[selected[i]];
	top_codes.resize(selected.size());
	set_names();

	return ((int)names.size());
}