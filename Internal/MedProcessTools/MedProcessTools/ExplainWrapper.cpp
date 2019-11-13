#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <cmath>
#include <random>
#include <omp.h>
#include "ExplainWrapper.h"
#include <MedAlgo/MedAlgo/MedXGB.h>
#include <MedStat/MedStat/MedStat.h>
#include "medial_utilities/medial_utilities/globalRNG.h"

#define LOCAL_SECTION LOG_MED_MODEL
#define LOCAL_LEVEL LOG_DEF_LEVEL

ExplainFilters::ExplainFilters() {
	max_count = 0;
	sum_ratio = 1;
	sort_mode = 0;
}

int ExplainFilters::init(map<string, string> &map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "max_count")
			max_count = med_stoi(it->second);
		else if (it->first == "sum_ratio")
			sum_ratio = med_stof(it->second);
		else if (it->first == "sort_mode")
			sort_mode = med_stoi(it->second);
		else
			MTHROW_AND_ERR("Error in ExplainFilters::init - Unknown param \"%s\"\n", it->first.c_str());
	}

	if (sum_ratio < 0 || sum_ratio > 1)
		MTHROW_AND_ERR("Error in ExplainFilters::init - sum_ratio should be in [0,1]\n");
	return 0;
}

void ExplainFilters::filter(map<string, float> &explain_list) const {
	vector<pair<string, float>> sorted;
	sorted.reserve(explain_list.size());
	for (const auto &it : explain_list)
		if (sort_mode == 0 || (sort_mode > 0 && it.second > 0) || (sort_mode < 0 && it.second < 0))
			sorted.push_back(pair<string, float>(it.first, it.second));
	if (sort_mode == 0)
		sort(sorted.begin(), sorted.end(), [](const pair<string, float>&pr1, const pair<string, float>&pr2)
	{ return abs(pr1.second) > abs(pr2.second); });
	else if (sort_mode > 0)
		sort(sorted.begin(), sorted.end(), [](const pair<string, float>&pr1, const pair<string, float>&pr2)
	{ return pr1.second > pr2.second; });
	else if (sort_mode < 0)
		sort(sorted.begin(), sorted.end(), [](const pair<string, float>&pr1, const pair<string, float>&pr2)
	{ return pr1.second < pr2.second; });

	//filter ratio:
	if (sum_ratio < 1) {
		float tot = 0, curr_sum = 0;
		for (const auto &it : sorted)
			tot += abs(it.second);
		int stop_at = 0;
		while (stop_at < sorted.size() && curr_sum / tot > sum_ratio) {
			curr_sum += abs(sorted[stop_at].second);
			++stop_at;
		}
		sorted.resize(stop_at);
	}
	//filter max_count:
	if (max_count > 0 && sorted.size() > max_count)
		sorted.resize(max_count);
	//commit selections:
	map<string, float> filterd;
	for (const auto &it : sorted)
		filterd[it.first] = it.second;
	explain_list = move(filterd);
}

ExplainProcessings::ExplainProcessings() {
	group_by_sum = false;
	learn_cov_matrix = false;
}

int ExplainProcessings::init(map<string, string> &map) {
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		if (it->first == "group_by_sum")
			group_by_sum = med_stoi(it->second) > 0;
		else if (it->first == "learn_cov_matrix")
			learn_cov_matrix = med_stoi(it->second) > 0;
		else if (it->first == "cov_features")
			cov_features.read_from_csv_file(it->second, 1);
		else if (it->first == "grouping")
			grouping = it->second;
		else if (it->first == "zero_missing")
			zero_missing = stoi(it->second);
		else if (it->first == "normalize_vals")
			normalize_vals = stoi(it->second);
		else
			MTHROW_AND_ERR("Error in ExplainProcessings::init - Unknown param \"%s\"\n", it->first.c_str());
	}

	return 0;
}


void ExplainProcessings::post_deserialization()
{
	abs_cov_features.clear();
	if (cov_features.m.size() > 0) {
		abs_cov_features.resize(cov_features.nrows, cov_features.ncols);
		for (int i = 0; i < cov_features.m.size(); i++)
			abs_cov_features.m[i] = abs(cov_features.m[i]);
	}

	groupName2Inds.clear();
	for (int i = 0; i < groupNames.size(); i++)
		groupName2Inds[groupNames[i]] = group2Inds[i];
}

void ExplainProcessings::learn(const MedFeatures &train_mat) {
	if (learn_cov_matrix) {
		//int feat_cnt = (int)train_mat.data.size();
		//cov_features.resize(feat_cnt, feat_cnt);
		MLOG("Calc Covariance mat\n");
		MedMat<float> x_mat;
		train_mat.get_as_matrix(x_mat);
		x_mat.normalize();
		//0 - no transpose, 1 - A_Transpose * B, 2 - A * B_Transpose, 3 - both transpose
		fast_multiply_medmat_transpose(x_mat, x_mat, cov_features, 1, 1.0 / x_mat.nrows);

		abs_cov_features = cov_features;
		for (auto &x : abs_cov_features.m) x = abs(x);

		//// debug
		//vector<string> f_names;
		//train_mat.get_feature_names(f_names);
		//for (int i = 0; i < cov_features.nrows; i++)
		//	for (int j = 0; j < cov_features.ncols; j++)
		//		MLOG("COV_DEBUG: (%d) %s (%d) %s :: %f\n", i, f_names[i].c_str(), j, f_names[j].c_str(), cov_features(i, j));
	}
}

float ExplainProcessings::get_group_normalized_contrib(const vector<int> &group_inds, vector<float> &contribs, float total_normalization_factor) const
{
	float group_normalization_factor = (float)1e-8;

	for (auto i : group_inds) group_normalization_factor += abs(contribs[i]);

	vector<int> group_mask(contribs.size(), 0);
	for (auto i : group_inds) group_mask[i] = 1;

	vector<float> alphas(contribs.size());

	for (int i = 0; i < group_mask.size(); i++) {
		if (group_mask[i])
			alphas[i] = 1.0f;
		else {
			alphas[i] = 0.0f;
			for (int j : group_inds)
				if (abs_cov_features(j, i) > alphas[i])
					alphas[i] = abs_cov_features(j, i);
			//alphas[i] += abs_cov_features(j, i) * abs(contribs[j]);
		//alphas[i] /= group_normalization_factor;
		}
	}

	float group_contrib = 0.0f;
	for (int i = 0; i < contribs.size(); i++)
		group_contrib += alphas[i] * contribs[i];

	group_contrib /= total_normalization_factor;
	return group_contrib;
}

void ExplainProcessings::process(map<string, float> &explain_list) const {

	if (cov_features.m.empty() && !group_by_sum && normalize_vals <= 0)
		return;

	unordered_set<string> skip_bias_names = { "b0", "Prior_Score" };
	for (auto &s : skip_bias_names) explain_list.erase(s);
	MedMat<float> orig_explain((int)explain_list.size(), 1);
	int k = 0;
	for (auto &e : explain_list) orig_explain(k++, 0) = e.second;


	float normalization_factor = 1.0;
	if (normalize_vals > 0) {

		normalization_factor = (float)1e-8; // starting with a small epsilon so that we never divide by 0 later
		for (auto &e : explain_list) normalization_factor += abs(e.second);
		//MLOG("====> DEBUG normalization_factor %f\n", normalization_factor);
	}

	//first do covarinace if has:
	if (!cov_features.m.empty()) {
		if (cov_features.ncols != explain_list.size() && cov_features.ncols != (int)explain_list.size() - 1)
			MTHROW_AND_ERR("Error in ExplainProcessings::process - processing covarince agg. wrong sizes. cov_features.ncols=%d, "
				"explain_list.size()=%zu\n", cov_features.ncols, explain_list.size());



		if (group_by_sum) {
			map<string, float> group_explain;
			for (int i = 0; i < group2Inds.size(); i++) {
				group_explain[groupNames[i]] = get_group_normalized_contrib(group2Inds[i], orig_explain.m, normalization_factor);
			}
			explain_list = move(group_explain);
		}
		else {
			MedMat<float> fixed_with_cov(cov_features.ncols, 1);

			fast_multiply_medmat(abs_cov_features, orig_explain, fixed_with_cov, (float)1.0 / normalization_factor);
			int k = 0;
			for (auto &e : explain_list) explain_list[e.first] = fixed_with_cov(k++, 0);
		}

		return; // ! -> since we treat group_by_sum differently in this case

	}

	//sum features in groups
	if (group_by_sum) {
		if (group2Inds.empty())
			MTHROW_AND_ERR("Error in ExplainProcessings::process - asked for group_by_sum but haven't provide groups in grouping\n");
		map<string, float> group_explain;
		for (size_t i = 0; i < group2Inds.size(); ++i)
		{
			const string &grp_name = groupNames[i];
			float contrib = 0.0f;
			for (int ind : group2Inds[i])
				contrib += orig_explain.m[ind];
			group_explain[grp_name] = contrib;

		}

		normalization_factor = 0;
		if (normalize_vals > 0) {
			for (const auto &e : group_explain) normalization_factor += abs(e.second);

			if (normalization_factor > 0)
				for (auto &e : group_explain) e.second /= normalization_factor;
		}

		explain_list = move(group_explain);
	}
	else {
		normalization_factor = 0;
		if (normalize_vals > 0) {
			for (const auto &e : explain_list) normalization_factor += abs(e.second);

			if (normalization_factor > 0)
				for (auto &e : explain_list) e.second /= normalization_factor;
		}
	}
}


void ExplainProcessings::process(map<string, float> &explain_list, unsigned char *missing_value_mask) const
{
	process(explain_list);
	if (zero_missing == 0 || missing_value_mask == NULL) 	return;


	if (!group_by_sum) {
		unordered_set<string> skip_bias_names = { "b0", "Prior_Score" };
		for (auto &s : skip_bias_names) explain_list.erase(s);

		// now zero all missing
		int k = 0;
		for (auto &e : explain_list) {
			//MLOG("feat[%d] : %s : %6.4f : mask = %x\n", k, e.first.c_str(), e.second, missing_value_mask[k]);
			if (missing_value_mask[k++] & MedFeatures::imputed_mask)
				e.second = 0;
		}

	}
	else {

		for (auto &g : groupName2Inds) {

			int is_empty = 1;
			for (auto v : g.second)
				if (!(missing_value_mask[v] & MedFeatures::imputed_mask)) {
					is_empty = 0;
					break;
				}

			if (is_empty)
				explain_list[g.first] = 0;
		}

	}



}

int ModelExplainer::init(map<string, string> &mapper) {
	map<string, string> left_to_parse;
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "processing")
			processing.init_from_string(it->second);
		else if (it->first == "filters")
			filters.init_from_string(it->second);
		else if (it->first == "attr_name")
			attr_name = it->second;
		else if (it->first == "use_split")
			use_split = stoi(it->second);
		else if (it->first == "use_p")
			use_p = stof(it->second);
		else if (it->first == "pp_type") {} //ignore
		else {
			left_to_parse[it->first] = it->second;
		}
	}
	_init(left_to_parse);

	return 0;
}

void ModelExplainer::explain(MedFeatures &matrix) const {
	vector<map<string, float>> explain_reasons; //for each sample, reasons and scores
	explain(matrix, explain_reasons);

	if (explain_reasons.size() != matrix.samples.size())
		MTHROW_AND_ERR("Error in ModelExplainer::explain - explain returned musmatch number of samples %zu, and requested %zu\n",
			explain_reasons.size(), matrix.samples.size());

	MedMat<unsigned char> masks_mat;
	if (processing.zero_missing)
		matrix.get_masks_as_mat(masks_mat);
	//process:
	for (size_t i = 0; i < explain_reasons.size(); ++i) {
		if (processing.zero_missing)
			processing.process(explain_reasons[i], &masks_mat.m[i*masks_mat.ncols]);
		else
			processing.process(explain_reasons[i]);
	}
	//filter:
	for (size_t i = 0; i < explain_reasons.size(); ++i)
		filters.filter(explain_reasons[i]);

	string group_name = attr_name;
	if (attr_name.empty()) //default name
		group_name = my_class_name();
#pragma omp critical
	{
		for (size_t i = 0; i < explain_reasons.size(); ++i)
			for (auto it = explain_reasons[i].begin(); it != explain_reasons[i].end(); ++it)
				matrix.samples[i].attributes[group_name + "::" + it->first] = it->second;
	}

}

///format TAB delim, 2 tokens: [Feature_name [TAB] group_name]
void read_feature_grouping(const string &file_name, const MedFeatures& data, vector<vector<int>>& group2index,
	vector<string>& group_names) {
	// Features
	vector<string> features;
	data.get_feature_names(features);
	int nftrs = (int)features.size();
	map<string, vector<int>> groups;
	vector<bool> grouped_ftrs(nftrs);

	if (file_name == "BY_SIGNAL") {
		for (int i = 0; i < nftrs; ++i)
		{
			vector<string> tokens;
			boost::split(tokens, features[i], boost::is_any_of("."));
			string word = tokens[0];
			if (tokens.size() > 1 && boost::starts_with(tokens[0], "FTR_"))
				word = tokens[1];

			groups[word].push_back(i);
			grouped_ftrs[i] = true;
		}
	}
	else if (file_name == "BY_SIGNAL_CATEG") {
		for (int i = 0; i < nftrs; ++i)
		{
			vector<string> tokens;
			boost::split(tokens, features[i], boost::is_any_of("."));
			string word = tokens[0];
			int idx = 0;
			if (tokens.size() > 1 && boost::starts_with(tokens[0], "FTR_")) {
				word = tokens[1];
				idx = 1;
			}
			if (idx + 1 < tokens.size()) {
				if (boost::starts_with(tokens[idx + 1], "category_")) {
					boost::replace_all(tokens[idx + 1], "category_set_count_", "");
					boost::replace_all(tokens[idx + 1], "category_set_sum_", "");
					boost::replace_all(tokens[idx + 1], "category_set_first_", "");
					boost::replace_all(tokens[idx + 1], "category_set_first_time_", "");
					boost::replace_all(tokens[idx + 1], "category_dep_set_", "");
					boost::replace_all(tokens[idx + 1], "category_set_", "");
					word += "." + tokens[idx + 1];
				}
			}

			groups[word].push_back(i);
			grouped_ftrs[i] = true;
		}
	}
	else {
		// Read Grouping
		ifstream inf(file_name);
		if (!inf.is_open())
			MTHROW_AND_ERR("Cannot open \'%s\' for reading\n", file_name.c_str());

		string curr_line;
		vector<string> fields;

		while (getline(inf, curr_line)) {
			boost::split(fields, curr_line, boost::is_any_of("\t"));
			if (fields.size() != 2)
				MTHROW_AND_ERR("Cannot parse line \'%s\' from %s\n", curr_line.c_str(), file_name.c_str());
			int feat_pos = find_in_feature_names(features, fields[0]);
			if (grouped_ftrs[feat_pos])
				MTHROW_AND_ERR("Features %s given twice\n", fields[0].c_str());

			grouped_ftrs[feat_pos] = true;
			groups[fields[1]].push_back(feat_pos);
		}
	}
	// Arrange
	for (auto& rec : groups) {
		group_names.push_back(rec.first);
		group2index.push_back(rec.second);
	}

	for (int i = 0; i < nftrs; i++) {
		if (!grouped_ftrs[i]) {
			group_names.push_back(features[i]);
			group2index.push_back({ i });
		}
	}

	MLOG("Grouping: %d features into %d groups\n", nftrs, (int)group_names.size());
}

void ModelExplainer::Learn(const MedFeatures &train_mat) {
	if (original_predictor == NULL)
		MTHROW_AND_ERR("Error ModelExplainer::Learn - please call init_post_processor before learn\n");
	if (!processing.grouping.empty())
		read_feature_grouping(processing.grouping, train_mat, processing.group2Inds, processing.groupNames);
	else {
		int icol = 0;
		for (auto& rec : train_mat.data) {
			processing.group2Inds.push_back({ icol++ });
			processing.groupNames.push_back(rec.first);
		}
	}
	processing.learn(train_mat);
	_learn(train_mat);
}

void ModelExplainer::dprint(const string &pref) const {
	string predictor_nm = "";
	if (original_predictor != NULL)
		predictor_nm = original_predictor->my_class_name();
	string filters_str = "", processing_str = "";
	char buffer[5000];
	snprintf(buffer, sizeof(buffer), "group_by_sum=%d, learn_cov_matrix=%d, normalize_vals=%d, zero_missing=%d, grouping=%s",
		int(processing.group_by_sum), int(processing.learn_cov_matrix), processing.normalize_vals
		, processing.zero_missing, processing.grouping.c_str());
	processing_str = string(buffer);
	snprintf(buffer, sizeof(buffer), "sort_mode=%d, max_count=%d, sum_ratio=%2.3f",
		filters.sort_mode, filters.max_count, filters.sum_ratio);
	filters_str = string(buffer);
	MLOG("%s :: ModelExplainer type %d(%s), original_predictor=%s, attr_name=%s, processing={%s}, filters={%s}\n",
		pref.c_str(), processor_type, my_class_name().c_str(), predictor_nm.c_str(), attr_name.c_str(),
		processing_str.c_str(), filters_str.c_str());
}

bool comp_score_str(const pair<string, float> &pr1, const pair<string, float> &pr2) {
	return abs(pr1.second) > abs(pr2.second); //bigger is better in absolute
}

void ModelExplainer::print_explain(MedSample &smp, int sort_mode) {
	vector<pair<string, float>> ranked;
	for (auto it = smp.attributes.begin(); it != smp.attributes.end(); ++it)
		if (boost::starts_with(it->first, "ModelExplainer::"))
			ranked.push_back(pair<string, float>(it->first, it->second));
	if (sort_mode == 0)
		sort(ranked.begin(), ranked.end(), comp_score_str);
	else if (sort_mode > 0)
		sort(ranked.begin(), ranked.end(), [](const pair<string, float>&pr1, const pair<string, float>&pr2) { return pr1.second > pr2.second; });
	else
		sort(ranked.begin(), ranked.end(), [](const pair<string, float>&pr1, const pair<string, float>&pr2) { return pr1.second < pr2.second; });

	for (size_t i = 0; i < ranked.size(); ++i)
		MLOG("%s = %f\n", ranked[i].first.c_str(), ranked[i].second);
	if (!smp.prediction.empty())
		MLOG("ModelExplainer::Prediction_Raw_Score = %f\n", smp.prediction[0]);
}

int get_max_rec(const vector<QRF_ResNode> &nodes, int idx) {
	const QRF_ResNode &curr = nodes[idx];
	if (curr.is_leaf)
		return 0; //reached leaf
	int max_d = get_max_rec(nodes, curr.left) + 1;
	int max_right = get_max_rec(nodes, curr.right) + 1;
	if (max_right > max_d)
		max_d = max_right;

	return max_d;
}

int get_tree_max_depth(const vector<QRF_ResNode> &nodes) {
	if (nodes.empty())
		return 0;
	int max_d = get_max_rec(nodes, 0) + 1;
	return max_d;
}

//all conversion functions
bool TreeExplainer::try_convert_trees() {
	//convert QRF, BART to generic_tree_model structure:
	if (original_predictor->classifier_type == MODEL_QRF) {
		MLOG("Converting QRF to generic ensemble trees\n");
		int num_outputs = original_predictor->n_preds_per_sample();
		const QRF_Forest &forest = static_cast<MedQRF *>(original_predictor)->qf;
		if (forest.n_categ == 2)
			num_outputs = 1;
		const vector<float> &all_vals = forest.sorted_values;
		if (all_vals.empty() && forest.n_categ != 2) {
			MWARN("Can't convert QRF. please retrain with keep_all_values to be able to convert categorical trees with n_categ > 2\n");
			return false;
		}
		int max_nodes = 0, max_depth = 0;
		for (size_t i = 0; i < forest.qtrees.size(); ++i)
		{
			if (max_nodes < forest.qtrees[i].qnodes.size())
				max_nodes = (int)forest.qtrees[i].qnodes.size();
			int mm = get_tree_max_depth(forest.qtrees[i].qnodes);
			if (mm > max_depth)
				max_depth = mm;
		}

		++max_nodes; //this is tree seperator index for model
		generic_tree_model.allocate((int)forest.qtrees.size(), max_nodes, num_outputs);
		int pos_in_model = 0;
		generic_tree_model.max_depth = max_depth;
		generic_tree_model.tree_limit = (int)forest.qtrees.size();
		generic_tree_model.max_nodes = max_nodes;
		generic_tree_model.num_outputs = num_outputs;
		generic_tree_model.base_offset = 0; //no bias
		for (size_t i = 0; i < forest.qtrees.size(); ++i) {
			//convert each tree:
			const vector<QRF_ResNode> &tr = forest.qtrees[i].qnodes;
			for (size_t j = 0; j < tr.size(); ++j)
			{
				generic_tree_model.children_left[pos_in_model + j] = tr[j].left;
				generic_tree_model.children_right[pos_in_model + j] = tr[j].right;
				generic_tree_model.children_default[pos_in_model + j] = tr[j].left; //smaller than is left
																					//generic_tree_model.children_default[pos_in_model + j] = -1; //no default - will fail

				generic_tree_model.features[pos_in_model + j] = int(tr[j].ifeat);
				if (tr[j].is_leaf)
				{
					generic_tree_model.children_right[pos_in_model + j] = -1;//mark leaf
					generic_tree_model.children_left[pos_in_model + j] = -1;//mark leaf
					generic_tree_model.children_default[pos_in_model + j] = -1;//mark leaf
				}
				generic_tree_model.thresholds[pos_in_model + j] = tr[j].split_val;
				generic_tree_model.node_sample_weights[pos_in_model + j] = float(tr[j].n_size); //no support for trained in weights for now
																								//if n_categ ==2:
				if (forest.n_categ == 2) {
					if (!tr[j].counts.empty()) {
						if (tr[j].n_size != 0)
							generic_tree_model.values[pos_in_model + j] = tr[j].counts[1] / float(tr[j].n_size);
						else
							MTHROW_AND_ERR("Error node leaf has 0 obs\n");
					}
					else
						generic_tree_model.values[pos_in_model + j] = tr[j].pred;
				}
				else {
					vector<float> scores(num_outputs);
					tr[j].get_scores(forest.mode, forest.get_counts_flag, forest.n_categ, scores);
					//convert values for each prediction:
					for (size_t k = 0; k < scores.size(); ++k)
						generic_tree_model.values[pos_in_model * num_outputs + j * num_outputs + k] = scores[k];
				}
			}

			pos_in_model += max_nodes;
		}

		return true;
	}

	return false;
}

TreeExplainerMode TreeExplainer::get_mode() const {
	if (proxy_predictor != NULL)
		return TreeExplainerMode::PROXY_IMPL;
	else if (generic_tree_model.is_allocate)
		return TreeExplainerMode::CONVERTED_TREES_IMPL;
	else if (original_predictor != NULL)
		return TreeExplainerMode::ORIGINAL_IMPL;
	else
		MTHROW_AND_ERR("Error TreeExplainer::get_mode() - unspecified mode. have you called init with predicotr first?");
}

void TreeExplainer::_init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "proxy_model_type")
			proxy_model_type = it->second;
		else if (it->first == "proxy_model_init")
			proxy_model_init = it->second;
		else if (it->first == "interaction_shap")
			interaction_shap = stoi(it->second) > 0;
		else if (it->first == "approximate")
			approximate = stoi(it->second) > 0;
		else if (it->first == "missing_value")
			missing_value = med_stof(it->second);
		else
			MTHROW_AND_ERR("Error in TreeExplainer::init - Unsupported parameter \"%s\"\n", it->first.c_str());
	}
}

void TreeExplainer::post_deserialization() {
	if (original_predictor != NULL) {
		if (original_predictor->classifier_type == MODEL_XGB) {
			const int PRED_CONTRIBS = 4, APPROX_CONTRIBS = 8, INTERACTION_SHAP = 16;
			if (interaction_shap)
				static_cast<MedXGB *>(original_predictor)->feat_contrib_flags = PRED_CONTRIBS | INTERACTION_SHAP;
			if (approximate)
				static_cast<MedXGB *>(original_predictor)->feat_contrib_flags |= APPROX_CONTRIBS;
		}
		try_convert_trees();
	}
}

void TreeExplainer::init_post_processor(MedModel& model) {
	ModelExplainer::init_post_processor(model);
	post_deserialization();
}

void TreeExplainer::_learn(const MedFeatures &train_mat) {
	if (processing.group2Inds.size() != train_mat.data.size() && processing.group_by_sum == 0) {
		processing.group_by_sum = 1;
		MWARN("Warning in TreeExplainer::Learn - no support for grouping in tree_shap not by sum. setting {group_by_sum:=1}\n");
	}

	if (original_predictor->classifier_type == MODEL_XGB) {
		const int PRED_CONTRIBS = 4, APPROX_CONTRIBS = 8, INTERACTION_SHAP = 16;
		if (interaction_shap)
			static_cast<MedXGB *>(original_predictor)->feat_contrib_flags = PRED_CONTRIBS | INTERACTION_SHAP;
		if (approximate)
			static_cast<MedXGB *>(original_predictor)->feat_contrib_flags |= APPROX_CONTRIBS;

		return; // no need to learn - will use XGB
	}
	//TODO: update LightGBM to support this
	if (original_predictor->classifier_type == MODEL_LIGHTGBM)
		return; // no need to learn - will use LigthGBM

	if (try_convert_trees())
		return; //success in convert to trees

				//Train XGboost model on model output.
	proxy_predictor = MedPredictor::make_predictor(proxy_model_type, proxy_model_init);
	//learn regression on input - TODO

	vector<float> labels_reg(train_mat.samples.size());
	MedMat<float> train_m;
	train_mat.get_as_matrix(train_m);
	if (proxy_predictor->transpose_for_predict != (train_m.transposed_flag > 0))
		train_m.transpose();
	original_predictor->predict(train_m, labels_reg);
	if (proxy_predictor->transpose_for_learn != (train_m.transposed_flag > 0))
		train_m.transpose();
	//proxy_predictor->prepare_x_mat(train_m);
	proxy_predictor->learn(train_m, labels_reg);
}

void conv_to_vec(const MedMat<float> &feat_contrib, vector<map<string, float>> &sample_explain_reasons) {
	sample_explain_reasons.resize(feat_contrib.nrows);
	for (int i = 0; i < sample_explain_reasons.size(); ++i)
	{
		map<string, float> &curr_sample_res = sample_explain_reasons[i];
		for (int j = 0; j < feat_contrib.ncols; ++j)
			curr_sample_res[feat_contrib.signals[j]] = feat_contrib(i, j);
	}
}

void TreeExplainer::explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const {
	TreeExplainerMode md = get_mode();
	MedMat<float> x_mat, feat_res;
	matrix.get_as_matrix(x_mat);

	vector<tfloat> shap_res; ///of size: Sample_Count, Features_count + 1(for bias/prior score), outputs_count
	ExplanationDataset data_set;
	vector<double> x, y, R;
	unique_ptr<bool> x_missing, R_missing;
	int M = x_mat.ncols;
	int num_outputs = original_predictor->n_preds_per_sample();
	int sel_output_channel = 0;
	if (num_outputs > 1)
		MWARN("Warning in TreeExplainer::explain - has several prediction channels(%d) - will explain only %d\n",
			num_outputs, sel_output_channel);

	string bias_name = "Prior_Score";
	if (md == CONVERTED_TREES_IMPL) {
		x.resize(x_mat.m.size());
		y.resize(matrix.samples.size());

		x_missing = unique_ptr<bool>(new bool[x_mat.m.size()]);
		for (size_t i = 0; i < x_mat.m.size(); ++i)
		{
			x[i] = (double)x_mat.m[i];
			//x_missing.get()[i] = x_mat.m[i] == missing_value;
			x_missing.get()[i] = false; //In QRF no missing values
		}
		for (size_t i = 0; i < matrix.samples.size(); ++i)
			y[i] = (double)matrix.samples[i].outcome;
		int num_X = (int)y.size();


		//TODO: init R
		//R.resize(x_mat.m.size());
		tfloat *R_p = NULL; // R.data()
		R_missing = NULL;
		//R_missing = unique_ptr<bool>(new bool[x_mat.m.size()]);
		int num_R = 0;

		data_set = ExplanationDataset(x.data(), x_missing.get(), y.data(), R_p, R_missing.get(), num_X, M, num_R);
		shap_res.resize(num_X * (M + 1)* num_outputs);
	}

	int tree_dep = FEATURE_DEPENDENCE::tree_path_dependent; //global is not supported in python - so not completed yet. indepent is usefull for complex transform, but can't be run with interaction
	int tranform = MODEL_TRANSFORM::identity; //this will explain raw score, the rest are use to explain loss/probability or some tranformation, based on model return function

	switch (md)
	{
	case ORIGINAL_IMPL:
		original_predictor->calc_feature_contribs(x_mat, feat_res);
		conv_to_vec(feat_res, sample_explain_reasons);
		break;
	case PROXY_IMPL:
		proxy_predictor->calc_feature_contribs(x_mat, feat_res);
		conv_to_vec(feat_res, sample_explain_reasons);
		break;
	case CONVERTED_TREES_IMPL:
		if (!approximate)
			dense_tree_shap(generic_tree_model, data_set, shap_res.data(), tree_dep, tranform, interaction_shap);
		else
			dense_tree_saabas(shap_res.data(), generic_tree_model, data_set);
		sample_explain_reasons.resize(matrix.samples.size());
		for (size_t i = 0; i < sample_explain_reasons.size(); ++i)
		{
			map<string, float> &curr_exp = sample_explain_reasons[i];
			tfloat *curr_res_exp = &shap_res[i * (M + 1) * num_outputs];
			//do only for sel_output_channel - the rest isn't supported yet
			for (size_t j = 0; j < M; ++j)
			{
				string &feat_name = x_mat.signals[j];
				curr_exp[feat_name] = (float)curr_res_exp[j * num_outputs + sel_output_channel];
			}
			curr_exp[bias_name] = (float)curr_res_exp[M * num_outputs + sel_output_channel];
		}
		break;
	default:
		MTHROW_AND_ERR("Error TreeExplainer::explain - Unsuppotrted mode %d\n", md);
	}
}

TreeExplainer::~TreeExplainer() {
	//TODO: use uniqe_ptr and than can remove those destructors..
	if (proxy_predictor != NULL) {
		delete proxy_predictor;
		proxy_predictor = NULL;
	}
	//generic_tree_model.free(); //points to existing memory in QRF, tree. not need to handle
}

MissingShapExplainer::~MissingShapExplainer() {
	//TODO: use uniqe_ptr and than can remove those destructors..
	if (retrain_predictor != NULL && !no_relearn) {
		delete retrain_predictor;
		retrain_predictor = NULL;
	}
}

template<typename T> int msn_count(const T *vals, int sz, T val) {
	int res = 0;
	for (size_t i = 0; i < sz; ++i)
		res += int(vals[i] == val);
	return res;
}

MissingShapExplainer::MissingShapExplainer() {
	processor_type = FTR_POSTPROCESS_MISSING_SHAP;
	max_test = 500;
	missing_value = MED_MAT_MISSING_VALUE;
	sample_masks_with_repeats = false;
	select_from_all = (float)0.8;
	uniform_rand = false;
	use_shuffle = false;
	add_new_data = 0;
	change_learn_args = "";
	verbose_learn = true;
	no_relearn = false;
	avg_bias_score = 0;
}

void MissingShapExplainer::_init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "missing_value")
			missing_value = med_stof(it->second);
		else if (it->first == "max_test")
			max_test = med_stoi(it->second);
		else if (it->first == "no_relearn")
			no_relearn = med_stoi(it->second) > 0;
		else if (it->first == "sample_masks_with_repeats")
			sample_masks_with_repeats = med_stoi(it->second) > 0;
		else if (it->first == "uniform_rand")
			uniform_rand = med_stoi(it->second) > 0;
		else if (it->first == "use_shuffle")
			use_shuffle = med_stoi(it->second) > 0;
		else if (it->first == "select_from_all")
			select_from_all = med_stof(it->second);
		else if (it->first == "add_new_data")
			add_new_data = med_stoi(it->second);
		else if (it->first == "change_learn_args")
			change_learn_args = it->second;
		else if (it->first == "verbose_learn")
			verbose_learn = stoi(it->second) > 0;
		else
			MTHROW_AND_ERR("Error SHAPExplainer::init - Unknown param \"%s\"\n", it->first.c_str());
	}
}

float get_avg_preds(const MedFeatures &train_mat, MedPredictor *original_predictor) {
	if (train_mat.samples.empty())
		MTHROW_AND_ERR("Error get_avg_preds learn matrix is empty\n");
	vector<float> preds_orig(train_mat.samples.size());
	float avg_bias_score = 0;
	if (train_mat.samples.front().prediction.empty()) {
		MedMat<float> mat_x;
		train_mat.get_as_matrix(mat_x);
		if (original_predictor->transpose_for_predict != (mat_x.transposed_flag > 0))
			mat_x.transpose();
		original_predictor->predict(mat_x, preds_orig);
	}
	else {
		for (size_t i = 0; i < preds_orig.size(); ++i)
			preds_orig[i] = train_mat.samples[i].prediction[0];
	}
	for (size_t i = 0; i < preds_orig.size(); ++i)
		avg_bias_score += preds_orig[i];
	avg_bias_score /= preds_orig.size();
	return avg_bias_score;
}

void MissingShapExplainer::_learn(const MedFeatures &train_mat) {
	avg_bias_score = get_avg_preds(train_mat, original_predictor);
	if (no_relearn) {
		retrain_predictor = original_predictor;
		return;
	}
	retrain_predictor = (MedPredictor *)medial::models::copyInfraModel(original_predictor, false);
	mt19937 gen(globalRNG::rand());
	MedMat<float> x_mat;
	train_mat.get_as_matrix(x_mat);
	int nftrs = x_mat.ncols;
	vector<float> labels(train_mat.samples.size()), weights(train_mat.samples.size() + add_new_data, 1);
	vector<int> miss_cnts(train_mat.samples.size() + add_new_data);
	vector<int> missing_hist(nftrs + 1), added_missing_hist(nftrs + 1);

	if (!train_mat.samples.front().prediction.empty())
		for (size_t i = 0; i < labels.size(); ++i)
			labels[i] = train_mat.samples[i].prediction[0];
	else {
		MedMat<float> tt;
		train_mat.get_as_matrix(tt);
		original_predictor->predict(tt, labels);
	}
	if (add_new_data > 0) {
		vector<float> rows_m(add_new_data * nftrs);
		unordered_set<vector<bool>> seen_mask;
		uniform_int_distribution<> rnd_row(0, (int)train_mat.samples.size() - 1);
		double log_max_opts = log(add_new_data) / log(2.0);
		if (log_max_opts >= nftrs) {
			if (!sample_masks_with_repeats)
				MWARN("Warning: you have request to sample masks without repeats, but it can't be done. setting sample with repeats\n");
			sample_masks_with_repeats = true;
		}
		if (verbose_learn)
			MLOG("Adding %d Data points\n", add_new_data);
		MedProgress add_progress("Add_Train_Data", add_new_data, 30, 1);
		for (size_t i = 0; i < add_new_data; ++i)
		{
			float *curr_row = &rows_m[i *  nftrs];
			//select row:
			int row_sel = rnd_row(gen);

			vector<bool> curr_mask; curr_mask.reserve(nftrs);
			for (int j = 0; j < nftrs; ++j)
				curr_mask.push_back(x_mat(row_sel, j) != missing_value);

			medial::shapley::generate_mask_(curr_mask, nftrs, gen, uniform_rand, use_shuffle);
			while (!sample_masks_with_repeats && seen_mask.find(curr_mask) != seen_mask.end())
				medial::shapley::generate_mask_(curr_mask, nftrs, gen, uniform_rand, use_shuffle);
			if (!sample_masks_with_repeats)
				seen_mask.insert(curr_mask);

			//commit mask to curr_row
			for (int j = 0; j < nftrs; ++j)
			{
				if (curr_mask[j])
					curr_row[j] = x_mat(row_sel, j);
				else
					curr_row[j] = missing_value;
			}
			labels.push_back(labels[row_sel]);
			add_progress.update();
		}
		x_mat.add_rows(rows_m);
	}

	// Add data with missing values according to sample masks
	for (size_t i = 0; i < x_mat.nrows; ++i) {
		miss_cnts[i] = msn_count<float>(&x_mat.m[i*nftrs], nftrs, missing_value);
		++missing_hist[miss_cnts[i]];
		if (i >= train_mat.samples.size())
			++added_missing_hist[miss_cnts[i]];
	}
	for (size_t i = 0; i < x_mat.nrows; ++i) {
		float curr_mask_w = x_mat.nrows / float(missing_hist[miss_cnts[i]]);
		weights[i] = curr_mask_w;
	}
	if (verbose_learn) {
		medial::print::print_hist_vec(miss_cnts, "missing_values hist", "%d");
		medial::print::print_hist_vec(added_missing_hist, "hist of added_missing_hist", "%d");
		medial::print::print_hist_vec(weights, "weights for learn", "%2.4f");
	}
	if (original_predictor->transpose_for_learn != (x_mat.transposed_flag > 0))
		x_mat.transpose();
	//reweight train_mat:
	retrain_predictor->init_from_string(change_learn_args);
	retrain_predictor->learn(x_mat, labels, weights);
	//test pref:
	if (verbose_learn) {
		vector<float> train_p;
		retrain_predictor->predict(x_mat, train_p);
		float rmse = medial::performance::rmse_without_cleaning(train_p, labels, &weights);
		MLOG("RMSE=%2.4f on train for model\n", rmse);
	}
}

void MissingShapExplainer::explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const {
	MedPredictor *predictor = retrain_predictor;
	if (no_relearn)
		predictor = original_predictor;
	const vector<vector<int>> *group_inds = &processing.group2Inds;
	const vector<string> *group_names = &processing.groupNames;
	vector<vector<int>> group_inds_loc;
	vector<string> group_names_loc;
	if (processing.group_by_sum) {
		int icol = 0;
		for (auto& rec : matrix.data) {
			group_inds_loc.push_back({ icol++ });
			group_names_loc.push_back(rec.first);
		}
		group_inds = &group_inds_loc;
		group_names = &group_names_loc;
	}

	sample_explain_reasons.resize(matrix.samples.size());
	string bias_name = "Prior_Score";
	vector<float> preds_orig(matrix.samples.size());
	if (matrix.samples.front().prediction.empty()) {
		MedMat<float> mat_x;
		matrix.get_as_matrix(mat_x);
		if (predictor->transpose_for_predict != (mat_x.transposed_flag > 0))
			mat_x.transpose();
		predictor->predict(mat_x, preds_orig);
	}
	else {
		for (size_t i = 0; i < preds_orig.size(); ++i)
			preds_orig[i] = matrix.samples[i].prediction[0];
	}
	int N_TOTAL_TH = omp_get_max_threads();

	MedProgress progress("MissingShapley", (int)matrix.samples.size(), 15);
#pragma omp parallel for if (matrix.samples.size() >= N_TOTAL_TH)
	for (int i = 0; i < matrix.samples.size(); ++i)
	{
		vector<float> features_coeff;
		float pred_shap = 0;
		medial::shapley::explain_shapley(matrix, (int)i, max_test, predictor, missing_value, *group_inds, *group_names, features_coeff,
			sample_masks_with_repeats, select_from_all, uniform_rand, use_shuffle, global_logger.levels[LOCAL_SECTION] < LOG_DEF_LEVEL);

		for (size_t j = 0; j < features_coeff.size(); ++j)
			pred_shap += features_coeff[j];

#pragma omp critical 
		{
			map<string, float> &curr_res = sample_explain_reasons[i];
			for (size_t j = 0; j < group_names->size(); ++j)
				curr_res[group_names->at(j)] = features_coeff[j];
			//Add prior to score:
			//curr_res[bias_name] = preds_orig[i] - pred_shap; //that will sum to current score
			curr_res[bias_name] = avg_bias_score; //that will sum to current score
		}

		progress.update();
	}
}

string GeneratorType_toStr(GeneratorType type) {
	switch (type)
	{
	case GIBBS:
		return "GIBBS";
	case GAN:
		return "GAN";
	case MISSING:
		return "MISSING";
	case RANDOM_DIST:
		return "RANDOM_DIST";
	default:
		MTHROW_AND_ERR("Unknown type %d\n", type);
	}
}
GeneratorType GeneratorType_fromStr(const string &type) {
	string tp = boost::to_upper_copy(type);
	if (tp == "GAN")
		return GeneratorType::GAN;
	else if (tp == "GIBBS")
		return GeneratorType::GIBBS;
	else if (tp == "MISSING")
		return GeneratorType::MISSING;
	else if (tp == "RANDOM_DIST")
		return GeneratorType::RANDOM_DIST;
	else
		MTHROW_AND_ERR("Unknown type %s\n", type.c_str());
}

void ShapleyExplainer::_init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "gen_type")
			gen_type = GeneratorType_fromStr(it->second);
		else if (it->first == "generator_args")
			generator_args = it->second;
		else if (it->first == "missing_value")
			missing_value = med_stof(it->second);
		else if (it->first == "n_masks")
			n_masks = med_stoi(it->second);
		else if (it->first == "sampling_args")
			sampling_args = it->second;
		else if (it->first == "use_random_sampling")
			use_random_sampling = med_stoi(it->second) > 0;
		else
			MTHROW_AND_ERR("Error in ShapleyExplainer::init - Unsupported param \"%s\"\n", it->first.c_str());
	}
	init_sampler(); //from args
}

void ShapleyExplainer::init_sampler(bool with_sampler) {
	switch (gen_type)
	{
	case GIBBS:
		if (with_sampler) {
			_gibbs.init_from_string(generator_args);
			_sampler = unique_ptr<SamplesGenerator<float>>(new GibbsSamplesGenerator<float>(_gibbs, true));
		}
		_gibbs_sample_params.init_from_string(sampling_args);
		sampler_sampling_args = &_gibbs_sample_params;
		break;
	case GAN:
		if (with_sampler) {
			_sampler = unique_ptr<SamplesGenerator<float>>(new MaskedGAN<float>);
			static_cast<MaskedGAN<float> *>(_sampler.get())->read_from_text_file(generator_args);
		}
		break;
	case MISSING:
		if (with_sampler)
			_sampler = unique_ptr<SamplesGenerator<float>>(new MissingsSamplesGenerator<float>(missing_value));
		break;
	case RANDOM_DIST:
		if (with_sampler)
			_sampler = unique_ptr<SamplesGenerator<float>>(new RandomSamplesGenerator<float>(0, 5));
		sampler_sampling_args = &n_masks;
		break;
	default:
		MTHROW_AND_ERR("Error in ShapleyExplainer::init_sampler() - Unsupported Type %d\n", gen_type);
	}
}

void ShapleyExplainer::_learn(const MedFeatures &train_mat) {
	_sampler->learn(train_mat.data);
	avg_bias_score = get_avg_preds(train_mat, original_predictor);
}

void ShapleyExplainer::explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const {
	sample_explain_reasons.resize(matrix.samples.size());
	string bias_name = "Prior_Score";
	vector<float> preds_orig(matrix.samples.size());
	if (matrix.samples.front().prediction.empty()) {
		MedMat<float> mat_x;
		matrix.get_as_matrix(mat_x);
		if (original_predictor->transpose_for_predict != (mat_x.transposed_flag > 0))
			mat_x.transpose();
		original_predictor->predict(mat_x, preds_orig);
	}
	else {
		for (size_t i = 0; i < preds_orig.size(); ++i)
			preds_orig[i] = matrix.samples[i].prediction[0];
	}

	const vector<vector<int>> *group_inds = &processing.group2Inds;
	const vector<string> *group_names = &processing.groupNames;
	vector<vector<int>> group_inds_loc;
	vector<string> group_names_loc;
	if (processing.group_by_sum) {
		int icol = 0;
		for (auto& rec : matrix.data) {
			group_inds_loc.push_back({ icol++ });
			group_names_loc.push_back(rec.first);
		}
		group_inds = &group_inds_loc;
		group_names = &group_names_loc;
	}

	int MAX_Threads = omp_get_max_threads();
	//copy sample for each thread:
	random_device rd;
	vector<mt19937> gen_thread(MAX_Threads);
	for (size_t i = 0; i < gen_thread.size(); ++i)
		gen_thread[i] = mt19937(globalRNG::rand());
	_sampler->prepare(sampler_sampling_args);

	MedProgress progress("ShapleyExplainer", (int)matrix.samples.size(), 15);
#pragma omp parallel for if (matrix.samples.size() >= 2)
	for (int i = 0; i < matrix.samples.size(); ++i)
	{
		int n_th = omp_get_thread_num();
		vector<float> features_coeff;
		float pred_shap = 0;
		medial::shapley::explain_shapley(matrix, (int)i, n_masks, original_predictor
			, *group_inds, *group_names, *_sampler, gen_thread[n_th], 1, sampler_sampling_args, features_coeff,
			use_random_sampling, global_logger.levels[LOCAL_SECTION] < LOCAL_LEVEL &&
			(!(matrix.samples.size() >= 2) || omp_get_thread_num() == 1));

		for (size_t j = 0; j < features_coeff.size(); ++j)
			pred_shap += features_coeff[j];

#pragma omp critical 
		{
			map<string, float> &curr_res = sample_explain_reasons[i];
			for (size_t j = 0; j < group_names->size(); ++j)
				curr_res[group_names->at(j)] = features_coeff[j];
			//Add prior to score:
			//curr_res[bias_name] = preds_orig[i] - pred_shap; //that will sum to current score
			curr_res[bias_name] = avg_bias_score;
		}

		progress.update();
	}
}

void ShapleyExplainer::post_deserialization() {
	init_sampler(false);
}

void ShapleyExplainer::load_GIBBS(MedPredictor *original_pred, const GibbsSampler<float> &gibbs, const GibbsSamplingParams &sampling_args) {
	this->original_predictor = original_pred;
	_gibbs = gibbs;
	_gibbs_sample_params = sampling_args;

	sampler_sampling_args = &_gibbs_sample_params;
	_sampler = unique_ptr<SamplesGenerator<float>>(new GibbsSamplesGenerator<float>(_gibbs, true));

	gen_type = GeneratorType::GIBBS;
}

void ShapleyExplainer::load_GAN(MedPredictor *original_pred, const string &gan_path) {
	this->original_predictor = original_pred;
	_sampler = unique_ptr<SamplesGenerator<float>>(new MaskedGAN<float>);
	static_cast<MaskedGAN<float> *>(_sampler.get())->read_from_text_file(gan_path);

	gen_type = GeneratorType::GAN;
}

void ShapleyExplainer::load_MISSING(MedPredictor *original_pred) {
	this->original_predictor = original_pred;
	_sampler = unique_ptr<SamplesGenerator<float>>(new MissingsSamplesGenerator<float>(missing_value));
	gen_type = GeneratorType::MISSING;
}

void ShapleyExplainer::load_sampler(MedPredictor *original_pred, unique_ptr<SamplesGenerator<float>> &&generator) {
	this->original_predictor = original_pred;
	_sampler = move(generator);
}

void ShapleyExplainer::dprint(const string &pref) const {
	string predictor_nm = "";
	if (original_predictor != NULL)
		predictor_nm = original_predictor->my_class_name();
	string filters_str = "", processing_str = "";
	char buffer[5000];
	snprintf(buffer, sizeof(buffer), "group_by_sum=%d, learn_cov_matrix=%d, normalize_vals=%d, zero_missing=%d, grouping=%s",
		int(processing.group_by_sum), int(processing.learn_cov_matrix), processing.normalize_vals
		, processing.zero_missing, processing.grouping.c_str());
	processing_str = string(buffer);
	snprintf(buffer, sizeof(buffer), "sort_mode=%d, max_count=%d, sum_ratio=%2.3f",
		filters.sort_mode, filters.max_count, filters.sum_ratio);
	filters_str = string(buffer);

	MLOG("%s :: ModelExplainer type %d(%s), original_predictor=%s, gen_type=%s, attr_name=%s, processing={%s}, filters={%s}\n",
		pref.c_str(), processor_type, my_class_name().c_str(), predictor_nm.c_str(),
		GeneratorType_toStr(gen_type).c_str(), attr_name.c_str(),
		processing_str.c_str(), filters_str.c_str());
}

void LimeExplainer::_init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "gen_type")
			gen_type = GeneratorType_fromStr(it->second);
		else if (it->first == "generator_args")
			generator_args = it->second;
		else if (it->first == "missing_value")
			missing_value = med_stof(it->second);
		else if (it->first == "sampling_args")
			sampling_args = it->second;
		else if (it->first == "p_mask")
			p_mask = med_stof(it->second);
		else if (it->first == "weight")
			weighting = get_weight_method(it->second);
		else if (it->first == "n_masks")
			n_masks = med_stoi(it->second);
		else
			MTHROW_AND_ERR("Error in LimeExplainer::init - Unsupported param \"%s\"\n", it->first.c_str());
	}
	init_sampler(); //from args
}

medial::shapley::LimeWeightMethod LimeExplainer::get_weight_method(string method_s) {

	boost::to_lower(method_s);
	if (method_s == "lime") return medial::shapley::LimeWeightLime;
	if (method_s == "unif" || method_s == "uniform") return medial::shapley::LimeWeightUniform;
	if (method_s == "shap" || method_s == "shapley") return medial::shapley::LimeWeightShap;
	if (method_s == "sum" || method_s == "shap_sum" || method_s == "shapley_sum") return medial::shapley::LimeWeightSum;

	MTHROW_AND_ERR("Unknown weighting method %s for LIME explainer\n", method_s.c_str());
	return medial::shapley::LimeWeightLast;
}

void LimeExplainer::init_sampler(bool with_sampler) {
	switch (gen_type)
	{
	case GIBBS:
		if (with_sampler) {
			_gibbs.init_from_string(generator_args);
			_sampler = unique_ptr<SamplesGenerator<float>>(new GibbsSamplesGenerator<float>(_gibbs, true));
		}
		_gibbs_sample_params.init_from_string(sampling_args);
		sampler_sampling_args = &_gibbs_sample_params;
		break;
	case GAN:
		if (with_sampler) {
			_sampler = unique_ptr<SamplesGenerator<float>>(new MaskedGAN<float>);
			static_cast<MaskedGAN<float> *>(_sampler.get())->read_from_text_file(generator_args);
		}
		break;
	case MISSING:
		if (with_sampler)
			_sampler = unique_ptr<SamplesGenerator<float>>(new MissingsSamplesGenerator<float>(missing_value));
		break;
	default:
		MTHROW_AND_ERR("Error in LimeExplainer::init_sampler() - Unsupported Type %d\n", gen_type);
	}
}

void LimeExplainer::load_GIBBS(MedPredictor *original_pred, const GibbsSampler<float> &gibbs, const GibbsSamplingParams &sampling_args) {
	this->original_predictor = original_pred;
	_gibbs = gibbs;
	_gibbs_sample_params = sampling_args;

	sampler_sampling_args = &_gibbs_sample_params;
	_sampler = unique_ptr<SamplesGenerator<float>>(new GibbsSamplesGenerator<float>(_gibbs, true));

	gen_type = GeneratorType::GIBBS;
}

void LimeExplainer::load_GAN(MedPredictor *original_pred, const string &gan_path) {
	this->original_predictor = original_pred;
	_sampler = unique_ptr<SamplesGenerator<float>>(new MaskedGAN<float>);
	static_cast<MaskedGAN<float> *>(_sampler.get())->read_from_text_file(gan_path);

	gen_type = GeneratorType::GAN;
}

void LimeExplainer::load_MISSING(MedPredictor *original_pred) {
	this->original_predictor = original_pred;
	_sampler = unique_ptr<SamplesGenerator<float>>(new MissingsSamplesGenerator<float>(missing_value));
	gen_type = GeneratorType::MISSING;
}

void LimeExplainer::load_sampler(MedPredictor *original_pred, unique_ptr<SamplesGenerator<float>> &&generator) {
	this->original_predictor = original_pred;
	_sampler = move(generator);
}

void LimeExplainer::post_deserialization() {
	init_sampler(false);
}

void LimeExplainer::_learn(const MedFeatures &train_mat) {
	_sampler->learn(train_mat.data);
}

void LimeExplainer::explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const {
	vector<vector<float>> alphas;

	const vector<vector<int>> *group_inds = &processing.group2Inds;
	const vector<string> *group_names = &processing.groupNames;
	vector<vector<int>> group_inds_loc;
	vector<string> group_names_loc;
	if (processing.group_by_sum) {
		int icol = 0;
		for (auto& rec : matrix.data) {
			group_inds_loc.push_back({ icol++ });
			group_names_loc.push_back(rec.first);
		}
		group_inds = &group_inds_loc;
		group_names = &group_names_loc;
	}

	medial::shapley::get_shapley_lime_params(matrix, original_predictor, _sampler.get(), p_mask, n_masks, weighting, missing_value,
		sampler_sampling_args, *group_inds, *group_names, alphas);

	sample_explain_reasons.resize(matrix.samples.size());

	for (size_t i = 0; i < sample_explain_reasons.size(); ++i)
	{
		map<string, float> &curr = sample_explain_reasons[i];
		const vector<float> &curr_res = alphas[i];
		for (size_t k = 0; k < group_names->size(); ++k)
			curr[group_names->at(k)] = curr_res[k];
	}
}

void LimeExplainer::dprint(const string &pref) const {
	string predictor_nm = "";
	if (original_predictor != NULL)
		predictor_nm = original_predictor->my_class_name();
	string filters_str = "", processing_str = "";
	char buffer[5000];
	snprintf(buffer, sizeof(buffer), "group_by_sum=%d, learn_cov_matrix=%d, normalize_vals=%d, zero_missing=%d, grouping=%s",
		int(processing.group_by_sum), int(processing.learn_cov_matrix), processing.normalize_vals
		, processing.zero_missing, processing.grouping.c_str());
	processing_str = string(buffer);
	snprintf(buffer, sizeof(buffer), "sort_mode=%d, max_count=%d, sum_ratio=%2.3f",
		filters.sort_mode, filters.max_count, filters.sum_ratio);
	filters_str = string(buffer);

	MLOG("%s :: ModelExplainer type %d(%s), original_predictor=%s, gen_type=%s, attr_name=%s, processing={%s}, filters={%s}\n",
		pref.c_str(), processor_type, my_class_name().c_str(), predictor_nm.c_str(),
		GeneratorType_toStr(gen_type).c_str(), attr_name.c_str(),
		processing_str.c_str(), filters_str.c_str());
}

void LinearExplainer::_init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//no arguments so far
		MTHROW_AND_ERR("Error in LinearExplainer::init - Unsupported param \"%s\"\n", it->first.c_str());
	}
}

void LinearExplainer::_learn(const MedFeatures &train_mat) {
	avg_bias_score = get_avg_preds(train_mat, original_predictor);
}

void LinearExplainer::explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const {
	sample_explain_reasons.resize(matrix.data.begin()->second.size());
	string bias_name = "Prior_Score";
	vector<float> preds_orig(matrix.samples.size());
	if (matrix.samples.front().prediction.empty()) {
		MedMat<float> mat_x;
		matrix.get_as_matrix(mat_x);
		if (original_predictor->transpose_for_predict != (mat_x.transposed_flag > 0))
			mat_x.transpose();
		original_predictor->predict(mat_x, preds_orig);
	}
	else {
		for (size_t i = 0; i < preds_orig.size(); ++i)
			preds_orig[i] = matrix.samples[i].prediction[0];
	}

	const vector<vector<int>> *group_inds = &processing.group2Inds;
	const vector<string> *group_names = &processing.groupNames;
	vector<vector<int>> group_inds_loc;
	vector<string> group_names_loc;
	if (processing.group_by_sum) {
		int icol = 0;
		for (auto& rec : matrix.data) {
			group_inds_loc.push_back({ icol++ });
			group_names_loc.push_back(rec.first);
		}
		group_inds = &group_inds_loc;
		group_names = &group_names_loc;
	}

	MedProgress progress("LinearExplainer", (int)matrix.samples.size(), 15);

	vector<float> x(matrix.samples.size() * matrix.data.size());
	for (int i = 0; i < matrix.samples.size(); ++i) {
		int j = 0;
		for (auto it = matrix.data.begin(); it != matrix.data.end(); ++it) {
			x[i* matrix.data.size() + j] = it->second[i];
			++j;
		}
	}


	vector<vector<float>> all_features_coeff(group_names->size());
	//no parallel - will happen in predict
	for (int i = 0; i < group_names->size(); ++i)
	{
		//put zeros in mask i
		vector<float> masked_x = x;
		for (int ind : group_inds->at(i))
			for (size_t j = 0; j < matrix.samples.size(); ++j)
				masked_x[j * matrix.data.size() + ind] = 0;
		vector<float> preds_masked;
		original_predictor->predict(masked_x, preds_masked, (int)matrix.samples.size(), (int)matrix.data.size());

		//commit:
#pragma omp critical
		{
			all_features_coeff[i].resize(matrix.samples.size());
			for (size_t j = 0; j < matrix.samples.size(); ++j)
				all_features_coeff[i][j] = preds_orig[j] - preds_masked[j];
		}

		progress.update();
	}

	//commit to memory:
	for (int i = 0; i < sample_explain_reasons.size(); ++i)
	{
		map<string, float> &curr_res = sample_explain_reasons[i];

		float pred_shap = 0;
		for (size_t j = 0; j < group_names->size(); ++j)
			pred_shap += all_features_coeff[j][i];

		for (size_t j = 0; j < group_names->size(); ++j)
			curr_res[group_names->at(j)] = all_features_coeff[j][i];
		//Add prior to score:
		//curr_res[bias_name] = preds_orig[i] - pred_shap; //that will sum to current score
		curr_res[bias_name] = avg_bias_score; //that will sum to current score
	}
}

void KNN_Explainer::_init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "fraction")
			fraction = med_stof(it->second);
		else if (it->first == "numClusters")
			numClusters = med_stoi(it->second);
		else if (it->first == "chosenThreshold")
			chosenThreshold = med_stof(it->second);
		else if (it->first == "thresholdQ")
			thresholdQ = med_stof(it->second);
		else
			MTHROW_AND_ERR("Error in KNN_Explainer::init - Unsupported param \"%s\"\n", it->first.c_str());
	}
}
void KNN_Explainer::_learn(const MedFeatures &train_mat) {
	if (numClusters == -1)numClusters = (int)train_mat.samples.size();
	if (numClusters > train_mat.samples.size()) {
		MWARN("Warning in KNN_Explainer::Learn - numClusters reduced to size of training \"%d>>%zu\"\n", numClusters, train_mat.samples.size());
		numClusters = (int)train_mat.samples.size();
	}

	MedMat<float> centers(numClusters, (int)train_mat.data.size());

	// get the features and normalize them
	MedFeatures normalizedFeatures = train_mat;
	MedMat<float> normalizedMatrix;
	normalizedFeatures.get_as_matrix(normalizedMatrix);

	normalizedMatrix.normalize();
	normalizedFeatures.set_as_matrix(normalizedMatrix);
	// keep normalization params for future use in apply
	average = normalizedMatrix.avg;
	std = normalizedMatrix.std;
	/* we will test on kmeans later because it is unstable
		//represent features by limitted number of clusters
		vector<int>clusters;
		MedMat<float>dists;
		KMeans(normalizedMatrix, numClusters, centers, clusters, dists);
		vector <float> weights(centers.nrows,0);
		for (int i=0;i<clusters.size();i++){
			weights[clusters[i]]+=1./clusters.size();
		}
		*/

		// random sample of space
	vector <int>krand;
	for (int k = 0; k < normalizedMatrix.nrows; k++)
		krand.push_back(k);
	shuffle(krand.begin(), krand.end(), default_random_engine(5246245));

	for (int i = 0; i < numClusters; i++)
		for (int col = 0; col < normalizedMatrix.ncols; col++)
			centers(i, col) = normalizedMatrix(krand[i], col);
	centers.signals = normalizedMatrix.signals;
	vector<float> weights(centers.nrows, 1);


	//keep the features for the apply phase
	trainingMap.set_as_matrix(centers);
	trainingMap.weights = weights;
	trainingMap.samples.resize(numClusters);
	trainingMap.init_pid_pos_len();

	// compute the thershold according to quantile
	MedFeatures myMat = train_mat;// train_mat is constant
	this->original_predictor->predict(myMat);

	//assign predictions to the sampled  features
	for (int i = 0; i < numClusters; i++)
		trainingMap.samples[i].prediction = vector <float>(1, myMat.samples[krand[i]].prediction[0]);

	// compute the thershold according to quantile
	if (chosenThreshold == MED_MAT_MISSING_VALUE) {
		if (thresholdQ != MED_MAT_MISSING_VALUE) {
			vector <float> predictions = {};
			vector <float> w(train_mat.samples.size(), 1);
			for (int k = 0; k < train_mat.samples.size(); k++)
				predictions.push_back(myMat.samples[k].prediction[0]);

			chosenThreshold = medial::stats::get_quantile(predictions, w, 1 - thresholdQ);
		}
	}
}

void KNN_Explainer::explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const
{
	MedFeatures explainedFeatures = matrix;
	MedMat<float> explainedMatrix, explainedMatrixCopy;

	vector <string> featureNames;
	trainingMap.get_feature_names(featureNames);

	// check if grouping required if not prepare knn groups from single features
	vector <vector<int>> knnGroups; // will hold the given groups from processing or single features group if not given
	vector<string> knnGroupNames;
	if ((processing.group2Inds.size() > 0) && (processing.group_by_sum == 0)) {
		knnGroups = processing.group2Inds;
		knnGroupNames = processing.groupNames;
	}
	else
		for (int col = 0; col < trainingMap.data.size(); col++) {
			knnGroups.push_back(vector<int>{col});
			knnGroupNames.push_back(featureNames[col]);
		}
	//normalize the explained features
	explainedFeatures.get_as_matrix(explainedMatrix);
	explainedFeatures.get_as_matrix(explainedMatrixCopy);// keep it to handle the missing
	MedMat <float> trainingCentersMatrix;
	trainingMap.get_as_matrix(trainingCentersMatrix);
	sample_explain_reasons = {};

	explainedMatrix.normalize(average, std, 1);
	for (int row = 0; row < explainedMatrix.nrows; row++)
		for (int col = 0; col < explainedMatrix.ncols; col++)
			if (explainedMatrixCopy.get(row, col) == MED_MAT_MISSING_VALUE)
				explainedMatrix.set(row, col) = MED_MAT_MISSING_VALUE;
	//for each sample compute the explanation
	vector <float> thisRow;
	for (int row = 0; row < explainedMatrix.nrows; row++) {
		explainedMatrix.get_row(row, thisRow);
		sample_explain_reasons.push_back({});
		computeExplanation(thisRow, sample_explain_reasons[row], knnGroups, knnGroupNames);
	}
}
void KNN_Explainer::computeExplanation(vector<float> thisRow, map<string, float> &sample_explain_reasons, vector <vector<int>> knnGroups, vector<string> knnGroupNames) const
// do the calculation for a single sample after normalization
{

	MedMat<float> centers; //matrix taken from features and holds the centers of clusters
	trainingMap.get_as_matrix(centers);
	MedMat<float> pDistance(centers.nrows, centers.ncols);//initialized to 0
	MedMat<float> gDistance(centers.nrows, (int)knnGroups.size());
	vector<float>totalDistance(centers.nrows, 0);
#define SQR(x)  ((x)*(x))
	for (int row = 0; row < centers.nrows; row++) {
		for (int col = 0; col < centers.ncols; col++)
			if (thisRow[col] != MED_MAT_MISSING_VALUE) {
				pDistance(row, col) = SQR(centers.get(row, col) - thisRow[col]);
				totalDistance[row] += pDistance(row, col);
			}
		for (int group = 0; group < knnGroupNames.size(); group++) {
			gDistance(row, group) = totalDistance[row];
			for (auto inGroup : knnGroups[group])
				if (thisRow[inGroup] != MED_MAT_MISSING_VALUE)
					gDistance(row, group) -= pDistance(row, inGroup);

		}
	}
	vector<float> thresholds(centers.nrows, 0);
	vector <float> colVector;
	float totalThreshold = medial::stats::get_quantile(totalDistance, trainingMap.weights, fraction);
	for (int group = 0; group < knnGroupNames.size(); group++) {
		gDistance.get_col(group, colVector);
		thresholds[group] = medial::stats::get_quantile(colVector, trainingMap.weights, fraction);
	}
	double sumWeights = 0;
	double pCol;
	double pTotal = 0;

	for (int row = 0; row < pDistance.nrows; row++)
		if (totalDistance[row] < totalThreshold) {
			float thisPred = trainingMap.samples[row].prediction[0];
			sumWeights += trainingMap.weights[row];
			if (chosenThreshold != MED_MAT_MISSING_VALUE)thisPred = thisPred > chosenThreshold;// threshold the predictions if needed
			pTotal += trainingMap.weights[row] * thisPred;
			//cout <<row<<" "<< trainingMap.samples[row].prediction[0] << "\n";
		}
	pTotal /= sumWeights;


	for (int group = 0; group < knnGroupNames.size(); group++) {
		pCol = 0;
		sumWeights = 0;
		for (int row = 0; row < gDistance.nrows; row++)
			if (gDistance.get(row, group) < thresholds[group]) {
				float thisPred = trainingMap.samples[row].prediction[0];
				if (chosenThreshold != MED_MAT_MISSING_VALUE)thisPred = thisPred > chosenThreshold;// threshold the predictions if needed
				pCol += trainingMap.weights[row] * thisPred;
				sumWeights += trainingMap.weights[row];
			}
		pCol /= sumWeights;
		sample_explain_reasons.insert(pair<string, float>(knnGroupNames[group], float(log((pTotal + 1e-10) / (pCol + 1e-10) / (1 - pTotal + 1e-10)*(1 - pCol + 1e-10)))));


	}




}