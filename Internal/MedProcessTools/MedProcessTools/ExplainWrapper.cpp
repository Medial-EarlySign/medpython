#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <cmath>
#include <random>
#include <omp.h>
#include "ExplainWrapper.h"
#include <MedAlgo/MedAlgo/MedXGB.h>

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL LOG_DEF_LEVEL

void ModelExplainer::explain(MedFeatures &matrix) const {
	vector<map<string, float>> explain_reasons; //for each sample, reasons and scores
	explain(matrix, explain_reasons);

	if (explain_reasons.size() != matrix.samples.size())
		MTHROW_AND_ERR("Error in ModelExplainer::explain - explain returned musmatch number of samples %zu, and requested %zu\n",
			explain_reasons.size(), matrix.samples.size());

	for (size_t i = 0; i < explain_reasons.size(); ++i)
		for (auto it = explain_reasons[i].begin(); it != explain_reasons[i].end(); ++it)
			matrix.samples[i].attributes["ModelExplainer::" + it->first] = it->second;

}

void ModelExplainer::Learn(MedModel &model, MedPidRepository& rep, const MedFeatures &train_mat) {
	Learn(model.predictor, train_mat);
}

bool comp_score_str(const pair<string, float> &pr1, const pair<string, float> &pr2) {
	return abs(pr1.second) > abs(pr2.second); //bigger is better in absolute
}

void ModelExplainer::print_explain(MedSample &smp) {
	vector<pair<string, float>> ranked;
	for (auto it = smp.attributes.begin(); it != smp.attributes.end(); ++it)
		if (boost::starts_with(it->first, "ModelExplainer::"))
			ranked.push_back(pair<string, float>(it->first, it->second));
	sort(ranked.begin(), ranked.end(), comp_score_str);

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

void read_feature_grouping(const string &file_name, const MedFeatures& data, vector<vector<int>>& group2index, vector<string>& group_names) {
	// Features
	vector<string> features;
	data.get_feature_names(features);
	int nftrs = (int)features.size();

	map<string, int> ftr2indx;
	for (int i = 0; i < nftrs; i++)
		ftr2indx[features[i]] = i;

	// Read Grouping
	ifstream inf(file_name);
	if (!inf.is_open())
		MTHROW_AND_ERR("Cannot open \'%s\' for reading\n", file_name.c_str());

	string curr_line;
	vector<string> fields;
	map<string, vector<int>> groups;
	unordered_set<string> grouped_ftrs;

	while (getline(inf, curr_line)) {
		boost::split(fields, curr_line, boost::is_any_of("\t"));
		if (fields.size() != 2 || ftr2indx.find(fields[0]) == ftr2indx.end())
			MTHROW_AND_ERR("Cannot parse line \'%s\' from %s\n", curr_line.c_str(), file_name.c_str());
		if (grouped_ftrs.find(fields[0]) != grouped_ftrs.end())
			MTHROW_AND_ERR("Features %s given twice\n", fields[0].c_str());

		grouped_ftrs.insert(fields[0]);
		groups[fields[1]].push_back(ftr2indx[fields[0]]);
	}

	// Arrange
	for (auto& rec : groups) {
		group_names.push_back(rec.first);
		group2index.push_back(rec.second);
	}

	for (int i = 0; i < nftrs; i++) {
		if (grouped_ftrs.find(features[i]) == grouped_ftrs.end()) {
			group_names.push_back(features[i]);
			group2index.push_back({ i });
		}
	}

	MLOG("Grouping: %d features into %d groups\n", nftrs, (int)group_names.size());
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

int TreeExplainer::init(map<string, string> &mapper) {
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
		else if (it->first == "pp_type") {}
		else
			MTHROW_AND_ERR("Error in TreeExplainer::init - Unsupported parameter \"%s\"\n", it->first.c_str());
	}
	return 0;
}

void TreeExplainer::post_deserialization() {
	if (original_predictor->classifier_type == MODEL_XGB) {
		const int PRED_CONTRIBS = 4, APPROX_CONTRIBS = 8, INTERACTION_SHAP = 16;
		if (interaction_shap)
			static_cast<MedXGB *>(original_predictor)->feat_contrib_flags = PRED_CONTRIBS | INTERACTION_SHAP;
		if (approximate)
			static_cast<MedXGB *>(original_predictor)->feat_contrib_flags |= APPROX_CONTRIBS;
	}
	try_convert_trees();
}

void TreeExplainer::Learn(MedPredictor *original_pred, const MedFeatures &train_mat) {
	this->original_predictor = original_pred;

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
	int tranform = MODEL_TRANSFORM::identity; //this will explain raw score, the rest are use to explain loss/probabilty or some tranformation, based on model return function

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
	if (retrain_predictor != NULL) {
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
}

int MissingShapExplainer::init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "missing_value")
			missing_value = med_stof(it->second);
		else if (it->first == "max_test")
			max_test = med_stoi(it->second);
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
		else if (it->first == "pp_type") {}
		else
			MTHROW_AND_ERR("Error SHAPExplainer::init - Unknown param \"%s\"\n", it->first.c_str());
	}
	return 0;
}

void MissingShapExplainer::Learn(MedPredictor *original_pred, const MedFeatures &train_mat) {
	this->original_predictor = original_pred;
	retrain_predictor = (MedPredictor *)medial::models::copyInfraModel(original_predictor, false);
	random_device rd;
	mt19937 gen(rd());
	MedMat<float> x_mat;
	train_mat.get_as_matrix(x_mat);
	int nftrs = x_mat.ncols;
	vector<float> labels(train_mat.samples.size()), weights(train_mat.samples.size() + add_new_data, 1);
	vector<int> miss_cnts(train_mat.samples.size() + add_new_data);
	vector<int> missing_hist(nftrs + 1);
	bool verbose_learn = true;

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
		}
		x_mat.add_rows(rows_m);
	}

	// Add data with missing values according to sample masks
	for (size_t i = 0; i < x_mat.nrows; ++i) {
		miss_cnts[i] = msn_count<float>(&x_mat.m[i*nftrs], nftrs, missing_value);
		++missing_hist[miss_cnts[i]];
	}
	for (size_t i = 0; i < x_mat.nrows; ++i) {
		float curr_mask_w = x_mat.nrows / float(missing_hist[miss_cnts[i]]);
		weights[i] = curr_mask_w;
	}
	if (verbose_learn) {
		medial::print::print_hist_vec(miss_cnts, "missing_values hist", "%d");
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
	sample_explain_reasons.resize(matrix.samples.size());
	vector<string> names;
	matrix.get_feature_names(names);
	string bias_name = "Prior_Score";
	vector<float> preds_orig(matrix.samples.size());
	if (matrix.samples.front().prediction.empty()) {
		MedMat<float> mat_x;
		matrix.get_as_matrix(mat_x);
		if (retrain_predictor->transpose_for_predict != (mat_x.transposed_flag > 0))
			mat_x.transpose();
		retrain_predictor->predict(mat_x, preds_orig);
	}
	else {
		for (size_t i = 0; i < preds_orig.size(); ++i)
			preds_orig[i] = matrix.samples[i].prediction[0];
	}

	MedTimer tm;
	tm.start();
	chrono::high_resolution_clock::time_point tm_prog = chrono::high_resolution_clock::now();
	int progress = 0;
	int max_loop = (int)matrix.samples.size();

#pragma omp parallel for
	for (int i = 0; i < matrix.samples.size(); ++i)
	{
		vector<float> features_coeff;
		float pred_shap = 0;
		medial::shapley::explain_shapley(matrix, (int)i, max_test, retrain_predictor, missing_value, features_coeff,
			sample_masks_with_repeats, select_from_all, uniform_rand, use_shuffle, false);

		for (size_t j = 0; j < names.size(); ++j)
			pred_shap += features_coeff[j];

#pragma omp critical 
		{
			map<string, float> &curr_res = sample_explain_reasons[i];
			for (size_t j = 0; j < names.size(); ++j)
				curr_res[names[j]] = features_coeff[j];
			//Add prior to score:
			curr_res[bias_name] = preds_orig[i] - pred_shap; //that will sum to current score
		}

#pragma omp atomic
		++progress;
		double duration = (unsigned long long)(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
			- tm_prog).count()) / 1000000.0;
		if (duration > 15 && progress % 50 == 0) {
#pragma omp critical
			tm_prog = chrono::high_resolution_clock::now();
			double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
				- tm.t[0]).count()) / 1000000.0;
			double estimate_time = int(double(max_loop - progress) / double(progress) * double(time_elapsed));
			MLOG("SHAPLEY Processed %d out of %d(%2.2f%%) time elapsed: %2.1f Minutes, "
				"estimate time to finish %2.1f Minutes\n", progress, max_loop, 100.0*(progress / float(max_loop)), time_elapsed / 60,
				estimate_time / 60.0);
		}
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
	else
		MTHROW_AND_ERR("Unknown type %s\n", type.c_str());
}

int ShapleyExplainer::init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "gen_type")
			gen_type = GeneratorType_fromStr(it->second);
		else if (it->first == "generator_args")
			generator_args = it->second;
		else if (it->first == "missing_value")
			missing_value = med_stof(it->second);
		else if (it->first == "max_test")
			max_test = med_stoi(it->second);
		else if (it->first == "sampling_args")
			sampling_args = it->second;
		else if (it->first == "grouping")
			grouping = it->second;
		else if (it->first == "pp_type") {}
		else
			MTHROW_AND_ERR("Error in ShapleyExplainer::init - Unsupported param \"%s\"\n", it->first.c_str());
	}

	init_sampler(); //from args

	return 0;
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
	default:
		MTHROW_AND_ERR("Error in ShapleyExplainer::init_sampler() - Unsupported Type %d\n", gen_type);
	}
}

void ShapleyExplainer::Learn(MedPredictor *original_pred, const MedFeatures &train_mat) {
	this->original_predictor = original_pred;
	_sampler->learn(train_mat.data);

	if (!grouping.empty())
		read_feature_grouping(grouping, train_mat, group2Ind, groupNames);
	else {
		int icol = 0;
		for (auto& rec : train_mat.data) {
			group2Ind.push_back({ icol++ });
			groupNames.push_back(rec.first);
		}
	}
}

void ShapleyExplainer::explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const {
	sample_explain_reasons.resize(matrix.samples.size());
	vector<string> names;
	matrix.get_feature_names(names);
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

	MedProgress progress("ShapleyExplainer", (int)matrix.samples.size(), 15);
#pragma omp parallel for
	for (int i = 0; i < matrix.samples.size(); ++i)
	{
		vector<float> features_coeff;
		float pred_shap = 0;
		medial::shapley::explain_shapley(matrix, (int)i, max_test, original_predictor, missing_value, *_sampler.get(), 1,
			sampler_sampling_args, features_coeff, false);

		for (size_t j = 0; j < names.size(); ++j)
			pred_shap += features_coeff[j];

#pragma omp critical 
		{
			map<string, float> &curr_res = sample_explain_reasons[i];
			for (size_t j = 0; j < names.size(); ++j)
				curr_res[names[j]] = features_coeff[j];
			//Add prior to score:
			curr_res[bias_name] = preds_orig[i] - pred_shap; //that will sum to current score
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


int LimeExplainer::init(map<string, string> &mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "gen_type")
			gen_type = GeneratorType_fromStr(it->second);
		else if (it->first == "generator_args")
			generator_args = it->second;
		else if (it->first == "missing_value")
			missing_value = med_stof(it->second);
		else if (it->first == "max_test")
			max_test = med_stoi(it->second);
		else if (it->first == "sampling_args")
			sampling_args = it->second;
		else if (it->first == "p_mask")
			p_mask = med_stof(it->second);
		else if (it->first == "n_masks")
			n_masks = med_stoi(it->second);
		else if (it->first == "grouping")
			grouping = it->second;
		else if (it->first == "pp_type") {}
		else
			MTHROW_AND_ERR("Error in LimeExplainer::init - Unsupported param \"%s\"\n", it->first.c_str());
	}
	init_sampler(); //from args

	return 0;
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

void LimeExplainer::Learn(MedPredictor *original_pred, const MedFeatures &train_mat) {
	this->original_predictor = original_pred;
	if (!grouping.empty())
		read_feature_grouping(grouping, train_mat, group2Ind, groupNames);
	else {
		int icol = 0;
		for (auto& rec : train_mat.data) {
			group2Ind.push_back({ icol++ });
			groupNames.push_back(rec.first);
		}
	}

	_sampler->learn(train_mat.data);
}

void LimeExplainer::explain(const MedFeatures &matrix, vector<map<string, float>> &sample_explain_reasons) const {
	vector<vector<float>> alphas;

	medial::shapley::get_shapley_lime_params(matrix, original_predictor, _sampler.get(), p_mask, n_masks, missing_value,
		sampler_sampling_args, group2Ind, groupNames, alphas);

	sample_explain_reasons.resize(matrix.samples.size());

	for (size_t i = 0; i < sample_explain_reasons.size(); ++i)
	{
		map<string, float> &curr = sample_explain_reasons[i];
		const vector<float> &curr_res = alphas[i];
		for (size_t k = 0; k < groupNames.size(); ++k)
			curr[groupNames[k]] = curr_res[k];
	}
}