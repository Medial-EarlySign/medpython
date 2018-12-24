/// \n
/// IterativeImputer\n
///\n
/// A general strong imputer that does the following:\n
/// (1) Runs a simple stratified imputer\n
/// (2) Runs iterations completing values (from the least missing to the max missing) where:\n
///     (a) continuous values are calculated using a regressor\n
///     (b) categorial values (less than some bound) are calculated using a multi category classifier\n
/// (3) Repeats the process several times until it converges or until max_iters is reached.\n
///
#include "IterativeImputer.h"
#include <regex>
#include <unordered_set>
#include <boost/regex.hpp>
#include <MedUtils/MedUtils/MedGlobalRNG.h>
#include <MedStat/MedStat/MedPerformance.h>

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//----------------------------------------------------------------------------------------------------------------------
int feature_info::prep_feats_for_pred(MedFeatures &mfd)
{
	feats_for_pred.clear();

	for (auto &feat : mfd.data) {
		if (feat.first != full_name)
			feats_for_pred.push_back(feat.first);
	}
	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
// assumes fi.data was already initialized
int feature_info::prep_indexes(const vector<int> &external_train_idx, const vector<int> &external_test_idx, float missing_value)
{
	train_idx.clear();
	test_idx.clear();
	pred_idx.clear();

	bool create_train = true;

	for (int i=0; i<data_len; i++)
		if (data[i] == missing_value)
			pred_idx.push_back(i);
	for (int i=0; i<external_train_idx.size(); i++)
		if (data[external_train_idx[i]] != missing_value)
			train_idx.push_back(external_train_idx[i]);
	for (int i=0; i<external_test_idx.size(); i++)
		if (data[external_test_idx[i]] != missing_value)
			test_idx.push_back(external_test_idx[i]);

	return 0;


}


//----------------------------------------------------------------------------------------------------------------------
int IterativeImputerParams::init(map<string, string>& mapper)
{
	features_to_impute.clear();
	for (auto entry : mapper) {
		string field = entry.first;
		//! [IterativeImputerParams::init]
		if (field == "regressor") regressor = entry.second;
		else if (field == "regressor_params") regressor_params = entry.second;
		else if (field == "multi_categ_classifier") multi_categ_classifier = entry.second;
		else if (field == "multi_categ_classifier_params") multi_categ_classifier_params = entry.second;
		else if (field == "add_ncateg_var_name") add_ncateg_var_name = entry.second;
		else if (field == "round1_strata") round1_strata = entry.second;
		else if (field == "do_round1") do_round1 = stoi(entry.second);
		else if (field == "round1_moment") round1_moment = stoi(entry.second);
		else if (field == "categorial_bound") categorial_bound = stoi(entry.second);
		else if (field == "max_iterations") max_iterations = stoi(entry.second);
		else if (field == "p_validation") p_validation = stof(entry.second);
		else if (field == "min_vals_for_training") min_vals_for_training = stoi(entry.second);
		else if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "missing_bound") missing_bound = stof(entry.second);
		else if (field == "round_to_resolution") round_to_resolution = stoi(entry.second);
		else if (field == "verbose") verbose = stoi(entry.second);
		else if (field == "features" || field == "names") {
			vector<string> f;
			boost::split(f, entry.second, boost::is_any_of(",:"));
			features_to_impute.insert(features_to_impute.end(), f.begin(), f.end());
		}
		//! [IterativeImputerParams::init]
	}

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::init_feature_info_update(MedFeatures &mfd, feature_info &fi)
{
	fi.data = &(mfd.data[fi.full_name][0]);
	fi.data_len = (int)mfd.data[fi.full_name].size();
	fi.is_missing.clear();

	// calculating several needed data descriptors
	fi.is_missing.resize(fi.data_len);
	for (int i=0; i<fi.data_len; i++) {
		if (fi.data[i] == params.missing_value) {
			fi.n_missing++;
			fi.is_missing[i] = 1;
		}
		else {
			fi.n_with_values++;
			if (fi.data[i] != 0) fi.n_with_non_zero_values++;
			fi.is_missing[i] = 0;
		}
	}
	//fi.prep_feats_for_pred(mfd);
	fi.prep_indexes({}, {}, params.missing_value);

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::init_feature_info(MedFeatures &mfd, string feat_name)
{
	// init train_idx if not initialized
	if (train_idx.size() == 0) {
		is_train.resize(mfd.samples.size(), 0);
		for (auto &pid_rec : mfd.pid_pos_len) {
			train_ids.insert(pid_rec.first);
			if (rand_1() >= params.p_validation) {
				for (int i=0; i<pid_rec.second.second; i++) {
					is_train[pid_rec.second.first + i] = 1;
					train_idx.push_back(pid_rec.second.first + i);
				}
			}
			else {
				test_ids.insert(pid_rec.first);
				for (int i=0; i<pid_rec.second.second; i++) {
					is_train[pid_rec.second.first + i] = 0;
					test_idx.push_back(pid_rec.second.first + i);
				}
			}
		}

	}

	FeatureProcessor fp;
	feature_info fi;
	// first test name uniqueness 
	fi.name = feat_name;
	fi.full_name = fp.resolve_feature_name(mfd, feat_name);
	//fi.name = regex_replace(fi.full_name, regex("FTR_[:digit]+\\."), "");
	boost::regex re("FTR_[[:digit:]]+\\.");
	fi.name = boost::regex_replace(fi.full_name, re, ".");

	fi.data = &(mfd.data[fi.full_name][0]);
	fi.data_len = (int)mfd.data[fi.full_name].size();

	// calculating several needed data descriptors
	unordered_map<float, int> val_cnt;
	fi.is_missing.resize(fi.data_len);
	for (int i=0; i<fi.data_len; i++) {
		if (fi.data[i] == params.missing_value) {
			fi.n_missing++;
			fi.is_missing[i] = 1;
		}
		else {
			fi.n_with_values++;
			if (fi.data[i] != 0) fi.n_with_non_zero_values++;
			fi.is_missing[i] = 0;
			if (fi.data[i] < fi.min) fi.min = fi.data[i];
			if (fi.data[i] > fi.max) fi.max = fi.data[i];
			if (val_cnt.find(fi.data[i]) == val_cnt.end())
				val_cnt[fi.data[i]] = 1;
			else
				val_cnt[fi.data[i]]++;

			//if (fi.name == "Proteinuria_State.last.win_0_10000") {
			//	if (fi.data[i] > 2) {
			//		MLOG("Found val %f in : i %d pid %d date %d\n", i, fi.data[i], mfd.samples[i].id, mfd.samples[i].time);
			//	}
			//}
		}
	}

	fi.n_diff_vals = (int)val_cnt.size();

	// find resolution
	float best_res = -1;
	int n_best = -1;
	float epsilon = (float)1e-4;
	for (float res = 1000; res >= (float) 0.001; res = res/10) {
		int n_in_res = 0;
		float res_above = res*10;
		for (auto &vc : val_cnt) {
			if (vc.first != 0) { // 0 never helps us to detect resolutions
				//if (params.verbose) MLOG("feature %s : resolution %f : n_in_res %d : val & count : %f %d : tests %f %f %f %f : %f %f\n",
				//	fi.name.c_str(), res, n_in_res, vc.first, vc.second, vc.first/res, (float)((int)(vc.first/res)),
				//	vc.first/res_above, (float)((int)(vc.first/res_above)), abs(vc.first/res - (float)((int)(vc.first/res))), abs(vc.first/res_above - (float)((int)(vc.first/res_above))));
				if ((res <= vc.first) && (abs(vc.first/res - (float)((int)(vc.first/res)))<epsilon) &&
					(abs(vc.first/res_above - (float)((int)(vc.first/res_above)))>epsilon))
					n_in_res += vc.second;
			}
		}
		if (params.verbose > 1) MLOG("feature %s : resolution %f : n_in_res %d \n", fi.name.c_str(), res, n_in_res);
		if (n_in_res >= n_best) {
			n_best = n_in_res;
			best_res = res;
		}
	}

	if ((float)n_best/(float)fi.n_with_non_zero_values < 0.25)
		best_res = 0; // means we will not round at all

	
	if (fi.n_diff_vals < params.categorial_bound)
		fi.is_categorial = 1;
	else
	{
		fi.is_categorial = 0;
	}

	fi.resolution = best_res;


	if (fi.n_missing == 0) fi.predictor_type = 0;
	else if (fi.is_categorial) fi.predictor_type = 2;
	else fi.predictor_type = 1;

	//fi.prep_feats_for_pred(mfd);
	fi.prep_indexes(train_idx, test_idx, params.missing_value);

	if (params.verbose) fi.print();

#pragma omp critical
///	if (fi.n_missing > 0 && fi.n_with_values > 0)
	feats.push_back(fi);

	return 0;

}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::init_internals(MedFeatures &mfd)
{
	// going over all features that were selected and initializing feature info for each
	
	// if no features to impute selected - 
	if (params.features_to_impute.size() == 0) {
		for (auto &elem : mfd.data)
			params.features_to_impute.push_back(elem.first);
	}

	int rc = 0;
	for (auto &feat : params.features_to_impute) {
		init_feature_info(mfd, feat);
	}

	return 0;

}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::round_arr(float *arr, int len, float resolution, float _min, float _max)
{
	if (resolution > 0) {
		for (int i=0; i<len; i++) {
			if (arr[i] != params.missing_value) {
				if (arr[i] > _max) arr[i] = _max;
				if (arr[i] < _min) arr[i] = _min;
				arr[i] = resolution * (float)((int)(arr[i]/resolution + (float)0.5));
			}
		}
	}
	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::round_to_resolution(MedFeatures &mfd)
{
	for (auto &fi : feats) {
		round_arr(fi.data, fi.data_len, fi.resolution, fi.min, fi.max);
/*
		if (fi.resolution > 0) {
			for (int i=0; i<fi.data_len; i++) {
				if (fi.data[i] != params.missing_value) {
					if (fi.data[i] > fi.max) fi.data[i] = fi.max;
					if (fi.data[i] < fi.min) fi.data[i] = fi.min;
					fi.data[i] = fi.resolution * (float)((int)(fi.data[i]/fi.resolution + (float)0.5));
				}
			}
		}
*/
	}
	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
// In some variants there should be a first round of imputations taking place
// We use the basic imputer with a strata given in the parameters
int IterativeImputer::learn_first_round(MedFeatures &mfd)
{
	for (auto &fi : feats) {
		if (fi.n_missing > 0) {
			if (params.verbose) MLOG("IterativeImputer learn_first_round: feature_name %s\n", fi.name.c_str());
			FeatureImputer basic_imputer;
			string init_str = "moment_type=" + to_string(params.round1_moment) + ";strata=" + params.round1_strata;
			basic_imputer.feature_name = fi.name;
			if (params.verbose) MLOG("IterativeImputer learn_first_round: feature %s init %s\n", fi.name.c_str(), init_str.c_str());
			if (basic_imputer.init_from_string(init_str) < 0)
				MTHROW_AND_ERR("Cannot init FeatureImputer  with init string \'%s\'\n", init_str.c_str());

			if (basic_imputer.learn(mfd, train_ids) < 0) {
				MERR("IterativeImputer : Failed 1st round basic imputer Learn on feature %s.... \n", fi.name.c_str());
				return -1;
			}
			if (params.verbose) MLOG("IterativeImputer learn_first_round: passed learn for feature %s\n", fi.name.c_str());

			first_round_imputers.push_back(basic_imputer);
		}
	}

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::apply_first_round(MedFeatures &mfd)
{
	for (auto &basic_imputer : first_round_imputers) {
		if (params.verbose) MLOG("IterativeImputer apply_first_round : feature %s\n", basic_imputer.feature_name.c_str());
		if (basic_imputer.apply(mfd) < 0) {
			MERR("IterativeImputer : Failed 1st round basic imputer Apply on feature %s.... \n", basic_imputer.feature_name.c_str());
			return -1;
		}

	}
	if (params.round_to_resolution && round_to_resolution(mfd) < 0) return -1;

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::Learn(MedFeatures &mfd) 
{
	MedFeatures &my_mfd = mfd;

	MedFeatures learn_mfd;
	if ((params.do_round1 && params.max_iterations) || params.max_iterations>1) {
		learn_mfd = mfd;
		my_mfd = learn_mfd;
	}

	if (params.verbose) MLOG("IterativeImputer Learn() init_internals\n");
	if (init_internals(mfd) < 0) return -1;
	if (params.verbose) MLOG("IterativeImputer Learn() learn_first_round (do %d)\n", params.do_round1);
	if (params.do_round1 && learn_first_round(my_mfd) < 0) return -1;
	if (params.do_round1 && apply_first_round(my_mfd) < 0) return -1;

	for (int iter=0; iter<params.max_iterations; iter++) {
		learn_iteration(mfd, iter);
		if (iter < params.max_iterations - 1) apply_iteration(mfd, iter);
	}
	return 0;
}


//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::Apply(MedFeatures &mfd)
{
	if (params.verbose) MLOG("IterativeImputer Apply() apply_first_round()\n");

	for (auto &fi : feats)
		init_feature_info_update(mfd, fi);

	if (params.do_round1 && apply_first_round(mfd) < 0) return -1;
	if (params.verbose) MLOG("IterativeImputer Apply() round_to_resolution()\n");
	for (int iter=0; iter<params.max_iterations; iter++) {
		apply_iteration(mfd, iter);
	}
	for (auto &fi : feats)
		mfd.attributes[fi.full_name].imputed = true;
	fprintf(stderr, "IterativeImputer::Apply() NOT IMPLEMENTED YET\n"); 
	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::learn_iteration(MedFeatures &mfd, int iter)
{
	vector<MedPredictor *> predictors_vec;
	for (int i=0; i<feats.size(); i++) {
		auto &fi = feats[i];
	//for (auto &fi : feats) {
		MedPredictor *predictor;
		find_feats_to_learn_from(i);
		feats_for_pred_inds_to_names(fi);
		if (params.verbose) MLOG("IterativeImputer::learn_iteration pred_feats %d %d\n", fi.feats_for_pred.size(), feats[i].feats_for_pred.size());
		if (fi.predictor_type == 0) {
			predictor = NULL;
		}
		else {

			// we need to train a regressor using the cases for which we have values
			// in order to do that we create 2 matrices:
			// (1) Train matrix : 1-p_validation of the data with values for training
			// (2) Test matrix : p_validation of the data with values for training
			// (3) Predict matrix : the cases for which we actually need to predict -> will be created and used in the apply stages
			MedMat<float> x_train, y_train, x_test, y_test;
			mfd.get_as_matrix(x_train, fi.feats_for_pred, fi.train_idx);
			x_train.normalized_flag = true;
			mfd.get_as_matrix(y_train, { fi.full_name }, fi.train_idx);
			y_train.normalized_flag = true;
			if (test_idx.size() > 0) {
				mfd.get_as_matrix(x_test, fi.feats_for_pred, fi.test_idx);
				x_test.normalized_flag = true;
				mfd.get_as_matrix(y_test, { fi.full_name }, fi.test_idx);
			}
			// now we are ready to train and test our regressor
			predictor = MedPredictor::make_predictor(params.regressor, params.regressor_params);
			predictor->learn(x_train, y_train);

			if (params.verbose) MLOG("IterativeImputer::learn_iteration :: iter %d :: feature %s :: train %d x %d , test %d x %d :: predictor is %d\n",
				iter, fi.name.c_str(), x_train.nrows, x_train.ncols, x_test.nrows, x_test.ncols, (int)predictor->classifier_type);


			if (x_test.nrows > 0) {
				vector<float> preds;
				predictor->predict(x_test, preds);
				round_arr(&preds[0], (int)preds.size(), fi.resolution, fi.min, fi.max);
				double corr = get_pearson_corr(y_test.m, preds);
				double d2 = get_squared_dist(y_test.m, preds);
				double dabs = get_abs_avg_dist(y_test.m, preds);
				double dabs_rel = get_abs_relative_avg_dist(y_test.m, preds);
				double acc = get_vecs_accuracy(y_test.m, preds, (double)fi.resolution);

				if (params.verbose) MLOG("IterativeImputer::learn_iteration :: iter %d :: feature %s :: corr %f d2 %f dabs %f dabs_rel %f acc %f\n", 
					iter, fi.name.c_str(), corr, d2, dabs, dabs_rel, acc);
			}
		}

		predictors_vec.push_back(predictor);
	}

	predictors.push_back(predictors_vec);

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::apply_iteration(MedFeatures &mfd, int iter)
{
	if (params.verbose) MLOG("IterativeImputer::apply_iteration :: iter %d :: feats size %d :: predictors %d\n", iter, feats.size(), predictors.size());
	MedFeatures interim;
	for (int j=0; j<feats.size(); j++) {
		auto &fi = feats[j];
		feats_for_pred_inds_to_names(fi);
		if (fi.predictor_type == 0) {
			if (params.verbose) MLOG("IterativeImputer::apply_iteration :: iter %d :: feature %s :: NULL predictor - nothing to do\n", iter, fi.name.c_str());
		}
		else {
			MedPredictor *predictor = predictors[iter][j];
			// we need to apply a regressor using the cases for the missing values
			MedMat<float> x_pred;
			mfd.get_as_matrix(x_pred, fi.feats_for_pred, fi.pred_idx);
			x_pred.normalized_flag = true;
			if (params.verbose) MLOG("IterativeImputer::apply_iteration :: iter %d :: feature %s :: pred mat %d x %d :: predictor is %d\n",
				iter, fi.name.c_str(), x_pred.nrows, x_pred.ncols, predictor->classifier_type);
			vector<float> preds;
			interim.data[fi.full_name] = preds;
			predictor->predict(x_pred, interim.data[fi.full_name]);
			//float *fdata = &(mfd.data[fi.full_name][0]);
			//for (int i=0; i<fi.pred_idx.size(); i++)
			//	fdata[fi.pred_idx[i]] = preds[i];
		}

	}

	// copy from interim to mfd
	for (int j=0; j<feats.size(); j++) {
		auto &fi = feats[j];
		if (fi.predictor_type != 0) {
			float *fdata = &(mfd.data[fi.full_name][0]);
			float *fdata_interim = &(interim.data[fi.full_name][0]);
			for (int i=0; i<fi.pred_idx.size(); i++)
				fdata[fi.pred_idx[i]] = fdata_interim[i];
		}
	}

	if (params.round_to_resolution && round_to_resolution(mfd) < 0) return -1;

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
// We have to deal with the problem that when some feature is missing it usually means other features are missing as well
// in a very correlated manner.
// The way to improve that is make sure that we don't train with features that are missing > min_bound (say 75%) of the times
// exactly when our feature is missing.
int IterativeImputer::find_feats_to_learn_from(int f_idx)
{
	int miss_size = (int)feats[f_idx].pred_idx.size();
	feats[f_idx].inds_for_pred.clear();
	for (int i=0; i<feats.size(); i++) {
		int cnt = 0;
		if (i != f_idx) {
			for (auto idx : feats[f_idx].pred_idx)
				if (feats[i].is_missing[idx])
					cnt++;
			float p_miss = (float)cnt/(float)(miss_size+1);
			if (p_miss > params.missing_bound) {
				if (params.verbose > 1) MLOG("IterativeImputer :: feat_to_learn for %s :: feat %s :: miss %f\n",
										      feats[f_idx].name.c_str(), feats[i].name.c_str(), p_miss);
			}
			else {
				feats[f_idx].inds_for_pred.push_back(i);
			}
		}
	}

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------
int IterativeImputer::feats_for_pred_inds_to_names(feature_info &fi)
{
	fi.feats_for_pred.clear();
	for (int i=0; i<fi.inds_for_pred.size(); i++)
		fi.feats_for_pred.push_back(feats[i].full_name);

	return 0;
}


/// check if a set of features is affected by the current processor
//.......................................................................................
bool FeatureIterativeImputer::are_features_affected(unordered_set<string>& out_req_features) {

	// If empty = all features are required
	if (out_req_features.empty())
		return true;

	if (imputer.params.features_to_impute.empty())
		return true;

	// Check intersections
	for (string ftr : imputer.params.features_to_impute) {
		if (out_req_features.find(ftr) != out_req_features.end())
			return true;
	}

	return false;
}

/// update sets of required as input according to set required as output to processor
//.......................................................................................
void FeatureIterativeImputer::update_req_features_vec(unordered_set<string>& out_req_features, unordered_set<string>& in_req_features) {

	if (are_features_affected(out_req_features))
		// If active, than everything before is required
		in_req_features.clear();
	else
		// If not, do nothing
		in_req_features = out_req_features;
}

/// Apply imputing model on subset of ids (TBI)
//.......................................................................................
int FeatureIterativeImputer::_apply(MedFeatures& features, unordered_set<int>& ids) {
	MERR("iterativeImputer on subset of ids is not implemented yet\n"); 
	return -1;
}

/// Learn imputing model on subset of ids (TBI)
//.......................................................................................
int FeatureIterativeImputer::Learn(MedFeatures& features, unordered_set<int>& ids) {
	MERR("iterativeImputer on subset of ids is not implemented yet\n"); 
	return -1; 
}