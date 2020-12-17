#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_FEATCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "FeatureProcess.h"
#include "MedStat/MedStat/MedBootstrap.h"
#include "MedStat/MedStat/bootstrap.h"
#include "ExplainWrapper.h"
#include <fstream>

// Selection top to bottom.
void IterativeFeatureSelector::doTop2BottomSelection(MedFeatures& features, map<string, vector<string>>& featureFamilies, MedBootstrapResult& bootstrapper) {
	ofstream f_progress;

	if (!progress_file_path.empty()) {
		f_progress.open(progress_file_path);
		if (!f_progress.good())
			MTHROW_AND_ERR("Error IterativeFeatureSelector::doTop2BottomSelection - can't open progress file\n");
	}

	if (verbose)
		MLOG("Running top - to - bottom feature selection on %d folds out of %d\n", folds.size(), nfolds);

	features.weights.resize(features.samples.size(), 1.0);
	int nRuns = (int)folds.size();

	// Prepare sets of families
	set<string> selectedFamilies;
	set<string> unSelectedFamilies;
	for (auto& familiy : featureFamilies)
		selectedFamilies.insert(familiy.first);

	// Build train+test list of indices + labels;
	vector<vector<int>> trainRows(nRuns), testRows(nRuns);
	vector<vector<float>> trainLabels(nRuns);
	vector<vector<MedSample>> testSamples(nRuns);
	vector<int> allTestRows;
	MedFeatures bootstrapFeatures;

	prepare_for_iterations(bootstrapper, features, folds, trainRows, testRows, trainLabels, testSamples, bootstrapFeatures);

	//Commit rew, ignore on slecetedFamilies
	//skip if all familiy is in required or ignored.
	set<string> selectedFamilies_filtered;
	for (const string &family : selectedFamilies) {
		//check candidate for filters:
		bool all_in_req = true, all_in_ignore = true;
		const vector<string> &family_feats = featureFamilies.at(family);
		for (int i = 0; i < family_feats.size() && all_in_req; ++i)
			all_in_req = resolved_required.find(family_feats[i]) != resolved_required.end();
		for (int i = 0; i < family_feats.size() && all_in_ignore; ++i)
			all_in_ignore = resolved_ignored.find(family_feats[i]) != resolved_ignored.end();
		if (!all_in_req && !all_in_ignore)
			selectedFamilies_filtered.insert(family);
		else {
			if (all_in_req && all_in_ignore)
				MWARN("WARN!!! Skip family %s - both in required and in ignored\n", family.c_str());
			else {
				if (all_in_req)
					MLOG("Skip family %s (in required)\n", family.c_str());
				if (all_in_ignore)
					MLOG("Skip family %s (in ignored)\n", family.c_str());
			}
		}
	}
	selectedFamilies = move(selectedFamilies_filtered);

	// Iterative addition of families
	MedTimer timer;
	while (selectedFamilies.size() > 0) {
		// Loop on families
		vector<pair<string, float> > scores;
		if (verbose)
			MLOG("Checking %d families to remove\n", selectedFamilies.size());

		set<string> families = selectedFamilies;
		timer.start();
		int counter = 0;
		for (const string &family : families) {
			selectedFamilies.erase(family);
			++counter;

			// names
			vector<string> selectedFeatures;
			//add required
			for (const string &ftr : resolved_required)
				selectedFeatures.push_back(ftr);
			//add current families (expect family - iteration loop)
			for (const string &addFamily : selectedFamilies)
				selectedFeatures.insert(selectedFeatures.end(), featureFamilies[addFamily].begin(), featureFamilies[addFamily].end());

			// Set parameters value:
			string predictorParamsString = predictor_params_vec.empty() ? predictor_params : predictor_params_vec[selectedFeatures.size()];

			// Cross Validation
			if (selectedFeatures.size() > 0) {
				int samplesIdx = 0;
				for (int iRun = 0; iRun < nRuns; iRun++) {

					// Get Train SubMatrix
					MedMat<float> matrix;
					features.get_as_matrix(matrix, selectedFeatures, trainRows[iRun]);

					// Initialize predictor
					MedPredictor *prd = MedPredictor::make_predictor(predictor, predictorParamsString);

					// Learn
					//	matrix.write_to_csv_file("tempMat");
					prd->learn(matrix, trainLabels[iRun]);

					// Get Test SubMatrix
					matrix.clear();
					features.get_as_matrix(matrix, selectedFeatures, testRows[iRun]);

					// Predict
					vector<float> preds;
					prd->predict(matrix, preds);
					delete prd;

					// Insert predictions to samples
					for (size_t i = 0; i < preds.size(); i++)
						bootstrapFeatures.samples[samplesIdx + i].prediction[0] = preds[i];
					samplesIdx += (int)preds.size();

				}

				// Performance
				bootstrapper.bootstrap(bootstrapFeatures);
				scores.push_back({ family,  bootstrapper.bootstrap_results["bs"][measurement_name] });
			}
			else
				scores.push_back({ family,(float)-1.0 });

			selectedFamilies.insert(family);

			string report_s = "checked backward : signal/family " + scores.back().first + " score : " + measurement_name + " = " + to_string(scores.back().second) + "\n";
			report.push_back(report_s);
			if (verbose) {
				medial::print::log_with_file(f_progress, "\tChecked removing family %s with %s = %f\n", scores.back().first.c_str(), measurement_name.c_str(), scores.back().second);

				timer.take_curr_time();
				double diff = timer.diff_sec() / 60.0;
				MLOG("\tCurrent round running for %2.2f min. Estimated time : %2.2f min.\n", diff, diff * ((float)families.size()) / counter);
			}
			else
				f_progress << "\tChecked removing family " << scores.back().first << " with " << measurement_name << " = " << scores.back().second << "\n";
		}

		// Top Performing families
		sort(scores.begin(), scores.end(), [](const pair<string, float>& left, const pair<string, float>& right) { return left.second > right.second; });
		for (int i = 0; i < rates_vec[selectedFamilies.size() + 1]; i++) {
			string report_s = "Removing family " + scores[i].first + " with " + measurement_name + " = " + to_string(scores[i].second);
			report.push_back(report_s);
			if (verbose)
				medial::print::log_with_file(f_progress, "REPORT: %s\n", report_s.c_str());
			else
				f_progress << report_s << "\n";
			unSelectedFamilies.insert(scores[i].first);
			selectedFamilies.erase(scores[i].first);
		}

		// Should we stop ?
		selected.clear();
		for (string ftr : resolved_required)
			selected.push_back(ftr);
		for (string family : selectedFamilies)
			selected.insert(selected.end(), featureFamilies[family].begin(), featureFamilies[family].end());
		if (selected.size() <= numToSelect)
			return;
	}

	if (!progress_file_path.empty())
		f_progress.close();
}

void IterativeFeatureSelector::doBottom2TopSelection(MedFeatures& features, map<string, vector<string>>& featureFamilies, MedBootstrapResult& bootstrapper) {
	ofstream f_progress;

	if (!progress_file_path.empty()) {
		f_progress.open(progress_file_path);
		if (!f_progress.good())
			MTHROW_AND_ERR("Error IterativeFeatureSelector::doBottom2TopSelection - can't open progress file\n");
	}

	if (verbose)
		MLOG("Running bottom - to - top feature selection on %d folds out of %d\n", folds.size(), nfolds);

	features.weights.resize(features.samples.size(), 1.0);

	int nRuns = (int)folds.size();

	// Prepare sets of families
	set<string> selectedFamilies;
	set<string> unSelectedFamilies;
	for (auto& familiy : featureFamilies)
		unSelectedFamilies.insert(familiy.first);
	int nFamilies = (int)unSelectedFamilies.size();

	// Build train+test list of indices + labels;
	vector<vector<int>> trainRows(nRuns), testRows(nRuns);
	vector<vector<float>> trainLabels(nRuns);
	vector<vector<MedSample>> testSamples(nRuns);

	MedFeatures bootstrapFeatures;
	prepare_for_iterations(bootstrapper, features, folds, trainRows, testRows, trainLabels, testSamples, bootstrapFeatures);

	// Iterative addition of families
	MedTimer timer;

	while (selectedFamilies.size() < featureFamilies.size()) {
		// Loop on families
		vector<pair<string, float> > scores;
		if (verbose)
			MLOG("Checking %d families to add\n", unSelectedFamilies.size());

		timer.start();
		int counter = 0;
		for (string family : unSelectedFamilies) {
			selectedFamilies.insert(family);
			counter++;

			// names
			vector<string> selectedFeatures;
			for (string ftr : resolved_required)
				selectedFeatures.push_back(ftr);

			for (string addFamily : selectedFamilies)
				selectedFeatures.insert(selectedFeatures.end(), featureFamilies[addFamily].begin(), featureFamilies[addFamily].end());

			// Set parameters value:
			string predictorParamsString = predictor_params_vec.empty() ? predictor_params : predictor_params_vec[selectedFeatures.size()];

			// Cross Validation
			int samplesIdx = 0;
			for (int iRun = 0; iRun < nRuns; iRun++) {
				// Get Train SubMatrix
				MedMat<float> matrix;
				features.get_as_matrix(matrix, selectedFeatures, trainRows[iRun]);

				// Initialize predictor
				MedPredictor *prd = MedPredictor::make_predictor(predictor, predictorParamsString);

				// Learn
				//	matrix.write_to_csv_file("tempMat");
				prd->learn(matrix, trainLabels[iRun]);

				// Get Test SubMatrix
				matrix.clear();
				features.get_as_matrix(matrix, selectedFeatures, testRows[iRun]);

				// Predict
				vector<float> preds;
				prd->predict(matrix, preds);
				delete prd;

				// Insert predictions to samples
				for (size_t i = 0; i < preds.size(); i++)
					bootstrapFeatures.samples[samplesIdx + i].prediction[0] = preds[i];
				samplesIdx += (int)preds.size();
			}

			// Performance
			bootstrapper.bootstrap(bootstrapFeatures);
			scores.push_back({ family,  bootstrapper.bootstrap_results["bs"][measurement_name] });

			selectedFamilies.erase(family);
			string report_s = "checked forward : signal/family " + scores.back().first + " score : " + measurement_name + " = " + to_string(scores.back().second) + "\n";
			report.push_back(report_s);

			if (verbose) {
				medial::print::log_with_file(f_progress, "\tChecked adding family %s with %s = %f\n", scores.back().first.c_str(), measurement_name.c_str(), scores.back().second);
				timer.take_curr_time();
				double diff = timer.diff_sec() / 60.0;
				MLOG("\tCurrent round: Adding to %d out of %d. Running for %2.2f min. Estimated time : %2.2f min.\n", selectedFamilies.size(), nFamilies, diff,
					diff * ((float)unSelectedFamilies.size()) / counter);
			}
			else
				f_progress << "\tChecked adding family " << scores.back().first << " with " << measurement_name << " = "
				<< scores.back().second << "\n";
		}

		// Top Performing families
		sort(scores.begin(), scores.end(), [](const pair<string, float>& left, const pair<string, float>& right) { return left.second > right.second; });
		for (int i = 0; i < rates_vec[selectedFamilies.size() + 1]; i++) {
			string report_s = "Adding family " + scores[i].first + " with " + measurement_name + " = " + to_string(scores[i].second);
			report.push_back(report_s);
			if (verbose)
				medial::print::log_with_file(f_progress, "REPORT: %s\n", report_s.c_str());
			else
				f_progress << report_s << "\n";
			selectedFamilies.insert(scores[i].first);
			unSelectedFamilies.erase(scores[i].first);
			if (unSelectedFamilies.empty())
				break;
		}

		// Should we stop ?
		selected.clear();
		for (string ftr : resolved_required)
			selected.push_back(ftr);
		for (string family : selectedFamilies)
			selected.insert(selected.end(), featureFamilies[family].begin(), featureFamilies[family].end());
		if (selected.size() <= numToSelect)
			return;
	}

	if (!progress_file_path.empty())
		f_progress.close();
}

// Selection top to bottom.
void IterativeFeatureSelector::retraceTop2BottomSelection(MedFeatures& features, map<string, vector<string>>& featureFamilies, MedBootstrapResult& bootstrapper, vector<string>& order, int start, int end) {

	if (verbose)
		MLOG("Running top - to - bottom tracing on %d folds out of %d\n", folds.size(), nfolds);

	features.weights.resize(features.samples.size(), 1.0);
	int nRuns = (int)folds.size();

	// Prepare sets of families
	set<string> selectedFamilies;
	set<string> unSelectedFamilies;
	for (auto& familiy : featureFamilies)
		selectedFamilies.insert(familiy.first);

	// Build train+test list of indices + labels;
	vector<vector<int>> trainRows(nRuns), testRows(nRuns);
	vector<vector<float>> trainLabels(nRuns);
	vector<vector<MedSample>> testSamples(nRuns);
	vector<int> allTestRows;
	MedFeatures bootstrapFeatures;

	prepare_for_iterations(bootstrapper, features, folds, trainRows, testRows, trainLabels, testSamples, bootstrapFeatures);

	// Remove prior to start
	for (int i = 0; i < start; i++)
		selectedFamilies.erase(order[i]);

	// Loop till end
	for (int i = start; i <= end; i++) {
		// Remove
		selectedFamilies.erase(order[i]);

		// names
		vector<string> selectedFeatures;
		for (string ftr : resolved_required)
			selectedFeatures.push_back(ftr);

		for (string addFamily : selectedFamilies)
			selectedFeatures.insert(selectedFeatures.end(), featureFamilies[addFamily].begin(), featureFamilies[addFamily].end());

		// Set parameters value:
		string predictorParamsString = predictor_params_vec.empty() ? predictor_params : predictor_params_vec[selectedFeatures.size()];

		// Cross Validation
		float score = -1.0;
		if (selectedFeatures.size() > 0) {
			int samplesIdx = 0;
			for (int iRun = 0; iRun < nRuns; iRun++) {

				// Get Train SubMatrix
				MedMat<float> matrix;
				features.get_as_matrix(matrix, selectedFeatures, trainRows[iRun]);

				// Initialize predictor
				MedPredictor *prd = MedPredictor::make_predictor(predictor, predictorParamsString);

				// Learn
				prd->learn(matrix, trainLabels[iRun]);

				// Get Test SubMatrix
				matrix.clear();
				features.get_as_matrix(matrix, selectedFeatures, testRows[iRun]);

				// Predict
				vector<float> preds;
				prd->predict(matrix, preds);
				delete prd;

				// Insert predictions to samples
				for (size_t i = 0; i < preds.size(); i++)
					bootstrapFeatures.samples[samplesIdx + i].prediction[0] = preds[i];
				samplesIdx += (int)preds.size();

			}

			// Performance
			bootstrapper.bootstrap(bootstrapFeatures);
			score = bootstrapper.bootstrap_results["bs"][measurement_name];
		}

		string report_s = "Removing family " + order[i] + " with " + measurement_name + " = " + to_string(score);
		report.push_back(report_s);
		if (verbose)
			MLOG("REPORT: %s\n", report_s.c_str());
	}
}

void IterativeFeatureSelector::retraceBottom2TopSelection(MedFeatures& features, map<string, vector<string>>& featureFamilies, MedBootstrapResult& bootstrapper, vector<string>& order, int start, int end) {

	if (verbose)
		MLOG("Running bottom - to - top tracing on %d folds out of %d\n", folds.size(), nfolds);

	features.weights.resize(features.samples.size(), 1.0);

	int nRuns = (int)folds.size();

	// Prepare sets of families
	set<string> selectedFamilies;

	// Build train+test list of indices + labels;
	vector<vector<int>> trainRows(nRuns), testRows(nRuns);
	vector<vector<float>> trainLabels(nRuns);
	vector<vector<MedSample>> testSamples(nRuns);

	MedFeatures bootstrapFeatures;
	prepare_for_iterations(bootstrapper, features, folds, trainRows, testRows, trainLabels, testSamples, bootstrapFeatures);

	// Add  prior to start
	for (int i = 0; i < start; i++)
		selectedFamilies.insert(order[i]);

	// From start to end
	for (int i = start; i <= end; i++) {
		selectedFamilies.insert(order[i]);

		// names
		vector<string> selectedFeatures;
		for (string ftr : resolved_required)
			selectedFeatures.push_back(ftr);

		for (string addFamily : selectedFamilies)
			selectedFeatures.insert(selectedFeatures.end(), featureFamilies[addFamily].begin(), featureFamilies[addFamily].end());

		// Set parameters value:
		string predictorParamsString = predictor_params_vec.empty() ? predictor_params : predictor_params_vec[selectedFeatures.size()];

		// Cross Validation
		int samplesIdx = 0;
		float score;
		for (int iRun = 0; iRun < nRuns; iRun++) {
			// Get Train SubMatrix
			MedMat<float> matrix;
			features.get_as_matrix(matrix, selectedFeatures, trainRows[iRun]);

			// Initialize predictor
			MedPredictor *prd = MedPredictor::make_predictor(predictor, predictorParamsString);

			// Learn
			//	matrix.write_to_csv_file("tempMat");
			prd->learn(matrix, trainLabels[iRun]);

			// Get Test SubMatrix
			matrix.clear();
			features.get_as_matrix(matrix, selectedFeatures, testRows[iRun]);

			// Predict
			vector<float> preds;
			prd->predict(matrix, preds);
			delete prd;

			// Insert predictions to samples
			for (size_t i = 0; i < preds.size(); i++)
				bootstrapFeatures.samples[samplesIdx + i].prediction[0] = preds[i];
			samplesIdx += (int)preds.size();
		}

		// Performance
		bootstrapper.bootstrap(bootstrapFeatures);
		score = bootstrapper.bootstrap_results["bs"][measurement_name];

		string report_s = "Adding family " + order[i] + " with " + measurement_name + " = " + to_string(score);
		report.push_back(report_s);
		if (verbose)
			MLOG("REPORT: %s\n", report_s.c_str());
	}
}

// Utilities
// Preparation
void IterativeFeatureSelector::prepare_for_iterations(MedBootstrapResult& bootstrapper, MedFeatures& features, vector<int>& folds, vector<vector<int>>& trainRows, vector<vector<int>>& testRows,
	vector<vector<float>>&trainLabels, vector<vector<MedSample>>&testSamples, MedFeatures& bootstrapFeatures) {

	int nRuns = (int)folds.size();
	int nSamples = (int)features.samples.size();

	// Select necessary features
	vector<string> required_features;
	for (Filter_Param& convert_param : bootstrapper.bootstrap_params.filter_cohort["bs"]) {
		if (convert_param.param_name != "Time-Window" && convert_param.param_name != "Label")
			required_features.push_back(resolve_feature_name(features, convert_param.param_name));
	}

	// Attributes
	bootstrapFeatures.time_unit = features.time_unit;
	for (string name : required_features) {
		bootstrapFeatures.attributes[name] = features.attributes[name];
		bootstrapFeatures.attributes[name].normalized = false;
	}

	// Denormalize
	map<string, vector<float>> denormalized;
	for (string name : required_features) {
		denormalized[name].resize(nSamples);
		for (int i = 0; i < nSamples; i++)
			denormalized[name][i] = features.data[name][i] * bootstrapFeatures.attributes[name].denorm_sdv + bootstrapFeatures.attributes[name].denorm_mean;
	}

	// Build train+test lists
	for (int iRun = 0; iRun < nRuns; iRun++) {
		int fold = folds[iRun];
		for (int i = 0; i < features.samples.size(); i++) {
			MedSample& sample = features.samples[i];
			if (sample.split == fold) {
				testRows[iRun].push_back(i);
				testSamples[iRun].push_back(sample);
				bootstrapFeatures.samples.push_back(sample);
				bootstrapFeatures.samples.back().prediction.assign(1, 0.0);
				for (string name : required_features)
					bootstrapFeatures.data[name].push_back(denormalized[name][i]);
			}
			else {
				trainRows[iRun].push_back(i);
				trainLabels[iRun].push_back(sample.outcome);
			}
		}
	}
}

// Get rates vector
void IterativeFeatureSelector::get_rates_vec() {

	rates_vec = { 1 };

	vector<string> fields, subFields;
	boost::split(fields, rates, boost::is_any_of(","));
	for (string& entry : fields) {
		boost::split(subFields, entry, boost::is_any_of(":"));
		if (subFields.size() != 2)
			MTHROW_AND_ERR("Cannot parse entry %s in rates\n", entry.c_str());
		int max = stoi(subFields[0]);
		int step = stoi(subFields[1]);
		if (max < rates_vec.size())
			MTHROW_AND_ERR("Problems parsing rate \'%s\'\n", rates.c_str());

		rates_vec.resize(max + 1, step);
	}
}

// Read parameters file
void IterativeFeatureSelector::read_params_vec()
{
	ifstream inf(predictor_params_file);

	MLOG("iterative_selector: reading %s\n", predictor_params_file.c_str());
	if (!inf)
		MTHROW_AND_ERR("iterative_selector: can't open file %s for read\n", predictor_params_file.c_str());

	string curr_line;
	vector<string> fields;

	// Read
	// Format: interval_min_ftrs interval_max_ftrs paramters_string 
	while (getline(inf, curr_line)) {
		if (curr_line[0] != '#') {
			boost::split(fields, curr_line, boost::is_any_of("\t"));
			if (fields.size() != 3)
				MTHROW_AND_ERR("iterative_seletor: can't parse line \'%s\' in %s\n", curr_line.c_str(), predictor_params_file.c_str());

			int from = stoi(fields[0]);
			int to = stoi(fields[1]);
			predictor_params_vec.resize(to + 1);
			for (int i = from; i <= to; i++)
				predictor_params_vec[i] = fields[2];
		}
	}
	inf.close();

	// Check
	for (size_t i = 0; i < predictor_params_vec.size(); i++) {
		if (predictor_params_vec[i].empty())
			MTHROW_AND_ERR("nFeatures-dependent predictor-params given, but missing for nFeatures=%zd. Currently, this is not allowed\n", i);
	}
}

// Get Families of signals
void IterativeFeatureSelector::get_features_families(MedFeatures& features, map<string, vector<string> >& featureFamilies) {



	vector<string> names;
	features.get_feature_names(names);
	if (!work_on_sets) {
		for (string name : names) {
			if (resolved_required.find(name) == resolved_required.end() && resolved_ignored.find(name) == resolved_ignored.end())
				featureFamilies[name].push_back(name);
		}
	}
	else { // Create sets

		vector<vector<int>> group_inds;
		vector<string> group_names;
		ExplainProcessings::read_feature_grouping(grouping_mode, names, group_inds, group_names);
		for (size_t i = 0; i < group_names.size(); ++i)
		{
			vector<string> &grp_childs = featureFamilies[group_names[i]];
			const vector<int> &curr_grp_idx = group_inds[i];
			grp_childs.resize(curr_grp_idx.size());
			for (size_t j = 0; j < curr_grp_idx.size(); ++j)
				grp_childs[j] = names[curr_grp_idx[j]];
		}
	}

	MLOG("Found %d sets for %d featurs\n", (int)featureFamilies.size(), (int)names.size());

}

// Bootstrapper initialization
// Cohort initialization. Example: Age:30,80/TimeWindow:0,360
void IterativeFeatureSelector::init_bootstrap_cohort(MedBootstrapResult& bootstrapper, string& init) {

	bootstrapper.bootstrap_params.filter_cohort.clear(); //clears "All" default cohort with no filters
	vector<Filter_Param> convert_params;
	if (!init.empty()) { //If empty - no filters, take all samples
		vector<string> params;
		boost::split(params, init, boost::is_any_of("/"));
		convert_params.resize(params.size());
		for (size_t i = 0; i < params.size(); ++i) {
			if (!params[i].empty())
				convert_params[i] = Filter_Param(params[i]);
		}
	}
	bootstrapper.bootstrap_params.filter_cohort["bs"] = convert_params;
}

// Measurement intialization. Examples:
// AUC
// SENS,FPR,0.1
void IterativeFeatureSelector::init_bootstrap_params(MedBootstrapResult& bootstrapper, string& init) {

	vector<string> params;
	boost::split(params, init, boost::is_any_of(","));

	bootstrapper.bootstrap_params.roc_Params.working_point_FPR.clear();
	bootstrapper.bootstrap_params.roc_Params.working_point_PR.clear();
	bootstrapper.bootstrap_params.roc_Params.working_point_Score.clear();
	bootstrapper.bootstrap_params.roc_Params.working_point_SENS.clear();

	if (params.size() == 1) {
		if (params[0] == "AUC")
			measurement_name = "AUC_Obs";
		else
			measurement_name = params[0];
		if ((measurement_name != "AUC_Obs") && (measurement_name != "AUC_Mean"))
			MTHROW_AND_ERR("Unknown single parameter \'%s\'\n", params[0].c_str())
		else
			bootstrapper.bootstrap_params.measurements_with_params = { pair<MeasurementFunctions, Measurement_Params *>(calc_only_auc, NULL) };
	}

	else if (params.size() == 3) {
		float target = stof(params[2]);
		if (params[1] == "PR" && (params[0] == "SENS" || params[0] == "SPEC"))
			bootstrapper.bootstrap_params.roc_Params.working_point_PR.push_back(100 * target);
		else if (params[1] == "FPR" && params[0] == "SENS")
			bootstrapper.bootstrap_params.roc_Params.working_point_FPR.push_back(100 * target);
		else if (params[1] == "SENS" && params[0] == "SPEC")
			bootstrapper.bootstrap_params.roc_Params.working_point_SENS.push_back(100 * target);
		else if (params[1] == "Score" && (params[0] == "SENS" || params[0] == "SPEC"))
			bootstrapper.bootstrap_params.roc_Params.working_point_Score.push_back(target);

		else
			MTHROW_AND_ERR("Cannot parse bootstrap-measuremnt initialization string \'%s\'\n", init.c_str());

		bootstrapper.bootstrap_params.measurements_with_params = { pair<MeasurementFunctions, Measurement_Params *>(calc_roc_measures_with_inc, &bootstrapper.bootstrap_params.roc_Params) };
		measurement_name = format_working_point(params[0] + "@" + params[1], target) + "_Obs";
	}

}