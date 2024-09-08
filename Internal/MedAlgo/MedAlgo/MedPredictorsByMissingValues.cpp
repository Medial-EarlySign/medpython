#include "MedPredictorsByMissingValues.h"
#include <MedProcessTools/MedProcessTools/Calibration.h>
#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL LOG_DEF_LEVEL

int MedPredictorsByMissingValues::init(map<string, string>& mapper) {
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		calibrate_predictions = "1"; // default - use calibration in predict
		if (it->first == "predictor_type")
			predictor_type = it->second;
		else if (it->first == "predictor_params")
			predictor_params = it->second;
		else if (it->first == "masks_params")
			masks_params = it->second;
		else if (it->first == "masks_tw")
			masks_tw = it->second;
		else if (it->first == "calibrate_predictions")
			calibrate_predictions = it->second;
		else
			MTHROW_AND_ERR("Unsupported argument \"%s\" for MedPredictorsByMissingValues\n", it->first.c_str());
	}
	vector<string> masks1, masks2;
	boost::split(masks1, masks_params, boost::is_any_of("|"));
	boost::split(masks2, masks_tw, boost::is_any_of(","));
	if (masks1.size() != masks2.size()) {
		MTHROW_AND_ERR("Number of masks and time windows for masks is not the same");
	}
	return 0;
}


vector<int> MedPredictorsByMissingValues::get_cols_of_mask(string &full_mask, vector<string> &signals) const {
	// for every signal in the mask, we find all features that should be removed

	vector<int> cols;
	int len_cols = 0;

	vector<string> signals_in_mask;
	boost::split(signals_in_mask, full_mask, boost::is_any_of(","));
	if (signals_in_mask.size() == 0) {
		MTHROW_AND_ERR("Check the masks - one of the masks is empty");
	}
	
	string sig, sig1;
	for (int m = 0; m < signals_in_mask.size(); m++) {
		len_cols = cols.size();
		if (signals_in_mask[m] == "Smoking_Intensity") {
			sig = "Smok_Pack_Years";
			sig1 = "Smoking_Intensity";
			for (int s = 0; s < signals.size(); s++) {
				bool isFound = signals[s].find(sig) != string::npos;
				bool isFound1 = signals[s].find(sig1) != string::npos;
				if (isFound | isFound1) {
					cols.insert(cols.end(), s);
				}
			}
		}
		else {
			sig = '.' + signals_in_mask[m] + '.'; 
			for (int s = 0; s < signals.size(); s++) {
				bool isFound = signals[s].find(sig) != string::npos;
				if (isFound) {
					cols.insert(cols.end(), s);
				}
			}
		}
		if (len_cols == cols.size()) {
			MTHROW_AND_ERR("Mask \"%s\" has no features", full_mask.c_str());
		}
	}
	return cols;
}

vector<int> MedPredictorsByMissingValues::get_cols_of_predictor(int i, vector<vector<int>> &cols_per_mask, vector<int> &all_cols) const {

	// what columns should be removed in predictor number i ?
	// -----------------------------------------------------------------------------
	// for every model we find the relevant masks by binary coding of model number i
	// i=0: full model
	// i=1=pow(2,0): no features included from masks[0]
	// i=2=pow(2,1) no features included from masks[1]
	// i=3=1+2=pow(2,0)+pow(2,1): no features included from masks[0] or masks[1]
	// ...
	// i=6=2+4=pow(2,1)+pow(2,2): no features included from masks[1] or masks[2]
	// -----------------------------------------------------------------------------

	vector<int> cols, cols_to_keep_p;
	for (int j = cols_per_mask.size() - 1; j >= 0; j--) {
		if (i >= pow(2, j)) {
			// mask j is part of model i
			// cout << "mask " << j << " is part of predictor " << i << " with " << cols_per_mask[j].size() << " features" << endl;
			i = i - pow(2, j);
			cols.insert(cols.end(), cols_per_mask[j].begin(), cols_per_mask[j].end());
		}
	}

	// keep only the columns that are not removed
	cols_to_keep_p = all_cols;
	// cout << "predictor " << i << " removing " << cols.size() << " features out of " << all_cols.size() << endl;
	for (int s = 0; s < cols.size(); s++)
		cols_to_keep_p.erase(remove(cols_to_keep_p.begin(), cols_to_keep_p.end(), cols[s]), cols_to_keep_p.end());

	return cols_to_keep_p;
}


int MedPredictorsByMissingValues::learn(MedMat<float> &x, MedMat<float> &y, const vector<float> &wgts) {

	MedMat<float> x_copy;
	vector<int> empty_is_all;
	vector<int> cols;
	vector<vector<int>> cols_per_mask;
	vector<int> all_cols, cols_to_keep_p;

	vector<string> masks;
	boost::split(masks, masks_params, boost::is_any_of("|"));
	if (masks.size() == 0) {
		MTHROW_AND_ERR("The predictor needs at least one mask");
	}

	int p = pow(2, masks.size()); // p - number of predictors to train

	for (int s = 0; s < x.signals.size(); s++)
		all_cols.insert(all_cols.end(), s);

	// preperation - what are the signals of each input mask
	for (int m = 0; m < masks.size(); m++) {
		cols = get_cols_of_mask(masks[m], x.signals);
		cols_per_mask.insert(cols_per_mask.end(), cols);
	}
	
	for (int i = 0; i < p; i++) {
		// what are the features of predictor i
		cols_to_keep_p = get_cols_of_predictor(i, cols_per_mask, all_cols);
	
		x_copy = x;
		x_copy.get_sub_mat(empty_is_all, cols_to_keep_p); // retrun sub mat with all the rows and just cols_to_keep columns
		MLOG("predictor %d has %d features\n", i + 1, cols_to_keep_p.size());

		// train model i
		predictors.push_back(MedPredictor::make_predictor(predictor_type, predictor_params));
		MLOG("Going to train predictor %d out of %d\n", i+1, p); 
		predictors[i]->learn(x_copy, y, wgts);
		//x_copy.write_to_csv_file("/nas1/Work/Users/Eitan/LungWithMask/predictions/before_" + to_string(i) + ".csv"); // write to csv
		
		// calibrate predictor i
		// 1. get predictions, and turn it to vector of MedSample
		// 2. add new calibrator to out list of calibrator, and prepare (learn) it using this vector
		
		// predict
		vector<float> scores;
		predictors[i]->predict(x_copy, scores);

		// prepare samples from the predictions
		vector<MedSample> vec(scores.size());
		for (size_t i = 0; i< vec.size(); ++i) {
			vec[i].id = i + 1;
			vec[i].outcome = x_copy.recordsMetadata[i].label;
			vec[i].prediction = { scores[i] };
		}

		// add and train the new calibrator
		calibrators.push_back(static_cast<Calibrator *>(PostProcessor::make_processor("calibrator", "calibration_type=isotonic_regression")));
		calibrators[i]->Learn(vec);

		// the rest was used for checking the code
		//MedSamples samples;
		//samples.import_from_sample_vec(vec);
		//samples.write_to_file("/nas1/Work/Users/Eitan/LungWithMask/predictions/before_" + to_string(i) + ".csv"); // write to csv
		//calibrators[i]->Apply(vec);
		//samples.import_from_sample_vec(vec);
		//samples.write_to_file("/nas1/Work/Users/Eitan/LungWithMask/predictions/after_" + to_string(i) + ".csv"); // write to csv

	}
	return 0;
}

int MedPredictorsByMissingValues::predict(MedFeatures& features) const {
	cout << "+++++++++++ enter predict +++++++++++++" << endl;
	vector<string> masks, masksTime, mask;
	boost::split(masks, masks_params, boost::is_any_of("|"));
	boost::split(masksTime, masks_tw, boost::is_any_of(","));

	vector<int> cols, all_cols, cols_to_keep_p;
	vector<vector<int>> cols_per_mask, cols_per_predictor;
	
	string sig;
	vector<vector<string>> indicator_cols;
	int s_found;
	vector<int> sample_number;
	vector<float> score;

	MedMat<float> x, x_copy;

	features.get_as_matrix(x);
	//cout << "write features matrix to /opt/medial_sign/Work/Eitan/AAA/output/predictions/features.csv" << endl;
	//x.write_to_csv_file("/opt/medial_sign/Work/Eitan/AAA/output/predictions/features.csv"); // write to csv

	//MedMat<unsigned char> x1;
	//cout << "going to prepare masks matrix" << endl;
	//features.get_masks_as_mat(x1);
	//cout << "write masks matrix to /opt/medial_sign/Work/Eitan/AAA/output/predictions/masks.csv" << endl;
	//x1.write_to_csv_file("/opt/medial_sign/Work/Eitan/AAA/output/predictions/masks.csv"); // write to csv

	for (int s = 0; s < x.signals.size(); s++)
		all_cols.insert(all_cols.end(), s);

	// prepare predictors
	int p = pow(2, masks.size()); // p - number of trained models
	for (int i = 0; i < p; i++) {
		// cout << "going to prepare predictor " << i + 1 << endl;
		predictors[i]->prepare_predict_single();
	}

	// preperation - what are the signals of each input mask
	for (int m = 0; m < masks.size(); m++) {
		cols = get_cols_of_mask(masks[m], x.signals);
		cols_per_mask.insert(cols_per_mask.end(), cols);
	}

	// preperation - what are the features of predictor i
	for (int i = 0; i < p; i++) {
		// what are the features of predictor i
		cols_to_keep_p = get_cols_of_predictor(i, cols_per_mask, all_cols);
		cols_per_predictor.insert(cols_per_predictor.end(), cols_to_keep_p);
	}

	// for every mask - find the indicator columns
	for (int i = 0; i < masks.size(); i++) {
		boost::split(mask, masks[i], boost::is_any_of(","));
		for (int m = 0; m < mask.size(); m++) {
			if (mask[m] == "Smoking_Intensity")
				sig = "Smoking_Intensity";
			else
				sig = "." + mask[m] + ".last.win_0_" + masksTime[i];
			s_found = 0;
			for (int s = 0; s < x.signals.size(); s++) {
				bool isFound = x.signals[s].find(sig) != string::npos;
				if (isFound) {
					mask[m] = x.signals[s];
					s = x.signals.size() + 1;
					s_found = 1;
				}
			}
			if (s_found == 0) {
				MTHROW_AND_ERR("Could not find feature \"%s\" to check if imputed", sig.c_str());
			}
		}
		indicator_cols.insert(indicator_cols.end(), mask);
	}
	vector<float> sample(x.signals.size());
	vector<float> sample4p(x.signals.size());

	vector<int> mask_stat(masks.size());
	vector<int> predictor_stat(p);

	// for every sample - find the right predictor
	cout << "+++++++++++ enter predict per sample ++++" << endl;
//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < x.get_nrows(); i++) {
		//cout << "++++++++ row " << i << endl;
		int p1 = 0; // default, no mask is needed
		// for every mask, check if value is missing in 50% or more of the indicator columns
		for (int s = 0; s < masks.size(); s++) {
			//cout << "++++++ mask " << s+1 << " out of " << masks.size() << endl;
			int missing = 0;
			for (int m = 0; m < indicator_cols[s].size(); m++) {
				//cout << "++++ indicator_cols+1 " << m << " out of " << indicator_cols[s].size() << endl;
				if (int(features.masks[indicator_cols[s][m]][i]) > 0)
					missing = missing + 1;
			}
			if (2 * missing >= indicator_cols[s].size()) {
				p1 = p1 + pow(2, s);
				mask_stat[s]++;
			}
		}
		predictor_stat[p1]++;

		// predict per sample
		x.get_row(i, sample);
		for (int j = 0; j < cols_per_predictor[p1].size(); j++) {
			sample4p[j] = sample[cols_per_predictor[p1][j]];
		}

		predictors[p1]->predict(sample4p, score, 1, cols_per_predictor[p1].size());
		
		// calibrate the score per sample
		vector<MedSample> sampleWprediction(1);
		sampleWprediction[0].id = 0;
		sampleWprediction[0].prediction = { score[0] };
		calibrators[p1]->Apply(sampleWprediction);
		features.samples[i].prediction = sampleWprediction[0].prediction; 

		features.samples[i].prediction.push_back(score[0]);
		features.samples[i].attributes["mask_number"] = p1;
		// cout << "sample: " << i << " predictor " << p1 << " score BEFORE calibration: " << score[0] << " and AFTER: " << features.samples[i].prediction[0] << endl;

	}

	// report some stat
	cout << endl << "MASKs frequency" << endl;
	for (int m = 0; m < masks.size(); m++) {
		cout << (float) 100 * mask_stat[m] / x.get_nrows() << "% " << masks[m] << endl;
	}
	cout << endl << "Predictors frequency - predictor 0 is no mask" << endl;
	for (int m = 0; m < p; m++) {
		cout << m << " " << (float)100 * predictor_stat[m] / x.get_nrows() << "%" << endl;
	}
	cout << endl;

	return 0;
}

