#include "commonHeader.h"

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL


// Read A matrix from csv/bin
void readMatrix(MedFeatures& features, po::variables_map& vm) { 
	readMatrix(features, vm, "inCsv", "inBin"); 
}

void readMatrix(MedFeatures& features, po::variables_map& vm, string csvName, string binName) {

	if (vm.count(csvName) + vm.count(binName) != 1)
		MTHROW_AND_ERR("Exactly one of %s and %s should be given\n",csvName.c_str(),binName.c_str());

	if (vm.count(csvName)) {
		string inFile = vm[csvName].as<string>();
		if (features.read_from_csv_mat(inFile) != 0)
			MTHROW_AND_ERR("Cannot read features from csv file %s\n", inFile.c_str());
	}
	else {
		string inFile = vm[binName].as<string>();
		if (features.read_from_file(inFile) != 0)
			MTHROW_AND_ERR("Cannot read features from bin file %s\n", inFile.c_str());

	}
}

// Write predictions to samples/features-bin/features-csv
void writePredictions(MedFeatures& features, po::variables_map& vm)
{

	if (vm.count("outSamples") + vm.count("outCsv") + vm.count("outBin") != 1)
		MTHROW_AND_ERR("Exactly one of ousSamples, outCsv, outBin should be given\n");

	if (vm.count("outSamples")) {
		MedSamples samples;
		features.get_samples(samples);
		string fileName = vm["outSamples"].as<string>();
		if (samples.write_to_file(fileName) != 0)
			MTHROW_AND_ERR("Cannot write samples to %s\n", fileName.c_str());
	}
	else if (vm.count("outCsv")) {
		string fileName = vm["outCsv"].as<string>();
		if (features.write_as_csv_mat(fileName) != 0)
			MTHROW_AND_ERR("Cannot write features to %s\n", fileName.c_str());
	}
	else if (vm.count("outBin")) {
		string fileName = vm["outBin"].as<string>();
		if (features.write_to_file(fileName) != 0)
			MTHROW_AND_ERR("Cannot write features to %s\n", fileName.c_str());
	}
}

// Read List from file
int readList(string& fname, vector<string>& list) {

	ifstream inf(fname);
	if (!inf) {
		cerr << "Cannot open " << fname << " for reading\n";
		return -1;
	}

	string line;
	list.clear();
	while (getline(inf, line))
		list.push_back(line);

	cerr << "Read " << list.size() << " members\n";
	inf.close();
	return 0;
}

// Read predictor from file
void read_predictor_from_file(MedPredictor* &pred, string& predictorFile) {
	unsigned char *blob;
	unsigned long long final_size;
	if (read_binary_data_alloc(predictorFile, blob, final_size) != 0)
		MTHROW_AND_ERR("Cannot read binary data file %s\n", predictorFile.c_str());

	MedPredictorTypes type;
	memcpy(&type, blob, sizeof(MedPredictorTypes));
	pred = MedPredictor::make_predictor(type);
	pred->deserialize(blob + sizeof(MedPredictorTypes));
}

// Get Folds To Take
void get_folds(po::variables_map& vm, vector<int>& folds, int nFolds) {

	if (vm.count("folds")) {
		string folds_s = vm["folds"].as<string>();
		vector<string> folds_v;
		boost::split(folds_v, folds_s, boost::is_any_of(","));

		set<int> _folds;
		for (unsigned int i = 0; i < folds_v.size(); i++) {
			int fold = stoi(folds_v[i]);
			if (fold >= nFolds || fold < 0)
				MTHROW_AND_ERR("Illegal fold \'%s\' in folds\n", folds_v[i].c_str());
			if (_folds.find(fold) != _folds.end())
				MTHROW_AND_ERR("Duplicate fold \'%s\' in folds \'%s\'\n", folds_v[i].c_str(),folds_s.c_str());
			folds.push_back(fold);
		}
	}
	else {
		for (int i = 0; i < nFolds; i++)
			folds.push_back(i);
	}

}


// Performance
void print_auc_performance(MedSamples& samples, int nfolds, string& outFile) {
	vector<int> folds;
	for (int i = 0; i < nfolds; i++)
		folds.push_back(i);
	print_auc_performance(samples, folds, outFile);
}

void print_auc_performance(MedSamples& samples, vector<int>& folds, string& outFile) {

	FILE *fp = fopen(outFile.c_str(), "w");
	if (fp == NULL)
		MTHROW_AND_ERR("Cannot open %s for writing\n", outFile.c_str());

	MedClassifierPerformance prf(samples);

	Measurement auc("AUC"), sensAt90("Sens", "Spec", 0.9F), specAt50("Spec", "Sens", 0.5F);
	prf.GetPerformanceParam(auc); prf.GetPerformanceParam(sensAt90);  prf.GetPerformanceParam(specAt50);

	vector<float> MeasureSums(3, 0.0);
	int nfolds = (int)folds.size();
	for (int i = 0; i <= nfolds; i++) {
		if (i > 0) {
			fprintf(fp, "Split %d: AUC = %f Sens@90 = %f Spec@50 = %f, \n", folds[i], prf.MeasurementValues[auc][i], prf.MeasurementValues[sensAt90][i], prf.MeasurementValues[specAt50][i]);
			MeasureSums[0] += prf.MeasurementValues[auc][i];
			MeasureSums[1] += prf.MeasurementValues[sensAt90][i];
			MeasureSums[2] += prf.MeasurementValues[specAt50][i];
		} else
			fprintf(fp, "Combined Used Splits: AUC = %f Sens@90 = %f Spec@50 = %f, \n",  prf.MeasurementValues[auc][i], prf.MeasurementValues[sensAt90][i], prf.MeasurementValues[specAt50][i]);
	}

	fprintf(fp, "Mean of  Used Splits: AUC = %f Sens@90 = %f Spec@50 = %f, \n", MeasureSums[0] / nfolds, MeasureSums[1] / nfolds, MeasureSums[2] / nfolds);


	// Sanity Check
	vector<int> nums(nfolds);
	vector<double> sums(nfolds);
	map<int, int> folds_m;
	for (int i = 0; i < nfolds; i++)
		folds_m[folds[i]] = i;

	for (MedIdSamples& idSample : samples.idSamples) {
		int fold = idSample.split;
		if (folds_m.find(fold) != folds_m.end()) {
			int split = folds_m[fold];
			for (MedSample& sample : idSample.samples) {
				nums[split] ++;
				sums[split] += sample.prediction[0];
			}
		}
	}
	for (int split = 0; split < nfolds; split++)
		fprintf(stderr, "Mean predictions for split %d = %.3f\n", folds[split], sums[split] / nums[split]);
}

void get_performance(MedSamples& samples, vector<Measurement>& msrs, vector<vector<float>>& prfs) {
	if (msrs.empty())
		return;

	MedClassifierPerformance prf(samples);
	for (auto& msr : msrs)
		prf.GetPerformanceParam(msr);

	size_t nMsrs = msrs.size();
	size_t nRes = prf.MeasurementValues[msrs[0]].size();
	prfs.resize(nMsrs, vector<float>(nRes, 0));
	for (size_t i = 0; i < nRes; i++) {
		for (size_t j = 0; j < nMsrs; j++)
			prfs[j][i] = prf.MeasurementValues[msrs[j]][i];
	}
}

void shuffleMatrix(MedFeatures& matrix) {

	vector<string> names;
	for (auto& rec : matrix.data)
		names.push_back(rec.first);

	cerr << "Shuffling " << names.size() << " vectors" << endl;

#pragma omp parallel for
	for (int i = 0; i <names.size(); i++) {
		shuffle(matrix.data[names[i]].begin(), matrix.data[names[i]].end(), globalRNG::get_engine());
	}

	cerr << "Shuffling label" << endl;

	vector<float> labels(matrix.samples.size());
	for (int i = 0; i < labels.size(); i++)
		labels[i] = matrix.samples[i].outcome;

	shuffle(labels.begin(), labels.end(), globalRNG::get_engine());

	for (int i = 0; i < labels.size(); i++)
		matrix.samples[i].outcome = labels[i];
}

// Functions for hyper-parameters optimization
// Get options list
void get_options(string& paramsFile, int nRuns, vector<map<string, string> >& options) {

	map<string, vector<string> > optimizationOptions;
	int nOptions = read_optimization_ranges(paramsFile, optimizationOptions);
	if (nRuns < 0 || nRuns > nOptions)
		nRuns = nOptions;
	float prob = (nRuns + 0.0) / nOptions;

	if (prob < 0.5) {
		// Select randomly
		unordered_set<string> selected;
		while (selected.size() < nRuns) {
			string option_s;
			map<string, string> option;

			for (auto& rec : optimizationOptions) {
				int index = rec.second.size() * (globalRNG::rand() / (globalRNG::max() + 0.0));
				option[rec.first] = rec.second[index];
				option_s += "," + rec.second[index];
			}

			if (selected.find(option_s) == selected.end()) {
				selected.insert(option_s);
				options.push_back(option);
			}
		}
	}
	else {
		// Select systematically
		vector<string> paramNames;
		for (auto& rec : optimizationOptions)
			paramNames.push_back(rec.first);
		int nParams = (int)paramNames.size();
		vector<int> pointers(nParams, 0);
		while (1) {
			map<string, string> option;
			for (int i = 0; i < nParams; i++)
				option[paramNames[i]] = optimizationOptions[paramNames[i]][pointers[i]];
			options.push_back(option);

			// Next
			int idx = (int)pointers.size() - 1;
			while (idx >= 0 && pointers[idx] == optimizationOptions[paramNames[idx]].size() - 1)
				idx--;

			if (idx < 0)
				break;

			pointers[idx] ++;
			for (int j = idx + 1; j < nParams; j++)
				pointers[j] = 0;
		}

		shuffle(options.begin(), options.end(), globalRNG::get_engine());
		options.resize(nRuns);
		options.shrink_to_fit();
	}

	cerr << "# of selected combinations = " << options.size() << "\n";
}

int read_optimization_ranges(string& optFile, map<string, vector<string> >& optimizationOptions) {

	ifstream inf(optFile);
	if (!inf)
		MTHROW_AND_ERR("Cannot open %s for reading\n", optFile.c_str());

	string line;
	int nOptions = 0;
	vector<string> fields;

	while (getline(inf, line)) {
		boost::split(fields, line, boost::is_any_of("\t"));
		if (fields.size() == 4) {
			float val = stof(fields[1]);
			float to = stof(fields[2]);
			float step = stof(fields[3]);
			while (val <= to) {
				// Keep int as int ...
				if (int(val) == val)
					optimizationOptions[fields[0]].push_back(to_string(int(val)));
				else
					optimizationOptions[fields[0]].push_back(to_string(val));
				val += step;
			}

			if (nOptions == 0)
				nOptions = 1;
			nOptions *= (int)optimizationOptions[fields[0]].size();
		}
		else if (fields.size() == 2) {
			vector<string> options;
			boost::split(options, fields[1], boost::is_any_of(","));
			for (string& val : options)
				optimizationOptions[fields[0]].push_back(val);

			if (nOptions == 0)
				nOptions = 1;
			nOptions *= (int)optimizationOptions[fields[0]].size();

		}
		else
			MTHROW_AND_ERR("Cannot parse input line \'%s\'\n", line.c_str());
	}

	if (nOptions == 0) {
		MTHROW_AND_ERR("Cannot optimzize without options\n");
	}
	else
		MLOG("Options file contains %d combinations\n", nOptions);

	return nOptions;
}

void print_performance(ofstream& of, vector<string>& predictorParams, vector<Measurement>& msrs, vector<map<string, string> >& predictorOptions, vector<int>& folds, vector<vector<vector<float>>>& all_prfs) {

	int nSplits = (int)folds.size();

	//Header
	of << "#";
	for (string& param : predictorParams)
		of << "\t" << param;
	for (auto& msr : msrs) {
		string msr_name = msr.name();
		for (int i = 0; i < nSplits; i++)
			of << "\t" << msr_name << "_Split" << folds[i];
		of << "\t" << msr_name << "_Mean";
	}
	of << "\n";

	for (int iOption = 0; iOption < predictorOptions.size(); iOption++) {
		of << iOption;

		for (string& param : predictorParams)
			of << "\t" << predictorOptions[iOption][param];

		for (int iMsr = 0; iMsr < msrs.size(); iMsr++) {
			float sum = 0;
			for (int iSplit = 0; iSplit < nSplits; iSplit++) {
				of << "\t" << all_prfs[iSplit][iOption][iMsr];
				sum += all_prfs[iSplit][iOption][iMsr];
			}
			of << "\t" << sum / nSplits;
		}
		of << "\n";
	}
}

// Get learning/test matrix
void get_features(MedFeatures& inMatrix, MedFeatures& outMatrix, int iFold, bool isLearning) {

	MedTimer timer;
	timer.start();
	vector<string> feature_names;
	inMatrix.get_feature_names(feature_names);
	for (string& name : feature_names)
		outMatrix.data[name].clear();
	outMatrix.samples.clear();

	for (auto& attr : inMatrix.attributes)
		outMatrix.attributes[attr.first] = attr.second;

	vector<int> inds;
	for (unsigned int i = 0; i < inMatrix.samples.size(); i++) {
		auto& sample = inMatrix.samples[i];
		if ((isLearning && sample.split != iFold) || ((!isLearning) && sample.split == iFold))
			inds.push_back(i);
	}

	for (string& name : feature_names) {
		outMatrix.data[name].resize(inds.size());
		for (unsigned int i = 0; i < inds.size(); i++)
			outMatrix.data[name][i] = inMatrix.data[name][inds[i]];
	}

	outMatrix.samples.resize(inds.size());
	for (unsigned int i = 0; i < inds.size(); i++)
		outMatrix.samples[i] = inMatrix.samples[inds[i]];

	/*
	for (unsigned int i = 0; i<inMatrix.samples.size(); i++) {
		auto& sample = inMatrix.samples[i];
		if ((isLearning && sample.split != iFold) || ((!isLearning) && sample.split == iFold)) {
			outMatrix.samples.push_back(sample);
			for (string& name : feature_names)
				outMatrix.data[name].push_back(inMatrix.data[name][i]);
		}
	}
	*/
	timer.take_curr_time();
	MLOG("get-features took %f seconds\n", timer.diff_sec());
}