#define _CRT_SECURE_NO_WARNINGS
#define _CRT_RAND_S

#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;


#include "Example.h"

int main(int argc, char *argv[])
{
	int rc = 0;
	po::variables_map vm;

	MedTimer timer;
	timer.start();

	// Running Parameters
	MLOG( "Reading params\n");
	rc = read_run_params(argc, argv, vm);
	assert(rc >= 0);

	// Read Signals
	MLOG( "Reading signals\n");

	vector<string> signals;
	rc = read_signals_list(vm, signals);
	assert(rc >= 0);

	// Read Samples
	MLOG( "Reading samples\n");

	MedSamples allSamples;
	get_samples(vm, allSamples);
	vector<int> ids;
	allSamples.get_ids(ids);

	// Read Repository
	MLOG("Initializing repository\n");

	MedPidRepository rep;
	rc = read_repository(vm, ids, signals, rep);
	assert(rc >= 0);

	timer.take_curr_time();
	MLOG("Reading params + rep time: %f sec\n", timer.diff_sec());

	// Define Model
	MedModel my_model;

	timer.start();

#define DIRECT_INIT 1
	if (DIRECT_INIT) {
		MLOG("Initializing RepCleaners and Features: nsignals: %d , n_ids: %d\n", signals.size(), allSamples.idSamples.size());
		for (auto sig : signals) {

			// cleaner for sig
			//my_model.add_process_to_set(0, "rp_type=nbrs_cln;take_log=0;signal="+sig);
			my_model.add_process_to_set(0, "rp_type=basic_cln;take_log=1;signal="+sig);

			// features for sig

			// basics
			my_model.add_process_to_set(0, "fg_type=basic; type=last; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=last; win_from=0; win_to=360; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=last; win_from=360; win_to=720; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=last; win_from=720; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=win_delta; win_from=0; win_to=180; d_win_from=360; d_win_to=720; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=win_delta; win_from=0; win_to=180; d_win_from=720; d_win_to=1080; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=first; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=last2; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=avg; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=max; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=min; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=std; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=last_time; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=last_time2; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);

			my_model.add_process_to_set(0, "fg_type=basic; type=slope; win_from=0; win_to=720; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=slope; win_from=0; win_to=10000; time_unit=Days; signal=" + sig);
			my_model.add_process_to_set(0, "fg_type=basic; type=slope; win_from=720; win_to=10000; time_unit=Days; signal=" + sig);

			//// binnedLM
			my_model.add_process_to_set(0, "fg_type=binnedLM; estimation_points=1440,720,360,180; signal=" + sig);
		}
		// Age/Gender features
		my_model.add_process_to_set(0, "fg_type=age");
		my_model.add_process_to_set(0, "fg_type=gender");
		//my_model.add_process_to_set(0, "fg_type=basic; type=category_set_sum; win_from=0; win_to=720; time_unit=Days; sets=ATC_A10_____; signal=Drug");

		// Add feature processors
		my_model.add_process_to_set(0, "fp_type=basic_cleaner");
		my_model.add_process_to_set(1, "fp_type=imputer;strata=Age,40,80,5;moment_type=0");
		//my_model.add_process_to_set(2, "fp_type=normalizer");

	}
	else {
		// Repository Cleaners
		MLOG("Initializing RepCleaners : nsignals: %d , n_ids: %d\n", signals.size(), allSamples.idSamples.size());
		my_model.add_rep_processors_set(REP_PROCESS_NBRS_OUTLIER_CLEANER, signals, vm["rep_cleaner_params"].as<string>());
		//	my_model.add_rep_processors_set(REP_PROCESS_BASIC_OUTLIER_CLEANER, signals, vm["rep_cleaner_params"].as<string>());

		// Signal-based feature generators
		MLOG("Initializing Features\n");

		vector<BasicFeatureTypes> sig_types ={ FTR_LAST_VALUE, FTR_FIRST_VALUE, FTR_LAST2_VALUE, FTR_AVG_VALUE, FTR_MAX_VALUE, FTR_MIN_VALUE, FTR_STD_VALUE, FTR_LAST_DELTA_VALUE , FTR_LAST_DAYS, FTR_LAST2_DAYS };
		//	vector<BasicFeatureTypes> sig_types = { FTR_LAST_VALUE };
		for (auto sig_type : sig_types) 
			my_model.add_feature_generators(FTR_GEN_BASIC, signals, "win_from=0; win_to=10000; type = " + std::to_string(sig_type));	
		my_model.add_feature_generators(FTR_GEN_BINNED_LM, signals, string("estimation_points=800,400,180"));

		// Age + Gender
		MLOG( "Initializing Extra Features\n");
		my_model.add_age();
		my_model.add_gender();

		// Add feature cleaners
		my_model.add_feature_processors_set(FTR_PROCESS_BASIC_OUTLIER_CLEANER, vm["feat_cleaner_params"].as<string>());

		// Add feature Imputers
		MLOG("Adding imputers\n");
		my_model.add_imputers("strata=Age,40,80,5");

		// Normalizers
		MLOG("Adding normalizers\n");
		my_model.add_normalizers();
	}



	// Predictor
	MLOG( "Initializing Predictor\n");
	my_model.set_predictor(vm["predictor"].as<string>(), vm["predictor_params"].as<string>());
	assert(my_model.predictor != NULL);

	timer.take_curr_time();
	MLOG("Init model time: %f sec\n", timer.diff_sec());

	if (vm.count("nfolds")) {
		// Cross Validator
		CrossValidator cv(&my_model);
		timer.start();
		MLOG( "Initializing Filters\n");

		BasicTrainFilter *trainFilter = new BasicTrainFilter;
		BasicTestFilter *testFilter = new BasicTestFilter;

		cv.add_learning_set_filter(trainFilter);
		cv.add_test_set_filter(testFilter);

		MLOG("Init filters time: %g sec\n", timer.diff_sec());

		int nfolds = vm["nfolds"].as<int>();
		MLOG( "Performing %d-fold cross-validation\n", nfolds);

		MedSamples cvOutSamples;

		if (cv.doCV(rep, allSamples, nfolds, cvOutSamples) < 0) {
			MLOG( "CV failed\n");
			return -1;
		}

		// analyze
		vector<float> y, preds;
		for (auto& idSample : cvOutSamples.idSamples) {
			for (auto& sample : idSample.samples) {
				//MLOG("Id=%d\t%f", idSample.id, sample.outcome);
				for (int i = 0; i < sample.prediction.size(); i++) {
					y.push_back(sample.outcome);
					preds.push_back(sample.prediction[i]);
					//MLOG("\t%f", sample.prediction[i]);
				}
				//MLOG("\n");
			}
		}

		float AUC = get_preds_auc(preds, y);
		MLOG("y size: %d , preds size: %d , cv AUC is : %f\n", y.size(), preds.size(), AUC);

	}
	else {
		if (vm.count("temp_file") == 0) {
			fprintf(stderr, "temp-file required if nfold not given\n");
			return -1;
		}
		string tempFile = vm["temp_file"].as<string>();

		// Learn on 50%; Predict on rest
		BasicTrainFilter trainFilter;
		BasicTestFilter testFilter;

		// Learning and test set
		MedSamples learningSamples, testSamples;
		for (auto& sample : allSamples.idSamples) {
			if (globalRNG::rand() % 2 == 0)
				learningSamples.idSamples.push_back(sample);
			else
				testSamples.idSamples.push_back(sample);
		}

		// Filter
		trainFilter.filter(rep, learningSamples);
		testFilter.filter(rep, testSamples);

		// Learn Model
		if (my_model.learn(rep, &learningSamples) < 0) {
			fprintf(stderr,"Learning model failed\n");
			return -1;
		}

		// Write to temporary file
		my_model.write_to_file(tempFile);
		fprintf(stderr, "Done writing to file %s\n", tempFile.c_str());

		// Read from temporary file
		MedModel new_model;
		new_model.read_from_file(tempFile);
		fprintf(stderr, "Done reading from file %s\n", tempFile.c_str());

		// Apply
		if (new_model.apply(rep, testSamples) < 0) {
			fprintf(stderr,"Applying model failed\n");
			return -1;
		}

		// analyze
		vector<float> y, preds;
		for (auto& idSample : testSamples.idSamples) {
			for (auto& sample : idSample.samples) {
				for (int i = 0; i < sample.prediction.size(); i++) {
					y.push_back(sample.outcome);
					preds.push_back(sample.prediction[i]);
				}
			}
		}

		float AUC = get_preds_auc(preds, y);
		MLOG("y size: %d , preds size: %d , cv AUC is : %f\n", y.size(), preds.size(), AUC);
	}
	
	return 0;

}

// Functions 
// Analyze Command Line
int read_run_params(int argc, char *argv[], po::variables_map& vm) {
	po::options_description desc("Program options");

	try {
		desc.add_options()
			("help", "produce help message")
			("config", po::value<string>()->required(), "repository file name")
			("ids",po::value<string>(),"file of ids to consider")
			("samples", po::value<string>()->required(), "samples file name")
			("features", po::value<string>()->required(), "file of signals to consider")
			("rep_cleaner", po::value<string>(), "repository cleaner")
			("rep_cleaner_params", po::value<string>()->default_value(""), "repository cleaner params")
			("feat_cleaner", po::value<string>(), "features cleaner")
			("feat_cleaner_params", po::value<string>()->default_value(""), "features cleaner params")
			("predictor", po::value<string>()->default_value("linear_model"), "predictor")
			("predictor_params", po::value<string>()->default_value(""), "predictor params")
			("temp_file",po::value<string>(), "temporary file for serialization")

			("nfolds", po::value<int>(), "number of cross-validation folds")
			;

		po::store(po::parse_command_line(argc, argv, desc), vm);
		if (vm.count("help")) {
			cerr << desc << "\n";
			exit(-1);

		}
		po::notify(vm);

		MLOG("=========================================================================================\n");
		MLOG("Command Line:");
		for (int i=0; i<argc; i++) MLOG(" %s", argv[i]);
		MLOG("\n");
		MLOG("..........................................................................................\n");
	}
	catch (exception& e) {
		cerr << "error: " << e.what() << "; run with --help for usage information\n";
		return -1;
	}
	catch (...) {
		cerr << "Exception of unknown type!\n";
		return -1;
	}

	return 0;
}

int read_repository(po::variables_map& vm, vector<int>& ids, vector<string>& signals, MedPidRepository& rep) {

	vector<string> sigs = signals;
	sigs.push_back("GENDER");
	sigs.push_back("BYEAR");
	sigs.push_back("TRAIN");
	MLOG("Before reading config file %s\n", vm["config"].as<string>().c_str());
	string config_file = vm["config"].as<string>();

	if (rep.read_all(config_file,ids,sigs) < 0) {
		MLOG("Cannot init repository %s\n", config_file.c_str());
		return -1;
	}

	size_t nids = rep.index.pids.size();
	size_t nsigs = rep.index.signals.size();
	MLOG("Read %d Ids and %d signals\n", (int)nids, (int)nsigs);

	return 0;
}

int read_signals_list(po::variables_map& vm, vector<string>& signals) {
	
	string file_name = vm["features"].as<string>();
	ifstream inf(file_name);

	if (!inf) {
		MLOG("Cannot open %s for reading\n", file_name.c_str());
		return -1;
	}

	string curr_line;
	while (getline(inf, curr_line)) {
		if (curr_line[curr_line.size() - 1] == '\r')
			curr_line.erase(curr_line.size() - 1);

		if (curr_line[0] != '#')
			signals.push_back(curr_line);
	}

	return 0;
}

int read_ids_list(po::variables_map& vm, vector<int>& ids) {

	string file_name = vm["ids"].as<string>();
	ifstream inf(file_name);

	if (!inf) {
		MLOG("Cannot open %s for reading\n", file_name.c_str());
		return -1;
	}

	string curr_line;
	while (getline(inf, curr_line)) {
		if (curr_line[curr_line.size() - 1] == '\r')
			curr_line.erase(curr_line.size() - 1);

		ids.push_back(stoi(curr_line));
	}

	return 0;
}

int get_samples(po::variables_map& vm, MedSamples& samples) {

	string file_name = vm["samples"].as<string>();
	return samples.read_from_file(file_name);
}
