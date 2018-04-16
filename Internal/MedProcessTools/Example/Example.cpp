#define _CRT_SECURE_NO_WARNINGS
#define _CRT_RAND_S

#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "Example.h"

#include <fenv.h>
#ifndef  __unix__
#pragma float_control( except, on )
#endif


void run_learn_apply(MedPidRepository &rep, MedSamples &allSamples, po::variables_map &vm, vector<string> signals)
{
	// Define Model
	MedModel my_model;
	string init_file = vm["model_init_file"].as<string>();
	vector<string> dummy;
	my_model.init_from_json_file_with_alterations(init_file, dummy);

	MedTimer timer;
	timer.start();
}

int main(int argc, char *argv[])
{

#if defined(__unix__)
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

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

	timer.take_curr_time();
	MLOG("Reading params time: %f sec\n", timer.diff_sec());

	MedPidRepository rep;

	run_learn_apply(rep, allSamples, vm, signals);

	return 0;
}

// Functions 
// Analyze Command Line
int read_run_params(int argc, char *argv[], po::variables_map& vm) {
	po::options_description desc("Program options");

	try {
		desc.add_options()
			("help", "produce help message")
			("config", po::value<string>()->default_value("W:\\CancerData\\Repositories\\THIN\\thin_jun2017\\thin.repository"), "repository file name")
			("ids",po::value<string>(),"file of ids to consider")
			("samples", po::value<string>()->required(), "samples file name")
			("sigs", po::value<string>()->default_value("NONE"), "list of signals to consider")
			("importance", "run importance when using qrf model")
			("drug_feats", "add drug based features to model")
			("csv_feat", po::value<string>()->default_value("NONE"), "file name to save features as csv (NONE = no saving)")
			("sigs_file", po::value<string>()->default_value("NONE"), "file of signals to consider")
			("rep_cleaner", po::value<string>(), "repository cleaner")
			("rep_cleaner_params", po::value<string>()->default_value(""), "repository cleaner params")
			("feat_cleaner", po::value<string>(), "features cleaner")
			("feat_cleaner_params", po::value<string>()->default_value(""), "features cleaner params")
			("predictor", po::value<string>()->default_value("linear_model"), "predictor")
			("predictor_params", po::value<string>()->default_value(""), "predictor params")
			("temp_file", po::value<string>(), "temporary file for serialization")
			("direct_init", po::value<int>()->default_value(0), "temporary file for serialization")
			("model_init_file", po::value<string>()->default_value("H:\\MR\\Libs\\Internal\\MedProcessTools\\Config_Examples\\model_json_version2_example.json"), "init json file for entire model")
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

int read_repository(string config_file, vector<int>& ids, vector<string>& signals, MedPidRepository& rep) {

	vector<string> sigs = signals;
	sigs.push_back("GENDER");
	sigs.push_back("BYEAR");
	sigs.push_back("TRAIN");
	MLOG("Before reading config file %s\n", config_file.c_str());


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
	
	string sigs = vm["sigs"].as<string>();
	if (sigs != "NONE") {
		signals.clear();
		boost::split(signals, sigs, boost::is_any_of(","));
		return 0;
	}

	string file_name = vm["sigs_file"].as<string>();
	if (file_name == "NONE")
		return 0;
	
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
