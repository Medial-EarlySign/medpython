//
// TARF_Tester
//


#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <InfraMed/InfraMed/InfraMed.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <TQRF/TQRF/TQRF.h>

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;
namespace po = boost::program_options;

//=========================================================================================================
int read_run_params(int argc, char *argv[], po::variables_map& vm) {
	po::options_description desc("Program options");

	try {
		desc.add_options()
			("help", "produce help message")
			("rep", po::value<string>()->default_value("/home/Repositories/THIN/thin_jun2017/thin.repository"), "repository file name")
			("samples_train", po::value<string>()->default_value(""), "samples file to train with")
			("samples_test", po::value<string>()->default_value(""), "samples file to test with")
			("model", po::value<string>()->default_value(""), "model file to generate features")
			("tqrf_params", po::value<string>()->default_value(""), "tqrf_params")
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


//========================================================================================
// MAIN
//========================================================================================

int main(int argc, char *argv[])
{
	int rc = 0;
	po::variables_map vm;

	// Reading run Parameters
	MLOG("Reading params\n");
	rc = read_run_params(argc, argv, vm);
	assert(rc >= 0);


	// Read train/test samples
	MedSamples samples_train, samples_test;
	if (samples_train.read_from_file(vm["samples_train"].as<string>()) < 0) {
		MERR("ERROR: failed reading samples file %s\n", vm["samples_train"].as<string>().c_str());
		return -1;
	}

	if (vm["samples_test"].as<string>() != "") {
		if (samples_train.read_from_file(vm["samples_test"].as<string>()) < 0) {
			MERR("ERROR: failed reading samples file %s\n", vm["samples_test"].as<string>().c_str());
			return -1;
		}
	}

	// Prepare a model with our feature generator
	MedModel model;
	if (model.read_from_file(vm["model"].as<string>()) < 0) {
		MERR("ERROR: failed reading model file %s\n", vm["model"].as<string>().c_str());
		return -1;
	}


	// get repository data for train
	vector<int> pids_train, pids_test, pids;
	samples_train.get_ids(pids_train);
	samples_test.get_ids(pids_test);
	pids = pids_train;
	pids.insert(pids.end(), pids_test.begin(), pids_test.end());
	vector<string> signals;
	model.get_required_signal_names(signals);

	MedPidRepository rep;
	if (rep.read_all(vm["rep"].as<string>(), pids_train, signals) < 0) {
		MERR("ERROR: failed reading repository %s\n", vm["rep"].as<string>().c_str());
		return -1;
	}


	MLOG("================================================================================================\n");
	MLOG(">>> Running Apply on model ....\n");

	// generate features using model
	model.apply(rep, samples_train, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_APPLY_FTR_PROCESSORS);

	MLOG(">>> After Apply on model ....\n");

	// initializing and training a tqrf model
	TQRF_Forest tqrf;

	MLOG(">>> Before tqrf init ....\n");
	tqrf.init_from_string(vm["tqrf_params"].as<string>());

	MLOG(">>> Before tqrf Train ....\n");
	tqrf.Train(model.features);

	return 0;
}


