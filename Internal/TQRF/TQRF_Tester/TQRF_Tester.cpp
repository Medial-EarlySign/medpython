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
#include <Medutils/MedUtils/MedGenUtils.h>
#include <TQRF/TQRF/TQRF.h>

#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;
namespace po = boost::program_options;

#include <External/Eigen/Core>
#include <External/Eigen/Dense>
using namespace Eigen;


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
			("preds_file", po::value<string>()->default_value(""), "preds file to write")
			("seed", po::value<int>()->default_value(-1), "random seed (-1 is time)")
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

	set_rand_seed(vm["seed"].as<int>());

	return 0;
}


//
// just playing with matrices to learn how to do stuff
void stam() {
	
	MatrixXf a(3, 2);
	MatrixXf b(2, 1);
	MatrixXf bsq(2, 1);
	MatrixXf c(3, 1);
	MatrixXf d(3, 2);
	MatrixXf e(1, 3);
	MatrixXf f(1, 2);


	a << 1, 2, 3, 4, 5, 6;
	b << 7, 5;

	cout << "a:" << endl << a << endl;
	cout << "b:" << endl << b << endl;

	c = a * b;

	cout << "c:" << endl << c << endl;

	d = a * b.asDiagonal();

	cout << "d:" << endl << d << endl;

	bsq = b.array()*b.array();

	cout << "bsq:" << endl << bsq << endl;

	MedMat<float> A(3, 2);
	A.m ={ (float)1.1,(float)2.2,(float)3.3,(float)4.4,(float)5.5,(float)6.6 };
	Map<MatrixXf> Af(&A.m[0], 2, 3);

	cout << "Af:" << endl << Af << endl;

	cout << "Af(0,1): " << Af(0, 1) << " a(0,1): " << a(0,1) << " A(0,1): " << A(0,1) <<  endl;

	e = b.transpose() * Af;
	cout << "e:" << endl << e << endl;

	f = a.colwise().sum();
	cout << "f=colsum a:" << endl << f << endl;


	exit(0);
}

//========================================================================================
// MAIN
//========================================================================================

int main(int argc, char *argv[])
{
	//stam();

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
		if (samples_test.read_from_file(vm["samples_test"].as<string>()) < 0) {
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
	MLOG("pids: train %d test: %d all: %d\n", pids_train.size(), pids_test.size(), pids.size());
	vector<string> signals;
	model.get_required_signal_names(signals);

	MedPidRepository rep;
	if (rep.read_all(vm["rep"].as<string>(), pids, signals) < 0) {
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

	MLOG(">>> Model built ... now predicting on test samples (if given)\n");
	if (vm["samples_test"].as<string>() != "") {

		model.features.clear();
		// run model
		model.apply(rep, samples_test, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_APPLY_FTR_PROCESSORS);

		// get a matrix
		MedMat<float> x;
		model.features.get_as_matrix(x);

		MLOG(">>> Got predict matrix x : %d x %d\n", x.nrows, x.ncols);

		// get tqrf predictions
		vector<float> preds;
		tqrf.Predict(x, preds);

		MLOG(">>> Got %d tqrf predictions : n_per_sample %d (%d ids)\n", preds.size(), tqrf.n_preds_per_sample(), preds.size()/tqrf.n_preds_per_sample());

		// add predictions to MedSamples
		int i = 0;
		int n_per_sample = tqrf.n_preds_per_sample();
		for (auto &Id : samples_test.idSamples)
			for (auto &s : Id.samples) {
				s.prediction.clear();
				for (int j=0; j<n_per_sample; j++)
					s.prediction.push_back(preds[i*n_per_sample + j]);
				i++;
			}

		// write to preds file
		samples_test.write_to_file(vm["preds_file"].as<string>());
		MLOG(">>> Wrote tqrf predictions to file %s\n", vm["preds_file"].as<string>().c_str());
	}

	return 0;
}


