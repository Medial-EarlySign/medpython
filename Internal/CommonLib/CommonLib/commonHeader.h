// A common header for all projects in solution

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_RAND_S
#define _SCL_SECURE_NO_WARNINGS

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "InfraMed/InfraMed/MedConvert.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/Utils.h"
#include <stdio.h>

#include <stdlib.h>
#include <stdarg.h>

#include <boost/program_options.hpp>

#include "MedStat/MedStat/MedStat.h"
#include "MedStat/MedStat/MedPerformance.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedUtils/MedUtils/MedCohort.h"
#include "MedProcessTools/MedProcessTools/MedModel.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"

#include <boost/serialization/map.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <omp.h>

#include <Logger/Logger/Logger.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedIO/MedIO/MedIO.h>

#include <algorithm>
#include <time.h>
#include <string.h>


#define DEFAULT_REP "/home/Repositories/THIN/thin_final/thin.repository"

using namespace std;
namespace po = boost::program_options; 

// Functions
void shuffleMatrix(MedFeatures& matrix);

// Read a matrix from csv/bin
void readMatrix(MedFeatures& features, po::variables_map& vm, string csvName, string binName);
void readMatrix(MedFeatures& features, po::variables_map& vm);

// Write predictions to csv/bin/samples
void writePredictions(MedFeatures& features, po::variables_map& vm);
// Read a list of names
int readList(string& fname, vector<string>& list);
// Read predictor from file
void read_predictor_from_file(MedPredictor* &pred, string& predictorFile);
void get_folds(po::variables_map& vm, vector<int>& folds, int nFolds);

// Get performance
void print_auc_performance(MedSamples& samples, int nfolds, string& outFile);
void print_auc_performance(MedSamples& samples, vector<int>& folds, string& outFile);
void get_performance(MedSamples& samples, vector<Measurement>& msr, vector<vector<float>>& prfs);

// For hyper-parameters optimization
void get_options(string& paramsFile, int nRuns, vector<map<string, string> >& options);
int read_optimization_ranges(string& optFile, map<string, vector<string> >& optimizationOptions);
void print_performance(ofstream& of, vector<string>& predictorParams, vector<Measurement>& msrs, vector<map<string, string> >& predictorOptions, vector<int>& folds, vector<vector<vector<float>>>& all_prfs);

// Split matrix
void get_features(MedFeatures& inMatrix, MedFeatures& outMatrix, int iFold, bool isLearning);

