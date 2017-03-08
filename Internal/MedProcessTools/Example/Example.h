
#include "MedProcessTools/MedProcessTools/MedModel.h"
#include "MedFeat/MedFeat/MedFeat.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedStat/MedStat/MedStat.h"
#include "MedStat/MedStat/MedPerformance.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/CrossValidator.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"

#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;


int read_run_params(int argc, char *argv[], po::variables_map& vm);
int read_repository(string config_file, vector<int>& ids, vector<string> &signals, MedPidRepository& rep);
int read_signals_list(po::variables_map& vm, vector<string>& signals);
int read_ids_list(po::variables_map& vm, vector<int>& ids);
int get_samples(po::variables_map& vm, MedSamples& samples);