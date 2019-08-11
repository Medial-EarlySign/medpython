#include <MedProcessTools/MedProcessTools/MedProcessWrappers.h>
#include <Logger/Logger/Logger.h>

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace std;

//---------------------------------------------------------------------------------------------------------------------------------------------
// gets f_model, f_samples, f_rep, reads them all, reads into rep only the required signals x needed pids.
int medial::model::prep_model_rep_samples(string f_model, string f_samples, string f_rep, MedPidRepository &rep, MedSamples &samples, MedModel &model)
{
	// read samples
	if (samples.read_from_file(f_samples) < 0) {
		MERR("ERROR: medial::model : Could not read samples file %s\n", f_samples.c_str());
		return -1;
	}

	// read model
	if (model.read_from_file(f_model) < 0) {
		MERR("ERROR: medial::model : Could not read model file %s\n", f_model.c_str());
		return -1;
	}

	// get pids and required sigs
	vector<int> pids;
	vector<string> required;

	samples.get_ids(pids);
	if (rep.init(f_rep) < 0)
		MTHROW_AND_ERR("ERROR could not read repository %s\n", f_rep.c_str());
	model.fit_for_repository(rep);

	model.get_required_signal_names(required);

	if (rep.read_all(f_rep, pids, required) < 0) {
		MERR("ERROR: medial::model : Could not read repository %s\n", f_rep.c_str());
		return -1;
	}

	return 0;
}
