#pragma once
//
// Containing a list of common operations done on top of MedProcessTools objects
//

#include <MedProcessTools/MedProcessTools/MedModel.h>
#include <InfraMed/InfraMed/InfraMed.h>
#include <InfraMed/InfraMed/MedPidRepository.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>

namespace medial {
	namespace model {
		// gets f_model, f_samples, f_rep, reads them all, reads into rep only the required signals x needed pids.
		int prep_model_rep_samples(string f_model, string f_samples, string f_rep, MedPidRepository &rep, MedSamples &samples, MedModel &model);
	}
}

