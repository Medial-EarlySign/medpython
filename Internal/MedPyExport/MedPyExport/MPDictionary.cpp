#include "MPDictionary.h"
#include "MPPidRepository.h"

#include <time.h>
#include <string>

#include "InfraMed/InfraMed/MedConvert.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/Utils.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedModel.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"


MPDictionary::MPDictionary(MPPidRepository* rep): o(rep->o) { }

int MPDictionary::section_id(string name) { return o->dict.section_id(name); }

int MPDictionary::id(string& signame) { return o->dict.id(signame); }

string MPDictionary::name(int id) { return o->dict.name(id); };

void MPDictionary::prep_sets_lookup_table(int section_id, const std::vector<string>& set_names, MEDPY_NP_OUTPUT(char** lut_array, int* lut_size)) {
	std::vector<char> lut;
	int retval = o->dict.prep_sets_lookup_table(section_id, set_names, lut);
	if (retval != 0) {
		throw runtime_error(str(boost::format("dict.prep_sets_lookup_table failed (retval= %1 )") % retval));
	}
	
	*lut_array = (char*)malloc(sizeof(char)*lut.size());
	*lut_size = (int)lut.size();
	memcpy(*lut_array, lut.data(), lut.size());
}

