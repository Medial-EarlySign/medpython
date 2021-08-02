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

int MPDictionary::id(int section_id, string& signame) { return o->dict.id(section_id, signame); }

string MPDictionary::name(int section_id, int id) { return o->dict.name(section_id, id); };

void MPDictionary::prep_sets_lookup_table(int section_id, const std::vector<string>& set_names, MEDPY_NP_OUTPUT(char** lut_array, unsigned long long* lut_size)) {
	std::vector<char> lut;
	int retval = o->dict.prep_sets_lookup_table(section_id, set_names, lut);
	if (retval != 0) {
		throw runtime_error(str(boost::format("dict.prep_sets_lookup_table failed (retval= %1 )") % retval));
	}
	
	*lut_array = (char*)malloc(sizeof(char)*lut.size());
	*lut_size = lut.size();
	memcpy(*lut_array, lut.data(), lut.size());
}

MPIntVecIntMapAdaptor MPDictionary::get_members_to_all_sets(int section_id, MEDPY_NP_INPUT(int* members_array, unsigned long long members_size)) {
	vector<int> members;
	buf_to_vector(members_array, members_size, members);
	unordered_map<int, std::vector<int>> Member2AllSets;
	o->dict.dict(section_id)->get_members_to_all_sets(members, Member2AllSets);
	MPIntVecIntMapAdaptor ret;
	for (auto& i : Member2AllSets) {
		ret.o->insert(i);
	}
	return ret;
}