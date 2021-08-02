#ifndef __MP_Dictionary_H
#define __MP_Dictionary_H

#include "MedPyCommon.h"

class MPPidRepository;
class MedPidRepository;

class MPDictionary {
public:
	MEDPY_IGNORE(MedPidRepository* o);

	MPDictionary(MPPidRepository* rep);

	int section_id(string name);
	int id(int section_id, string& signame);
	string name(int section_id, int id);
	void prep_sets_lookup_table(int section_id, const std::vector<string>& set_names, MEDPY_NP_OUTPUT(char** lut_array, unsigned long long* lut_size));

	MPIntVecIntMapAdaptor get_members_to_all_sets(int section_id, MEDPY_NP_INPUT(int* members_array, unsigned long long members_size));
};




#endif // !__MP_Dictionary_H
