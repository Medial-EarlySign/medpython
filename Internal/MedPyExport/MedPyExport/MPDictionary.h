#ifndef __MP_Dictionary_H
#define __MP_Dictionary_H

#include "MedPyCommon.h"

class MPPidRepository;
class MedPidRepository;

class MPDictionary {
public:
	MedPidRepository* o;

	MPDictionary(MPPidRepository* rep);

	int section_id(string name);
	int id(string& signame);
	string name(int id);
	void prep_sets_lookup_table(int section_id, const std::vector<string>& set_names, MEDPY_NP_OUTPUT(char** lut_array, int* lut_size));
};




#endif // !__MP_Dictionary_H
