//
// MedRepositoryType.h
//
// A class used to handle the differences between General-Practice repostiories (THIN, MHS), Hospital repositories (MIMIC) and potentially other
//
// The library defines a global initiated instance called med_rep_type
//
//

#ifndef __MED__REP_TYPE__H__
#define __MED__REP_TYPE__H__

#include <string>
#include <MedTime/MedTime/MedTime.h>
#include "Logger/Logger/Logger.h"

using namespace std;

typedef enum {
	REP_TYPE_GP,
	REP_TYPE_HOSPITAL,
	REP_TYPE_LAST
} MedRepTypes;

class MedRepositoryType {

public:
	// type
	MedRepTypes repType = REP_TYPE_GP;

	// Differences
	string genderSignalName;
	bool ageDirectlyGiven;
	int basicTimeUnit;

	// Constructor
	MedRepositoryType() { setRepositoryType(REP_TYPE_GP); }
	int setRepositoryType(MedRepTypes type);
};

extern MedRepositoryType med_rep_type;


#endif
