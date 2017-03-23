//
// MedRepositoryType.cpp
//

#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "MedRepositoryType.h"

#define LOCAL_SECTION LOG_REPTYPE
#define LOCAL_LEVEL	LOG_DEF_LEVEL

MedRepositoryType med_rep_type;

int MedRepositoryType::setRepositoryType(MedRepTypes type) {
	if (type == REP_TYPE_GP) {
		genderSignalName = "GENDER";
		ageDirectlyGiven = false;
		basicTimeUnit = MedTime::Date;
		windowTimeUnit = MedTime::Days;
	}
	else if (type == REP_TYPE_HOSPITAL) {
		genderSignalName = "Gender";
		ageDirectlyGiven = true;
		basicTimeUnit = MedTime::Minutes;
		windowTimeUnit = MedTime::Minutes;
	}
	else {
		MERR("Unknown Repository Type %d\n", type);
		return -1;
	}

	return 0;
}