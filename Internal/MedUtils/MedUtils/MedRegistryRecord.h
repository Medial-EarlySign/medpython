#ifndef __MED_REGISTRY_RECORD_H__
#define __MED_REGISTRY_RECORD_H__
#include <SerializableObject/SerializableObject/SerializableObject.h>
/**
* A class which represnt a registry record of patient in time range from start_date to end_date
* It has min_allowed and max_allowed dates for sampling and it has optional mark for
* current sample age
*/
class MedRegistryRecord : public SerializableObject
{
public:
	int pid; ///< patient ID
			 //defines the registry value apply date range
	int start_date; ///< the start_date range for the record
	int end_date; ///< the end_date range for the record

	float registry_value; ///< the registry value/state

	ADD_CLASS_NAME(MedRegistryRecord)
	ADD_SERIALIZATION_FUNCS(pid, start_date, end_date, registry_value)
};

MEDSERIALIZE_SUPPORT(MedRegistryRecord)

#endif
