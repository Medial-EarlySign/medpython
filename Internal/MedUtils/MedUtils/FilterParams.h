#ifndef __FILTER_PARAMS_H__
#define __FILTER_PARAMS_H__

#include <map>
#include <string>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include "LabelParams.h"
using namespace std;

/**
* A Class to constraint sampling options for creating medsamples prior to sample process.
* Gender,Train constraints can be applied later with FilterSamples, no need to be prior sampling constraints.
*/
class FilterParams : public SerializableObject {
public:
	int min_age; ///< Filter on min age
	int max_age; ///< Filter on max age
	int min_time; ///< Filter on min time
	int max_time; ///< Filter on max time
	TimeWindowMode interaction_mode; ///< the interaction mode in filtering

	FilterParams();
	/// init function
	int init(map<string, string>& map);
};

#endif