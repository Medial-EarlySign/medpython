/// @file
/// registry methods over MedRegistry Object
#ifndef __MED_ENUM_H__
#define __MED_ENUM_H__
#include <vector>
#include <string>

using namespace std;

/// @enum
/// Conflicting options
enum class ConflictMode {
	All = 0,///< "all" - take all
	Drop = 1, ///< "drop" - drop when conflcit
	Max = 2, ///< "max" - take max on conflict 
	Last = 3, ///< "last" - take last value
	Bitwise_Max = 4 ///< "bitwise_max" - max on each multilabel entry
};

extern vector<string> ConflictMode_to_name;
ConflictMode ConflictMode_name_to_type(const string& ConflictMode_name);

/// @enum
/// Sampling options
enum class TimeWindowMode {
	Before_End = 0, ///< "before_end" - need to be before end_time of registry
	Before_Start = 1, ///< "before_start" - need to be before start_time of registry
	After_Start = 2, ///< "after_start" - need to be after start_time of registry
	Within = 3, ///< "within" - need to be within start_time and end_time - contained fully time window. 
	All_ = 4 ///< "all" - takes all not testing for anything
			 //no None, after_end - useless for now
};
extern vector<string> TimeWindow_to_name;
TimeWindowMode TimeWindow_name_to_type(const string& TimeWindow_name);

#endif
