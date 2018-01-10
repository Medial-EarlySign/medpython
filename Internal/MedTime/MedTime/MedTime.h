/// @file
/// MedTime.h
///
/// time conversion routines (years, date, days, hours, minutes)
///
/// The library defines a global initiated instance called med_time_converter which can be used to convert.
///
///

#ifndef __MED__TIME__H__
#define __MED__TIME__H__

#include <vector>
#include <string>
using namespace std;

/// The following is a class to handle time conversions
/// There's a global instance of it declared below (med_time_convert) to be used in conversions.
class MedTime {

public:

	// names for timing options
	const static int Undefined = 0;			///< undefined time unit
	const static int Date = 1;				///< dates are in full regular format YYYYMMDD
	const static int Years = 2;				///< years since 1900 (not since 0!)
	const static int Months = 3;			///< months since 1900/01/01
	const static int Days = 4;				///< days since 1900/01/01
	const static int Hours = 5;				///< hours since 1900/01/01
	const static int Minutes = 6;			///< minutes since 1900/01/01

	vector<int> YearsMonths2Days;
	vector<int> Years2Days;
	vector<int> Days2Years;
	vector<int> Days2Months;
	vector<int> Days2Date;
	vector<int> days2month ={ 0,31,59,90,120,151,181,212,243,273,304,334,365 };

	MedTime() { init_time_tables(); }

	// general converters
	// when times are given/returned in int they are rounded to the floor
	// then times are given/returned in double they are in fractions of the given/requested unit

	/// <summary>
	/// Converts time formats. usefull for converting dates to days from 1900 for example to
	/// calculate days diff.
	/// </summary>
	/// <returns>
	/// returns the time in new format
	/// </returns>
	int convert_times(int from_type, int to_type, int in_time);
	int convert_times(int from_type, int to_type, double in_time);
	double convert_times_D(int from_type, int to_type, int in_time);
	double convert_times_D(int from_type, int to_type, double in_time);


	void init_time_tables();

	int convert_date(int to_type, int in_time);
	int convert_years(int to_type, int in_time);
	int convert_months(int to_type, int in_time);
	int convert_days(int to_type, int in_time);
	int convert_hours(int to_type, int in_time);
	int convert_minutes(int to_type, int in_time);
	int add_subtruct_days(int in_time, int delta_days);

	int string_to_type(const string &str);

};

extern MedTime med_time_converter;


#endif