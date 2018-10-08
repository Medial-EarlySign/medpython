#include "MPTime.h"
#include "MedTime/MedTime/MedTime.h"

const int MPTime::Undefined = MedTime::Undefined;			///< undefined time unit
const int MPTime::Date = MedTime::Date;						///< dates are in full regular format YYYYMMDD
const int MPTime::Years = MedTime::Years;					///< years since 1900 (not since 0!)
const int MPTime::Months = MedTime::Months;					///< months since 1900/01/01
const int MPTime::Days = MedTime::Days;						///< days since 1900/01/01
const int MPTime::Hours = MedTime::Hours;					///< hours since 1900/01/01
const int MPTime::Minutes = MedTime::Minutes;				///< minutes since 1900/01/01
const int MPTime::DateTimeString = MedTime::DateTimeString;	///< string only format "YYYYMMDDHHMI"
