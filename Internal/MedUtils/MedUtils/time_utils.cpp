#define _CRT_SECURE_NO_WARNINGS
#include "MedGenUtils.h"

int days2month[] = {0,31,59,90,120,151,181,212,243,273,304,334,365} ;

int get_day_approximate(int val) {
	if (val <= 19000101 || val >= 21000101)
		return -1;
	if (get_day(val) != -1)
			return get_day(val);
	if (get_day((val/100)*100 + 1) != -1)
		return get_day((val / 100) * 100 + 1);
	if (get_day((val / 10000) * 10000 + 101) != -1)
		return get_day((val / 10000) * 10000 + 101);
	fprintf(stderr, "unable to convert %d to days \n", val);
	throw exception();
}
int get_day(int val) {
	int year = val/100/100 ;
	int month = (val/100)%100 ;
	int day = val%100 ;

	if (month < 1 || month > 12 || year < 1900 || day < 0)
		return -1;

	// Full years
	int days = 365 * (year-1900) ;
	days += (year-1897)/4 ;
	days -= (year-1801)/100 ;
	days += (year-1601)/400 ;

	// Full Months

	int days2month[] = {0,31,59,90,120,151,181,212,243,273,304,334,365} ;
	days += days2month[month-1] ;
	if (month>2 && (year%4)==0 && ((year%100)!=0 || (year%400)==0))
		days ++ ;
	days += (day-1) ;

	return days;
}

// Days -> Date
int get_date(int days) {

	// Full Years
	int year = 1900 + days / 365;
	days %= 365;

	days -= (year - 1897) / 4;
	days += (year - 1801) / 100;
	days -= (year - 1601) / 400;

	if (days < 0) {
		year--;
		days += 365;
		if ((year % 4) == 0 && ((year % 100) != 0 || (year % 400) == 0)) {
			days++;
			if (days == 366) {
				days = 0;
				year++;
			}
		}
	}

	// Full Months
	bool leap_year = ((year % 4) == 0 && ((year % 100) != 0 || (year % 400) == 0));
	int month;
	for (int i = 1; i <= 12; i++) {
		int mdays = days2month[i] + ((leap_year && i > 1) ? 1 : 0);
		if (days < mdays) {
			month = i;
			days -= (days2month[i - 1] + ((leap_year && i > 2) ? 1 : 0));
			break;
		}
	}

	days++;
	return days + 100 * month + 10000 * year;
}
