#include "MedTime.h"

MedTime med_time_converter;

// implementations

//.....................................................................................................
void MedTime::init_time_tables()
{

	// date to days tables

	YearsMonths2Days.resize(3001*100, -1);
	Years2Days.resize(3001, -1);

	for (int year = 1900; year<=3000; year++) {
		// Full years
		int days = 365 * (year-1900);
		days += (year-1897)/4;
		days -= (year-1801)/100;
		days += (year-1601)/400;

		Years2Days[year-1900] = days;
		//if (year<1905) fprintf(stderr, "y2d[%d] = %d\n", year-1900, days);
		YearsMonths2Days[year*100 + 0] = days; // month 0

		// month 0 - for lazy people !!
		YearsMonths2Days[year*100] = days;

		// months 1-12
		for (int month = 1; month<=12; month++) {
			int ym = year*100 + month;

			// Full Months
			int d = days + days2month[month-1];
			if (month>2 && (year%4)==0 && ((year%100)!=0 || (year%400)==0))
				d++;
			YearsMonths2Days[ym] = d;
			//if (year<1905) fprintf(stderr, "ym2d[%d] = %d\n", ym, d);
		}
	}

	// days to dates tables
	Days2Years.resize(1100*365, -1); // covering dates up to 3000
	Days2Months.resize(1100*365, -1); // covering dates up to 3000
	Days2Date.resize(1100*365, -1); // covering dates up to 3000
	for (int d=0; d<1100*365; d++) {
		// Full Years
		int year = 1900 + d / 365;
		int days = d % 365;

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

		Days2Years[d] = year - 1900;

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

		Days2Months[d] = (year-1900)*12 + month;

		days++;

		Days2Date[d] = days + 100 * month + 10000 * year;

	}

}

//.....................................................................................................
int MedTime::convert_days(int to_type, int in_time)
{
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Days) return in_time;
	if (to_type == MedTime::Hours) return in_time*24;
	if (to_type == MedTime::Minutes) return in_time*24*60;
	if (to_type == MedTime::Date) return Days2Date[in_time];
	if (to_type == MedTime::Months) return Days2Months[in_time];
	if (to_type == MedTime::Years) return Days2Years[in_time];

	return -1;
}

//.....................................................................................................
int MedTime::convert_date(int to_type, int in_time)
{
	if (to_type == MedTime::Date) return in_time;
	if (to_type == MedTime::Years) return in_time/10000 - 1900;
	if (to_type == MedTime::Months) return ((in_time/10000)-1900)*12 + (in_time%10000)/100 - 1;

	int ym = in_time/100;
	int days = (in_time % 100) - 1;

	days += YearsMonths2Days[ym];
	//fprintf(stderr, "it %d ym %d days %d ym2d %d\n", in_time, ym, days, YearsMonths2Days[ym]);

	return convert_days(to_type, days);
}

//.....................................................................................................
int MedTime::convert_years(int to_type, int in_time)
{
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Date) return ((in_time+1900)*10000 + 101);
	if (to_type == MedTime::Years) return in_time;
	if (to_type == MedTime::Months) return in_time*12;

	return convert_days(to_type, Years2Days[in_time]);
}

//.....................................................................................................
int MedTime::convert_months(int to_type, int in_time)
{
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Months) return in_time;

	int year = 1900 + (in_time/12);
	if (to_type == MedTime::Years) return year;

	int month = 1 + (in_time % 12);
	int ym = year * 100 + month;

	if (to_type == MedTime::Date) return (ym*100 + 1);

	int days = YearsMonths2Days[ym];
	return convert_days(to_type, days);
}

//.....................................................................................................
int MedTime::convert_hours(int to_type, int in_time)
{
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Hours) return in_time;
	if (to_type == MedTime::Minutes) return in_time*60;
	int days = in_time/24;
	return convert_days(to_type, days);
}

//.....................................................................................................
int MedTime::convert_minutes(int to_type, int in_time)
{
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Minutes) return in_time;
	int hours = in_time/60;
	return convert_hours(to_type, hours);
}

//.....................................................................................................
int MedTime::convert_times(int from_type, int to_type, int in_time)
{
	if (from_type == MedTime::Date) return convert_date(to_type, in_time);
	if (from_type == MedTime::Days) return convert_days(to_type, in_time);
	if (from_type == MedTime::Minutes) return convert_minutes(to_type, in_time);
	if (from_type == MedTime::Hours) return convert_hours(to_type, in_time);
	if (from_type == MedTime::Months) return convert_months(to_type, in_time);
	if (from_type == MedTime::Years) return convert_years(to_type, in_time);
	return -1;
}

//.....................................................................................................
int MedTime::convert_times(int from_type, int to_type, double in_time)
{
	if (from_type == MedTime::Date) return convert_date(to_type, (int)in_time);

	if (from_type == MedTime::Minutes) return convert_minutes(to_type, (int)in_time);
	if (from_type == MedTime::Hours) return convert_minutes(to_type, (int)(60.0*in_time));
	if (from_type == MedTime::Days) return convert_minutes(to_type, (int)(60.0*24.0*in_time));
	if (from_type == MedTime::Months) return convert_minutes(to_type, (int)((365.0/12.0)*24.0*60.0*in_time));
	if (from_type == MedTime::Years) return convert_minutes(to_type, (int)(365.0*24.0*60.0*in_time));

	return -1;
}


//.....................................................................................................
double MedTime::convert_times_D(int from_type, int to_type, int in_time)
{
	if (to_type == MedTime::Date || (from_type!=MedTime::Date && to_type>=from_type)) return (double)convert_times(from_type, to_type, in_time);
	if (from_type == MedTime::Date && (to_type >= MedTime::Days)) return (double)convert_times(from_type, to_type, in_time);

	int minutes1 = convert_times(from_type, MedTime::Minutes, in_time);
	int int_time = convert_times(from_type, to_type, in_time);
	int minutes2 = convert_times(to_type, MedTime::Minutes, int_time);

	double res = (double)int_time;

	//fprintf(stderr, "m1 %d it %d m2 %d\n", minutes1, int_time, minutes2);

	if (to_type == MedTime::Years) { res += (double)(minutes1 - minutes2)/(365.0*24.0*60.0); return res; }
	if (to_type == MedTime::Months) { res += (double)(minutes1 - minutes2)/((365.0/12.0)*24.0*60.0); return res; }
	if (to_type == MedTime::Days) { res += (double)(minutes1 - minutes2)/(24.0*60.0); return res; }
	if (to_type == MedTime::Hours) { res += (double)(minutes1 - minutes2)/60.0; return res; }

	return -1;
}

//.....................................................................................................
double MedTime::convert_times_D(int from_type, int to_type, double in_time)
{
	if (from_type == MedTime::Date) return (double)convert_date(to_type, (int)in_time);
	if (from_type == MedTime::Minutes) return (double)convert_minutes(to_type, (int)in_time);

	double it = 0;

	if (from_type == MedTime::Hours) it = in_time * 60.0;
	else if (from_type == MedTime::Days) it = in_time * 60.0*24.0;
	else if (from_type == MedTime::Months) it = in_time * ((365.0/12.0)*24.0*60.0);
	else if (from_type == MedTime::Years) it = in_time * 365.0*24.0*60.0;

	return (double)convert_minutes(to_type, (int)it);
}

//.....................................................................................................
int MedTime::string_to_type(const string &str)
{
	if (str == "Date" || str == "date") return MedTime::Date;
	if (str == "Years" || str == "years" || str == "Year" || str == "year") return MedTime::Years;
	if (str == "Months" || str == "months" || str == "Month" || str == "month") return MedTime::Months;
	if (str == "Days" || str == "days" || str == "Day" || str == "day") return MedTime::Days;
	if (str == "Hours" || str == "hours" || str == "Hour" || str == "hour") return MedTime::Hours;
	if (str == "Minutes" || str == "minutes" || str == "Minute" || str == "minute") return MedTime::Minutes;
	return -1;
}
//.....................................................................................................
int MedTime::add_subtruct_days(int in_time, int delta_days) {
	int in_time_n = convert_times(Date, Days, in_time);
	in_time_n += delta_days;
	int out_time = convert_times(Days, Date, in_time_n);
	return out_time;
}