#include "MedTime.h"
#include <Logger/Logger/Logger.h>
#include <boost/algorithm/string.hpp>
#include <SerializableObject/SerializableObject/SerializableObject.h>
//#include <MedUtils/MedUtils/MedUtils.h>

#define LOCAL_SECTION LOG_MED_UTILS
#define LOCAL_LEVEL	LOG_DEF_LEVEL

MedTime med_time_converter;
int global_default_time_unit = MedTime::Date;
int global_default_windows_time_unit = MedTime::Days;

// implementations

bool MedTime::is_valid_date(int date) const {
	int year = int(date / 10000);
	bool valid = year >= 1900 && year <= 3000;

	if (valid) {
		//check month:
		int month = int(date / 100) % 100;
		valid = month >= 1 && month <= 12;
		if (valid) { //check day:
			int day = date % 100;
			int days_in_month = days2month[month] - days2month[month - 1];
			if ((month == 2) &&  //check years with 29 days
				(year % 4) == 0 && ((year % 100) != 0 || (year % 400) == 0))
					++days_in_month;
			
			valid = day >= 1 && day <= days_in_month;
		}
	}

	return valid;
}

//.....................................................................................................
void MedTime::init_time_tables()
{

	// date to days tables

	YearsMonths2Days.resize(3001 * 100, -1);
	Years2Days.resize(3001, -1);

	for (int year = 1900; year <= 3000; year++) {
		// Full years
		int days = 365 * (year - 1900);
		days += (year - 1897) / 4;
		days -= (year - 1801) / 100;
		days += (year - 1601) / 400;

		Years2Days[year - 1900] = days;
		//if (year<1905) fprintf(stderr, "y2d[%d] = %d\n", year-1900, days);
		YearsMonths2Days[year * 100 + 0] = days; // month 0

		// month 0 - for lazy people !!
		YearsMonths2Days[year * 100] = days;

		// months 1-12
		for (int month = 1; month <= 12; month++) {
			int ym = year * 100 + month;

			// Full Months
			int d = days + days2month[month - 1];
			if (month > 2 && (year % 4) == 0 && ((year % 100) != 0 || (year % 400) == 0))
				d++;
			YearsMonths2Days[ym] = d;
			//if (year<1905) fprintf(stderr, "ym2d[%d] = %d\n", ym, d);
		}
	}

	// days to dates tables
	Days2Years.resize(1100 * 365, -1); // covering dates up to 3000
	Days2Months.resize(1100 * 365, -1); // covering dates up to 3000
	Days2Date.resize(1100 * 365, -1); // covering dates up to 3000
	for (int d = 0; d < 1100 * 365; d++) {
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
		int month = 0;
		for (int i = 1; i <= 12; i++) {
			int mdays = days2month[i] + ((leap_year && i > 1) ? 1 : 0);
			if (days < mdays) {
				month = i;
				days -= (days2month[i - 1] + ((leap_year && i > 2) ? 1 : 0));
				break;
			}
		}

		Days2Months[d] = (year - 1900) * 12 + month - 1;

		days++;

		Days2Date[d] = days + 100 * month + 10000 * year;

	}

}

//.....................................................................................................
int MedTime::convert_days(int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_days - to_time_type is MedTime::Undefined");
	if (to_type == MedTime::Days) return in_time;
	if (to_type == MedTime::Hours) return in_time * 24;
	if (to_type == MedTime::Minutes) return in_time * 24 * 60;
	if (in_time < 0) MTHROW_AND_ERR("Error in MedTime::convert_days - tried to convert negative date %d\n", in_time);
	if (to_type == MedTime::Date) return Days2Date[in_time];
	if (to_type == MedTime::Months) return Days2Months[in_time];
	if (to_type == MedTime::Years) return Days2Years[in_time];

	MTHROW_AND_ERR("Error in MedTime::convert_days - unsupported type value %d\n", to_type);
	return -1;
}

int MedTime::convert_datetime_safe(int to_type, string in_time, char handle_ilegal_date) {
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_datetime_safe - to_time_type is MedTime::Undefined");
	static int warning_count = 0;
	string out_t = in_time;

	boost::replace_all(out_t, " ", "");
	boost::replace_all(out_t, "-", "");
	boost::replace_all(out_t, ":", "");
	if (out_t.size() > 12)
		out_t = out_t.substr(0, 12);

	int date_part = med_stoi(out_t.substr(0, 8));
	int year_part = int(date_part / 10000);
	if (year_part < 1900) {
		if (handle_ilegal_date > 1)
			MTHROW_AND_ERR("Error in MedTime:convert_datetime_safe - recieved date before 1900. date_part=%d\n",
				date_part);
		if (handle_ilegal_date > 0 && warning_count < 5)
			MWARN("Warning in MedTime:convert_datetime_safe - recieved date before 1900."
				" date_part=%d, truncating to 19000101\n", date_part);
		++warning_count;
		date_part = 19000101;
	}

	if (to_type == MedTime::Minutes) {
		int minutes = convert_date(to_type, date_part);
		if (out_t.size() >= 10) minutes += med_stoi(out_t.substr(8, 2)) * 60;
		if (out_t.size() >= 12) minutes += med_stoi(out_t.substr(10, 2));
		return minutes;
	}
	else return convert_date(to_type, date_part);
}

/// handles YYYYMMDD and YYYY-MM-DD HH:MI:SS formats
int MedTime::convert_datetime(int to_type, string in_time) {
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_datetime - to_time_type is MedTime::Undefined");
	string out_t = in_time;
	boost::replace_all(out_t, " ", "");
	boost::replace_all(out_t, "-", "");
	boost::replace_all(out_t, ":", "");
	if (out_t.size() > 12)
		out_t = out_t.substr(0, 12);

	int date_part = med_stoi(out_t.substr(0, 8));
	if (to_type == MedTime::Minutes) {
		int minutes = convert_date(to_type, date_part);
		if (out_t.size() >= 10) minutes += med_stoi(out_t.substr(8, 2)) * 60;
		if (out_t.size() >= 12) minutes += med_stoi(out_t.substr(10, 2));
		return minutes;
	}
	else return convert_date(to_type, date_part);
}

//.....................................................................................................
int MedTime::convert_date(int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_date - to_time_type is MedTime::Undefined");
	if (to_type == MedTime::Date) return in_time;
	if (to_type == MedTime::Years) return in_time / 10000 - 1900;
	if (to_type == MedTime::Months) return ((in_time / 10000) - 1900) * 12 + (in_time % 10000) / 100 - 1;

	// ihadanny - removing this obscure code that tries to guess that you actually meant days instead of date
	//if (in_time >= 30000000)
	//	return convert_days(to_type, 1100 * 365);

	int ym = in_time / 100;
	int days = (in_time % 100) - 1;

	days += YearsMonths2Days[ym];
	//fprintf(stderr, "it %d ym %d days %d ym2d %d\n", in_time, ym, days, YearsMonths2Days[ym]);

	return convert_days(to_type, days);
}

//.....................................................................................................
int MedTime::convert_years(int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_years - to_time_type is MedTime::Undefined");
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Date) return ((in_time + 1900) * 10000 + 101);
	if (to_type == MedTime::Years) return in_time;
	if (to_type == MedTime::Months) return in_time * 12;

	return convert_days(to_type, Years2Days[in_time]);
}

//.....................................................................................................
int MedTime::convert_months(int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_months - to_time_type is MedTime::Undefined");
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Months) return in_time;

	int year = 1900 + (in_time / 12);
	if (to_type == MedTime::Years) return year;

	int month = 1 + (in_time % 12);
	int ym = year * 100 + month;

	if (to_type == MedTime::Date) return (ym * 100 + 1);

	int days = YearsMonths2Days[ym];
	return convert_days(to_type, days);
}

//.....................................................................................................
int MedTime::convert_hours(int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_hours - to_time_type is MedTime::Undefined");
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Hours) return in_time;
	if (to_type == MedTime::Minutes) return in_time * 60;
	int days = in_time / 24;
	return convert_days(to_type, days);
}

//.....................................................................................................
int MedTime::convert_minutes(int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_minutes - to_time_type is MedTime::Undefined");
	if (in_time < 0) in_time = 0;
	if (to_type == MedTime::Minutes) return in_time;
	int hours = in_time / 60;
	return convert_hours(to_type, hours);
}

//.....................................................................................................
string MedTime::convert_times_S(int from_type, int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_times_S - to_time_type is MedTime::Undefined");
	if (to_type == MedTime::DateTimeString) {
		int d = convert_times(from_type, MedTime::Date, in_time);
		if (from_type == MedTime::Minutes || from_type == MedTime::Hours) {
			int total_m = convert_times(from_type, MedTime::Minutes, in_time);
			char buff[20];
			snprintf(buff, sizeof(buff), "%d-%02d:%02d", d, (total_m / 60) % 24, total_m % 60);
			return buff;
		}
		else return to_string(d);
	}
	else
		return to_string(convert_times(from_type, to_type, in_time));
}

//.....................................................................................................
int MedTime::convert_times(int from_type, int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_times - to_time_type is MedTime::Undefined");
	if (from_type == MedTime::Date) return convert_date(to_type, in_time);
	if (from_type == MedTime::Days) return convert_days(to_type, in_time);
	if (from_type == MedTime::Minutes) return convert_minutes(to_type, in_time);
	if (from_type == MedTime::Hours) return convert_hours(to_type, in_time);
	if (from_type == MedTime::Months) return convert_months(to_type, in_time);
	if (from_type == MedTime::Years) return convert_years(to_type, in_time);

	MTHROW_AND_ERR("Error in MedTime::convert_times - unsupported from_type value %d to_type %d\n", 
		from_type, to_type);
	return -1;
}

//.....................................................................................................
int MedTime::convert_times(int from_type, int to_type, double in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_times - to_time_type is MedTime::Undefined");
	if (from_type == MedTime::Date) return convert_date(to_type, (int)in_time);

	if (from_type == MedTime::Minutes) return convert_minutes(to_type, (int)in_time);
	if (from_type == MedTime::Hours) return convert_minutes(to_type, (int)(60.0*in_time));
	if (from_type == MedTime::Days) return convert_minutes(to_type, (int)(60.0*24.0*in_time));
	if (from_type == MedTime::Months) return convert_minutes(to_type, (int)((365.0 / 12.0)*24.0*60.0*in_time));
	if (from_type == MedTime::Years) return convert_minutes(to_type, (int)(365.0*24.0*60.0*in_time));

	MTHROW_AND_ERR("Error in MedTime::convert_times - unsupported from_type value %d\n", from_type);
	return -1;
}


//.....................................................................................................
double MedTime::convert_times_D(int from_type, int to_type, int in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_times_D - to_time_type is MedTime::Undefined");
	if (to_type == MedTime::Date || (from_type != MedTime::Date && to_type >= from_type)) return (double)convert_times(from_type, to_type, in_time);
	if (from_type == MedTime::Date && (to_type >= MedTime::Days)) return (double)convert_times(from_type, to_type, in_time);

	int minutes1 = convert_times(from_type, MedTime::Minutes, in_time);
	int int_time = convert_times(from_type, to_type, in_time);
	int minutes2 = convert_times(to_type, MedTime::Minutes, int_time);

	double res = (double)int_time;

	//fprintf(stderr, "m1 %d it %d m2 %d\n", minutes1, int_time, minutes2);

	if (to_type == MedTime::Years) { res += (double)(minutes1 - minutes2) / (365.0*24.0*60.0); return res; }
	if (to_type == MedTime::Months) { res += (double)(minutes1 - minutes2) / ((365.0 / 12.0)*24.0*60.0); return res; }
	if (to_type == MedTime::Days) { res += (double)(minutes1 - minutes2) / (24.0*60.0); return res; }
	if (to_type == MedTime::Hours) { res += (double)(minutes1 - minutes2) / 60.0; return res; }

	MTHROW_AND_ERR("Error in MedTime::convert_times_D - unsupported to_type value %d\n", to_type);
	return -1;
}

//.....................................................................................................
double MedTime::convert_times_D(int from_type, int to_type, double in_time)
{
	if (to_type == MedTime::Undefined) MTHROW_AND_ERR("Error in MedTime::convert_times_D - to_time_type is MedTime::Undefined");
	if (from_type == MedTime::Date) return (double)convert_date(to_type, (int)in_time);
	if (from_type == MedTime::Minutes) return (double)convert_minutes(to_type, (int)in_time);

	double it = 0;

	if (from_type == MedTime::Hours) it = in_time * 60.0;
	else if (from_type == MedTime::Days) it = in_time * 60.0*24.0;
	else if (from_type == MedTime::Months) it = in_time * ((365.0 / 12.0)*24.0*60.0);
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
string MedTime::type_to_string(int type)
{
	string res = "Unknown";
	switch (type)
	{
	case MedTime::Date:
		res = "Date";
		break;
	case MedTime::Years:
		res = "Years";
		break;

	case MedTime::Months:
		res = "Months";
		break;
	case MedTime::Days:
		res = "Days";
		break;
	case MedTime::Hours:
		res = "Hours";
		break;
	case MedTime::Minutes:
		res = "Minutes";
		break;
	case MedTime::Undefined:
		res = "Undefined";
		break;
	case MedTime::DateTimeString:
		res = "DateTimeString";
		break;
	default:
		MTHROW_AND_ERR("Error in MedTime::type_to_string - not implement name for %d\n", type);
	}
	return res;
}
//.....................................................................................................
int MedTime::add_subtruct_days(int in_time, int delta_days) {
	int in_time_n = convert_times(Date, Days, in_time);
	in_time_n += delta_days;
	int out_time = convert_times(Days, Date, in_time_n);
	return out_time;
}


//.....................................................................................................
int MedTime::diff_times(int d1, int d2, int in_type, int out_type)
{
	int t1 = convert_times(in_type, out_type, d1);
	int t2 = convert_times(in_type, out_type, d2);
	return (t1 - t2);
}

//.....................................................................................................
double MedTime::diff_times_D(int d1, int type1, int d2, int type2, int out_type)
{
	double t1 = convert_times_D(type1, type2, d1);
	double d = t1 - (double)d2;
	double t = convert_times_D(type2, out_type, d);

	return t;
}

//.....................................................................................................
int MedTime::diff_times(int d1, int type1, int d2, int type2, int out_type)
{
	int t1 = convert_times(type1, type2, d1);
	int d = t1 - d2;
	int t = convert_times(type2, out_type, d);

	return t;
}

//.....................................................................................................
double MedTime::get_age(int t, int type_t, int byear)
{
	double years1 = convert_times_D(type_t, MedTime::Years, (double)t);
	double years2 = (double)byear - 1900 + 0.5;
	return (years1 - years2);
}

//.....................................................................................................
double MedTime::get_age_from_bdate(int t, int type_t, int bdate)
{
	double years1 = convert_times_D(type_t, MedTime::Years, t);
	double years2 = convert_times_D(MedTime::Date, MedTime::Years, bdate);
	return (years1 - years2);
}

//.....................................................................................................
int MedTime::add_subtract_time(int in_time, int in_type, int delta_time, int delta_type)
{
	int conved = convert_times(in_type, delta_type, in_time);
	conved += delta_time;
	
	if (conved < 0) {
		MWARN("WARNING: Found negative time in add_subtract_time. Truncating to 0 (%d)\n", in_time);
		conved = 0;
	}

	return convert_times(delta_type, in_type, conved);
}