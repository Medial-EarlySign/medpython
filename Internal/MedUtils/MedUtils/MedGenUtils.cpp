#include "MedGenUtils.h"
#include <Logger/Logger/Logger.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_MED_UTILS
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;


void set_rand_seed(int seed)
{
	if (seed != -1)
		globalRNG::srand(seed);
	else {
		MedTimer t;
		unsigned long long curr = t.get_clock_micro();
		unsigned int s = (unsigned int)(curr & 0xffffffff);
		globalRNG::srand(s);
	}
}

void get_rand_vector_no_repetitions(vector<int> &v, int N, int size)
{
	vector<int> w(N);

	for (int i = 0; i<N; i++)
		w[i] = i;

	shuffle(w.begin(), w.end(), globalRNG::get_engine());

	v.resize(size);

	for (int i = 0; i<size; i++)
		v[i] = w[i];
}

void get_rand_vector_with_repetitions(vector<int> &v, int N, int size)
{
	v.resize(size);
	for (int i = 0; i<size; i++)
		v[i] = rand_N(N);
}


void get_rand_splits(vector<int> &split, int nsplits, int size)
{
	split.resize(size);
	for (int i = 0; i<size; i++)
		split[i] = i % nsplits;
	shuffle(split.begin(), split.end(),globalRNG::get_engine());
}

void categorize_vec(vector<float> &in, vector<float> &bounds, vector<float> &out)
{
	out.resize(in.size());

	for (int i = 0; i<in.size(); i++) {
		int j = 0;
		while (j<bounds.size() && in[i]>bounds[j + 1]) j++;
		out[i] = (float)j;
	}
}

void get_probs_vec(vector<float> &v)
{
	float sum = 0;

	for (int i = 0; i<v.size(); i++) {
		sum += v[i];
	}

	if (sum > 0) {
		for (int i = 0; i<v.size(); i++)
			v[i] = v[i] / sum;
	}
}

bool is_windows_os(void)
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32)
	return true;
#else
	return false;
#endif
}

int get_day_approximate(int val) {
	if (val <= 19000101 || val >= 21000101)
		return -1;
	if (get_day(val) != -1)
		return get_day(val);
	if (get_day((val / 100) * 100 + 1) != -1)
		return get_day((val / 100) * 100 + 1);
	if (get_day((val / 10000) * 10000 + 101) != -1)
		return get_day((val / 10000) * 10000 + 101);
	fprintf(stderr, "unable to convert %d to days \n", val);
	throw exception();
}

int get_day(int val) {
	int year = val / 100 / 100;
	int month = (val / 100) % 100;
	int day = val % 100;

	if (month < 1 || month > 12 || year < 1900 || day < 0)
		return -1;

	// Full years
	int days = 365 * (year - 1900);
	days += (year - 1897) / 4;
	days -= (year - 1801) / 100;
	days += (year - 1601) / 400;

	// Full Months

	int days2month[] = { 0,31,59,90,120,151,181,212,243,273,304,334,365 };
	days += days2month[month - 1];
	if (month>2 && (year % 4) == 0 && ((year % 100) != 0 || (year % 400) == 0))
		days++;
	days += (day - 1);

	return days;
}

// Days -> Date
int get_date(int days) {

	int days2month[] = { 0,31,59,90,120,151,181,212,243,273,304,334,365 };
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
	int month = 0;
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
