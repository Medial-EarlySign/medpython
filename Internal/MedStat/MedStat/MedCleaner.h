#ifndef _MED_CLEANER_H_
#define _MED_CLEANER_H_

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include <assert.h>
#include <math.h>

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <SerializableObject/SerializableObject/SerializableObject.h>

//======================================================================================
// MedCleaner - class to handle cleaning/normalization of data
//======================================================================================

#define MED_CLEANER_MAX_Z 15
#define MED_CLEANER_EPSILON 0.0001

// Cleaner class : Normalizing and cleaning of outliers
class MedCleaner : SerializableObject {
public:

	float missing_value;
	float min_trim;

	int n, nvals, most_common_count, nzeros;
	float median, q1, q3, iqr, mean, sdv, skew, min, max;
	float most_common_value;

	bool trim_flag, remove_flag, normalize_flag, replace_missing_to_mean_flag;
	float trim_min, trim_max;
	float remove_min, remove_max;

	float sk;
	int skew_sign;

	MedCleaner();

	void print(const string& prefix);
	void print_short(const string& prefix);
	void calculate(vector<float> &values);
	void get_mean_and_sdv(vector<float> &values, bool take_missing_into_account = false);
	void get_cleaning_range(vector<float>& values, float& min_val, float& max_val, float std_mult = MED_CLEANER_MAX_Z);
	void get_limits_iteratively(vector<float> values, float std_mult = MED_CLEANER_MAX_Z);
	void get_cleaning_params(vector<float> values);
	int clear(vector<float>& values);
	int clean(vector<float>& values) { return clear(values); };
	void remove_trim_replace(vector<float> &values);

	void normalize(vector<float>& values);

	bool is_valid(float value);
	float get_trimmed(float value);
	float get_value(float value);
	int trim(float& value);
	void single_remove_trim_replace(float &val);
	void single_normalize(float &val);

	ADD_CLASS_NAME(MedCleaner)
		size_t get_size();
	size_t serialize(unsigned char *buffer);
	size_t deserialize(unsigned char *buffer);
	string object_json() const;
};

MEDSERIALIZE_SUPPORT(MedCleaner)


#endif