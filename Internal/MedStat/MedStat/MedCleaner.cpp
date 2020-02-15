#ifndef __MED_CLEANER_CPP__
#define __MED_CLEANER_CPP__

#include "MedCleaner.h"
#include "MedStat.h"
#include "MedUtils/MedUtils/MedUtils.h"


#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

// Constructor
MedCleaner::MedCleaner() {
	missing_value = MED_DEFAULT_MISSING_VALUE; remove_flag = trim_flag = normalize_flag = replace_missing_to_mean_flag = true; n = 0; nvals = 0;
	most_common_count = 0; nzeros = 0; median = q1 = q3 = iqr = MED_DEFAULT_MISSING_VALUE; median = MED_DEFAULT_MISSING_VALUE;mean = MED_DEFAULT_MISSING_VALUE; sdv = 0; skew = 0;
	most_common_value = MED_DEFAULT_MISSING_VALUE; trim_min = MED_DEFAULT_MISSING_VALUE;  min_trim = MED_DEFAULT_MIN_TRIM;
	trim_max = MED_DEFAULT_MISSING_VALUE; remove_min = MED_DEFAULT_MISSING_VALUE; remove_max = -1; sk = 0; skew_sign = 0;
}

// Return the trimmed value, or 'missing-value' if missing-value/invalid (outside removal limits)
float MedCleaner::get_value(float value) {

	if (value == missing_value || ! is_valid(value))
		return missing_value ;
	else
		return get_trimmed(value) ;

}

// Check if value is valid : 'missing-value' or value withing removal limits
bool MedCleaner::is_valid(float value) {
	return (value == missing_value || (value <= remove_max && value >= remove_min)) ;
}


// Trim value in-place : missing-value, or value within trimming limits
int MedCleaner::trim(float& value) {

	if (value == missing_value)
		return 0 ;
	else if (value > trim_max) {
		value = trim_max ;
		return 1 ;
	} else if (value < trim_min) {
		value = trim_min;
		return 1;
	} else 
		return 0 ;
}

// Return trimmed value : missing value, or value within trimming limits
float MedCleaner::get_trimmed(float value) {
	if (value == missing_value)
		return value ;
	else if (value > trim_max)
		return trim_max ;
	else if (value < trim_min)
		return trim_min;
	else
		return value ;
}

// Clear data - remove and trim
int MedCleaner::clear(vector<float>& values) {

	int changed = 0 ;

	for (unsigned int i=0; i<values.size(); i++) {
		if (remove_flag && ! is_valid(values[i])) {
			values[i] = missing_value ;
			changed ++ ;
		}
		if (trim_flag && values[i] != missing_value)
			changed += trim(values[i]) ;
	}

	return changed ;
}

void MedCleaner::single_remove_trim_replace(float &val)
{
	if (remove_flag && !is_valid(val))
		val = missing_value;
	if (replace_missing_to_mean_flag && val == missing_value)
		val = mean;
}

void MedCleaner::remove_trim_replace(vector<float> &values)
{
	for (int i=0; i<values.size(); i++) {
		single_remove_trim_replace(values[i]);
	}
}

void MedCleaner::single_normalize(float &val) 
{
	if (replace_missing_to_mean_flag && val == missing_value)
		val = (normalize_flag ? 0 : mean);
	if (normalize_flag && val != missing_value)
		val = (val - mean)/sdv ;
}

// Normalize data according to mean + sdv
void MedCleaner::normalize(vector<float>& values) {
	
	for (unsigned int i=0; i<values.size(); i++) {
		single_normalize(values[i]);
	}
}


// Calculate trimming min/max iteratively - first version
// std_mult - allow values within std_mult std's of mean. should be positive.
// iteratively trim and remove outliers
void MedCleaner::get_limits_iteratively(vector<float> values, float std_mult) {
	
	int changed = 1 ;
	//int iter = 0;
	
	float min_val,max_val ;
	
	//MLOG("get_limits_iteratively iter %d\n", iter++);
	get_cleaning_range(values,min_val,max_val,std_mult) ;
	trim_min = min_val ; remove_min = trim_min - 1 ;
	trim_max = max_val ; remove_max = trim_max + 1 ;

	vector<float> orig_values = values;
	while (changed) {
		changed = clear(values);
		if (changed) {
			//MLOG("get_limits_iteratively iter %d\n", iter++);
			get_cleaning_range(values,min_val,max_val,std_mult) ;
			trim_min = min_val ; remove_min = trim_min - 1 ;
			trim_max = max_val ; remove_max = trim_max + 1 ;
		}

	}
	// once we got our final limits, recalc sdv and take into account missing values
	get_mean_and_sdv(values, true);
}

// Calculate trimming min/max iteratively - second version
// Trim min/max are minimal/maximal values within remove min/max range
void MedCleaner::get_cleaning_params(vector<float> values) {

	int changed = 1 ;
	int total_change = 0 ;
	
	float min_val,max_val ;
	
	// Remove min/max
	get_cleaning_range(values,min_val,max_val) ;
	trim_min = remove_min = min_val ;
	trim_max = remove_max = max_val ;

	while (changed) {
		changed = clear(values) ;

		if (changed) {
			total_change += changed ;

			get_cleaning_range(values,min_val,max_val) ;
			trim_min = remove_min = min_val ;
			trim_max = remove_max = max_val ;
			//fprintf(stderr, "min_val %f trim_min %f remove_min %f \n", min_val, trim_min, remove_min);
			//fprintf(stderr, "max_val %f trim_max %f remove_max %f \n", max_val, trim_max, remove_max);

		}
	}

	// after removing all outliers, we set trim to the min/max in the remaining values
	int init = 1 ;
	for (unsigned int i=0; i<values.size(); i++) {
		if (values[i] != missing_value) {
			if (init == 1) {
				init = 0 ;
				trim_min = trim_max = values[i] ;
			} else if (values[i] > trim_max) {
				trim_max = values[i] ;
			} else if (values[i] < trim_min) {
				trim_min = values[i] ;
			}
		}
	}
}

float calc_quartile(map<float, int> counter, int q) {
	int sub_count = 0;
	for (auto it = counter.begin(); it != counter.end(); it++) {
		sub_count += it->second;
		if (sub_count >= q)
			return it->first;		
	}
	throw exception();
}

// Calculate all moments, ignores missing values
void MedCleaner::calculate(vector<float> &_values) {

	vector<float> values ;
	for (unsigned int i=0; i<_values.size(); i++) {
		if (_values[i] != missing_value)
			values.push_back(_values[i]) ;
	}

	// Size
	n = (int) values.size() ;

	// Number of different values
	map<float,int> counter ;
	for (int i=0; i<n; i++)
		counter[values[i]]++ ;

	nvals = (int) counter.size() ;
	nzeros = counter[0] ;

	// Most common value
	most_common_count= 0 ;
	for (auto it = counter.begin(); it != counter.end(); it++) {
		if (it->second > most_common_count) {
			most_common_count = it->second ;
			most_common_value = it->first ;
		}
	}

	// Median
	median = calc_quartile(counter, n / 2);
	q1 = calc_quartile(counter, n / 4);
	q3 = calc_quartile(counter, n*3 / 4);
	max = calc_quartile(counter, n);
	iqr = q3 - q1;


	get_mean_and_sdv(values);

	// MLOG("cleaner : avg %f std %f\n", mean, sdv);
	// Skew
	float s = 0 ;
	for (auto it = counter.begin(); it != counter.end(); it++) {
		float v = (it->first-mean)/sdv ;
		s += v*v*v ;
	}
	skew = (sdv != 0) ? s/n : (float) -99999.99;

	skew_sign = (skew > 0) ? 1 : -1 ;
	sk = ((skew > 0) ? 1 : -1) * exp(log(fabs(skew)/3)) ;

	// Min/Max
	trim_max = mean + MED_CLEANER_MAX_Z * sdv ;
	trim_min = mean - MED_CLEANER_MAX_Z * sdv ;
	if (trim_min < min_trim) trim_min = min_trim;

	remove_min = trim_min - 1 ;
	remove_max = trim_max - 1 ;

}

// Calculate mean + sdv
// take_missing_into_account: sdv should be adjusted if there are lots of missing values that are to be replaced with mean
void MedCleaner::get_mean_and_sdv(vector<float>& values, bool take_missing_into_account) {

	// Mean
	double s = 0, s2 = 0  ;
	double n = 0 ;

	int size = (int)values.size();

	for (auto val : values)
		if (val != missing_value) {
			n++;
			s += val;
		}

	double mean_without_missing;
	if (n > 0)
		mean_without_missing = s/n;
	else
		mean_without_missing = 0;
	
	if (take_missing_into_account) {
		s += (size-n)*mean_without_missing;
		n = size;
	}

	if (n > 0)
		mean = (float)(s/n);
	else
		mean = 0;


	for (auto val : values)
		if (val != missing_value) 
			s2 += (double)(val - mean) * (val - mean);

	if (n > 0)
		sdv = (float)(sqrt(s2/n));
	else
		sdv = 1;

	if (sdv < MED_CLEANER_EPSILON) sdv = 1;

	//MLOG("cleaner : size %d mean %f std %f n %d adjust %d\n", size, mean, sdv, (int)n, (int)take_missing_into_account);

	return;
}

// Calculate Min/Max
// Note - std_mult must be positive
void MedCleaner::get_cleaning_range (vector<float>& values, float& min_val, float& max_val, float std_mult) {
	get_mean_and_sdv(values) ;
	min_val = mean - std_mult * sdv;
	max_val = mean + std_mult * sdv;
	return ;
}

// Print
void MedCleaner::print(const string& prefix) {

	MOUT("%s : N = %d\n",prefix.c_str(),n) ;
	MOUT("%s : N-Vals = %d\n",prefix.c_str(),nvals) ;
	MOUT("%s : No. of zeros = %d\n",prefix.c_str(),nzeros) ;
	MOUT("%s : Most common Value %f appears %d times\n",prefix.c_str(),most_common_value,most_common_count) ;
	MOUT("%s : Median = %f\n",prefix.c_str(),median) ;
	MOUT("%s : Mean = %f\n",prefix.c_str(),mean) ;
	MOUT("%s : Standard Deviation = %f\n",prefix.c_str(),sdv) ;
	MOUT("%s : Skewness = %f\n",prefix.c_str(),skew) ;
	if (trim_flag)
		MOUT("%s : Trimming range = [%f,%f]\n",prefix.c_str(),trim_min,trim_max) ;
	if (remove_flag)
		MOUT("%s : Removing range = [%f,%f]\n",prefix.c_str(),remove_min,remove_max) ;
}


void MedCleaner::print_short(const string& prefix) {

//	MOUT("%s Median %6.2f Mean %6.2f Std %6.2f N %8d\n",prefix.c_str(),median,mean,sdv,n) ;
	MOUT("%s N %d non-zeros %d Median %6.2f Mean %6.2f Std %6.2f Max %6.2f\n",prefix.c_str(),n, n-nzeros, median, mean, sdv, max) ;
}

// Get size of object
size_t MedCleaner::get_size() {
	 return  (11 * sizeof(float) + 5 * sizeof(int) + 2 * sizeof(bool))  ;
}

// Serialize
size_t MedCleaner::serialize(unsigned char *buffer) {

	size_t size = get_size() ;	

	size_t new_size = 0 ;
	memcpy(buffer+new_size,&missing_value,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&median,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&mean,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&sdv,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&skew,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&most_common_count,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&trim_min,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&trim_max,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&remove_min,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&remove_max,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(buffer+new_size,&sk,sizeof(float)) ; new_size += sizeof(float) ;

	memcpy(buffer+new_size,&n,sizeof(int)) ; new_size += sizeof(int) ;
	memcpy(buffer+new_size,&nvals,sizeof(int)) ; new_size += sizeof(int) ;
	memcpy(buffer+new_size,&most_common_count,sizeof(int)) ; new_size += sizeof(int) ;
	memcpy(buffer+new_size,&nzeros,sizeof(int)) ; new_size += sizeof(int) ;
	memcpy(buffer+new_size,&skew_sign,sizeof(int)) ; new_size += sizeof(int) ;

	memcpy(buffer+new_size,&trim_flag,sizeof(bool)) ; new_size += sizeof(bool) ;
	memcpy(buffer+new_size,&remove_flag,sizeof(bool)) ; new_size += sizeof(bool) ;

	assert(new_size == size) ;

	return size ;
}

// Deserialize
size_t MedCleaner::deserialize(unsigned char *buffer) {

	int new_size = 0 ;
	memcpy(&missing_value,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&median,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&mean,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&sdv,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&skew,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&most_common_count,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&trim_min,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&trim_max,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&remove_min,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&remove_max,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;
	memcpy(&sk,buffer+new_size,sizeof(float)) ; new_size += sizeof(float) ;

	memcpy(&n,buffer+new_size,sizeof(int)) ; new_size += sizeof(int) ;
	memcpy(&nvals,buffer+new_size,sizeof(int)) ; new_size += sizeof(int) ;
	memcpy(&most_common_count,buffer+new_size,sizeof(int)) ; new_size += sizeof(int) ;
	memcpy(&nzeros,buffer+new_size,sizeof(int)) ; new_size += sizeof(int) ;
	memcpy(&skew_sign,buffer+new_size,sizeof(int)) ; new_size += sizeof(int) ;

	memcpy(&trim_flag,buffer+new_size,sizeof(bool)) ; new_size += sizeof(bool) ;
	memcpy(&remove_flag,buffer+new_size,sizeof(bool)) ; new_size += sizeof(bool) ;

	assert(new_size == get_size()) ;
	return new_size ;
}

string MedCleaner::object_json() const {
	return "MedValueCleaner - not implemented yet";
}

#endif
