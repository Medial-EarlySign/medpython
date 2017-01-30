#define _CRT_SECURE_NO_WARNINGS

#include "FeatureGenerator.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// FeatGenerator
//=======================================================================================
// Generator types
FeatureGeneratorTypes ftr_generator_name_to_type(const string& generator_name) {

	if (generator_name == "basic")
		return FTR_GEN_BASIC;
	else if (generator_name == "age")
		return FTR_GEN_AGE;
	else if (generator_name == "gender")
		return FTR_GEN_GENDER;
	else if (generator_name == "binnedLmEstimates")
		return FTR_GEN_BINNED_LM;
	else
		return FTR_GEN_LAST;
}

// Initialize featurse
//.......................................................................................
void FeatureGenerator::init(MedFeatures &features) {

	for (auto& name : names)
		features.attributes[name].normalized = false;
}

// Initialization
//.......................................................................................
FeatureGenerator *FeatureGenerator::make_generator(string generator_name) {

	return make_generator(ftr_generator_name_to_type(generator_name));
}

//.......................................................................................
FeatureGenerator *FeatureGenerator::make_generator(string generator_name, string init_string) {

	return make_generator(ftr_generator_name_to_type(generator_name), init_string);
}

//.......................................................................................
FeatureGenerator *FeatureGenerator::make_generator(FeatureGeneratorTypes generator_type) {

	if (generator_type == FTR_GEN_BASIC)
		return new BasicFeatGenerator;
	else if (generator_type == FTR_GEN_AGE)
		return new AgeGenerator;
	else if (generator_type == FTR_GEN_GENDER)
		return new GenderGenerator;
	else if (generator_type == FTR_GEN_BINNED_LM)
		return new BinnedLmEstimates;
	else
		return NULL;

}

//.......................................................................................
FeatureGenerator * FeatureGenerator::make_generator(FeatureGeneratorTypes generator_type, string init_string) {

	FeatureGenerator *newFtrGenerator = make_generator(generator_type);
	newFtrGenerator->init_from_string(init_string);
	return newFtrGenerator;
}

//.......................................................................................
// Add at end of feature vector
int FeatureGenerator::generate(PidDynamicRec& in_rep, MedFeatures& features) {

	return Generate(in_rep, features, features.get_pid_pos(in_rep.pid), features.get_pid_len(in_rep.pid));

}

//.......................................................................................
// Add uncleaned data at end of feature vector
int FeatureGenerator::generate(MedPidRepository& rep, int id, MedFeatures& features) {

	int samples_size = (int)features.samples.size();
	int data_size = (int)features.data[names[0]].size();

	if (data_size > samples_size) {
		MERR("Data (%d) is longer than Samples (%d) for %s. Cannot generate features \n", data_size, samples_size, names[0].c_str());
		return -1;
	}

	features.data[names[0]].resize(samples_size);
	return generate(rep, id, features, data_size, (int)(samples_size - data_size));
}

//.......................................................................................
// Add uncleaned data
int FeatureGenerator::generate(MedPidRepository& rep, int id, MedFeatures& features, int index, int num) {

	PidDynamicRec rec;
	rec.prealloc(DYNAMIC_REC_SIZE);
	
	if (req_signal_ids.empty())
		get_required_signal_ids(rec.my_base_rep->dict);

	rec.init_from_rep(std::addressof(rep), id, req_signal_ids, num);

	return Generate(rec, features, index, num);
}

// (De)Serialize
//.......................................................................................
size_t FeatureGenerator::get_generator_size() {
	return sizeof(generator_type) + get_size();
}

//.......................................................................................
size_t FeatureGenerator::generator_serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &generator_type, sizeof(FeatureGeneratorTypes)); ptr += sizeof(FeatureGeneratorTypes);
	ptr += serialize(blob + ptr);

	return ptr;
}

// Required signals
//.......................................................................................
void FeatureGenerator::get_required_signal_ids(MedDictionarySections& dict){

	req_signal_ids.resize(req_signals.size());

	for (unsigned int i = 0; i < req_signals.size(); i++)
		req_signal_ids[i] = dict.id(req_signals[i]);
}

//.......................................................................................
void FeatureGenerator::get_required_signal_ids(unordered_set<int>& signalIds, MedDictionarySections& dict) {

	if (req_signal_ids.empty())
		get_required_signal_ids(dict);

	for (int signalId : req_signal_ids)
		signalIds.insert(signalId);
}

//=======================================================================================
// Single signal features that do not require learning(e.g. last hemoglobin)
//=======================================================================================
// 
//.......................................................................................
void BasicFeatGenerator::set_names() {
	
	if (names.empty()) {
		string name = signalName + ".";
		switch (type) {
		case FTR_LAST_VALUE:	name += "last"; break;
		case FTR_FIRST_VALUE:	name += "first"; break;
		case FTR_LAST2_VALUE:	name += "last2"; break;
		case FTR_AVG_VALUE:		name += "avg"; break;
		case FTR_MAX_VALUE:		name += "max"; break;
		case FTR_MIN_VALUE:		name += "min"; break;
		case FTR_STD_VALUE:		name += "std"; break;
		case FTR_LAST_DELTA_VALUE:		name += "last_delta"; break;
		case FTR_LAST_DAYS:		name += "last_days"; break;
		case FTR_LAST2_DAYS:		name += "last2_days"; break;
		default: name += "ERROR";
		}

		name += ".win_" + std::to_string(win_from) + "_" + std::to_string(win_to);
		names.push_back(name);
	}
}

// Generate
//.......................................................................................
int BasicFeatGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	string& name = names[0];

	float *p_feat = &(features.data[name][index]);
	for (int i = 0; i < num; i++)
//		features.data[name][index + i] = get_value(rec, i, features.samples[i].date);
		p_feat[i] = get_value(rec, i, features.samples[i].date);

	return 0;
}

//.......................................................................................
float BasicFeatGenerator::get_value(PidDynamicRec& rec, int idx, int date) {

	// signalId
	if (signalId == -1)
		signalId = rec.my_base_rep->dict.id(signalName);

	int len;
	SDateVal *signal = (SDateVal *)rec.get(signalId, idx, len);

	switch (type) {
		case FTR_LAST_VALUE:	return get_last_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_FIRST_VALUE:	return get_first_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_LAST2_VALUE:	return get_last2_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_AVG_VALUE:		return get_avg_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_MAX_VALUE:		return get_max_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_MIN_VALUE:		return get_min_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_STD_VALUE:		return get_std_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_LAST_DELTA_VALUE:		return get_last_delta_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_LAST_DAYS:			return get_last_days_value(signal, len, date, win_from, win_to, missing_val);
		case FTR_LAST2_DAYS:		return get_last2_days_value(signal, len, date, win_from, win_to, missing_val);

		default:	return missing_val;
	}

	return missing_val;
}

// Init
//.......................................................................................
int BasicFeatGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "type") type = (BasicFeatureTypes) stoi(entry.second);
		else if (field == "win_from") win_from = stoi(entry.second);
		else if (field == "win_to") win_to = stoi(entry.second);
		else if (field == "signalName") signalName = entry.second;
		else MLOG("Unknonw parameter \'%s\' for FeatureNormalizer\n", field.c_str());
	}

	names.clear();
	set_names();
	return 0;
}

// (De)Serialization
//.......................................................................................
size_t BasicFeatGenerator::get_size() {

	size_t size = 0;

	size += sizeof(BasicFeatureTypes); //  BasicFeatureTypes type;
	size += 2 * sizeof(int); // win from-to

	// signalName
	size += sizeof(size_t); 
	size += signalName.length()+1 ;

	return size;
}

extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t BasicFeatGenerator::serialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(blob + ptr, &type, sizeof(BasicFeatureTypes)); ptr += sizeof(BasicFeatureTypes);
	memcpy(blob + ptr, &win_from, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &win_to, sizeof(int)); ptr += sizeof(int);

	// SignalName
	size_t nameLen = signalName.length();
	assert(nameLen < MAX_NAME_LEN);

	strcpy(signalName_c, signalName.c_str());

	memcpy(blob + ptr, &nameLen, sizeof(size_t)); ptr += sizeof(size_t);
	memcpy(blob + ptr, signalName_c, nameLen + 1); ptr += nameLen + 1;

	return ptr;
}

//.......................................................................................
size_t BasicFeatGenerator::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(&type, blob + ptr, sizeof(BasicFeatureTypes)); ptr += sizeof(BasicFeatureTypes);
	memcpy(&win_from, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&win_to, blob + ptr, sizeof(int)); ptr += sizeof(int);

	// SignalName
	size_t nameLen;
	memcpy(&nameLen, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);
	assert(nameLen < MAX_NAME_LEN);

	memcpy(signalName_c, blob + ptr,nameLen+1); ptr += nameLen + 1;
	signalName = signalName_c;

	set_names();

	return ptr;
}


//=======================================================================================
// Age
//=======================================================================================
int AgeGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	// signalId
	if (byearId == -1)
		byearId = rec.my_base_rep->dict.id("BYEAR");

	int len;
	SVal *bYearSignal = (SVal *)rec.get(byearId, len);
	if (len != 1) { MERR("id %d , got len %d for signal %d (BYEAR)...\n", rec.pid, len, byearId); }
	assert(len == 1);
	int byear = (int)(bYearSignal[0].val);

	for (int i = 0; i < num; i++)
		features.data[names[0]][index + i] = (float) (features.samples[i].date / 10000 - byear); 

	return 0;
}

//=======================================================================================
// Gender
//=======================================================================================
int GenderGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	// signalId
	if (genderId == -1)
		genderId = rec.my_base_rep->dict.id("GENDER");

	int len;
	SVal *genderSignal = (SVal *)rec.get(genderId, len);
	assert(len == 1);
	int gender = (int)(genderSignal[0].val);

	for (int i = 0; i < num; i++)
		features.data[names[0]][index + i] = (float) gender;

	return 0;
}

//=======================================================================================
// Utilities
//=======================================================================================
// in all the following timing is: .... -> (date-win_to) -> .... -> (date-win_from) -> ... -> date -> ...

//.......................................................................................
// get the last value in the window [win_to, win_from) before date
float get_last_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val) {

	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	float val = missing_val;
	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) {
			if (i>0 && signal[i - 1].date >= min_date)
			{
				val = signal[i - 1].val; break;
			}
			else {
				val = missing_val; break;
			}
		}
	}
//	MLOG("len %d date %d win_from %d win_to %d missing %f val %f min_date %d max_date %d\n", len, date, win_from, win_to, missing_val, val, min_date, max_date);
//		return missing_val;
	return val;
}

//.......................................................................................
// get the first value in the window [win_to, win_from) before date
float get_first_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val) {

	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= min_date) {
			if (signal[i].date >= max_date)
				return missing_val;
			else
				return signal[i].val;
		}
	}
	return missing_val;
}

//.......................................................................................
// get the last2 value (the one before the last) in the window [win_to, win_from) before date
float get_last2_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val)
{
	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) {
			if (i>1 && signal[i - 2].date >= min_date)
				return signal[i - 2].val;
			else
				return missing_val;
		}
	}

	return missing_val;
}

//.......................................................................................
// get the average value in the window [win_to, win_from) before date
float get_avg_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val)
{
	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	double sum = 0, nvals = 0;

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) break;
		if (signal[i].date >= min_date) {
			sum += signal[i].val;
			nvals++;
		}
	}

	if (nvals > 0)
		return (float)(sum/nvals);

	return missing_val;
}

//.......................................................................................
// get the max value in the window [win_to, win_from) before date
float get_max_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val)
{
	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	float max_val = -1e10;

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) break;
		if (signal[i].date >= min_date && signal[i].val > max_val)	max_val = signal[i].val;
	}

	if (max_val > -1e10)
		return max_val;

	return missing_val;
}

//.......................................................................................
float get_min_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val)
// get the min value in the window [win_to, win_from) before date
{
	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	float min_val = 1e10;

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) break;
		if (signal[i].date >= min_date && signal[i].val < min_val) 	min_val = signal[i].val;
	}

	if (min_val < 1e10)
		return min_val;

	return missing_val;
}

//.......................................................................................
// get the std in the window [win_to, win_from) before date
float get_std_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val)
{
	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	double sum = 0, sum_sq = 0, nvals = 0;

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) break;
		if (signal[i].date >= min_date) {
			sum += signal[i].val;
			sum_sq += signal[i].val * signal[i].val;
			nvals++;
		}
	}

	if (nvals > 1) {
		double avg = sum/nvals;
		double var = sum_sq/nvals - avg*avg;
		return (float)sqrt(var);
	}

	return missing_val;
}

//.......................................................................................
float get_last_delta_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val)
{
	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) {
			if (i>0 && signal[i - 1].date >= min_date)
				return (signal[i].val - signal[i - 1].val);
			else
				return missing_val;
		}
	}
	return missing_val;
}

//.......................................................................................
float get_last_days_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val)
{
	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) {
			if (i>0 && signal[i - 1].date >= min_date)
				return (float)(get_day(date) - get_day(signal[i-1].date));
			else
				return missing_val;
		}
	}
	return missing_val;
}

//.......................................................................................
float get_last2_days_value(SDateVal *signal, int len, int date, int win_from, int win_to, float missing_val)
{
	int min_date = get_date(get_day(date) - win_to);
	int max_date = get_date(get_day(date) - win_from);

	for (int i = 0; i < len; i++) {
		if (signal[i].date >= max_date) {
			if (i>1 && signal[i - 2].date >= min_date)
				return (float)(get_day(date) - get_day(signal[i-2].date));
			else
				return missing_val;
		}
	}

	return missing_val;
}