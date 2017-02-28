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
	else if (generator_name == "binnedLmEstimates" || generator_name == "binnedLm"  || generator_name == "binnedLM")
		return FTR_GEN_BINNED_LM;
	else
		return FTR_GEN_LAST;
}

// Initialize featurse
//.......................................................................................
void FeatureGenerator::init(MedFeatures &features) {

	if (names.size() == 0)
		set_names();
	for (auto& name : names) 
		features.attributes[name].normalized = false;
}

//.......................................................................................
FeatureGenerator *FeatureGenerator::create_generator(string &params)
{
	string fg_type;
	get_single_val_from_init_string(params, "fg_type", fg_type);
	return (make_generator(fg_type, params));
}

// Initialization
//.......................................................................................
FeatureGenerator *FeatureGenerator::make_generator(string generator_name) {

	return make_generator(ftr_generator_name_to_type(generator_name));
}

//.......................................................................................
FeatureGenerator *FeatureGenerator::make_generator(string generator_name, string init_string) {

	//MLOG("making generator %s , %s\n", generator_name.c_str(), init_string.c_str());
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

	//MLOG("making generator %d , %s\n", (int)generator_type, init_string.c_str());
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
BasicFeatureTypes BasicFeatGenerator::name_to_type(const string &name)
{

	if (name == "last")				return FTR_LAST_VALUE;
	if (name == "first")			return FTR_FIRST_VALUE;
	if (name == "last2")			return FTR_LAST2_VALUE;
	if (name == "avg")				return FTR_AVG_VALUE;
	if (name == "max")				return FTR_MAX_VALUE;
	if (name == "min")				return FTR_MIN_VALUE;
	if (name == "std")				return FTR_STD_VALUE;
	if (name == "last_delta")		return FTR_LAST_DELTA_VALUE;
	if (name == "last_time")		return FTR_LAST_DAYS;
	if (name == "last_time2")		return FTR_LAST2_DAYS;
	if (name == "slope")			return FTR_SLOPE_VALUE;
	if (name == "win_delta")				return FTR_WIN_DELTA_VALUE;
	if (name == "category_set")				return FTR_CATEGORY_SET;
	if (name == "category_set_count")		return FTR_CATEGORY_SET_COUNT;
	if (name == "category_set_sum")			return FTR_CATEGORY_SET_SUM;

	return (BasicFeatureTypes)stoi(name);
}

//.......................................................................................
void BasicFeatGenerator::set_names() {
	
	names.clear();

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
		case FTR_LAST_DAYS:				name += "last_time"; break;
		case FTR_LAST2_DAYS:			name += "last2_time"; break;
		case FTR_SLOPE_VALUE:			name += "slope"; break;
		case FTR_WIN_DELTA_VALUE:		name += "win_delta"; break;
		case FTR_CATEGORY_SET:			name += "category_set"; break;
		case FTR_CATEGORY_SET_COUNT:	name += "category_set_count"; break;
		case FTR_CATEGORY_SET_SUM:		name += "category_set_sum"; break;
		default: name += "ERROR";
		}

		name += ".win_" + std::to_string(win_from) + "_" + std::to_string(win_to);
		if (time_channel!=0 || val_channel != 0)
			name += ".t" + std::to_string(time_channel) + "v" + std::to_string(val_channel);
		names.push_back(name);
	}

	//time_unit_sig = rep.sigs.Sid2Info[sid].time_unit; !! this is an issue to SOLVE !!
}

// Generate
//.......................................................................................
int BasicFeatGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	string& name = names[0];
	if (time_unit_sig == MedTime::Undefined)	time_unit_sig = rec.my_base_rep->sigs.Sid2Info[signalId].time_unit;

	float *p_feat = &(features.data[name][index]);
	for (int i = 0; i < num; i++)
		//		features.data[name][index + i] = get_value(rec, i, features.samples[i].date);
		p_feat[i] = get_value(rec, i, med_time_converter.convert_times(features.time_unit, time_unit_win, features.samples[index+i].time));
	
	return 0;
}

//.......................................................................................
float BasicFeatGenerator::get_value(PidDynamicRec& rec, int idx, int time) {

	rec.uget(signalId, idx);


	switch (type) {
	case FTR_LAST_VALUE:	return uget_last(rec.usv, time, win_from, win_to);
	case FTR_FIRST_VALUE:	return uget_first(rec.usv, time);
	case FTR_LAST2_VALUE:	return uget_last2(rec.usv, time);
	case FTR_AVG_VALUE:		return uget_avg(rec.usv, time);
	case FTR_MAX_VALUE:		return uget_max(rec.usv, time);
	case FTR_MIN_VALUE:		return uget_min(rec.usv, time);
	case FTR_STD_VALUE:		return uget_std(rec.usv, time);
	case FTR_LAST_DELTA_VALUE:	return uget_last_delta(rec.usv, time);
	case FTR_LAST_DAYS:			return uget_last_time(rec.usv, time);
	case FTR_LAST2_DAYS:		return uget_last2_time(rec.usv, time);
	case FTR_SLOPE_VALUE:		return uget_slope(rec.usv, time);
	case FTR_WIN_DELTA_VALUE:	return uget_win_delta(rec.usv, time);
	case FTR_CATEGORY_SET:				return uget_category_set(rec, rec.usv, time);
	case FTR_CATEGORY_SET_COUNT:		return uget_category_set_count(rec, rec.usv, time);
	case FTR_CATEGORY_SET_SUM:			return uget_category_set_sum(rec, rec.usv, time);

	default:	return missing_val;
	}

	return missing_val;
}

// Init
//.......................................................................................
int BasicFeatGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "type") { type = name_to_type(entry.second); }
		else if (field == "win_from") win_from = stoi(entry.second);
		else if (field == "win_to") win_to = stoi(entry.second);
		else if (field == "d_win_from") d_win_from = stoi(entry.second);
		else if (field == "d_win_to") d_win_to = stoi(entry.second);
		else if (field == "signalName" || field == "signal") signalName = entry.second;
		else if (field == "time_unit") time_unit_win = med_time_converter.string_to_type(entry.second);
		else if (field == "time_channel") time_channel = stoi(entry.second);
		else if (field == "val_channel") val_channel = stoi(entry.second);
		else if (field == "sum_channel") sum_channel = stoi(entry.second);
		else if (field == "sets") boost::split(sets, entry.second, boost::is_any_of(","));
		else if (field != "fg_type")
				MLOG("Unknown parameter \'%s\' for BasicFeatGenerator\n", field.c_str());
	}

	// names for BasicFeatGenerator are set as a first step in the Learn call as we must have access to the MedRepository
	names.clear();

	set_names();

	req_signals.assign(1,signalName);

	return 0;
}

// (De)Serialization
//.......................................................................................
size_t BasicFeatGenerator::get_size() {

	size_t size = 0;

	size += sizeof(BasicFeatureTypes); //  BasicFeatureTypes type;
	size += 4 * sizeof(int); // win from-to d_win from-to
	size += 4 * sizeof(int); // time_unit_win, time_channel, val_channel, sum_channel

	// signalName
	size += MedSerialize::get_size(signalName);

	// sets
	size += MedSerialize::get_size(sets);

	return size;
}

extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t BasicFeatGenerator::serialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(blob + ptr, &type, sizeof(BasicFeatureTypes)); ptr += sizeof(BasicFeatureTypes);
	memcpy(blob + ptr, &win_from, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &win_to, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &d_win_from, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &d_win_to, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &time_unit_win, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &time_channel, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &val_channel, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &sum_channel, sizeof(int)); ptr += sizeof(int);

	// SignalName
	ptr += MedSerialize::serialize(blob + ptr, signalName);

	// sets
	ptr += MedSerialize::serialize(blob + ptr, sets);

	return ptr;
}

//.......................................................................................
size_t BasicFeatGenerator::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(&type, blob + ptr, sizeof(BasicFeatureTypes)); ptr += sizeof(BasicFeatureTypes);
	memcpy(&win_from, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&win_to, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&d_win_from, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&d_win_to, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&time_unit_win, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&time_channel, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&val_channel, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&sum_channel, blob + ptr, sizeof(int)); ptr += sizeof(int);

	//// SignalName
	ptr += MedSerialize::deserialize(blob + ptr, signalName);
	
	// sets
	ptr += MedSerialize::deserialize(blob + ptr, sets);

	req_signals.assign(1,signalName);
	
	names.clear();
	set_names();

	return ptr;
}


//=======================================================================================
// Age
//=======================================================================================
int AgeGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	// Sanity check
	if (byearId == -1) {
		MERR("Uninitialized byearId\n");
		return -1;
	}

	int len;
	SVal *bYearSignal = (SVal *)rec.get(byearId, len);
	if (len != 1) { MERR("id %d , got len %d for signal %d (BYEAR)...\n", rec.pid, len, byearId); }
	assert(len == 1);
	int byear = (int)(bYearSignal[0].val);

	for (int i = 0; i < num; i++)
		features.data[names[0]][index + i] = (float) (med_time_converter.convert_times(features.time_unit, MedTime::Years, features.samples[i].time) - byear);

	return 0;
}

//=======================================================================================
// Gender
//=======================================================================================
int GenderGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	// Sanity check
	if (genderId == -1) {
		MERR("Uninitialized genderId\n");
		return -1;
	}

	int len;
	SVal *genderSignal = (SVal *)rec.get(genderId, len);
	assert(len == 1);
	int gender = (int)(genderSignal[0].val);

	for (int i = 0; i < num; i++)
		features.data[names[0]][index + i] = (float) gender;

	return 0;
}



//................................................................................................................
// in all following uget funcs the relevant time window is [min_time, max_time] and time is given in time_unit_win
//................................................................................................................

void BasicFeatGenerator::get_window_in_sig_time(int _win_from, int _win_to, int _time_unit_win, int _time_unit_sig, int _win_time, int &_min_time, int &_max_time)
{
	_min_time = med_time_converter.convert_times(_time_unit_win, _time_unit_sig, _win_time -_win_to);
	_max_time = med_time_converter.convert_times(_time_unit_win, _time_unit_sig, _win_time -_win_from);
}

// get the last value in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_last(UniversalSigVec &usv, int time, int _win_from, int _win_to) 
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	for (int i=usv.len-1; i>=0; i--) {
		int itime = usv.Time(i, time_channel);
		if (itime <= max_time) {
			if (itime >= min_time)
				return usv.Val(i, val_channel);
			else
				return missing_val;
		}
	}

	return missing_val;
}

//.......................................................................................
// get the first value in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_first(UniversalSigVec &usv, int time) 
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime >= min_time) {
			if (itime > max_time)
				return missing_val;
			else
				return usv.Val(i, val_channel);
		}
	}
	return missing_val;
}

//.......................................................................................
// get the last2 value (the one before the last) in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_last2(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	for (int i=usv.len-1; i>=0; i--) {
		if (usv.Time(i, time_channel) <= max_time) {
			if (i>0 && usv.Time(i-1, time_channel) >= min_time)
				return usv.Val(i-1, val_channel);
			else
				return missing_val;
		}
	}

	return missing_val;
}

//.......................................................................................
// get the average value in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_avg(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	double sum = 0, nvals = 0;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time) {
			sum += usv.Val(i, val_channel);
			nvals++;
		}
	}

	if (nvals > 0)
		return (float)(sum/nvals);

	return missing_val;
}

//.......................................................................................
// get the max value in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_max(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	float max_val = -1e10;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time && usv.Val(i, val_channel) > max_val)	max_val = usv.Val(i, val_channel);
	}

	if (max_val > -1e10)
		return max_val;

	return missing_val;
}

//.......................................................................................
// get the min value in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_min(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	float min_val = (float)1e20;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time && usv.Val(i, val_channel) < min_val) 	min_val = usv.Val(i, val_channel);
	}

	if (min_val < (float)1e20)
		return min_val;

	return missing_val;
}

//.......................................................................................
// get the std in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_std(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	double sum = 0, sum_sq = 0, nvals = 0;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time) {
			float ival = usv.Val(i, val_channel);
			sum += ival;
			sum_sq += ival * ival;
			nvals++;
		}
	}

	if (nvals > 1) {
		double avg = sum/nvals;
		double var = sum_sq/nvals - avg*avg;
		if (var < 0.0001) var = 0.0001;
		return (float)sqrt(var);
	}

	return missing_val;
}

//.......................................................................................
float BasicFeatGenerator::uget_last_delta(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	for (int i=usv.len-1; i>=0; i--) {
		if (usv.Time(i, time_channel) <= max_time) {
			if (i>0 && usv.Time(i-1, time_channel) >= min_time)
				return (usv.Val(i, val_channel) - usv.Val(i-1, val_channel));
			else
				return missing_val;
		}
	}
	return missing_val;
}

//.......................................................................................
float BasicFeatGenerator::uget_last_time(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	for (int i=usv.len-1; i>=0; i--) {
		int itime = usv.Time(i, time_channel);
		if (itime <= max_time)
			if (itime >= min_time)
				return (float)(time - usv.TimeU(i, time_channel, time_unit_win));
			else
				return missing_val;
	}

	return missing_val;
}

//.......................................................................................
float BasicFeatGenerator::uget_last2_time(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	for (int i=usv.len-1; i>=0; i--) {
		if (usv.Time(i, time_channel) <= max_time) {
			if (i>0 && usv.Time(i-1, time_channel) >= min_time)
				return (float)(time - usv.TimeU(i-1, time_channel, time_unit_win));
			else
				return missing_val;
		}
	}

	return missing_val;
}

//.......................................................................................
// get the slope in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_slope(UniversalSigVec &usv, int time)
{

	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	double sx = 0, sy = 0, sxx = 0, sxy = 0, n = 0;
	double t_start = -1;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time) {
			if (t_start < 0) t_start = usv.TimeU(i, time_channel, time_unit_win);
			double t_curr = usv.TimeU(i, time_channel, time_unit_win);
			double x = (t_curr - t_start)/500.0;
			double ival = usv.Val(i, val_channel);
			sx += x;
			sy += ival;
			sxx += x*x;
			sxy += x*ival;
			n++;
		}
	}

	if (n < 2) return missing_val;

	double cov = sxy - sx*(sy/n);
	double var = sxx - sx*(sx/n);

	if (var < 0.1)		return 0;

	return ((float)(cov/var));

}

//.......................................................................................
float BasicFeatGenerator::uget_win_delta(UniversalSigVec &usv, int time)
{
	float val1 = uget_last(usv, time, win_from, win_to);
	if (val1 == missing_val) return missing_val;

	float val2 = uget_last(usv, time, d_win_from, d_win_to);
	if (val2 == missing_val) return missing_val;

	return (val1 - val2);
}

//.......................................................................................
float BasicFeatGenerator::uget_category_set(PidDynamicRec &rec, UniversalSigVec &usv, int time)
{
#pragma omp critical
	if (lut.size() == 0) {

		int section_id = rec.my_base_rep->dict.section_id(signalName);
		//MLOG("signalName %s section_id %d sets size %d sets[0] %s\n", signalName.c_str(), section_id, sets.size(), sets[0].c_str());
		rec.my_base_rep->dict.prep_sets_lookup_table(section_id, sets, lut);
		int n1=0; for (auto l : lut) n1 += l;
		//MLOG("size of lut %d , n1 %d\n", lut.size(), n1);
	}

	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time && lut[(int)usv.Val(i, val_channel)]) 	return 1;
	}

	return 0;
}

//.......................................................................................
float BasicFeatGenerator::uget_category_set_count(PidDynamicRec &rec, UniversalSigVec &usv, int time)
{
#pragma omp critical
	if (lut.size() == 0) {
		int section_id = rec.my_base_rep->dict.section_id(signalName);
		rec.my_base_rep->dict.prep_sets_lookup_table(section_id, sets, lut);
	}

	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	int cnt = 0;
	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time && lut[(int)usv.Val(i, val_channel)]) 	cnt++;
	}

	return (float)cnt;
}

//.......................................................................................
float BasicFeatGenerator::uget_category_set_sum(PidDynamicRec &rec, UniversalSigVec &usv, int time)
{
#pragma omp critical
	if (lut.size() == 0) {
		int section_id = rec.my_base_rep->dict.section_id(signalName);
		rec.my_base_rep->dict.prep_sets_lookup_table(section_id, sets, lut);
	}

	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	float sum = 0;
	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time && lut[(int)usv.Val(i, val_channel)]) 	sum += usv.Val(i, sum_channel);
	}

	return sum;
}
