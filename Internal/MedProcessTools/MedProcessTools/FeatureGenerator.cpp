#define _CRT_SECURE_NO_WARNINGS

#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>

#include "FeatureGenerator.h"
#include "SmokingGenerator.h"
#include "KpSmokingGenerator.h"
#include "DrugIntakeGenerator.h"
#include "AlcoholGenerator.h"

#include "MedProcessTools/MedProcessTools/MedModel.h"

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
	else if (generator_name == "singleton")
		return FTR_GEN_SINGLETON;
	else if (generator_name == "binnedLmEstimates" || generator_name == "binnedLm" || generator_name == "binnedLM")
		return FTR_GEN_BINNED_LM;
	else if (generator_name == "smoking")
		return FTR_GEN_SMOKING;
	else if (generator_name == "kp_smoking")
		return FTR_GEN_KP_SMOKING;
	else if (generator_name == "alcohol")
		return FTR_GEN_ALCOHOL;
	else if (generator_name == "range")
		return FTR_GEN_RANGE;
	else if (generator_name == "drugIntake")
		return FTR_GEN_DRG_INTAKE;
	else if (generator_name == "model")
		return FTR_GEN_MODEL;
	else MTHROW_AND_ERR("unknown generator name [%s]", generator_name.c_str());
}

// Prepare for feature Generation
//.......................................................................................
void FeatureGenerator::prepare(MedFeatures &features, MedPidRepository& rep, MedSamples& samples) {

	if (!iGenerateWeights) {
		//MLOG("FeatureGenerator::init _features\n");
		if (names.size() == 0)
			set_names();

		// Attributes and data
		for (auto& name : names) {
			features.attributes[name].normalized = false;
			features.data[name].resize(0, 0);
		}

		// Tags
		for (auto& name : names) {
			for (string& tag : tags)
				features.tags[name].insert(tag);
		}
	}
	else
		features.weights.resize(0, 0);
}

// Get pointers to data vectors
//.......................................................................................
void FeatureGenerator::get_p_data(MedFeatures &features) {

	p_data.clear();
	if (iGenerateWeights)
		p_data.push_back(&(features.weights[0]));
	else {
		for (string& name : names)
			p_data.push_back(&(features.data[name][0]));
	}

	return;
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
	else if (generator_type == FTR_GEN_SINGLETON)
		return new SingletonGenerator;
	else if (generator_type == FTR_GEN_BINNED_LM)
		return new BinnedLmEstimates;
	else if (generator_type == FTR_GEN_SMOKING)
		return new SmokingGenerator;
	else if (generator_type == FTR_GEN_KP_SMOKING)
		return new KpSmokingGenerator;
	else if (generator_type == FTR_GEN_ALCOHOL)
		return new AlcoholGenerator;
	else if (generator_type == FTR_GEN_RANGE)
		return new RangeFeatGenerator;
	else if (generator_type == FTR_GEN_DRG_INTAKE)
		return new DrugIntakeGenerator;
	else if (generator_type == FTR_GEN_MODEL)
		return new ModelFeatGenerator;

	else MTHROW_AND_ERR("dont know how to make_generator for [%s]", to_string(generator_type).c_str());
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
	//MLOG("gen [%s]\n", this->names[0].c_str());
	return _generate(in_rep, features, features.get_pid_pos(in_rep.pid), features.get_pid_len(in_rep.pid));
}

//.......................................................................................
// Add uncleaned data at end of feature vector
int FeatureGenerator::generate(MedPidRepository& rep, int id, MedFeatures& features) {

	int samples_size = (int)features.samples.size();
	int data_size;

	if (iGenerateWeights) {
		data_size = (int)features.weights.size();

		if (data_size > samples_size) {
			MERR("Data (%d) is longer than Samples (%d) for %s. Cannot generate weights \n", data_size, samples_size);
			return -1;
		}
		features.weights.resize(samples_size);
		p_data[0] = &(features.weights[0]);
	}
	else {
		data_size = (int)features.data[names[0]].size();

		if (data_size > samples_size) {
			MERR("Data (%d) is longer than Samples (%d) for %s. Cannot generate feature \n", data_size, samples_size, names[0].c_str());
			return -1;
		}
		for (string& name : names)
			features.data[name].resize(samples_size);

		get_p_data(features);
	}
	return generate(rep, id, features, data_size, (int)(samples_size - data_size));
}

//.......................................................................................
// Add uncleaned data
int FeatureGenerator::generate(MedPidRepository& rep, int id, MedFeatures& features, int index, int num) {

	PidDynamicRec rec;
	rec.prealloc(DYNAMIC_REC_SIZE);

	rec.init_from_rep(std::addressof(rep), id, req_signal_ids, num);

	return _generate(rec, features, index, num);
}

//.......................................................................................
// Init 
int FeatureGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "weights_generator")
			iGenerateWeights = med_stoi(entry.second);
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for FeatureGenerator\n", field.c_str());
	}
	set_names();
	return 0;
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

// Set required signal ids
//.......................................................................................
void FeatureGenerator::set_required_signal_ids(MedDictionarySections& dict) {
	if (req_signals.empty()) {
		dprint("", 1);
		MTHROW_AND_ERR("FeatureGenerator::set_required_signal_ids got empty req_signals\n");
	}
	req_signal_ids.resize(req_signals.size());

	for (unsigned int i = 0; i < req_signals.size(); i++)
		req_signal_ids[i] = dict.id(req_signals[i]);
}


// Get Required Signals
//.......................................................................................
void FeatureGenerator::get_required_signal_names(unordered_set<string>& signalNames) {
	for (auto sig : req_signals)
		signalNames.insert(sig);
}

//.......................................................................................
void FeatureGenerator::get_required_signal_ids(unordered_set<int>& signalIds) {
	for (auto sig : req_signal_ids)
		signalIds.insert(sig);
}


//.......................................................................................
// Filter generated features according to a set. return number of valid features (does not affect single-feature genertors, just returns 1/0 if feature name in set)
int FeatureGenerator::filter_features(unordered_set<string>& validFeatures) {

	int idx = 0;
	for (int i = 0; i < names.size(); i++) {
		if (validFeatures.find(names[i]) != validFeatures.end())
			names[idx++] = names[i];
	}

	names.resize(idx);

	return ((int)names.size());
}

inline bool isInteger(const std::string & s)
{
	if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;

	char * p;
	strtol(s.c_str(), &p, 10);

	return (*p == 0);
}


//.......................................................................................
void FeatureGenerator::dprint(const string &pref, int fg_flag)
{
	if (fg_flag > 0) {
		MLOG("%s :: FG type %d : serial_id %d : ", pref.c_str(), generator_type, serial_id);
		MLOG("names(%d) : ", names.size());
		if (fg_flag > 1) for (auto &name : names) MLOG("%s,", name.c_str());
		MLOG(" tags(%d) : ", tags.size());
		if (fg_flag > 1) for (auto &t : tags) MLOG("%s,", t.c_str());
		MLOG(" req_signals(%d) : ", req_signals.size());
		if (fg_flag > 1) for (auto &rsig : req_signals) MLOG("%s,", rsig.c_str());
		MLOG("\n");
	}
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
	if (name == "last2_time")		return FTR_LAST2_DAYS;
	if (name == "slope")			return FTR_SLOPE_VALUE;
	if (name == "win_delta")				return FTR_WIN_DELTA_VALUE;
	if (name == "category_set")				return FTR_CATEGORY_SET;
	if (name == "category_set_count")		return FTR_CATEGORY_SET_COUNT;
	if (name == "category_set_sum")			return FTR_CATEGORY_SET_SUM;
	if (name == "nsamples")			return FTR_NSAMPLES;
	if (name == "exists")			return FTR_EXISTS;
	if (name == "max_diff")			return FTR_MAX_DIFF;
	if (name == "first_time")		return FTR_FIRST_DAYS;
	if (name == "category_set_first")				return FTR_CATEGORY_SET_FIRST;


	if (isInteger(name))
		return (BasicFeatureTypes)med_stoi(name);
	else
		MTHROW_AND_ERR("unknown name [%s]\n", name.c_str());
}

//.......................................................................................
void BasicFeatGenerator::set_names() {

	names.clear();
	string name = signalName + ".";
	//string name = signalName + ".";
	string set_names = in_set_name;
	if (set_names == "" && this->sets.size() > 0)
		set_names = boost::algorithm::join(this->sets, "_");
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
	case FTR_CATEGORY_SET:			name += "category_set_" + set_names; break;
	case FTR_CATEGORY_SET_COUNT:	name += "category_set_count_" + set_names; break;
	case FTR_CATEGORY_SET_SUM:		name += "category_set_sum_" + set_names; break;
	case FTR_CATEGORY_SET_FIRST:	name += "category_set_first_" + set_names; break;
	case FTR_NSAMPLES:			name += "nsamples"; break;
	case FTR_EXISTS:			name += "exists"; break;
	case FTR_MAX_DIFF:			name += "max_diff"; break;
	case FTR_FIRST_DAYS:		name += "first_time"; break;

	default: {
		name += "ERROR";
		MTHROW_AND_ERR("Got a wrong type in basic feature generator %d\n", type);
	}
	}

	name += ".win_" + std::to_string(win_from) + "_" + std::to_string(win_to);
	if (type == FTR_WIN_DELTA_VALUE)
		name += "_" + std::to_string(d_win_from) + "_" + std::to_string(d_win_to);
	if (time_channel != 0 || val_channel != 0)
		name += ".t" + std::to_string(time_channel) + "v" + std::to_string(val_channel);
	names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + name);
	// add the undecorated feature name as a tag, so we can later remove/select it with TagFeatureSelector
	tags.push_back(name);
	//MLOG("Created %s\n", name.c_str());

	//time_unit_sig = rep.sigs.Sid2Info[sid].time_unit; !! this is an issue to SOLVE !!
}

// Init
//.......................................................................................
void BasicFeatGenerator::init_defaults() {
	generator_type = FTR_GEN_BASIC;
	signalId = -1;
	time_unit_sig = MedTime::Undefined;
	time_unit_win = global_default_time_unit;
	string _signalName = "";
	bound_outcomeTime = false;
	//set(_signalName, FTR_LAST, 0, 360000);
};

// Generate
//.......................................................................................
int BasicFeatGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {


	if (time_unit_sig == MedTime::Undefined)	time_unit_sig = rec.my_base_rep->sigs.Sid2Info[signalId].time_unit;

	float *p_feat = p_data[0] + index;
	MedSample *p_samples = &(features.samples[index]);

	for (int i = 0; i < num; i++)
		p_feat[i] = get_value(rec, i, med_time_converter.convert_times(features.time_unit, time_unit_win, p_samples[i].time),
			med_time_converter.convert_times(features.time_unit, time_unit_sig, p_samples[i].outcomeTime));

	return 0;
}

// Init look-up table
//.......................................................................................
void BasicFeatGenerator::init_tables(MedDictionarySections& dict) {

	if (type == FTR_CATEGORY_SET || type == FTR_CATEGORY_SET_COUNT || type == FTR_CATEGORY_SET_SUM || type == FTR_CATEGORY_SET_FIRST) {
		if (lut.size() == 0) {
			int section_id = dict.section_id(signalName);
			//MLOG("BEFORE_LEARN:: signalName %s section_id %d sets size %d sets[0] %s\n", signalName.c_str(), section_id, sets.size(), sets[0].c_str());
			dict.prep_sets_lookup_table(section_id, sets, lut);
			//MLOG("AFTER_LEARN:: signalName %s section_id %d sets size %d sets[0] %s LUT %d\n", signalName.c_str(), section_id, sets.size(), sets[0].c_str(), lut.size());
		}
	}
	else
		lut.clear();

	return;
}

//.......................................................................................
float BasicFeatGenerator::get_value(PidDynamicRec& rec, int idx, int time, int outcomeTime) {

	rec.uget(signalId, idx);

	switch (type) {
	case FTR_LAST_VALUE:	return uget_last(rec.usv, time, win_from, win_to, outcomeTime);
	case FTR_FIRST_VALUE:	return uget_first(rec.usv, time, outcomeTime);
	case FTR_LAST2_VALUE:	return uget_last2(rec.usv, time, outcomeTime);
	case FTR_AVG_VALUE:		return uget_avg(rec.usv, time, outcomeTime);
	case FTR_MAX_VALUE:		return uget_max(rec.usv, time, outcomeTime);
	case FTR_MIN_VALUE:		return uget_min(rec.usv, time, outcomeTime);
	case FTR_STD_VALUE:		return uget_std(rec.usv, time, outcomeTime);
	case FTR_LAST_DELTA_VALUE:	return uget_last_delta(rec.usv, time, outcomeTime);
	case FTR_LAST_DAYS:			return uget_last_time(rec.usv, time, outcomeTime);
	case FTR_LAST2_DAYS:		return uget_last2_time(rec.usv, time, outcomeTime);
	case FTR_SLOPE_VALUE:		return uget_slope(rec.usv, time, outcomeTime);
	case FTR_WIN_DELTA_VALUE:	return uget_win_delta(rec.usv, time, outcomeTime);
	case FTR_CATEGORY_SET:				return uget_category_set(rec, rec.usv, time, outcomeTime);
	case FTR_CATEGORY_SET_COUNT:		return uget_category_set_count(rec, rec.usv, time, outcomeTime);
	case FTR_CATEGORY_SET_SUM:			return uget_category_set_sum(rec, rec.usv, time, outcomeTime);
	case FTR_NSAMPLES:			return uget_nsamples(rec.usv, time, win_from, win_to, outcomeTime);
	case FTR_EXISTS:			return uget_exists(rec.usv, time, win_from, win_to, outcomeTime);
	case FTR_MAX_DIFF:			return uget_max_diff(rec.usv, time, outcomeTime);
	case FTR_FIRST_DAYS:		return uget_first_time(rec.usv, time, outcomeTime);
	case FTR_CATEGORY_SET_FIRST:		return uget_category_set_first(rec, rec.usv, time, outcomeTime);

	default:	return missing_val;
	}

	return missing_val;
}

// Init
//.......................................................................................
int BasicFeatGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [BasicFeatGenerator::init]
		if (field == "type") { type = name_to_type(entry.second); }
		else if (field == "win_from") win_from = med_stoi(entry.second);
		else if (field == "win_to") win_to = med_stoi(entry.second);
		else if (field == "d_win_from") d_win_from = med_stoi(entry.second);
		else if (field == "d_win_to") d_win_to = med_stoi(entry.second);
		else if (field == "signalName" || field == "signal") signalName = entry.second;
		else if (field == "time_unit") time_unit_win = med_time_converter.string_to_type(entry.second);
		else if (field == "time_channel") time_channel = med_stoi(entry.second);
		else if (field == "val_channel") val_channel = med_stoi(entry.second);
		else if (field == "sum_channel") sum_channel = med_stoi(entry.second);
		else if (field == "sets") boost::split(sets, entry.second, boost::is_any_of(","));
		else if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "in_set_name") in_set_name = entry.second;
		else if (field == "bound_outcomeTime") bound_outcomeTime = stoi(entry.second) > 0;
		else if (field == "weights_generator") iGenerateWeights = med_stoi(entry.second);
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for BasicFeatGenerator\n", field.c_str());
		//! [BasicFeatGenerator::init]
	}

	// names for BasicFeatGenerator are set as a first step in the Learn call as we must have access to the MedRepository
	names.clear();

	set_names();

	req_signals.assign(1, signalName);

	return 0;
}

//=======================================================================================
// Age
//=======================================================================================
int AgeGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	// Sanity check
	if (signalId == -1)
		MTHROW_AND_ERR("Uninitialized signalId in age generation\n");

	float *p_feat = p_data[0] + index;

	int len;

	SVal *sig = (SVal *)rec.get(signalId, len);
	if (len != 1) { MTHROW_AND_ERR("id %d , got len %d for signal %d [%s])...\n", rec.pid, len, signalId, signalName.c_str()); }
	if (len == 0) throw MED_EXCEPTION_NO_BYEAR_GIVEN;
	if (signalName == "BYEAR") {
		int byear = (int)(sig[0].val);
		for (int i = 0; i < num; i++)
			p_feat[i] = (float)(med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time) / 10000 - byear);
	}
	else if (signalName == "BDATE") {
		int bdate = (int)(sig[0].val);
		for (int i = 0; i < num; i++) {
			int time = med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time);
			int days_since_birth = get_day_approximate(time) - get_day_approximate(bdate);
			p_feat[i] = (float)(1.0 * days_since_birth) / 365;

		}
	}
	else MTHROW_AND_ERR("Unknown age signal [%s] \n", signalName.c_str());

	return 0;
}

int AgeGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "signal")
			signalName = entry.second;
		else if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter [%s] for AgeGenerator\n", field.c_str());
	}
	set_names();

	req_signals.clear();
	req_signals.push_back(signalName);
	return 0;
}

//=======================================================================================
// Gender
//=======================================================================================
int GenderGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	// Sanity check
	if (genderId == -1) {
		MERR("Uninitialized genderId\n");
		return -1;
	}

	rec.uget(genderId, 0);
	if (rec.usv.len == 0) throw MED_EXCEPTION_NO_GENDER_GIVEN;
	int gender = (int)(rec.usv.Val(0));

	float *p_feat = p_data[0] + index;
	for (int i = 0; i < num; i++)
		p_feat[i] = (float)gender;

	return 0;
}


// Init
//.......................................................................................
int GenderGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [SingletonGenerator::init]
		if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "weights_generator") iGenerateWeights = med_stoi(entry.second);
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for SingletonGenerator\n", field.c_str());
		//! [SingletonGenerator::init]
	}

	// naming and required signals

	names.clear();
	set_names();

	req_signals.assign(1, "GENDER");

	return 0;
}


//=======================================================================================
// Singleton
//=======================================================================================
int SingletonGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	rec.uget(signalId, 0);
	float value;
	if (rec.usv.len == 0)
		value = missing_val;
	else {
		if (sets.size() == 0)
		{
			// Normal Singleton, just return value
			value = (float)((int)(rec.usv.Val(0)));
		}
		else
		{
			// Categorial Variable - check whether exists in LUT. Return 0/1
			value = (float)lut[((int)(rec.usv.Val(0)))];
		}
	}
	float *p_feat = p_data[0] + index;
	for (int i = 0; i < num; i++)
		p_feat[i] = value;

	return 0;
}

void SingletonGenerator::set_names()
{
	if (names.empty()) {
		string name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + signalName ; 
		//string name = signalName + ".";
		string set_names = in_set_name;
		if (set_names == "" && this->sets.size() > 0)
			set_names = boost::algorithm::join(this->sets, "_");

		if (set_names != "")
			name += ".category_set_" + set_names;

		names.push_back(name);
		//MLOG("Created %s\n", name.c_str());
	}
}

void SingletonGenerator::init_tables(MedDictionarySections& dict) {
	MLOG("sets size = %d \n", lut.size());
	if (sets.size() > 0) {
		// This is a categorial variable.
		if (lut.size() == 0) {
			int section_id = dict.section_id(signalName);
			//MLOG("BEFORE_LEARN:: signalName %s section_id %d sets size %d sets[0] %s\n", signalName.c_str(), section_id, sets.size(), sets[0].c_str());
			dict.prep_sets_lookup_table(section_id, sets, lut);
			//MLOG("AFTER_LEARN:: signalName %s section_id %d sets size %d sets[0] %s LUT %d\n", signalName.c_str(), section_id, sets.size(), sets[0].c_str(), lut.size());
		}
	}
	else
		lut.clear();

	return;
}

// Init
//.......................................................................................
int SingletonGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [SingletonGenerator::init]
		if (field == "signalName" || field == "signal") signalName = entry.second;
		else if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "weights_generator") iGenerateWeights = med_stoi(entry.second);
		else if (field == "sets") boost::split(sets, entry.second, boost::is_any_of(","));
		else if (field == "in_set_name") in_set_name = entry.second;
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for SingletonGenerator\n", field.c_str());
		//! [SingletonGenerator::init]
	}

	// naming and required signals

	names.clear();
	set_names();

	req_signals.assign(1, signalName);

	return 0;
}

//=======================================================================================
// ComorbidityGenerator
//=======================================================================================

//................................................................................................................
void RangeFeatGenerator::set_names() {

	names.clear();

	string name = signalName + ".";

	switch (type) {
	case FTR_RANGE_CURRENT:	name += "current"; break;
	case FTR_RANGE_LATEST:	name += "latest"; break;
	case FTR_RANGE_MIN:		name += "min"; break;
	case FTR_RANGE_MAX:		name += "max"; break;
	case FTR_RANGE_EVER:	name += "ever_" + sets[0]; break;
	case FTR_RANGE_TIME_DIFF: name += "time_diff_" + to_string(check_first) + sets[0]; break;
	default: {
		name += "ERROR";
		MTHROW_AND_ERR("Got a wrong type in range feature generator %d\n", type);
	}
	}

	name += ".win_" + std::to_string(win_from) + "_" + std::to_string(win_to);
	if (val_channel != 0)
		name += ".v" + std::to_string(val_channel);
	names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + name);
	// add the undecorated feature name as a tag, so we can later remove/select it with TagFeatureSelector
	tags.push_back(name);
}

// Init
//.......................................................................................
int RangeFeatGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RangeFeatGenerator::init]
		if (field == "type") { type = name_to_type(entry.second); }
		else if (field == "win_from") win_from = med_stoi(entry.second);
		else if (field == "win_to") win_to = med_stoi(entry.second);
		else if (field == "signalName" || field == "signal") signalName = entry.second;
		else if (field == "time_unit") time_unit_win = med_time_converter.string_to_type(entry.second);
		else if (field == "val_channel") val_channel = med_stoi(entry.second);
		else if (field == "sets") boost::split(sets, entry.second, boost::is_any_of(","));
		else if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "weights_generator") iGenerateWeights = med_stoi(entry.second);
		else if (field == "check_first") check_first = med_stoi(entry.second);
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for RangeFeatGenerator\n", field.c_str());
		//! [RangeFeatGenerator::init]
	}

	// set names and required signals
	set_names();

	req_signals.assign(1, signalName);

	return 0;
}

void RangeFeatGenerator::init_tables(MedDictionarySections& dict) {

	if (type == FTR_RANGE_EVER || type == FTR_RANGE_TIME_DIFF) {
		if (lut.size() == 0) {
			int section_id = dict.section_id(signalName);
			dict.prep_sets_lookup_table(section_id, sets, lut);
		}
	}
	else
		lut.clear();
	return;
}

// Init
//.......................................................................................
void RangeFeatGenerator::init_defaults() {
	generator_type = FTR_GEN_RANGE;
	signalId = -1;
	sets.clear();
	time_unit_sig = MedTime::Undefined;
	time_unit_win = global_default_time_unit;
	string _signalName = "";
	set(_signalName, FTR_RANGE_CURRENT, 0, 360000);
};

// Get type from name
//.......................................................................................
RangeFeatureTypes RangeFeatGenerator::name_to_type(const string &name)
{

	if (name == "current")				return FTR_RANGE_CURRENT;
	if (name == "latest")			return FTR_RANGE_LATEST;
	if (name == "max")			return FTR_RANGE_MAX;
	if (name == "min")			return FTR_RANGE_MIN;
	if (name == "ever")			return FTR_RANGE_EVER;
	if (name == "time_diff")  return FTR_RANGE_TIME_DIFF;

	return (RangeFeatureTypes)med_stoi(name);
}

// Generate
//.......................................................................................
int RangeFeatGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	if (time_unit_sig == MedTime::Undefined)	time_unit_sig = rec.my_base_rep->sigs.Sid2Info[signalId].time_unit;

	float *p_feat = p_data[0] + index;
	MedSample *p_samples = &(features.samples[index]);

	for (int i = 0; i < num; i++)
		p_feat[i] = get_value(rec, i, med_time_converter.convert_times(features.time_unit, time_unit_win, p_samples[i].time));

	return 0;
}

//.......................................................................................
float RangeFeatGenerator::get_value(PidDynamicRec& rec, int idx, int time) {

	rec.uget(signalId, idx);

	switch (type) {
	case FTR_RANGE_CURRENT:	return uget_range_current(rec.usv, time);
	case FTR_RANGE_LATEST:	return uget_range_latest(rec.usv, time);
	case FTR_RANGE_MIN:	return uget_range_min(rec.usv, time);
	case FTR_RANGE_MAX:		return uget_range_max(rec.usv, time);
	case FTR_RANGE_EVER:		return uget_range_ever(rec.usv, time);
	case FTR_RANGE_TIME_DIFF: 	return uget_range_time_diff(rec.usv, time);

	default:	return missing_val;
	}

	return missing_val;
}


//................................................................................................................
// in all following uget funcs the relevant time window is [min_time, max_time] and time is given in time_unit_win
//................................................................................................................

// get the last value in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_last(UniversalSigVec &usv, int time, int _win_from, int _win_to, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(_win_from, _win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	for (int i = usv.len - 1; i >= 0; i--) {
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
float BasicFeatGenerator::uget_first(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

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
float BasicFeatGenerator::uget_last2(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	for (int i = usv.len - 1; i >= 0; i--) {
		if (usv.Time(i, time_channel) <= max_time) {
			if (i > 0 && usv.Time(i - 1, time_channel) >= min_time)
				return usv.Val(i - 1, val_channel);
			else
				return missing_val;
		}
	}

	return missing_val;
}

//.......................................................................................
// get the average value in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_avg(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

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
		return (float)(sum / nvals);

	return missing_val;
}

//.......................................................................................
// get the max value in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_max(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

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
float BasicFeatGenerator::uget_min(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

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
float BasicFeatGenerator::uget_std(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

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
		double avg = sum / nvals;
		double var = sum_sq / nvals - avg*avg;
		if (var < 0.0001) var = 0.0001;
		return (float)sqrt(var);
	}

	return missing_val;
}

//.......................................................................................
float BasicFeatGenerator::uget_last_delta(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	for (int i = usv.len - 1; i >= 0; i--) {
		if (usv.Time(i, time_channel) <= max_time) {
			if (i > 0 && usv.Time(i - 1, time_channel) >= min_time)
				return (usv.Val(i, val_channel) - usv.Val(i - 1, val_channel));
			else
				return missing_val;
		}
	}
	return missing_val;
}

//.......................................................................................
float BasicFeatGenerator::uget_last_time(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	for (int i = usv.len - 1; i >= 0; i--) {
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
float BasicFeatGenerator::uget_first_time(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime >= min_time) {
			if (itime > max_time)
				return missing_val;
			else
				return (float)(time - usv.TimeU(i, time_channel, time_unit_win));
		}
	}
	return missing_val;
}

//.......................................................................................
float BasicFeatGenerator::uget_last2_time(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	for (int i = usv.len - 1; i >= 0; i--) {
		if (usv.Time(i, time_channel) <= max_time) {
			if (i > 0 && usv.Time(i - 1, time_channel) >= min_time)
				return (float)(time - usv.TimeU(i - 1, time_channel, time_unit_win));
			else
				return missing_val;
		}
	}

	return missing_val;
}

//.......................................................................................
// get the slope in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_slope(UniversalSigVec &usv, int time, int outcomeTime)
{

	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	double sx = 0, sy = 0, sxx = 0, sxy = 0, n = 0;
	double t_start = -1;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time) {
			if (t_start < 0) t_start = usv.TimeU(i, time_channel, time_unit_win);
			double t_curr = usv.TimeU(i, time_channel, time_unit_win);
			double x = (t_curr - t_start) / 500.0;
			double ival = usv.Val(i, val_channel);
			sx += x;
			sy += ival;
			sxx += x*x;
			sxy += x*ival;
			n++;
		}
	}

	if (n < 2) return missing_val;

	double cov = sxy - sx*(sy / n);
	double var = sxx - sx*(sx / n);

	if (var < 0.1)		return 0;

	return ((float)(cov / var));

}

//.......................................................................................
float BasicFeatGenerator::uget_win_delta(UniversalSigVec &usv, int time, int outcomeTime)
{
	float val1 = uget_last(usv, time, win_from, win_to, outcomeTime);
	if (val1 == missing_val) return missing_val;

	float val2 = uget_last(usv, time, d_win_from, d_win_to, outcomeTime);
	if (val2 == missing_val) return missing_val;


	return (val1 - val2);
}

//.......................................................................................
float BasicFeatGenerator::uget_category_set(PidDynamicRec &rec, UniversalSigVec &usv, int time, int outcomeTime)
{

	//#pragma omp critical
	//if (lut.size() == 0) {
	//		int section_id = rec.my_base_rep->dict.section_id(signalName);
	//		rec.my_base_rep->dict.prep_sets_lookup_table(section_id, sets, temp_lut);
	//	}

	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time && lut[(int)usv.Val(i, val_channel)]) 	return 1;
	}

	return 0;
}

float BasicFeatGenerator::uget_category_set_first(PidDynamicRec &rec, UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) return 0; // passed window
		if (lut[(int)usv.Val(i, val_channel)]) // what we look for
			return itime >= min_time; // inside window
	}

	return 0;
}

//.......................................................................................
float BasicFeatGenerator::uget_category_set_count(PidDynamicRec &rec, UniversalSigVec &usv, int time, int outcomeTime)
{
	//#pragma omp critical
	//	if (lut.size() == 0) {
	//		int section_id = rec.my_base_rep->dict.section_id(signalName);
	//		rec.my_base_rep->dict.prep_sets_lookup_table(section_id, sets, lut);
	//	}

	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	int cnt = 0;
	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time && lut[(int)usv.Val(i, val_channel)]) 	cnt++;
	}

	return (float)cnt;
}

//.......................................................................................
float BasicFeatGenerator::uget_category_set_sum(PidDynamicRec &rec, UniversalSigVec &usv, int time, int outcomeTime)
{
	//#pragma omp critical
	//	if (lut.size() == 0) {
	//		int section_id = rec.my_base_rep->dict.section_id(signalName);
	//		rec.my_base_rep->dict.prep_sets_lookup_table(section_id, sets, lut);
	//	}

	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	float sum = 0;
	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime > max_time) break;
		if (itime >= min_time && lut[(int)usv.Val(i, val_channel)]) 	sum += usv.Val(i, sum_channel);
	}

	return sum;
}

//.......................................................................................
// get the number of samples in [win_to, win_from] before time
float BasicFeatGenerator::uget_nsamples(UniversalSigVec &usv, int time, int _win_from, int _win_to, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;
	int i, j;
	for (i = usv.len - 1; i >= 0 && usv.Time(i, time_channel) > max_time; i--);
	for (j = 0; j < usv.len && usv.Time(j, time_channel) < min_time; j++);
	if (usv.Time(i, time_channel) <= max_time && usv.Time(j, time_channel) >= min_time) return (float)i - j + 1;
	return 0;
}

//.......................................................................................
// get 1.0 if there were any samples in [win_to, win_from] before time, else 0.0
float BasicFeatGenerator::uget_exists(UniversalSigVec &usv, int time, int _win_from, int _win_to, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;
	int i, j;
	for (i = usv.len - 1; i >= 0 && usv.Time(i, time_channel) > max_time; i--);
	for (j = 0; j < usv.len && usv.Time(j, time_channel) < min_time; j++);
	if (usv.Time(i, time_channel) <= max_time && usv.Time(j, time_channel) >= min_time && i - j >= 0)
		return 1.0;
	else return 0.0;
}

//.......................................................................................
// get values for RangeFeatGenerator
//.......................................................................................
// get the value in a range that includes time - win_from, if available
float RangeFeatGenerator::uget_range_current(UniversalSigVec &usv, int time)
{
	int dummy_time, time_to_check;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, dummy_time, time_to_check);

	for (int i = 0; i < usv.len; i++) {
		int fromTime = usv.Time(i, 0);
		int toTime = usv.Time(i, 1);

		if (fromTime > time_to_check)
			break;
		else if (toTime >= time_to_check)
			return usv.Val(i, val_channel);
	}

	return missing_val;
}

//.......................................................................................
// get the value in the latest range that intersets with time-win_to to time-win_from
float RangeFeatGenerator::uget_range_latest(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	float val = missing_val;
	for (int i = 0; i < usv.len; i++) {
		int fromTime = usv.Time(i, 0);
		int toTime = usv.Time(i, 1);

		if (fromTime > max_time)
			break;
		else if (toTime < min_time)
			continue;
		else if (fromTime >= min_time || toTime <= max_time)
			val = usv.Val(i, val_channel);
	}

	return val;
}

//.......................................................................................
// get the minimal value in a range that intersets with time-win_to to time-win_from
float RangeFeatGenerator::uget_range_min(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	float min_val = (float)1e20;

	for (int i = 0; i < usv.len; i++) {
		int fromTime = usv.Time(i, 0);
		int toTime = usv.Time(i, 1);

		if (fromTime > max_time)
			break;
		else if (toTime < min_time)
			continue;
		else if ((fromTime >= min_time || toTime <= max_time) && usv.Val(i, val_channel) < min_val)
			min_val = usv.Val(i, val_channel);
	}

	if (min_val < (float)1e20)
		return min_val;
	else
		return missing_val;
}

//.......................................................................................
// get the maximal value in a range that intersets with time-win_to to time-win_from
float RangeFeatGenerator::uget_range_max(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	float max_val = (float)-1e10;

	for (int i = 0; i < usv.len; i++) {
		int fromTime = usv.Time(i, 0);
		int toTime = usv.Time(i, 1);

		if (fromTime > max_time)
			break;
		else if (toTime < min_time)
			continue;
		else if ((fromTime >= min_time || toTime <= max_time) && usv.Val(i, val_channel) > max_val)
			max_val = usv.Val(i, val_channel);
	}

	if (max_val > (float)-1e10)
		return max_val;
	else
		return missing_val;
}

//.......................................................................................
// returns 1 if the range ever (up to time) had the value signalValue
float RangeFeatGenerator::uget_range_ever(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	for (int i = 0; i < usv.len; i++) {
		int fromTime = usv.Time(i, 0);
		int toTime = usv.Time(i, 1);

		if (fromTime > max_time)
			break;
		else if (toTime < min_time)
			continue;
		else if (lut[(int)usv.Val(i, val_channel)]) 	return 1;
	}
	return 0.0;
}

//.......................................................................................
// returns time diff if the range ever (up to time) had the value signalValue
float RangeFeatGenerator::uget_range_time_diff(UniversalSigVec &usv, int time)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);

	int no_lut_ind = 0;
	float time_diff = missing_val;
	for (int i = 0; i < usv.len; i++) {
		int fromTime = usv.Time(i, 0);
		int toTime = usv.Time(i, 1);
		if (fromTime > max_time)
			break;
		else if (toTime < min_time)
			continue;
		else if (lut[(int)usv.Val(i, val_channel)]) {
			if (check_first == 1) {
				//in case of first range
				int max_time = fromTime;
				if (win_from > max_time) max_time = win_from;

				time_diff = (float)time - med_time_converter.convert_times(MedTime::Date, MedTime::Days, max_time);
				//fprintf(stderr, "max_time: %i time :%i from_time:%i win_from:%i time_diff:%i\n", max_time, time, fromTime, win_from, time_diff);
				return time_diff;
			}
			else {
				//in case of last range
				int time_to_diff = toTime;
				if (win_to < toTime) time_to_diff = win_to;
				time_diff = (float)+time - med_time_converter.convert_times(MedTime::Date, MedTime::Days, time_to_diff);
			}
		}
		else
			no_lut_ind = 1;
	}

	//in case of last range
	if (check_first == 0 && time_diff != missing_val)
		return time_diff;

	//in case of range exists but no lut
	if (no_lut_ind == 1) {
		time_diff = -1.0F* win_to;
		return time_diff;
	}
	//in case of no range in the time window
	else {
		return missing_val;
	}
}

// ModelFeatureGenerator
//=======================================================================================

//................................................................................................................
void ModelFeatGenerator::set_names() {
	names.clear();

	string name;
	if (modelName != "")
		name = modelName;
	else if (modelFile != "")
		name = modelFile;
	else
		name = "ModelPred";

	if (impute_existing_feature)
		names.push_back(modelName);
	else {
		for (int i = 0; i < n_preds; i++)
			names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + name + "." + to_string(i + 1));
	}
}

// Init from map
//.......................................................................................
int ModelFeatGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [ModelFeatGenerator::init]
		if (field == "name") modelName = entry.second;
		else if (field == "file") modelFile = entry.second;
		else if (field == "impute_existing_feature") impute_existing_feature = med_stoi(entry.second);
		else if (field == "n_preds") n_preds = med_stoi(entry.second);
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter \'%s\' for ModelFeatureGenerator\n", field.c_str());
		//! [ModelFeatGenerator::init]
	}

	// set names
	set_names();
	// {name} is magical
	boost::replace_all(modelFile, "{name}", modelName);
	// Read Model and get required signal
	MedModel *_model = new MedModel;
	if (_model->read_from_file(modelFile) != 0)
		MTHROW_AND_ERR("Cannot read model from binary file %s\n", modelFile.c_str());

	init_from_model(_model);
	return 0;
}

// Copy Model and get required signal
//.......................................................................................
int ModelFeatGenerator::init_from_model(MedModel *_model) {

	generator_type = FTR_GEN_MODEL;
	model = _model;

	unordered_set<string> required;
	model->get_required_signal_names(required);
	for (string signal : required)
		req_signals.push_back(signal);

	set_names();
	return 0;
}

/// Load predictions from a MedSamples object. Compare to the models MedSamples (unless empty)
//.......................................................................................
void ModelFeatGenerator::override_predictions(MedSamples& inSamples, MedSamples& modelSamples) {

	// Sanity check ...
	if (modelSamples.idSamples.size() && !inSamples.same_as(modelSamples, 0))
		MTHROW_AND_ERR("inSamples is not identical to model samples in ModelFeatGenerator::load\n");

	preds.resize(inSamples.nSamples()*n_preds);

	int idx = 0;
	for (auto& idSamples : inSamples.idSamples) {
		for (auto& sample : idSamples.samples) {
			if (sample.prediction.size() < n_preds)
				MTHROW_AND_ERR("Cannot extract %d predictions from sample in ModelFeatGenerator::load\n", n_preds);

			for (int i = 0; i < n_preds; i++)
				preds[idx++] = sample.prediction[i];
		}
	}

	set_names();

	use_overriden_predictions = 1;

}

// Do the actual prediction prior to feature generation , only if vector is empty
//.......................................................................................
void ModelFeatGenerator::prepare(MedFeatures & features, MedPidRepository& rep, MedSamples& samples) {
	if (!use_overriden_predictions) {
		// Predict
		if (model->apply(rep, samples, MED_MDL_APPLY_FTR_GENERATORS, MED_MDL_APPLY_PREDICTOR) != 0)
			MTHROW_AND_ERR("ModelFeatGenerator::prepare feature %s failed to apply model\n", modelName.c_str());
		// Extract predictions
		if (model->features.samples[0].prediction.size() < n_preds)
			MTHROW_AND_ERR("ModelFeatGenerator::prepare cannot generate feature %s\n", modelName.c_str());

		preds.resize(n_preds*model->features.samples.size());

		for (int i = 0; i < model->features.samples.size(); i++) {
			for (int j = 0; j < n_preds; j++) {
				float new_val = model->features.samples[i].prediction[j];
				if (!isfinite(new_val))
					MTHROW_AND_ERR("ModelFeatGenerator::prepare feature %s nan in row %d\n", modelName.c_str(), i);
				preds[i*n_preds + j] = new_val;
			}
		}

	}
	if (impute_existing_feature) {
		FeatureProcessor p;
		string res = p.resolve_feature_name(features, modelName);
		names.clear();
		names.push_back(res);
	}
	else {
		FeatureGenerator::prepare(features, rep, samples);
	}
}

// Put relevant predictions in place
//.......................................................................................
int ModelFeatGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {
	float *p_feat = p_data[0] + index;
	for (int i = 0; i < num; i++) {
		if (!impute_existing_feature || p_feat[i] == missing_val) {
			float new_val = preds[index*n_preds + i];
			if (!isfinite(new_val))
				MTHROW_AND_ERR("ModelFeatGenerator::_generate nan in row %d\n", index*n_preds + i);
			p_feat[i] = new_val;
		}
	}

	return 0;

}

// (De)Serialize
//.......................................................................................
size_t ModelFeatGenerator::get_size() {
	size_t size = MedSerialize::get_size(generator_type, tags, modelFile, modelName, n_preds, names, req_signals, impute_existing_feature);
	return size + model->get_size();
}

size_t ModelFeatGenerator::serialize(unsigned char *blob) {
	size_t ptr1 = MedSerialize::serialize(blob, generator_type, tags, modelFile, modelName, n_preds, names, req_signals, impute_existing_feature);
	size_t ptr2 = model->serialize(blob + ptr1);
	return ptr2+ptr1;
}

size_t ModelFeatGenerator::deserialize(unsigned char *blob) {

	size_t ptr1 = MedSerialize::deserialize(blob, generator_type, tags, modelFile, modelName, n_preds, names, req_signals, impute_existing_feature);
	model = new MedModel;
	size_t ptr2 = model->deserialize(blob + ptr1);
	return ptr2+ptr1;
}

//................................................................................................................
// Helper function for time conversion
//................................................................................................................
void get_window_in_sig_time(int _win_from, int _win_to, int _time_unit_win, int _time_unit_sig, int _win_time, int &_min_time, int &_max_time)
{
	_min_time = med_time_converter.convert_times(_time_unit_win, _time_unit_sig, _win_time - _win_to);
	_max_time = med_time_converter.convert_times(_time_unit_win, _time_unit_sig, _win_time - _win_from);
}


//.......................................................................................
// get the max diiference in values in the window [win_to, win_from] before time
float BasicFeatGenerator::uget_max_diff(UniversalSigVec &usv, int time, int outcomeTime)
{
	int min_time, max_time;
	get_window_in_sig_time(win_from, win_to, time_unit_win, time_unit_sig, time, min_time, max_time);
	if (bound_outcomeTime && outcomeTime < max_time)
		max_time = outcomeTime;

	float max_diff = missing_val;
	vector<float> _vals_vec;
	for (int i = 0; i < usv.len; i++) {
		int itime = usv.Time(i, time_channel);
		if (itime >= min_time) {
			if (itime > max_time)
				break;
			else {
				if (_vals_vec.size() > 0) {
					nth_element(_vals_vec.begin(), _vals_vec.begin() + _vals_vec.size() / 2, _vals_vec.end());
					float median_prev_val = _vals_vec[_vals_vec.size() / 2];
					//float prev_val = median_prev_val;
					float prev_val = _vals_vec.back();
					float diff = usv.Val(i, val_channel) - prev_val;
					if (diff > max_diff || max_diff == missing_val)
						max_diff = diff;
				}
				_vals_vec.push_back(usv.Val(i, val_channel));
			}
		}
	}
	return max_diff;
}