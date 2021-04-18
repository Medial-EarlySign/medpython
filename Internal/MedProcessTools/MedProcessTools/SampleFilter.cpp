#include "SampleFilter.h"
#include "Logger/Logger/Logger.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedProcessTools/MedProcessTools/FeatureProcess.h"
#include <MedUtils/MedUtils/MedUtils.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

#define LOCAL_SECTION LOG_SMPL_FILTER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// SampleFilter
//=======================================================================================
// Filter Types from names
//.......................................................................................
SampleFilterTypes sample_filter_name_to_type(const string& filter_name) {

	if (filter_name == "train")
		return SMPL_FILTER_TRN;
	else if (filter_name == "test")
		return SMPL_FILTER_TST;
	else if (filter_name == "outliers")
		return SMPL_FILTER_OUTLIERS;
	else if (filter_name == "match")
		return SMPL_FILTER_MATCH;
	else if (filter_name == "required")
		return SMPL_FILTER_REQ_SIGNAL;
	else if (filter_name == "basic")
		return SMPL_FILTER_BASIC;
	else
		return SMPL_FILTER_LAST;
}

// Initialization : given filter name
//.......................................................................................
SampleFilter* SampleFilter::make_filter(string filter_name) {

	return make_filter(sample_filter_name_to_type(filter_name));
}

// Initialization : given filter name and intialization string
//.......................................................................................
SampleFilter * SampleFilter::make_filter(string filter_name, string init_string) {

	return make_filter(sample_filter_name_to_type(filter_name), init_string);
}

//.......................................................................................
void *SampleFilter::new_polymorphic(string dname)
{
	CONDITIONAL_NEW_CLASS(dname, BasicTrainFilter);
	CONDITIONAL_NEW_CLASS(dname, BasicTestFilter);
	CONDITIONAL_NEW_CLASS(dname, OutlierSampleFilter);
	CONDITIONAL_NEW_CLASS(dname, MatchingSampleFilter);
	CONDITIONAL_NEW_CLASS(dname, RequiredSignalFilter);
	CONDITIONAL_NEW_CLASS(dname, BasicSampleFilter);
	MWARN("Warning in SampleFilter::new_polymorphic - Unsupported class %s\n", dname.c_str());
	return NULL;
}


// Initialization : given filter type
//.......................................................................................
SampleFilter * SampleFilter::make_filter(SampleFilterTypes filter_type) {

	if (filter_type == SMPL_FILTER_TRN)
		return new BasicTrainFilter;
	else if (filter_type == SMPL_FILTER_TST)
		return new BasicTestFilter;
	else if (filter_type == SMPL_FILTER_OUTLIERS)
		return new OutlierSampleFilter;
	else if (filter_type == SMPL_FILTER_MATCH)
		return new MatchingSampleFilter;
	else if (filter_type == SMPL_FILTER_REQ_SIGNAL)
		return new RequiredSignalFilter;
	else if (filter_type == SMPL_FILTER_BASIC)
		return new BasicSampleFilter;
	else
		return NULL;

}

// Initialization : given filter type and intialization string
//.......................................................................................
SampleFilter * SampleFilter::make_filter(SampleFilterTypes filter_type, string init_string) {

	SampleFilter *newSampleFilter = make_filter(filter_type);
	newSampleFilter->init_from_string(init_string);
	if (newSampleFilter->init_from_string(init_string) < 0)
		MTHROW_AND_ERR("Cannot init SampleFilter of type %d with init string \'%s\'\n", filter_type, init_string.c_str());
	return newSampleFilter;
}

//.......................................................................................

// (De)Serialize
// Add filter-type to (De)Serialization of inheriting class
//.......................................................................................
size_t SampleFilter::get_filter_size() {
	return sizeof(filter_type) + get_size();
}

//.......................................................................................
size_t SampleFilter::filter_serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &filter_type, sizeof(SampleFilterTypes)); ptr += sizeof(SampleFilterTypes);
	ptr += serialize(blob + ptr);

	return ptr;
}

// In-place filtering : apply filtering and copy
//.......................................................................................
int SampleFilter::filter(MedSamples& samples) {

	MedSamples out_samples;
	int rc = filter(samples, out_samples);

	if (rc == 0)
		samples = out_samples;

	return rc;
}

// In-place filtering with repository : apply filtering and copy
//.......................................................................................
int SampleFilter::filter(MedRepository& rep, MedSamples& samples) {

	MedSamples out_samples;

	int rc = filter(rep, samples, out_samples);

	if (rc == 0)
		samples = out_samples;

	return rc;
}

//=======================================================================================
// BasicTrainFilter : 
//			take all controls samples (outcome=0) and all cases before outcomeTime
//=======================================================================================
// Filter
//.......................................................................................
int BasicTrainFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {

	outSamples.time_unit = inSamples.time_unit;

	// Take only samples before outcome
	for (MedIdSamples& idSamples : inSamples.idSamples) {

		MedIdSamples outIdSamples(idSamples.id);

		for (MedSample& sample : idSamples.samples) {
			// Negative or pre-outcome
			if (sample.outcome == 0 || sample.outcomeTime > sample.time)
				outIdSamples.samples.push_back(sample);
		}

		if (outIdSamples.samples.size() > 0)
			outSamples.idSamples.push_back(outIdSamples);
	}

	return 0;
}

//=======================================================================================
// BasicTestFilter :
//			dummy filter - take everything
//=======================================================================================
// Filter
//.......................................................................................
int BasicTestFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {

	// Take them all
	outSamples = inSamples;

	return 0;
}

//=======================================================================================
// OutlierSampleFilter :
//	- A filter that remove samples with outlier-outcomes (suitable for regression)
//	  Outliers detection is done using MedValueCleaner's methods (through inheritance)
//=======================================================================================
// Filter
//.......................................................................................
int OutlierSampleFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {

	outSamples.time_unit = inSamples.time_unit;

	// Filter by value of outcome
	for (MedIdSamples& idSample : inSamples.idSamples) {
		MedIdSamples newIdSample(idSample.id);

		for (MedSample& sample : idSample.samples) {
			if (sample.outcome >= removeMin - NUMERICAL_CORRECTION_EPS && sample.outcome <= removeMax + NUMERICAL_CORRECTION_EPS)
				newIdSample.samples.push_back(sample);
		}

		if (newIdSample.samples.size() > 0)
			outSamples.idSamples.push_back(newIdSample);
	}

	return 0;

}

// Learning : check outlier-detection method and call appropriate learner (iterative/quantile)
//.......................................................................................
int OutlierSampleFilter::_learn(MedSamples& samples) {

	if (params.type == VAL_CLNR_ITERATIVE)
		return iterativeLearn(samples);
	else if (params.type == VAL_CLNR_QUANTILE)
		return quantileLearn(samples);
	else {
		MERR("Unknown cleaning type %d\n", params.type);
		return -1;
	}
}

// Learning : learn outliers using MedValueCleaner's iterative approximation of moments
//.......................................................................................
int OutlierSampleFilter::iterativeLearn(MedSamples& samples) {

	// Get all values
	vector<float> values;
	get_values(samples, values);

	return get_iterative_min_max(values);
}

// Learning : learn outliers using MedValueCleaner's quantile appeoximation of moments
//.......................................................................................
int OutlierSampleFilter::quantileLearn(MedSamples& samples) {

	// Get all values
	vector<float> values;
	get_values(samples, values);

	return get_quantile_min_max(values);
}

// Helper for learning - extract all outcomes from samples.
//.......................................................................................
void OutlierSampleFilter::get_values(MedSamples& samples, vector<float>& values) {

	for (MedIdSamples& idSample : samples.idSamples) {
		for (MedSample& sample : idSample.samples)
			values.push_back(sample.outcome);
	}
}

//=======================================================================================
// MatchingSampleFilter
//=======================================================================================

// Init from map
//.......................................................................................
int MatchingSampleFilter::init(map<string, string>& mapper) {

	vector<string> strata;

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "priceRatio") eventToControlPriceRatio = med_stof(entry.second);
		else if (field == "maxRatio") maxControlToEventRatio = med_stof(entry.second);
		else if (field == "verbose") verbose = med_stoi(entry.second);
		else if (field == "minGroup") min_group_size = med_stoi(entry.second);
		else if (field == "strata") {
			boost::split(strata, entry.second, boost::is_any_of(":"));
			for (string& stratum : strata) addMatchingStrata(stratum);
		}
		else if (field == "match_to_prior")
			match_to_prior = med_stof(entry.second);
		else
			MLOG("Unknonw parameter \'%s\' for MatchingSampleFilter\n", field.c_str());
	}

	return 0;
}

// Add a matching stratum defined by a string
// Possibilities are :
//	"age",resolution-in-years
//	"time",time-unit(string, from MedTime),resolution
//	"signal",name,resolution,time-window,time-unit(string, from MedTime)
//  "gender"
//.......................................................................................
int MatchingSampleFilter::addMatchingStrata(string& init_string) {

	vector<string> fields;

	boost::split(fields, init_string, boost::is_any_of(","));

	matchingParams newStrata;
	if (fields[0] == "age") {
		if (fields.size() > 2) {
			MERR("Wrong number of features for matching strata\n");
			return -1;
		}

		newStrata.match_type = SMPL_MATCH_AGE;
		newStrata.resolution = (fields.size() > 1) ? stof(fields[1]) : (float)1.0;

	}
	else if (fields[0] == "time") {
		if (fields.size() > 3 || fields.size() < 2) {
			MERR("Wrong number of features for matching strata\n");
			return -1;
		}

		newStrata.match_type = SMPL_MATCH_TIME;
		newStrata.matchingTimeUnit = med_time_converter.string_to_type(fields[1]);
		newStrata.resolution = (fields.size() > 2) ? stof(fields[2]) : (float)1.0;
	}
	else if (fields[0] == "signal") {
		if (fields.size() > 5 || fields.size() < 2) {
			MERR("Wrong number of features for matching strata\n");
			return -1;
		}

		newStrata.match_type = SMPL_MATCH_SIGNAL;
		newStrata.signalName = fields[1];
		newStrata.resolution = (fields.size() > 2) ? stof(fields[2]) : (float)1.0;
		newStrata.timeWindow = (fields.size() > 3) ? (int)stof(fields[3]) : (int)1.0;
		newStrata.windowTimeUnit = (fields.size() > 4) ? med_time_converter.string_to_type(fields[4]) : global_default_windows_time_unit;
	}
	else if (fields[0] == "feature") {
		if (fields.size() > 3 || fields.size() < 2) {
			MERR("Wrong number of features for matching strata\n");
			return -1;
		}

		newStrata.match_type = SMPL_MATCH_FTR;
		newStrata.featureName = fields[1];
		newStrata.resolution = (fields.size() > 2) ? stof(fields[2]) : (float)1.0;
	}
	else if (fields[0] == "gender") {
		if (fields.size() != 1) {
			MERR("Wrong number of features for matching strata\n");
			return -1;
		}

		newStrata.match_type = SMPL_MATCH_SIGNAL;
		newStrata.signalName = "GENDER";
		newStrata.resolution = 1.0;
		newStrata.timeWindow = 99999999;
		newStrata.windowTimeUnit = global_default_windows_time_unit;
	}
	else {
		MERR("Unknown matching strata type %s\n", fields[0].c_str());
		return -1;
	}

	matchingStrata.push_back(newStrata);
	return 0;
}

// Filter with repository
// Find signature of each sample according to matching strata and then use match_by_general()
//...........................................................................................
int MatchingSampleFilter::_filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {

	outSamples.time_unit = inSamples.time_unit;

	// Create a MedFatures object ..
	MedFeatures features;
	inSamples.export_to_sample_vec(features.samples);

	// Init helpers
	if (initHelpers(inSamples, features, rep) < 0)
		return -1;

	// Mark samples according to strata
	vector<string> signatures(features.samples.size());

	for (unsigned int i = 0; i < signatures.size(); i++)
		getSampleSignature(features.samples[i], features, i, rep, signatures[i]);

	// Do the filtering
	vector<int> filtered;
	if (match_to_prior <= 0 || match_to_prior >= 1)
		medial::process::match_by_general(features, signatures, filtered, eventToControlPriceRatio, maxControlToEventRatio, min_group_size, (verbose > 0));
	else
		medial::process::match_to_prior(features, signatures, match_to_prior, filtered);
	outSamples.import_from_sample_vec(features.samples);

	return 0;

}

// Filter without repository (return -1 if repository is required)
//.......................................................................................
int MatchingSampleFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {

	if (isRepRequired()) {
		MERR("Cannot perform required matching without repository\n");
		return -1;
	}
	else {

		MedRepository dummyRep;
		return filter(dummyRep, inSamples, outSamples);
	}
}

// Filter with matrix (return -1 if repository is required)
//.......................................................................................
int MatchingSampleFilter::_filter(MedFeatures& features, MedSamples& inSamples, MedSamples& outSamples) {


	if (isRepRequired()) {
		MERR("Cannot perform required matching without repository\n");
		return -1;
	}
	else {

		MedRepository dummyRep;
		return _filter(features, dummyRep, inSamples, outSamples);
	}
}

// Filter with matrix + repository
//.......................................................................................
int MatchingSampleFilter::_filter(MedFeatures& features, MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {

	// Init helpers
	if (initHelpers(inSamples, features, rep) < 0)
		return -1;

	// Mark samples according to strata
	vector<string> signatures(features.samples.size());

	for (unsigned int i = 0; i < signatures.size(); i++)
		getSampleSignature(features.samples[i], features, i, rep, signatures[i]);

	// Do the filtering
	vector<int> filtered;
	if (match_to_prior <= 0 || match_to_prior >= 1)
		medial::process::match_by_general(features, signatures, filtered, eventToControlPriceRatio, maxControlToEventRatio, min_group_size, (verbose > 0));
	else
		medial::process::match_to_prior(features, signatures, match_to_prior, filtered);
	outSamples.import_from_sample_vec(features.samples);


	return 0;
}

// Utilities
// Check if repository is needed for matching (strata includes signal/age)
//.......................................................................................
bool MatchingSampleFilter::isRepRequired() {

	for (auto& strata : matchingStrata) {
		if (strata.match_type == SMPL_MATCH_AGE || strata.match_type == SMPL_MATCH_SIGNAL)
			return true;
	}

	return false;
}

// Check if age is needed for matching
//.......................................................................................
bool MatchingSampleFilter::isAgeRequired() {

	for (auto& strata : matchingStrata) {
		if (strata.match_type == SMPL_MATCH_AGE)
			return true;
	}

	return false;
}

// initialize values of helpers
// return 0 upon success, -1 if any of the required signals does not appear in the dictionary or any of the required features is not given </returns>
//.......................................................................................
int MatchingSampleFilter::initHelpers(MedSamples& inSamples, MedFeatures& features, MedRepository& rep) {

	// Helpers
	// Time unit from samples
	samplesTimeUnit = inSamples.time_unit;

	if (rep.dict.read_state == 0) {
		MWARN("WARNING Rep dictionary is empty\n");
	}

	// Age : either as a signal or using BDATE
	if (isAgeRequired()) {
		bdateId = rep.dict.id("BDATE");
		if (bdateId == -1) {
			MERR("Cannot find signalId for BDATE\n");
			return -1;
		}
	}

	// Check time dependence of matching signals
	for (matchingParams& stratum : matchingStrata) {
		if (stratum.match_type == SMPL_MATCH_SIGNAL) {
			stratum.signalId = rep.dict.id(stratum.signalName);
			if (stratum.signalId == -1) {
				MERR("Cannot find signalId for %s\n", stratum.signalName.c_str());
				return -1;
			}

			stratum.isTimeDependent = (rep.sigs.Sid2Info[stratum.signalId].n_time_channels > 0);
			stratum.signalTimeUnit = rep.sigs.Sid2Info[stratum.signalId].time_unit;
		}
	}

	// Check features
	FeatureProcessor dummy;
	for (matchingParams& stratum : matchingStrata) {
		if (stratum.match_type == SMPL_MATCH_FTR)
			stratum.resolvedFeatureName = dummy.resolve_feature_name(features, stratum.featureName);
	}

	return 0;
}

// Indexing of a single sample according to strata
// an index is a colon-separated string of bins per stratum
//.......................................................................................
int MatchingSampleFilter::getSampleSignature(MedSample& sample, MedFeatures& features, int i, MedRepository& rep, string& signature) {

	signature = "";
	for (auto& stratum : matchingStrata) {
		if (addToSampleSignature(sample, stratum, features, i, rep, signature) < 0)
			return -1;
	}

	return 0;
}

// add indexing of a single sample according to a single stratum to sample's index
//.......................................................................................
int MatchingSampleFilter::addToSampleSignature(MedSample& sample, matchingParams& stratum, MedFeatures& features, int i, MedRepository& rep, string& signature) {

	int age;
	UniversalSigVec usv;
	int bin;

	if (stratum.match_type == SMPL_MATCH_TIME) {
		// Take binned time in 'matchingTimeUnit'
		int time = med_time_converter.convert_times(samplesTimeUnit, stratum.matchingTimeUnit, sample.time);
		int bin = (int)(time / stratum.resolution);
		signature += to_string(bin) + ":";
	}
	else if (stratum.match_type == SMPL_MATCH_AGE) {
		// Take binned age
		int bdate = medial::repository::get_value(rep, sample.id, bdateId);
		int byear = int(bdate / 10000);
		age = med_time_converter.convert_times(samplesTimeUnit, MedTime::Date, sample.time) / 10000 - byear;

		bin = (int)((float)age / stratum.resolution);
		signature += to_string(bin) + ":";
	}
	else if (stratum.match_type == SMPL_MATCH_FTR) {
		bin = (int)(features.data[stratum.resolvedFeatureName][i] / stratum.resolution);
		signature += to_string(bin) + ":";
	}
	else if (stratum.match_type == SMPL_MATCH_SIGNAL) {
		rep.uget(sample.id, stratum.signalId, usv);
		if (!stratum.isTimeDependent) {
			// Signal is not time dependent - take binned value
			bin = (int)(usv.Val(0) / stratum.resolution);
			signature += to_string(bin) + ":";
		}
		else {
			int target = med_time_converter.convert_times(samplesTimeUnit, stratum.windowTimeUnit, sample.time);
			int maxTime = med_time_converter.convert_times(samplesTimeUnit, stratum.signalTimeUnit, sample.time);
			int minTime = med_time_converter.convert_times(stratum.windowTimeUnit, stratum.signalTimeUnit, target - stratum.timeWindow);
			//			MLOG("units = %d/%d/%d time = %d Target = %d min = %d\n", samplesTimeUnit, stratum.signalTimeUnit, stratum.windowTimeUnit, sample.time, maxTime, minTime);

			string tempSignature = "NULL";
			// Find first value after maxTime and check previous value
			for (int idx = 0; idx < usv.len; idx++) {
				if (usv.Time(idx) > maxTime) {
					if (idx > 0 && usv.Time(idx - 1) >= minTime)
						tempSignature = to_string((int)(0.001 + usv.Val(idx - 1) / stratum.resolution));
					break;
				}
			}

			// Is last value between minTime and maxTime ? (missed by previous check)
			if (usv.len > 0 && usv.Time(usv.len - 1) <= maxTime && usv.Time(usv.len - 1) >= minTime)
				tempSignature = to_string((int)(0.001 + usv.Val(usv.len - 1) / stratum.resolution));

			signature += tempSignature + ":";
		}
	}
	else {
		MERR("Unknown matching type %d\n", stratum.match_type);
		return -1;
	}

	return 0;
}

//Get all signals required  for matching
//.......................................................................................
void MatchingSampleFilter::get_required_signals(vector<string>& req_sigs)
{
	req_sigs.clear();
	if (isAgeRequired())
		req_sigs.push_back("BDATE");

	for (auto &s : matchingStrata) {
		if (s.signalName != "")
			req_sigs.push_back(s.signalName);
	}

	return;

}

//=======================================================================================
// Required Signal Filter 
//	- Keep only samples with a required signal appearing in a time-window.
//	- OBSOLETE - REPLACED BY BasicSampleFilter. KEPT HERE FOR BACKWARD COMPETABILITY
//=======================================================================================

// Init
//.......................................................................................
int RequiredSignalFilter::init(map<string, string>& mapper) {

	vector<string> strata;

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "signalName") signalName = entry.second;
		else if (field == "timeWindow") timeWindow = med_stoi(entry.second);
		else if (field == "timeUnit") windowTimeUnit = med_time_converter.string_to_type(entry.second);
		else
			MLOG("Unknonw parameter \'%s\' for RequiredSampleFilter\n", field.c_str());
	}

	return 0;
}

// Init to defaults
//.......................................................................................
void RequiredSignalFilter::init_defaults() {

	timeWindow = 0;
	windowTimeUnit = global_default_windows_time_unit;
}

// Filter
//.......................................................................................
int RequiredSignalFilter::_filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {

	outSamples.time_unit = inSamples.time_unit;

	int signalId = rep.dict.id(signalName);
	int signalTimeUnit = rep.sigs.Sid2Info[signalId].time_unit;

	UniversalSigVec usv;
	for (auto& idSamples : inSamples.idSamples) {
		MedIdSamples outIdSamples(idSamples.id);

		rep.uget(idSamples.id, signalId, usv);
		int idx = 0;
		for (auto& sample : idSamples.samples) {

			int target = med_time_converter.convert_times(inSamples.time_unit, windowTimeUnit, sample.time);
			int maxTime = med_time_converter.convert_times(inSamples.time_unit, signalTimeUnit, sample.time);
			int minTime = med_time_converter.convert_times(windowTimeUnit, signalTimeUnit, target - timeWindow);
			//	MLOG("units = %d/%d/%d time = %d Target = %d min = %d\n", samplesTimeUnit, stratum.signalTimeUnit, stratum.windowTimeUnit, sample.time, maxTime, minTime);

			while (idx < usv.len) {
				if (usv.Time(idx) == maxTime) {
					outIdSamples.samples.push_back(sample);
					break;
				}
				else if (usv.Time(idx) > maxTime) {
					if (idx > 0 && usv.Time(idx - 1) >= minTime)
						outIdSamples.samples.push_back(sample);
					break;
				}
				idx++;
			}
		}

		if (!outIdSamples.samples.empty())
			outSamples.idSamples.push_back(outIdSamples);
	}

	return 0;
}

// Filter without repository : Return an error
//.......................................................................................
int RequiredSignalFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {
	MERR("A repository is required for Required-Signal Filter\n");
	return -1;
}

//=======================================================================================
// A general filter to allow the following basics:
// (1) min and max time of outcomeTime
// (2) option to allow for a signal to be in some given range (if there is a sample in some given window)
// (3) option to force N samples of a specific signal within the given range (in the given time window)
//
// examples:
// sig:TRAIN,min_val:1,max_val:1,min_Nvals:1
// sig:Creatinine,win_from:0,win_to:720,min_Nvals:2
//=======================================================================================

// BasicFilteringParams Initialization from string
//.......................................................................................
int BasicFilteringParams::init_from_string(const string &init_str)
{
	vector<string> fields;

	boost::split(fields, init_str, boost::is_any_of(":=,"));

	for (int i = 0; i < fields.size(); i++) {
		if (fields[i] == "sig") { sig_name = fields[++i]; }
		if (fields[i] == "min_val") { min_val = stof(fields[++i]); }
		if (fields[i] == "max_val") { max_val = stof(fields[++i]); }
		if (fields[i] == "win_from") { win_from = stoi(fields[++i]); }
		if (fields[i] == "win_to") { win_to = stoi(fields[++i]); }
		if (fields[i] == "min_Nvals") { min_Nvals = stoi(fields[++i]); }
		if (fields[i] == "time_ch") { time_channel = stoi(fields[++i]); }
		if (fields[i] == "val_ch") { val_channel = stoi(fields[++i]); }
	}


	return 0;
}

// Test filtering criteria
// Return 1 if passing and 0 otherwise
//.......................................................................................
int BasicFilteringParams::test_filter(MedSample &sample, MedRepository &rep, int win_time_unit)
{
	//MLOG("id %d sig_id %d %s time %d\n", sample.id, sig_id, sig_name.c_str(), sample.time);
	if (sig_id < 0) {
		if (sig_name == "Age") {
			use_byear = 1;
			sig_id = rep.sigs.sid("BDATE");
		}
		else
			sig_id = rep.sigs.sid(sig_name);
	}

	UniversalSigVec usv;

	rep.uget(sample.id, sig_id, usv);
	//MLOG("id %d sig_id %d len %d %f\n", sample.id, sig_id, usv.len, usv.Val(0));
	//MLOG("id %d sig_id %d len %d\n", sample.id, sig_id, usv.len);

	if (usv.len == 0 && min_Nvals > 0) return 0;
	//if (min_Nvals <= 0) return 1;

	// Special handling of age through byear
	if (use_byear) {
		int year = 1900 + med_time_converter.convert_times(global_default_time_unit, MedTime::Years, sample.time);
		int bdate_v = (int)usv.Val(0);
		int age = year - int(bdate_v / 10000);
		if (age < min_val || age > max_val)
			return 0;
		return 1;
	}

	// Otherwise ...
	if (usv.n_time_channels() == 0) {
		// timeless signal - checking the first
		//MLOG("id %d val %f\n", sample.id, usv.Val(0));
		if (usv.Val(0) < min_val || usv.Val(0) > max_val)
			return 0;
		return 1;
	}
	else {

		int ref_time = med_time_converter.convert_times(usv.time_unit(), win_time_unit, sample.time);

		// go over all values
		int nvals = 0;
		int nvals_in_range = 0;
		for (int i = 0; i < usv.len; i++) {

			// check if in relevant window
			int i_time = usv.Time(i, time_channel);
			int i_time_converted = med_time_converter.convert_times(usv.time_unit(), win_time_unit, i_time);
			int dtime = ref_time - i_time_converted;
			//MLOG("id %d i_time %d %f %d time %d %d dtime %d win %d %d\n", sample.id, i_time, usv.Val(i, val_channel), i_time_converted, sample.time, ref_time, dtime, win_from, win_to);
			if (dtime < win_from) break;
			nvals++;
			if (dtime <= win_to) {
				// in relevant time window, checking the value range
				float i_val = usv.Val(i, val_channel);
				//MLOG("i %d id %d i_val %f min %f max %f minNvals %d\n", i, sample.id, i_val, min_val, max_val, min_Nvals);
				if (i_val >= min_val && i_val <= max_val) {
					nvals_in_range++;
					if (min_Nvals > 0 && nvals_in_range >= min_Nvals)
						return 1;
				}
			}
		}

		if (min_Nvals < 0 && nvals == nvals_in_range)
			return 1;
	}

	return 0;
}

// Initialize from map
//.......................................................................................
int BasicSampleFilter::init(map<string, string>& mapper)
{
	req_sigs.clear();
	for (auto &m : mapper) {
		if (m.first == "min_sample_time") min_sample_time = stoi(m.second);
		if (m.first == "max_sample_time") max_sample_time = stoi(m.second);
		if (m.first == "win_time_unit") winsTimeUnit = med_time_converter.string_to_type(m.second);
		if (m.first == "min_bfilter") min_bfilter = med_time_converter.string_to_type(m.second);
		if (m.first == "bfilter") {
			vector<string> fields;
			boost::split(fields, m.second, boost::is_any_of("+"));
			for (auto &f : fields) {
				BasicFilteringParams bfp;
				bfp.init_from_string(f);
				bfilters.push_back(bfp);
				if (bfp.sig_name == "Age")
					req_sigs.push_back("BDATE");
				else
					req_sigs.push_back(bfp.sig_name);
			}
		}

	}

	if (bfilters.size() > 0 && (min_bfilter < 0 || min_bfilter > bfilters.size())) min_bfilter = (int)bfilters.size();

	return 0;
}

//.......................................................................................
void BasicSampleFilter::get_required_signals(vector<string> &reqs)
{
	if (req_sigs.size() == 0) {
		req_sigs.clear();
		for (auto &bf : bfilters) {
			if (bf.sig_name == "Age")
				req_sigs.push_back("BDATE");
			else
				req_sigs.push_back(bf.sig_name);
		}
	}

	reqs = req_sigs;
	return;
}

// Filter with repository
//.......................................................................................
int BasicSampleFilter::_filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples)
{
	// assumes rep is already loaded with relevant signals

	outSamples = inSamples;
	outSamples.idSamples.clear();

	for (auto &id_s : inSamples.idSamples) {

		MedIdSamples id_out;
		id_out = id_s;
		id_out.samples.clear();

		for (auto &in_s : id_s.samples) {
			int take_it = 1;
			if (in_s.time < min_sample_time || in_s.time > max_sample_time) take_it = 0;
			//MLOG("id %d time %d min %d max %d take_it %d\n", in_s.id, in_s.time, min_sample_time, max_sample_time,take_it);
			if (take_it) {
				take_it = 0;
				for (auto &bf : bfilters) {
					if (bf.test_filter(in_s, rep, winsTimeUnit)) {
						take_it++;
					}
				}

				if (take_it >= min_bfilter) id_out.samples.push_back(in_s);
			}
		}

		if (id_out.samples.size() > 0) {
			// id passed with some samples, keeping them
			outSamples.idSamples.push_back(id_out);
		}

	}

	outSamples.sort_by_id_date();
	return 0;
}

// Filter without repository
// relevant only if bfilters is empty. Otherwise, return -1
//.......................................................................................
int BasicSampleFilter::_filter(MedSamples& inSamples, MedSamples& outSamples)
{

	if (bfilters.empty()) {
		MedRepository dummy;
		return filter(dummy, inSamples, outSamples);
	}
	else {
		MERR("A repository is required for Required-Signal Filter\n");
		return -1;
	}
}

// SanitySimpleFilter Initialization from string
//.......................................................................................
int SanitySimpleFilter::init_from_string(const string &init_str)
{
	vector<string> fields;

	boost::split(fields, init_str, boost::is_any_of(":=;"));

	for (int i = 0; i < fields.size(); i++) {
		//MLOG("INIT: %s -> %s\n", fields[i].c_str(), fields[i+1].c_str());
		if (fields[i] == "sig") { sig_name = fields[++i]; }
		if (fields[i] == "min_val") { min_val = stof(fields[++i]); }
		if (fields[i] == "max_val") { max_val = stof(fields[++i]); }
		if (fields[i] == "win_from") { win_from = stoi(fields[++i]); }
		if (fields[i] == "win_to") { win_to = stoi(fields[++i]); }
		if (fields[i] == "min_Nvals") { min_Nvals = stoi(fields[++i]); }
		if (fields[i] == "max_Nvals") { max_Nvals = stoi(fields[++i]); }
		if (fields[i] == "min_left") { min_left = stoi(fields[++i]); }
		if (fields[i] == "max_outliers") { max_outliers = stoi(fields[++i]); }
		if (fields[i] == "time_ch") { time_channel = stoi(fields[++i]); }
		if (fields[i] == "val_ch") { val_channel = stoi(fields[++i]); }
		if (fields[i] == "win_time_unit") { win_time_unit = med_time_converter.string_to_type(fields[++i]); }
		if (fields[i] == "samples_time_unit") { samples_time_unit = med_time_converter.string_to_type(fields[++i]); }
		if (fields[i] == "values_in_dictionary") { values_in_dictionary = stoi(fields[++i]); }
		if (fields[i] == "allowed_values") {
			vector<string> svals;
			boost::split(svals, fields[++i], boost::is_any_of(","));
			for (auto &s : svals) allowed_values.insert(stof(s));
		}
	}


	return 0;
}

#define SANITY_FILTER_DBG 0

// Test filtering criteria
// Returns one of the codes defined as static in the h file
//.......................................................................................
int SanitySimpleFilter::test_filter(MedSample &sample, MedRepository &rep, int &nvals, int &noutliers)
{

#if SANITY_FILTER_DBG
	MLOG("SanitySimpleFilter::test_filter(1) ==> id %d sig_id %d %s time %d\n", sample.id, sig_id, sig_name.c_str(), sample.time);
#endif

	if (sig_id < 0) {
		if (boost::iequals(sig_name, "Age")) {
			sig_id = 0;
			bdate_id = rep.sigs.sid("BDATE");
			if (bdate_id < 0) {
				used_byear = true;
				bdate_id = rep.sigs.sid("BYEAR");
				if (bdate_id < 0)
					MWARN("WARNING: !!!! ===> Using SanitySimpleFilter for age but without BDATE/BYEAR... Are you using a repository with an AGE signal??\n");
			}
		}
		else
			sig_id = rep.sigs.sid(sig_name);
#if SANITY_FILTER_DBG
		MLOG("SanitySimpleFilter::test_filter(2) ==> id %d sig_id %d %s time %d\n", sample.id, sig_id, sig_name.c_str(), sample.time);
#endif
		}
	if (sig_id < 0)
		return SanitySimpleFilter::Signal_Not_Valid;

	// Age case
	if (sig_id == 0) {
		// TBD: Must make this work also for the cases in which Age is given as a signal
		if (bdate_id > 0) {
			// calculate using byear
			float y = 1900 + (float)med_time_converter.convert_times(samples_time_unit, MedTime::Years, sample.time);
			int bdate_val = medial::repository::get_value(rep, sample.id, bdate_id);
			int byear = bdate_val;
			if (byear < 0)
				return SanitySimpleFilter::Failed_Age_No_Byear;

			if (!used_byear)
				byear = int(bdate_val / 10000);
#if SANITY_FILTER_DBG
			MLOG("SanitySimpleFilter::test_filter(3) ====> AGE : id %d byear %f y %f time %d : age %f min_val %f max_val %f\n", sample.id, byear, y, sample.time, y - byear, min_val, max_val);
#endif
			float age = y - byear;

			if (age < min_val || age > max_val)
				return SanitySimpleFilter::Failed_Age;
		}
		else
			return SanitySimpleFilter::Failed_Age_No_Byear;
	}

	if (sig_id > 0) {
		// regular signal case
		if (section_id < 0 && sig_id > 0) {
			section_id = rep.dict.section_id(sig_name);
		}

		UniversalSigVec usv;

		rep.uget(sample.id, sig_id, usv);
#if SANITY_FILTER_DBG
		MLOG("SanitySimpleFilter::test_filter(3.5) id %d sig %s sig_id %d\n", sample.id, sig_name.c_str(), sig_id);
		MLOG("SanitySimpleFilter::test_filter(4) id %d sig_id %d len %d\n", sample.id, sig_id, usv.len);
		MLOG("SanitySimpleFilter::test_filter(5) id %d sig_id %d len %d min_Nvals %d max_Nvals %d\n", sample.id, sig_id, usv.len, min_Nvals, max_Nvals);
#endif
		nvals = 0;
		noutliers = 0;
		int n_not_in_dict = 0;
		int n_not_allowed = 0;
		int n_left = 0;
		if (usv.len == 0 && min_Nvals > 0) return SanitySimpleFilter::Failed_Min_Nvals;

		if (usv.n_time_channels() == 0) {
			// timeless signal

			nvals = usv.len;
			for (int i = 0; i < usv.len; i++) {
				float i_val = usv.Val(i);
				if (i_val < min_val || i_val > max_val)
					noutliers++;
				else
					n_left++;
				if (values_in_dictionary && section_id > 0) {
					if (rep.dict.dicts[section_id].Id2Name.find((int)i_val) == rep.dict.dicts[section_id].Id2Name.end())
						n_not_in_dict++;
				}
				if (allowed_values.size() > 0) {
					if (allowed_values.find(i_val) == allowed_values.end())
						n_not_allowed++;
				}

			}

		}
		else {

			int ref_time = med_time_converter.convert_times(usv.time_unit(), win_time_unit, sample.time);

			// go over all values
			for (int i = 0; i < usv.len; i++) {

				// check if in relevant window
				int i_time = usv.Time(i, time_channel);
				int i_time_converted = med_time_converter.convert_times(usv.time_unit(), win_time_unit, i_time);
				int dtime = ref_time - i_time_converted;
#if SANITY_FILTER_DBG
				MLOG("SanitySimpleFilter::test_filter(6) id %d i_time %d %f %d time %d %d dtime %d win %d %d\n", sample.id, i_time, usv.Val(i, val_channel), i_time_converted, sample.time, ref_time, dtime, win_from, win_to);
#endif
				if (dtime < win_from) break;
				if (dtime <= win_to) {
					nvals++;
					// in relevant time window, checking the value range
					float i_val = usv.Val(i, val_channel);
					if (i_val < min_val || i_val > max_val) noutliers++;
					else n_left++;
#if SANITY_FILTER_DBG
					MLOG("SanitySimpleFilter::test_filter(7) i %d id %d i_val %f min %f max %f minNvals %d nvals %d noutliers %d\n", i, sample.id, i_val, min_val, max_val, min_Nvals, nvals, noutliers);
#endif

					if (values_in_dictionary && section_id > 0) {
						//MLOG("dictionary test: section_id %d sig_name %s i_val %f dsize %d\n", section_id, sig_name.c_str(), i_val, rep.dict.dicts[section_id].Id2Name.size());
						if (rep.dict.dicts[section_id].Id2Name.find((int)i_val) == rep.dict.dicts[section_id].Id2Name.end())
							n_not_in_dict++;
					}

					if (allowed_values.size() > 0) {
						if (allowed_values.find(i_val) == allowed_values.end())
							n_not_allowed++;
					}
				}

			}
		}

#if SANITY_FILTER_DBG
		MLOG("SanitySimpleFilter::test_filter(8) ###>>> id %d time %d sig %s (len %d) : min %d max %d maxout %d : nvals %d noutliers %d not_in_dict %d not_allowed %d\n",
			sample.id, sample.time, sig_name.c_str(), usv.len, min_Nvals, max_Nvals, max_outliers, nvals, noutliers, n_not_in_dict, n_not_allowed);
#endif

		if (min_Nvals >= 0 && nvals < min_Nvals) return SanitySimpleFilter::Failed_Min_Nvals;
		if (max_Nvals >= 0 && nvals > max_Nvals) return SanitySimpleFilter::Failed_Max_Nvals;
		if (max_outliers >= 0 && noutliers > max_outliers) return SanitySimpleFilter::Failed_Outliers;
		if (values_in_dictionary && ((n_not_in_dict > 0) || section_id < 0)) return SanitySimpleFilter::Failed_Dictionary_Test;
		if ((allowed_values.size() > 0) && (n_not_allowed > 0)) return SanitySimpleFilter::Failed_Allowed_Values;
		if (min_left >= 0 && n_left < min_left) return SanitySimpleFilter::Failed_Not_Enough_Non_Outliers_Left;
			}

	return SanitySimpleFilter::Passed;

	return 0;
		}
