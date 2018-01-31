#include "SampleFilter.h"
#include "Logger/Logger/Logger.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
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
int BasicTrainFilter::_filter(MedSamples& inSamples,MedSamples& outSamples) {

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
	get_values(samples,values);

	return get_iterative_min_max(values);
}

// Learning : learn outliers using MedValueCleaner's quantile appeoximation of moments
//.......................................................................................
int OutlierSampleFilter::quantileLearn(MedSamples& samples) {

	// Get all values
	vector<float> values;
	get_values(samples,values);

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

// (De)Serialization
//.......................................................................................
size_t OutlierSampleFilter::get_size() {

	size_t size = 0;

	size += sizeof(int); // int take_log
	size += 2 * sizeof(float); // float removeMax, removeMin;

	return size;
}

//.......................................................................................
size_t OutlierSampleFilter::serialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(blob + ptr, &params.take_log, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &removeMax, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &removeMin, sizeof(float)); ptr += sizeof(float);

	return ptr;
}

//.......................................................................................
size_t OutlierSampleFilter::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(&params.take_log, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&removeMax, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&removeMin, blob + ptr, sizeof(float)); ptr += sizeof(float);

	return ptr;
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

		if (field == "priceRatio") eventToControlPriceRatio = stof(entry.second);
		else if (field == "maxRatio") matchMaxRatio = stof(entry.second);
		else if (field == "verbose") verbose = stoi(entry.second);
		else if (field == "strata") {
			boost::split(strata, entry.second, boost::is_any_of(":"));
			for (string& stratum : strata) addMatchingStrata(stratum);
		} else
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
		newStrata.windowTimeUnit = (fields.size() > 4) ? med_time_converter.string_to_type(fields[4]) : med_rep_type.windowTimeUnit;
	}
	else if (fields[0] == "gender") {
		if (fields.size() != 1) {
			MERR("Wrong number of features for matching strata\n");
			return -1;
		}

		newStrata.match_type = SMPL_MATCH_SIGNAL;
		newStrata.signalName = med_rep_type.genderSignalName;
		newStrata.resolution = 1.0;
		newStrata.timeWindow = 99999999;
		newStrata.windowTimeUnit =  med_rep_type.windowTimeUnit;
	}
	else {
		MERR("Unknown matching strata type %s\n", fields[0].c_str());
		return -1;
	}

	matchingStrata.push_back(newStrata);
	return 0;
}

// Filter with repository
// Find signature of each sample according to matching strata, thus assigning samples (cases+controls) to bins
// Find optimal case-control matching ratio according to bins and relative puhishment for
// removing cases and controls (eventToControlPriceRatio). Adjust ratio to maximal allowed 
// ratio (matchMaxRatio) and then sample randomly from each bin.
//.......................................................................................
int MatchingSampleFilter::_filter(MedRepository& rep, MedSamples& inSamples, MedSamples& outSamples) {

	outSamples.time_unit = inSamples.time_unit;

	// Init helpers
	if (initHelpers(inSamples, rep) < 0)
		return -1;

	// Mark samples according to strata
	map<string, pair<int, int>> cnts;
	map<string, vector<pair<int,int>>> control_ids;
	map<string, vector<pair<int, int>>> event_ids;

	string signature;
	pair<int, int> p0(0, 0);
	vector<pair<int, int>> empty;

	int idx = 0;
	for (unsigned int idIdx = 0; idIdx < inSamples.idSamples.size(); idIdx++) {
		for (unsigned int sampleIdx = 0; sampleIdx < inSamples.idSamples[idIdx].samples.size(); sampleIdx++) {

			MedSample& sample = inSamples.idSamples[idIdx].samples[sampleIdx];
			getSampleSignature(sample, rep, signature);

			if (cnts.find(signature) == cnts.end()) {
				cnts[signature] = p0;
				control_ids[signature] = empty;
				event_ids[signature] = empty;
			}

			if (sample.outcome > 0) {
				cnts[signature].first++;
				event_ids[signature].push_back({ idIdx,sampleIdx });
			}
			else {
				cnts[signature].second++;
				control_ids[signature].push_back({ idIdx,sampleIdx });
			}
		}
	}

	// Identify pairing ratio
	float opt_factor = get_pairing_ratio(cnts, eventToControlPriceRatio);
	float factor = opt_factor;
	if (factor > matchMaxRatio) {
		if (verbose) MLOG("updating factor {%8.3f} to maxFactor=%.3f\n", factor, matchMaxRatio);
		factor = matchMaxRatio;
	}
	else if (factor < 1 / matchMaxRatio) {
		if (verbose) MLOG("updating factor {%8.3f} to 1/(maxFactor=%.3f)\n", factor, matchMaxRatio);
		factor = 1/matchMaxRatio;
	}
	if (verbose) MLOG("opt ratio is %8.3f (effective factor = %8.3f)\n", opt_factor, factor);


	// sample controls and events
	int n_events = 0, n_ctrl = 0;
	vector<pair<int, int>> selected;
	for (auto it = cnts.begin(); it != cnts.end(); it++) {
		signature = it->first;
		int ev_cnt = it->second.first;
		int ctrl_cnt = it->second.second;
		if (ev_cnt == 0 || ctrl_cnt == 0)
			continue;
		int ntake_ctrl = (int)((float)ev_cnt * factor);
		//if (ntake == 0) ntake = 1 + (int)factor/2;
		if (ntake_ctrl > (int)control_ids[signature].size()) ntake_ctrl = (int)control_ids[signature].size();

		int ntake_ev = (int)((float)ctrl_cnt / factor);
		//if (ntake == 0) ntake = 1 + (int)factor/2;
		if (ntake_ev > (int)event_ids[signature].size()) ntake_ev = (int)event_ids[signature].size();

		random_shuffle(control_ids[signature].begin(), control_ids[signature].end());
		random_shuffle(event_ids[signature].begin(), event_ids[signature].end());
		for (int i = 0; i<ntake_ctrl; i++)
			selected.push_back(control_ids[signature][i]);
		for (int i = 0; i<ntake_ev; i++)
			selected.push_back(event_ids[signature][i]);
		n_ctrl += ntake_ctrl;
		n_events += ntake_ev;
		if (verbose) MLOG("%s : ntake_ctrl: %d/%d ntake_ev: %d/%d total size is %d\n", signature.c_str(), ntake_ctrl, control_ids[signature].size(),ntake_ev, event_ids[signature].size(), n_ctrl+n_events);
	}

	MLOG("Added %d controls, %d events, with a factor of %f\n", n_ctrl, n_events, n_events, factor);

	// Fill outSamples
	sort(selected.begin(), selected.end(), [](const pair<int, int> &v1, const pair<int, int> &v2) {return (v1.first < v2.first || (v1.first == v2.first && v1.second < v2.second)); });
	for (unsigned int i = 0; i < selected.size(); i++) {
		if (i == 0 || selected[i].first != selected[i - 1].first)
			outSamples.idSamples.push_back(MedIdSamples(inSamples.idSamples[selected[i].first].id));
		outSamples.idSamples.back().samples.push_back(inSamples.idSamples[selected[i].first].samples[selected[i].second]);
	}

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
//.......................................................................................
int MatchingSampleFilter::initHelpers(MedSamples& inSamples, MedRepository& rep) {

	// Helpers
	// Time unit from samples
	samplesTimeUnit = inSamples.time_unit;

	if (rep.dict.read_state == 0) {
		MWARN("WARNING Rep dictionary is empty\n");
	}

	// Age : either as a signal or using BYEAR
	if (isAgeRequired()) {
		if (med_rep_type.ageDirectlyGiven) {
			ageId = rep.dict.id("Age");
			if (ageId == -1) {
				MERR("Cannot find signalId for Age\n");
				return -1;
			}
		}
		else {
			byearId = rep.dict.id("BYEAR");
			if (byearId == -1) {
				MERR("Cannot find signalId for BYEAR\n");
				return -1;
			}
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

	return 0;
}

// Indexing of a single sample according to strata
// an index is a colon-separated string of bins per stratum
//.......................................................................................
int MatchingSampleFilter::getSampleSignature(MedSample& sample, MedRepository& rep, string& signature) {

	signature = "";
	for (auto& stratum : matchingStrata) {
		if (addToSampleSignature(sample, stratum, rep, signature) < 0)
			return -1;
	}

	return 0;
}

// add indexing of a single sample according to a single stratum to sample's index
//.......................................................................................
int MatchingSampleFilter::addToSampleSignature(MedSample& sample, matchingParams& stratum, MedRepository& rep, string& signature) {

	int len, age;
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
		if (med_rep_type.ageDirectlyGiven) {
			rep.uget(sample.id, ageId, usv);
			age = (int)usv.Val(0, 0);
		}
		else {
			int byear = (int)((SVal *)rep.get(sample.id, byearId, len))[0].val;
			age = med_time_converter.convert_times(samplesTimeUnit, MedTime::Date, sample.time) / 10000 - byear;
		}
		bin = (int)((float)age / stratum.resolution);
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
			if (usv.len>0 && usv.Time(usv.len-1) <= maxTime && usv.Time(usv.len-1) >= minTime)
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

// search for the optimal ratio between control/case samples
// the price of giving up 1 control is 1.0, the price of giving up 1 case is w 
//.......................................................................................
float MatchingSampleFilter::get_pairing_ratio(map<string, pair<int, int>> cnts, float w) {

	// Get all ratios
	vector<float> ratios;
	for (auto& rec : cnts) {
		int ev_cnt = rec.second.first;
		int ctrl_cnt = rec.second.second;
		if (ev_cnt > 0 && ctrl_cnt > 0)
			ratios.push_back((float)(ctrl_cnt) / (float)(ev_cnt));
	}
	sort(ratios.begin(), ratios.end());

	MLOG("min ratio %8.3f max ratio %8.3f\n", ratios[0], ratios[ratios.size() - 1]);

	// Find Optimal ratio - cnt1 and cnt2 are the overall number of samples we're giving up on 
	int opt_cnt1 = -1, opt_cnt2 = -1;
	float opt_r = 0;
	for (float r : ratios) {

		int cnt1 = 0, cnt2 = 0;
		for (auto it = cnts.begin(); it != cnts.end(); it++) {
			int ev_cnt = it->second.first;
			int ctrl_cnt = it->second.second;
			if (ev_cnt == 0.0)
				cnt1 += ctrl_cnt;
			else if (ctrl_cnt == 0.0)
				cnt2 += ev_cnt;
			else {
				double iratio = (float)(ctrl_cnt) / (float)(ev_cnt);
				if (iratio < r)
					// more events than we want, have to give up on some:
					cnt2 += (ev_cnt - (int)(ctrl_cnt / r + 0.5));
				else
					cnt1 += (ctrl_cnt - (int)(ev_cnt * r + 0.5));
			}
		}
		double current_price = cnt1 + w*cnt2;
		double opt_price = opt_cnt1 + w*opt_cnt2;

		if (opt_r == 0 || current_price < opt_price) {
			opt_r = r;
			opt_cnt1 = cnt1;
			opt_cnt2 = cnt2;
		}
		if (verbose) MLOG("ratio %8.3f price %8.3f (lose %d events, %d controls) \t opt ratio %8.3f lose %d events and %d controls = price of %8.3f \n", r, current_price, cnt2, cnt1, opt_r, opt_cnt2, opt_cnt1, opt_price);

	}

	return opt_r;
}

//Get all signals required  for matching
//.......................................................................................
void MatchingSampleFilter::get_required_signals(vector<string>& req_sigs)
{
	req_sigs.clear();
	if (isAgeRequired()) {
		if (med_rep_type.ageDirectlyGiven)
			req_sigs.push_back("Age");
		else
			req_sigs.push_back("BYEAR");
	}

	for (auto &s : matchingStrata) {
		if (s.signalName != "")
			req_sigs.push_back(s.signalName);
	}

	return ;

}

// (De)Serialization
//.......................................................................................
size_t MatchingSampleFilter::get_size() {
	return MedSerialize::get_size(matchingStrata, eventToControlPriceRatio, matchMaxRatio);
}

//.......................................................................................
size_t MatchingSampleFilter::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, matchingStrata, eventToControlPriceRatio, matchMaxRatio);
}

//.......................................................................................
size_t MatchingSampleFilter::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, matchingStrata, eventToControlPriceRatio, matchMaxRatio);
}

// (De)Serialization of matchingParams
//.......................................................................................
size_t matchingParams::get_size() {
	return MedSerialize::get_size(match_type,signalName,timeWindow,matchingTimeUnit,resolution);
}

//.......................................................................................
size_t matchingParams::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, match_type, signalName, timeWindow, matchingTimeUnit, resolution);
}

//.......................................................................................
size_t matchingParams::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, match_type, signalName, timeWindow, matchingTimeUnit, resolution);
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
		else if (field == "timeWindow") timeWindow = stoi(entry.second);
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
	windowTimeUnit = med_rep_type.windowTimeUnit;
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

// (De)Serialization
//.......................................................................................
size_t RequiredSignalFilter::get_size() {
	return MedSerialize::get_size(signalName,timeWindow,windowTimeUnit);
}

//.......................................................................................
size_t RequiredSignalFilter::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, signalName, timeWindow, windowTimeUnit);
}

//.......................................................................................
size_t RequiredSignalFilter::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, signalName, timeWindow, windowTimeUnit);
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

	for (int i=0; i<fields.size(); i++) {
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
		if (sig_name == "Age" && !med_rep_type.ageDirectlyGiven) {
			use_byear = 1;
			sig_id = rep.sigs.sid("BYEAR");
		}
		else
			sig_id = rep.sigs.sid(sig_name);
	}

	UniversalSigVec usv;

	rep.uget(sample.id, sig_id, usv);
	//MLOG("id %d sig_id %d len %d %f\n", sample.id, sig_id, usv.len, usv.Val(0));
	//MLOG("id %d sig_id %d len %d\n", sample.id, sig_id, usv.len);

	if (usv.len == 0 && min_Nvals > 0) return 0;
	if (min_Nvals <= 0) return 1;

	// Special handling of age through byear
	if (use_byear) {
		int year = 1900 + med_time_converter.convert_times(med_rep_type.basicTimeUnit, MedTime::Years, sample.time);
		int age = year - (int)usv.Val(0);
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
		for (int i=0; i<usv.len; i++) {

			// check if in relevant window
			int i_time = usv.Time(i, time_channel);
			int i_time_converted = med_time_converter.convert_times(usv.time_unit(), win_time_unit, i_time);
			int dtime = ref_time - i_time_converted;
			//MLOG("id %d i_time %d %f %d time %d %d dtime %d win %d %d\n", sample.id, i_time, usv.Val(i, val_channel), i_time_converted, sample.time, ref_time, dtime, win_from, win_to);
			if (dtime < win_from) break;
			if (dtime <= win_to) {
				// in relevant time window, checking the value range
				float i_val = usv.Val(i, val_channel);
				//MLOG("i %d id %d i_val %f min %f max %f minNvals %d\n", i, sample.id, i_val, min_val, max_val, min_Nvals);
				if (i_val >= min_val && i_val <= max_val) {
					nvals++;
					if (nvals >= min_Nvals)
						return 1;
				}
			}
		}
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
		if (m.first == "bfilter") {
			vector<string> fields;
			boost::split(fields, m.second, boost::is_any_of("+"));
			for (auto &f : fields) {
				BasicFilteringParams bfp;
				bfp.init_from_string(f);
				bfilters.push_back(bfp);
				if (bfp.sig_name == "Age" && !med_rep_type.ageDirectlyGiven)
					req_sigs.push_back("BYEAR");
				else
					req_sigs.push_back(bfp.sig_name);
			}
		}
		
	}

	return 0;
}

//.......................................................................................
void BasicSampleFilter::get_required_signals(vector<string> &reqs)
{
	if (req_sigs.size() == 0) {
		req_sigs.clear();
		for (auto &bf : bfilters) {
			if (bf.sig_name == "Age" && !med_rep_type.ageDirectlyGiven)
				req_sigs.push_back("BYEAR");
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
				for (auto &bf : bfilters) {
					if (!bf.test_filter(in_s, rep, winsTimeUnit)) {
						take_it = 0;
						break;
					}
				}

				if (take_it) id_out.samples.push_back(in_s);
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

	for (int i=0; i<fields.size(); i++) {
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
			boost:split(svals, fields[++i], boost::is_any_of(","));	
			for (auto &s : svals) allowed_values.insert(stof(s));
		}
	}


	return 0;
}

// Test filtering criteria
// Returns one of the codes defined as static in the h file
//.......................................................................................
int SanitySimpleFilter::test_filter(MedSample &sample, MedRepository &rep, int &nvals, int &noutliers)
{
	//MLOG("SanitySimpleFilter::test_filter() ==> id %d sig_id %d %s time %d\n", sample.id, sig_id, sig_name.c_str(), sample.time);
	if (sig_id < 0) {
		if (boost::iequals(sig_name,"Age")) {
			sig_id = 0;
			byear_id = rep.sigs.sid("BYEAR");
			if (byear_id < 0) {
				MWARN("WARNING: !!!! ===> Using SanitySimpleFilter for age but without BYEAR... Are you using a repository with an AGE signal??\n");
			}
		} else
			sig_id = rep.sigs.sid(sig_name);
	}
	if (sig_id < 0)
		return SanitySimpleFilter::Signal_Not_Valid;

	// Age case
	if (sig_id == 0) {
		// TBD: Must make this work also for the cases in which Age is given as a signal
		if (byear_id > 0) {
			// calculate using byear
			int len;
			float y = 1900 + (float)med_time_converter.convert_times(samples_time_unit, MedTime::Years, sample.time);
			SVal *sv = (SVal *)rep.get(sample.id, byear_id, len);
			if (len > 0) {
				float age = y - sv[0].val;
				//MLOG("====> AGE : byear %f y %f time %d : age %f min_val %f max_val %f\n", sv[0].val, y, sample.time, age, min_val, max_val);
				if (age < min_val || age > max_val)
					return SanitySimpleFilter::Failed_Age;
			}
			else
				return SanitySimpleFilter::Failed_Age_No_Byear;
		} else
			return SanitySimpleFilter::Failed_Age_No_Byear;
	}

	if (sig_id > 0) {
		// regular signal case
		if (section_id < 0 && sig_id > 0) {
			section_id = rep.dict.section_id(sig_name);
		}

		UniversalSigVec usv;

		rep.uget(sample.id, sig_id, usv);
		//MLOG("id %d sig_id %d len %d %f\n", sample.id, sig_id, usv.len, usv.Val(0));
		//MLOG("id %d sig_id %d len %d\n", sample.id, sig_id, usv.len);

		nvals = 0;
		noutliers = 0;
		int n_not_in_dict = 0;
		int n_not_allowed = 0;
		int n_left = 0;
		if (usv.len == 0 && min_Nvals > 0) return SanitySimpleFilter::Failed_Min_Nvals;

		if (usv.n_time_channels() == 0) {
			// timeless signal

			nvals = usv.len;
			for (int i=0; i<usv.len; i++) {
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
			for (int i=0; i<usv.len; i++) {

				// check if in relevant window
				int i_time = usv.Time(i, time_channel);
				int i_time_converted = med_time_converter.convert_times(usv.time_unit(), win_time_unit, i_time);
				int dtime = ref_time - i_time_converted;
				//MLOG("id %d i_time %d %f %d time %d %d dtime %d win %d %d\n", sample.id, i_time, usv.Val(i, val_channel), i_time_converted, sample.time, ref_time, dtime, win_from, win_to);
				if (dtime < win_from) break;
				if (dtime <= win_to) {
					nvals++;
					// in relevant time window, checking the value range
					float i_val = usv.Val(i, val_channel);
					if (i_val < min_val || i_val > max_val) noutliers++;
					else n_left++;
					//MLOG("i %d id %d i_val %f min %f max %f minNvals %d nvals %d noutliers %d\n", i, sample.id, i_val, min_val, max_val, min_Nvals, nvals, noutliers);

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
		}

		//	MLOG("###>>> id %d time %d sig %s (len %d) : %f : min %d max %d maxout %d : nvals %d noutliers %d not_in_dict %d not_allowed %d\n",
		//		sample.id, sample.time, sig_name.c_str(), usv.len, usv.Val(0), min_Nvals, max_Nvals, max_outliers, nvals, noutliers, n_not_in_dict , n_not_allowed);

		if (min_Nvals >= 0 && nvals < min_Nvals) return SanitySimpleFilter::Failed_Min_Nvals;
		if (max_Nvals >= 0 && nvals > max_Nvals) return SanitySimpleFilter::Failed_Max_Nvals;
		if (max_outliers >= 0 && noutliers > max_outliers) return SanitySimpleFilter::Failed_Outliers;
		if (values_in_dictionary && ((n_not_in_dict > 0) || section_id < 0)) return SanitySimpleFilter::Failed_Dictionary_Test;
		if ((allowed_values.size() > 0) && (n_not_allowed > 0)) return SanitySimpleFilter::Failed_Allowed_Values;
		if (min_left >=0 && n_left < min_left) return SanitySimpleFilter::Failed_Not_Enough_Non_Outliers_Left;
	}

	return SanitySimpleFilter::Passed;

	return 0;
}
