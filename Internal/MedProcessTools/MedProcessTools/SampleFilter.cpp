#include "SampleFilter.h"
#include "Logger/Logger/Logger.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/MedPidRepository.h"

#define LOCAL_SECTION LOG_SMPL_FILTER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// SampleFilter
//=======================================================================================
// Filter Types
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
	else
		return SMPL_FILTER_LAST;
}

// Initialization
//.......................................................................................
SampleFilter* SampleFilter::make_filter(string filter_name) {

	return make_filter(sample_filter_name_to_type(filter_name));
}

//.......................................................................................
SampleFilter * SampleFilter::make_filter(string filter_name, string init_string) {

	return make_filter(sample_filter_name_to_type(filter_name), init_string);
}

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
	else
		return NULL;

}

//.......................................................................................
SampleFilter * SampleFilter::make_filter(SampleFilterTypes filter_type, string init_string) {

	SampleFilter *newSampleFilter = make_filter(filter_type);
	newSampleFilter->init_from_string(init_string);
	return newSampleFilter;
}

//.......................................................................................

// (De)Serialize
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

// Filter
//.......................................................................................
int SampleFilter::filter(MedSamples& samples) {

	MedSamples out_samples;
	int rc = filter(samples, out_samples);

	if (rc == 0)
		samples = out_samples;

	return rc;
}

// Filter
//.......................................................................................
int SampleFilter::filter(MedRepository& rep, MedSamples& samples) {

	MedSamples out_samples;
	int rc = filter(rep, samples, out_samples);

	if (rc == 0)	
		samples = out_samples;
	
	return rc;
}


//=======================================================================================
// BasicTrainFilter
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
// BasicTestFilter
//=======================================================================================
// Filter
//.......................................................................................
int BasicTestFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {

	// Take them all
	outSamples = inSamples;

	return 0;
}

//=======================================================================================
// OutlierSampleFilter
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

// Learn
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

// Learn
//.......................................................................................
int OutlierSampleFilter::iterativeLearn(MedSamples& samples) {

	// Get all values
	vector<float> values;
	get_values(samples,values);

	return get_iterative_min_max(values);
}

// Learn
//.......................................................................................
int OutlierSampleFilter::quantileLearn(MedSamples& samples) {

	// Get all values
	vector<float> values;
	get_values(samples,values);

	return get_quantile_min_max(values);
}

// Utility for learning
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

// Init
//.......................................................................................
int MatchingSampleFilter::init(map<string, string>& mapper) {

	vector<string> strata;

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "priceRatio") eventToCasePriceRatio = stof(entry.second);
		else if (field == "maxRatio") matchMaxRatio = stof(entry.second);
		else if (field == "strata") {
			boost::split(strata, entry.second, boost::is_any_of(":"));
			for (string& stratum : strata) addMatchingStrata(stratum);
		} else
			MLOG("Unknonw parameter \'%s\' for MatchingSampleFilter\n", field.c_str());
	}

	return 0;
}

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
		newStrata.resolution = (fields.size() > 1) ? stof(fields[1]) : 1.0;

	}
	else if (fields[0] == "time") {
		if (fields.size() > 3 || fields.size() < 2) {
			MERR("Wrong number of features for matching strata\n");
			return -1;
		}

		newStrata.match_type = SMPL_MATCH_TIME;
		newStrata.matchingTimeUnit = med_time_converter.string_to_type(fields[1]);
		newStrata.resolution = (fields.size() > 2) ? stof(fields[2]) : 1.0;
	}
	else if (fields[0] == "signal") {
		if (fields.size() > 5 || fields.size() < 2) {
			MERR("Wrong number of features for matching strata\n");
			return -1;
		}

		newStrata.match_type = SMPL_MATCH_SIGNAL;
		newStrata.signalName = fields[1];
		newStrata.resolution = (fields.size() > 2) ? stof(fields[2]) : 1.0;
		newStrata.timeWindow = (fields.size() > 3) ? stof(fields[3]) : 1.0;
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

// Filter
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
	float opt_factor = get_pairing_ratio(cnts, eventToCasePriceRatio);
	float factor = opt_factor;
	if (factor > matchMaxRatio) {
		MLOG("updating factor {%8.3f} to maxFactor=%.3f\n", factor, matchMaxRatio);
		factor = matchMaxRatio;
	}
	else if (factor < 1 / matchMaxRatio) {
		MLOG("updating factor {%8.3f} to 1/(maxFactor=%.3f)\n", factor, matchMaxRatio);
		factor = 1/matchMaxRatio;
	}
//	MLOG("opt ratio is %8.3f (effective factor = %8.3f)\n", opt_factor, factor);


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
//		MLOG("%s : ntake_ctrl: %d/%d ntake_ev: %d/%d total size is %d\n", signature.c_str(), ntake_ctrl, control_ids[signature].size(),ntake_ev, event_ids[signature].size(), n_ctrl+n_events);
	}

//	MLOG("Added %d controls, %d events, with a factor of %f\n", n_ctrl, n_events, n_events, factor);

	// Fill outSamples
	sort(selected.begin(), selected.end(), [](const pair<int, int> &v1, const pair<int, int> &v2) {return (v1.first < v2.first || (v1.first == v2.first && v1.second < v2.second)); });
	for (unsigned int i = 0; i < selected.size(); i++) {
		if (i == 0 || selected[i].first != selected[i - 1].first)
			outSamples.idSamples.push_back(MedIdSamples(inSamples.idSamples[selected[i].first].id));
		outSamples.idSamples.back().samples.push_back(inSamples.idSamples[selected[i].first].samples[selected[i].second]);
	}

	return 0;

}

//.......................................................................................
int MatchingSampleFilter::_filter(MedSamples& inSamples, MedSamples& outSamples) {

	if (isRepRequired()) {
		MERR("Cannot perform required matching without repository\n");
		return -1;
	}
	else {
		MedRepository dummyRep;
		return _filter(dummyRep, inSamples, outSamples);
	}
}

// Utilities
//.......................................................................................
bool MatchingSampleFilter::isRepRequired() {

	for (auto& strata : matchingStrata) {
		if (strata.match_type == SMPL_MATCH_AGE || strata.match_type == SMPL_MATCH_SIGNAL)
			return true;
	}

	return false;
}

//.......................................................................................
bool MatchingSampleFilter::isAgeRequired() {

	for (auto& strata : matchingStrata) {
		if (strata.match_type == SMPL_MATCH_AGE)
			return true;
	}

	return false;
}

//.......................................................................................
int MatchingSampleFilter::initHelpers(MedSamples& inSamples, MedRepository& rep) {

	// Helpers
	samplesTimeUnit = inSamples.time_unit;

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
//.......................................................................................
int MatchingSampleFilter::getSampleSignature(MedSample& sample, MedRepository& rep, string& signature) {

	signature = "";
	for (auto& stratum : matchingStrata) {
		if (addToSampleSignature(sample, stratum, rep, signature) < 0)
			return -1;
	}

	return 0;
}

//.......................................................................................
int MatchingSampleFilter::addToSampleSignature(MedSample& sample, matchingParams& stratum, MedRepository& rep, string& signature) {

	int len, age;
	UniversalSigVec usv;
	int bin;

	if (stratum.match_type == SMPL_MATCH_TIME) {
		int time = med_time_converter.convert_times(samplesTimeUnit, stratum.matchingTimeUnit, sample.time);
		int bin = (int)(time / stratum.resolution);
		signature += to_string(bin) + ":";
	}
	else if (stratum.match_type == SMPL_MATCH_AGE) {
		if (med_rep_type.ageDirectlyGiven) {
			rep.uget(sample.id, ageId, usv);
			age = (int)usv.Val(0, 0);
		}
		else {
			int byear = (int)((SVal *)rep.get(sample.id, byearId, len))[0].val;
			age = med_time_converter.convert_times(samplesTimeUnit, MedTime::Date, sample.time) / 10000 - byear;
		}
		bin = age / stratum.resolution;
		signature += to_string(bin) + ":";
	}
	else if (stratum.match_type == SMPL_MATCH_SIGNAL) {
		rep.uget(sample.id, stratum.signalId, usv);
		if (!stratum.isTimeDependent) {
			bin = (int)(usv.Val(0) / stratum.resolution);
			signature += to_string(bin) + ":";
		}
		else {
			int target = med_time_converter.convert_times(samplesTimeUnit, stratum.windowTimeUnit, sample.time);
			int maxTime = med_time_converter.convert_times(samplesTimeUnit, stratum.signalTimeUnit, sample.time);
			int minTime = med_time_converter.convert_times(stratum.windowTimeUnit, stratum.signalTimeUnit, target - stratum.timeWindow);
//			MLOG("units = %d/%d/%d time = %d Target = %d min = %d\n", samplesTimeUnit, stratum.signalTimeUnit, stratum.windowTimeUnit, sample.time, maxTime, minTime);

			string tempSignature = "NULL";
			for (int idx = 0; idx < usv.len; idx++) {
				if (usv.Time(idx) > maxTime) {
					if (idx > 0 && usv.Time(idx - 1) >= minTime)
						tempSignature = to_string((int)(0.001 + usv.Val(idx - 1) / stratum.resolution));
					break;
				}
			}

			if (usv.len>0 && usv.Time(usv.len-1) <= maxTime && usv.Time(usv.len-1) >= minTime)
				tempSignature = to_string((int)(0.001 + usv.Val(usv.len - 1) / stratum.resolution));

			signature += tempSignature + ":";  
		}
	}

	return 0;
}

// search for the optimal ratio between control/event samples
// the price of giving up 1 control is 1.0, the price of giving up 1 event is w 
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

//	MLOG("min ratio %8.3f max ratio %8.3f\n", ratios[0], ratios[ratios.size() - 1]);

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
//		MLOG("ratio %8.3f price %8.3f \t opt ratio %8.3f lose %d events and %d controls = price of %8.3f \n", r, current_price, opt_r, opt_cnt2, opt_cnt1, opt_price);

	}

	return opt_r;
}


// (De)Serialization
//.......................................................................................
size_t MatchingSampleFilter::get_size() {
	return MedSerialize::get_size(matchingStrata, eventToCasePriceRatio, matchMaxRatio);
}

//.......................................................................................
size_t MatchingSampleFilter::serialize(unsigned char *blob) {
	return MedSerialize::serialize(blob, matchingStrata, eventToCasePriceRatio, matchMaxRatio);
}

//.......................................................................................
size_t MatchingSampleFilter::deserialize(unsigned char *blob) {
	return MedSerialize::deserialize(blob, matchingStrata, eventToCasePriceRatio, matchMaxRatio);
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
// RequiredSignalFilter
//=======================================================================================

// Init
//.......................................................................................
int RequiredSignalFilter::init(map<string, string>& mapper) {

	vector<string> strata;

	for (auto entry : mapper) {
		string field = entry.first;

		if (field == "signalName") signalName = entry.second;
		else if (field == "timeWindow") timeWindow = stof(entry.second);
		else if (field == "timeUnit") windowTimeUnit = med_time_converter.string_to_type(entry.second);
		else
			MLOG("Unknonw parameter \'%s\' for RequiredSampleFilter\n", field.c_str());
	}

	return 0;
}

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