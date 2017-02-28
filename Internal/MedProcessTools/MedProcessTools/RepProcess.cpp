#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "RepProcess.h"

//=======================================================================================
// RepProcessors
//=======================================================================================
// Processors types
RepProcessorTypes rep_processor_name_to_type(const string& processor_name) {

	if (processor_name == "multi_processor" || processor_name == "multi")
		return REP_PROCESS_MULTI;
	else if (processor_name == "basic_outlier_cleaner" || processor_name == "basic_cln")
		return REP_PROCESS_BASIC_OUTLIER_CLEANER;
	else if (processor_name == "nbrs_outlier_cleaner" || processor_name == "nbrs_cln")
		return REP_PROCESS_NBRS_OUTLIER_CLEANER;
	else
		return REP_PROCESS_LAST;
}

//.......................................................................................
// create from params : must have type parameter
RepProcessor *RepProcessor::create_processor(string &params)
{
	string rp_type;
	get_single_val_from_init_string(params, "rp_type", rp_type);
	return (make_processor(rp_type, params));
}

// Initialization
//.......................................................................................
RepProcessor * RepProcessor::make_processor(string processor_name) {

	return make_processor(rep_processor_name_to_type(processor_name));
}

//.......................................................................................
RepProcessor * RepProcessor::make_processor(string processor_name, string init_string) {

	return make_processor(rep_processor_name_to_type(processor_name), init_string);
}

//.......................................................................................
RepProcessor * RepProcessor::make_processor(RepProcessorTypes processor_type) {

	if (processor_type == REP_PROCESS_MULTI)
		return new RepMultiProcessor;
	else if (processor_type == REP_PROCESS_BASIC_OUTLIER_CLEANER)
		return new RepBasicOutlierCleaner;
	else if (processor_type == REP_PROCESS_NBRS_OUTLIER_CLEANER)
		return new RepNbrsOutlierCleaner;
	else
		return NULL;

}

//.......................................................................................
RepProcessor * RepProcessor::make_processor(RepProcessorTypes processor_type, string init_string) {

	//MLOG("Processor type is %d\n", (int)processor_type);
	RepProcessor *newRepProcessor = make_processor(processor_type);
	newRepProcessor->init_from_string(init_string);
	return newRepProcessor;
}

// Applying
//.......................................................................................
int RepProcessor::apply(PidDynamicRec& rec, MedIdSamples& samples) {

	vector<int> time_points(samples.samples.size());
	for (unsigned int i = 0; i < time_points.size(); i++)
		time_points[i] = samples.samples[i].date;

	return apply(rec, time_points);
}

//.......................................................................................
int RepProcessor::apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignalIds) {

	for (int signalId : neededSignalIds) {
		if (is_signal_affected(signalId))
			return apply(rec, time_points);
	}

	return 0;
}

//.......................................................................................
int RepProcessor::apply(PidDynamicRec& rec, MedIdSamples& samples, vector<int>& neededSignalIds) {

	vector<int> time_points(samples.samples.size());
	for (unsigned int i = 0; i < time_points.size(); i++)
		time_points[i] = samples.samples[i].date;

	return apply(rec, time_points, neededSignalIds);
}

// Required signals
//.......................................................................................
void RepProcessor::get_required_signal_ids(MedDictionarySections& dict) {

	req_signal_ids.resize(req_signals.size());

	for (unsigned int i = 0; i < req_signals.size(); i++)
		req_signal_ids[i] = dict.id(req_signals[i]);
}

//.......................................................................................
void RepProcessor::get_required_signal_ids(unordered_set<int>& signalIds, MedDictionarySections& dict) {

	if (req_signal_ids.empty())
		get_required_signal_ids(dict);

	for (int signalId : req_signal_ids)
		signalIds.insert(signalId);
}

// Affected signals
//.......................................................................................
void RepProcessor::get_affected_signal_ids(MedDictionarySections& dict) {

	for (string signalName : aff_signals)
		aff_signal_ids.insert(dict.id(signalName));
}

// (De)Serialize
//.......................................................................................
size_t RepProcessor::get_processor_size() {
	return sizeof(processor_type) + get_size();
}

//.......................................................................................
size_t RepProcessor::processor_serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &processor_type, sizeof(RepProcessorTypes)); ptr += sizeof(RepProcessorTypes);
	ptr += serialize(blob + ptr);

	return ptr;
}

//=======================================================================================
// RepMultiProcessor
//=======================================================================================
// Required Signals
//.......................................................................................
void RepMultiProcessor::get_required_signal_ids(MedDictionarySections& dict) {

	req_signal_ids.clear();
	for (auto& processor : processors) {
		req_signal_ids.clear();
		processor->get_required_signal_ids(dict);
		req_signal_ids.insert(req_signal_ids.end(), processor->req_signal_ids.begin(), processor->req_signal_ids.end());
	}
}

// Affected Signals
//.......................................................................................
void RepMultiProcessor::get_affected_signal_ids(MedDictionarySections& dict) {

	aff_signal_ids.clear();
	for (auto& processor : processors) {
		processor->aff_signal_ids.clear();
		processor->get_affected_signal_ids(dict);

		processor->aff_signal_ids.clear();
		for (int signalId : processor->aff_signal_ids)
			aff_signal_ids.insert(signalId);
	}
}

// Signal ids
//.......................................................................................
void RepMultiProcessor::get_signal_ids(MedDictionarySections& dict) {

	for (auto& processor : processors)
		processor->get_signal_ids(dict);

}

// Learning
//.......................................................................................
int RepMultiProcessor::Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) {

	vector<int> rc(processors.size(), 0);
#pragma omp parallel for
	for (int j=0; j<processors.size(); j++) {
	//for (auto& cleaner : cleaners) {
		rc[j] = processors[j]->Learn(rep, ids, prev_processors);
	}

	for (int r : rc) if (r<0) return -1;
	return 0;
}

// Apply
//.......................................................................................
int RepMultiProcessor::apply(PidDynamicRec& rec, vector<int>& time_points) {

	vector<int> rc(processors.size(),0);
#pragma omp parallel for
	for (int j=0; j<processors.size(); j++) {
		rc[j] = processors[j]->apply(rec, time_points);
	}

	for (int r : rc) if (r<0) return -1;
	return 0;
}

//.......................................................................................
int RepMultiProcessor::apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignalIds) {

	vector<int> rc(processors.size(), 0);

#pragma omp parallel for
	for (int j = 0; j<processors.size(); j++) {
		rc[j] = processors[j]->apply(rec, time_points, neededSignalIds);
	}

	for (int r : rc) if (r<0) return -1;
	return 0;
}

// Add processors
//.......................................................................................
void  RepMultiProcessor::add_processors_set(RepProcessorTypes type, vector<string>& signals) {

	for (string& signal : signals) {
		RepProcessor *processor = RepProcessor::make_processor(type); 
		processor->set_signal(signal);
		processors.push_back(processor);
	}

}

//.......................................................................................
void  RepMultiProcessor::add_processors_set(RepProcessorTypes type, vector<string>& signals, string init_string) {

	for (string& signal : signals) {
		RepProcessor *processor = RepProcessor::make_processor(type, init_string);
		processor->set_signal(signal);
		processors.push_back(processor);
	}

}

// (De)Serialization
//.......................................................................................
size_t RepMultiProcessor::get_size() {

	size_t size = sizeof(int); // Number of cleaners
	for (auto& processor : processors)
		size += processor->get_processor_size();

	return size;
}

//.......................................................................................
size_t RepMultiProcessor::serialize(unsigned char *blob) {

	size_t ptr = 0;
	int nProcessors = (int)processors.size();
	memcpy(blob + ptr, &nProcessors, sizeof(int)); ptr += sizeof(int);

	for (auto& processor : processors)
		ptr += processor->processor_serialize(blob + ptr);

	return ptr;
}

//.......................................................................................
size_t RepMultiProcessor::deserialize(unsigned char *blob) {

	size_t ptr = 0;
	int nProcessors;
	memcpy(&nProcessors, blob + ptr, sizeof(int)); ptr += sizeof(int);

	processors.resize(nProcessors);
	for (int i = 0; i < nProcessors; i++) {
		RepProcessorTypes type; 
		memcpy(&type, blob + ptr, sizeof(RepProcessorTypes)); ptr += sizeof(RepProcessorTypes);
		processors[i] = RepProcessor::make_processor(type);
		ptr += processors[i]->deserialize(blob + ptr);
	}

	return ptr;
}

//=======================================================================================
// BasicOutlierCleaner
//=======================================================================================
//.......................................................................................
int RepBasicOutlierCleaner::init(map<string, string>& mapper) 
{ 
	init_defaults(); 

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "signal") { signalName = entry.second; req_signals.push_back(signalName); }
		else if (field == "time_channel") time_channel = stoi(entry.second);
		else if (field == "val_channel") time_channel = stoi(entry.second);
	}

	return MedValueCleaner::init(mapper); 
}

// Learn bounds
//.......................................................................................
int RepBasicOutlierCleaner::Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_cleaners) {

	if (params.type == VAL_CLNR_ITERATIVE)
		return iterativeLearn(rep, ids, prev_cleaners);
	else if (params.type == VAL_CLNR_QUANTILE)
		return quantileLearn(rep, ids, prev_cleaners);
	else {
		MERR("Unknown cleaning type %d\n",params.type);
		return -1;
	}
}

//.......................................................................................
int RepBasicOutlierCleaner::iterativeLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_cleaners) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Get all values
	vector<float> values;
	get_values(rep, ids, signalId, time_channel, val_channel, params.range_min, params.range_max, values, prev_cleaners);

	int rc =  get_iterative_min_max(values);
	//print();
	return rc;
}

//.......................................................................................
int RepBasicOutlierCleaner::quantileLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_cleaners) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Get all values
	vector<float> values;
	get_values(rep, ids, signalId, time_channel, val_channel, params.range_min, params.range_max, values, prev_cleaners);

	return get_quantile_min_max(values);
}

// Clean
//.......................................................................................
int  RepBasicOutlierCleaner::apply(PidDynamicRec& rec, vector<int>& time_points) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int len;
	int iver = (int)time_points.size() - 1;

	while (iver >= 0) {
		// Check all versions that points to the same location
		int jver = iver - 1;
		while (jver >= 0 && rec.versions_are_the_same(signalId, iver, jver))
			jver--;
		jver++;

		// we now know that versions [jver ... iver] are all pointing to the same version
		// hence we need to clean version jver, and then make sure all these versions point to it

		// Clean 
		rec.uget(signalId, jver, rec.usv); // get into the internal usv obeject - this statistically saves init time

		len = rec.usv.len; 
		vector<int> remove(len);
		vector<pair<int, float>> change(len);
		int nRemove = 0, nChange = 0;

		// Collect
		for (int i = 0; i < len; i++) {
			int itime = rec.usv.Time(i, time_channel);
			float ival = rec.usv.Val(i, val_channel);

			if (itime > time_points[jver])	break;

			if (params.doRemove && (ival < removeMin - NUMERICAL_CORRECTION_EPS || ival > removeMax + NUMERICAL_CORRECTION_EPS))
				remove[nRemove++] = i;
			else if (params.doTrim) {
				if (ival < trimMin)
					change[nChange++] = pair<int, float>(i, trimMin);
				else if (ival > trimMax)
					change[nChange++] = pair<int, float>(i, trimMax);
			}
		}

		// Apply removals + changes
		
		change.resize(nChange);
		remove.resize(nRemove);
		if (rec.update(signalId, jver, val_channel, change, remove) < 0)
			return -1;

		while (iver > jver) {
			rec.point_version_to(signalId, jver, iver);
			iver--;
		}
		iver--;
	}

	return 0;

}

// (De)Serialization
//.......................................................................................
size_t RepBasicOutlierCleaner::get_size() {

	size_t size = 0;

	// signalName
	size += sizeof(size_t);
	size += signalName.length() + 1;
	
	size += 2*sizeof(int); // time_channel, val_channel

	size += sizeof(int); // int take_log
	size += sizeof(float); // float missing value

	size += sizeof(bool); //  bool doTrim;
	size += 2 * sizeof(float); // float trimMax, trimMin;

	size += sizeof(bool); //  bool doRemove;
	size += 2 * sizeof(float); // float removeMax, removeMin;

	return size;
}

extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t RepBasicOutlierCleaner::serialize(unsigned char *blob) {

	size_t ptr = 0;
	// SignalName
	size_t nameLen = signalName.length();
	assert(nameLen < MAX_NAME_LEN);

	strcpy(signalName_c, signalName.c_str());

	memcpy(blob + ptr, &nameLen, sizeof(size_t)); ptr += sizeof(size_t);
	memcpy(blob + ptr, signalName_c, nameLen + 1); ptr += nameLen + 1;
	memcpy(blob + ptr, &time_channel, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &val_channel, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &params.take_log, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &params.missing_value, sizeof(float)); ptr += sizeof(int);
	memcpy(blob + ptr, &params.doTrim, sizeof(bool)); ptr += sizeof(bool);
	memcpy(blob + ptr, &trimMax, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &trimMin, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &params.doRemove, sizeof(bool)); ptr += sizeof(bool);
	memcpy(blob + ptr, &removeMax, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &removeMin, sizeof(float)); ptr += sizeof(float);

	return ptr;
}

//.......................................................................................
size_t RepBasicOutlierCleaner::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	// FeatureName
	size_t nameLen;
	memcpy(&nameLen, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);
	assert(nameLen < MAX_NAME_LEN);

	memcpy(signalName_c, blob + ptr, nameLen + 1); ptr += nameLen + 1;
	signalName = signalName_c;
	memcpy(&time_channel, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&val_channel, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&params.take_log, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&params.missing_value, blob + ptr, sizeof(float)); ptr += sizeof(int);
	memcpy(&params.doTrim, blob + ptr, sizeof(bool)); ptr += sizeof(bool);
	memcpy(&trimMax, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&trimMin, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&params.doRemove, blob + ptr, sizeof(bool)); ptr += sizeof(bool);
	memcpy(&removeMax, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&removeMin, blob + ptr, sizeof(float)); ptr += sizeof(float);

	return ptr;
}

//.......................................................................................
void RepBasicOutlierCleaner::print()
{
	MLOG("BasicOutlierCleaner: signal: %d : doTrim %d trimMax %f trimMin %f : doRemove %d : removeMax %f removeMin %f\n",
		signalId, params.doTrim, trimMax, trimMin, params.doRemove, removeMax, removeMin);
}

//=======================================================================================
// NbrsOutlierCleaner
//=======================================================================================
//.......................................................................................
int RepNbrsOutlierCleaner::init(map<string, string>& mapper)
{
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "signal") { signalName = entry.second; req_signals.push_back(signalName); }
		else if (field == "time_channel") time_channel = stoi(entry.second);
		else if (field == "val_channel") time_channel = stoi(entry.second);
		else if (field == "nbr_time_unit") nbr_time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "nbr_time_width") nbr_time_width = stoi(entry.second);

	}

	return MedValueCleaner::init(mapper);
}

// Learn bounds
//.......................................................................................
int RepNbrsOutlierCleaner::Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_cleaners) {

	if (params.type == VAL_CLNR_ITERATIVE)
		return iterativeLearn(rep, ids, prev_cleaners);
	else if (params.type == VAL_CLNR_QUANTILE)
		return quantileLearn(rep, ids, prev_cleaners);
	else {
		MERR("Unknown cleaning type %d\n", params.type);
		return -1;
	}
}

//.......................................................................................
int RepNbrsOutlierCleaner::iterativeLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_cleaners) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Get all values
	vector<float> values;
	get_values(rep, ids, signalId, time_channel, val_channel, params.range_min, params.range_max, values, prev_cleaners);

	return get_iterative_min_max(values);
}

//.......................................................................................
int RepNbrsOutlierCleaner::quantileLearn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_cleaners) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Get all values
	vector<float> values;
	get_values(rep, ids, signalId, time_channel, val_channel, params.range_min, params.range_max, values, prev_cleaners);

	return get_quantile_min_max(values);
}

// Clean
//.......................................................................................
int  RepNbrsOutlierCleaner::apply(PidDynamicRec& rec, vector<int>& time_points) {
	
	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int len;

	for (int iver = 0; iver < time_points.size(); iver++) {

		rec.uget(signalId, iver, rec.usv);

		len = rec.usv.len;
		vector<int> remove(len);
		vector<pair<int, float>> change(len);
		int nRemove = 0, nChange = 0;

		// Clean 
		int verLen = 0;
		vector<int> candidates(len, 0);
		vector<int> removed(len, 0);
		for (int i = 0; i < len; i++) {
			int itime = rec.usv.Time(i, time_channel);

			if (itime > time_points[iver]) {
				verLen = i - 1;
				break;
			}

			float ival = rec.usv.Val(i, val_channel);

			// Remove ?
			if (params.doRemove && (ival < removeMin - NUMERICAL_CORRECTION_EPS || ival > removeMax + NUMERICAL_CORRECTION_EPS)) {
				remove[nRemove++] = i;
				removed[i] = 1;
			}
			else if (params.doTrim) {
				if (ival < trimMin - NUMERICAL_CORRECTION_EPS) {
					candidates[i] = -1;
				}
				else if (ival > trimMax + NUMERICAL_CORRECTION_EPS) {
					candidates[i] = 1;
				}
			}
		}

		// Check candidates
		for (int i = 0; i < verLen; i++) {
			if (candidates[i] != 0) {
				int dir = candidates[i];

				// Get weighted values from neighbours
				double sum = 0, norm = 0;
				double priorSum = 0, priorNorm = 0;
				double postSum = 0, postNorm = 0;

				int time_i = rec.usv.TimeU(i, nbr_time_unit);

				for (int j = 0; j < verLen; j++) {

					if (j != i && !removed[j]) {
						int diff = abs(rec.usv.TimeU(j, nbr_time_unit) - time_i) / nbr_time_width;
						double w = 1.0 / (diff + 1);

						float jval = rec.usv.Val(j, val_channel);

						sum += w * jval;
						norm += w;

						if (j > i) {
							postSum += w * jval;
							postNorm += w;
						}
						else {
							priorSum += w * jval;
							priorNorm += w;
						}
					}
				}

				// Check it up
				int found_nbr = 0;
				if (norm > 0) {
					double win_val = sum / norm;
					if ((dir == 1 && win_val > nbrsMax) || (dir == -1 && win_val < nbrsMin))
						found_nbr = 1;
				}

				if (!found_nbr && priorNorm > 0) {
					double win_val = priorSum / priorNorm;
					if ((dir == 1 && win_val > nbrsMax) || (dir == -1 && win_val < nbrsMin))
						found_nbr = 1;
				}

				if (!found_nbr && postNorm > 0) {
					double win_val = postSum / postNorm;
					if ((dir == 1 && win_val > nbrsMax) || (dir == -1 && win_val < nbrsMin))
						found_nbr = 1;
				}

				// Should we clip ?
				if (!found_nbr) {
					float cval = (dir == 1) ? trimMax : trimMin;
					change[nChange++] = pair<int, float>(i, cval);
				}
			}
		}


		// Apply removals + changes
		change.resize(nChange);
		remove.resize(nRemove);
		if (rec.update(signalId, iver, val_channel, change, remove) < 0)
			return -1;
	}

	return 0;

}


// (De)Serialization
//.......................................................................................
size_t RepNbrsOutlierCleaner::get_size() {

	size_t size = 0;

	// signalName
	size += sizeof(size_t);
	size += signalName.length() + 1;

	size += 2*sizeof(int); // time_channel, val_channel
	size += 2*sizeof(int); // nbr_time_unit, nbr_time_width

	size += sizeof(int); // int take_log
	size += sizeof(float); // float missing value

	size += sizeof(bool); //  bool doTrim;
	size += 2 * sizeof(float); // float trimMax, trimMin;

	size += sizeof(bool); //  bool doRemove;
	size += 2 * sizeof(float); // float removeMax, removeMin;

	size += 2 * sizeof(float); // float nbrsMax, nbrsMin;

	return size;
}

extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t RepNbrsOutlierCleaner::serialize(unsigned char *blob) {

	size_t ptr = 0;
	// SignalName
	size_t nameLen = signalName.length();
	assert(nameLen < MAX_NAME_LEN);

	strcpy(signalName_c, signalName.c_str());

	memcpy(blob + ptr, &nameLen, sizeof(size_t)); ptr += sizeof(size_t); 
	memcpy(blob + ptr, signalName_c, nameLen + 1); ptr += nameLen + 1;
	memcpy(blob + ptr, &time_channel, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &val_channel, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &nbr_time_unit, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &nbr_time_width, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &params.take_log, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &params.missing_value, sizeof(float)); ptr += sizeof(int);
	memcpy(blob + ptr, &params.doTrim, sizeof(bool)); ptr += sizeof(bool);
	memcpy(blob + ptr, &trimMax, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &trimMin, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &params.doRemove, sizeof(bool)); ptr += sizeof(bool);
	memcpy(blob + ptr, &removeMax, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &removeMin, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &nbrsMax, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &nbrsMin, sizeof(float)); ptr += sizeof(float);

	return ptr;
}

//.......................................................................................
size_t RepNbrsOutlierCleaner::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	// FeatureName
	size_t nameLen;
	memcpy(&nameLen, blob + ptr, sizeof(size_t)); ptr += sizeof(size_t);
	assert(nameLen < MAX_NAME_LEN);

	memcpy(signalName_c, blob + ptr, nameLen + 1); ptr += nameLen + 1;
	signalName = signalName_c;
	memcpy(&time_channel, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&val_channel, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&nbr_time_unit, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&nbr_time_width, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&params.take_log, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&params.missing_value, blob + ptr, sizeof(float)); ptr += sizeof(int);
	memcpy(&params.doTrim, blob + ptr, sizeof(bool)); ptr += sizeof(bool);
	memcpy(&trimMax, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&trimMin, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&params.doRemove, blob + ptr, sizeof(bool)); ptr += sizeof(bool);
	memcpy(&removeMax, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&removeMin, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&nbrsMax, blob + ptr, sizeof(float)); ptr += sizeof(float);
	memcpy(&nbrsMin, blob + ptr, sizeof(float)); ptr += sizeof(float);

	return ptr;
}

//.......................................................................................
void RepNbrsOutlierCleaner::print()
{
	MLOG("RepNbrsOutlierCleaner: signal: %d : doTrim %d trimMax %f trimMin %f : doRemove %d : removeMax %f : removeMin %f : nbrsMax %f : nbrsMin %f\n",
		signalId, params.doTrim, trimMax, trimMin, params.doRemove, removeMax, removeMin,nbrsMax, nbrsMin);
}

//=======================================================================================
// Utility Functions
//=======================================================================================
//.......................................................................................
// Get values of a signal from a set of ids
int get_values(MedRepository& rep, vector<int>& ids, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values, vector<RepProcessor *>& prev_processors) 
{
	vector<int> neededSignalIds = { signalId };
	PidDynamicRec rec;
	vector<int> req_signal_ids(1, signalId);

	UniversalSigVec usv;

	for (int id : ids) {

		// Get signal
		rep.uget(id, signalId, usv);

		if (prev_processors.size()) {

			// Get all time points
			vector<int> time_points(usv.len);
			for (int i = 0; i < usv.len; i++)
				time_points[i] = usv.Time(i, time_channel);

			// Init Dynamic Rec
			rec.init_from_rep(std::addressof(rep), id, req_signal_ids, (int)time_points.size());
			
			// Clean at all time-points
			for (auto& processor : prev_processors)
				processor->apply(rec, time_points, neededSignalIds);

			// Collect 
			for (int i = 0; i < usv.len; i++) {
				rec.uget(signalId, i, rec.usv);
				float ival = rec.usv.Val(i, val_channel);
				if (ival >= range_min && ival <= range_max)
					values.push_back(ival);
			}
		}
		else {
			// Collect 
			for (int i = 0; i < usv.len; i++) {
				float ival = usv.Val(i, val_channel);
				if (ival >= range_min && ival <= range_max)
					values.push_back(ival);
			}
		}
	}
	
	return 0;
}

//.......................................................................................
int get_values(MedRepository& rep, vector<int>& ids, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values) {
	vector<RepProcessor *> temp;
	return get_values(rep, ids, signalId, time_channel, val_channel, range_min, range_max, values, temp); 
}
