#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "RepProcess.h"

//=======================================================================================
// RepProcessors
//=======================================================================================
// Processors types from names
//.......................................................................................
RepProcessorTypes rep_processor_name_to_type(const string& processor_name) {

	if (processor_name == "multi_processor" || processor_name == "multi")
		return REP_PROCESS_MULTI;
	else if (processor_name == "basic_outlier_cleaner" || processor_name == "basic_cln")
		return REP_PROCESS_BASIC_OUTLIER_CLEANER;
	else if (processor_name == "nbrs_outlier_cleaner" || processor_name == "nbrs_cln")
		return REP_PROCESS_NBRS_OUTLIER_CLEANER;
	else if (processor_name == "configured_outlier_cleaner" || processor_name == "conf_cln")
		return REP_PROCESS_CONFIGURED_OUTLIER_CLEANER;
	else
		return REP_PROCESS_LAST;
}

// Create processor from params string (type must be given within string)
//.......................................................................................
RepProcessor *RepProcessor::create_processor(string &params)
{
	string rp_type;
	get_single_val_from_init_string(params, "rp_type", rp_type);
	return (make_processor(rp_type, params));
}

// Initialization : given processor name
//.......................................................................................
RepProcessor * RepProcessor::make_processor(string processor_name) {

	return make_processor(rep_processor_name_to_type(processor_name));
}

// Initialization : given processor name and intialization string
//.......................................................................................
RepProcessor * RepProcessor::make_processor(string processor_name, string init_string) {

	return make_processor(rep_processor_name_to_type(processor_name), init_string);
}

// Initialization : given processor type
//.......................................................................................
RepProcessor * RepProcessor::make_processor(RepProcessorTypes processor_type) {

	if (processor_type == REP_PROCESS_MULTI)
		return new RepMultiProcessor;
	else if (processor_type == REP_PROCESS_BASIC_OUTLIER_CLEANER)
		return new RepBasicOutlierCleaner;
	else if (processor_type == REP_PROCESS_NBRS_OUTLIER_CLEANER)
		return new RepNbrsOutlierCleaner;
	else if (processor_type == REP_PROCESS_CONFIGURED_OUTLIER_CLEANER)
		return new RepConfiguredOutlierCleaner;
	else
		return NULL;

}

// Initialization : given processor type and intialization string
//.......................................................................................
RepProcessor * RepProcessor::make_processor(RepProcessorTypes processor_type, string init_string) {

	//MLOG("Processor type is %d\n", (int)processor_type);
	RepProcessor *newRepProcessor = make_processor(processor_type);
	newRepProcessor->init_from_string(init_string);
	return newRepProcessor;
}

// Apply processing on a single PidDynamicRec at a set of time-points given by samples
//.......................................................................................
int RepProcessor::apply(PidDynamicRec& rec, MedIdSamples& samples) {

	vector<int> time_points(samples.samples.size());
	for (unsigned int i = 0; i < time_points.size(); i++)
		time_points[i] = samples.samples[i].time;

	return apply(rec, time_points);
}

// Apply processing on a single PidDynamicRec at a set of time-points given by time-points,
// only if affecting any of the signals given in neededSignalIds
//.......................................................................................
int RepProcessor::apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignalIds) {

	for (int signalId : neededSignalIds) {
		if (is_signal_affected(signalId))
			return apply(rec, time_points);
	}

	return 0;
}

// Apply processing on a single PidDynamicRec at a set of time-points given by samples,
// only if affecting any of the signals given in neededSignalIds
//.......................................................................................
int RepProcessor::apply(PidDynamicRec& rec, MedIdSamples& samples, vector<int>& neededSignalIds) {

	vector<int> time_points(samples.samples.size());
	for (unsigned int i = 0; i < time_points.size(); i++)
		time_points[i] = samples.samples[i].time;

	return apply(rec, time_points, neededSignalIds);
}

// Fill req_signal_ids from req_signals
//.......................................................................................
void RepProcessor::set_required_signal_ids(MedDictionarySections& dict) {

	req_signal_ids.resize(req_signals.size());

	for (unsigned int i = 0; i < req_signals.size(); i++)
		req_signal_ids[i] = dict.id(req_signals[i]);
}

// Append req_signal_ids to vector (fill req_signal_ids from req_signals if empty)
//.......................................................................................
void RepProcessor::get_required_signal_ids(unordered_set<int>& signalIds, MedDictionarySections& dict) {

	if (req_signal_ids.empty())
		set_required_signal_ids(dict);

	for (int signalId : req_signal_ids)
		signalIds.insert(signalId);
}

// Append req_signals to vector
//.......................................................................................
void RepProcessor::get_required_signal_names(unordered_set<string>& signalNames) {
	for (auto sig : req_signals)
		signalNames.insert(sig);
}

// Affected signals - set aff_signal_ids from aff_signals (id->name)
//.......................................................................................
void RepProcessor::set_affected_signal_ids(MedDictionarySections& dict) {

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
// Required Signals ids : Fill the member vector - req_signal_ids
//.......................................................................................
void RepMultiProcessor::set_required_signal_ids(MedDictionarySections& dict) {

	req_signal_ids.clear();
	for (auto& processor : processors) {
		req_signal_ids.clear();
		processor->set_required_signal_ids(dict);
		req_signal_ids.insert(req_signal_ids.end(), processor->req_signal_ids.begin(), processor->req_signal_ids.end());
	}
}

// Affected Signals : Fill the member vector aff_signal_ids
//.......................................................................................
void RepMultiProcessor::set_affected_signal_ids(MedDictionarySections& dict) {

	aff_signal_ids.clear();
	for (auto& processor : processors) {
		processor->aff_signal_ids.clear();
		processor->set_affected_signal_ids(dict);

		processor->aff_signal_ids.clear();
		for (int signalId : processor->aff_signal_ids)
			aff_signal_ids.insert(signalId);
	}
}

// Set signal-ids for all linked signals
//.......................................................................................
void RepMultiProcessor::set_signal_ids(MedDictionarySections& dict) {

	for (auto& processor : processors)
		processor->set_signal_ids(dict);

}

// Required Signals names : Fill the unordered set signalNames
//.......................................................................................
void RepMultiProcessor::get_required_signal_names(unordered_set<string>& signalNames) {
	for (auto& processor : processors)
		processor->get_required_signal_names(signalNames);
}

// Learn processors
//.......................................................................................
int RepMultiProcessor::Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_processors) {

	vector<int> rc(processors.size(), 0);

#pragma omp parallel for schedule(dynamic)
	for (int j=0; j<processors.size(); j++) {
		rc[j] = processors[j]->Learn(rep, ids, prev_processors);
	}

	for (int r : rc) if (r<0) return -1;
	return 0;
}

// Apply processors
//.......................................................................................
int RepMultiProcessor::apply(PidDynamicRec& rec, vector<int>& time_points) {

	vector<int> rc(processors.size(),0);

// ??? chances are this next parallelization is not needed, as we parallel before on recs...
#pragma omp parallel for schedule(dynamic)
	for (int j=0; j<processors.size(); j++) {
		rc[j] = processors[j]->apply(rec, time_points);
	}

	for (int r : rc) if (r<0) return -1;
	return 0;
}

// Apply processors that affect any of the needed signals
//.......................................................................................
int RepMultiProcessor::apply(PidDynamicRec& rec, vector<int>& time_points, vector<int>& neededSignalIds) {

	vector<int> rc(processors.size(), 0);

// ??? chances are this next parallelization is not needed, as we parallel before on recs...
#pragma omp parallel for schedule(dynamic)
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

// Add processors with initialization string
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

	for (auto& processor : processors) {
		ptr += processor->processor_serialize(blob + ptr);
	}

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
// Fill req- and aff-signals vectors
//.......................................................................................
void RepBasicOutlierCleaner::init_lists() {

	req_signals.push_back(signalName); 
	aff_signals.insert(signalName);
}

// Init from map
//.......................................................................................
 int RepBasicOutlierCleaner::init(map<string, string>& mapper) 
{ 
	init_defaults(); 

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "signal") { signalName = entry.second; req_signals.push_back(signalName); }
		else if (field == "time_channel") time_channel = stoi(entry.second);
		else if (field == "val_channel") val_channel = stoi(entry.second);
	}

	init_lists();
	return MedValueCleaner::init(mapper); 
}

 // Learn cleaning boundaries
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

// Learning : learn cleaning boundaries using MedValueCleaner's iterative approximation of moments
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

	// Iterative approximation of moments
	int rc =  get_iterative_min_max(values);

	return rc;
}

// Learning : learn cleaning boundaries using MedValueCleaner's quantile approximation of moments
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

	// Quantile approximation of moments
	return get_quantile_min_max(values);
}

// Apply cleaning model
//.......................................................................................
int  RepBasicOutlierCleaner::apply(PidDynamicRec& rec, vector<int>& time_points) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Check that we have the correct number of dynamic-versions : one per time-point
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

			// No need to clean past the latest relevant time-point
			if (itime > time_points[iver])	break;

			// Identify values to change or remove
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

		// Point versions jver+1..iver to jver
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

	size += MedSerialize::get_size(processor_type, signalName, time_channel, val_channel);
	size += MedSerialize::get_size(params.take_log, params.missing_value, params.doTrim, params.doRemove);
	size += MedSerialize::get_size(trimMax, trimMin, removeMax, removeMin);

	return size;
}

//.......................................................................................
size_t RepBasicOutlierCleaner::serialize(unsigned char *blob) {

	size_t ptr = 0;
	ptr += MedSerialize::serialize(blob + ptr, processor_type, signalName, time_channel, val_channel);
	ptr += MedSerialize::serialize(blob + ptr, params.take_log, params.missing_value, params.doTrim, params.doRemove);
	ptr += MedSerialize::serialize(blob + ptr, trimMax, trimMin, removeMax, removeMin);


	return ptr;
}

//.......................................................................................
size_t RepBasicOutlierCleaner::deserialize(unsigned char *blob) {

	size_t ptr = 0;
	ptr += MedSerialize::deserialize(blob + ptr, processor_type, signalName, time_channel, val_channel);
	ptr += MedSerialize::deserialize(blob + ptr, params.take_log, params.missing_value, params.doTrim, params.doRemove);
	ptr += MedSerialize::deserialize(blob + ptr, trimMax, trimMin, removeMax, removeMin);


	return ptr;
}

//.......................................................................................
void RepBasicOutlierCleaner::print()
{
	MLOG("BasicOutlierCleaner: signal: %d : doTrim %d trimMax %f trimMin %f : doRemove %d : removeMax %f removeMin %f\n",
		signalId, params.doTrim, trimMax, trimMin, params.doRemove, removeMax, removeMin);
}


//=======================================================================================
// ConfiguredOutlierCleaner  
//=======================================================================================
//.......................................................................................
int readConfFile(string confFileName, map<string, confRecord>& outlierParams)
// read from outlierParamFile into outlierParams map
{
	ifstream infile;
	confRecord thisRecord;
	string thisLine;
	infile.open(confFileName.c_str(),ifstream::in);
	if(!infile.is_open()){
		fprintf(stderr, "Cannot open %s for reading\n", confFileName.c_str());
		return -1;
	}
	getline(infile,thisLine);//consume title line.
	while (!infile.eof()) {
		getline(infile,thisLine);
		vector<string> f;
		boost::split(f, thisLine, boost::is_any_of(","));
		if(f.size()!=7){
			fprintf(stderr, "Wrong field count in  %s \n", confFileName.c_str());
			infile.close();
			return -1;
		}
		thisRecord.confirmedLow=thisRecord.logicalLow =atof( f[1].c_str());
		thisRecord.confirmedHigh=thisRecord.logicalHigh = atof(f[2].c_str());
		
			
		thisRecord.distLow = f[4];
		thisRecord.distHigh = f[6];
		if (thisRecord.distLow != "none")thisRecord.confirmedLow = atof(f[3].c_str());
		if (thisRecord.distHigh != "none")thisRecord.confirmedHigh = atof(f[5].c_str());
		outlierParams[f[0]] = thisRecord;
	}
	infile.close();
	return(0);

}
int RepConfiguredOutlierCleaner::init(map<string, string>& mapper)
{
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "signal") { signalName = entry.second; req_signals.push_back(signalName); }
		else if (field == "time_channel") time_channel = stoi(entry.second);
		else if (field == "val_channel") val_channel = stoi(entry.second);
		else if (field == "conf_file") {
			confFileName = entry.second; if (int res = readConfFile(confFileName, outlierParams))return(res);
		}
		else if (field == "clean_method")cleanMethod = entry.second;
	}

	init_lists();
	return MedValueCleaner::init(mapper);
}

// Learn bounds
//.......................................................................................
int RepConfiguredOutlierCleaner::Learn(MedPidRepository& rep, vector<int>& ids, vector<RepProcessor *>& prev_cleaners) {
	if (outlierParams.find(signalName) == outlierParams.end()) {
		MERR("MedModel learn() : ERROR: Signal not supported by conf_cln()\n");
		return -1;
	}
	trimMax = 1e+98;
	trimMin = -1e+98;

	if (cleanMethod == "logical") {
		removeMax = outlierParams[signalName].logicalHigh;
		removeMin = outlierParams[signalName].logicalLow;
		return(0);
	}
	else if (cleanMethod == "confirmed") {
		removeMax = outlierParams[signalName].confirmedHigh;
		removeMin = outlierParams[signalName].confirmedLow;
		return(0);
	}
	else if (cleanMethod == "learned") {
		removeMax = outlierParams[signalName].logicalHigh;
		removeMin = outlierParams[signalName].logicalLow;
		string thisDistHi = outlierParams[signalName].distHigh;
		string thisDistLo = outlierParams[signalName].distLow;
		if (thisDistHi == "none" && thisDistLo == "none") return(0);//nothing to learn

		else {
			
			vector<float> values, filteredValues;

			double borderHi, borderLo, logBorderHi, logBorderLo;
			get_values(rep, ids, signalId, time_channel, val_channel, removeMin, removeMax, values, prev_cleaners);
			for (auto& el : values)if (el != 0)filteredValues.push_back(el);
			sort(filteredValues.begin(), filteredValues.end());
			if (thisDistHi == "norm" || thisDistLo == "norm")
				learnDistributionBorders(borderHi, borderLo, filteredValues);
			if (thisDistHi == "lognorm" || thisDistLo == "lognorm") {
			/*	ofstream dFile;
				dFile.open("DFILE");
				for (auto& el : filteredValues)dFile << el << "\n";
				dFile.close();
				*/


				for (auto& el : filteredValues)
					if (el > 0)el = log(el);
					else return(-1);
				
				learnDistributionBorders(logBorderHi, logBorderLo, filteredValues);
			}
			if (thisDistHi == "norm")removeMax = borderHi;
			else if (thisDistHi == "lognorm")removeMax = exp(logBorderHi);
			else if (thisDistHi == "manual")removeMax = outlierParams[signalName].confirmedHigh;
			if (thisDistLo == "norm")removeMin = borderLo;
			else if (thisDistLo == "lognorm")removeMin = exp(logBorderLo);
			else if (thisDistLo == "manual")removeMin = outlierParams[signalName].confirmedLow;
			
			return(0);
		}

	}

	else {
		MERR("Unknown cleaning method %s\n", cleanMethod.c_str());
		return -1;
	}
}


void learnDistributionBorders(double& borderHi, double& borderLo, vector<float> filteredValues)
// a function that takes sorted vector of filtered values and estimates the +- 7 sd borders based on the center of distribution
// predefined calibration constants are used for estimation of the borders. 
{
	double sum = 0;
	double sumsq = 0;
	const float margin[] = {(float) 0.01, (float)0.99 };// avoid tails of distribution
	const float varianceFactor = 0.8585;
	const float meanShift = 0; // has value when margins are asymetric
	const float sdNums = 7; // how many standard deviation on each side of the mean.

	int start = round(filteredValues.size()*margin[0]);
	int stop = round(filteredValues.size()*margin[1]);
	for (vector<float>::iterator el = filteredValues.begin()+start; el < filteredValues.begin()+stop; el++) {
		
		sum += *el;
		sumsq += *el* *el;
	}
	 double mean= sum/(stop-start);
	 double var =sumsq/(stop-start) - mean*mean;
	 printf("sum %f sumsq %f  stop %d start %d\n", sum, sumsq, stop, start);
	 var = var / varianceFactor;
	 mean=  mean - meanShift*sqrt(var);
	 borderHi = mean + sdNums*sqrt(var);
	 borderLo = mean - sdNums*sqrt(var);



}


//.......................................................................................
size_t RepConfiguredOutlierCleaner::get_size() {

	size_t size = 0;

	size += MedSerialize::get_size(processor_type, signalName, time_channel, val_channel);
	size += MedSerialize::get_size(params.take_log, params.missing_value, params.doTrim, params.doRemove);
	size += MedSerialize::get_size(trimMax, trimMin, removeMax, removeMin);
	size += MedSerialize::get_size(confFileName, cleanMethod, outlierParams);

	return size;
}

//.......................................................................................
size_t RepConfiguredOutlierCleaner::serialize(unsigned char *blob) {

	size_t ptr = 0;
	ptr += MedSerialize::serialize(blob + ptr, processor_type, signalName, time_channel, val_channel);
	ptr += MedSerialize::serialize(blob + ptr, params.take_log, params.missing_value, params.doTrim, params.doRemove);
	ptr += MedSerialize::serialize(blob + ptr, trimMax, trimMin, removeMax, removeMin);
	ptr += MedSerialize::serialize(blob + ptr, confFileName, cleanMethod, outlierParams);


	return ptr;
}

//.......................................................................................
size_t RepConfiguredOutlierCleaner::deserialize(unsigned char *blob) {

	size_t ptr = 0;
	ptr += MedSerialize::deserialize(blob + ptr, processor_type, signalName, time_channel, val_channel);
	ptr += MedSerialize::deserialize(blob + ptr, params.take_log, params.missing_value, params.doTrim, params.doRemove);
	ptr += MedSerialize::deserialize(blob + ptr, trimMax, trimMin, removeMax, removeMin);
	ptr += MedSerialize::deserialize(blob + ptr, confFileName, cleanMethod, outlierParams);


	return ptr;
}

//.......................................................................................
void RepConfiguredOutlierCleaner::print()
{
	MLOG("BasicOutlierCleaner: signal: %d : doTrim %d trimMax %f trimMin %f : doRemove %d : removeMax %f removeMin %f\n",
		signalId, params.doTrim, trimMax, trimMin, params.doRemove, removeMax, removeMin, confFileName.c_str(), cleanMethod.c_str());
}

//=======================================================================================
// NbrsOutlierCleaner
//=======================================================================================
void RepNbrsOutlierCleaner::init_lists() {

	req_signals.push_back(signalName);
	aff_signals.insert(signalName);
}

//.......................................................................................
int RepNbrsOutlierCleaner::init(map<string, string>& mapper)
{
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "signal") { signalName = entry.second; req_signals.push_back(signalName); }
		else if (field == "time_channel") time_channel = stoi(entry.second);
		else if (field == "val_channel") val_channel = stoi(entry.second);
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

	int rc = get_iterative_min_max(values);
	return rc;
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
		int verLen = len;
		vector<int> candidates(len, 0);
		vector<int> removed(len, 0);
		
		for (int i = 0; i < len; i++) {
			int itime = rec.usv.Time(i, time_channel);

			if (itime > time_points[iver]) {
				verLen = i ;
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

	size += MedSerialize::get_size(processor_type, signalName, time_channel, val_channel, nbr_time_unit, nbr_time_width);
	size += MedSerialize::get_size(params.take_log, params.missing_value, params.doTrim, params.doRemove);
	size += MedSerialize::get_size(trimMax, trimMin, removeMax, removeMin, nbrsMax, nbrsMin);

	return size;
}

extern char signalName_c[MAX_NAME_LEN + 1];

//.......................................................................................
size_t RepNbrsOutlierCleaner::serialize(unsigned char *blob) {

	size_t ptr = 0;

	ptr += MedSerialize::serialize(blob + ptr, processor_type, signalName, time_channel, val_channel, nbr_time_unit, nbr_time_width);
	ptr += MedSerialize::serialize(blob + ptr, params.take_log, params.missing_value, params.doTrim, params.doRemove);
	ptr += MedSerialize::serialize(blob + ptr, trimMax, trimMin, removeMax, removeMin, nbrsMax, nbrsMin);

	return ptr;
}

//.......................................................................................
size_t RepNbrsOutlierCleaner::deserialize(unsigned char *blob) {

	size_t ptr = 0;

	ptr += MedSerialize::deserialize(blob + ptr, processor_type, signalName, time_channel, val_channel, nbr_time_unit, nbr_time_width);
	ptr += MedSerialize::deserialize(blob + ptr, params.take_log, params.missing_value, params.doTrim, params.doRemove);
	ptr += MedSerialize::deserialize(blob + ptr, trimMax, trimMin, removeMax, removeMin, nbrsMax, nbrsMin);

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
		
		// Nothing to do if empty ...
		if (usv.len == 0)
			continue;

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
