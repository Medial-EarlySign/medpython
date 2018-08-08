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
	else if (processor_name == "configured_outlier_cleaner" || processor_name == "conf_cln")
		return REP_PROCESS_CONFIGURED_OUTLIER_CLEANER;
	else if (processor_name == "rulebased_outlier_cleaner" || processor_name == "rule_cln")
		return REP_PROCESS_RULEBASED_OUTLIER_CLEANER;
	else if (processor_name == "calc_signals" || processor_name == "calculator")
		return REP_PROCESS_CALC_SIGNALS;
	else if (processor_name == "complete")
		return REP_PROCESS_COMPLETE;
	else if (processor_name == "req" || processor_name == "requirements")
		return REP_PROCESS_CHECK_REQ;
	else if (processor_name == "sim_val" || processor_name == "sim_val_handler")
		return REP_PROCESS_SIM_VAL;
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
	else if (processor_type == REP_PROCESS_RULEBASED_OUTLIER_CLEANER)
		return new RepRuleBasedOutlierCleaner;
	else if (processor_type == REP_PROCESS_CALC_SIGNALS)
		return new RepCalcSimpleSignals;
	else if (processor_type == REP_PROCESS_COMPLETE)
		return new RepPanelCompleter;
	else if (processor_type == REP_PROCESS_CHECK_REQ)
		return new RepCheckReq;
	else if (processor_type == REP_PROCESS_SIM_VAL)
		return new RepSimValHandler;
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

// learn on all pids in repository, using fake samples - works only for repProcessors that ignore sample dates
//.......................................................................................
int RepProcessor::learn(MedPidRepository& rep) {
	MedSamples fakeSamples;
	for (int pid : rep.pids)
		fakeSamples.insertRec(pid, 0);
	this->learn(rep, fakeSamples);
	return 0;
}

// Learn processing parameters only if affecting any of the signals given in neededSignalIds
//.......................................................................................
int RepProcessor::_conditional_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors, unordered_set<int>& neededSignalIds) {
	for (int signalId : neededSignalIds) {
		if (is_signal_affected(signalId)) 
			return _learn(rep, samples, prev_processors);
	}
	return 0;
}

// Check if processor can be filtered
//...................................
bool RepProcessor::filter(unordered_set<string>& neededSignals) {

	if (unconditional)
		return false;

	for (string signal : neededSignals) {
		if (is_signal_affected(signal))
			return false;
	}

	MLOG("RepProcessor::filter filtering out processor of type %d, affected signals: ", processor_type);
	for (string signal : aff_signals)
		MLOG("[%s] ", signal.c_str());
	MLOG("\n");
	return true;

}

// Apply processing on a single PidDynamicRec at a set of time-points given by samples
//.......................................................................................
int RepProcessor::apply(PidDynamicRec& rec, MedIdSamples& samples) {

	vector<int> time_points(samples.samples.size());
	for (unsigned int i = 0; i < time_points.size(); i++)
		time_points[i] = samples.samples[i].time;

	vector<vector<float>> attributes_mat(time_points.size(), vector<float>(attributes.size(),0)) ;
	int rc = apply(rec, time_points,attributes_mat);

	if (rc == 0) {
		for (unsigned int i = 0; i < time_points.size(); i++) {
			for (int j = 0; j < attributes.size(); j++)
				samples.samples[i].attributes[attributes[j]] += attributes_mat[i][j];
		}
	}

	return rc;
}

// Apply processing on a single PidDynamicRec at a set of time-points given by samples,
// only if affecting any of the signals given in neededSignalIds
//.......................................................................................
int RepProcessor::conditional_apply(PidDynamicRec& rec, MedIdSamples& samples, unordered_set<int>& neededSignalIds) {

	vector<int> time_points(samples.samples.size());
	for (unsigned int i = 0; i < time_points.size(); i++)
		time_points[i] = samples.samples[i].time;

	vector<vector<float>> attributes_mat(time_points.size(), vector<float>(attributes.size(), 0));
	int rc = conditional_apply(rec, time_points, neededSignalIds, attributes_mat);

	if (rc == 0) {
		for (unsigned int i = 0; i < time_points.size(); i++) {
			for (int j = 0; j < attributes.size(); j++) {
				samples.samples[i].attributes[attributes[j]] += attributes_mat[i][j];
			}
		}
	}

	return rc;
}

// Apply processing on a single PidDynamicRec at a set of time-points given by time-points,
// only if affecting any of the signals given in neededSignalIds
//.......................................................................................
int RepProcessor::_conditional_apply(PidDynamicRec& rec, vector<int>& time_points, unordered_set<int>& neededSignalIds, vector<vector<float>>& attributes_mat) {

	if (unconditional)
		return apply(rec, time_points, attributes_mat);

	for (int signalId : neededSignalIds) {
		if (is_signal_affected(signalId)) 
			return apply(rec, time_points, attributes_mat);
	}

	return 0;
}

// Fill req_signal_ids from req_signals
//.......................................................................................
void RepProcessor::set_required_signal_ids(MedDictionarySections& dict) {

	for (string signal : req_signals)
		req_signal_ids.insert(dict.id(signal));

}

// Add req_signals to set
//.......................................................................................
void RepProcessor::get_required_signal_names(unordered_set<string>& signalNames) {

	for (auto sig : req_signals)
		signalNames.insert(sig);
}

// Add req_signals to set only if processor is required for any of preReqSignalNames
//.......................................................................................
void RepProcessor::get_required_signal_names(unordered_set<string>& signalNames, unordered_set<string> preReqSignalNames) {
	
	if (unconditional)
		get_required_signal_names(signalNames);
	else {
		for (string signal : preReqSignalNames) {
			if (is_signal_affected(signal)) {
				get_required_signal_names(signalNames);
				return;
			}
		}
	}

}

//.......................................................................................
void RepProcessor::get_required_signal_ids(unordered_set<int>& signalIds) {

	for (auto sig : req_signal_ids)
		signalIds.insert(sig);
}

// Add req_signals to set only if processor is required for any of preReqSignalNames
//.......................................................................................
void RepProcessor::get_required_signal_ids(unordered_set<int>& signalIds, unordered_set<int> preReqSignals) {

	for (int signal : preReqSignals) {
		if (is_signal_affected(signal)) {
			get_required_signal_ids(signalIds);
			return;
		}
	}
}

// Affected signals - set aff_signal_ids from aff_signals (id->name)
//.......................................................................................
void RepProcessor::set_affected_signal_ids(MedDictionarySections& dict) {

	for (string signalName : aff_signals)
		aff_signal_ids.insert(dict.id(signalName));
}

//.......................................................................................
void RepProcessor::dprint(const string &pref, int rp_flag)
{
	if (rp_flag > 0) {
		MLOG("%s :: RP type %d : required(%d): ", pref.c_str(), processor_type, req_signals.size());
		if (rp_flag > 1) for (auto &rsig : req_signals) MLOG("%s,", rsig.c_str());
		MLOG(" affected(%d): ", aff_signals.size());
		if (rp_flag > 1) for (auto &asig : aff_signals) MLOG("%s, ", asig.c_str());
		MLOG(" virtual(%d): ", virtual_signals.size());
		if (rp_flag > 1) for (auto &vsig : virtual_signals) MLOG("%s ", vsig.first.c_str());
		MLOG("\n");
	}
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
//.......................................................................................
void RepMultiProcessor::clear()
{
	for (auto p : processors) {
		if (p != NULL) {
			delete p;
			p = NULL;
		}
	}
	processors.clear();
}

// Required Signals ids : Fill the member vector - req_signal_ids
//.......................................................................................
void RepMultiProcessor::set_required_signal_ids(MedDictionarySections& dict) {

	req_signal_ids.clear();
	for (auto& processor : processors) {
		req_signal_ids.clear();
		processor->set_required_signal_ids(dict);

		for (int signalId : processor->req_signal_ids)
			req_signal_ids.insert(signalId);
	}
}

// Affected Signals : Fill the member vector aff_signal_ids
//.......................................................................................
void RepMultiProcessor::set_affected_signal_ids(MedDictionarySections& dict) {

	aff_signal_ids.clear();
	for (auto& processor : processors) {
		processor->aff_signal_ids.clear();
		processor->set_affected_signal_ids(dict);

		for (int signalId : processor->aff_signal_ids)
			aff_signal_ids.insert(signalId);
	}
}

// Check if processor can be filtered
//...................................
bool RepMultiProcessor::filter(unordered_set<string>& neededSignals) {

	vector<RepProcessor *> filtered;
	for (auto& processor : processors) {
		if (!processor->filter(neededSignals))
			filtered.push_back(processor);
	}

	if (filtered.empty()) {
		MLOG("RepMultiProcessor::filter filtering out processor of type %d\n", processor_type);
		return true;
	}
	else {
		processors = filtered;
		return false;
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

// Add req_signals to set only if processor is required for any of preReqSignalNames
// Note that preReq is copied so it is not affected by enlarging signalNames
//.......................................................................................
void RepMultiProcessor::get_required_signal_names(unordered_set<string>& signalNames, unordered_set<string> preReqSignalNames) {
	for (auto& processor : processors)
		processor->get_required_signal_names(signalNames, preReqSignalNames);
}

// Required Signals names : Fill the unordered set signalIds
//.......................................................................................
void RepMultiProcessor::get_required_signal_ids(unordered_set<int>& signalIds) {
	for (auto& processor : processors)
		processor->get_required_signal_ids(signalIds);
}

// Virtual Signals names : Get the virtual signals map
//.......................................................................................
void RepMultiProcessor::add_virtual_signals(map<string, int> &_virtual_signals)
{
	for (auto& processor : processors)
		processor->add_virtual_signals(_virtual_signals);
}

// Add req_signals to set only if processor is required for any of preReqSignalNames
// Note that preReq is copied so it is not affected by enlarging signalNames
//.......................................................................................
void RepMultiProcessor::get_required_signal_ids(unordered_set<int>& signalIds, unordered_set<int> preReqSignals) {
	for (auto& processor : processors)
		processor->get_required_signal_ids(signalIds, preReqSignals);
}

// Learn processors
//.......................................................................................
int RepMultiProcessor::_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors) {

	vector<int> rc(processors.size(), 0);

#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < processors.size(); j++) {
		rc[j] = processors[j]->learn(rep, samples, prev_processors);
	}

	for (int r : rc) if (r < 0) return -1;
	return 0;
}

// Learn processing parameters only if affecting any of the signals given in neededSignalIds
//.......................................................................................
int RepMultiProcessor::_conditional_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_processors, unordered_set<int>& neededSignalIds) {

	vector<int> rc(processors.size(), 0);

#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < processors.size(); j++) {
		rc[j] = processors[j]->conditional_learn(rep, samples, prev_processors,neededSignalIds);
	}

	for (int r : rc) if (r < 0) return -1;
	return 0;
}

// Apply processors
//.......................................................................................
int RepMultiProcessor::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {

	vector<int> rc(processors.size(), 0);

	// If attributes are required, prepare space for collecting them
	vector<vector<vector<float> > > all_attributes_mats(processors.size());
	if (!attributes_mat.empty()) {
		all_attributes_mats.resize(processors.size());
		for (int j = 0; j < processors.size(); j++) {
			all_attributes_mats[j].resize(time_points.size());
			for (int i = 0; i < time_points.size(); i++)
				all_attributes_mats[j][i].resize(processors[j]->attributes.size(), 0);
		}
	}

	// ??? chances are this next parallelization is not needed, as we parallel before on recs...
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < processors.size(); j++) {
		rc[j] = processors[j]->apply(rec, time_points, all_attributes_mats[j]);
	}

	for (int r : rc) if (r < 0) return -1;

	// If attributes are required, collect them 
	if (!attributes_mat.empty()) {
		for (int j = 0; j < processors.size(); j++) {
			for (int i = 0; i < time_points.size(); i++) {
				for (int k = 0; k < processors[j]->attributes.size(); k++)
					attributes_mat[i][attributes_map[j][k]] += all_attributes_mats[j][i][k];
			}
		}
	}

	return 0;
}

// Apply processors that affect any of the needed signals
//.......................................................................................
int RepMultiProcessor::_conditional_apply(PidDynamicRec& rec, vector<int>& time_points, unordered_set<int>& neededSignalIds, vector<vector<float>>& attributes_mat) {

	vector<int> rc(processors.size(), 0);

	// If attributes are required, prepare space for collecting them
	vector<vector<vector<float> > > all_attributes_mats(processors.size());
	if (!attributes_mat.empty()) {
		all_attributes_mats.resize(processors.size());
		for (int j = 0; j < processors.size(); j++) {
			all_attributes_mats[j].resize(time_points.size());
			for (int i = 0; i < time_points.size(); i++) 
				all_attributes_mats[j][i].resize(processors[j]->attributes.size(), 0);
		}
	}

	// ??? chances are this next parallelization is not needed, as we parallel before on recs...
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < processors.size(); j++) {
		if (unconditional)
			rc[j] = processors[j]->apply(rec, time_points, all_attributes_mats[j]);
		else
			rc[j] = processors[j]->conditional_apply(rec, time_points, neededSignalIds, all_attributes_mats[j]);
	}

	for (int r : rc) if (r < 0) return -1;

	// If attributes are required, collect them 
	if (!attributes_mat.empty()) {
		for (int j = 0; j < processors.size(); j++) {
			for (int i = 0; i < time_points.size(); i++) {
				for (int k = 0; k < processors[j]->attributes.size(); k++) 
					attributes_mat[i][attributes_map[j][k]] += all_attributes_mats[j][i][k];
			}
		}
	}

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

// init attributes list and attributes map
//.......................................................................................
void RepMultiProcessor::init_attributes() {

	attributes.clear();
	map<string, int> attributes_pos;

	attributes_map.resize(processors.size());
	for (int i = 0; i < processors.size(); i++) {
		processors[i]->init_attributes();
		attributes_map[i].resize(processors[i]->attributes.size());
		for (int j = 0; j < processors[i]->attributes.size(); j++) {
			if (attributes_pos.find(processors[i]->attributes[j]) == attributes_pos.end()) {
				attributes.push_back(processors[i]->attributes[j]);
				attributes_pos[attributes.back()] = (int) attributes.size() - 1;
			}
			attributes_map[i][j] = attributes_pos[processors[i]->attributes[j]];
		}
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

//.......................................................................................
void RepMultiProcessor::dprint(const string &pref, int rp_flag)
{
	if (rp_flag > 0) {
		MLOG("%s :: RP MULTI(%d) -->\n", pref.c_str(), processors.size());
		for (auto& proc : processors) {
			proc->dprint(pref+"->Multi", rp_flag);
		}
	}
}

//=======================================================================================
// BasicOutlierCleaner
//=======================================================================================
// Fill req- and aff-signals vectors
//.......................................................................................
void RepBasicOutlierCleaner::init_lists() {

	req_signals.insert(signalName);
	aff_signals.insert(signalName);
}

// Init from map
//.......................................................................................
int RepBasicOutlierCleaner::init(map<string, string>& mapper) 
{ 
	init_defaults(); 

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepBasicOutlierCleaner::init]
		if (field == "signal") { signalName = entry.second; }
		else if (field == "time_channel") time_channel = med_stoi(entry.second);
		else if (field == "val_channel") val_channel = med_stoi(entry.second);
		else if (field == "nrem_attr") nRem_attr = entry.second;
		else if (field == "ntrim_attr") nTrim_attr = entry.second;
		else if (field == "nrem_suff") nRem_attr_suffix = entry.second;
		else if (field == "ntrim_suff") nTrim_attr_suffix = entry.second;
		//! [RepBasicOutlierCleaner::init]
	}

	init_lists();
	return MedValueCleaner::init(mapper);
}

// init attributes list
//.......................................................................................
void RepBasicOutlierCleaner::init_attributes() {

	string _signal_name = signalName;
	if (val_channel != 0) _signal_name += "_" + to_string(val_channel);

	attributes.clear();
	if (!nRem_attr.empty()) attributes.push_back(nRem_attr);
	if (!nRem_attr_suffix.empty()) attributes.push_back(_signal_name + "_" + nRem_attr_suffix);

	if (!nTrim_attr.empty()) attributes.push_back(nTrim_attr);
	if (!nTrim_attr_suffix.empty()) attributes.push_back(_signal_name + "_" + nTrim_attr_suffix);
}

// Learn cleaning boundaries
//.......................................................................................
int RepBasicOutlierCleaner::_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_cleaners) {

	if (params.type == VAL_CLNR_ITERATIVE)
		return iterativeLearn(rep, samples, prev_cleaners);
	else if (params.type == VAL_CLNR_QUANTILE)
		return quantileLearn(rep, samples, prev_cleaners);
	else {
		MERR("Unknown cleaning type %d\n", params.type);
		return -1;
	}
}

// Learning : learn cleaning boundaries using MedValueCleaner's iterative approximation of moments
//.......................................................................................
int RepBasicOutlierCleaner::iterativeLearn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_cleaners) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Get all values
	vector<float> values;
	get_values(rep, samples, signalId, time_channel, val_channel, params.range_min, params.range_max, values, prev_cleaners);
	//MLOG("basic Iterative clean Learn: signalName %s signalId %d :: got %d values()\n", signalName.c_str(), signalId, values.size());

	// Iterative approximation of moments
	int rc =  get_iterative_min_max(values);

	return rc;
}

// Learning : learn cleaning boundaries using MedValueCleaner's quantile approximation of moments
//.......................................................................................
int RepBasicOutlierCleaner::quantileLearn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_cleaners) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Get all values
	vector<float> values;
	get_values(rep, samples, signalId, time_channel, val_channel, params.range_min, params.range_max, values, prev_cleaners);

	if (values.empty()) {
		MWARN("RepBasicOutlierCleaner::quantileLearn WARNING signal [%d] = [%s] is empty, will not clean outliers\n", signalId,
			this->signalName.c_str());
		return 0;
	}
	// Quantile approximation of moments
	return get_quantile_min_max(values);
}

// Apply cleaning model
//.......................................................................................
int  RepBasicOutlierCleaner::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float> >& attributes_mat) {

	//MLOG("basic cleaner _apply: signalName %s signalId %d\n", signalName.c_str(), signalId);

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Check that we have the correct number of dynamic-versions : one per time-point (if given)
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int len;

	differentVersionsIterator vit(rec, signalId);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		// Clean 
		rec.uget(signalId, iver, rec.usv); // get into the internal usv obeject - this statistically saves init time

		len = rec.usv.len;
		vector<int> remove(len);
		vector<pair<int, float>> change(len);
		int nRemove = 0, nChange = 0;

		// Collect
		for (int i = 0; i < len; i++) {
			int itime = rec.usv.Time(i, time_channel);
			float ival = rec.usv.Val(i, val_channel);

			// No need to clean past the latest relevant time-point
			if (time_points.size() != 0 && itime > time_points[iver])	break;

			// Identify values to change or remove
			if (params.doRemove && (ival < removeMin - NUMERICAL_CORRECTION_EPS || ival > removeMax + NUMERICAL_CORRECTION_EPS)) {
//				MLOG("pid %d ver %d time %d %s %f removed\n", rec.pid, iver, itime, signalName.c_str(), ival);
				remove[nRemove++] = i;
			}
			else if (params.doTrim) {
				if (ival < trimMin) {
//					MLOG("pid %d ver %d time %d %s %f trimmed\n", rec.pid, iver, itime, signalName.c_str(), ival);
					change[nChange++] = pair<int, float>(i, trimMin);
				}
				else if (ival > trimMax) {
//					MLOG("pid %d ver %d time %d %s %f trimmed\n", rec.pid, iver, itime, signalName.c_str(), ival);
					change[nChange++] = pair<int, float>(i, trimMax);
				}
			}
		}

		// Apply removals + changes
		change.resize(nChange);
		remove.resize(nRemove);
		if (rec.update(signalId, iver, val_channel, change, remove) < 0)
			return -1;

		// Collect atttibutes
		int idx = 0;
		if (!nRem_attr.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++) 
				attributes_mat[pVersion][idx] = (float) nRemove;
			idx++;
		}

		if (!nRem_attr_suffix.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float)nRemove;
			idx++;
		}

		if (!nTrim_attr.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++) 
				attributes_mat[pVersion][idx] = (float) nChange;
			idx++;
		}

		if (!nTrim_attr_suffix.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float) nChange;
			idx++;
		}
	}

	return 0;

}

//.......................................................................................
void RepBasicOutlierCleaner::print()
{
	MLOG("BasicOutlierCleaner: signal: %d %s : v_channel %d : doTrim %d trimMax %f trimMin %f : doRemove %d : removeMax %f removeMin %f\n",
		signalId, signalName.c_str(), val_channel, params.doTrim, trimMax, trimMin, params.doRemove, removeMax, removeMin);
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
	infile.open(confFileName.c_str(), ifstream::in);
	if (!infile.is_open()) {
		fprintf(stderr, "Cannot open %s for reading\n", confFileName.c_str());
		return -1;
	}
	getline(infile, thisLine);//consume title line.
	while (!infile.eof()) {
		getline(infile, thisLine);
		if (thisLine.back() == '\r') thisLine.pop_back();
	
		vector<string> f;
		boost::split(f, thisLine, boost::is_any_of(","));
		if (f.size() != 7) {
			fprintf(stderr, "Wrong field count in  %s (%s : %zd) \n", confFileName.c_str(),thisLine.c_str(),f.size());
			infile.close();
			return -1;
		}
	
		thisRecord.confirmedLow = thisRecord.logicalLow = (float)atof(f[1].c_str());
		thisRecord.confirmedHigh = thisRecord.logicalHigh = (float)atof(f[2].c_str());

		thisRecord.distLow = f[4];
		thisRecord.distHigh = f[6];
		if (thisRecord.distLow != "none")thisRecord.confirmedLow = (float)atof(f[3].c_str());
		if (thisRecord.distHigh != "none")thisRecord.confirmedHigh = (float)atof(f[5].c_str());
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
		//! [RepConfiguredOutlierCleaner::init]
		if (field == "signal") { signalName = entry.second; req_signals.insert(signalName); }
		else if (field == "time_channel") time_channel = med_stoi(entry.second);
		else if (field == "val_channel") val_channel = med_stoi(entry.second);
		else if (field == "nrem_attr") nRem_attr = entry.second;
		else if (field == "ntrim_attr") nTrim_attr == entry.second;
		else if (field == "nrem_suff") nRem_attr_suffix = entry.second;
		else if (field == "ntrim_suff") nTrim_attr_suffix = entry.second;
		else if (field == "conf_file") {
			confFileName = entry.second; if (int res = readConfFile(confFileName, outlierParams))return(res);
		}
		else if (field == "clean_method")cleanMethod = entry.second;
		//! [RepConfiguredOutlierCleaner::init]
	}

	init_lists();
	return MedValueCleaner::init(mapper);
}

// Learn bounds
//.......................................................................................
int RepConfiguredOutlierCleaner::_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_cleaners) {
	if (outlierParams.find(signalName) == outlierParams.end()) {
		MERR("MedModel learn() : ERROR: Signal %s not supported by conf_cln()\n", signalName.c_str());
		return -1;
	}

	trimMax = 1e30F;
	trimMin = -1e+30F;

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

			float borderHi, borderLo, logBorderHi, logBorderLo;
			get_values(rep, samples, signalId, time_channel, val_channel, removeMin, removeMax, values, prev_cleaners);
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
			else if (thisDistHi == "lognorm")removeMax = expf(logBorderHi);
			else if (thisDistHi == "manual")removeMax = outlierParams[signalName].confirmedHigh;
			if (thisDistLo == "norm")removeMin = borderLo;
			else if (thisDistLo == "lognorm")removeMin = expf(logBorderLo);
			else if (thisDistLo == "manual")removeMin = outlierParams[signalName].confirmedLow;

			return(0);
		}

	}

	else {
		MERR("Unknown cleaning method %s\n", cleanMethod.c_str());
		return -1;
	}
}

void learnDistributionBorders(float& borderHi, float& borderLo, vector<float> filteredValues)
// a function that takes sorted vector of filtered values and estimates the +- 7 sd borders based on the center of distribution
// predefined calibration constants are used for estimation of the borders. 
{
	double sum = 0;
	double sumsq = 0;
	const float margin[] = { 0.01F, 0.99F };// avoid tails of distribution
	const float varianceFactor = 0.8585F;
	const float meanShift = 0; // has value when margins are asymetric
	const float sdNums = 7; // how many standard deviation on each side of the mean.

	int start = (int) round(filteredValues.size()*margin[0]);
	int stop = (int) round(filteredValues.size()*margin[1]);
	for (vector<float>::iterator el = filteredValues.begin() + start; el < filteredValues.begin() + stop; el++) {

		sum += *el;
		sumsq += *el* *el;
	}
	double mean = sum / (stop - start);
	double var = sumsq / (stop - start) - mean*mean;
	//printf("sum %f sumsq %f  stop %d start %d\n", sum, sumsq, stop, start);
	var = var / varianceFactor;
	mean = mean - meanShift*sqrt(var);
	borderHi = (float) (mean + sdNums*sqrt(var));
	borderLo = (float) (mean - sdNums*sqrt(var));



}

//.......................................................................................

void RepConfiguredOutlierCleaner::print()
{
	MLOG("BasicOutlierCleaner: signal: %d : doTrim %d trimMax %f trimMin %f : doRemove %d : removeMax %f removeMin %f\n",
		signalId, params.doTrim, trimMax, trimMin, params.doRemove, removeMax, removeMin, confFileName.c_str(), cleanMethod.c_str());
}
//=======================================================================================
// RuleBasedOutlierCleaner
//=======================================================================================

int RepRuleBasedOutlierCleaner::init(map<string, string>& mapper)
{

	init_defaults();
	set<string> rulesStrings;

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepRuleBasedOutlierCleaner::init]
		if (field == "signals") {
			boost::split(aff_signals, entry.second, boost::is_any_of(",")); // build list of  affected signals 
			for (auto sig : aff_signals)  // all affected are of course required
				req_signals.insert(sig);
		}
		else if (field == "time_channel") time_channel = med_stoi(entry.second);
		else if (field == "val_channel") val_channel = med_stoi(entry.second);
		else if (field == "addRequiredSignals")addRequiredSignals = med_stoi(entry.second)!=0;
		else if (field == "nrem_attr") nRem_attr = entry.second;
		else if (field == "nrem_suff") nRem_attr_suffix = entry.second;
		else if (field == "tolerance") tolerance = med_stof(entry.second);
		else if (field == "consideredRules") {
			boost::split(rulesStrings, entry.second, boost::is_any_of(","));
			for (auto& rule : rulesStrings) {
				int ruleNum = med_stoi(rule);
				consideredRules.push_back(ruleNum);
				if (ruleNum == 0)break;
			}
		}
		//! [RepRuleBasedOutlierCleaner::init]

	}


	for (auto& rule : rules2Signals) {
		if (std::find(consideredRules.begin(), consideredRules.end(), 0) != consideredRules.end() ||
			std::find(consideredRules.begin(), consideredRules.end(), rule.first) != consideredRules.end())
			continue;// rule remains
		else 
			rules2Signals.erase(rule.first);// rule removed
	}

	// add required signals according to rules that apply to affected signals
	for (auto& rule : rules2Signals) {
		for (auto& sig : aff_signals) {
			if (std::find(rule.second.begin(), rule.second.end(), sig) != rule.second.end()) {

				rulesToApply.push_back(rule.first);
				bool loopBreak = false;
				for (auto& reqSig : rule.second) {					
					bool found = false;
					for (auto& existReqSig : req_signals) {
						if (reqSig == existReqSig) {
							found = true;
							break;  //already there
						}
					}

					if (!found)
						if (addRequiredSignals)
							req_signals.insert(reqSig);    //add required signal
						else {
							rulesToApply.pop_back(); //We were asked not to load additional signals so ignore this rule
							loopBreak = true;
							break;
						}
				}
				if (loopBreak)break;
			}
		}
	}

	return MedValueCleaner::init(mapper);
}

void RepRuleBasedOutlierCleaner::init_attributes() {

	attributes.clear();
	if (!nRem_attr_suffix.empty()) {
		for (string signalName : aff_signals)
			attributes.push_back(signalName + "_" + nRem_attr_suffix);
	}

	if (!nRem_attr.empty()) attributes.push_back(nRem_attr);
}

void RepRuleBasedOutlierCleaner::set_signal_ids(MedDictionarySections& dict) {
	for (auto reqSig : req_signals)reqSignalIds.insert(dict.id(reqSig));
	for (auto affSig : aff_signals)affSignalIds.insert(dict.id(affSig));
}

void RepRuleBasedOutlierCleaner::init_tables(MedDictionarySections& dict, MedSignals& sigs) {

	rules_sids.resize(rulesToApply.size());
	affected_by_rules.resize(rulesToApply.size());

	for (int i = 0; i < rulesToApply.size(); i++) {
		// build set of the participating signals

		for (auto& sname : rules2Signals[rulesToApply[i]]) {
			int thisSid = dict.id(sname);
			rules_sids[i].push_back(thisSid);
			affected_by_rules[i].push_back(affSignalIds.find(thisSid) != affSignalIds.end());
		}
	}
}
	

int RepRuleBasedOutlierCleaner::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float> >& attributes_mat) {
	
	// get the signals
	map <int, UniversalSigVec> usvs;// from signal to its USV
	map <int, vector <int>> removePoints; // from signal id to its remove points
	
	// Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	differentVersionsIterator vit(rec, reqSignalIds);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		map <int, set<int>> removePoints; // from sid to indices to be removed

		// Clean 
		for (auto reqSigId : reqSignalIds) {
			rec.uget(reqSigId, iver, usvs[reqSigId]);
			removePoints[reqSigId] = {};
		}

		//Now loop on rules
		//printf("removePointsSize:%d\n", removePoints.size());

		for (int iRule = 0; iRule < rulesToApply.size(); iRule++) {
			int rule = rulesToApply[iRule];
			vector <UniversalSigVec>ruleUsvs;
			vector<int>& mySids = rules_sids[iRule];

			// build set of the participating signals
			for (int sid : mySids)
				ruleUsvs.push_back(usvs[sid]);

			bool signalEmpty = false;
			for (auto& thisUsv : ruleUsvs) // look for empty signals and skip the rule
				if (thisUsv.len == 0)signalEmpty = true;
			if (signalEmpty) continue; //skip this rule. one of the signals is empty (maybe was cleaned in earlier stage ) 

			// loop and find times where you have all signals
			vector <int>sPointer(mySids.size(), 0);
			int thisTime;
			for (sPointer[0] = 0; sPointer[0] < ruleUsvs[0].len; sPointer[0]++) {
				//printf("start loop %d %d \n", sPointer[0], ruleUsvs[0].len);
				thisTime = ruleUsvs[0].Time(sPointer[0], time_channel);
				if (time_points.size() != 0 && thisTime > time_points[iver])break;
				bool ok = true;
				for (int i = 1; i < mySids.size(); i++) {
					while (ruleUsvs[i].Time(sPointer[i], time_channel) < thisTime && sPointer[i] < ruleUsvs[i].len - 1)sPointer[i]++;
					//printf("before ok_check: %d %d %d %d %d %d\n", i, sPointer[0], sPointer[1], sPointer[2],thisTime, ruleUsvs[i].Time(sPointer[i], time_channel));
					if (ruleUsvs[i].Time(sPointer[i], time_channel) != thisTime ) {
						//printf("before ok_0: %d %d %d %d %d\n", rule, sPointer[0], sPointer[1], sPointer[2]);
						ok = 0;
						break;
					}
				}
				if (ok) {
					// if found all signals from same date eliminate doubles and take the last one for comparison
					for (int i = 0; i < mySids.size(); i++) 
						while(sPointer[i] < ruleUsvs[i].len - 1)
							if (ruleUsvs[i].Time(sPointer[i], time_channel) == ruleUsvs[i].Time(sPointer[i] + 1, time_channel)) {
								if(affected_by_rules[iRule][i])
										removePoints[mySids[i]].insert(sPointer[i]);
								sPointer[i]++;
							}
							else break;
					// check rule and mark for removement
					//printf("before apply: %d %d %d %d\n", rule, sPointer[0],sPointer[1],sPointer[2]);
					bool ruleFlagged = applyRule(rule, ruleUsvs, sPointer);
					/*
					printf("%d R: %d P: %d t: %d   ",ruleFlagged, rule, rec.pid, thisTime);
					for (int k = 0; k < sPointer.size(); k++)printf(" %f", ruleUsvs[k].Val(sPointer[k]));
					printf("\n");
					*/
					if (ruleFlagged) {
						
						for (int sIndex = 0; sIndex < mySids.size(); sIndex++)
							if (affected_by_rules[iRule][sIndex])
								removePoints[mySids[sIndex]].insert(sPointer[sIndex]);
					}
				}
			}
		}

		// Apply removals
		size_t nRemove = 0;
		int idx = 0;
		for (auto sig : affSignalIds) {
			vector <int> toRemove(removePoints[sig].begin(), removePoints[sig].end());
			vector <pair<int, float>>noChange;
			if (rec.update(sig, iver, val_channel, noChange, toRemove) < 0)
				return -1;
			if (!nRem_attr_suffix.empty() && !attributes_mat.empty()) {
				for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
					attributes_mat[pVersion][idx] = (float) toRemove.size();
				idx++;
			}
			nRemove += toRemove.size();
		}

		// Collect atttibutes
		if (!nRem_attr.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float) nRemove;
		}
	}

	return 0;


}

bool  RepRuleBasedOutlierCleaner::applyRule(int rule, vector <UniversalSigVec> ruleUsvs, vector<int> sPointer)
// apply the rule and return true if data is consistent with the rule
//ruleUsvs hold the signals in the order they appear in the rule in the rules2Signals above
{

	float left, right; // sides of the equality or inequality of the rule

	switch (rule) {
	case 1://BMI=Weight/Height^2*1e4
		if (ruleUsvs[2].Val(sPointer[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0]);
		right = ruleUsvs[1].Val(sPointer[1]) / ruleUsvs[2].Val(sPointer[2]) / ruleUsvs[2].Val(sPointer[2]) * (float)1e4;
		//printf("inputs %f %f\n", ruleUsvs[1].Val(sPointer[1]), ruleUsvs[2].Val(sPointer[2]));
		return (abs(left / right - 1) > tolerance);

	case 2://MCH=Hemoglobin/RBC*10
	case 3://MCV=Hematocrit/RBC*10
		if (ruleUsvs[2].Val(sPointer[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0]);
		right = ruleUsvs[1].Val(sPointer[1]) / ruleUsvs[2].Val(sPointer[2]) * 10;
		return(abs(left / right - 1) > tolerance);

	case 4://MCHC-M=MCH/MCV*100
		if (ruleUsvs[2].Val(sPointer[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0]);
		right = ruleUsvs[1].Val(sPointer[1]) / ruleUsvs[2].Val(sPointer[2]) * 100;
		return(abs(left / right - 1) > tolerance);

	case 11://HDL_over_nonHDL=HDL/NonHDLCholesterol
	case 12://HDL_over_Cholesterol=HDL/Cholesterol
		if (ruleUsvs[2].Val(sPointer[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0]);
		right =round( ruleUsvs[1].Val(sPointer[1]) / ruleUsvs[2].Val(sPointer[2])*10)/(float)10.; //resolution in THIN is 0.1
		return(abs(left / right - 1) > tolerance);

	case 6://MPV=Platelets_Hematocrit/Platelets
		if (ruleUsvs[2].Val(sPointer[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0]);
		right = ruleUsvs[1].Val(sPointer[1]) / ruleUsvs[2].Val(sPointer[2]);
		return(abs(left / right - 1) > tolerance);

	case 8://UrineAlbumin_over_Creatinine = UrineAlbumin / UrineCreatinine
		if (ruleUsvs[2].Val(sPointer[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0]);
		right = round(ruleUsvs[1].Val(sPointer[1]) / ruleUsvs[2].Val(sPointer[2])*10)/10;//resolution in THIN is 0.1
		return(abs(left / right - 1) > tolerance);

	case 13://HDL_over_LDL=HDL/LDL
	case 15://Cholesterol_over_HDL=Cholesterol/HDL
	case 18://LDL_over_HDL=LDL/HDL
		if (ruleUsvs[2].Val(sPointer[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0]);
		right = ruleUsvs[1].Val(sPointer[1]) / ruleUsvs[2].Val(sPointer[2]);
		return(abs(left / right - 1) > tolerance);

	case 5://Eosinophils#+Monocytes#+Basophils#+Lymphocytes#+Neutrophils#<=WBC
		left = ruleUsvs[0].Val(sPointer[0]) + ruleUsvs[1].Val(sPointer[1]) + ruleUsvs[2].Val(sPointer[2]) + ruleUsvs[3].Val(sPointer[3]) + ruleUsvs[4].Val(sPointer[4]);
		right = ruleUsvs[5].Val(sPointer[5]);
		return (left*(1 - tolerance) >= right);
	
	case 19://Albumin<=Protein_Total	
	case 21://NRBC<=RBC
	case 22://CHADS2<=CHADS2_VASC
		left = ruleUsvs[0].Val(sPointer[0]);
		right = ruleUsvs[1].Val(sPointer[1]);
		return(left*(1 - tolerance) >= right);

	case 7://UrineAlbumin <= UrineTotalProtein
	case 20://FreeT4<=T4
		left = ruleUsvs[0].Val(sPointer[0]);
		right = ruleUsvs[1].Val(sPointer[1]);
		return(left*(1 - tolerance) >= right*1000); // T4 is nmol/L free T4 is pmol/L ;  Albumin mg/L versus protein g/L

	case 9://LDL+HDL<=Cholesterol
		left = ruleUsvs[0].Val(sPointer[0]) + ruleUsvs[1].Val(sPointer[1]);
		right = ruleUsvs[2].Val(sPointer[2]);
		return (left*(1 - tolerance) > right);

	case 10://NonHDLCholesterol + HDL = Cholesterol
		left = ruleUsvs[0].Val(sPointer[0]) + ruleUsvs[1].Val(sPointer[1]);
		right = ruleUsvs[2].Val(sPointer[2]);
		return (abs(left / right - 1) > tolerance);

	case 14://HDL_over_LDL=1/LDL_over_HDL
	case 17://Cholesterol_over_HDL = 1 / HDL_over_Cholestrol
		if (ruleUsvs[2].Val(sPointer[1]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0]);
		right =(float) 1. / ruleUsvs[1].Val(sPointer[1]);
		return (abs(left / right - 1) > tolerance);

	default: assert(0); return false; // return is never executed but eliminates warning
	}
}

//=======================================================================================
// NbrsOutlierCleaner
//=======================================================================================
void RepNbrsOutlierCleaner::init_lists() {

	req_signals.insert(signalName);
	aff_signals.insert(signalName);
}

//.......................................................................................
int RepNbrsOutlierCleaner::init(map<string, string>& mapper)
{
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepNbrsOutlierCleaner::init]
		if (field == "signal") { signalName = entry.second; }
		else if (field == "time_channel") time_channel = med_stoi(entry.second);
		else if (field == "val_channel") val_channel = med_stoi(entry.second);
		else if (field == "nbr_time_unit") nbr_time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "nbr_time_width") nbr_time_width = med_stoi(entry.second);
		else if (field == "nrem_attr") nRem_attr = entry.second;
		else if (field == "ntrim_attr") nTrim_attr == entry.second;
		else if (field == "nrem_suff") nRem_attr_suffix = entry.second;
		else if (field == "ntrim_suff") nTrim_attr_suffix = entry.second;
		//! [RepNbrsOutlierCleaner::init]
	}

	init_lists();
	return MedValueCleaner::init(mapper);
}

// init attributes list
//.......................................................................................
void RepNbrsOutlierCleaner::init_attributes() {

	string _signal_name = signalName;
	if (val_channel != 0) _signal_name += "_" + to_string(val_channel);

	attributes.clear();
	if (!nRem_attr.empty()) attributes.push_back(nRem_attr);
	if (!nRem_attr_suffix.empty()) attributes.push_back(_signal_name + "_" + nRem_attr_suffix);

	if (!nTrim_attr.empty()) attributes.push_back(nTrim_attr);
	if (!nTrim_attr_suffix.empty()) attributes.push_back(_signal_name + "_" + nTrim_attr_suffix);
}

// Learn bounds
//.......................................................................................
int RepNbrsOutlierCleaner::_learn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_cleaners) {

	if (params.type == VAL_CLNR_ITERATIVE)
		return iterativeLearn(rep, samples, prev_cleaners);
	else if (params.type == VAL_CLNR_QUANTILE)
		return quantileLearn(rep, samples, prev_cleaners);
	else {
		MERR("Unknown cleaning type %d\n", params.type);
		return -1;
	}
}

//.......................................................................................
int RepNbrsOutlierCleaner::iterativeLearn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_cleaners) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Get all values
	vector<float> values;
	get_values(rep, samples, signalId, time_channel, val_channel, params.range_min, params.range_max, values, prev_cleaners);

	int rc = get_iterative_min_max(values);
	return rc;
}

//.......................................................................................
int RepNbrsOutlierCleaner::quantileLearn(MedPidRepository& rep, MedSamples& samples, vector<RepProcessor *>& prev_cleaners) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Get all values
	vector<float> values;
	get_values(rep, samples, signalId, time_channel, val_channel, params.range_min, params.range_max, values, prev_cleaners);

	return get_quantile_min_max(values);
}

// Clean
//.......................................................................................
int  RepNbrsOutlierCleaner::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float> >& attributes_mat) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int len;
	allVersionsIterator vit(rec, signalId);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

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

			if (time_points.size() != 0 && itime > time_points[iver]) {
				verLen = i;
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

		// Collect atttibutes
		int idx = 0;
		if (!nRem_attr.empty() && !attributes_mat.empty())
			attributes_mat[iver][idx++] = (float) nRemove;
		if (!nRem_attr_suffix.empty() && !attributes_mat.empty())
			attributes_mat[iver][idx++] = (float)nRemove;
		if (!nTrim_attr.empty() && !attributes_mat.empty())
			attributes_mat[iver][idx++] = (float)nChange;
		if (!nTrim_attr_suffix.empty() && !attributes_mat.empty())
			attributes_mat[iver][idx++] = (float)nChange;
	}

	return 0;

}

//.......................................................................................
void RepNbrsOutlierCleaner::print()
{
MLOG("RepNbrsOutlierCleaner: signal: %d : doTrim %d trimMax %f trimMin %f : doRemove %d : removeMax %f : removeMin %f : nbrsMax %f : nbrsMin %f\n",
	signalId, params.doTrim, trimMax, trimMin, params.doRemove, removeMax, removeMin,nbrsMax, nbrsMin);
}

//=======================================================================================
// RepMinimalReq - check requirement and set attributes accordingly
//=======================================================================================
//.......................................................................................
int RepCheckReq::init(map<string, string>& mapper)
{

	time_channels.clear();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepMinimalReq::init]
		if (field == "signals") 
			boost::split(signalNames, entry.second, boost::is_any_of(","));
		else if (field == "time_channels") {
			vector<string> channels;
			boost::split(channels, entry.second, boost::is_any_of(","));
			for (string& channel : channels)
				time_channels.push_back(stoi(channel));
		}
		else if (field == "time_unit") 
			window_time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "win_from") 
			win_from = med_stoi(entry.second);
		else if (field == "win_to") 
			win_to = med_stoi(entry.second);
		else if (field == "attr")
			attrName = entry.second;
		//! [RepMinimalReq::init]
	}

	// Take care of time-channels
	if (time_channels.size() == 0)
		time_channels.push_back(0);
	if (time_channels.size() == 1) {
		int channel = time_channels[0];
		time_channels.resize(signalNames.size(), channel);
	}

	init_lists();

	return 0;
}

//.......................................................................................
void RepCheckReq::set_signal_ids(MedDictionarySections& dict) {

	signalIds.resize(signalNames.size());
	for (int i = 0; i < signalNames.size(); i++)
		signalIds[i] = (dict.id(signalNames[i]));
}

//.......................................................................................
void RepCheckReq::init_lists() {
	req_signals.insert(signalNames.begin(), signalNames.end());
}

//.......................................................................................
void RepCheckReq::init_tables(MedDictionarySections& dict, MedSignals& sigs) {

	sig_time_units.resize(signalIds.size());
	for (int i = 0; i < signalIds.size(); i++)
		sig_time_units[i] = sigs.Sid2Info[signalIds[i]].time_unit;

}

//.......................................................................................
int RepCheckReq::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {

	// Should we do anything ?
	if (attributes_mat.empty())
		return 0;

	// Sanity checks
	if (signalIds.size() != signalNames.size()) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	for (int signalId : signalIds) {
		if (signalId == -1) {
			MERR("Uninitialized signalId\n");
			return -1;
		}
	}
	
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	// Loop on versions
	set<int> _signalIds(signalIds.begin(), signalIds.end());
	allVersionsIterator vit(rec, _signalIds);
	vector<UniversalSigVec> usvs(signalIds.size());


	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		// NOTE : This is not perfect : we assume that the samples' time-unit is the default time unit
		int time_point = med_time_converter.convert_times(global_default_time_unit, window_time_unit, time_points[iver]);
		
		int nMissing = 0;
		for (int i = 0; i < signalIds.size(); i++) {
			rec.uget(signalIds[i], iver, usvs[i]);

			bool found = false;
			for (int j = 0; j < usvs[i].len; j++) {
				if (usvs[i].Time(j,time_channels[i]) > med_time_converter.convert_times(window_time_unit, sig_time_units[i], time_point - win_from))
					break;
				if (usvs[i].Time(j, time_channels[i]) >= med_time_converter.convert_times(window_time_unit, sig_time_units[i], time_point - win_to)) {
					found = true;
					break;
				}
			}
			
			if (!found)
				nMissing++;
		}

		// Set attribute
		attributes_mat[iver][0] = (float) nMissing;
	}

	return 0;
}

//=======================================================================================
// RepSimValHandler - handle multiple simultanous values
//=======================================================================================
// Fill req- and aff-signals vectors
//.......................................................................................
void RepSimValHandler::init_lists() {

	req_signals.insert(signalName);
	aff_signals.insert(signalName);
}

// Init from map
//.......................................................................................
int RepSimValHandler::init(map<string, string>& mapper)
{
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepSimValHandler::init]
		if (field == "signal") { signalName = entry.second; }
		else if (field == "type") handler_type = get_sim_val_handle_type(entry.second);
		else if (field == "time_channels") {
			vector<string> channels;
			boost::split(channels, entry.second, boost::is_any_of(","));
			for (auto& channel : channels)
				time_channels.push_back(stoi(channel));
		}
		else if (field == "attr") nHandle_attr = entry.second;
		else if (field == "suff") nHandle_attr_suffix = entry.second;
		//! [RepSimValHandler::init]
	}

	init_lists();
	return 0;
}

// init attributes list
//.......................................................................................
void RepSimValHandler::init_attributes() {

	string _signal_name = signalName;

	if (time_channels.size() > 1 || (time_channels.size()==1 && time_channels[0] != 0)) {
		vector<string> channels_s(time_channels.size());
		for (unsigned int i = 0; i < channels_s.size(); i++)
			channels_s[i] = to_string(time_channels[i]);
		_signal_name += "_" + boost::join(channels_s, "_");
	}
	

	attributes.clear();
	if (!nHandle_attr.empty()) attributes.push_back(nHandle_attr);
	if (!nHandle_attr_suffix.empty()) attributes.push_back(_signal_name + "_" + nHandle_attr_suffix);
}

// name to SimValHandleTypes
//.......................................................................................
SimValHandleTypes RepSimValHandler::get_sim_val_handle_type(string& name) {

	if (name == "first" || name == "first_val")
		return SIM_VAL_FIRST_VAL;
	else if (name == "last" || name == "last_val")
		return SIM_VAL_LAST_VAL;
	else if (name == "mean" || name == "avg")
		return SIM_VAL_MEAN;
	else if (name == "rem" || name == "remvoe")
		return SIM_VAL_REM;
	else if (name == "rem_diff" || name == "remove_diff")
		return SIM_VAL_REM_DIFF;
	else
		MTHROW_AND_ERR("Unkwnon sim_val_hand_type \'%s\'\n", name.c_str());

}

// Get time-channels (if empty)
//.......................................................................................
void RepSimValHandler::init_tables(MedDictionarySections& dict, MedSignals& sigs) {

	if (time_channels.empty()) {
		int n = sigs.Sid2Info[signalId].n_time_channels;
		time_channels.resize(n);
		for (int i = 0; i < n; i++)
			time_channels[i] = i;
	}

	nValChannels = sigs.Sid2Info[signalId].n_val_channels;
}

// Apply
//.......................................................................................
int  RepSimValHandler::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float> >& attributes_mat) {

	// Sanity check
	if (signalId == -1) {
		MERR("Uninitialized signalId\n");
		return -1;
	}

	// Check that we have the correct number of dynamic-versions : one per time-point (if given)
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int len;
	differentVersionsIterator vit(rec, signalId);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		// Do it 
		rec.uget(signalId, iver, rec.usv); // get into the internal usv obeject - this statistically saves init time

		len = rec.usv.len;
		vector<int> remove(len);
		vector<pair<int, vector<float>>> change(len);
		int nRemove = 0, nChange = 0, nTimes = 0;

		// Collect
		int start = 0, end = 0;
		for (int i = 1; i < len; i++) {

			// No need to clean past the latest relevant time-point (valid only when using a single time-channel == 0)
			if (time_points.size() != 0 && time_channels.size() == 1 && time_channels[0] == 0 && rec.usv.Time(i) > time_points[iver]) break;

			// Are we simultanous ?
			bool sim = true;
			for (int channel : time_channels) {
				if (rec.usv.Time(i, channel) != rec.usv.Time(i - 1, channel)) {
					sim = false;
					break;
				}
			}

			if (!sim) {
				if (end > start)
					handle_block(start, end, rec.usv, remove, nRemove, change, nChange, nTimes);
				start = end = i;
			}
			else
				end++;
		}

		// Handle last block
		if (end > start)
			handle_block(start, end, rec.usv, remove, nRemove, change, nChange, nTimes);
	

		// Apply removals + changes
		change.resize(nChange);
		remove.resize(nRemove);
		if (rec.update(signalId, iver, change, remove) < 0)
			return -1;

		// Collect atttibutes
		int idx = 0;
		if (!nHandle_attr.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float)nTimes;
			idx++;
		}
		
		if (!nHandle_attr_suffix.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float)nTimes;
			idx++;
		}
	}

	return 0;

}

// Utility : handle a block
//.......................................................................................
void RepSimValHandler::handle_block(int start, int end, UniversalSigVec& usv, vector<int>& remove, int& nRemove, vector<pair<int, vector<float>>>& change, int& nChange, int& nTimes) {

	if (handler_type == SIM_VAL_FIRST_VAL) {
		for (int j = start + 1; j <= end; j++)
			remove[nRemove++] = j;
		nTimes++;
	}
	else if (handler_type == SIM_VAL_LAST_VAL) {
		for (int j = start; j < end; j++)
			remove[nRemove++] = j;
		nTimes++;
	}
	else if (handler_type == SIM_VAL_MEAN) {
		vector<float> sums(nValChannels, 0);
		for (int j = start; j < end; j++) {
			for (int iChannel=0; iChannel < nValChannels; iChannel++)
				sums[iChannel] += usv.Val(j,iChannel);
			remove[nRemove++] = j;
		}
		pair<int, vector<float>> newChange;
		newChange.first = end;
		newChange.second.resize(nValChannels);
		for (int iChannel=0; iChannel < nValChannels; iChannel++)
			newChange.second[iChannel] = (sums[iChannel] + usv.Val(end, iChannel)) / (end + 1 - start);
		change[nChange++] = newChange;
		nTimes++;
	}
	else if (handler_type == SIM_VAL_REM) {
		for (int j = start; j <= end; j++)
			remove[nRemove++] = j;
		nTimes++;
	}
	else if (handler_type == SIM_VAL_REM_DIFF) {
		bool rem = false;
		for (int j = start + 1; j <= end; j++) {
			for (int channel = 0; channel < nValChannels; channel++) {
				if (usv.Val(j, channel) != usv.Val(j - 1, channel)) {
					rem = true;
					break;
				}
			}
			if (rem) break;
		}

		if (rem) {
			for (int j = start; j <= end; j++)
				remove[nRemove++] = j;
			nTimes++;
		}
	}
}

//=======================================================================================
// RepCalcSimpleSignals - calculators with no learning stage, can be parametric.
//=======================================================================================
//.......................................................................................
int RepCalcSimpleSignals::init(map<string, string>& mapper)
{
	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "calculator") calculator = entry.second;
		else if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "coeffs") {
			vector<string> fields;
			boost::split(fields, entry.second, boost::is_any_of(",:"));
			coeff.clear();
			for (auto &f : fields) coeff.push_back(stof(f));
		}
		else if (field == "names") {
			boost::split(V_names, entry.second, boost::is_any_of(",:"));
		}
		else if (field == "signals") {
			boost::split(signals, entry.second, boost::is_any_of(",:"));
		}
		else if (field == "timer") timer_signal = entry.second;
		else if (field == "time_step") time_step_str = entry.second;
	}

	calc_type = get_calculator_type(calculator);

	if (calc_type != CALC_TYPE_UNDEF) {

		//time_step
		if (!time_step_str.empty()) {
			time_step = stoi(time_step_str);
			if (time_step <= 0) {
				MERR("Non-positive time_step: %d\n", time_step);
				return -1;
			}
		}

		// Timer
		// Given when not required
		if (!timer_signal.empty() && calc_without_timers.find(calculator) != calc_without_timers.end()) {
			MERR("Caclulator %s must not be given a timer signal\n", calculator.c_str());
			return -1;
		}

		// for CALC_TYPE_HOSP_PROCESSOR default timer is signal
		if (calc_type == CALC_TYPE_HOSP_PROCESSOR) {
			if (!timer_signal.empty() && timer_signal != signals[0]) {
				MERR("Calculator %s timer must be the same as signal\n", calculator.c_str());
				return -1;
			}
			timer_signal = signals[0];
		}

		// add required signals depending on the actual calculator we run
		// might be overidden from json
		req_signals.clear();
		if (signals.size() == 0)
			signals = calc2req_sigs.find(calculator)->second;

		for (auto & req_s : signals)
			req_signals.insert(req_s);

		if ((!timer_signal.empty()) && req_signals.find(timer_signal) == req_signals.end())
			req_signals.insert(timer_signal);
						
		// add coefficients if needed
		if (coeff.size() == 0)
			coeff = calc2coeffs.find(calculator)->second;
		size_t c_size = calc2coeffs.find(calculator)->second.size();
		if (coeff.size() != c_size) {
			MERR("RepCalcSigs: ERROR: calculator %s , expecting %d coefficients and got only %d\n", calculator.c_str(), c_size, coeff.size());
			return -1;
		}

		// add V_names
		if (V_names.size() == 0) {
			for (auto &vsig : calc2virtual.find(calculator)->second)
				V_names.push_back(vsig.first);
		}
		size_t v_size = calc2virtual.find(calculator)->second.size();
		if (V_names.size() != v_size) {
			MERR("RepCalcSigs: ERROR: calculator %s , expecting %d virtual names and got only %d\n", calculator.c_str(), v_size, V_names.size());
			return -1;
		}

		// add V_types
		for (auto &vsig : calc2virtual.find(calculator)->second)
			V_types.push_back(vsig.second);

		// add names to required, affected and virtual_signals
		aff_signals.clear();
		virtual_signals.clear();
		for (int i=0; i<V_names.size(); i++) {
			aff_signals.insert(V_names[i]);
			virtual_signals.push_back({ V_names[i], V_types[i] });
		}
		for (int i = 0; i < signals.size(); i++)
			req_signals.insert(signals[i]);

		return 0;
	}

	MERR("RepCalcSigs: ERROR: calculator %s not defined\n", calculator.c_str());
	return -1;
}

mutex RepCalcSimpleSignals_init_tables_mutex;
//.......................................................................................
void RepCalcSimpleSignals::init_tables(MedDictionarySections& dict, MedSignals& sigs)
{
	lock_guard<mutex> guard(RepCalcSimpleSignals_init_tables_mutex);

	V_ids.clear();
	sigs_ids.clear();
	for (auto &vsig : V_names)
		V_ids.push_back(dict.id(vsig));
	// In the next loop it is VERY important to go over items in the ORDER they are given in calc2req
	// This is since we create a vector of sids (sigs_ids) that matches it exactly, and enables a much
	// more efficient code without going to this map for every pid. (See for example the egfr calc function)
	for (auto &rsig : signals) 
		sigs_ids.push_back(dict.id(rsig));

	if (!timer_signal.empty()) {
		timer_signal_id = dict.id(timer_signal);
		signals_time_unit = sigs.Sid2Info[timer_signal_id].time_unit;
	}

	//hack for saving run-time: code for african american is added on the fly as a coeff,
	//to be used by the calculator
	if (calc_type == CALC_TYPE_HOSP_IS_AFRICAN_AMERICAN) {
		int african_american_dict_id = dict.id("BLACK/AFRICAN AMERICAN");
		coeff = {(float)african_american_dict_id };
	}

}

//.......................................................................................
void RepCalcSimpleSignals::add_virtual_signals(map<string, int> &_virtual_signals)
{
	for (int i=0; i<V_names.size(); i++)
		_virtual_signals[V_names[i]] = V_types[i];
	//for (auto &vsig : calc2virtual.find(calculator)->second)
	//	_virtual_signals[vsig.first] = vsig.second;
}

//.......................................................................................
int RepCalcSimpleSignals::get_calculator_type(const string &calc_name)
{
	if (calc2type.find(calc_name) != calc2type.end())
		return calc2type.find(calc_name)->second;
	return CALC_TYPE_UNDEF;
}

//.......................................................................................
int RepCalcSimpleSignals::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat)
{
	//handle special calculations
	if (calc_type == CALC_TYPE_EGFR)
		return _apply_calc_eGFR(rec, time_points);

	if (calc_type == CALC_TYPE_DEBUG) 
		return _apply_calc_debug(rec, time_points);

	if (calc_type == CALC_TYPE_HOSP_24H_URINE_OUTPUT) 
		return _apply_calc_24h_urine_output(rec, time_points);

	if (calc_type == CALC_TYPE_LOG)
		return _apply_calc_log(rec, time_points);

	//handle calculation that are done by first extrapolating each signal to the required 
	//time points and then performing a pointwise calculation for the values in each time point.
	float(*calcFunc)(const vector<float>&, const vector<float>&) = NULL;

	switch (calc_type) {
		case CALC_TYPE_HOSP_PROCESSOR: 
			calcFunc = identity; 
			if (signals.size() != 1) {
				//MERR("calc_hosp_processor calculator requires exactly one input signal. Found %d\n", (int)(signals.size()));
				return -1;
			}
			break;
		case CALC_TYPE_HOSP_MELD: calcFunc = calc_hosp_MELD; break;
		case CALC_TYPE_HOSP_BMI: calcFunc = calc_hosp_BMI; break;
		case CALC_TYPE_HOSP_APRI: calcFunc = calc_hosp_APRI; break;
		case CALC_TYPE_HOSP_SIDA: calcFunc = calc_hosp_SIDA; break;
		case CALC_TYPE_HOSP_PaO2_FiO2_RATIO: calcFunc = calc_hosp_PaO2_FiO2_ratio; break;
		case CALC_TYPE_HOSP_IS_AFRICAN_AMERICAN: calcFunc = calc_hosp_is_african_american; break;
		case CALC_TYPE_HOSP_SOFA_NERVOUS: calcFunc = calc_hosp_SOFA_nervous; break;
		case CALC_TYPE_HOSP_SOFA_LIVER: calcFunc = calc_hosp_SOFA_liver; break;
		case CALC_TYPE_HOSP_SOFA_COAGULATION: calcFunc = calc_hosp_SOFA_coagulation; break;
		case CALC_TYPE_HOSP_DOPAMINE_PER_KG: calcFunc = calc_hosp_dopamine_per_kg; break;
		case CALC_TYPE_HOSP_EPINEPHRINE_PER_KG: calcFunc = calc_hosp_epinephrine_per_kg; break;
		case CALC_TYPE_HOSP_NOREPINEPHRINE_PER_KG: calcFunc = calc_hosp_norepinephrine_per_kg; break;
		case CALC_TYPE_HOSP_DOBUTAMINE_PER_KG: calcFunc = calc_hosp_dobutamine_per_kg; break;
		case CALC_TYPE_HOSP_QSOFA: calcFunc = calc_hosp_qSOFA; break;
		case CALC_TYPE_HOSP_SIRS: calcFunc = calc_hosp_SIRS; break;
		case CALC_TYPE_HOSP_PRESSURE_ADJUSTED_HR: calcFunc = calc_hosp_pressure_adjusted_hr; break;
		case CALC_TYPE_HOSP_MODS: calcFunc = calc_hosp_MODS; break;
		case CALC_TYPE_HOSP_SHOCK_INDEX: calcFunc = calc_hosp_shock_index; break;
		case CALC_TYPE_HOSP_PULSE_PRESSURE: calcFunc = calc_hosp_pulse_pressure; break;
		case CALC_TYPE_HOSP_EFGR: calcFunc = calc_hosp_eGFR; break;
		case CALC_TYPE_HOSP_SOFA_RESPIRATORY: calcFunc = calc_hosp_SOFA_respiratory; break;
		case CALC_TYPE_HOSP_SOFA_RENAL: calcFunc = calc_hosp_SOFA_renal; break;
		case CALC_TYPE_HOSP_SOFA_CARDIO: calcFunc = calc_hosp_SOFA_cardio; break;
		case CALC_TYPE_HOSP_SOFA: calcFunc = calc_hosp_SOFA; break;		
	}

	if (calcFunc) {
		int rv = _apply_calc_hosp_pointwise(rec, time_points, calcFunc);
		return rv;
	}

	//handle time-dependent calculations
	//calculation are done by first extrapolating each signal to the required 
	//time points and then performing a pointwise calculation for the values in each time point.
	//in contrast to previous signals, here the calculation does depend on the orignal times 
	//of the signals and their reference to the requested times
	float(*calcTimeFunc)(const vector<pair<int, float> >&, int, const vector<float>&) = NULL;

	switch (calc_type) {	
	case CALC_TYPE_HOSP_BP_SYS: calcTimeFunc = interleave; break;
	case CALC_TYPE_HOSP_BP_DIA: calcTimeFunc = interleave; break;
	}

	if (calcTimeFunc)
		return _apply_calc_hosp_time_dependent_pointwise(rec, time_points, calcTimeFunc, false); //use past/future observations

	switch (calc_type) {
	case CALC_TYPE_HOSP_IS_MECHANICALLY_VENTILATED: calcTimeFunc = anySeenRecently; break; //special function
	}

	if (calcTimeFunc)
		return _apply_calc_hosp_time_dependent_pointwise(rec, time_points, calcTimeFunc, true); //use only past obeservations

	return -1;
}

//=======================================================================================
// Utility Functions
//=======================================================================================
//.......................................................................................
// Get values of a signal from a set of ids
int get_values(MedRepository& rep, MedSamples& samples, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values, vector<RepProcessor *>& prev_processors)
{

	// Required signals
	vector<int> req_signal_ids_v;
	vector<unordered_set<int> > current_required_signal_ids(prev_processors.size());
	vector<FeatureGenerator *> noGenerators ;
	unordered_set<int> extra_req_signal_ids = { signalId };
	handle_required_signals(prev_processors, noGenerators, extra_req_signal_ids, req_signal_ids_v, current_required_signal_ids);

	PidDynamicRec rec;
	UniversalSigVec usv;

	bool signalIsVirtual = (bool)(rep.sigs.Sid2Info[signalId].virtual_sig != 0);

	for (MedIdSamples& idSamples : samples.idSamples) {

		int id = idSamples.id;
	
		vector<int> time_points;
		// Special care for virtual signals - use samples 
		if (signalIsVirtual) {
			time_points.resize(idSamples.samples.size());
			for (size_t i = 0; i < time_points.size(); i++)
				time_points[i] = idSamples.samples[i].time;
		}
		else {
			// Get signal
			rep.uget(id, signalId, usv);

			time_points.resize(usv.len);
			for (int i = 0; i < usv.len; i++)
				time_points[i] = usv.Time(i, time_channel);
		}

		// Nothing to do if empty ...
		if (time_points.empty())
			continue;

		if (prev_processors.size()) {

			// Init Dynamic Rec
			rec.init_from_rep(std::addressof(rep), id, req_signal_ids_v, (int)time_points.size());

			// Process at all time-points
			vector<vector<float>> dummy_attributes_mat;
			for (size_t i = 0; i < prev_processors.size(); i++)
				prev_processors[i]->conditional_apply(rec, time_points, current_required_signal_ids[i], dummy_attributes_mat);
		
			// If virtual - we need to get the signal now
			if (signalIsVirtual)
				rec.uget(signalId, 0, usv);

			// Collect
			int iVersion = 0;
			rec.uget(signalId, iVersion, rec.usv);

			for (int i = 0; i < usv.len; i++) {
				// Get a new version if we past the current one
				if (usv.Time(i) > time_points[iVersion]) {
					iVersion++;
					if (iVersion == rec.get_n_versions())
						break;
					rec.uget(signalId, iVersion, rec.usv);
				}

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
int get_values(MedRepository& rep, MedSamples& samples, int signalId, int time_channel, int val_channel, float range_min, float range_max, vector<float>& values) {
	vector<RepProcessor *> temp;
	return get_values(rep, samples, signalId, time_channel, val_channel, range_min, range_max, values, temp);
}
