#define _CRT_SECURE_NO_WARNINGS

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#include "RepProcess.h"
#include <MedUtils/MedUtils/MedUtils.h>
#include <cmath>

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
	else if (processor_name == "signal_rate")
		return REP_PROCESS_SIGNAL_RATE;
	else if (processor_name == "combine")
		return REP_PROCESS_COMBINE;
	else if (processor_name == "split")
		return REP_PROCESS_SPLIT;
	else if (processor_name == "aggregation_period")
		return REP_PROCESS_AGGREGATION_PERIOD;
	else if (processor_name == "basic_range_cleaner" || processor_name == "range_cln")
		return REP_PROCESS_BASIC_RANGE_CLEANER;
	else if (processor_name == "aggregate")
		return REP_PROCESS_AGGREGATE;
	else if (processor_name == "limit_history" || processor_name == "history_limit")
		return REP_PROCESS_HISTORY_LIMIT;
	else if (processor_name == "create_registry")
		return REP_PROCESS_CREATE_REGISTRY;
	else
		return REP_PROCESS_LAST;
}

// rep processors get a new derived class
//.......................................................................................
void *RepProcessor::new_polymorphic(string dname)
{
	CONDITIONAL_NEW_CLASS(dname, RepMultiProcessor);
	CONDITIONAL_NEW_CLASS(dname, RepBasicOutlierCleaner);
	CONDITIONAL_NEW_CLASS(dname, RepNbrsOutlierCleaner);
	CONDITIONAL_NEW_CLASS(dname, RepConfiguredOutlierCleaner);
	CONDITIONAL_NEW_CLASS(dname, RepRuleBasedOutlierCleaner);
	CONDITIONAL_NEW_CLASS(dname, RepCalcSimpleSignals);
	CONDITIONAL_NEW_CLASS(dname, RepPanelCompleter);
	CONDITIONAL_NEW_CLASS(dname, RepCheckReq);
	CONDITIONAL_NEW_CLASS(dname, RepSimValHandler);
	CONDITIONAL_NEW_CLASS(dname, RepSignalRate);
	CONDITIONAL_NEW_CLASS(dname, RepCombineSignals);
	CONDITIONAL_NEW_CLASS(dname, RepSplitSignal);
	CONDITIONAL_NEW_CLASS(dname, RepAggregationPeriod);
	CONDITIONAL_NEW_CLASS(dname, RepBasicRangeCleaner);
	CONDITIONAL_NEW_CLASS(dname, RepAggregateSignal);
	CONDITIONAL_NEW_CLASS(dname, RepHistoryLimit);
	CONDITIONAL_NEW_CLASS(dname, RepCreateRegistry);
	return NULL;
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
	else if (processor_type == REP_PROCESS_SIGNAL_RATE)
		return new RepSignalRate;
	else if (processor_type == REP_PROCESS_COMBINE)
		return new RepCombineSignals;
	else if (processor_type == REP_PROCESS_SPLIT)
		return new RepSplitSignal;
	else if (processor_type == REP_PROCESS_AGGREGATION_PERIOD)
		return new RepAggregationPeriod;
	else if (processor_type == REP_PROCESS_BASIC_RANGE_CLEANER)
		return new RepBasicRangeCleaner;
	else if (processor_type == REP_PROCESS_AGGREGATE)
		return new RepAggregateSignal;
	else if (processor_type == REP_PROCESS_HISTORY_LIMIT)
		return new RepHistoryLimit;
	else if (processor_type == REP_PROCESS_CREATE_REGISTRY)
		return new RepCreateRegistry;
	else
		return NULL;

}

// Initialization : given processor type and intialization string
//.......................................................................................
RepProcessor * RepProcessor::make_processor(RepProcessorTypes processor_type, string init_string) {

	//MLOG("Processor type is %d\n", (int)processor_type);
	RepProcessor *newRepProcessor = make_processor(processor_type);
	if (newRepProcessor->init_from_string(init_string) < 0)
		MTHROW_AND_ERR("Cannot init RepProcessor of type %d with init string \'%s\'\n", processor_type, init_string.c_str());

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

	MLOG_D("RepProcessor::filter filtering out processor of type %d, affected signals: ", processor_type);
	for (string signal : aff_signals)
		MLOG_D("[%s] ", signal.c_str());
	MLOG_D("\n");
	return true;

}

// Apply processing on a single PidDynamicRec at a set of time-points given by samples
//.......................................................................................
int RepProcessor::apply(PidDynamicRec& rec, MedIdSamples& samples) {

	vector<int> time_points(samples.samples.size());
	for (unsigned int i = 0; i < time_points.size(); i++)
		time_points[i] = samples.samples[i].time;

	vector<vector<float>> attributes_mat(time_points.size(), vector<float>(attributes.size(), 0));
	int rc = apply(rec, time_points, attributes_mat);

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

	vector<int> time_points;
	samples.get_times(time_points);

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
		MLOG("%s :: RP type %d(%s) : required(%d): ", pref.c_str(), processor_type, my_class_name().c_str(), req_signals.size());
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
	bool did_something = false;
	for (auto& processor : processors) {
		if (!processor->filter(neededSignals))
			filtered.push_back(processor);
		else {
			delete processor;
			processor = NULL;
			did_something = true;
		}
	}
	if (did_something)
		MLOG_D("Filtering uneeded rep_processors in RepMultiProcessor. left with %zu processors out of %zu\n",
			filtered.size(), processors.size());

	if (filtered.empty()) {
		MLOG_D("RepMultiProcessor::filter filtering out processor of type %d\n", processor_type);
		processors.clear();
		return true;
	}
	else {
		processors.swap(filtered);
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
		rc[j] = processors[j]->conditional_learn(rep, samples, prev_processors, neededSignalIds);
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

//.......................................................................................
int RepMultiProcessor::_apply_simple(PidDynamicRec& rec, vector<int>& time_points)
{
	for (auto p : processors) {
		if ((p->_apply_simple(rec, time_points)) < 0)
			return -1;
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
				attributes_pos[attributes.back()] = (int)attributes.size() - 1;
			}
			attributes_map[i][j] = attributes_pos[processors[i]->attributes[j]];
		}
	}
}

//.......................................................................................
void RepMultiProcessor::dprint(const string &pref, int rp_flag)
{
	if (rp_flag > 0) {
		MLOG("%s :: RP MULTI(%d) -->\n", pref.c_str(), processors.size());
		for (auto& proc : processors) {
			proc->dprint(pref + "->Multi", rp_flag);
		}
	}
}

void RepMultiProcessor::register_virtual_section_name_id(MedDictionarySections& dict) {
	for (size_t i = 0; i < processors.size(); ++i)
		processors[i]->register_virtual_section_name_id(dict);
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
	int rc = get_iterative_min_max(values);

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
				attributes_mat[pVersion][idx] = (float)nRemove;
			idx++;
		}

		if (!nRem_attr_suffix.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float)nRemove;
			idx++;
		}

		if (!nTrim_attr.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float)nChange;
			idx++;
		}

		if (!nTrim_attr_suffix.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float)nChange;
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
	while (getline(infile, thisLine)) {
		boost::trim(thisLine);
		if (thisLine.empty() || thisLine.at(0) == '#')
			continue; //skip empty line

		vector<string> f;
		boost::split(f, thisLine, boost::is_any_of(","));
		if (f.size() != 8) {
			fprintf(stderr, "Wrong field count in  %s (%s : %zd) \n", confFileName.c_str(), thisLine.c_str(), f.size());
			infile.close();
			return -1;
		}

		thisRecord.confirmedLow = thisRecord.logicalLow = (float)atof(f[1].c_str());
		thisRecord.confirmedHigh = thisRecord.logicalHigh = (float)atof(f[2].c_str());

		thisRecord.distLow = f[4];
		thisRecord.distHigh = f[6];
		thisRecord.val_channel = stoi(f[7]);
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

void RepConfiguredOutlierCleaner::set_signal_ids(MedDictionarySections& dict) {
	RepBasicOutlierCleaner::set_signal_ids(dict); //call base class init
	//fetch val_channel from file
	if (outlierParams.find(signalName) == outlierParams.end())
		MTHROW_AND_ERR("RepConfiguredOutlierCleaner : ERROR: Signal %s not supported by conf_cln\n", signalName.c_str());
	val_channel = outlierParams.at(signalName).val_channel;
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

	int start = (int)round(filteredValues.size()*margin[0]);
	int stop = (int)round(filteredValues.size()*margin[1]);
	for (vector<float>::iterator el = filteredValues.begin() + start; el < filteredValues.begin() + stop; el++) {

		sum += *el;
		sumsq += *el* *el;
	}
	double mean = sum / (stop - start);
	double var = sumsq / (stop - start) - mean * mean;
	//printf("sum %f sumsq %f  stop %d start %d\n", sum, sumsq, stop, start);
	var = var / varianceFactor;
	mean = mean - meanShift * sqrt(var);
	borderHi = (float)(mean + sdNums * sqrt(var));
	borderLo = (float)(mean - sdNums * sqrt(var));



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

void RepRuleBasedOutlierCleaner::parse_rules_signals(const string &path) {
	ifstream fr(path);
	if (!fr.good())
		MTHROW_AND_ERR("Error RepRuleBasedOutlierCleaner::parse_rules_signals - can't read file %s\n", path.c_str());
	string line;
	while (getline(fr, line)) {
		boost::trim(line);
		if (line.empty() || line[0] == '#')
			continue;
		vector<string> tokens, list_of_sigs;
		boost::split(tokens, line, boost::is_any_of("\t"));
		if (tokens.size() != 2)
			MTHROW_AND_ERR("Error RepRuleBasedOutlierCleaner::parse_rules_signals - line should contain 2 tokens with TAB. got line:\n%s\n",
				line.c_str());
		int rule_id = med_stoi(tokens[0]);
		boost::split(list_of_sigs, tokens[1], boost::is_any_of(","));
		if (rules2Signals[rule_id].size() != list_of_sigs.size())
			MTHROW_AND_ERR("Error RepRuleBasedOutlierCleaner::parse_rules_signals - rule %d contains %zu signals, got %zu signals\n",
				rule_id, rules2Signals[rule_id].size(), list_of_sigs.size());
		rules2Signals[rule_id] = list_of_sigs;
	}
	fr.close();
}

void RepRuleBasedOutlierCleaner::parse_sig_channels(const string &path) {
	ifstream fr(path);
	if (!fr.good())
		MTHROW_AND_ERR("Error RepRuleBasedOutlierCleaner::parse_sig_channels - can't read file %s\n", path.c_str());
	string line;
	while (getline(fr, line)) {
		boost::trim(line);
		if (line.empty() || line[0] == '#')
			continue;
		vector<string> tokens, list_of_sigs;
		boost::split(tokens, line, boost::is_any_of("\t"));
		if (tokens.size() != 3)
			MTHROW_AND_ERR("Error RepRuleBasedOutlierCleaner::parse_sig_channels - line should contain 3 tokens with TAB (signal name, time channel, val channel). got line:\n%s\n",
				line.c_str());
		string sigName = tokens[0];
		int time_channel = med_stoi(tokens[1]);
		int val_channel = med_stoi(tokens[1]);
		signal_channels[sigName].first = time_channel;
		signal_channels[sigName].second = val_channel;
	}
	fr.close();
}


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
		else if (field == "addRequiredSignals")addRequiredSignals = med_stoi(entry.second) != 0;
		else if (field == "rules2Signals") parse_rules_signals(entry.second); //each line is rule_id [TAB] list of signals with "," 
		else if (field == "signal_channels") parse_sig_channels(entry.second); //each line is signal_name [TAB] time_channel [TAB] val_channel 
		else if (field == "time_window") time_window = med_stoi(entry.second);
		else if (field == "nrem_attr") nRem_attr = entry.second;
		else if (field == "verbose_file") verbose_file = entry.second;
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

	if (consideredRules.empty()) {
		//init deafault to use all:
		for (const auto &rule : rules2Signals)
			consideredRules.push_back(rule.first);
	}

	for (auto& rule : rules2Signals) {
		if (std::find(consideredRules.begin(), consideredRules.end(), 0) != consideredRules.end() ||
			std::find(consideredRules.begin(), consideredRules.end(), rule.first) != consideredRules.end())
			continue;// rule remains
		else
			rules2Signals.erase(rule.first);// rule removed
	}

	// add required signals according to rules that apply to affected signals
	unordered_set<int> seen_rule;
	for (auto& rule : rules2Signals) {
		for (auto& sig : aff_signals) {
			if (std::find(rule.second.begin(), rule.second.end(), sig) != rule.second.end()) {

				if (seen_rule.find(rule.first) == seen_rule.end()) {
					rulesToApply.push_back(rule.first);
					seen_rule.insert(rule.first);
				}
				bool loopBreak = false;
				for (auto& reqSig : rule.second) {
					bool found = false;
					for (auto& existReqSig : req_signals) {
						if (reqSig == existReqSig) {
							found = true;
							break;  //already there
						}
					}

					if (!found) {
						if (addRequiredSignals)
							req_signals.insert(reqSig);    //add required signal
						else {
							rulesToApply.pop_back(); //We were asked not to load additional signals so ignore this rule
							loopBreak = true;
							break;
						}
					}
				}
				if (loopBreak)break;
			}
		}
	}
	if (!verbose_file.empty())
	{
		ofstream fw(verbose_file);
		fw.close(); //rewrite empty file
	}

	return 0;
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
	for (const auto &reqSig : req_signals)reqSignalIds.insert(dict.id(reqSig));
	for (const auto &affSig : aff_signals)affSignalIds.insert(dict.id(affSig));
	for (int affSig_id : affSignalIds)
		affected_ids_to_name[affSig_id] = dict.name(affSig_id);
	if (!verbose_file.empty() && !log_file.is_open()) {
		log_file.open(verbose_file, ios::app);
		if (!log_file.good())
			MWARN("Warnning in RepRuleBasedOutlierCleaner - verbose_file %s can't be opened\n", verbose_file.c_str());
	}
}

void RepRuleBasedOutlierCleaner::init_tables(MedDictionarySections& dict, MedSignals& sigs) {

	//rules_sids.resize(rulesToApply.size());
	//affected_by_rules.resize(rulesToApply.size());
	rules_sids.clear();
	affected_by_rules.clear();

	for (int i = 0; i < rulesToApply.size(); i++) {
		// build set of the participating signals

		for (auto& sname : rules2Signals[rulesToApply[i]]) {
			int thisSid = dict.id(sname);
			rules_sids[rulesToApply[i]].push_back(thisSid);
			affected_by_rules[rulesToApply[i]].push_back(affSignalIds.find(thisSid) != affSignalIds.end());
			if (signal_channels.find(sname) != signal_channels.end()) {
				signal_id_channels[thisSid] = signal_channels[sname];
				//check channels exists:
				if (signal_id_channels[thisSid].first >= sigs.Sid2Info.at(thisSid).n_time_channels ||
					signal_id_channels[thisSid].second >= sigs.Sid2Info.at(thisSid).n_val_channels
					)
					MTHROW_AND_ERR("Error in RepRuleBasedOutlierCleaner::init_tables - signal %s reffer to channel that not exists\n"
						"existed time_channels %d, requested %d, existed val_channels %d, request %d\n",
						sname.c_str(), sigs.Sid2Info.at(thisSid).n_time_channels, signal_id_channels[thisSid].first,
						sigs.Sid2Info.at(thisSid).n_val_channels, signal_id_channels[thisSid].second);
			}
		}
	}


}


int RepRuleBasedOutlierCleaner::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float> >& attributes_mat) {

	// get the signals
	map <int, UniversalSigVec> usvs;// from signal to its USV
	//map <int, vector <int>> removePoints; // from signal id to its remove points


	// Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	differentVersionsIterator vit(rec, reqSignalIds);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		map<int, set<int>> removePoints; // from sid to indices to be removed
		unordered_map<int, vector<int>> removePoints_Time; // from sid to Time to be removed - for printings

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
			vector<int>& mySids = rules_sids[rule];
			const vector<string> &rule_signals = rules2Signals[rule];

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
			pair<int, int> first_chan(0, 0);
			if (signal_channels.find(rule_signals.front()) != signal_channels.end())
				first_chan = signal_channels[rule_signals.front()];
			for (sPointer[0] = 0; sPointer[0] < ruleUsvs[0].len; sPointer[0]++) {
				//printf("start loop %d %d \n", sPointer[0], ruleUsvs[0].len);

				thisTime = ruleUsvs[0].Time(sPointer[0], first_chan.first);
				if (time_points.size() != 0 && thisTime > time_points[iver])break;
				bool ok = true;
				for (int i = 1; i < mySids.size(); i++) {
					pair<int, int> sig_channels(0, 0);
					if (signal_channels.find(rule_signals[i]) != signal_channels.end())
						sig_channels = signal_channels[rule_signals[i]];
					while (ruleUsvs[i].Time(sPointer[i], sig_channels.first) < thisTime - time_window && sPointer[i] < ruleUsvs[i].len - 1)
						++sPointer[i];
					//find closest (or exact):
					int try_more = 0;
					while (ruleUsvs[i].Time(sPointer[i], sig_channels.first) < thisTime && sPointer[i] + try_more < ruleUsvs[i].len - 1)
						++try_more;
					if (try_more > 0) {
						if (sPointer[i] + try_more < ruleUsvs[i].len)
							sPointer[i] += try_more; //still good pointer so use it
						else
							sPointer[i] += (try_more - 1); //still better pointer so use it's best
					}

					//printf("before ok_check: %d %d %d %d %d %d\n", i, sPointer[0], sPointer[1], sPointer[2],thisTime, ruleUsvs[i].Time(sPointer[i], time_channel));
					int time_diff = abs(ruleUsvs[i].Time(sPointer[i], sig_channels.first) - thisTime);
					//if (ruleUsvs[i].Time(sPointer[i], sig_channels.first) != thisTime) {
					if (time_diff > time_window) { //not found any candidate
						//printf("before ok_0: %d %d %d %d %d\n", rule, sPointer[0], sPointer[1], sPointer[2]);
						ok = 0;
						break;
					}
				}
				if (ok) {
					// if found all signals from same date eliminate doubles and take the last one for comparison
					vector<int> rule_val_channels;
					for (int i = 0; i < mySids.size(); i++) {
						pair<int, int> sig_channels(0, 0);
						if (signal_channels.find(rule_signals[i]) != signal_channels.end())
							sig_channels = signal_channels[rule_signals[i]];
						rule_val_channels.push_back(sig_channels.second);
					}
					for (int i = 0; i < mySids.size(); i++) {
						pair<int, int> sig_channels(0, 0);
						if (signal_channels.find(rule_signals[i]) != signal_channels.end())
							sig_channels = signal_channels[rule_signals[i]];
						int remove_same_time = 1; //try remove for signal with same time value - Can use SimValHandler before, it's better
						while (sPointer[i] - remove_same_time >= 0 &&
							ruleUsvs[i].Time(sPointer[i], sig_channels.first) == ruleUsvs[i].Time(sPointer[i] - remove_same_time, sig_channels.first)) {
							if (affected_by_rules[rulesToApply[iRule]][i]) {
								removePoints[mySids[i]].insert(sPointer[i] - remove_same_time);
								if (!verbose_file.empty())
									removePoints_Time[mySids[i]].push_back(ruleUsvs[i].Time(sPointer[i] - remove_same_time, sig_channels.first));
							}
							++remove_same_time;
						}

						// check rule and mark for removement
						//printf("before apply: %d %d %d %d\n", rule, sPointer[0],sPointer[1],sPointer[2]);
						bool ruleFlagged = applyRule(rule, ruleUsvs, rule_val_channels, sPointer);
						/*
						printf("%d R: %d P: %d t: %d   ",ruleFlagged, rule, rec.pid, thisTime);
						for (int k = 0; k < sPointer.size(); k++)printf(" %f", ruleUsvs[k].Val(sPointer[k]));
						printf("\n");
						*/
						if (ruleFlagged) {

							for (int sIndex = 0; sIndex < mySids.size(); sIndex++)
								if (affected_by_rules[rulesToApply[iRule]][sIndex]) {
									removePoints[mySids[sIndex]].insert(sPointer[sIndex]);
									if (!verbose_file.empty())
										removePoints_Time[mySids[i]].push_back(ruleUsvs[i].Time(sPointer[sIndex], sig_channels.first));
								}
						}
					}
				}
			}
		}

		// Apply removals
		size_t nRemove = 0;
		int idx = 0;
		for (auto sig : affSignalIds) {
			vector <int> toRemove(removePoints[sig].begin(), removePoints[sig].end());
			if (!verbose_file.empty()) {
				string sig_name = affected_ids_to_name[sig];
				string time_points = "";
				for (size_t i = 0; i < removePoints_Time[sig].size(); ++i)
					time_points += "," + to_string(removePoints_Time[sig][i]);
				log_file << "signal " << sig_name << " pid " << rec.pid << " removed "
					<< toRemove.size() << " int_times " << time_points << "\n";
			}
			vector <pair<int, float>>noChange;
			pair<int, int> sig_channels(0, 0);
			if (signal_id_channels.find(sig) != signal_id_channels.end())
				sig_channels = signal_id_channels[sig];
			if (rec.update(sig, iver, sig_channels.second, noChange, toRemove) < 0)
				return -1;
			if (!nRem_attr_suffix.empty() && !attributes_mat.empty()) {
				for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
					attributes_mat[pVersion][idx] = (float)toRemove.size();
				idx++;
			}
			nRemove += toRemove.size();
		}

		// Collect atttibutes
		if (!nRem_attr.empty() && !attributes_mat.empty()) {
			for (int pVersion = vit.block_first(); pVersion <= vit.block_last(); pVersion++)
				attributes_mat[pVersion][idx] = (float)nRemove;
		}
	}

	return 0;


}

bool  RepRuleBasedOutlierCleaner::applyRule(int rule, const  vector<UniversalSigVec> &ruleUsvs,
	const vector<int> &val_channels, const vector<int> &sPointer)
	// apply the rule and return true if data is consistent with the rule
	//ruleUsvs hold the signals in the order they appear in the rule in the rules2Signals above
{

	float left, right; // sides of the equality or inequality of the rule

	switch (rule) {
	case 1://BMI=Weight/Height^2*1e4
		if (ruleUsvs[2].Val(sPointer[2], val_channels[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = ruleUsvs[1].Val(sPointer[1], val_channels[1]) / ruleUsvs[2].Val(sPointer[2], val_channels[2]) / ruleUsvs[2].Val(sPointer[2], val_channels[2]) * (float)1e4;
		//printf("inputs %f %f\n", ruleUsvs[1].Val(sPointer[1]), ruleUsvs[2].Val(sPointer[2]));
		return (abs(left / right - 1) > tolerance);

	case 2://MCH=Hemoglobin/RBC*10
	case 3://MCV=Hematocrit/RBC*10
		if (ruleUsvs[2].Val(sPointer[2], val_channels[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = ruleUsvs[1].Val(sPointer[1], val_channels[1]) / ruleUsvs[2].Val(sPointer[2], val_channels[2]) * 10;
		return(abs(left / right - 1) > tolerance);

	case 4://MCHC-M=MCH/MCV*100
		if (ruleUsvs[2].Val(sPointer[2], val_channels[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = ruleUsvs[1].Val(sPointer[1], val_channels[1]) / ruleUsvs[2].Val(sPointer[2], val_channels[2]) * 100;
		return(abs(left / right - 1) > tolerance);

	case 11://HDL_over_nonHDL=HDL/NonHDLCholesterol
	case 12://HDL_over_Cholesterol=HDL/Cholesterol
		if (ruleUsvs[2].Val(sPointer[2], val_channels[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = round(ruleUsvs[1].Val(sPointer[1], val_channels[1]) / ruleUsvs[2].Val(sPointer[2], val_channels[2]) * 10) / (float)10.; //resolution in THIN is 0.1
		return(abs(left / right - 1) > tolerance);

	case 6://MPV=Platelets_Hematocrit/Platelets
		if (ruleUsvs[2].Val(sPointer[2], val_channels[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = ruleUsvs[1].Val(sPointer[1], val_channels[1]) / ruleUsvs[2].Val(sPointer[2], val_channels[2]);
		return(abs(left / right - 1) > tolerance);

	case 8://UrineAlbumin_over_Creatinine = UrineAlbumin / UrineCreatinine
		if (ruleUsvs[2].Val(sPointer[2], val_channels[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = round(ruleUsvs[1].Val(sPointer[1], val_channels[1]) / ruleUsvs[2].Val(sPointer[2], val_channels[2]) * 10) / 10;//resolution in THIN is 0.1
		return(abs(left / right - 1) > tolerance);

	case 13://HDL_over_LDL=HDL/LDL
	case 15://Cholesterol_over_HDL=Cholesterol/HDL
	case 18://LDL_over_HDL=LDL/HDL
		if (ruleUsvs[2].Val(sPointer[2], val_channels[2]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = ruleUsvs[1].Val(sPointer[1], val_channels[1]) / ruleUsvs[2].Val(sPointer[2], val_channels[2]);
		return(abs(left / right - 1) > tolerance);

	case 5://Eosinophils#+Monocytes#+Basophils#+Lymphocytes#+Neutrophils#<=WBC
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]) + ruleUsvs[1].Val(sPointer[1], val_channels[1]) + ruleUsvs[2].Val(sPointer[2], val_channels[2]) + ruleUsvs[3].Val(sPointer[3], val_channels[3]) + ruleUsvs[4].Val(sPointer[4], val_channels[4]);
		right = ruleUsvs[5].Val(sPointer[5], val_channels[5]);
		return (left*(1 - tolerance) >= right);

	case 19://Albumin<=Protein_Total	
	case 21://NRBC<=RBC
	case 22://CHADS2<=CHADS2_VASC
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = ruleUsvs[1].Val(sPointer[1], val_channels[1]);
		return(left*(1 - tolerance) >= right);

	case 7://UrineAlbumin <= UrineTotalProtein
	case 20://FreeT4<=T4
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = ruleUsvs[1].Val(sPointer[1], val_channels[1]);
		return(left*(1 - tolerance) >= right * 1000); // T4 is nmol/L free T4 is pmol/L ;  Albumin mg/L versus protein g/L

	case 9://LDL+HDL<=Cholesterol
		left = ruleUsvs[0].Val(sPointer[0]) + ruleUsvs[1].Val(sPointer[1]);
		right = ruleUsvs[2].Val(sPointer[2]);
		return (left*(1 - tolerance) > right);

	case 10://NonHDLCholesterol + HDL = Cholesterol
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]) + ruleUsvs[1].Val(sPointer[1], val_channels[1]);
		right = ruleUsvs[2].Val(sPointer[2], val_channels[2]);
		return (abs(left / right - 1) > tolerance);

	case 14://HDL_over_LDL=1/LDL_over_HDL
	case 17://Cholesterol_over_HDL = 1 / HDL_over_Cholestrol
		if (ruleUsvs[2].Val(sPointer[1], val_channels[1]) == 0)return(true);
		left = ruleUsvs[0].Val(sPointer[0], val_channels[0]);
		right = (float) 1. / ruleUsvs[1].Val(sPointer[1], val_channels[1]);
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
			attributes_mat[iver][idx++] = (float)nRemove;
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
		signalId, params.doTrim, trimMax, trimMin, params.doRemove, removeMax, removeMin, nbrsMax, nbrsMin);
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
		//! [RepCheckReq::init]
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
		//! [RepCheckReq::init]
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
				if (usvs[i].Time(j, time_channels[i]) > med_time_converter.convert_times(window_time_unit, sig_time_units[i], time_point - win_from))
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
		attributes_mat[iver][0] = (float)nMissing;
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
		else if (field == "debug") debug = stoi(entry.second) > 0;
		else if (field == "unconditional") unconditional = stoi(entry.second) > 0;
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

	if (time_channels.size() > 1 || (time_channels.size() == 1 && time_channels[0] != 0)) {
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
	//! [RepSimValHandler::get_sim_val_handle_type]
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
	else if (name == "min")
		return SIM_VAL_MIN;
	else if (name == "max")
		return SIM_VAL_MAX;
	else
		MTHROW_AND_ERR("Unkwnon sim_val_hand_type \'%s\'\n", name.c_str());
	//! [RepSimValHandler::get_sim_val_handle_type]
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
	int total_nTimes = 0;
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
		total_nTimes += nTimes;
	}

	if (debug && total_nTimes > 0 && verbose_cnt < 1) {
		MLOG("RepSimValHandler for %s - patient %d handled %d samples\n",
			signalName.c_str(), rec.pid, total_nTimes);
		++verbose_cnt;
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
			for (int iChannel = 0; iChannel < nValChannels; iChannel++)
				sums[iChannel] += usv.Val(j, iChannel);
			remove[nRemove++] = j;
		}
		pair<int, vector<float>> newChange;
		newChange.first = end;
		newChange.second.resize(nValChannels);
		for (int iChannel = 0; iChannel < nValChannels; iChannel++)
			newChange.second[iChannel] = (sums[iChannel] + usv.Val(end, iChannel)) / (end + 1 - start);
		change[nChange++] = newChange;
		nTimes++;
	}
	else if (handler_type == SIM_VAL_MIN) {
		vector<float> mins(nValChannels);
		for (int iChannel = 0; iChannel < nValChannels; iChannel++)
			mins[iChannel] = usv.Val(start, iChannel);
		remove[nRemove++] = start;

		for (int j = start + 1; j <= end; j++) {
			for (int iChannel = 0; iChannel < nValChannels; iChannel++) {
				if (usv.Val(j, iChannel) < mins[iChannel])
					mins[iChannel] = usv.Val(j, iChannel);
			}
			if (j != end) remove[nRemove++] = j;
		}

		pair<int, vector<float>> newChange;
		newChange.first = end;
		newChange.second.resize(nValChannels);
		for (int iChannel = 0; iChannel < nValChannels; iChannel++)
			newChange.second[iChannel] = mins[iChannel];
		change[nChange++] = newChange;
		nTimes++;
	}
	else if (handler_type == SIM_VAL_MAX) {
		vector<float> maxs(nValChannels);
		for (int iChannel = 0; iChannel < nValChannels; iChannel++)
			maxs[iChannel] = usv.Val(start, iChannel);
		remove[nRemove++] = start;

		for (int j = start + 1; j <= end; j++) {
			for (int iChannel = 0; iChannel < nValChannels; iChannel++) {
				if (usv.Val(j, iChannel) > maxs[iChannel])
					maxs[iChannel] = usv.Val(j, iChannel);
			}
			if (j != end) remove[nRemove++] = j;
		}
		pair<int, vector<float>> newChange;
		newChange.first = end;
		newChange.second.resize(nValChannels);
		for (int iChannel = 0; iChannel < nValChannels; iChannel++)
			newChange.second[iChannel] = maxs[iChannel];
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
		//! [RepCalcSimpleSignals::init]
		if (field == "calculator") calculator = entry.second;
		else if (field == "missing_value") missing_value = stof(entry.second);
		else if (field == "work_channel") work_channel = stoi(entry.second);
		else if (field == "max_time_search_range") max_time_search_range = stoi(entry.second);
		else if (field == "calculator_init_params") calculator_init_params = entry.second;
		else if (field == "names") boost::split(V_names, entry.second, boost::is_any_of(",:"));
		else if (field == "signals") boost::split(signals, entry.second, boost::is_any_of(",:"));
		else if (field == "signals_time_unit") signals_time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "unconditional") unconditional = stoi(entry.second) > 0;
		else if (field == "rp_type") {}
		else MTHROW_AND_ERR("Error in RepCalcSimpleSignals::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepCalcSimpleSignals::init]
	}
	if (signals_time_unit == -1 || signals_time_unit == MedTime::Undefined) {
		MWARN("Warning in RepCalcSimpleSignals::init - using signals_time_unit = Days as defualt time unit\n");
		signals_time_unit = MedTime::Days;
	}

	MLOG_D("DBG===> in RepCalcSimpleSignals init: calculator %s , time %d\n", calculator.c_str(), signals_time_unit);
	calculator_logic = SimpleCalculator::make_calculator(calculator);

	if (!calculator_init_params.empty()) {
		if (calculator_logic->init_from_string(calculator_init_params) < 0)
			return -1;
	}

	calculator_logic->missing_value = missing_value;
	calculator_logic->work_channel = work_channel;

	req_signals.clear();
	if (signals.empty() && calc2req_sigs.find(calculator) != calc2req_sigs.end())
		signals = calc2req_sigs.at(calculator);
	if (signals.empty())
		MTHROW_AND_ERR("Error in RepCalcSimpleSignals::init please provide input signals for \"%s\" calculator. no defaluts\n",
			calculator.c_str());

	for (auto & req_s : signals)
		req_signals.insert(req_s);

	// add V_names
	vector<pair<string, int>> default_virtual_signals;
	calculator_logic->list_output_signals(signals, default_virtual_signals);
	if (V_names.size() == 0) {
		//fetch from default:
		V_names.resize(default_virtual_signals.size());
		for (size_t i = 0; i < default_virtual_signals.size(); ++i)
			V_names[i] = default_virtual_signals[i].first;
	}

	// add names to required, affected and virtual_signals
	aff_signals.clear();
	virtual_signals.clear();
	for (int i = 0; i < V_names.size(); i++) {
		aff_signals.insert(V_names[i]);
		virtual_signals.push_back({ V_names[i], default_virtual_signals[i].second });
	}
	for (int i = 0; i < signals.size(); i++)
		req_signals.insert(signals[i]);

	calculator_logic->validate_arguments(signals, V_names);

	return 0;


	MERR("RepCalcSigs: ERROR: calculator %s not defined\n", calculator.c_str());
	return -1;
}

SimpleCalculator *SimpleCalculator::make_calculator(const string &calc_type) {
	//! [SimpleCalculator::make_calculator]
	if (calc_type == "ratio" || calc_type == "calc_ratio")
		return new RatioCalculator();
	else if (calc_type == "eGFR" || calc_type == "calc_eGFR")
		return new eGFRCalculator();
	else if (calc_type == "log" || calc_type == "calc_log")
		return new logCalculator();
	else if (calc_type == "sum" || calc_type == "calc_sum")
		return new SumCalculator();
	else if (calc_type == "range" || calc_type == "calc_range")
		return new RangeCalculator();
	else if (calc_type == "multiply" || calc_type == "calc_multiply")
		return new MultiplyCalculator();
	else if (calc_type == "set" || calc_type == "calc_set")
		return new SetCalculator();
	else
		HMTHROW_AND_ERR("Error: SimpleCalculator::make_calculator - unsupported calculator: %s\n",
			calc_type.c_str());
	//! [SimpleCalculator::make_calculator]
}

RepCalcSimpleSignals::~RepCalcSimpleSignals() {
	if (calculator_logic != NULL) {
		delete calculator_logic;
		calculator_logic = NULL;
	}
}

mutex RepCalcSimpleSignals_init_tables_mutex;
//.......................................................................................
void RepCalcSimpleSignals::init_tables(MedDictionarySections& dict, MedSignals& sigs)
{
	lock_guard<mutex> guard(RepCalcSimpleSignals_init_tables_mutex);
	static_input_signals.resize(signals.size());

	V_ids.clear();
	sigs_ids.clear();
	for (auto &vsig : V_names)
		V_ids.push_back(sigs.sid(vsig));
	aff_signal_ids.clear();
	aff_signal_ids.insert(V_ids.begin(), V_ids.end());

	// In the next loop it is VERY important to go over items in the ORDER they are given in calc2req
	// This is since we create a vector of sids (sigs_ids) that matches it exactly, and enables a much
	// more efficient code without going to this map for every pid. (See for example the egfr calc function)
	for (auto &rsig : signals)
		sigs_ids.push_back(sigs.sid(rsig));
	req_signal_ids.clear();
	req_signal_ids.insert(sigs_ids.begin(), sigs_ids.end());
	vector<bool> all_sigs_static(T_Last);
	all_sigs_static[T_TimeStamp] = true;
	all_sigs_static[T_Value] = true;
	all_sigs_static[T_ValShort2] = true;
	all_sigs_static[T_ValShort4] = true;
	for (size_t i = 0; i < signals.size(); ++i)
		static_input_signals[i] = all_sigs_static[sigs.Sid2Info[sigs_ids[i]].type];
	if (calculator_logic == NULL) { //recover from serialization
		calculator_logic = SimpleCalculator::make_calculator(calculator);

		if (!calculator_init_params.empty()) {
			if (calculator_logic->init_from_string(calculator_init_params) < 0)
				MTHROW_AND_ERR("Cannot init calculator from \'%s\'\n", calculator_init_params.c_str());
		}
		calculator_logic->missing_value = missing_value;
	}
	calculator_logic->init_tables(dict, sigs, signals);
}

//.......................................................................................

bool is_in_time_range(vector<UniversalSigVec> &usvs, vector<int> idx, int active_id,
	int time_range, int time_unit, int &sum_diff) {
	int time = usvs[active_id].TimeU(idx[active_id] - 1, time_unit);
	sum_diff = 0; //if not found
	for (size_t i = 0; i < idx.size(); ++i)
	{
		if (i == active_id)
			continue;//skip current
		if (idx[i] == 0) //one signal is not yet started - will be happen in future, so waiting for it!
			return false;

		int ref_time = usvs[i].TimeU(idx[i] - 1, time_unit);
		if (time - ref_time > time_range)
			return false;
		sum_diff += time - ref_time;
	}
	return true;
}

bool no_missings(const vector<float> &vals, float missing_value) {
	for (size_t i = 0; i < vals.size(); ++i)
		if (vals[i] == missing_value)
			return false;
	return true;
}
int RepCalcSimpleSignals::apply_calc_in_time(PidDynamicRec& rec, vector<int>& time_points) {
	// Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("RepCalcSimpleSignals::apply_calc_in_time nversions mismatch\n");
		return -1;
	}
	int v_out_sid = V_ids[0];
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepCalcSimpleSignals::apply_calc_in_time - V_ids is not initialized - bad call\n");
	int n_vals = work_channel + 1;
	//first lets fetch "static" signals without Time field:

	//MLOG("DBG3===>: apply_calc_in_time: pid %d\n", rec.pid);

	set<int> iteratorSignalIds;
	vector<int> timed_sigs;
	for (size_t i = 0; i < sigs_ids.size(); ++i)
		if (!static_input_signals[i]) {
			iteratorSignalIds.insert(sigs_ids[i]);
			timed_sigs.push_back(sigs_ids[i]);
		}

	allVersionsIterator vit(rec, iteratorSignalIds);
	rec.usvs.resize(timed_sigs.size());

	int first_ver = vit.init();
	vector<float> static_signals_values(sigs_ids.size(), missing_value);
	for (size_t i = 0; i < static_signals_values.size(); ++i)
		if (static_input_signals[i]) {
			UniversalSigVec usv;
			rec.uget(sigs_ids[i], first_ver, usv);
			static_signals_values[i] = usv.Val(0);
		}

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		for (size_t i = 0; i < timed_sigs.size(); ++i)
			rec.uget(timed_sigs[i], iver, rec.usvs[i]);
		bool all_non_empty = true;
		for (size_t i = 0; i < rec.usvs.size() && all_non_empty; ++i)
			all_non_empty = rec.usvs[i].len > 0;
		int last_time = -1;
		if (all_non_empty) {
			vector<int> idx(timed_sigs.size());
			int active_id = medial::repository::fetch_next_date(rec.usvs, idx);
			int final_size = 0;
			vector<float> v_vals;
			vector<int> v_times;
			int max_diff = -1;
			while (active_id >= 0) {
				//iterate on time ordered of signals - Let's try to calc signal:
				bool can_calc = is_in_time_range(rec.usvs, idx, active_id, max_time_search_range, signals_time_unit, max_diff);
				if (can_calc) {
					vector<float> collected_vals(sigs_ids.size());
					int time_idx = 0;
					for (size_t i = 0; i < sigs_ids.size(); ++i) {
						if (static_input_signals[i])
							collected_vals[i] = static_signals_values[i];
						else {
							collected_vals[i] = rec.usvs[time_idx].Val(idx[time_idx] - 1, work_channel);
							++time_idx;
						}
					}
					if (no_missings(collected_vals, missing_value)) {
						float prev_val = missing_value;
						if (last_time == rec.usvs[active_id].Time(idx[active_id] - 1)) {
							--final_size; //override last value
							prev_val = v_vals[final_size];
						}
						if (v_times.size() < final_size + 1) {
							v_times.resize(final_size + 1);
							v_vals.resize(n_vals * final_size + n_vals);
						}
						v_times[final_size] = rec.usvs[active_id].Time(idx[active_id] - 1);
						v_vals[n_vals * final_size + n_vals - 1] = calculator_logic->do_calc(collected_vals);
						for (int kk = 0; kk < n_vals - 1; ++kk)
							v_vals[n_vals * final_size + kk] = rec.usvs[active_id].Val(idx[active_id] - 1, kk);

						if (v_vals[n_vals * final_size + n_vals - 1] != missing_value) { //insert only legal values (missing_value when ilegal)!
							++final_size;
							last_time = rec.usvs[active_id].Time(idx[active_id] - 1);
						}
						else if (last_time == rec.usvs[active_id].Time(idx[active_id] - 1)) {
							v_vals[n_vals * final_size + n_vals - 1] = prev_val; //return previous val that was not missing
							//Pay attention it still update rest value channels
							++final_size;
						}
					}
				}

				active_id = medial::repository::fetch_next_date(rec.usvs, idx);
			}
			// pushing virtual data into rec (into orig version)
			rec.set_version_universal_data(v_out_sid, iver, &v_times[0], &v_vals[0], final_size);
		}
	}

	return 0;
}

void RepCalcSimpleSignals::print() {
	MLOG("RepCalcSimpleSignals: calculator: %s : calculator_init_params %s : max_time_search_range %d signals_time_unit %s signals: %s, V_names: %s, req_signals: %s, aff_signals: %s, work_channel: %d\n",
		calculator.c_str(), calculator_init_params.c_str(), max_time_search_range, med_time_converter.type_to_string(signals_time_unit).c_str(),
		medial::io::get_list(signals).c_str(), medial::io::get_list(V_names).c_str(),
		medial::io::get_list(req_signals).c_str(), medial::io::get_list(aff_signals).c_str(),
		work_channel);
}

//.......................................................................................
int RepCalcSimpleSignals::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat)
{
	//handle special calculations
	apply_calc_in_time(rec, time_points);


	return 0;
}

int RepCombineSignals::init(map<string, string> &mapper) {
	vector<string> tokens;
	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepCombineSignals::init]
		if (field == "names") output_name = entry.second;
		else if (field == "signals") signals = boost::split(signals, entry.second, boost::is_any_of(","));
		else if (field == "unconditional") unconditional = stoi(entry.second) > 0;
		else if (field == "rp_type") {}
		else if (field == "factors") {
			boost::split(tokens, entry.second, boost::is_any_of(","));
			factors.resize(tokens.size());
			for (size_t i = 0; i < tokens.size(); ++i)
				factors[i] = stof(tokens[i]);
		}
		else MTHROW_AND_ERR("Error in RepCombineSignals::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepCombineSignals::init]
	}
	if (signals.empty())
		MTHROW_AND_ERR("Error in RepCombineSignals::init - parameter \"signals\" should be passed.\n");
	factors.resize(signals.size(), 1);
	if (output_name.empty()) {
		output_name = "COMBO_" + signals[0];
		for (size_t i = 1; i < signals.size(); ++i)
			output_name += "_" + signals[i];
	}

	aff_signals.clear();
	aff_signals.insert(output_name);
	req_signals.clear();
	req_signals.insert(signals.begin(), signals.end());
	virtual_signals.clear();
	virtual_signals.push_back(pair<string, int>(output_name, T_DateRangeVal2));

	return 0;
}

int RepCombineSignals::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {
	//uses each time point - If have only drug amount  (2nd signal) so using second signal value
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepCombineSignals::_apply - v_out_sid is not initialized - bad call\n");
	//first lets fetch "static" signals without Time field:

	set<int> set_ids(sigs_ids.begin(), sigs_ids.end());
	allVersionsIterator vit(rec, set_ids);
	rec.usvs.resize(sigs_ids.size());

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		for (size_t i = 0; i < sigs_ids.size(); ++i)
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);


		vector<int> idx(sigs_ids.size());
		int active_id = medial::repository::fetch_next_date(rec.usvs, idx);
		int final_size = 0;
		vector<float> v_vals;
		vector<int> v_times;
		int last_time = -1;
		int last_val_ch1 = -1;
		while (active_id >= 0) {
			if (last_time == rec.usvs[active_id].Time(idx[active_id] - 1) &&
				last_val_ch1 == rec.usvs[active_id].Val(idx[active_id] - 1)) {
				active_id = medial::repository::fetch_next_date(rec.usvs, idx);
				continue; //skip same time, and same first channel value
			}

			if (v_times.size() < 2 * final_size + 1) {
				v_times.resize(2 * final_size + 2);
				v_vals.resize(2 * final_size + 2);
			}
			v_times[2 * final_size] = rec.usvs[active_id].Time(idx[active_id] - 1);
			v_vals[2 * final_size] = rec.usvs[active_id].Val(idx[active_id] - 1);
			v_vals[2 * final_size + 1] = factors[active_id] * rec.usvs[active_id].Val(idx[active_id] - 1, 1);
			++final_size;

			last_time = rec.usvs[active_id].Time(idx[active_id] - 1);
			last_val_ch1 = rec.usvs[active_id].Val(idx[active_id] - 1);
			active_id = medial::repository::fetch_next_date(rec.usvs, idx);
		}
		// pushing virtual data into rec (into orig version)
		rec.set_version_universal_data(v_out_sid, iver, &v_times[0], &v_vals[0], final_size);

	}

	return 0;
}

void RepCombineSignals::init_tables(MedDictionarySections& dict, MedSignals& sigs) {
	v_out_sid = sigs.sid(output_name);
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepCombineSignals::init_tables - virtual output signal %s not found\n",
			output_name.c_str());
	aff_signal_ids.clear();
	aff_signal_ids.insert(v_out_sid);
	sigs_ids.resize(signals.size());
	for (size_t i = 0; i < signals.size(); ++i)
	{
		sigs_ids[i] = sigs.sid(signals[i]);
		if (sigs_ids[i] < 0)
			MTHROW_AND_ERR("Error in RepCombineSignals::init_tables - input signal %s not found\n",
				signals[i].c_str());
		SignalInfo &si = sigs.Sid2Info[sigs_ids[i]];

		if (si.n_val_channels < 2)
			MTHROW_AND_ERR("ERROR in RepCombineSignals::init_tables - input signal %s should contain 2 val channels\n",
				signals[i].c_str());
	}
	req_signal_ids.clear();
	req_signal_ids.insert(sigs_ids.begin(), sigs_ids.end());
}

void RepCombineSignals::register_virtual_section_name_id(MedDictionarySections& dict) {
	dict.SectionName2Id[output_name] = dict.section_id(signals.front());
}

void RepCombineSignals::print() {
	MLOG("RepCombineSignals: output_name: %s : signals %s : req_signals %s aff_signals %s\n",
		output_name.c_str(), medial::io::get_list(signals).c_str(), medial::io::get_list(req_signals).c_str(), medial::io::get_list(aff_signals).c_str());
}

int RepSignalRate::init(map<string, string> &mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepSignalRate::init]
		if (field == "names") output_name = entry.second;
		else if (field == "input_name") input_name = entry.second;
		else if (field == "unconditional") unconditional = stoi(entry.second) > 0;
		else if (field == "work_channel") work_channel = stoi(entry.second);
		else if (field == "factor") factor = stof(entry.second);
		else if (field == "rp_type") {}
		else MTHROW_AND_ERR("Error in RepSignalRate::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepSignalRate::init]
	}
	if (input_name.empty())
		MTHROW_AND_ERR("Error in RepSignalRate::init - input signal should be passed\n");
	if (work_channel > 1)
		MTHROW_AND_ERR("Error in RepSignalRate::init - unsupported work_channel > 1\n");
	aff_signals.clear();
	aff_signals.insert(output_name);
	req_signals.clear();
	req_signals.insert(input_name);
	virtual_signals.clear();
	if (work_channel == 0)
		virtual_signals.push_back(pair<string, int>(output_name, T_DateRangeVal));
	else
		virtual_signals.push_back(pair<string, int>(output_name, T_DateRangeVal2));

	return 0;
}

void RepSignalRate::add_virtual_signals(map<string, int> &_virtual_signals) {
	_virtual_signals[output_name] = virtual_signals.front().second;
}

void RepSignalRate::init_tables(MedDictionarySections& dict, MedSignals& sigs) {
	v_out_sid = sigs.sid(output_name);
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepSignalRate::init_tables - virtual output signal %s not found\n",
			output_name.c_str());
	aff_signal_ids.clear();
	aff_signal_ids.insert(v_out_sid);
	in_sid = sigs.sid(input_name);
	if (in_sid < 0)
		MTHROW_AND_ERR("Error in RepSignalRate::init_tables - input signal %s not found\n",
			input_name.c_str());
	req_signal_ids.clear();
	req_signal_ids.insert(in_sid);

	SignalInfo &si = sigs.Sid2Info[in_sid];

	if (si.n_time_channels < 2)
		MTHROW_AND_ERR("ERROR in RepSignalRate::init_tables - input signal %s should contain 2 time channels\n",
			input_name.c_str());
	if (si.n_val_channels < work_channel + 1)
		MTHROW_AND_ERR("ERROR in RepSignalRate::init_tables - input signal %s should contain %d val channels\n",
			input_name.c_str(), work_channel + 1);
}

void RepSignalRate::register_virtual_section_name_id(MedDictionarySections& dict) {
	dict.SectionName2Id[output_name] = dict.section_id(input_name);
}

int RepSignalRate::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {
	//uses each time point - I have only signal value need to tranform into signal_rate divide by time unit
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepCombineSignals::_apply - v_out_sid is not initialized - bad call\n");
	//first lets fetch "static" signals without Time field:

	set<int> set_ids;
	set_ids.insert(in_sid);
	allVersionsIterator vit(rec, set_ids);
	rec.usvs.resize(1);

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		rec.uget(in_sid, iver, rec.usvs[0]);

		vector<float> v_vals;
		vector<int> v_times;
		for (int i = 0; i < rec.usvs[0].len; ++i) {
			int start_time = rec.usvs[0].Time(i);
			int end_time = rec.usvs[0].Time(i, 1);
			int diff_time = med_time_converter.diff_times(end_time, start_time, rec.usvs[0].time_unit(),
				global_default_time_unit);
			if (diff_time == 0 || start_time == 0 || end_time == 0)
				continue;

			v_times.push_back(start_time);
			v_times.push_back(end_time);

			//add previous channels
			for (int k = 0; k < work_channel; ++k)
				v_vals.push_back(rec.usvs[0].Val(i, k));
			//update current channel to be rate
			v_vals.push_back(factor * rec.usvs[0].Val(i, work_channel) / diff_time);
		}
		// pushing virtual data into rec (into orig version)
		if (rec.usvs[0].len > 0)
			rec.set_version_universal_data(v_out_sid, iver, &v_times[0], &v_vals[0], (int)v_vals.size() / (work_channel + 1));

	}

	return 0;
}

void RepSignalRate::print() {
	MLOG("RepSignalRate: input_name: %s, output_name: %s : factor: %f, work_channel: %d, req_signals %s aff_signals %s\n",
		input_name.c_str(), output_name.c_str(), factor, work_channel, medial::io::get_list(req_signals).c_str(), medial::io::get_list(aff_signals).c_str());
}

int RepSplitSignal::init(map<string, string>& mapper) {
	vector<string> tokens;
	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepSplitSignal::init]
		if (field == "input_name") input_name = entry.second;
		else if (field == "names") boost::split(names, entry.second, boost::is_any_of(","));
		else if (field == "sets") boost::split(sets, entry.second, boost::is_any_of(","));
		else if (field == "factors") {
			boost::split(tokens, entry.second, boost::is_any_of(","));
			factors.resize(tokens.size());
			for (size_t i = 0; i < tokens.size(); ++i)
				factors[i] = stof(tokens[i]);
		}
		else if (field == "unconditional") unconditional = stoi(entry.second) > 0;
		else if (field == "rp_type") {}
		else MTHROW_AND_ERR("Error in RepSplitSignal::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepSplitSignal::init]
	}
	if (input_name.empty())
		MTHROW_AND_ERR("ERROR in RepSplitSignal::init - must provide input_name\n");
	if (names.size() != 2)
		MTHROW_AND_ERR("ERROR in RepSplitSignal::init - must provide only 2 names for output\n");
	if (sets.empty())
		MTHROW_AND_ERR("ERROR in RepSplitSignal::init - must provide sets\n");
	factors.resize(names.size(), 1);

	aff_signals.clear();
	aff_signals.insert(names.begin(), names.end());
	req_signals.clear();
	req_signals.insert(input_name);
	virtual_signals.clear();
	for (size_t i = 0; i < names.size(); ++i)
		virtual_signals.push_back(pair<string, int>(names[i], T_DateRangeVal2));

	return 0;
}

void RepSplitSignal::init_tables(MedDictionarySections& dict, MedSignals& sigs) {
	//init Flags:
	int section_id = dict.section_id(input_name);
	dict.prep_sets_lookup_table(section_id, sets, Flags);

	in_sid = sigs.sid(input_name);
	if (in_sid < 0)
		MTHROW_AND_ERR("Error in RepSplitSignal::init_tables - input signal %s not found\n",
			input_name.c_str());
	req_signal_ids.clear();
	req_signal_ids.insert(in_sid);

	V_ids.resize(names.size());
	for (size_t i = 0; i < V_ids.size(); ++i)
	{
		V_ids[i] = sigs.sid(names[i]);
		if (V_ids[i] < 0)
			MTHROW_AND_ERR("Error in RepSplitSignal::init_tables - virtual output signal %s not found\n",
				names[i].c_str());
	}
	aff_signal_ids.clear();
	aff_signal_ids.insert(V_ids.begin(), V_ids.end());

	SignalInfo &si = sigs.Sid2Info[in_sid];

	if (si.n_val_channels < 2)
		MTHROW_AND_ERR("ERROR in RepSplitSignal::init_tables - input signal %s should contain 2 val channels\n",
			input_name.c_str());
}

void RepSplitSignal::add_virtual_signals(map<string, int> &_virtual_signals) {
	for (size_t i = 0; i < names.size(); ++i)
		_virtual_signals[names[i]] = virtual_signals[i].second;
}

void RepSplitSignal::register_virtual_section_name_id(MedDictionarySections& dict) {
	for (size_t i = 0; i < names.size(); ++i)
		dict.SectionName2Id[names[i]] = dict.section_id(input_name);
}

int RepSplitSignal::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}
	for (size_t i = 0; i < V_ids.size(); ++i)
		if (V_ids[i] < 0)
			MTHROW_AND_ERR("Error in RepSplitSignal::_apply - V_ids is not initialized - bad call\n");


	set<int> set_ids;
	set_ids.insert(in_sid);
	allVersionsIterator vit(rec, set_ids);
	rec.usvs.resize(1);

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		rec.uget(in_sid, iver, rec.usvs[0]);

		vector<vector<int>> v_times(names.size());
		vector<vector<float>> v_vals(names.size());
		for (int i = 0; i < rec.usvs[0].len; ++i) {
			int time = rec.usvs[0].Time(i);
			float orig_val = rec.usvs[0].Val(i);
			int idx = Flags[orig_val];

			v_times[idx].push_back(time);
			v_times[idx].push_back(0);
			v_vals[idx].push_back(orig_val); //chan 0
			v_vals[idx].push_back(rec.usvs[0].Val(i, 1) * factors[idx]); //chan 1
		}
		// pushing virtual data into rec (into orig version)
		for (size_t i = 0; i < names.size(); ++i)
			if (!v_times[i].empty())
				rec.set_version_universal_data(V_ids[i], iver, &v_times[i][0], &v_vals[i][0], (int)v_times[i].size() / 2);
	}

	return 0;
}

void RepSplitSignal::print() {
	MLOG("RepSplitSignal: input_name: %s, names: %s, req_signals %s aff_signals %s\n",
		input_name.c_str(), medial::io::get_list(names).c_str(), medial::io::get_list(req_signals).c_str(), medial::io::get_list(aff_signals).c_str());
}

//=======================================================================================
// RepAggregationPeriod
//=======================================================================================

int RepAggregationPeriod::init(map<string, string>& mapper) {
	vector<string> tokens;
	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepAggregationPeriod::init]
		if (field == "input_name") input_name = entry.second;
		else if (field == "output_name") output_name = entry.second;
		else if (field == "sets") boost::split(sets, entry.second, boost::is_any_of(","));
		else if (field == "rp_type") {}
		else if (field == "period") period = med_stoi(entry.second);
		else if (field == "time_unit_sig") time_unit_sig = med_time_converter.string_to_type(entry.second);
		else if (field == "time_unit_win") time_unit_win = med_time_converter.string_to_type(entry.second);
		else MTHROW_AND_ERR("Error in RepAggregationPeriod::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepAggregationPeriod::init]
	}
	if (input_name.empty())
		MTHROW_AND_ERR("ERROR in RepAggregationPeriod::init - must provide input_name\n");
	if (output_name.empty())
		MTHROW_AND_ERR("ERROR in RepAggregationPeriod::init - must provide output_name\n");
	if (sets.empty())
		MTHROW_AND_ERR("ERROR in RepAggregationPeriod::init - must provide sets\n");
	if (period == 0)
		MLOG("WARNING in RepAggregationPeriod::init  - period set to default value: %d\n", period);

	aff_signals.insert(output_name);
	req_signals.insert(input_name);
	virtual_signals.push_back(pair<string, int>(output_name, T_TimeRange));

	return 0;
}

void RepAggregationPeriod::init_tables(MedDictionarySections& dict, MedSignals& sigs) {
	//init:
	int section_id = dict.section_id(input_name);
	dict.prep_sets_lookup_table(section_id, sets, lut);

	in_sid = sigs.sid(input_name);
	if (in_sid < 0)
		MTHROW_AND_ERR("Error in RepAggregationPeriod::init_tables - input signal %s not found\n",
			input_name.c_str());
	req_signal_ids.insert(in_sid);

	V_ids.resize(1);
	V_ids[0] = sigs.sid(output_name);
	if (V_ids[0] < 0)
		MTHROW_AND_ERR("Error in RepAggregationPeriod::init_tables - virtual output signal %s not found\n",
			output_name.c_str());

	aff_signal_ids.insert(V_ids.begin(), V_ids.end());
}

void RepAggregationPeriod::add_virtual_signals(map<string, int> &_virtual_signals) {
	_virtual_signals[output_name] = virtual_signals[0].second;
}

int RepAggregationPeriod::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}
	for (size_t i = 0; i < V_ids.size(); ++i)
		if (V_ids[i] < 0)
			MTHROW_AND_ERR("Error in RepAggregationPeriod::_apply - V_ids is not initialized - bad call\n");

	set<int> set_ids;
	set_ids.insert(in_sid);
	allVersionsIterator vit(rec, set_ids);
	rec.usvs.resize(1);
	int sig_period = med_time_converter.convert_times(time_unit_win, time_unit_sig, period);

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		rec.uget(in_sid, iver, rec.usvs[0]);

		vector<int> v_times;
		vector<float> v_vals;
		if (rec.usvs[0].len < 1) { //in case this version of the signal is empty
			continue;
		}
		int start_time = 0, end_time = 0;
		bool first = true;
		for (int i = 0; i < rec.usvs[0].len; ++i) { // find remaining valid values
			if (lut[rec.usvs[0].Val(i)] == 0) { // value not in set
				continue;
			}
			int time = rec.usvs[0].Time(i);
			if (first) {
				start_time = time;
				end_time = start_time + sig_period;
				first = false;
			}
			else {
				if (med_time_converter.diff_times(time, end_time, time_unit_sig, time_unit_sig) <= 0) {
					end_time = max(end_time, time + sig_period);
				}
				else { // found a signal that is not included in the current period, close old period and open new one
					v_times.push_back(start_time);
					v_times.push_back(end_time);

					start_time = time;
					end_time = time + sig_period;
				}
			}


		}
		if (!first) { // else - no valid set values were found
			v_times.push_back(start_time);
			v_times.push_back(end_time);
			// pushing virtual data into rec (into orig version)
			rec.set_version_universal_data(V_ids[0], iver, &v_times[0], &v_vals[0], (int)v_times.size() / 2);
		}
	}

	return 0;
}

void RepAggregationPeriod::print() {
	MLOG("RepAggregationPeriod: input_name: %s, output_name: %s, req_signals %s aff_signals %s\n",
		input_name.c_str(), output_name.c_str(), medial::io::get_list(req_signals).c_str(), medial::io::get_list(aff_signals).c_str());
}


//=======================================================================================
// BasicRangeCleaner
//=======================================================================================

int RepBasicRangeCleaner::init(map<string, string>& mapper)
{
	MLOG("In RepBasicRangeCleaner init\n");
	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepBasicRangeCleaner::init]
		if (field == "signal_name") { signal_name = entry.second; }
		else if (field == "rp_type") {}
		else if (field == "ranges_sig_name") { ranges_name = entry.second; }
		else if (field == "output_name") { output_name = entry.second; }
		else if (field == "time_channel") time_channel = med_stoi(entry.second);
		else if (field == "output_type") output_type = med_stoi(entry.second); // needs to match the input signal type! defaults to range-value signal (3)
		else MTHROW_AND_ERR("Error in RepBasicRangeCleaner::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepBasicRangeCleaner::init]
	}
	if (signal_name.empty())
		MTHROW_AND_ERR("ERROR in RepBasicRangeCleaner::init - must provide signal_name\n");
	if (ranges_name.empty())
		MTHROW_AND_ERR("ERROR in RepBasicRangeCleaner::init - must provide ranges_sig_name\n");
	if (output_name.empty()) {
		output_name = signal_name + "_" + ranges_name;
		MLOG("WARNING in RepBasicRangeCleaner::init - no output_name provided, using input signal combination: %s", output_name.c_str());
	}


	req_signals.insert(signal_name);
	req_signals.insert(ranges_name);
	aff_signals.insert(output_name);

	virtual_signals.push_back(pair<string, int>(output_name, output_type));
	return 0;
}

void RepBasicRangeCleaner::init_tables(MedDictionarySections& dict, MedSignals& sigs) {
	signal_id = sigs.sid(signal_name);
	ranges_id = sigs.sid(ranges_name);
	output_id = sigs.sid(output_name);
	req_signal_ids.insert(signal_id);
	req_signal_ids.insert(ranges_id);
	aff_signal_ids.insert(output_id);
}

void RepBasicRangeCleaner::add_virtual_signals(map<string, int> &_virtual_signals) {
	_virtual_signals[output_name] = virtual_signals[0].second;
}

int  RepBasicRangeCleaner::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float> >& attributes_mat) {

	if (signal_id == -1) {
		MERR("Uninitialized signal_id\n");
		return -1;
	}
	if (ranges_id == -1) {
		MERR("Uninitialized ranges_id\n");
		return -1;
	}
	if (output_id == -1) {
		MERR("Uninitialized output_id\n");
		return -1;
	}
	// Check that we have the correct number of dynamic-versions : one per time-point (if given)
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}
	int len;
	set<int> set_ids;
	set_ids.insert(signal_id);
	set_ids.insert(ranges_id);
	set_ids.insert(output_id);
	allVersionsIterator vit(rec, set_ids);
	rec.usvs.resize(3);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		// setup
		rec.uget(signal_id, iver, rec.usvs[0]); // original signal
		rec.uget(ranges_id, iver, rec.usvs[1]); // range signal
		rec.uget(output_id, iver, rec.usvs[2]); // output - virtual signal
		int time_channels = rec.usvs[0].n_time_channels();
		int val_channels = rec.usvs[0].n_val_channels();
		len = rec.usvs[0].len;
		vector<int> v_times(len * time_channels); // initialize size to avoid multiple resizings for long signals
		vector<float> v_vals(len * val_channels);

		// Collect elements to keep
		int nKeep = 0;
		int j = 0;
		for (int i = 0; i < len; i++) { //iterate over input signal
			int time = rec.usvs[0].Time(i, time_channel);
			// remove element only if it doesn't appear in any range
			bool doRemove = true;
			//for (int j = 0; j < rec.usvs[1].len; j++) { // slower version
			//	if (time >= rec.usvs[1].Time(j, 0) && time <= rec.usvs[1].Time(j, 1)) {
			//		doRemove = false;
			//		break;
			//	}
			//}
			for (; j < rec.usvs[1].len; j++) { // iterate over range signal
				if (time > rec.usvs[1].Time(j, 1)) continue;
				if (time >= rec.usvs[1].Time(j, 0) && time <= rec.usvs[1].Time(j, 1)) doRemove = false;
				break;
			}
			if (!doRemove) {
				for (int t = 0; t < time_channels; t++) v_times[nKeep * time_channels + t] = rec.usvs[0].Time(i, t);
				for (int v = 0; v < val_channels; v++) v_vals[nKeep * val_channels + v] = rec.usvs[0].Val(i, v);
				nKeep++;
			}
		}
		// v_times and v_vals are likely longer than necessary, it's ok because nKeep defines which part of the vector is used.
		rec.set_version_universal_data(output_id, iver, &v_times[0], &v_vals[0], nKeep);
	}

	return 0;
}

void RepBasicRangeCleaner::print()
{
	MLOG("RepBasicRangeCleaner: signal: %d %s : t_channel %d : ranges_signal: %d %s : output_signal: %d %s\n",
		signal_id, signal_name.c_str(), time_channel, ranges_id, ranges_name.c_str(), output_id, output_name.c_str());
}

int RepAggregateSignal::init(map<string, string> &mapper) {
	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepAggregateSignal::init]
		if (field == "output_name") output_name = entry.second;
		else if (field == "signal") signalName = entry.second;
		else if (field == "unconditional") unconditional = med_stoi(entry.second) > 0;
		else if (field == "work_channel") work_channel = med_stoi(entry.second);
		else if (field == "start_time_channel") start_time_channel = med_stoi(entry.second);
		else if (field == "end_time_channel") end_time_channel = med_stoi(entry.second);
		else if (field == "time_window") time_window = med_stoi(entry.second);
		else if (field == "time_unit") time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "factor") factor = med_stof(entry.second);
		else if (field == "drop_missing_rate") drop_missing_rate = med_stof(entry.second);
		else if (field == "buffer_first") buffer_first = med_stoi(entry.second) > 0;
		else if (field == "rp_type") {}
		else MTHROW_AND_ERR("Error in RepAggregateSignal::init - Unsupported param \"%s\"\n", field.c_str());
		//! [RepAggregateSignal::init]
	}
	if (signalName.empty())
		MTHROW_AND_ERR("Error in RepAggregateSignal::init - signal should be passed\n");
	if (time_window == 0)
		MTHROW_AND_ERR("Error in RepAggregateSignal::init - time_window should be passed\n");
	if (output_name.empty())
		MTHROW_AND_ERR("Error in RepAggregateSignal::init - output_name should be passed\n");

	aff_signals.clear();
	aff_signals.insert(output_name);
	req_signals.clear();
	req_signals.insert(signalName);
	virtual_signals.clear();
	virtual_signals.push_back(pair<string, int>(output_name, T_DateVal));

	return 0;
}

void RepAggregateSignal::add_virtual_signals(map<string, int> &_virtual_signals) {
	_virtual_signals[output_name] = virtual_signals.front().second;
}

void RepAggregateSignal::init_tables(MedDictionarySections& dict, MedSignals& sigs) {
	v_out_sid = sigs.sid(output_name);
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepAggregateSignal::init_tables - virtual output signal %s not found\n",
			output_name.c_str());
	aff_signal_ids.clear();
	aff_signal_ids.insert(v_out_sid);
	in_sid = sigs.sid(signalName);
	if (in_sid < 0)
		MTHROW_AND_ERR("Error in RepAggregateSignal::init_tables - input signal %s not found\n",
			signalName.c_str());
	req_signal_ids.clear();
	req_signal_ids.insert(in_sid);

	SignalInfo &si = sigs.Sid2Info[in_sid];

	if (si.n_time_channels <= max(end_time_channel, start_time_channel))
		MTHROW_AND_ERR("ERROR in RepAggregateSignal::init_tables - input signal %s should contain [%d, %d] time channels\n",
			signalName.c_str(), start_time_channel, end_time_channel);
	if (si.n_val_channels < work_channel + 1)
		MTHROW_AND_ERR("ERROR in RepAggregateSignal::init_tables - input signal %s should contain %d val channels\n",
			signalName.c_str(), work_channel + 1);
}

void RepAggregateSignal::register_virtual_section_name_id(MedDictionarySections& dict) {
	dict.SectionName2Id[output_name] = dict.section_id(signalName);
}

void update_collected(vector<float> &collected, vector<int> collected_times[], int start_time, int end_time) {
	//iterate throght collected and remove indexes with no intersect with start_time->end_time
	vector<float> sel_vals;
	vector<int> sel_times[2];
	for (int i = 0; i < collected.size(); ++i)
	{
		int start = collected_times[0][i];
		int end = collected_times[1][i];
		if (start <= end_time && end > start_time) {
			sel_vals.push_back(collected[i]);
			sel_times[0].push_back(start);
			sel_times[1].push_back(end);
		}
	}
	//comit selection:
	collected.swap(sel_vals);
	collected_times[0].swap(sel_times[0]);
	collected_times[1].swap(sel_times[1]);
}

float calc_value(const vector<int> collected_times[], const vector<float> &collected,
	int start_time, int end_time, float threshold) {
	float res = 0;

	int window_len = end_time - start_time;
	//asuume sorted by collected_times[0] which is start_time
	int coverage = 0;
	int prev_end = 0;
	for (int i = 0; i < collected.size(); ++i)
	{
		int start = collected_times[0][i];
		int end = collected_times[1][i];
		int curr_len = end - start;
		float v = collected[i];

		int real_end = end_time;
		int real_start = start_time;
		if (end < end_time)
			real_end = end;
		if (start > start_time)
			real_start = start;
		int peroid_len = real_end - real_start;
		float ratio_weight = 1;
		if (curr_len > 0)
			ratio_weight = peroid_len / float(curr_len);

		res += ratio_weight * v;
		if (i == 0)
			coverage += peroid_len;
		else {
			if (real_start >= prev_end)
				coverage += peroid_len;
			else if (real_end >= prev_end)
				coverage += real_end - prev_end;
		}
		prev_end = real_end;
	}
	float missing_rate = 0;
	if (window_len > 0)
		missing_rate = 1 - coverage / float(window_len);
	if (missing_rate > threshold)
		return MED_MAT_MISSING_VALUE;
	return res;
}

int RepAggregateSignal::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {

	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepAggregateSignal::_apply - v_out_sid is not initialized - bad call\n");
	//first lets fetch "static" signals without Time field:

	set<int> set_ids;
	set_ids.insert(in_sid);
	allVersionsIterator vit(rec, set_ids);
	rec.usvs.resize(1);

	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		rec.uget(in_sid, iver, rec.usvs[0]);

		vector<float> v_vals;
		vector<int> v_times;
		vector<float> collected;
		vector<int> collected_times[2];
		int first_time = 0;
		for (int i = 0; i < rec.usvs[0].len; ++i) {
			int end_time = rec.usvs[0].Time(i, end_time_channel);
			int start_time = rec.usvs[0].Time(i, start_time_channel);
			if (start_time <= 0 || end_time <= 0)
				continue;
			if (first_time == 0)
				first_time = med_time_converter.convert_times(global_default_time_unit, time_unit, rec.usvs[0].Time(i, start_time_channel));

			int end_win_time = med_time_converter.convert_times(global_default_time_unit, time_unit, end_time);
			int start_win_time = end_win_time - time_window;

			collected.push_back(rec.usvs[0].Val(i, work_channel));
			collected_times[0].push_back(med_time_converter.convert_times(global_default_time_unit, time_unit, start_time));
			collected_times[1].push_back(end_win_time);

			update_collected(collected, collected_times, start_win_time, end_win_time);
			if (buffer_first && collected_times[0].back() - first_time < time_window)
				continue; //do not add - wait till buffer filled
			//get value:
			float val = calc_value(collected_times, collected, start_win_time, end_win_time, drop_missing_rate);
			if (val != MED_MAT_MISSING_VALUE) {
				val *= factor;
				v_times.push_back(end_time);
				v_vals.push_back(val);
			}
		}
		// pushing virtual data into rec (into orig version)
		if (rec.usvs[0].len > 0)
			rec.set_version_universal_data(v_out_sid, iver, &v_times[0], &v_vals[0], (int)v_vals.size());

	}

	return 0;
}

void RepAggregateSignal::print() {
	MLOG("RepAggregateSignal:: signal:%s, output_name:%s, work_channel=%d, factor=%2.4f, time_window=%d, time_unit=%s"
		", start_time_channel=%d, end_time_channel=%d, drop_missing_rate=%2.4f, buffer_first=%d\n", signalName.c_str(), output_name.c_str(),
		work_channel, factor, time_window, med_time_converter.type_to_string(time_unit).c_str(),
		start_time_channel, end_time_channel, drop_missing_rate, buffer_first);
}


//----------------------------------------------------------------------------------------
// RepHistoryLimit : given a signal : chomps history to be at a given window relative
//                   to prediction points
//----------------------------------------------------------------------------------------

// Fill req- and aff-signals vectors
//.......................................................................................
void RepHistoryLimit::init_lists() {

	req_signals.insert(signalName);
	aff_signals.insert(signalName);
}

// Init from map
//.......................................................................................
int RepHistoryLimit::init(map<string, string>& mapper)
{
	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		//! [RepBasicOutlierCleaner::init]
		if (field == "signal") { signalName = entry.second; }
		else if (field == "time_channel") time_channel = med_stoi(entry.second);
		else if (field == "win_from") win_from = med_stoi(entry.second);
		else if (field == "win_to") win_to = med_stoi(entry.second);
		else if (field == "delete_sig") delete_sig = med_stoi(entry.second);
		else if (field == "rep_time_unit") rep_time_unit = med_time_converter.string_to_type(entry.second);
		else if (field == "win_time_unit") win_time_unit = med_time_converter.string_to_type(entry.second);
	}

	init_lists();

	return 0;
}

int RepHistoryLimit::get_sub_usv_data(UniversalSigVec &usv, int from_time, int to_time, vector<char> &data, int &len)
{
	data.clear();
	len = 0;
	char *udata = (char *)usv.data;
	int element_size = (int)usv.size();
	for (int i = 0; i < usv.len; i++) {
		int i_time = usv.Time(i, time_channel);
		if (i_time > from_time && i_time <= to_time) {
			for (int j = element_size * i; j < element_size*(i + 1); j++)
				data.push_back(udata[j]);
			len++;
		}
	}
	return 0;
}

//---------------------------------------------------------------------------------------------------------------
int RepHistoryLimit::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat)
{

	// goal : for each time_points[i] generate the signal i version to contain only data within the given time window

	int len = 0;
	UniversalSigVec usv;
	vector<char> data;

	if (delete_sig == 0) {
		for (int ver = 0; ver < time_points.size(); ver++) {
			rec.uget(signalId, ver, usv);
			int curr_time = med_time_converter.convert_times(rep_time_unit, win_time_unit, time_points[ver]);
			int from_time = med_time_converter.convert_times(win_time_unit, rep_time_unit, curr_time - win_to);
			int to_time = med_time_converter.convert_times(win_time_unit, rep_time_unit, curr_time - win_from);
			get_sub_usv_data(usv, from_time, to_time, data, len);
			if (len < usv.len) {
				rec.set_version_data(signalId, ver, &data[0], len);
			}
		}
	}
	else {
		// simply delete signal and point all versions to the deleted signal
		rec.uget(signalId, 0, usv);
		rec.set_version_data(signalId, 0, &data[0], 0);
		for (int ver = 1; ver < time_points.size(); ver++)
			rec.point_version_to(signalId, 0, ver);
	}

	return 0;
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
	vector<FeatureGenerator *> noGenerators;
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
