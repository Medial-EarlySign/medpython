#include "ElixhauserGenerator.h"
#include <boost/algorithm/string/predicate.hpp>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

using namespace boost;

// Parse sets name string
// Expected format "type=set1,set2,set3;type=...."
void ElixhauserGenerator::parseSets(string& init_string, vector<vector<string>>& sets) {

	// Prepare map  - type -> index
	map<string, int> type2index;
	for (int i = 0; i < types.size(); i++)
		type2index[types[i]] = i;

	// Parse
	init_string.erase(remove(init_string.begin(), init_string.end(), ' '), init_string.end());

	vector<string> tokens, keyVal;
	boost::split(tokens, init_string, boost::is_any_of(";"));
	for (string& token : tokens) {
		boost::split(keyVal, token, boost::is_any_of("="));
		if (keyVal.size() != 2)
			MTHROW_AND_ERR("Cannot parse token %s from init-string %s\n", token.c_str(), init_string.c_str());
		if (type2index.find(keyVal[0]) == type2index.end())
			MTHROW_AND_ERR("Unknown type %s in token %s from init_string %s\n", keyVal[0].c_str(), token.c_str(), init_string.c_str());

		boost::split(sets[type2index[keyVal[0]]], keyVal[1], boost::is_any_of(","));
	}
}

// Init default sets
void ElixhauserGenerator::initSets() {

}

// Init from String
int ElixhauserGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [SmokingGenerator::init]
		if (field == "drg_signal")
			drgSignalName = entry.second;
		else if (field == "diagnosis_signal")
			diagSignalName = entry.second;
		else if (field == "drg_sets")
			parseSets(entry.second, drgSets);
		else if (field == "diagnosis_sets")
			parseSets(entry.second, diagSets);
		else if (field == "type") {
			type = entry.second;
			boost::algorithm::to_lower(type);
			if (_weights.find(type) == _weights.end())
				MTHROW_AND_ERR("Unknown Elixhauser type \'%s\'\n", type.c_str());
		}
		else if (field == "name" || field == "names")
			names = { entry.second };
		else if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "weights_generator")
			iGenerateWeights = med_stoi(entry.second);
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter \'%s\' for SmokingGenerator\n", field.c_str());
		//! [SmokingGenerator::init]

	}
	set_names();
	req_signals.clear();
	req_signals.push_back(drgSignalName);
	req_signals.push_back(diagSignalName);
	
	return 0;
}

// Init look-up tables
//.......................................................................................
void ElixhauserGenerator::init_tables(MedDictionarySections& dict) {

	// Sections
	int drgSection = dict.section_id(drgSignalName);
	int diagSection = dict.section_id(diagSignalName);

	drgLuts.resize(types.size());
	diagLuts.resize(types.size());

	for (size_t i = 0; i < types.size(); i++) {
		dict.prep_sets_indexed_lookup_table(drgSection, drgSets[i], drgLuts[i]);
		dict.prep_sets_indexed_lookup_table(diagSection, diagSets[i], diagLuts[i]);
	}

	weights = _weights[type];
}

// Learn a generator - just time units
//.......................................................................................
int ElixhauserGenerator::_learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) {
	drg_time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(drgSignalName)].time_unit;
	diag_time_unit_sig = rep.sigs.Sid2Info[rep.sigs.sid(diagSignalName)].time_unit;
	
	return 0;
}

int ElixhauserGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {

	if (drg_time_unit_sig == MedTime::Undefined)	drg_time_unit_sig = rec.my_base_rep->sigs.Sid2Info[drgSignalId].time_unit;
	if (diag_time_unit_sig == MedTime::Undefined)	diag_time_unit_sig = rec.my_base_rep->sigs.Sid2Info[diagSignalId].time_unit;

	float *p_feat = _p_data[0] + index;
	MedSample *p_samples = &(features.samples[index]);

	UniversalSigVec drg, diag;
	int min_time, max_time, time, outcomeTime;
	for (int i = 0; i < num; i++) {
		time = med_time_converter.convert_times(features.time_unit, time_unit_win, p_samples[i].time);
		outcomeTime = med_time_converter.convert_times(features.time_unit, time_unit_win, p_samples[i].outcomeTime);

		// DRG
		get_window_in_sig_time(win_from, win_to, time_unit_win, drg_time_unit_sig, time, min_time, max_time, true, outcomeTime);
		rec.uget(drgSignalId, i, drg);

		vector<int> drgIndices(types.size(), 0);
		for (int i = 0; i < drg.len; i++) {
			int itime = drg.Time(i);
			if (itime > max_time) break;
			if (itime >= min_time) {
				for (size_t itype = 0; itype < types.size(); itype++) {
					if (drgLuts[itype][drg.Val<int>(i, 0)])
						drgIndices[itype] = 1;
				}
			}
		}

		// Diagnosis
		get_window_in_sig_time(win_from, win_to, time_unit_win, diag_time_unit_sig, time, min_time, max_time, true, outcomeTime);
		rec.uget(diagSignalId, i, diag);

		vector<int> diagIndices(types.size(), 0);
		for (int i = 0; i < diag.len; i++) {
			int itime = diag.Time(i);
			if (itime > max_time) break;
			if (itime >= min_time) {
				for (size_t itype = 0; itype < types.size(); itype++) {
					if (diagLuts[itype][diag.Val<int>(i, 0)])
						diagIndices[itype] = 1;
				}
			}
		}

		// Elixhauser
		p_feat[i] = 0;
		for (size_t itype = 0; itype < types.size(); itype++) {
			if (!drgIndices[itype] && diagIndices[itype])
				p_feat[i] += weights[itype];
		}

	}
	return 0;
}
