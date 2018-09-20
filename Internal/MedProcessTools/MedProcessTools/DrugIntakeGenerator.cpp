#include "DrugIntakeGenerator.h"

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//=======================================================================================
// FutureDrugsIntake : Look at drug intake within time-window (IN THE FUTURE !)
//=======================================================================================
// 

void DrugIntakeGenerator::set_names() {

	names.clear();

	if (names.empty()) {
		string name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + signalName + ".";
		string set_names = in_set_name;
		if (set_names == "" && this->sets.size() > 0)
			set_names = boost::algorithm::join(this->sets, "_");

		name += set_names + ".win_" + std::to_string(win_from) + "_" + std::to_string(win_to);
		names.push_back(name);
		//MLOG("Created %s\n", name.c_str());
	}

}

// Init
//.......................................................................................
void DrugIntakeGenerator::init_defaults() {
	generator_type = FTR_GEN_DRG_INTAKE;
	signalId = -1;
	time_unit_sig = MedTime::Undefined;
	time_unit_win = global_default_windows_time_unit;
	signalName = "";
};

// Generate
//.......................................................................................
int DrugIntakeGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {

	string& name = names[0];

	if (time_unit_sig == MedTime::Undefined)	time_unit_sig = rec.my_base_rep->sigs.Sid2Info[signalId].time_unit;

	float *p_feat = iGenerateWeights ? &(features.weights[index]) : &(features.data[name][index]);
	for (int i = 0; i < num; i++) {
		rec.uget(signalId, i);
		p_feat[i] = get_value(rec, rec.usv, med_time_converter.convert_times(features.time_unit, time_unit_win, features.samples[index + i].time),
			med_time_converter.convert_times(features.time_unit, time_unit_win, features.samples[index + i].outcomeTime));
	}

	return 0;
}

// Init look-up table
//.......................................................................................
void DrugIntakeGenerator::init_tables(MedDictionarySections& dict) {

	if (lut.size() == 0) {
		int section_id = dict.section_id(signalName);
		assert(sets.size() < 255); // Make sure we're fine with unsigned chars.
		//MLOG("BEFORE_LEARN:: signalName %s section_id %d sets size %d sets[0] %s\n", signalName.c_str(), section_id, sets.size(), sets[0].c_str());
		dict.prep_sets_indexed_lookup_table(section_id, sets, lut);
		//MLOG("AFTER_LEARN:: signalName %s section_id %d sets size %d sets[0] %s LUT %d\n", signalName.c_str(), section_id, sets.size(), sets[0].c_str(), lut.size());
	}

	return;
}

// Get Value
//.......................................................................................
float DrugIntakeGenerator::get_value(PidDynamicRec &rec, UniversalSigVec &usv, int time, int sig_outcomeTime)
{

	int min_time = time - win_to;
	int max_time = time - win_from;
	if (sig_outcomeTime < max_time)
		max_time = sig_outcomeTime;

	// Check Drugs
	vector<vector<pair<int, int> > > drugIntakePeriods(sets.size());

	// Collect period of administration per drug
	for (int i = 0; i < usv.len; i++) {
		if (lut[(int)usv.Val(i,0)] > 0) {
			int setId = lut[(int)usv.Val(i, 0)] - 1;
			int currentTime = med_time_converter.convert_times(time_unit_sig,time_unit_win,usv.Time(i));
		
			int period = (int)usv.Val(i,1);

			if (drugIntakePeriods[setId].empty() || currentTime > drugIntakePeriods[setId].back().second)
				drugIntakePeriods[setId].push_back({ currentTime,currentTime + period - 1 });
			else
				drugIntakePeriods[setId].back().second = currentTime + period - 1;

			if (currentTime > max_time)
				break;
		}
	}

	// Collect total period of administrating anti-coagulants
	vector<pair<int, int> > periods;
	for (int i = 0; i < drugIntakePeriods.size(); i++)
		periods.insert(periods.end(), drugIntakePeriods[i].begin(), drugIntakePeriods[i].end());

	sort(periods.begin(), periods.end(), [](const pair<int, int> &v1, const pair<int, int> &v2) {return (v1.first < v2.first); });

	int adminTime = 0;
	int lastCovered = -1;
	for (int i = 0; i < periods.size(); i++) {
 		if (periods[i].second < min_time)
			continue;

		if (periods[i].first < min_time)
			periods[i].first = min_time;

		if (periods[i].second > max_time)
			periods[i].second = max_time;

		if (periods[i].first > periods[i].second)
			continue;

		if (lastCovered == -1 || periods[i].first > lastCovered)
			adminTime += periods[i].second - periods[i].first;
		else if (periods[i].second > lastCovered)
			adminTime += periods[i].second - lastCovered;

		if (periods[i].second > lastCovered)
			lastCovered = periods[i].second;
	}

	float coverage = ((float)adminTime) / (float)(max_time + 1 - min_time);
	return coverage;
}

// Init
//.......................................................................................
int DrugIntakeGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [DrugIntakeGenerator::init]
		if (field == "win_from") win_from = med_stoi(entry.second);
		else if (field == "win_to") win_to = med_stoi(entry.second);
		else if (field == "signalName" || field == "signal") signalName = entry.second;
		else if (field == "time_unit") time_unit_win = med_time_converter.string_to_type(entry.second);
		else if (field == "sets") boost::split(sets, entry.second, boost::is_any_of(","));
		else if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "in_set_name") in_set_name = entry.second;
		else if (field == "weights_generator") iGenerateWeights = med_stoi(entry.second);
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for DrugIntakeGenerator\n", field.c_str());
		//! [DrugIntakeGenerator::init]
	}

	names.clear();
	set_names();

	req_signals.assign(1, signalName);

	return 0;
}


