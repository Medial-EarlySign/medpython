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
		string name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + signalName + ".category_intake_";
		string set_names = in_set_name;
		if (set_names == "" && this->sets.size() > 0)
			set_names = boost::algorithm::join(this->sets, "_");

		name += set_names + ".win_" + std::to_string(win_from) + "_" + std::to_string(win_to);
		if (timeRangeSignalName != "")
			name += ".time_range_" + timeRangeSignalName + "_" + time_range_type_to_name(timeRangeType);
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
	timeRangeSignalId = -1;
};

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

void DrugIntakeGenerator::get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const {
	signal_categories_in_use[signalName] = sets;
}

// Generate
//.......................................................................................
int DrugIntakeGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {

	string& name = names[0];

	if (time_unit_sig == MedTime::Undefined)	time_unit_sig = rec.my_base_rep->sigs.Sid2Info[signalId].time_unit;
	if (timeRangeSignalName != "" && time_unit_range_sig == MedTime::Undefined)
		time_unit_range_sig = rec.my_base_rep->sigs.Sid2Info[timeRangeSignalId].time_unit;

	float *p_feat = iGenerateWeights ? &(features.weights[index]) : &(features.data[name][index]);
	for (int i = 0; i < num; i++) {
		p_feat[i] = get_value(rec, i, med_time_converter.convert_times(features.time_unit, time_unit_win, features.samples[index + i].time),
			bound_outcomeTime ?
			med_time_converter.convert_times(features.time_unit, time_unit_win, features.samples[index + i].outcomeTime) : -1);
	}

	return 0;
}

// Get Value
//.......................................................................................
float DrugIntakeGenerator::get_value(PidDynamicRec &rec, int idx, int time, int sig_outcomeTime)
{

	rec.uget(signalId, idx);

	int updated_win_from = win_from, updated_win_to = win_to;
	int updated_d_win_from, updated_d_win_to;
	if (timeRangeSignalId != -1) {
		UniversalSigVec time_range_usv;
		rec.uget(timeRangeSignalId, idx, time_range_usv);
		get_updated_time_window(time_range_usv, timeRangeType, time_unit_range_sig, time_unit_win, time_unit_sig, time,
			win_from, updated_win_from, win_to, updated_win_to, false, 0, updated_d_win_from, 0, updated_d_win_to);
	}

	int min_time = time - updated_win_to;
	int max_time = time - updated_win_from;
	if (bound_outcomeTime && sig_outcomeTime < max_time)
		max_time = sig_outcomeTime;

	//MLOG("rugIntakeGenerator::get_value -> pid %d time %d min_time %d max_time %d sig_outcomeTime %d\n", rec.pid, time, min_time, max_time, sig_outcomeTime);

	if (max_time < min_time)
		return MED_MAT_MISSING_VALUE;

	// Check Drugs
	vector<vector<pair<int, int> > > drugIntakePeriods(sets.size());

	// Collect periods of administration per drug
	for (int i = 0; i < rec.usv.len; i++) {
		if (lut[(int)rec.usv.Val(i,0)] > 0) {
			int setId = lut[(int)rec.usv.Val(i, 0)] - 1;
			int currentTime = med_time_converter.convert_times(time_unit_sig,time_unit_win,rec.usv.Time(i));
		
			int period = (int)rec.usv.Val(i,1);

			if (drugIntakePeriods[setId].empty() || currentTime > drugIntakePeriods[setId].back().second)
				drugIntakePeriods[setId].push_back({ currentTime,currentTime + period - 1 });
			else
				drugIntakePeriods[setId].back().second = currentTime + period - 1;

			if (currentTime > max_time)
				break;
		}
	}

	// Collect total period of administrating drugs
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
		else if (field == "time_range_signal") timeRangeSignalName = entry.second;
		else if (field == "time_range_signal_type") timeRangeType = time_range_name_to_type(entry.second);
		else if (field == "bound_outcomeTime") bound_outcomeTime = stoi(entry.second) > 0;
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for DrugIntakeGenerator\n", field.c_str());
		//! [DrugIntakeGenerator::init]
	}

	names.clear();
	set_names();

	req_signals.assign(1, signalName);
	if (timeRangeSignalName != "")
		req_signals.push_back(timeRangeSignalName);

	return 0;
}


