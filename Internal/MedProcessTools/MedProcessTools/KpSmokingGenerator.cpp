#include "KpSmokingGenerator.h"
#include <boost/algorithm/string/predicate.hpp>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace boost;

void KpSmokingGenerator::set_names() {
	names.clear();
	unordered_set<string> legal_features({ "Current_Smoker", "Ex_Smoker", "Unknown_Smoker","Never_Smoker", "Passive_Smoker","Smoke_Days_Since_Quitting", "Smoke_Pack_Years_Max", "Smoke_Pack_Years_Last","NLST_Criterion" });

	if (raw_feature_names.size() == 0)
		MTHROW_AND_ERR("KpSmokingGenerator got no smoking_features");
	for (string s : raw_feature_names) {
		if (legal_features.find(s) == legal_features.end())
			MTHROW_AND_ERR("KpSmokingGenerator does not know how to generate [%s]", s.c_str());
		if ((s == "NLST_Criterion") && (nonDefaultNlstCriterion == true))
			names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".Screening_Criterion" + s + "_min_age_" + to_string((int)nlstMinAge) + "_max_age_" + to_string((int)nlstMaxAge) + "_pack_years_" + to_string((int)nlstPackYears));
		else
			names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);
	}
}

int KpSmokingGenerator::init(map<string, string>& mapper) {

	// Set NLST default values:
	nlstMinAge = 55;
	nlstMaxAge = 80;
	nlstPackYears = 30;
	nlstQuitTimeYears = 15;
	nonDefaultNlstCriterion = false;

	for (auto entry : mapper) {
		string field = entry.first;
		//! [SmokingGenerator::init]
		if (field == "smoking_features")
			boost::split(raw_feature_names, entry.second, boost::is_any_of(","));
		else if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
		else if (entry.first == "min_age") {
			nlstMinAge = stof(entry.second);
			nonDefaultNlstCriterion = true;
		}
		else if (entry.first == "max_age") {
			nlstMaxAge = stof(entry.second);
			nonDefaultNlstCriterion = true;
		}
		else if (entry.first == "pack_years") {
			nlstPackYears = stof(entry.second);
			nonDefaultNlstCriterion = true;
		}
		else if (entry.first == "quit_time_years") {
			nlstQuitTimeYears = stof(entry.second);
			nonDefaultNlstCriterion = true;
		}
		else if (field == "weights_generator")
			iGenerateWeights = stoi(entry.second);
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter \'%s\' for KpSmokingGenerator\n", field.c_str());
		//! [SmokingGenerator::init]

	}

	set_names();
	req_signals.clear();
	req_signals.push_back("Smoking_Status");
	req_signals.push_back("Smoking_Quit_Date");
	req_signals.push_back("Pack_Years");
	req_signals.push_back("BDATE");
	return 0;
}

int KpSmokingGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data)
{
	int quitTime = (int)missing_val, unknownSmoker = 1, neverSmoker = 0, passiveSmoker = 0, formerSmoker = 0, currentSmoker = 0;
	float lastPackYears = missing_val, daysSinceQuitting = missing_val;
	UniversalSigVec smokingStatusUsv, quitTimeUsv, SmokingPackYearsUsv, bdateUsv;
	MedTime &tm = med_time_converter;
	string prevStatus = "Never";
	int prevStatusDate = 19000101;

	for (int i = 0; i < num; i++) {
		quitTime = (int)missing_val;
		unknownSmoker = 1, neverSmoker = 0, passiveSmoker = 0, formerSmoker = 0, currentSmoker = 0;
		lastPackYears = missing_val, daysSinceQuitting = missing_val;
		prevStatus = "Never";
		prevStatusDate = 19000101;

		// get signals:
		rec.uget("Smoking_Status", i, smokingStatusUsv);
		rec.uget("Smoking_Quit_Date", i, quitTimeUsv);
		rec.uget("Pack_Years", i, SmokingPackYearsUsv);
		rec.uget("BDATE", i, bdateUsv);
		int birthDate;
		if (bdateUsv.n_val_channels() > 0)
			birthDate = bdateUsv.Val(0);
		else
			birthDate = bdateUsv.Time(0);

		int testDate = med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time);

		// calc age
		int age = (int)((float)tm.diff_times(testDate, birthDate, MedTime::Date, MedTime::Days) / 365.0);
		// If Smoking Status vec exists, then status is known.
		if (smokingStatusUsv.len > 0)
		{
			if (smokingStatusUsv.Time(0) <= testDate)
			{
				unknownSmoker = 0;
				neverSmoker = 1;
				passiveSmoker = 0;
				formerSmoker = 0;
				currentSmoker = 0;
			}
		}
		// Get last quit time before test date.
		for (int timeInd = 0; timeInd < quitTimeUsv.len; timeInd++)
		{
			if (quitTimeUsv.Time(timeInd) > testDate) { break; }
			quitTime = (int)quitTimeUsv.Val(timeInd);
			neverSmoker = 0;
			formerSmoker = 1;
		}

		// Calculate smoking status   
		int smokingStatusSid = rec.my_base_rep->dict.section_id("Smoking_Status");

		for (int timeInd = 0; timeInd < smokingStatusUsv.len; timeInd++)
		{
			if (smokingStatusUsv.Time(timeInd) > testDate) { break; }
			string sigVal = rec.my_base_rep->dict.name(smokingStatusSid, (int)smokingStatusUsv.Val(timeInd));

			// If has Quit or Yes, must be Smoker or Ex smoker. will be set later according to last value
			if ((sigVal == "Quit") | (sigVal == "Yes"))
			{
				neverSmoker = 0;
				formerSmoker = 1;
			}
			// Check if also Passive
			if (sigVal == "Passive")
			{
				passiveSmoker = 1;
			}

			// Check Whether Smoked (Yest) after last quitTime - if so, set quit time to the quitting day.
			if ((prevStatus == "Yes") & (sigVal == "Quit"))
			{
				if (prevStatusDate > quitTime)
				{
					quitTime = (int)smokingStatusUsv.Time(timeInd);
				}
			}
			prevStatus = sigVal;
			prevStatusDate = smokingStatusUsv.Time(timeInd);
		}

		// if last value is "Yes" then this is a current smoker.
		if (prevStatus == "Yes")
		{
			formerSmoker = 0;
			currentSmoker = 1;
			quitTime = testDate;
		}

		// Pack Years
		float maxPackYears = 0;
		for (int timeInd = 0; timeInd < SmokingPackYearsUsv.len; timeInd++)
		{
			if (SmokingPackYearsUsv.Time(timeInd) > testDate) { break; }
			if (SmokingPackYearsUsv.Val(timeInd) > 0)
			{
				lastPackYears = SmokingPackYearsUsv.Val(timeInd);
				if (SmokingPackYearsUsv.Val(timeInd) > maxPackYears)
				{
					maxPackYears = SmokingPackYearsUsv.Val(timeInd);
				}
			}
			neverSmoker = 0;
			if (currentSmoker == 0) {
				formerSmoker = 1;
			}
		}
		// This means that there wasn't any value.
		if (lastPackYears == missing_val) { maxPackYears = missing_val; }

		if (neverSmoker == 1)
		{
			quitTime = (int)KP_NEVER_SMOKER_QUIT_TIME;
			maxPackYears = 0;
			lastPackYears = 0;
		}

		// Calculate time since quitting
		if (quitTime != (int)missing_val)
		{
			//	MLOG("quit time: %d %d \n", testDate, (int)quitTime);
			daysSinceQuitting = (float)med_time_converter.diff_times(testDate, quitTime, MedTime::Date, MedTime::Days);
		}
		else
		{
			daysSinceQuitting = missing_val;
		}

		// Add data to matrix:
		// Current_Smoker
		if (_p_data[SMX_KP_CURRENT_SMOKER] != NULL) _p_data[SMX_KP_CURRENT_SMOKER][index + i] = (float)currentSmoker;
		// Ex_Smoker
		if (_p_data[SMX_KP_EX_SMOKER] != NULL) _p_data[SMX_KP_EX_SMOKER][index + i] = (float)formerSmoker;
		// Smoke_Days_Since_Quitting
		if (_p_data[SMX_KP_DAYS_SINCE_QUITTING] != NULL) _p_data[SMX_KP_DAYS_SINCE_QUITTING][index + i] = (float)daysSinceQuitting;
		// Smok_Pack_Years_max
		if (_p_data[SMX_KP_SMOK_PACK_YEARS_MAX] != NULL) _p_data[SMX_KP_SMOK_PACK_YEARS_MAX][index + i] = (float)maxPackYears;
		// last pack years
		if (_p_data[SMX_KP_SMOK_PACK_YEARS_LAST] != NULL) _p_data[SMX_KP_SMOK_PACK_YEARS_LAST][index + i] = (float)lastPackYears;
		// Never_Smoker
		if (_p_data[SMX_KP_NEVER_SMOKER] != NULL) _p_data[SMX_KP_NEVER_SMOKER][index + i] = (float)neverSmoker;
		// Unknown_Smoker
		if (_p_data[SMX_KP_UNKNOWN_SMOKER] != NULL) _p_data[SMX_KP_UNKNOWN_SMOKER][index + i] = (float)unknownSmoker;
		// Passive_Smoker
		if (_p_data[SMX_KP_PASSIVE_SMOKER] != NULL) _p_data[SMX_KP_PASSIVE_SMOKER][index + i] = (float)passiveSmoker;
		//NLST_Criterion (only after everything was calculated, we can calc. the NLST criterion)
		if (_p_data[NLST_CRITERION] != NULL) _p_data[NLST_CRITERION][index + i] = (float)calcNlst(age, unknownSmoker, daysSinceQuitting, lastPackYears);
	}
	return 0;
}

int KpSmokingGenerator::calcNlst(int age, int unknownSmoker, int daysSinceQuitting, float lastPackYears)
{
	if ((unknownSmoker == 1) || (daysSinceQuitting == (int)missing_val) || (lastPackYears == missing_val)) {
		return missing_val;
	}
	return ((age >= nlstMinAge) && (age <= nlstMaxAge) && (lastPackYears >= nlstPackYears) && (daysSinceQuitting <= nlstQuitTimeYears * 365.0));
}

void KpSmokingGenerator::get_p_data(MedFeatures& features, vector<float *> &_p_data) {
	p_data.resize(SMX_KP_LAST, NULL);

	if (iGenerateWeights) {
		if (names.size() != 1)
			MTHROW_AND_ERR("Cannot generate weights using a multi-feature generator (type %d generates %d features)\n", generator_type, (int)names.size())
		else
			p_data[0] = &(features.weights[0]);
	}

	for (string &name : names) {
		if (algorithm::ends_with(name, "Current_Smoker"))
			_p_data[SMX_KP_CURRENT_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Ex_Smoker"))
			_p_data[SMX_KP_EX_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smoke_Days_Since_Quitting"))
			_p_data[SMX_KP_DAYS_SINCE_QUITTING] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smoke_Pack_Years_Max"))
			_p_data[SMX_KP_SMOK_PACK_YEARS_MAX] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smoke_Pack_Years_Last"))
			_p_data[SMX_KP_SMOK_PACK_YEARS_LAST] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Never_Smoker"))
			_p_data[SMX_KP_NEVER_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Unknown_Smoker"))
			_p_data[SMX_KP_UNKNOWN_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Passive_Smoker"))
			_p_data[SMX_KP_PASSIVE_SMOKER] = &(features.data[name][0]);
		else if (algorithm::contains(name, "NLST_Criterion"))
			_p_data[NLST_CRITERION] = &(features.data[name][0]);
		else
			MTHROW_AND_ERR("unknown feature name [%s]", name.c_str());
	}
}
