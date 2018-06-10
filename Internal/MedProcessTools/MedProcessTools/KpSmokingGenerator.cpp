#include "KpSmokingGenerator.h"
#include <boost/algorithm/string/predicate.hpp>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace boost;

void KpSmokingGenerator::set_names() {
	names.clear();
	unordered_set<string> legal_features({ "Current_Smoker", "Ex_Smoker", "Unknown_Smoker","Never_Smoker", "Passive_Smoker","Smoke_Days_Since_Quitting", "Smoke_Pack_Years_Max", "Smoke_Pack_Years_Last" });

	if (raw_feature_names.size() == 0)
		MTHROW_AND_ERR("KpSmokingGenerator got no smoking_features");
	for (string s : raw_feature_names) {
		if (legal_features.find(s) == legal_features.end())
			MTHROW_AND_ERR("KpSmokingGenerator does not know how to generate [%s]", s.c_str());
		names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);
	}
}

int KpSmokingGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [SmokingGenerator::init]
		if (field == "smoking_features")
			boost::split(raw_feature_names, entry.second, boost::is_any_of(","));
		else if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
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
	return 0;
}

int KpSmokingGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num)
{
	int quitTime = (int)missing_val, unknownSmoker = 1, neverSmoker = 0, passiveSmoker = 0, formerSmoker = 0, currentSmoker = 0;
	float smokingStatus = missing_val, lastPackYears = missing_val, maxPackYears = missing_val, daysSinceQuitting = missing_val;
	UniversalSigVec smokingStatusUsv, quitTimeUsv, SmokingPackYearsUsv;
	int version = 0;
	string prevStatus = "Never";
	int prevStatusDate = 19000101;

	for (int i = 0; i < num; i++) {
		quitTime = (int)missing_val;
		unknownSmoker = 1, neverSmoker = 0, passiveSmoker = 0, formerSmoker = 0, currentSmoker = 0;
		smokingStatus = missing_val, lastPackYears = missing_val, maxPackYears = missing_val, daysSinceQuitting = missing_val;
		prevStatus = "Never";
		prevStatusDate = 19000101;

		// get signals:
		rec.uget("Smoking_Status", i, smokingStatusUsv);
		rec.uget("Smoking_Quit_Date", i, quitTimeUsv);
		rec.uget("Pack_Years", i, SmokingPackYearsUsv);

		int testDate = med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time);
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
		if (lastPackYears == missing_val) {	maxPackYears = missing_val;}

		if (neverSmoker == 1)
		{
			quitTime = (int)KP_NEVER_SMOKER_QUIT_TIME;
			maxPackYears = 0;
			lastPackYears = 0;
		}

		if (neverSmoker == 1)
			smokingStatus = 0;
		if (passiveSmoker)
			smokingStatus = 1;
		if (formerSmoker)
			smokingStatus = 2;
		if (currentSmoker)
			smokingStatus = 3;

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
		if (p_data[SMX_KP_CURRENT_SMOKER] != NULL) p_data[SMX_KP_CURRENT_SMOKER][index + i] = (float)currentSmoker;
		// Ex_Smoker
		if (p_data[SMX_KP_EX_SMOKER] != NULL) p_data[SMX_KP_EX_SMOKER][index + i] = (float)formerSmoker;
		// Smoke_Days_Since_Quitting
		if (p_data[SMX_KP_DAYS_SINCE_QUITTING] != NULL) p_data[SMX_KP_DAYS_SINCE_QUITTING][index + i] = (float)daysSinceQuitting;
		// Smok_Pack_Years_max
		if (p_data[SMX_KP_SMOK_PACK_YEARS_MAX] != NULL) p_data[SMX_KP_SMOK_PACK_YEARS_MAX][index + i] = (float)maxPackYears;
		// last pack years
		if (p_data[SMX_KP_SMOK_PACK_YEARS_LAST] != NULL) p_data[SMX_KP_SMOK_PACK_YEARS_LAST][index + i] = (float)lastPackYears;
		// Never_Smoker
		if (p_data[SMX_KP_NEVER_SMOKER] != NULL) p_data[SMX_KP_NEVER_SMOKER][index + i] = (float)neverSmoker;
		// Unknown_Smoker
		if (p_data[SMX_KP_UNKNOWN_SMOKER] != NULL) p_data[SMX_KP_UNKNOWN_SMOKER][index + i] = (float)unknownSmoker;
		// Passive_Smoker
		if (p_data[SMX_KP_PASSIVE_SMOKER] != NULL) p_data[SMX_KP_PASSIVE_SMOKER][index + i] = (float)passiveSmoker;

	}
	return 0;
}

void KpSmokingGenerator::get_p_data(MedFeatures& features) {
	p_data.resize(SMX_KP_LAST, NULL);

	if (iGenerateWeights) {
		if (names.size() != 1)
			MTHROW_AND_ERR("Cannot generate weights using a multi-feature generator (type %d generates %d features)\n", generator_type, (int)names.size())
		else
			p_data[0] = &(features.weights[0]);
	}

	for (string &name : names) {
		if (algorithm::ends_with(name, "Current_Smoker"))
			p_data[SMX_KP_CURRENT_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Ex_Smoker"))
			p_data[SMX_KP_EX_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smoke_Days_Since_Quitting"))
			p_data[SMX_KP_DAYS_SINCE_QUITTING] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smoke_Pack_Years_Max"))
			p_data[SMX_KP_SMOK_PACK_YEARS_MAX] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smoke_Pack_Years_Last"))
			p_data[SMX_KP_SMOK_PACK_YEARS_LAST] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Never_Smoker"))
			p_data[SMX_KP_NEVER_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Unknown_Smoker"))
			p_data[SMX_KP_UNKNOWN_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Passive_Smoker"))
			p_data[SMX_KP_PASSIVE_SMOKER] = &(features.data[name][0]);
		else
			MTHROW_AND_ERR("unknown feature name [%s]", name.c_str());
	}
}

