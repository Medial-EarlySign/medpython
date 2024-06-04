#include "UnifiedSmokingGenerator.h"
#include <boost/algorithm/string/predicate.hpp>
#include <algorithm>
#include <omp.h>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace boost;

void UnifiedSmokingGenerator::set_names() {
	names.clear();
	unordered_set<string> legal_features({ "debug_file", "Current_Smoker", "Ex_Smoker", "Unknown_Smoker","Never_Smoker", "Passive_Smoker","Smok_Days_Since_Quitting","Smok_Years_Since_Quitting", "Smok_Pack_Years_Max", "Smok_Pack_Years", "Smok_Pack_Years_Last","NLST_Criterion","Smoking_Intensity", "Smoking_Years" });

	if (raw_feature_names.size() == 0)
		MTHROW_AND_ERR("UnifiedSmokingGenerator got no smoking_features");
	for (string s : raw_feature_names) {
		if (legal_features.find(s) == legal_features.end())
			MTHROW_AND_ERR("UnifiedSmokingGenerator does not know how to generate [%s]", s.c_str());
		if ((s == "NLST_Criterion") && (nonDefaultNlstCriterion == true))
			names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".Smoking." + s + "_min_age_" + to_string((int)nlstMinAge) + "_max_age_" + to_string((int)nlstMaxAge) + "_pack_years_" + to_string((int)nlstPackYears) + "_quit_years_" + to_string((int)nlstQuitTimeYears));
		else
			names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + ".Smoking." + s);
	}
}

int UnifiedSmokingGenerator::init(map<string, string>& mapper) {

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
		else if (entry.first == "use_data_complition") {
			useDataComplition = stof(entry.second) > 0;
			cout << endl << "######################### useDataComplition was changed to " << useDataComplition << endl << endl;
		}
		else if (field == "timeSinceQuittingModelSlope")
			timeSinceQuittingModelSlope = stof(entry.second);
		else if (field == "timeSinceQuittingModelConst")
			timeSinceQuittingModelConst = stof(entry.second);
		else if (field == "weights_generator")
			iGenerateWeights = stoi(entry.second);
		else if (field == "debug_file")
			debug_file = entry.second;
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter \'%s\' for UnifiedSmokingGenerator\n", field.c_str());
		//! [SmokingGenerator::init]
	}

	set_names();
	req_signals.clear();
	req_signals.push_back("BDATE");
	req_signals.push_back("Smoking_Status");

	req_signals.push_back("Smoking_Quit_Date");
	req_signals.push_back("Pack_Years");
	req_signals.push_back("Smoking_Intensity");
	req_signals.push_back("Smoking_Duration");

	//char *filename = "W:/Users/Ron/Projects/LungCancer/results/unified_smoking/tmp_output4.tsv";
	if (debug_file != "")
	{
		fp = fopen(debug_file.c_str(), "w");
		if (fp == NULL)
			MTHROW_AND_ERR("Unable to open file: %s \n", debug_file.c_str());
	}
	return 0;
}

int UnifiedSmokingGenerator::update(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "tags")
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
		else if (entry.first == "use_data_complition") {
			useDataComplition = stof(entry.second) > 0;
			cout << endl << "#########################(u) useDataComplition was changed to " << useDataComplition << endl << endl;
		}
		else if (field == "weights_generator")
			iGenerateWeights = stoi(entry.second);
		else if (field == "debug_file")
			debug_file = entry.second;
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter \'%s\' for UnifiedSmokingGenerator\n", field.c_str());
	}

	if (debug_file != "")
	{
		fp = fopen(debug_file.c_str(), "w");
		if (fp == NULL)
			MTHROW_AND_ERR("Unable to open file: %s \n", debug_file.c_str());
	}

	return 0;
}


void UnifiedSmokingGenerator::init_tables(MedDictionarySections& dict) {
	//init luts:
	smoke_status_sec_id = dict.section_id("Smoking_Status");
	unordered_map<string, vector<string>> smoke_stat_categs;
	get_required_signal_categories(smoke_stat_categs);
	vector<string> smoke_st_vals = smoke_stat_categs.at("Smoking_Status");
	smoke_status_luts.resize(smoke_st_vals.size());
	for (size_t i = 0; i < smoke_status_luts.size(); ++i)
	{
		vector<string> set_vals = { smoke_st_vals[i] };
		if (dict.id(smoke_status_sec_id, smoke_st_vals[i]) < 0) {
			//MWARN("WARN UnifiedSmokingGenerator::init_tables - Smoking_Status \"%s\" is undefined\n",
			//	smoke_st_vals[i].c_str());
			set_vals.clear(); //Allow Unknown smoking status to be undefined
		}
		dict.prep_sets_lookup_table(smoke_status_sec_id, set_vals, smoke_status_luts[i]);
	}
}

int UnifiedSmokingGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data)
{
	int unknownSmoker = 1, neverSmoker = 0, passiveSmoker = 0, formerSmoker = 0, currentSmoker = 0;
	float lastPackYears = missing_val, daysSinceQuitting = missing_val, smokingDuration = missing_val, smokingIntensity = missing_val, yearsSinceQuitting = missing_val;
	float daysSinceQuittingOriginal = missing_val;
	UniversalSigVec smokingStatusUsv, quitTimeUsv, SmokingPackYearsUsv, bdateUsv, SmokingIntensityUsv, SmokingDurationUsv;
	//cout << endl << "NUM " << num << endl << endl;
	int birthDate, qa_print = 1;

	for (int i = 0; i < num; i++) {
		qa_print = 1;
		unknownSmoker = 1, neverSmoker = 0, passiveSmoker = 0, formerSmoker = 0, currentSmoker = 0;
		lastPackYears = daysSinceQuitting = yearsSinceQuitting = smokingIntensity = missing_val;
		daysSinceQuittingOriginal = missing_val;

		// get signals:
		rec.uget(smoking_status_id, i, smokingStatusUsv);
		if (smoking_quit_date_id > 0)
			rec.uget(smoking_quit_date_id, i, quitTimeUsv);
		if (smoking_pack_years_id > 0)
			rec.uget(smoking_pack_years_id, i, SmokingPackYearsUsv);
		if (smoking_intensity_id > 0)
			rec.uget(smoking_intensity_id, i, SmokingIntensityUsv);
		if (smoking_duration_id > 0)
			rec.uget(smoking_duration_id, i, SmokingDurationUsv);

		rec.uget(bdate_sid, i, bdateUsv);

		if (bdateUsv.n_val_channels() > 0)
			birthDate = bdateUsv.Val(0);
		else
			birthDate = bdateUsv.Time(0);

		int testDate = med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time);

		// calc age
		int age = (int)((float)med_time_converter.diff_times(testDate, birthDate, MedTime::Date, MedTime::Days) / 365.0);

		// Generate First and Last Dates and dates vector
		// Map between status to <First Date, Last Date>
		map<SMOKING_STATUS, pair<int, int>> smokingStatusDates = { { NEVER_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } } ,{ PASSIVE_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } },{ EX_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } },{ CURRENT_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } },{ NEVER_OR_EX_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } } };
		vector<int> dates = {};
		genFirstLastSmokingDates(rec, smokingStatusUsv, quitTimeUsv, testDate, smokingStatusDates, dates, birthDate);

		// Determine Smoking Status per date
		vector<pair<SMOKING_STATUS, int>> smokingStatusVec = { { UNKNOWN_SMOKER,  (int)bdateUsv.Val(0) } };

		//genSmokingVec(rec, smokingStatusUsv, smokingStatusVec, testDate, unknownSmoker, neverSmoker, passiveSmoker, formerSmoker, currentSmoker); KP logic
		genSmokingStatus(smokingStatusDates, dates, testDate, birthDate, smokingStatusVec);

		// Create Ranges:
		vector<RangeStatus> smokeRanges = {};
		genSmokingRanges(smokingStatusVec, testDate, birthDate, smokeRanges);

		// Extract Last Status
		genLastStatus(smokeRanges, unknownSmoker, neverSmoker, formerSmoker, currentSmoker, passiveSmoker);

		calcQuitTime(smokeRanges, formerSmoker, neverSmoker, currentSmoker, testDate, birthDate, daysSinceQuitting, yearsSinceQuitting);
		// years since quitting just from original data 
		calcQuitTimeOriginalData(rec, smokingStatusUsv, quitTimeUsv, testDate, formerSmoker, neverSmoker, currentSmoker, daysSinceQuittingOriginal);

		calcSmokingIntensity(SmokingIntensityUsv, testDate, neverSmoker, smokingIntensity);

		int lastPackYearsDate = missing_val;
		float maxPackYears = missing_val;
		calcPackYears(SmokingPackYearsUsv, testDate, neverSmoker, currentSmoker, formerSmoker, lastPackYearsDate, lastPackYears, maxPackYears);
		// pack years just from original data
		float lastPackYearsOriginal = lastPackYears;
		calcPackYearsOriginalData(testDate, lastPackYearsDate, lastPackYears, lastPackYearsOriginal, SmokingIntensityUsv, SmokingDurationUsv);

		//cout << i << " age " << age << endl;
		float smokingDurationBeforePackYears = missing_val;
		calcSmokingDuration(neverSmoker, unknownSmoker, smokeRanges, birthDate, lastPackYearsDate, SmokingDurationUsv, testDate, smokingDurationBeforePackYears, smokingDuration);

		fixPackYearsSmokingIntensity(smokingDurationBeforePackYears, smokingIntensity, smokingDuration, lastPackYears, maxPackYears);
		printDebug(smokeRanges, qa_print, smokingStatusUsv, SmokingIntensityUsv, birthDate, testDate, smokingStatusVec, rec, quitTimeUsv, SmokingPackYearsUsv, smokingIntensity, smokingDuration, yearsSinceQuitting, maxPackYears);

		addDataToMat(_p_data, index, i, age, currentSmoker, formerSmoker, daysSinceQuitting, daysSinceQuittingOriginal, maxPackYears, lastPackYears, lastPackYearsOriginal, neverSmoker, unknownSmoker, passiveSmoker, yearsSinceQuitting, smokingIntensity, smokingDuration);

	}
	return 0;
}

int UnifiedSmokingGenerator::calcNlst(int age, int unknownSmoker, int daysSinceQuitting, float lastPackYears)
{
	//cout << "nlstMinAge " << nlstMinAge << " nlstMaxAge " << nlstMaxAge << " age " << age << endl;
	//cout << "lastPackYears " << lastPackYears << endl << endl;
	if ((unknownSmoker == 1) || (daysSinceQuitting == (int)missing_val) || (lastPackYears == missing_val)) {
		return missing_val;
	}

	return ((age >= nlstMinAge) && (age <= nlstMaxAge) && (lastPackYears >= nlstPackYears) && (daysSinceQuitting <= nlstQuitTimeYears * 365.0));
}

void UnifiedSmokingGenerator::calcQuitTimeOriginalData(PidDynamicRec& rec, UniversalSigVec &smokingStatusUsv, UniversalSigVec &quitTimeUsv, int testDate, int formerSmoker, int neverSmoker, int currentSmoker, float &daysSinceQuittingOriginal)
// this function calculates daysSinceQuittingOriginal - daysSinceQuitting based only on original data without any model
// we use it in nlst elgibility check, when we do not allow complition by model
// 0 for current_smoker, large numeber for never_smoker
// for former_smoker, take the max between last reported quitDate and last report currentSmoker+1 
{
	int lastQuitTimeDate = missing_val;
	int lastCurrentDate = missing_val;
	int dateQuittingOriginal;
	string sigVal, inVal;

	daysSinceQuittingOriginal = -missing_val; // large number, default for any smokingstatus but formerSmoker

	if (currentSmoker == 1) { daysSinceQuittingOriginal = 0; }

	if (formerSmoker == 1) {
		// get last quit_time reported (might be missing)
		for (int timeInd = 0; timeInd < quitTimeUsv.len; timeInd++)
		{
			if (quitTimeUsv.Time(timeInd) > testDate) { break; }
			if (quitTimeUsv.Val(timeInd) > 0)
			{
				lastQuitTimeDate = quitTimeUsv.Val(timeInd);
			}
		}
		// get Last time current was reported
		for (int timeInd = 0; timeInd < smokingStatusUsv.len; timeInd++)
		{
			if (smokingStatusUsv.Time(timeInd) > testDate) { break; }
			if (smokingStatusUsv.Val(timeInd) > 0)
			{
				int sigVal_num = (int)smokingStatusUsv.Val(timeInd);
				string smoking_status = rec.my_base_rep->dict.name(smoke_status_sec_id, sigVal_num);
				if (smoking_status == "Current") {
					lastCurrentDate = smokingStatusUsv.Time(timeInd);
				}
				// cout << "lastQuitTimeDate " << lastQuitTimeDate << " timeInd " << smokingStatusUsv.Time(timeInd) << " lastCurrentDate " << lastCurrentDate << endl;
			}
		}
		dateQuittingOriginal = max(lastCurrentDate + 1, lastQuitTimeDate);
		// cout << "CALCULATING lastQuitTimeDate " << lastQuitTimeDate << " lastCurrentDate " << lastCurrentDate << " dateQuittingOriginal " << dateQuittingOriginal << endl;
		if (dateQuittingOriginal > 0) {
			daysSinceQuittingOriginal = (float)med_time_converter.diff_times(testDate, dateQuittingOriginal, MedTime::Date, MedTime::Days);
		}
	}
}

void UnifiedSmokingGenerator::calcPackYearsOriginalData(int testDate, int lastPackYearsDate, float lastPackYears, float &lastPackYearsOriginal, UniversalSigVec SmokingIntensityUsv, UniversalSigVec SmokingDurationUsv)
// this function calculates lastPackYearsOriginal - lastPackYears based only on original data without any model
// we use it in nlst elgibility check, when we do not allow complition by model
// lastPackYears is the last packYears signal
// we replace it with lastIntensity * lastDuration / PACK_SIZE when we have both signals and thay are more updated 
{
	int lastDuration = missing_val;
	int lastDurationDate = missing_val;
	int lastIntensity = missing_val;
	int lastIntensityDate = missing_val;

	for (int timeInd = 0; timeInd < SmokingDurationUsv.len; timeInd++)
	{
		if (SmokingDurationUsv.Time(timeInd) > testDate) { break; }
		if (SmokingDurationUsv.Val(timeInd) > 0)
		{
			lastDuration = SmokingDurationUsv.Val(timeInd);
			lastDurationDate = SmokingDurationUsv.Time(timeInd);
		}
	}

	for (int timeInd = 0; timeInd < SmokingIntensityUsv.len; timeInd++)
	{
		if (SmokingIntensityUsv.Time(timeInd) > testDate) { break; }
		if (SmokingIntensityUsv.Val(timeInd) > 0)
		{
			lastIntensity = SmokingIntensityUsv.Val(timeInd);
			lastIntensityDate = SmokingIntensityUsv.Time(timeInd);
		}
	}

	if ((lastDuration != missing_val) && (lastIntensity != missing_val)) {
		if ((lastPackYears == missing_val) || ((lastIntensityDate > lastPackYearsDate) && (lastDurationDate > lastPackYearsDate))) {
			lastPackYearsOriginal = lastIntensity * lastDuration / PACK_SIZE;
		}
	}
}

void UnifiedSmokingGenerator::get_p_data(MedFeatures& features, vector<float *> &_p_data) {
	p_data.resize(SMX_UNIFIED_LAST, NULL);

	if (iGenerateWeights) {
		if (names.size() != 1)
			MTHROW_AND_ERR("Cannot generate weights using a multi-feature generator (type %d generates %d features)\n", generator_type, (int)names.size())
		else
			p_data[0] = &(features.weights[0]);
	}

	for (string &name : names) {
		if (algorithm::ends_with(name, "Current_Smoker"))
			_p_data[SMX_UNIFIED_CURRENT_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Ex_Smoker"))
			_p_data[SMX_UNIFIED_EX_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smok_Days_Since_Quitting"))
			_p_data[SMX_UNIFIED_DAYS_SINCE_QUITTING] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smok_Pack_Years_Max"))
			_p_data[SMX_UNIFIED_SMOK_PACK_YEARS_MAX] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smok_Pack_Years_Last"))
			_p_data[SMX_UNIFIED_SMOK_PACK_YEARS_LAST] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smok_Pack_Years"))
			_p_data[SMX_UNIFIED_SMOK_PACK_YEARS] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Never_Smoker"))
			_p_data[SMX_UNIFIED_NEVER_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Unknown_Smoker"))
			_p_data[SMX_UNIFIED_UNKNOWN_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Passive_Smoker"))
			_p_data[SMX_UNIFIED_PASSIVE_SMOKER] = &(features.data[name][0]);
		else if (algorithm::contains(name, "NLST_Criterion"))
			_p_data[UNIFIED_NLST_CRITERION] = &(features.data[name][0]);
		else if (algorithm::contains(name, "Smok_Years_Since_Quitting"))
			_p_data[SMX_UNIFIED_YEARS_SINCE_QUITTING] = &(features.data[name][0]);
		else if (algorithm::contains(name, "Smoking_Intensity"))
			_p_data[SMX_UNIFIED_SMOKING_INTENSITY] = &(features.data[name][0]);
		else if (algorithm::contains(name, "Smoking_Years"))
			_p_data[SMX_UNIFIED_SMOKING_YEARS] = &(features.data[name][0]);
		else
			MTHROW_AND_ERR("unknown feature name [%s]", name.c_str());
	}
}

SMOKING_STATUS UnifiedSmokingGenerator::val2SmokingStatus(int sigVal, int smokingStatusSid, PidRec &rec)
{
	for (size_t i = 0; i < smoke_status_luts.size(); ++i)
		if (smoke_status_luts[i][sigVal])
		{
			if (SMOKING_STATUS(i) == PASSIVE_SMOKER) // convert Passive smoker to Never Smoker (currently there is no treatment for passive smoker)
				return NEVER_SMOKER;
			return SMOKING_STATUS(i);
		}


	string sigVal_name = rec.my_base_rep->dict.name(smokingStatusSid, sigVal);
	MTHROW_AND_ERR("unknown smoking status name [%s]", sigVal_name.c_str());
}

void UnifiedSmokingGenerator::genSmokingVec(PidDynamicRec& rec, UniversalSigVec &smokingStatusUsv, vector<pair<SMOKING_STATUS, int>> &smokingStatusVec, int testDate, int &unknownSmoker, int &neverSmoker, int &passiveSmoker, int &formerSmoker, int &currentSmoker)
{
	SMOKING_STATUS inVal, outVal = UNKNOWN_SMOKER;
	if (smokingStatusUsv.len > 0)
	{
		if (smokingStatusUsv.Time(0) <= testDate)
		{
			outVal = NEVER_SMOKER;
			unknownSmoker = 0;
			neverSmoker = 1;
			passiveSmoker = 0;
			formerSmoker = 0;
			currentSmoker = 0;
		}
	}

	// Calculate smoking status   

	string sigVal;
	int timeInd = 0;
	while (timeInd < smokingStatusUsv.len)
	{
		int currTime = smokingStatusUsv.Time(timeInd);
		if (currTime > testDate) { break; }
		int sigVal_num = (int)smokingStatusUsv.Val(timeInd);
		//cout << "Curr time: " << testDate << ". Sig time: " << smokingStatusUsv.Time(timeInd) << ". Sig Val: " << sigVal << endl;
		// If has Quit or Current, must be Smoker or Ex smoker. will be set later according to last value
		inVal = val2SmokingStatus(sigVal_num, smoke_status_sec_id, rec);
		timeInd++;

		// If there are several indications per time, pick maximal
		while (timeInd < smokingStatusUsv.len)
		{
			if (currTime == smokingStatusUsv.Time(timeInd))
			{
				inVal = max(inVal, val2SmokingStatus((int)smokingStatusUsv.Val(timeInd), smoke_status_sec_id, rec));
				timeInd++;
			}
			else
				break;
		}

		if ((inVal == EX_SMOKER) | (inVal == CURRENT_SMOKER))
		{
			outVal = EX_SMOKER;
			neverSmoker = 0;
			formerSmoker = 1;
		}
		if (inVal == CURRENT_SMOKER)
		{
			smokingStatusVec.push_back({ CURRENT_SMOKER, currTime });
			continue;
		}

		smokingStatusVec.push_back({ outVal, currTime });
	}

	// if last value is "Current" then this is a current smoker.
	if (outVal == CURRENT_SMOKER)
	{
		formerSmoker = 0;
		currentSmoker = 1;
	}

}

void UnifiedSmokingGenerator::genFirstLastSmokingDates(PidDynamicRec& rec, UniversalSigVec &smokingStatusUsv, UniversalSigVec &quitTimeUsv, int testDate, map<SMOKING_STATUS, pair<int, int>> &smokingStatusDates, vector<int> &dates, int birth_date)
{
	for (int timeInd = 0; timeInd < smokingStatusUsv.len; timeInd++)
	{
		int currTime = smokingStatusUsv.Time(timeInd);
		if (currTime > testDate) { break; }

		SMOKING_STATUS inVal = val2SmokingStatus((int)smokingStatusUsv.Val(timeInd), smoke_status_sec_id, rec);

		// convert Passive to Never Smoker  
		if (inVal == UNKNOWN_SMOKER)
			continue; //skip unknown smoker - later assume it's only current,ex or never and it can cause problems when patient has only unknown status
		dates.push_back(currTime);

		if (smokingStatusDates[inVal].first == NA_SMOKING_DATE)
			smokingStatusDates[inVal].first = smokingStatusUsv.Time(timeInd);
		smokingStatusDates[inVal].second = smokingStatusUsv.Time(timeInd);

	}

	// Add values according to quittime - Assume that a day before quittime, the patient was current smoker, and after that he was Former.
	// Get last quit time before test date.
	unordered_set<int> prevQuitTimes = {};
	for (int timeInd = 0; timeInd < quitTimeUsv.len; timeInd++)
	{
		int currTime = quitTimeUsv.Time(timeInd);
		if (currTime > testDate) { break; }
		int quitTime = (int)quitTimeUsv.Val(timeInd);
		if (quitTime < MedTime::MIN_DATE_SUPPORT || quitTime < birth_date) { continue; }
		if (prevQuitTimes.find(quitTime) != prevQuitTimes.end())
			continue;
		prevQuitTimes.insert(quitTime);
		int currSmokingtime = (int)med_time_converter.add_subtruct_days(quitTime, -1);
		int exSmokingTime = quitTime;
		dates.push_back(currSmokingtime);
		dates.push_back(exSmokingTime);
		if ((smokingStatusDates[CURRENT_SMOKER].first == NA_SMOKING_DATE) || (currSmokingtime < smokingStatusDates[CURRENT_SMOKER].first))
			smokingStatusDates[CURRENT_SMOKER].first = currSmokingtime;
		if ((smokingStatusDates[CURRENT_SMOKER].second == NA_SMOKING_DATE) || (currSmokingtime > smokingStatusDates[CURRENT_SMOKER].second))
			smokingStatusDates[CURRENT_SMOKER].second = currSmokingtime;
		if ((smokingStatusDates[EX_SMOKER].first == NA_SMOKING_DATE) || (exSmokingTime < smokingStatusDates[EX_SMOKER].first))
			smokingStatusDates[EX_SMOKER].first = exSmokingTime;
		if ((smokingStatusDates[EX_SMOKER].second == NA_SMOKING_DATE) || (exSmokingTime > smokingStatusDates[EX_SMOKER].second))
			smokingStatusDates[EX_SMOKER].second = exSmokingTime;
	}

	// set Ex or Never Smoking Values:
	if (smokingStatusDates[NEVER_OR_EX_SMOKER].first != NA_SMOKING_DATE)
	{
		if ((smokingStatusDates[CURRENT_SMOKER].first == NA_SMOKING_DATE) && (smokingStatusDates[EX_SMOKER].first == NA_SMOKING_DATE))
		{
			// not ex or current
			smokingStatusDates[NEVER_SMOKER].first = smokingStatusDates[NEVER_SMOKER].first != NA_SMOKING_DATE ? min(smokingStatusDates[NEVER_OR_EX_SMOKER].first, smokingStatusDates[NEVER_SMOKER].first) : smokingStatusDates[NEVER_OR_EX_SMOKER].first;
			smokingStatusDates[NEVER_SMOKER].second = smokingStatusDates[NEVER_SMOKER].first != NA_SMOKING_DATE ? max(smokingStatusDates[NEVER_OR_EX_SMOKER].second, smokingStatusDates[NEVER_SMOKER].second) : smokingStatusDates[NEVER_OR_EX_SMOKER].second;
			//cout << smokingStatusDates[NEVER_SMOKER].first << " " << smokingStatusDates[NEVER_SMOKER].second << endl;
		}
		else
		{
			// There is a current or ex status - meaning it's prob. Ex
			smokingStatusDates[EX_SMOKER].first = smokingStatusDates[EX_SMOKER].first != NA_SMOKING_DATE ? min(smokingStatusDates[NEVER_OR_EX_SMOKER].first, smokingStatusDates[EX_SMOKER].first) : smokingStatusDates[NEVER_OR_EX_SMOKER].first;
			smokingStatusDates[EX_SMOKER].second = smokingStatusDates[EX_SMOKER].first != NA_SMOKING_DATE ? max(smokingStatusDates[NEVER_OR_EX_SMOKER].second, smokingStatusDates[EX_SMOKER].second) : smokingStatusDates[NEVER_OR_EX_SMOKER].second;
		}
	}

	std::sort(dates.begin(), dates.end());
}

void UnifiedSmokingGenerator::genSmokingStatus(map<SMOKING_STATUS, pair<int, int>> &smokingStatusDates, vector<int> &dates, int testDate, int birthDate, vector<pair<SMOKING_STATUS, int>> &smokingStatusVec)
{
	// build tree:
	for (int currDate : dates)
	{
		if (currDate > testDate)
			break;

		// later we will loop on it
		SMOKING_STATUS outStatus;
		if ((smokingStatusDates[EX_SMOKER].first != NA_SMOKING_DATE) || (smokingStatusDates[NEVER_SMOKER].first != NA_SMOKING_DATE)) // Ex or Never smokers exist - node 2
			if ((smokingStatusDates[CURRENT_SMOKER].first <= currDate) && (smokingStatusDates[CURRENT_SMOKER].first != NA_SMOKING_DATE)) // First Smoker <= current Date -  node 4
				if ((smokingStatusDates[CURRENT_SMOKER].second <= currDate)) // node 21 if ex and current date together  - currently prefex ex. current is more permissivie
					if ((smokingStatusDates[NEVER_SMOKER].first <= currDate) && (smokingStatusDates[NEVER_SMOKER].first != NA_SMOKING_DATE))// node 8 
						if (smokingStatusDates[NEVER_SMOKER].second <= smokingStatusDates[CURRENT_SMOKER].second) // node14 
							if ((smokingStatusDates[EX_SMOKER].first <= currDate) && (smokingStatusDates[EX_SMOKER].first != NA_SMOKING_DATE)) // node 16
								if (smokingStatusDates[EX_SMOKER].second < smokingStatusDates[CURRENT_SMOKER].second) // node 18 - needs to reconsider
									outStatus = CURRENT_SMOKER;
								else
									outStatus = EX_SMOKER;
							else // node 16
								outStatus = CURRENT_SMOKER;
						else // node 14 
							outStatus = EX_SMOKER;
					else // node 8
						if ((smokingStatusDates[EX_SMOKER].first <= currDate) && (smokingStatusDates[EX_SMOKER].second != NA_SMOKING_DATE)) //node 9
							if ((smokingStatusDates[EX_SMOKER].second <= smokingStatusDates[CURRENT_SMOKER].second) && (smokingStatusDates[CURRENT_SMOKER].second != NA_SMOKING_DATE))//node 11
								outStatus = CURRENT_SMOKER;
							else  //node 11
								outStatus = EX_SMOKER;
						else //node 9
							outStatus = CURRENT_SMOKER;
				else // node 21 
					outStatus = CURRENT_SMOKER;
			else // node 4
				if ((smokingStatusDates[EX_SMOKER].first <= currDate) && (smokingStatusDates[EX_SMOKER].first != NA_SMOKING_DATE))//node 5
					outStatus = EX_SMOKER;
				else
					outStatus = NEVER_SMOKER;
		else // node 0
			outStatus = CURRENT_SMOKER;

		if (currDate < birthDate)
			currDate = birthDate;

		smokingStatusVec.push_back({ outStatus, currDate });
	}
}

void UnifiedSmokingGenerator::genSmokingRanges(vector<pair<SMOKING_STATUS, int>> &smokingStatusVec, int testDate, int birthDate, vector<RangeStatus> &smokeRanges)
{
	int groupStartDate;
	SMOKING_STATUS groupStatus;
	RangeStatus group;
	int estimatedSmokingStart = (birthDate / 10000 + AGE_AT_START_SMOKING) * 10000 + 101;
	int prevDate = birthDate;

	groupStatus = smokingStatusVec[0].first;
	groupStartDate = smokingStatusVec[0].second;

	for (int k = 1; k < smokingStatusVec.size(); ++k)
	{
		if (birthDate < MedTime::MIN_DATE_SUPPORT)
			break;
		SMOKING_STATUS currStatus = smokingStatusVec[k].first;
		int currDate = smokingStatusVec[k].second;
		if (currDate < MedTime::MIN_DATE_SUPPORT) { continue; }
		if (currStatus != groupStatus)
		{
			if (((groupStatus == NEVER_SMOKER) || (groupStatus == UNKNOWN_SMOKER)) && (currStatus == CURRENT_SMOKER))
			{
				int startSmokeTime = max(estimatedSmokingStart, groupStatus == NEVER_SMOKER ? groupStartDate : estimatedSmokingStart);
				startSmokeTime = min(startSmokeTime, currDate);
				group.startDate = groupStartDate;
				group.endDate = (int)med_time_converter.add_subtruct_days(startSmokeTime, -1);
				group.smokingStatus = groupStatus;
				groupStatus = currStatus;
				groupStartDate = startSmokeTime;
			}
			else if (((groupStatus == CURRENT_SMOKER) && (currStatus == EX_SMOKER)) || ((groupStatus == EX_SMOKER) && (currStatus == CURRENT_SMOKER)))
			{
				int diff = med_time_converter.diff_times(currDate, prevDate, MedTime::Date, MedTime::Days) / 2;
				group.startDate = groupStartDate;
				group.endDate = (int)med_time_converter.add_subtruct_days(prevDate, diff - 1);
				group.smokingStatus = groupStatus;

				groupStartDate = (int)med_time_converter.add_subtruct_days(prevDate, diff);
				groupStatus = currStatus;
			}
			else if (((groupStatus == NEVER_SMOKER) || (groupStatus == UNKNOWN_SMOKER)) && (currStatus == EX_SMOKER))
			{
				int ageAtEx = (int)((float)med_time_converter.diff_times(currDate, birthDate, MedTime::Date, MedTime::Days) / 365.0);
				if ((timeSinceQuittingModelSlope == missing_val) || (timeSinceQuittingModelConst == missing_val))
					MTHROW_AND_ERR("UnifiedSmokingGenerator : Ex smokers estimating Model is unknown");

				int estQuitTimeSinceQuitting = timeSinceQuittingModelSlope * ageAtEx + timeSinceQuittingModelConst;
				int max_date_back = med_time_converter.diff_times(currDate, estimatedSmokingStart, MedTime::Date, MedTime::Days);
				if (estQuitTimeSinceQuitting > max_date_back)
					estQuitTimeSinceQuitting = max_date_back;
				if (estQuitTimeSinceQuitting < 0) estQuitTimeSinceQuitting = 0;
				int quitTime = (int)med_time_converter.add_subtruct_days(currDate, -1 * estQuitTimeSinceQuitting);
				estQuitTimeSinceQuitting = max(estQuitTimeSinceQuitting, groupStatus == NEVER_SMOKER ? prevDate : estQuitTimeSinceQuitting); // if never cannot happen before prev state which is unknown or never

				int startSmokeTime = max(estimatedSmokingStart, groupStatus == NEVER_SMOKER ? prevDate : estimatedSmokingStart);
				startSmokeTime = min(startSmokeTime, currDate);
				quitTime = max(quitTime, startSmokeTime);

				//close  unknown/never smoking group
				group.endDate = startSmokeTime;
				group.smokingStatus = groupStatus;
				group.startDate = groupStartDate;
				smokeRanges.push_back(group);

				// close current smoking group
				group.startDate = startSmokeTime;
				group.smokingStatus = CURRENT_SMOKER;
				group.endDate = quitTime;

				// set Ex smoking group
				groupStatus = currStatus;
				groupStartDate = quitTime;
			}
			else if ((groupStatus == UNKNOWN_SMOKER) && (currStatus == NEVER_SMOKER))
			{
				group.startDate = groupStartDate;
				group.endDate = max(group.startDate, (int)med_time_converter.add_subtruct_days(currDate, -1));
				group.smokingStatus = groupStatus;
				groupStatus = currStatus;
				groupStartDate = currDate;
			}
			else
			{
				MTHROW_AND_ERR("UnifiedSmokingGenerator : Tranisition not covered %d -> %d", groupStatus, currStatus);
			}
			smokeRanges.push_back(group);
		}
		prevDate = currDate;
	}
	group.startDate = groupStartDate;
	group.endDate = testDate;
	group.smokingStatus = groupStatus;
	smokeRanges.push_back(group);
}

void UnifiedSmokingGenerator::genLastStatus(vector<RangeStatus> &smokeRanges, int &unknownSmoker, int &neverSmoker, int &formerSmoker, int &currentSmoker, int &passiveSmoker)
{
	unknownSmoker = neverSmoker = formerSmoker = currentSmoker = passiveSmoker = 0;
	SMOKING_STATUS lastStatus = smokeRanges[smokeRanges.size() - 1].smokingStatus;
	switch (lastStatus) {
	case UNKNOWN_SMOKER:
		unknownSmoker = 1;
		break;
	case NEVER_SMOKER:
		neverSmoker = 1;
		break;
	case EX_SMOKER:
		formerSmoker = 1;
		break;
	case CURRENT_SMOKER:
		currentSmoker = 1;
		break;
	default:
		MTHROW_AND_ERR("UnifiedSmokingGenerator : lastStatus not covered");
	}
}

void UnifiedSmokingGenerator::calcQuitTime(vector<RangeStatus> &smokeRanges, int formerSmoker, int neverSmoker, int currentSmoker, int testDate, int birthDate, float &daysSinceQuitting, float &yearsSinceQuitting)
{
	////////////////////////////////////////
	// Determine Quittime
	// find the last current -> ex transition, pick the middle of it - see genSmokingRanges 
	int quitTime = (int)missing_val;
	if (formerSmoker)
	{
		quitTime = smokeRanges.back().startDate;
	}
	else if (neverSmoker)
		quitTime = birthDate;
	else if (currentSmoker)
		quitTime = testDate;

	// Calculate time since quitting
	if (quitTime != (int)missing_val)
	{
		//MLOG("quit time: %d %d \n", testDate, (int)quitTime);
		daysSinceQuitting = (float)med_time_converter.diff_times(testDate, quitTime, MedTime::Date, MedTime::Days);
		yearsSinceQuitting = daysSinceQuitting / 365.0;
	}
	else
	{
		daysSinceQuitting = missing_val;
		yearsSinceQuitting = missing_val;
	}
}

void UnifiedSmokingGenerator::calcSmokingIntensity(UniversalSigVec &SmokingIntensityUsv, int testDate, int neverSmoker, float &smokingIntensity)
{
	////////////////////////////////////////////
	// Calc intensity

	vector<float> validSmokingIntensityValues = {};
	for (int timeInd = 0; timeInd < SmokingIntensityUsv.len; timeInd++)
	{
		if (SmokingIntensityUsv.Time(timeInd) > testDate) { break; }
		if ((SmokingIntensityUsv.Val(timeInd) > MAX_INTENSITY_TO_REMOVE) || (SmokingIntensityUsv.Val(timeInd) <= 0))
			continue;

		float value = SmokingIntensityUsv.Val(timeInd);
		value = min((float)MAX_INTENSITY_TO_TRIM, value);
		validSmokingIntensityValues.push_back(value);
	}
	if (validSmokingIntensityValues.size() > 0)
		smokingIntensity = medial::stats::mean_without_cleaning(validSmokingIntensityValues);

	if (neverSmoker)
		smokingIntensity = 0;
}

void UnifiedSmokingGenerator::getQuitAge(PidDynamicRec& rec, int lastDate, float &ageAtEx, float &deltaTime)
{
	deltaTime = ageAtEx = missing_val;
	UniversalSigVec smokingStatusUsv, quitTimeUsv, bdateUsv;

	// get signals:
	rec.uget(smoking_status_id, 0, smokingStatusUsv);
	if (smoking_quit_date_id > 0)
		rec.uget(smoking_quit_date_id, 0, quitTimeUsv);
	rec.uget(bdate_sid, 0, bdateUsv);
	int birthDate = bdateUsv.Val(0);

	// Generate First and Last Dates and dates vector
	// Map between status to <First Date, Last Date>
	map<SMOKING_STATUS, pair<int, int>> smokingStatusDates = { { NEVER_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } } ,{ PASSIVE_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } },{ EX_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } },{ CURRENT_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } },{ NEVER_OR_EX_SMOKER,{ NA_SMOKING_DATE,NA_SMOKING_DATE } } };
	vector<int> dates = {};
	genFirstLastSmokingDates(rec, smokingStatusUsv, quitTimeUsv, lastDate, smokingStatusDates, dates, birthDate);

	// Determine Smoking Status per date
	vector<pair<SMOKING_STATUS, int>> smokingStatusVec = { { UNKNOWN_SMOKER,  (int)bdateUsv.Val(0) } };

	//genSmokingVec(rec, smokingStatusUsv, smokingStatusVec, testDate, unknownSmoker, neverSmoker, passiveSmoker, formerSmoker, currentSmoker); KP logic
	genSmokingStatus(smokingStatusDates, dates, lastDate, birthDate, smokingStatusVec);

	// calc age of ex, and quit time (middle of curret and ex)
	if (smokingStatusVec.size() > 1)
	{
		for (int i = 0; i < smokingStatusVec.size() - 1; i++)
		{
			if ((smokingStatusVec[i].first == CURRENT_SMOKER) && (smokingStatusVec[i + 1].first == EX_SMOKER))
			{
				float delta = (float)med_time_converter.diff_times(smokingStatusVec[i + 1].second, smokingStatusVec[i].second, MedTime::Date, MedTime::Days) / 2;
				if (delta != 0.5)
				{
					deltaTime = delta;
					ageAtEx = ((float)med_time_converter.diff_times(smokingStatusVec[i + 1].second, birthDate, MedTime::Date, MedTime::Days) / 365.0);
				}

				else
					// check next status:
					if (i + 2 < smokingStatusVec.size())
						if (smokingStatusVec[i + 2].first == EX_SMOKER)
						{
							deltaTime = (float)med_time_converter.diff_times(smokingStatusVec[i + 2].second, smokingStatusVec[i].second, MedTime::Date, MedTime::Days);
							ageAtEx = ((float)med_time_converter.diff_times(smokingStatusVec[i + 2].second, birthDate, MedTime::Date, MedTime::Days) / 365.0);
						}
				break;
			}
		}
	}
}

void UnifiedSmokingGenerator::get_required_signal_categories(unordered_map<string, vector<string>> &signal_categories_in_use) const {
	vector<string> &in_use = signal_categories_in_use["Smoking_Status"];
	//IMPORTANT - FETCH the names in the same order as SMOKING_STATUS enum - will use it to map values
	vector<string> status = { "Unknown" ,"Never", "Passive", "Former", "Current", "Never_or_Former" };
	in_use.insert(in_use.end(), status.begin(), status.end());
}

void UnifiedSmokingGenerator::set_signal_ids(MedSignals& sigs) {
	smoking_quit_date_id = sigs.sid("Smoking_Quit_Date");
	smoking_status_id = sigs.sid("Smoking_Status");
	if (smoking_status_id < 0)
		MTHROW_AND_ERR("Error in UnifiedSmokingGenerator::set_signal_ids - repository must have Smoking_Status\n");
	smoking_intensity_id = sigs.sid("Smoking_Intensity");
	smoking_duration_id = sigs.sid("Smoking_Duration");
	smoking_pack_years_id = sigs.sid("Pack_Years");
	bdate_sid = sigs.sid("BDATE");
	if (bdate_sid < 0)
		MTHROW_AND_ERR("Error in UnifiedSmokingGenerator::set_signal_ids - repository must have BDATE\n");

}

void UnifiedSmokingGenerator::fit_for_repository(MedPidRepository &rep) {
	vector<bool> exists_sigs(4); //Smoking_Quit_Date, Smoking_Intensity, Smoking_Duration, Pack_Years
	exists_sigs[0] = rep.sigs.sid("Smoking_Quit_Date") > 0;
	exists_sigs[1] = rep.sigs.sid("Smoking_Intensity") > 0;
	exists_sigs[2] = rep.sigs.sid("Smoking_Duration") > 0;
	exists_sigs[3] = rep.sigs.sid("Pack_Years") > 0;

	req_signals.clear();
	req_signals.push_back("BDATE");
	req_signals.push_back("Smoking_Status");
	if (exists_sigs[0])
		req_signals.push_back("Smoking_Quit_Date");
	else
		MWARN("WARN:: UnifiedSmokingGenerator - repository doesn't have Smoking_Quit_Date\n");
	if (exists_sigs[1])
		req_signals.push_back("Smoking_Intensity");
	else
		MWARN("WARN:: UnifiedSmokingGenerator - repository doesn't have Smoking_Intensity\n");
	if (exists_sigs[2])
		req_signals.push_back("Smoking_Duration");
	else
		MWARN("WARN:: UnifiedSmokingGenerator - repository doesn't have Smoking_Duration\n");
	if (exists_sigs[3])
		req_signals.push_back("Pack_Years");
	else
		MWARN("WARN:: UnifiedSmokingGenerator - repository doesn't have Pack_Years\n");
}

int UnifiedSmokingGenerator::_learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors) {
	// preparing records and features for threading
	int N_tot_threads = omp_get_max_threads();
	//	MLOG("MedModel::learn/apply() : feature generation with %d threads\n", N_tot_threads);
	vector<PidDynamicRec> idRec(N_tot_threads);
	vector<int> req_signals_ids;
	for (string signal : req_signals)
		req_signals_ids.push_back(rep.dict.id(signal));
	vector<float> ageAtEx(samples.idSamples.size());
	vector<float> deltaTime(samples.idSamples.size());
	vector<float> pids(samples.idSamples.size());
#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < samples.idSamples.size(); j++) {

		int n_th = omp_get_thread_num();

		// Generate DynamicRec with all relevant signals
		int recSize = 1;
		if (idRec[n_th].init_from_rep(std::addressof(rep), samples.idSamples[j].id, req_signals_ids, recSize) < 0)
			MTHROW_AND_ERR("UnifiedSmokingGenerator : init Failed");

		//No need to  apply rep-processing
		int lastDate = (samples.idSamples[j].samples.end() - 1)->time;
		getQuitAge(idRec[n_th], lastDate, ageAtEx[j], deltaTime[j]);
		pids[j] = idRec[n_th].pid;
	}
	// clean missing values:
	vector<float> ageAtExCleaned, deltaTimeCleaned = {};

	for (int i = 0; i < ageAtEx.size(); i++)
		if (ageAtEx[i] != missing_val)
		{
			ageAtExCleaned.push_back(ageAtEx[i]);
			deltaTimeCleaned.push_back(deltaTime[i]);
		}

	MedLM lm;
	lm.init_from_string("rfactor=1;niter=10000;eiter=1e-4");
	if (!ageAtExCleaned.empty()) {
		lm.Learn(&ageAtExCleaned[0], &deltaTimeCleaned[0], (int)ageAtExCleaned.size(), 1);
		timeSinceQuittingModelSlope = lm.b[0];
		timeSinceQuittingModelConst = lm.b0;
	}
	else {
		MWARN("WARNING, no samples for learning timeSinceQuit. Will fail in apply unless you will pass default params in timeSinceQuittingModelSlope,timeSinceQuittingModelConst in json init");
	}
	return 0;
}

void UnifiedSmokingGenerator::calcSmokingDuration(int neverSmoker, int unknownSmoker, vector<RangeStatus> &smokeRanges, int birthDate, int lastPackYearsDate, UniversalSigVec &SmokingDurationUsv, int testDate, float &smokingDurationBeforeLastPackYears, float &smokingDuration)
{
	if (neverSmoker)
		smokingDuration = 0;
	else if (!unknownSmoker)
	{
		int lastDurationDateBeforTestDate;
		float lastDurationValueBeforeTestDate;
		getLastSmokingDuration(birthDate, SmokingDurationUsv, testDate, lastDurationDateBeforTestDate, lastDurationValueBeforeTestDate);
		float smokingDurationDays = lastDurationValueBeforeTestDate * 365.0;

		float smokingDurationBeforeLastPackDays = missing_val;
		if (lastPackYearsDate != missing_val)
		{
			int lastDurationDateBeforePackYears;
			float lastDurationValueBeforePackYears;
			getLastSmokingDuration(birthDate, SmokingDurationUsv, lastPackYearsDate, lastDurationDateBeforePackYears, lastDurationValueBeforePackYears);
			smokingDurationBeforeLastPackDays = lastDurationValueBeforePackYears * 365.0;
		}

		for (int k = 0; k < smokeRanges.size(); k++)
		{
			int startDate = smokeRanges[k].startDate;
			int endDate = smokeRanges[k].endDate;
			SMOKING_STATUS currStatus = smokeRanges[k].smokingStatus;

			if (currStatus == CURRENT_SMOKER)
			{
				int addedSmokingTimeFoDuration = (int)med_time_converter.diff_times(endDate, max(lastDurationDateBeforTestDate, startDate), MedTime::Date, MedTime::Days);
				smokingDurationDays += addedSmokingTimeFoDuration > 0 ? addedSmokingTimeFoDuration : 0;

				if (lastPackYearsDate != missing_val)
				{

					int addedSmokingTimeForPack = (int)med_time_converter.diff_times(min(endDate, lastPackYearsDate), max(lastDurationDateBeforTestDate, startDate), MedTime::Date, MedTime::Days);
					smokingDurationBeforeLastPackDays += addedSmokingTimeForPack > 0 ? addedSmokingTimeForPack : 0;
				}
			}
		}

		if (smokingDurationBeforeLastPackDays != missing_val)
		{
			smokingDurationBeforeLastPackYears = smokingDurationBeforeLastPackDays / 365.0;
		}
		smokingDuration = smokingDurationDays / 365.0;
	}
}

void UnifiedSmokingGenerator::getLastSmokingDuration(int birthDate, UniversalSigVec &SmokingDurationUsv, int testDate, int &lastDurationDate, float &lastDurationValue)
{
	// Keep last value of duration (but the date should  be the first that it was found  as it is usually duplicated). 
	lastDurationDate = birthDate;
	lastDurationValue = 0;

	for (int k = 0; k < SmokingDurationUsv.len; k++)
	{
		if (SmokingDurationUsv.Time(k) > testDate)
			break;
		if (SmokingDurationUsv.Val(k) != lastDurationValue)
		{
			lastDurationDate = SmokingDurationUsv.Time(k);
			lastDurationValue = SmokingDurationUsv.Val(k);
		}
	}
}

void UnifiedSmokingGenerator::calcPackYears(UniversalSigVec &SmokingPackYearsUsv, int testDate, int &neverSmoker, int &currentSmoker, int &formerSmoker, int &lastPackYearsDate, float &lastPackYears, float &maxPackYears)
{
	// Pack Years
	maxPackYears = 0;
	lastPackYearsDate = missing_val;

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

			neverSmoker = 0;
			if (currentSmoker == 0) {
				formerSmoker = 1;
			}
			lastPackYearsDate = SmokingPackYearsUsv.Time(timeInd);
		}
	}
	// This means that there wasn't any value.
	if (lastPackYears == missing_val) { maxPackYears = missing_val; }

	if (neverSmoker == 1)
	{
		maxPackYears = 0;
		lastPackYears = 0;
	}
}

void UnifiedSmokingGenerator::fixPackYearsSmokingIntensity(float smokingDurationBeforePackYears, float &smokingIntensity, float smokingDuration, float &lastPackYears, float &maxPackYears)
{
	if (lastPackYears == missing_val)
	{
		if ((smokingDuration != missing_val) && (smokingIntensity != missing_val))
			lastPackYears = maxPackYears = smokingIntensity / PACK_SIZE * smokingDuration;
	}
	else
	{
		if ((smokingDurationBeforePackYears != missing_val) && (smokingDurationBeforePackYears != 0))
		{
			if (smokingIntensity == missing_val)
			{
				if (smokingDurationBeforePackYears != 0)
				{
					smokingIntensity = maxPackYears / smokingDurationBeforePackYears * PACK_SIZE;
				}
			}
			if (smokingIntensity != missing_val)
			{
				float smokingDurationSinceLastPackYears = smokingDuration - smokingDurationBeforePackYears;
				if (smokingDurationSinceLastPackYears > 0)
				{
					lastPackYears += smokingIntensity / PACK_SIZE * smokingDurationSinceLastPackYears;
					maxPackYears += smokingIntensity / PACK_SIZE * smokingDurationSinceLastPackYears;
				}
			}
		}
	}
}

void UnifiedSmokingGenerator::printDebug(vector<RangeStatus> &smokeRanges, int qa_print, UniversalSigVec &smokingStatusUsv,
	UniversalSigVec & SmokingIntensityUsv, int birthDate, int testDate, vector<pair<SMOKING_STATUS, int>> & smokingStatusVec, PidDynamicRec& rec,
	UniversalSigVec &quitTimeUsv, UniversalSigVec &SmokingPackYearsUsv, float smokingIntensity, float smokingDuration, float yearsSinceQuitting, float maxPackYears)
{
	// debug
	vector <SMOKING_STATUS> currStatsVec = {};
	for (auto currRange : smokeRanges)
		currStatsVec.push_back(currRange.smokingStatus);
	int oldLen = (int)possibleCombinations.size();

	//possibleCombinations.insert(currStatsVec);
	if (possibleCombinations.size() != oldLen)
	{
		qa_print = 0;
		char *smokingStatusDesc[] = { "UNKNOWN_SMOKER", "NEVER_SMOKER", "PASSIVE_SMOKER", "EX_SMOKER", "CURRENT_SMOKER" };
		for (auto stat : currStatsVec)
			cout << smokingStatusDesc[stat] << " ";
		cout << endl;
	}

	if (debug_file != "")
#pragma omp critical
	{
		fprintf(fp, "********** Unified Smoking  pid: %i, birthdate %i, date %i  ************\n", rec.pid, birthDate, testDate);
		char *smokingStatusDesc[] = { "UNKNOWN_SMOKER", "NEVER_SMOKER", "PASSIVE_SMOKER", "EX_SMOKER", "CURRENT_SMOKER" };
		//fprintf(fp, "%d\t%d\t\n", rec.pid, testDate);
		fprintf(fp, "Smoking Status\t");
		for (int timeInd = 0; timeInd < smokingStatusUsv.len; timeInd++)
		{
			int currTime = smokingStatusUsv.Time(timeInd);
			if (currTime > testDate) { break; }
			string val = rec.my_base_rep->dict.name(smoke_status_sec_id, (int)smokingStatusUsv.Val(timeInd));
			fprintf(fp, "%d %s\t", currTime, val.c_str());
		}
		fprintf(fp, "\n");
		fprintf(fp, "Smoking Status Processed\t");
		for (auto currStat : smokingStatusVec)
		{
			fprintf(fp, "%d %s\t", currStat.second, smokingStatusDesc[currStat.first]);
		}
		fprintf(fp, "\n");


		fprintf(fp, "\n");
		fprintf(fp, "Smoking Intensity\t");
		for (int timeInd = 0; timeInd < SmokingIntensityUsv.len; timeInd++)
		{
			int currTime = SmokingIntensityUsv.Time(timeInd);
			if (currTime > testDate) { break; }
			float val = (float)SmokingIntensityUsv.Val(timeInd);
			fprintf(fp, "%d %f\t", currTime, val);
		}
		fprintf(fp, "\n");
		fprintf(fp, "Quit time vector\t");
		for (int timeInd = 0; timeInd < quitTimeUsv.len; timeInd++)
		{
			int currTime = quitTimeUsv.Time(timeInd);
			if (currTime > testDate) { break; }
			int val = (int)quitTimeUsv.Val(timeInd);
			fprintf(fp, "%d %d\t", currTime, val);
		}
		fprintf(fp, "\n");
		fprintf(fp, "Pack years vector\t");
		for (int timeInd = 0; timeInd < SmokingPackYearsUsv.len; timeInd++)
		{
			int currTime = SmokingPackYearsUsv.Time(timeInd);
			if (currTime > testDate) { break; }
			float val = (float)SmokingPackYearsUsv.Val(timeInd);
			fprintf(fp, "%d %f\t", currTime, val);
		}
		fprintf(fp, "\n");

		for (int kk = 0; kk < smokeRanges.size(); kk++) {
			fprintf(fp, "%i-%i %s\t", smokeRanges[kk].startDate, smokeRanges[kk].endDate, smokingStatusDesc[smokeRanges[kk].smokingStatus]);
		}
		fprintf(fp, "\n");
		fprintf(fp, "Intensity Out:\t%f\n", smokingIntensity);
		fprintf(fp, "Duration Out:\t%f\n", smokingDuration);
		fprintf(fp, "Quit time:\t%f\n", yearsSinceQuitting);
		fprintf(fp, "Pack years:\t%f\n", maxPackYears);
		//			cin.get();
	}
}

void UnifiedSmokingGenerator::addDataToMat(vector<float *> &_p_data, int index, int i, int age, int currentSmoker, int formerSmoker, float daysSinceQuitting, float daysSinceQuittingOriginal, float maxPackYears, float lastPackYears, float lastPackYearsOriginal, int neverSmoker, int unknownSmoker, int passiveSmoker, float yearsSinceQuitting, float smokingIntensity, float smokingDuration)
{
	// Add data to matrix:
	// Current_Smoker
	if (_p_data[SMX_UNIFIED_CURRENT_SMOKER] != NULL) _p_data[SMX_UNIFIED_CURRENT_SMOKER][index + i] = (float)currentSmoker;
	// Ex_Smoker
	if (_p_data[SMX_UNIFIED_EX_SMOKER] != NULL) _p_data[SMX_UNIFIED_EX_SMOKER][index + i] = (float)formerSmoker;
	// Smoke_Days_Since_Quitting
	if (_p_data[SMX_UNIFIED_DAYS_SINCE_QUITTING] != NULL) _p_data[SMX_UNIFIED_DAYS_SINCE_QUITTING][index + i] = (float)daysSinceQuitting;
	// Smok_Pack_Years_max
	if (_p_data[SMX_UNIFIED_SMOK_PACK_YEARS_MAX] != NULL) _p_data[SMX_UNIFIED_SMOK_PACK_YEARS_MAX][index + i] = (float)maxPackYears;
	// last pack years
	if (_p_data[SMX_UNIFIED_SMOK_PACK_YEARS_LAST] != NULL) _p_data[SMX_UNIFIED_SMOK_PACK_YEARS_LAST][index + i] = (float)lastPackYears;
	// Never_Smoker
	if (_p_data[SMX_UNIFIED_NEVER_SMOKER] != NULL) _p_data[SMX_UNIFIED_NEVER_SMOKER][index + i] = (float)neverSmoker;
	// Unknown_Smoker
	if (_p_data[SMX_UNIFIED_UNKNOWN_SMOKER] != NULL) _p_data[SMX_UNIFIED_UNKNOWN_SMOKER][index + i] = (float)unknownSmoker;
	// Passive_Smoker
	if (_p_data[SMX_UNIFIED_PASSIVE_SMOKER] != NULL) _p_data[SMX_UNIFIED_PASSIVE_SMOKER][index + i] = (float)passiveSmoker;
	//NLST_Criterion (only after everything was calculated, we can calc. the NLST criterion)
	// cout << "WRITING useDataComplition " << useDataComplition << " daysSinceQuitting " << daysSinceQuitting << " daysSinceQuittingOriginal " << daysSinceQuittingOriginal << endl;
	if (useDataComplition == false) {
		lastPackYears = lastPackYearsOriginal;
		daysSinceQuitting = daysSinceQuittingOriginal;
	}
	if (_p_data[UNIFIED_NLST_CRITERION] != NULL) _p_data[UNIFIED_NLST_CRITERION][index + i] = (float)calcNlst(age, unknownSmoker, daysSinceQuitting, lastPackYears);
	// Smoke_Years_Since_Quitting
	if (_p_data[SMX_UNIFIED_YEARS_SINCE_QUITTING] != NULL) _p_data[SMX_UNIFIED_YEARS_SINCE_QUITTING][index + i] = (float)yearsSinceQuitting;
	// Smok_Pack_Years
	if (_p_data[SMX_UNIFIED_SMOK_PACK_YEARS] != NULL) _p_data[SMX_UNIFIED_SMOK_PACK_YEARS][index + i] = (float)maxPackYears;
	// Smoking_Intensity
	if (_p_data[SMX_UNIFIED_SMOKING_INTENSITY] != NULL) _p_data[SMX_UNIFIED_SMOKING_INTENSITY][index + i] = (float)smokingIntensity;
	// Smoking_Years
	if (_p_data[SMX_UNIFIED_SMOKING_YEARS] != NULL) _p_data[SMX_UNIFIED_SMOKING_YEARS][index + i] = (float)smokingDuration;
}
