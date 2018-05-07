//
// This file holds calculator implementations 
// Currently for the class RepCalcSimpleSignals

#include <MedProcessTools/MedProcessTools/RepProcess.h>

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL


//.......................................................................................
float RepCalcSimpleSignals::get_age(int byear, int date)
{
	float y = (float)med_time_converter.convert_date(MedTime::Years, date);
	return (y - (float)(byear - 1900));
}

//.......................................................................................
//mode 1 (https://en.wikipedia.org/wiki/Model_for_End-Stage_Liver_Disease): 3.78ªln[serum bilirubin (mg/dL)] + 11.2ªln[INR] + 9.57ªln[serum creatinine (mg/dL)] + 6.43
//mode 2 (http://www.thecalculator.co/health/MELD-Na-Score-Calculator-846.html): MELD_mode1 ? Na ? [0.025 ª MELD_mode1 ª (140 ? Na)] + 140
float RepCalcSimpleSignals::calc_hosp_MELD(float bilirubin, float inr, float creatinine, float na, int mode)
{
	if (isAny({ bilirubin , inr, creatinine }, isMissingOrNonPositive))
		return badValue();

	if (mode == 2 && isAny({na}, isMissingOrNonPositive))
		return badValue();

	float meld1 = (float)(3.78 * log(bilirubin) + 11.2 * log(inr) + 9.57 * log(creatinine) + 6.43);
	
	if (mode == 1) {
		return meld1;
	}
	else {//no check that actually mode = 2...
		return meld1 - na - 0.025F * meld1 * (140.0F - na) + 140.0F;
	}
}

//.......................................................................................
//BMI - defined as weight / height^2. expects weight and height
float RepCalcSimpleSignals::calc_hosp_BMI(float weight, float height) {
	if (isMissingValue(weight) || isMissingValue(height) || height <= 0 || weight <= 0)
		return badValue();
	return (weight * 10000 / (height * height));
};

//.......................................................................................
//APRI - defined as ast / platelets * 2.5, expects ast, platelets.see http ://www.hepatitisc.uw.edu/page/clinical-calculators/apri
float RepCalcSimpleSignals::calc_hosp_APRI(float ast, float plt) {
	if (isMissingValue(ast) || isMissingValue(plt) || plt <= 0)
		return badValue();
	else
		return (ast * 2.5F / plt);
};

//.......................................................................................
//SIDa - defined as [Na+] - [Cl-] (mode 1) or [Na+] + [K+] - [Cl-] see Amland paper on sepsis
//expects Na, K, and Cl
float RepCalcSimpleSignals::calc_hosp_SIDA(float na, float k, float cl, int mode) {
	if (mode == 1) {
		if (isAny({ na, cl }, isMissingOrNegative))
			return badValue();
		else
			return (na - cl);
	}
	else {//no check that actually mode = 2...
		if (isAny({ na, cl, k }, isMissingOrNegative))
			return badValue();
		else
			return (na + k - cl);
	}
}
//.......................................................................................
float RepCalcSimpleSignals::calc_hosp_PaO2_FiO2_ratio(float paO2, float fiO2) {
	if (isMissingValue(paO2) || isMissingValue(fiO2) || fiO2 <= 0 || paO2 < 0)
		return badValue();
	else
		return (100.0F * paO2 / fiO2);
};

//.......................................................................................
float RepCalcSimpleSignals::calc_hosp_is_african_american(float ethnicity, float african_american_dict_id) {
	if (isMissingValue(ethnicity))
		return badValue();
	else
		//return (ethnicity == 17612 ? 1.0F : 0.0F);
		return (ethnicity == african_american_dict_id ? 1.0F : 0.0F);
};

//.......................................................................................
//return nervous SOFA
//x = Glasgow Coma Scale score
float RepCalcSimpleSignals::calc_hosp_SOFA_nervous(float x, bool optimistic) {
	if (optimistic) {
		if (isMissingValue(x)) x = 15.0;
	}
	else {
		if (isMissingValue(x)) x = 0.0;
	}

	int s = 0;
	if (x < 6) s = 4; else if (x <= 9) s = 3; else if (x <= 12) s = 2; else if (x <= 14) s = 1;
	return (float)s;
}

//.......................................................................................
//return liver SOFA
//x = serum bilirubin level mg/dl
float RepCalcSimpleSignals::calc_hosp_SOFA_liver(float x, bool optimistic) {
	if (optimistic) {
		if (isMissingValue(x)) x = 0.0;
	}
	else {
		if (isMissingValue(x)) x = 13.0;
	}

	int s = 0;
	if (x > 12.0) s = 4; else if (x >= 6.0) s = 3; else if (x >= 2.0) s = 2; else if (x >= 1.2) s = 1;
	return (float)s;
}

//.......................................................................................
//return coagulation SOFA
//x = plateletCount 10^9/L
float RepCalcSimpleSignals::calc_hosp_SOFA_coagulation(float x, bool optimistic) {
	if (optimistic) {
		if (isMissingValue(x)) x = 150.0;
	}
	else {
		if (isMissingValue(x)) x = 0.0;
	}

	int s = 0;
	if (x < 20) s = 4; else if (x < 50) s = 3; else if (x < 100) s = 2; else if (x < 150) s = 1;
	return (float)s;
}

//.......................................................................................
//functions to calculate medication per kg from existing source signals on MIMIC
//Note that some medications have only a per kg signal and some also have an absolute quantity signal,
//hence the need for combining the two sources. This is very MIMIC-specific
float RepCalcSimpleSignals::calc_hosp_dopamine_per_kg(float medPerKg) { return (isMissingOrNegative(medPerKg) ? badValue() : medPerKg); }
float RepCalcSimpleSignals::calc_hosp_epinephrine_per_kg(float medPerKg, float med, float kg) { return med2KgPlusMed(medPerKg, med, kg); }
float RepCalcSimpleSignals::calc_hosp_norepinephrine_per_kg(float medPerKg, float med, float kg) { return med2KgPlusMed(medPerKg, med, kg); }
float RepCalcSimpleSignals::calc_hosp_dobutamine_per_kg(float medPerKg) { return (isMissingOrNegative(medPerKg) ? badValue() : medPerKg); }

//.......................................................................................
//calculate SIRS score. 2 out of the following 4 - Temp > 38 or < 36; Heart Rate > 90; Repoiratory Rate > 20 per min or PaCO2 < 32 mm Hg; 
//WBC > 12 or WBC < 4 or > 10% bands
//expects Temperature, Heart_Rate, Resp_Rate, Art_PaCO2, WBC, Neutrophils_Bands. 
//mode 1: output is 1 if certainly >=2, 0 if certainly < 2, badValue() if undecided
//mode 2: raw score if all constituents present, badValue() otherwise.
//mode 3: optimistic
//mode 4: pessimistic
float RepCalcSimpleSignals::calc_hosp_SIRS(float temp, float hr, float resp, float paco2, float wbc, float bands, int mode) {
	int maxSIRS = 0, minSIRS = 0;
	// Check Temperature
	if (isMissingValue(temp)) maxSIRS++; else if (temp < 36 || temp > 38) { maxSIRS++; minSIRS++; }

	// Check Heart Rate
	if (isMissingValue(hr)) maxSIRS++; else if (hr > 90) { maxSIRS++; minSIRS++; }

	// Check Respiratory Rate and PaCO2
	bool c1 = (!isMissingValue(resp) && resp > 20);
	bool c2 = (!isMissingValue(paco2) && paco2 < 32);

	if (isMissingValue(resp) || isMissingValue(paco2) || c1 || c2)
		maxSIRS++;

	if (c1 || c2)
		minSIRS++;

	// Check WBC and bands
	c1 = (!isMissingValue(wbc) && (wbc > 12 || wbc < 4));
	c2 = (!isMissingValue(bands) && bands > 10);

	if (isMissingValue(wbc) || isMissingValue(bands) || c1 || c2)
		maxSIRS++;

	if (c1 || c2)
		minSIRS++;

	if (mode == 1) {
		if (minSIRS >= 2)
			return 1.0F;
		else if (maxSIRS < 2)
			return 0.0F;
		else
			return badValue();
	}
	else if (mode == 2) {
		if (maxSIRS == minSIRS)
			return (float)maxSIRS;
		else
			return badValue();
	}
	else if (mode == 3) {
		return (float)minSIRS;
	}
	else {//no check if actually 4...
		return (float)maxSIRS;
	}
}

//.......................................................................................

//pressure adjusted hr: HR * CVP / MAP
//expects hr, cvp, map
float RepCalcSimpleSignals::calc_hosp_pressure_adjusted_hr(float hr, float cvp, float map) {
	if (isAny({ cvp, hr, map }, isMissingValue) || map <= 0 || cvp < 0 || hr < 0)
		return badValue();
	else
		return (hr * cvp / map);
};

//.......................................................................................
//MODS score. expects paO2, FiO2, platelets, serum bilirubin, hr, cvp, map, gcs, serum creatinine
//http://reference.medscape.com/calculator/mods-score-multiple-organ-dysfunction
//see also Marshall et al. "Multiple organ dysfunction score : a reliable descriptor of a complex clinical outcome"
float RepCalcSimpleSignals::calc_hosp_MODS(float paO2, float fiO2, float plt, float bili, float hr, float cvp, float map, float gcs, float cre) {
	if (isAny({ paO2, fiO2, plt, bili, hr, cvp, map, gcs, cre }, isMissingOrNegative) || map == 0.0F || fiO2 == 0.0F)
		return badValue();

	float x;
	int s = 0;

	x = 100.0F * paO2 / fiO2; if (x < 76) s += 4; else if (x < 151) s += 3; else if (x < 226) s += 2; else if (x < 301) s += 1;
	x = plt; if (x < 21) s += 4; else if (x < 51) s += 3; else if (x < 81) s += 2; else if (x < 121) s += 1;
	x = bili; if (x > 14) s += 4; else if (x > 7) s += 3; else if (x > 3.5) s += 2; else if (x > 1.2) s += 1;
	x = hr * cvp / map; if (x > 30) s += 4; else if (x > 20) s += 3; else if (x > 15) s += 2; else if (x > 10) s += 1;
	x = gcs; if (x < 7) s += 4; else if (x < 10) s += 3; else if (x < 13) s += 2; else if (x < 15) s += 1;
	x = cre; if (x > 5.7) s += 4; else if (x > 4.0) s += 3; else if (x > 2.3) s += 2; else if (x > 1.1) s += 1;

	return (float)s;
}

//.......................................................................................
//shockIndex - defined as heart rate divided by systolic blood pressure (mode 1) or mean pressure (mode=2). 
//expects heart rate, systolic bp, diastolic bp
float RepCalcSimpleSignals::calc_hosp_shock_index(float hr, float sysBp, float diaBp, int mode) {
	if (mode == 1) {
		if (isMissingValue(hr) || isMissingValue(sysBp) || sysBp <= 0)
			return badValue();
		else
			return hr / sysBp;
	}
	else {//mode == 2, no check!
		float map = (2.0F*diaBp + sysBp) / 3.0F;
		if (isMissingValue(hr) || isMissingValue(sysBp) || isMissingValue(diaBp) || map == 0 || diaBp < 0 || sysBp < 0)
			return badValue();
		else
			return hr / map;
	}
}

//.......................................................................................
//pulse pressure: systolic bp - diastolic bp
//expects systolic bp, diastolic bp
float RepCalcSimpleSignals::calc_hosp_pulse_pressure(float sysBp, float diaBp) {
	if (isMissingValue(sysBp) || isMissingValue(diaBp) || sysBp < 0 || diaBp < 0)
		return badValue();
	else
		return (sysBp - diaBp);
};

//.......................................................................................
//feature for various types of eGFR. mode 1: MDRD, mode 2: CKD EPI, mode 3: KeGFR (not implemented)
//useEthnicity - a boolean flag
//expects: creatinine, age, gender (1-male/2-female/missing), isAfricanAmerican (0/1/missing), 
//important: KeGFR looks 72 hours back. 
float RepCalcSimpleSignals::calc_hosp_eGFR(float cr, float age, float gender, float isAfricanAmerican, bool useEthnicity, int mode) {
	if (isAny({ cr, age, gender, isAfricanAmerican }, isMissingValue))
		return badValue();

	if (cr <= 0 || age <= 0 || (gender != 1.0 && gender != 2.0) || (isAfricanAmerican != 0.0 && isAfricanAmerican != 1.0))
		return badValue();

	bool isAfAm = (isAfricanAmerican == 1.0);
	bool isFemale = (gender == 2.0);

	//double MDRD = log(175.0) - 1.154*log((double)cr) - 0.203*log(age) + (isFemale ? log(0.742) : 0.0) + ((isAfAm && useEthnicity) ? log(1.212) : 0.0);
	//kdigo says 186, not 175
	double MDRD = log(186.0) - 1.154*log((double)cr) - 0.203*log(age) + (isFemale ? -0.298406 : 0.0) + ((isAfAm && useEthnicity) ? 0.1922719 : 0.0);
	MDRD = exp(MDRD);

	if (mode == 1) {
		return (float)MDRD;
	}
	else if (mode == 2) {
		double a = (isFemale ? -0.329 : -0.411);
		double k = (isFemale ? 0.7 : 0.9);

		double ckdepi = log(141.0) + a * log(min((double)cr / k, 1.0)) - 1.209 * log(max((double)cr / k, 1.0)) + age * log(0.993) + (isFemale ? log(1.018) : 0.0) +
			((isAfAm && useEthnicity) ? log(1.159) : 0.0);

		ckdepi = exp(ckdepi);
		return (float)ckdepi;
	}

	return -1;
}

//.......................................................................................
//x - mmHg PaO2 / FiO2
//y - is mechanically ventilated, 1 if true, 0 otherwise			
//optimistic - if true, treats missing data as optimistically as possible, otherwise treats it pessimistically.
float RepCalcSimpleSignals::calc_hosp_SOFA_respiratory(float x, float y, bool optimistic) {
	if (optimistic) {
		if (isMissingValue(x)) x = 1000.0;
		if (isMissingValue(y)) y = 0.0;
	}
	else {
		if (isMissingValue(x)) x = 0.0;
		if (isMissingValue(y)) y = 1.0;
	}

	int s = 0;
	if (x < 100.0 && y > 0.0) s = 4; else if (x < 200.0 && y > 0.0) s = 3; else if (x < 300.0) s = 2; else if (x < 400.0) s = 1;
	return (float)s;
}

//.......................................................................................
//return renal SOFA
//x = serum creatinine level mg/dl
//y = urine output, ml/min
float RepCalcSimpleSignals::calc_hosp_SOFA_renal(float x, float y, bool optimistic) {
	if (optimistic) {
		if (isMissingValue(x)) x = 0.0;
		if (isMissingValue(y)) y = 1000.0;
	}
	else {
		if (isMissingValue(x)) x = 10.0;
		if (isMissingValue(y)) y = 0.0;
	}

	int s = 0;
	if (x > 5.0 || y < 200.0) s = 4; else if (x >= 3.5 || y < 500.0) s = 3; else if (x >= 2.0) s = 2; else if (x >= 1.2) s = 1;
	return (float)s;
}

//.......................................................................................
//return cardiovascular SOFA
//x - mean arterial pressure (mmHg)
//y - dopamine, vasopressor µg/kg/min
//z - epinephrine, vasopressor µg/kg/min
//u - norepinephrine, vasopressor µg/kg/min
//v - dobutamine, vasopressor µg/kg/min
float RepCalcSimpleSignals::calc_hosp_SOFA_cardio(float x, float y, float z, float u, float v, bool optimistic) {
	if (optimistic) {
		if (isMissingValue(x)) x = 100.0;
		if (isMissingValue(y)) y = 0.0;
		if (isMissingValue(z)) z = 0.0;
		if (isMissingValue(u)) u = 0.0;
		if (isMissingValue(v)) v = 0.0;
	}
	else {
		if (isMissingValue(x)) x = 0.0;
		if (isMissingValue(y)) y = 100.0;
		if (isMissingValue(z)) z = 100.0;
		if (isMissingValue(u)) u = 100.0;
		if (isMissingValue(v)) v = 100.0;
	}

	int s = 0;
	if (y > 15.0 || z > 0.1 || u > 0.1) s = 4;
	else if (y > 5.0 || z > 0.0 || u > 0.0) s = 3;
	else if (y > 0.0 || v > 0.0) s = 2;
	else if (x < 70.0) s = 1;

	return (float)s;
}

//.......................................................................................
//SOFA score: mode 1 is sum of constituents, mode 2 is max.
//it is assumed that all constituents were calculated similarly - all optimistic or all pessimistic
//as a result, there should be no missing values
float RepCalcSimpleSignals::calc_hosp_SOFA(float snerve, float sliver, float scoag, float sresp, float srenal, float scardio, int mode) {
	return (mode == 1 ? (snerve + sliver + scoag + sresp + srenal + scardio) : max(max(max(max(max(snerve,sliver),scoag),sresp),srenal),scardio));
}


//.......................................................................................
//return qSOFA
//gcs - glasgow coma score
//sBp - systolic bp
//resp - respiration
//mode 1: if at least 2 of 3 conditions are fulfilled, return 1, if not, return 0, otherwise (undecided) return badValue
//mode 2 - return the raw number of conditions met, bad value if not clear
float RepCalcSimpleSignals::calc_hosp_qSOFA(float gcs, float sBp, float resp, int mode) {
	int sMin = 0, sMax = 0;
	if (isMissingValue(gcs))
		sMax++;
	else if (gcs <= 13) {
		sMin++; sMax++;
	}

	if (isMissingValue(sBp))
		sMax++;
	else if (sBp <= 100) {
		sMin++; sMax++;
	}

	if (isMissingValue(resp))
		sMax++;
	else if (resp >= 22) {
		sMin++; sMax++;
	}

	if (mode == 1) {
		if (sMax < 2)
			return (float)0;
		else if (sMin >= 2)
			return (float)1;
		else
			return badValue();
	}
	else { //no check if actually 2...
		if (sMax == sMin)
			return (float)sMin;
		else
			return badValue();
	}
}

//.......................................................................................
/*
int RepCalcSimpleSignals::_apply_calc_hosp_time_dependent_pointwise(PidDynamicRec& rec, vector<int>& time_points,
	float(*calcFunc)(const vector<pair<int, float> >&, int, const vector<float>&)) {

	//Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	//sid for the calculated signal
	int v_sid = V_ids[0];

	//length of calculated data
	size_t len = time_points.size();

	//handle abnormal case
	if (len == 0)
		return -1;

	//sigs_ids is a vector containing the sids of the components in the right order
	size_t nComponents = sigs_ids.size();

	//rec should contain space for the required sids
	if (rec.usvs.size() < nComponents)
		rec.usvs.resize(nComponents);

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end()); iteratorSignalIds.insert(v_sid);
	versionIterator vit(rec, iteratorSignalIds);

	//auto it = calc2req_sigs.find(calculator);
	//if (it == calc2req_sigs.end()) {//unexpected
	//	cout << "cannot find name " << calculator << " in calc2req_sigs" << endl;
	//	return -1;
	//}

	for (int iver = vit.init(); iver >= 0; iver = vit.next_different()) {
		//componentData will hold the data of component signals after alignment to required time_points
		//for each time point we take the most recent data or missing_data otherwise.
		vector<vector<pair<int, float> > > componentData(nComponents);
		vector<int> curTimes; //times of data for current component
		vector<float> curVals; //(processed) vals of data for current component	

		vector<int> iverTimes(time_points.begin(), time_points.begin() + iver + 1);

		//get component data
		for (size_t i = 0; i < nComponents; ++i) {
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);
			curTimes.clear();
			curVals.clear();

			process_hosp_signal(signals[i], rec.usvs[i], curTimes, curVals);

			int curLen = (int)curTimes.size();

			vector<size_t> indices;
			find_sorted_vec_in_sorted_vec(iverTimes, curTimes, indices);

			//get the values from the correct indices
			componentData[i].reserve(len);

			if (curTimes.empty())
				componentData[i].insert(componentData[i].begin(), iver + 1, pair<int, float>(beforeEverything(), badValue()));
			else {
				for (int j = 0; j <= iver; ++j) {
					size_t curInd = indices[j];
					//curInd is the place we would add time_points[j] in curTimes; 0 means at the beginning
					if (curInd == 0)
						componentData[i].push_back(pair<int, float>(beforeEverything(), badValue()));
					else
						componentData[i].push_back(pair<int, float>(curTimes[(int)curInd - 1], curVals[(int)curInd - 1]));
				}
			}
		}

		//we now have a matrix of data in componentData. For a given j, componentData[i][j] represent values in time_points[j]
		vector<pair<int, float> > inValues(nComponents); //for a given time, the inputs
		vector<float> calcValues(iver + 1); //for all times, the calculated values

		for (int j = 0; j <= iver; ++j) {
			for (int i = 0; i < nComponents; ++i)
				inValues[i] = componentData[i][j];

			float val = calcFunc(inValues, time_points[j], coeff);
			calcValues[j] = val;
		}

		// pushing virtual data into rec (into orig version)
		rec.set_version_universal_data(v_sid, iver, &time_points[0], &calcValues[0], (int)(iver + 1));
	}

	return 0;
}
*/

//.......................................................................................
/*
int RepCalcSimpleSignals::_apply_calc_hosp_time_dependent_pointwise(PidDynamicRec& rec, vector<int>& time_points,
	float(*calcFunc)(const vector<pair<int,float> >&, int, const vector<float>&)) {

	//sid for the calculated signal
	int v_sid = V_ids[0];

	//length of calculated data
	size_t len = time_points.size();

	//handle abnormal case
	if (len == 0)
		return -1;

	//sigs_ids is a vector containing the sids of the components in the right order
	size_t nComponents = sigs_ids.size();

	//rec should contain space for the required sids
	if (rec.usvs.size() < nComponents)
		rec.usvs.resize(nComponents);

	//componentData will hold the data (actual time+val) of component signals after alignment to required time_points
	//for each time point we take the most recent data or missing_data otherwise.
	vector<vector<pair<int,float> > > componentData(nComponents);
	vector<int> curTimes; //times of data for current component
	vector<float> curVals; //(processed) vals of data for current component

	auto it = calc2req_sigs.find(V_names[0]);
	if (it == calc2req_sigs.end()) //unexpected
		return -1;

	//get component data
	for (size_t i = 0; i < nComponents; ++i) {
		rec.uget(sigs_ids[i], 0, rec.usvs[i]);
		curTimes.clear();
		curVals.clear();

		process_hosp_signal((it->second)[i], rec.usvs[i], curTimes, curVals);

		int curLen = (int)curTimes.size();
		
		//find the correct indices in the signal to extrapolate from 

		vector<size_t> indices;
		find_sorted_vec_in_sorted_vec(time_points, curTimes, indices);

		//get the values from the correct indices
		componentData[i].reserve(len);
		
		if (curTimes.empty())
			//IMPORTANT: we not only put a bad value, but time it to a minimal time, making the data maximally useless for the user
			componentData[i].insert(componentData[i].begin(), time_points.size(), pair<int,float>(beforeEverything(), badValue()));
		else {
			for (int j = 0; j < len; ++j) {
				size_t curInd = indices[j];
				//curInd is the place we would add time_points[j] in curTimes; 0 means at the beginning
				if (curInd == 0)
					componentData[i].push_back(pair<int, float>(beforeEverything(), badValue()));
				else
					componentData[i].push_back(pair<int, float>(curTimes[(int)curInd - 1], curVals[(int)curInd - 1]));
			}
		}
	}

	//we now have a matrix of data in componentData. For a given j, componentData[i][j] represent values in time_points[j]
	vector<pair<int, float> > inValues(nComponents); //for a given time, the inputs
	vector<float> calcValues(len); //for all times, the calculated values

	for (int j = 0; j < len; ++j) {
		for (int i = 0; i < nComponents; ++i)
			inValues[i] = componentData[i][j];

		float val = calcFunc(inValues, time_points[j], coeff);
		calcValues[j] = val;
	}

	// pushing virtual data into rec (into orig version)
	rec.set_version_universal_data(v_sid, 0, &time_points[0], &calcValues[0], (int)len);

	// pointing all versions to the 0 one
	for (int ver = 1; ver<rec.get_n_versions(); ver++)
		rec.point_version_to(v_sid, 0, ver);

	return 0;
}
*/
int RepCalcSimpleSignals::_apply_calc_24h_urine_output(PidDynamicRec& rec, vector<int>& time_points) {

	//Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int window = 1440; //24 hours

	//sid for the calculated signal
	int v_sid = V_ids[0];

	//length of calculated data
	size_t len = time_points.size();

	//handle abnormal case
	if (len == 0)
		return -1;

	if (rec.usvs.size() < 1)
		rec.usvs.resize(1);

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end()); iteratorSignalIds.insert(v_sid);
	versionIterator vit(rec, iteratorSignalIds);

	//auto it = calc2req_sigs.find(V_names[0]);
	//if (it == calc2req_sigs.end()) //unexpected
	//	return -1;

	for (int iver = vit.init(); iver >= 0; iver = vit.next_different()) {
		rec.uget(sigs_ids[0], iver, rec.usvs[0]);
		int origLen = rec.usvs[0].len;

		vector<int> iverTimes(time_points.begin(), time_points.begin() + iver + 1);

		//data will hold pairs of end time and val.
		vector<pair<int, float> > data; data.reserve(origLen + 1);

		//get signal data
		//if end time < start time, we ignore the point. end time = start time is ok
		//negative values are turned to 0
		for (int i = 0; i < origLen; ++i) {
			int startTime = rec.usvs[0].Time(i, 0);
			int endTime = rec.usvs[0].Time(i, 1);

			if (endTime < startTime)
				continue;

			float val = max(0.0F, rec.usvs[0].Val(i, 0));

			data.push_back(pair<int, float>(endTime, val));
		}

		//the raw signal is understood as sum total over a time period; however, "spikes" are allowed, namely, 
		//a time range of size 0. In addition, overlaps are allowed. Important: only the end times are originally provided, so only cumulative
		//quantities can be calculated.(signals with the same end time have their values summed.)
		//we don't expect the times to be sorted, so sorting is required.
		//Note that the first end time records an amount that accumulated since nobody knows when. 

		//sort by end time
		std::sort(data.begin(), data.end());

		//skip first time (<data> may also be empty, or length 1, so loop could be ignored)
		int index = -1;
		for (int i = 1; i < data.size(); ++i) {
			if (data[i].first != data[i - 1].first) {
				index = i;
				break;
			}
		}

		//handle the case where all have the same, first, time (including when data.size()<=1)
		if (index == -1) {
			vector<float> calcValues(iverTimes.size(), badValue());
			// pushing virtual data into rec (into orig version)
			rec.set_version_universal_data(v_sid, iver, &iverTimes[0], &calcValues[0], (int)(iver + 1));

			continue;
		}

		//at this point, data.size()>1, index points to the second unique end time

		//vector to hold the cumulative amount for each end time. 
		vector<float> cumVals; cumVals.reserve(data.size());
		vector<int> times; times.reserve(data.size());

		float curSum = 0; //the actual amount in the first time is immaterial

		for (int i = index; i < data.size(); ++i) {		
			if (data[i].first != data[i - 1].first) {
				cumVals.push_back(curSum);
				times.push_back(data[i - 1].first);
			}
			curSum += data[i].second;
		}

		//handle last point 
		cumVals.push_back(curSum);
		times.push_back(data.back().first);

		//to estimate the cumulative amount over a window we do the following:
		//for the second point of the window, we take the latest point and not interpolate so as not to have a leak.
		//the first point is <window> before the point we took instead of the second. 
		//if the first point is between points in <times>, we will interpolate the cumulative amount.		
		
		
		vector<size_t> indicesWindowStart, indicesWindowEnd;

		//to estimate the end value for iverTimes[i] we will use times[indicesWindowEnd[i] - 1] if it's valid
		find_sorted_vec_in_sorted_vec(iverTimes, times, indicesWindowEnd);

		//find were the beginning points of windows would be
		//if the end is before all times, then so will the beginning be, and that window will be discarded later
		//otherwise, <window> time before the last seen actual time
		vector<int> timePointsMinusWindow; timePointsMinusWindow.reserve(iverTimes.size());
		for (int i = 0; i < iverTimes.size(); ++i)
			timePointsMinusWindow.push_back(indicesWindowEnd[i] == 0 ? beforeEverything() : times[indicesWindowEnd[i] - 1] - window);

		find_sorted_vec_in_sorted_vec(timePointsMinusWindow, times, indicesWindowStart);

		//to hold the calcluated values
		vector<float> calcValues; calcValues.reserve(iverTimes.size());

		for (int i = 0; i < indicesWindowStart.size(); ++i) {
			//if the window starts before all times - bad value
			if (indicesWindowStart[i] == 0)
				calcValues.push_back(badValue());
			else {//the window starts between two points 
				//interpolate 
				float coef = (float)(timePointsMinusWindow[i] - times[indicesWindowStart[i] - 1]) / 
					(float)(times[indicesWindowStart[i]] - times[indicesWindowStart[i] - 1]);

				//val is the cumulative value at the end of the window minus the interpolated val <window> before it
				float val = cumVals[indicesWindowEnd[i] - 1] - (1.0F - coef)*cumVals[indicesWindowStart[i] - 1]
					- coef*cumVals[indicesWindowStart[i]];

				//if (!isfinite(val))
				//	cout << "timePointsMinusWindow[i] = " << timePointsMinusWindow[i] << " times[indicesWindowStart[i] - 1] = " << times[indicesWindowStart[i] - 1] << " times[indicesWindowEnd[i] - 1] = " << times[indicesWindowEnd[i] - 1] << " times[indicesWindowStart[i] - 1] = " << times[indicesWindowStart[i] - 1] << " coef = " << coef << " cumVals[indicesWindowStart[i] - 1] = " << cumVals[indicesWindowStart[i] - 1] << " cumVals[indicesWindowEnd[i] - 1] = " << cumVals[indicesWindowEnd[i] - 1] << " val = " << val << endl;
				
				calcValues.push_back(val);
			}
		}

		// pushing virtual data into rec (into orig version)
		rec.set_version_universal_data(v_sid, iver, &iverTimes[0], &calcValues[0], (int)(iver + 1));

	}

	return 0;
}

//.......................................................................................


int RepCalcSimpleSignals::_apply_calc_log(PidDynamicRec& rec, vector<int>& time_points) {

	//Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	//sid for the calculated signal
	int v_sid = V_ids[0]; // Output signal
	int i_sid = sigs_ids[0]; // Input signal

	if (rec.usvs.size() < 1) rec.usvs.resize(1);

	// Loop on versions of signal
	set<int> iteratorSignalIds = { i_sid, v_sid };
	versionIterator vit(rec, iteratorSignalIds);

	for (int iver = vit.init(); iver >= 0; iver = vit.next_different()) {

		rec.uget(i_sid, iver, rec.usvs[0]); // getting the requested input signal into usv[0]

		float epsilon = (float)0.00001;
		if (rec.usvs[0].len) {

			vector<float> v_vals;
			vector<int> v_times;

			// calculating for each creatinine point up to the relevant time-point
			for (int i = 0; i<rec.usvs[0].len; i++) {
				if (time_points.size() != 0 && rec.usvs[0].Time(i) > time_points[iver])
					break;

				if (rec.usvs[0].Val(i) > epsilon) {
					v_times.push_back(rec.usvs[0].Time(i));
					v_vals.push_back(log(rec.usvs[0].Val(i)));
				}
			}

			// pushing virtual data into rec (into orig version)
			rec.set_version_universal_data(v_sid, iver, &v_times[0], &v_vals[0], (int)v_times.size());
		}
	}


	return 0;

}


//.......................................................................................

// Given two strictly increasing vectors : 'target' and 'given', create a vector of indices of the same length of 'target',
// Such that the value of given at index[i] is the closest to the value of target[i] while requiring - 
//		1. abs(given[index[i]] - target[i]) < max_diff or max_diff < 0
//		2. given[index[i]] < max
// Insert -1 if none can be found
/*void RepCalcSimpleSignals::index_targets_in_given_vector(const vector<int> &target, const vector<int> &given, int& max_diff, int signals_time_unit, int diff_time_unit, int& max, vector<size_t>& indices, bool onlyPast) {

	indices.resize(target.size(), -1);

	size_t iGiven = 0;
	for (size_t iTarget = 0; iTarget < target.size(); iTarget++) {

		int targetTime = med_time_converter.convert_times(signals_time_unit, diff_time_unit, target[iTarget]);

		int iMinDiff = -1;
		int diff;
		while (iGiven < given.size() && given[iGiven] <= max) {
			// Passing/Hitting the target
			int givenTime = med_time_converter.convert_times(signals_time_unit, diff_time_unit, given[iGiven]);
			if (given[iGiven] >= target[iTarget]) {
				if (iGiven == 0) {
					iMinDiff = (int)iGiven;
					diff = givenTime;
				}
				else {
					int prevGivenTime = med_time_converter.convert_times(signals_time_unit, diff_time_unit, given[iGiven - 1]);
					if (givenTime - targetTime < targetTime - prevGivenTime) {
						iMinDiff = (int)iGiven;
						diff = givenTime - targetTime;
					}
					else {
						iMinDiff = (int)(--iGiven);
						diff = targetTime - prevGivenTime;
					}
				}
				break;
			}
			iGiven++;
		}

		// Have we reached the end ?
		if (iGiven > 0 && (iGiven == given.size() || given[iGiven] >= max)) {
			iMinDiff = (int)(--iGiven);
			int givenTime = med_time_converter.convert_times(signals_time_unit, diff_time_unit, given[iGiven]);
			diff = targetTime - givenTime;
		}

		// Are we ok ?
		if (iMinDiff != -1 && (max_diff < 0 || fabs((float)diff) <= max_diff))
			indices[iTarget] = iMinDiff;

	}
}
*/

// Given two strictly increasing vectors : 'target' and 'given', create a vector of indices of the same length of 'target',
// Such that the value of given at index[i] is the closest to the value of target[i] while requiring - 
//		1. abs(given[index[i]] - target[i]) < max_diff or max_diff < 0
//		2. given[index[i]] < max
// Insert -1 if none can be found
// if onlyPast is true, we require also that given[index[i]] >= target[i] 
void RepCalcSimpleSignals::index_targets_in_given_vector(const vector<int> &target, const vector<int> &given, int& max_diff, int signals_time_unit, int diff_time_unit, int& max, vector<size_t>& indices, bool onlyPast) {
	indices.resize(target.size(), -1);

	if (given.size() == 0 || target.size() == 0)
		return;

	size_t iGiven = 0;
	for (size_t iTarget = 0; iTarget < target.size(); iTarget++) {

		int targetTime = med_time_converter.convert_times(signals_time_unit, diff_time_unit, target[iTarget]);

		int iMinDiff = -1; //the best given index for the current target
		int diff; //the absolute time difference to the best given time

		//increment the given candidate, aiming to pass the target, 
		//while having a valid candidate (before max)
		while (iGiven < given.size() && given[iGiven] <= max) {
			int givenTime = med_time_converter.convert_times(signals_time_unit, diff_time_unit, given[iGiven]);
			
			// passing or hitting the target
			if (given[iGiven] >= target[iTarget]) {
				//hitting the target
				if (given[iGiven] == target[iTarget]) {
					iMinDiff = (int)iGiven;
					diff = 0;
				}
				else { //passing the target
					if (iGiven == 0) {
						if (!onlyPast) {							
							iMinDiff = (int)iGiven;
							diff = givenTime - targetTime;
						}

						//else it's in the future, and invalid (iMinDiff remains -1)
					}
					else {//iGiven > 0, there is a previous candidate before the target
						int prevGivenTime = med_time_converter.convert_times(signals_time_unit, diff_time_unit, given[iGiven - 1]);
						if (givenTime - targetTime < targetTime - prevGivenTime && !onlyPast) {//future time is closer than past time
							iMinDiff = (int)iGiven;
							diff = givenTime - targetTime;
						}
						else {//use the past time
							iMinDiff = (int)(--iGiven);
							diff = targetTime - prevGivenTime;
						}
					}
				}
				break;
			}

			iGiven++;
		}

		// at this point, either iGiven = 0, and given[0] >= max (no good given time)
		// or it's the first given that is >= max or iGiven = given.size() (should try the previous given)
		// in any case, the given did not hit or pass the target
		if (iGiven > 0 && (iGiven == given.size() || given[iGiven] >= max)) {
			iMinDiff = (int)(--iGiven);
			int givenTime = med_time_converter.convert_times(signals_time_unit, diff_time_unit, given[iGiven]);
			diff = targetTime - givenTime;
		}

		// check if the candidate (if there is one) is closer than max_diff (if defined)
		// if we don't enter the condition, the indices[iTarget]  will remain -1
		if (iMinDiff != -1 && (max_diff < 0 || fabs((float)diff) <= max_diff))
			indices[iTarget] = iMinDiff;

	}
}


//.......................................................................................

int RepCalcSimpleSignals::_apply_calc_hosp_pointwise(PidDynamicRec& rec, vector<int>& time_points,
	float(*calcFunc)(const vector<float>&, const vector<float>&), bool onlyPast) {

	//we either use a timer signal times or the union of component times
	//that don't belong to missing values
	bool useTimer = (timer_signal_id != -1);

	//Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	//sid for the calculated signal
	int v_sid = V_ids[0];

	//sigs_ids is a vector containing the sids of the components in the right order
	size_t nComponents = sigs_ids.size();

	//rec should contain space for the required sids, and also, if used, the timer signal (last one)
	if (rec.usvs.size() < nComponents + 1)
		rec.usvs.resize(nComponents + 1);

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end());

	iteratorSignalIds.insert(v_sid);

	if (useTimer)
		iteratorSignalIds.insert(timer_signal_id);

	versionIterator vit(rec, iteratorSignalIds);

	//auto it = calc2req_sigs.find(calculator);
	//if (it == calc2req_sigs.end()) {//unexpected
	//	cout << "cannot find name " << calculator << " in calc2req_sigs" << endl;
	//	return -1;
	//}

	for (int iver = vit.init(); iver >= 0; iver = vit.next()) {
		//componentData will hold the data of component signals after alignment to relevant times
		//for each time we take the most recent valid data or missing_data otherwise.

		vector<vector<float> > componentData(nComponents + 1);
		vector<vector<int> > times(nComponents + 1); //times of data for components
		vector<vector<float> > vals(nComponents + 1); //(processed) vals of data for components

		for (size_t i = 0; i < nComponents; ++i) {
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);
			process_hosp_signal(signals[i], rec.usvs[i], times[i], vals[i]);
		}

		if (useTimer) {
			rec.uget(timer_signal_id, iver, rec.usvs[nComponents]);
			process_hosp_signal(timer_signal, rec.usvs[nComponents], times[nComponents], vals[nComponents]);
		}

		//if there is a timer_id, we use those times. Otherwise, we use the union
		//of times of all component times which don't belong to missing values
		vector<int> iverTimes;

		if (useTimer) {
			const vector<int>& curTimes = times[nComponents];
			iverTimes.reserve(curTimes.size());

			for (auto t : curTimes) {
				if (t > time_points[iver])
					break;
				else
					iverTimes.push_back(t);
			}
		}
		else {
			set<int> timeSet;
			for (int k = 0; k < nComponents; ++k) {
				const vector<int>& curTimes = times[k];
				for (auto t : curTimes) {
					if (t > time_points[iver])
						break;
					else
						timeSet.insert(t);
				}
			}
			iverTimes.reserve(timeSet.size());
			for (auto x : timeSet)
				iverTimes.push_back(x);
		}

		//if time_step > 0 we add a time before everything, at time_points[iver] and a time every time_step before the first time
		if (time_step > 0) {
			if (iverTimes.empty()) {
				iverTimes.push_back(beforeEverything());
				iverTimes.push_back(time_points[iver]);
			}
			else {
				int firstTime = iverTimes[0];
				iverTimes.reserve(iverTimes.size() + (time_points[iver] - firstTime) / time_step + 2);
				while (firstTime <= time_points[iver]) {
					iverTimes.push_back(firstTime);
					firstTime += time_step;
				}
				set<int> timeSet(iverTimes.begin(), iverTimes.end());
				timeSet.insert(beforeEverything());
				timeSet.insert(time_points[iver]);
				iverTimes = vector<int>(timeSet.begin(), timeSet.end());
			}
		}

		int nPoints = (int)iverTimes.size();

		//get component data
		for (size_t i = 0; i < nComponents; ++i) {
			vector<size_t> indices;
			int max_diff = -1; 
			int diff_time_unit = med_rep_type.windowTimeUnit;
			index_targets_in_given_vector(iverTimes, times[i], max_diff, signals_time_unit, diff_time_unit, time_points[iver], indices, onlyPast);
			const vector<float> & curVals = vals[i];

			//get the values from the correct indices
			componentData[i].resize(nPoints);

			for (int j = 0; j < nPoints; j++) {
				if (indices[j] == -1)
					componentData[i][j] = badValue();
				else
					componentData[i][j] = curVals.at(indices[j]);
			}
		}

		//we now have a matrix of data in componentData. 
		vector<int> finalTimes(nPoints); 
		vector<float> finalValues(nPoints);
		vector<float> inValues(nComponents); //for a given time, the inputs
		int nGoodPoints = 0;

		for (int j = 0; j < nPoints; j++) {
			for (int i = 0; i < nComponents; ++i)
				inValues[i] = componentData[i][j];

			float val = calcFunc(inValues, coeff);

			if (!isMissingValue(val)) {
				finalTimes[nGoodPoints] = iverTimes[j];
				finalValues[nGoodPoints] = val;
				nGoodPoints++;
			}
		}

		// replace bad value markers with a default value given in the missing_value field
		if (missing_value != badValue()) {
			for (auto& x : finalValues)
				if (x == badValue())
					x = missing_value;
		}

		// pushing virtual data into rec (into orig version), if there is data to push
		if (nGoodPoints > 0)
			rec.set_version_universal_data(v_sid, iver, &(finalTimes[0]), &(finalValues[0]), nGoodPoints);
	}

	return 0;
}

int RepCalcSimpleSignals::_apply_calc_hosp_time_dependent_pointwise(PidDynamicRec& rec, vector<int>& time_points,
	float(*calcFunc)(const vector<pair<int, float> >&, int, const vector<float>&), bool onlyPast) {

	//we either use a timer signal times or the union of component times
	//that don't belong to missing values
	bool useTimer = (timer_signal_id != -1);

	//Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	//sid for the calculated signal
	int v_sid = V_ids[0];

	//sigs_ids is a vector containing the sids of the components in the right order
	size_t nComponents = sigs_ids.size();

	//rec should contain space for the required sids, and also, if used, the timer signal (last one)
	if (rec.usvs.size() < nComponents + 1)
		rec.usvs.resize(nComponents + 1);

	// Loop on versions
	set<int> iteratorSignalIds(sigs_ids.begin(), sigs_ids.end());

	iteratorSignalIds.insert(v_sid);

	if (useTimer)
		iteratorSignalIds.insert(timer_signal_id);

	versionIterator vit(rec, iteratorSignalIds);

	//auto it = calc2req_sigs.find(calculator);
	//if (it == calc2req_sigs.end()) {//unexpected
	//	cout << "cannot find name " << calculator << " in calc2req_sigs" << endl;
	//	return -1;
	//}

	for (int iver = vit.init(); iver >= 0; iver = vit.next()) {
		//componentData will hold the data of component signals after alignment to relevant times
		//for each time we take the most recent valid data or missing_data otherwise.

		vector<vector<pair<int, float> > > componentData(nComponents + 1);
		vector<vector<int> > times(nComponents + 1); //times of data for components
		vector<vector<float> > vals(nComponents + 1); //(processed) vals of data for components

		for (size_t i = 0; i < nComponents; ++i) {
			rec.uget(sigs_ids[i], iver, rec.usvs[i]);
			process_hosp_signal(signals[i], rec.usvs[i], times[i], vals[i]);
		}

		if (useTimer) {
			rec.uget(timer_signal_id, iver, rec.usvs[nComponents]);
			process_hosp_signal(timer_signal, rec.usvs[nComponents], times[nComponents], vals[nComponents]);
		}

		//if there is a timer_id, we use those times. Otherwise, we use the union
		//of times of all component times which don't belong to missing values
		vector<int> iverTimes;

		if (useTimer) {
			const vector<int>& curTimes = times[nComponents];
			iverTimes.reserve(curTimes.size());

			for (auto t : curTimes) {
				if (t > time_points[iver])
					break;
				else
					iverTimes.push_back(t);
			}
		}
		else {
			set<int> timeSet;
			for (int k = 0; k < nComponents; ++k) {
				const vector<int>& curTimes = times[k];
				for (auto t : curTimes) {
					if (t > time_points[iver])
						break;
					else
						timeSet.insert(t);
				}
			}
			iverTimes.reserve(timeSet.size());
			for (auto x : timeSet)
				iverTimes.push_back(x);
		}

		//if time_step > 0 we add a time before everything, at time_points[iver] and a time every time_step before the first time
		if (time_step > 0) {
			if (iverTimes.empty()) {
				iverTimes.push_back(beforeEverything());
				iverTimes.push_back(time_points[iver]);
			}
			else {
				int firstTime = iverTimes[0];
				iverTimes.reserve(iverTimes.size() + (time_points[iver] - firstTime) / time_step + 2);
				while (firstTime <= time_points[iver]) {
					iverTimes.push_back(firstTime);
					firstTime += time_step;
				}
				set<int> timeSet(iverTimes.begin(), iverTimes.end());
				timeSet.insert(beforeEverything());
				timeSet.insert(time_points[iver]);
				iverTimes = vector<int>(timeSet.begin(), timeSet.end());
			}
		}

		int nPoints = (int)iverTimes.size();

		//get component data
		for (size_t i = 0; i < nComponents; ++i) {
			vector<size_t> indices;
			int max_diff = -1;
			int diff_time_unit = med_rep_type.windowTimeUnit;
			index_targets_in_given_vector(iverTimes, times[i], max_diff, signals_time_unit, diff_time_unit, time_points[iver], indices, onlyPast);
			const vector<float> & curVals = vals[i];
			const vector<int>& curTimes = times[i];

			//get the values from the correct indices
			componentData[i].resize(nPoints);

			for (int j = 0; j < nPoints; j++) {
				if (indices[j] == -1)
					componentData[i][j] = pair<int, float>(beforeEverything(), badValue());
				else
					componentData[i][j] = pair<int, float>(curTimes.at(indices[j]), curVals.at(indices[j]));
			}
		}

		//we now have a matrix of data in componentData. 
		vector<int> finalTimes(nPoints);
		vector<float> finalValues(nPoints);
		vector<pair<int,float> > inValues(nComponents); //for a given time, the inputs
		int nGoodPoints = 0;

		for (int j = 0; j < nPoints; j++) {
			for (int i = 0; i < nComponents; ++i)
				inValues[i] = componentData[i][j];

			float val = calcFunc(inValues, iverTimes[j], coeff);

			if (!isMissingValue(val)) {
				finalTimes[nGoodPoints] = iverTimes[j];
				finalValues[nGoodPoints] = val;
				nGoodPoints++;
			}
		}

		// replace bad value markers with a default value given in the missing_value field
		if (missing_value != badValue()) {
			for (auto& x : finalValues)
				if (x == badValue())
					x = missing_value;
		}

		// pushing virtual data into rec (into orig version)
		if (nGoodPoints > 0)
			rec.set_version_universal_data(v_sid, iver, &(finalTimes[0]), &(finalValues[0]), nGoodPoints);
	}

	return 0;
}

//.......................................................................................

float RepCalcSimpleSignals::calc_egfr_ckd_epi(float creatinine, int gender, float age, int ethnicity)
{
	double eGFR_CKD_EPI = pow(0.993, (double)age);

	if (creatinine <= 0)
		return -1;

	if (ethnicity == 1)
		eGFR_CKD_EPI *= 1.159;

	if (gender == 1) {
		// Male
		eGFR_CKD_EPI *= 141.0;
		if (creatinine <= 0.9)
			eGFR_CKD_EPI *= pow(creatinine/0.9, -0.441);
		else
			eGFR_CKD_EPI *= pow(creatinine/0.9, -1.209);
	}
	else {
		// Female
		eGFR_CKD_EPI *= 144.0;
		if (creatinine <= 0.7)
			eGFR_CKD_EPI *= pow(creatinine/0.7, -0.329);
		else
			eGFR_CKD_EPI *= pow(creatinine/0.7, -1.209);
	}

	return (float)eGFR_CKD_EPI;
}

//.......................................................................................
int RepCalcSimpleSignals::_apply_calc_eGFR(PidDynamicRec& rec, vector<int>& time_points)
{

	// Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int v_eGFR_sid = V_ids[0];

	// next must be in the same order of required signals in the calc2req_sigs map
	int Creatinine_sid = sigs_ids[0];
	int GENDER_sid = sigs_ids[1];
	int BYEAR_sid = sigs_ids[2];
	
	if (rec.usvs.size() < 3) rec.usvs.resize(3);
	
	// Loop on versions of Creatinine
	set<int> iteratorSignalIds ={ Creatinine_sid,v_eGFR_sid };
	versionIterator vit(rec, iteratorSignalIds);
	int gender = -1, byear = -1;
	int ethnicity = 0; // currently assuming 0, not having a signal for it

	for (int iver = vit.init(); iver>=0; iver = vit.next_different()) {

		rec.uget(Creatinine_sid, iver, rec.usvs[0]); // getting Creatinine into usv[0]
	//	MLOG("###>>>> pid %d version %d Creatinine %d , %s :: len %d\n", rec.pid, iver, Creatinine_sid, rec.my_base_rep->sigs.name(Creatinine_sid).c_str(), rec.usvs[0].len);

		int idx = 0;
		if (rec.usvs[0].len) {

			if (gender == -1) { // Assuming that birth-year and gender are version-independent

				rec.uget(GENDER_sid, iver, rec.usvs[1]); // getting GENDER into usv[1]
				rec.uget(BYEAR_sid, iver, rec.usvs[2]); // getting BYEAR into usv[2]

				gender = (int)rec.usvs[1].Val(0);
				byear = (int)rec.usvs[2].Val(0);
			}

			vector<float> v_vals(rec.usvs[0].len);
			vector<int> v_times(rec.usvs[0].len);

			// calculating for each creatinine point up to the relevant time-point
			for (int i = 0; i<rec.usvs[0].len; i++) {
				if (rec.usvs[0].Time(i) > time_points[iver])
					break;

				float age = RepCalcSimpleSignals::get_age(byear, rec.usvs[0].Time(i));
				float creatinine = rec.usvs[0].Val(i);
				v_times[idx] = rec.usvs[0].Time(i);
				v_vals[idx] = RepCalcSimpleSignals::calc_egfr_ckd_epi(creatinine, gender, age, ethnicity);
//				MLOG("pid %d ver %d time %d i %d creatinine %f gender %d age %f ethnicity %d : efgr = %f\n", rec.pid, iver, v_times[i], i, creatinine, gender, age, ethnicity, v_vals[i]);

				// Insert only legal values
				if (v_vals[i] >= 0)
					idx++;
			}

			// pushing virtual data into rec (into orig version)
			rec.set_version_universal_data(v_eGFR_sid, iver, &v_times[0], &v_vals[0], idx);
		}
	}

	
	return 0;
}


//.......................................................................................
int RepCalcSimpleSignals::_apply_calc_debug(PidDynamicRec& rec, vector<int>& time_points)
{

	// Check that we have the correct number of dynamic-versions : one per time-point
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	int v_debug_sid = V_ids[0];

	// next must be in the same order of required signals in the calc2req_sigs map
	int sid = sigs_ids[0];

	UniversalSigVec usv;

	// Loop on versions of Creatinine
	set<int> iteratorSignalIds ={ sid , v_debug_sid };
	versionIterator vit(rec, iteratorSignalIds);
	int gender = -1, byear = -1;
	int ethnicity = 0; // currently assuming 0, not having a signal for it

	for (int iver = vit.init(); iver>=0; iver = vit.next_different()) {

		rec.uget(sid, iver, usv); // getting sid into usv[0]
		//MLOG("###>>>> calc_debug : v_sig %d : pid %d version %d sig  %d , %s :: len %d\n", v_debug_sid, rec.pid, iver, sid, rec.my_base_rep->sigs.name(sid).c_str(), usv.len);

		int idx = 0;
		if (usv.len > 1) {

			vector<float> v_vals(usv.len-1);
			vector<int> v_times(usv.len-1);

			// calculating for each creatinine point up to the relevant time-point
			for (int i = 1; i<usv.len; i++) {
				if (rec.usvs[0].Time(i) > time_points[iver])
					break;

				v_times[idx] = usv.Time(i);
				v_vals[idx] = usv.Val(i) - usv.Val(i-1);

				idx++;
			}

			// pushing virtual data into rec (into orig version)
			rec.set_version_universal_data(v_debug_sid, iver, &v_times[0], &v_vals[0], idx);
		}
	}

	//rec.print_sigs({ "Creatinine","Hemoglobin","calc_eGFR","calc_debug" });
	return 0;
}

//.......................................................................................

//combine two sources of medication: first (med2kg) is normalized by kg, the other (med) is not, and weight (kg) in kg is provided
//of course, the dosages should be otherwise in the same units
float RepCalcSimpleSignals::med2KgPlusMed(float med2kg, float med, float kg) {
	if (isMissingValue(med2kg) && (isMissingValue(med) || isMissingValue(kg)))
		return badValue();
	else {
		if (isMissingValue(med2kg)) med2kg = 0.0F;
		if (isMissingValue(med) || kg <= 0 || isMissingValue(kg)) { med = 0.0F; kg = 1.0F; }
		return med2kg + med / kg;
	}
}

//.......................................................................................

//given several signals (<data>, a vector of time + val) at a given time (<time>), 
//interleave data as follows: consider non-missing values, and take the value with the most recent
//time. If all are "never", say badValue(). if there are multiple values for the same time, take the source with the smaller index 
//(in the vector of input features). Example of use: invasive and non-invasive bp.
float RepCalcSimpleSignals::interleave(const vector<pair<int, float> >& data, int time, const vector<float>& params) {
	int len = (int)data.size();
	float val = badValue();
	int bestTime = beforeEverything();

	for (int j = 0; j < len; ++j) {
		pair<int, float> x = data[j];				
		if (!isMissingValue(x.second) && x.first > bestTime) {
			val = x.second;
			bestTime = x.first;					
		}
	}

	return val;
}

//.......................................................................................

//given several signals (<data>, a vector of time + val) at a given time (<time>), 
//return 1.0 if at least one was not missing and with actual time at most params[0] before <time>, 0 otherwise
float RepCalcSimpleSignals::anySeenRecently(const vector<pair<int, float> >& data, int time, const vector<float>& params) {
	int allowedDelay = (int)(params[0]);
	int len = (int)data.size();

	for (int i = 0; i < len; ++i) {
		pair<int, float> x = data[i];
		if (!isMissingValue(x.second) && (time - x.first <= allowedDelay))
			return 1.0F;
	}

	return 0.0F;
}

//.......................................................................................

int RepCalcSimpleSignals::process_hosp_signal(const string& name, UniversalSigVec& usv, vector<int>& times, vector<float>& vals) {
	int len = usv.len;
	times.clear(); times.reserve(len + 1);
	vals.clear(); vals.reserve(len + 1);

	//get type
	SigType sigType = usv.get_type();

	//handle fields with Amount Per Time Range semantic
	if (sigType == T_TimeRangeVal && (name.substr(0, 6) == "INPUT_" || name.substr(0, 7) == "OUTPUT_" || name.substr(0, 12) == "PRESCRIPTION")) {

		//the signal is understood as sum total over a time period; however, "spikes" are allowed, namely, 
		//a time range of size 0. In addition, overlaps are allowed. Important: only the end times are originally provided, so only cumulative
		//quantities can be calculated.
		//issues: 
		// - there are signals with the same end time. They are summed
		// - the first end time records an amount that accumulated since nobody knows when. It cannot be used to calculate rate. 	
		// note: start_time = -1 signifies missing data (immaterial here). end_time cannot be -1 that way
		// negative vals are turned to 0.
		// if end time smaller than start time, data point is ignored
		// in general, the signal is sorted by the time end points and, then replaced by delta / delta_t (discrete derivative)

		if (len == 0)
			return 0;

		vector<pair<int, float> > v1; v1.reserve(len + 1);

		for (int i = 0; i < len; ++i) {
			int startTime = usv.Time(i, 0);
			int endTime = usv.Time(i, 1);

			if (endTime < startTime)
				continue;

			float val = max(0.0F, usv.Val(i, 0));

			if (endTime < startTime)
				continue;

			v1.push_back(pair<int, float>(endTime, val));
		}

		sort(v1.begin(), v1.end());

		//skip first time (v1 may also be empty, or length 1, so loop could be ignored)
		int index = -1;
		for (int i = 1; i < v1.size(); ++i) {
			if (v1[i].first != v1[i - 1].first) {
				index = i;
				break;
			}
		}

		if (index == -1) //all have the same, first, time, or v1.size()<=1
			return 0;

		//at this point, v1.size()>1, index points to the second unique end time

		//append a fictitious element to v1 with time higher than all entries to simplify handling
		v1.push_back(pair<int, float>(v1.back().first + 1, 0.0F));

		float curSum = 0;
		int curCount = 0;
		long long prevTime = v1[index - 1].first; //the first unique end time

		//note that index+1 is valid in v1

		for (int i = index + 1; i < v1.size(); ++i) {
			curSum += v1[i - 1].second;
			curCount++;
			if (v1[i].first != v1[i - 1].first) {
				float value = curSum / (v1[i - 1].first - prevTime);
				prevTime = v1[i - 1].first;
				times.push_back(v1[i - 1].first);
				vals.push_back(value);
				curSum = 0;
				curCount = 0;
			}
		}

		return 0;
	}

	//handle all other signals by simply taking the value and time from the 0-th channels

	for (int i = 0; i < len; ++i) {
		times.push_back(usv.Time(i, 0));
		vals.push_back(usv.Val(i, 0));
	}
	
	return 0;
}

