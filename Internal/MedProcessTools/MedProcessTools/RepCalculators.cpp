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