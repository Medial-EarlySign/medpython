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
	// we ignore time_points in this case, as we simply calculate the eGFR for
	// each time point in which we did a Creatinine test

	int v_eGFR_sid = V_ids[0];

	// next must be in the same order of required signals in the calc2req_sigs map
	int Creatinine_sid = sigs_ids[0];
	int GENDER_sid = sigs_ids[1];
	int BYEAR_sid = sigs_ids[2];


	if (rec.usvs.size() < 3) rec.usvs.resize(3);

	rec.uget(Creatinine_sid, 0, rec.usvs[0]); // getting Creatinine into usv[0]
//	MLOG("###>>>> pid %d Creatinine %d , %s :: len %d\n", rec.pid, Creatinine_sid, rec.my_base_rep->sigs.name(Creatinine_sid).c_str(), rec.usvs[0].len);

	if (rec.usvs[0].len > 0) {

		rec.uget(GENDER_sid, 0, rec.usvs[1]); // getting GENDER into usv[1]
		rec.uget(BYEAR_sid, 0, rec.usvs[2]); // getting BYEAR into usv[2]

		int gender = (int)rec.usvs[1].Val(0);
		int byear = (int)rec.usvs[2].Val(0);
		int ethnicity = 0; // currently assuming 0, not having a signal for it

		vector<float> v_vals(rec.usvs[0].len);
		vector<int> v_times(rec.usvs[0].len);

		// calculating for each creatinine point
		for (int i=0; i<rec.usvs[0].len; i++) {
			float age = RepCalcSimpleSignals::get_age(byear, rec.usvs[0].Time(i));
			float creatinine = rec.usvs[0].Val(i);
			v_times[i] = rec.usvs[0].Time(i);
			v_vals[i] = RepCalcSimpleSignals::calc_egfr_ckd_epi(creatinine, gender, age, ethnicity);
			//MLOG("pid %d time %d i %d creatinine %f gender %d age %f ethnicity %d\n", rec.pid, v_times[i], i, creatinine, gender, age, ethnicity);
			if (v_vals[i] < 0) v_vals[i] = missing_value;
		}

		// pushing virtual data into rec (into orig version)
		rec.set_version_universal_data(v_eGFR_sid, 0, &v_times[0], &v_vals[0], rec.usvs[0].len);

		// pointing all versions to the 0 one
		for (int ver=1; ver<rec.get_n_versions(); ver++)
			rec.point_version_to(v_eGFR_sid, 0, ver);
	}
	return 0;
}