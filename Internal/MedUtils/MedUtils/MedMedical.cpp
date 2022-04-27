#include "MedMedical.h"
#include "MedGenUtils.h"

#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_MED_UTILS
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//=================================================================================
// Calculated sigs
//=================================================================================
// unless otherwise stated gender is 1 for males and 2 for females
//

//---------------------------------------------------------------------------------------------------------------------------
float get_KFRE_Model_2( float age,	int gender,	float eGFR)
{
	vector <float> X(3);
	
	// unless otherwise stated gender is 1 for males and 2 for females
	if (gender == 1)
		X[0] = 1.;
	else
		X[0] = 0.;

	X[1] = age / 10;
	X[2] = eGFR / 5;

#ifdef KFRE_DEBUG
	for (int i = 0; i<X.size(); i++)
		cout << "X[" << i << "] = " << X[i] << endl;
#endif

	vector <float> Coeff = {
		(float)0.37548,
		(float)-0.29351,
		(float)-0.61217,
	};

	vector <float> Xbar = {
		(float)0.5642,
		(float)7.0355,
		(float)7.2216,
	};

	vector <float> betaXbar;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaXbar.push_back(Coeff[i] * Xbar[i]);
	}

	vector <float> betaX;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaX.push_back(Coeff[i] * X[i]);
	}

#ifdef KFRE_DEBUG
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		cout << "betaXbar[" << i << "] = " << betaXbar[i] << endl;
	}

	for (unsigned int i = 0; i < Coeff.size(); i++) {
		cout << "BetaX[" << i << "] = " << betaX[i] << endl;
	}
#endif

	float betaXbar_sum = std::accumulate(
		betaXbar.begin(),
		betaXbar.end(),
		decltype(betaXbar)::value_type(0)
	);

	float betaX_sum = std::accumulate(
		betaX.begin(),
		betaX.end(),
		decltype(betaX)::value_type(0)
	);

	float baseline = (float)0.916;
	float risk = 1 - pow(baseline, exp(betaX_sum - betaXbar_sum));

	return risk;
}

float get_KFRE_Model_3(
	float age,
	int gender,
	float eGFR,
	float UACR)
{
	//	Validate ranges, e.g. UACR>0
	if (UACR <= 0)
		return false;

	vector <float> X(4);

	if (gender > 0)
		X[0] = 1.;
	else
		X[0] = 0.;

	X[1] = age / 10;
	X[2] = eGFR / 5;
	X[3] = log(UACR);

#ifdef KFRE_DEBUG
	for (int i = 0; i<X.size(); i++)
		cout << "X[" << i << "] = " << X[i] << endl;
#endif

	vector <float> Coeff = {
		(float)0.26940,
		(float)-0.21670,
		(float)-0.55418,
		(float)0.45608,
	};

	vector <float> Xbar = {
		(float)0.5642,
		(float)7.0355,
		(float)7.2216,
		(float)5.2774,
	};

	vector <float> betaXbar;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaXbar.push_back(Coeff[i] * Xbar[i]);
	}

	vector <float> betaX;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaX.push_back(Coeff[i] * X[i]);
	}

#ifdef KFRE_DEBUG
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		cout << "betaXbar[" << i << "] = " << betaXbar[i] << endl;
	}

	for (unsigned int i = 0; i < Coeff.size(); i++) {
		cout << "BetaX[" << i << "] = " << betaX[i] << endl;
	}
#endif

	float betaXbar_sum = std::accumulate(
		betaXbar.begin(),
		betaXbar.end(),
		decltype(betaXbar)::value_type(0)
	);

	float betaX_sum = std::accumulate(
		betaX.begin(),
		betaX.end(),
		decltype(betaX)::value_type(0)
	);

	float baseline = (float)0.924;
	float risk = 1 - pow(baseline, exp(betaX_sum - betaXbar_sum));

	return risk;

}

float get_KFRE_Model_6(
	float age,
	int gender,
	float eGFR,
	float UACR,
	float Calcium,
	float Phosphorus,
	float Albumin,
	float Bicarbonate)
{
	//	Validate ranges, e.g. UACR>0
	if (UACR <= 0)
		return false;

	vector <float> X(8);

	if (gender > 0)
		X[0] = 1.;
	else
		X[0] = 0.;

	X[1] = age / 10;
	X[2] = eGFR / 5;
	X[3] = log(UACR);
	X[4] = Calcium;
	X[5] = Phosphorus;
	X[6] = Albumin;
	X[7] = Bicarbonate;

#ifdef KFRE_DEBUG
	for (int i = 0; i<X.size(); i++)
		cout << "X[" << i << "] = " << X[i] << endl;
#endif

	vector <float> Coeff = {
		(float)0.16117,
		(float)-0.19883,
		(float)-0.49360,
		(float)0.35066,
		(float)-0.22129,
		(float)0.24197,
		(float)-0.33867,
		(float)-0.07429,
	};

	vector <float> Xbar = {
		(float)0.5642,
		(float)7.0355,
		(float)7.2216,
		(float)5.2774,
		(float)9.3510,
		(float)3.9221,
		(float)3.9925,
		(float)25.5441,
	};

	vector <float> betaXbar;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaXbar.push_back(Coeff[i] * Xbar[i]);
	}

	vector <float> betaX;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaX.push_back(Coeff[i] * X[i]);
	}

#ifdef KFRE_DEBUG
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		cout << "betaXbar[" << i << "] = " << betaXbar[i] << endl;
	}

	for (unsigned int i = 0; i < Coeff.size(); i++) {
		cout << "BetaX[" << i << "] = " << betaX[i] << endl;
	}
#endif

	float betaXbar_sum = std::accumulate(
		betaXbar.begin(),
		betaXbar.end(),
		decltype(betaXbar)::value_type(0)
	);

	float betaX_sum = std::accumulate(
		betaX.begin(),
		betaX.end(),
		decltype(betaX)::value_type(0)
	);

	float baseline = (float)0.929;
	float delta = betaX_sum - betaXbar_sum;
	
	errno = 0;
	float risk = 1 - pow(baseline, exp(delta));

	if (errno == ERANGE) {
		printf("exp(%f) overflows\n", delta);

		cout << "age " << age << endl;
		cout << "age " << age << endl;
		cout << "age " << age << endl;
		cout << "age " << age << endl;
		cout << "age " << age << endl;
		cout << "age " << age << endl;2w

			int gender,
			float eGFR,
			float UACR,
			float Calcium,
			float Phosphorus,
			float Albumin,
			float Bicarbonate

		return -1.;
	}
	return risk;

}



//---------------------------------------------------------------------------------------------------------------------------
float get_eGFR_CKD_EPI(float age, float creatinine, int gender, int ethnicity)
{
	double eGFR_CKD_EPI = pow(0.993, (double)age);

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

//---------------------------------------------------------------------------------------------------------------------------
float get_eGFR_MDRD(float age, float creatinine, int gender, int ethnicity)
{
	if (age <= 1 || creatinine <= 0.1) return -1;

	double eGFR_MDRD = 175.0 * pow((double)creatinine, -1.154) * pow((double)age, -0.203);
	if (gender == 2) eGFR_MDRD *= 0.742;
	if (ethnicity == 1) eGFR_MDRD *= 1.212;

	return ((float)eGFR_MDRD);
}

//---------------------------------------------------------------------------------------------------------------------------
float get_Framingham(float age, float total_cholesterol, float hdl, float bp_systolic, int smoking, int gender)
{
	//int framingham = 0;

	//if (gender == 1) {

	//	// Age
	//	if (age <= 34) framingham -= 7;
	//	else if (age <= 39) framingham -= 3;
	//	else if (age <= 44) framingham += 0;
	//	else if (age <= 49) framingham += 3;
	//	else if (age <= 54) framingham += 6;
	//	else if (age <= 59) framingham += 8;
	//	else if (age <= 64) framingham += 10;
	//	else if (age <= 69) framingham += 12;
	//	else if (age <= 74) framingham += 14;
	//	else framingham += 16;

	//	// Cholesterol 

	//}
	MERR("Framingham score not implemented yet !!!\n");
	return -1;
}


//=================================================================================
// Registries helpers
//=================================================================================

// data_mode can be mhs or thin (if left empty it will be detected automatically using the type of the Drug )
int get_diabetes_dates(MedRepository &rep, int pid, string data_mode, int &last_healthy_date, int &first_pre_diabetes_date, int &first_diabetes_date)
{
	last_healthy_date = 0;
	first_pre_diabetes_date = 0;
	first_diabetes_date = 0;

	int evidence = 0;

	if (rep.sigs.sid("Drug") < 0 || rep.sigs.sid("Glucose") < 0 || rep.sigs.sid("HbA1C") < 0)
		return 1;
	// assumes rep was already loaded with Glucose , HbA1C and Drugs
	if (data_mode == "") {
		int dsid = rep.sigs.sid("Drug");
		int dtype = rep.sigs.Sid2Info[dsid].type;
		if (dtype == T_DateVal2) data_mode = "mhs";
		if (dtype == T_DateShort2) data_mode = "thin";
	}

	int glu_len, hba1c_len;
	SDateVal *glu_sdv = (SDateVal *)rep.get(pid, "Glucose", glu_len);
	SDateVal *hba1c_sdv = (SDateVal *)rep.get(pid, "HbA1C", hba1c_len);

	if (glu_len == 0 && hba1c_len == 0) return 1;

	vector<vector<int>> events; // vector of triplets <date> <test 0: glu 1: hba1c 2: drugs> <type: 0 - healthy 1 - pre diabetic 2 - diabetic (but need 2 of those) , 3 diabetic even with a signle test>

	int type;

	// glucose and hba1c
	for (int i=0; i<glu_len; i++) {
		if (glu_sdv[i].val <= 100) type = 0;
		else if (glu_sdv[i].val <= 125) type = 1;
		else if (glu_sdv[i].val <= 300) type = 2;
		else type = 3;
		if (type > 0) evidence |= 10;
		events.push_back({ glu_sdv[i].date, 0, type, date_to_days(glu_sdv[i].date) });
	}


	for (int i=0; i<hba1c_len; i++) {
		if (hba1c_sdv[i].val <= 5.7) type = 0;
		else if (hba1c_sdv[i].val <= 6.4) type = 1;
		else if (hba1c_sdv[i].val <= 8.5) type = 2;
		else type = 3;
		if (type > 0) evidence |= 100;
		events.push_back({ hba1c_sdv[i].date, 1, type, date_to_days(hba1c_sdv[i].date) });
	}


	// Drugs
	int min_days = 30;
	int first_date = 0;
	int sum_days = 0;

	if (data_mode == "mhs") {
		int drug_len;
		SDateVal2 *drug_sdv2 = (SDateVal2 *)rep.get(pid, "Drug", drug_len);
		int section_id = rep.dict.section_id("Drug");
		int drug_set = rep.dict.id(section_id, "ATC_A10_____");
		int is_in = 0;
		for (int i=0; i<drug_len; i++) {
			if ((is_in = rep.dict.is_in_set(section_id, (int)drug_sdv2[i].val, drug_set))) {
				if (first_date == 0) first_date = drug_sdv2[i].date;
				sum_days += drug_sdv2[i].val2;
				if (sum_days > min_days) {
					evidence += 1;
					break;
				}
				MLOG("is_in %d :: %s\n", is_in, rep.dict.name((int)drug_sdv2[i].val).c_str());
			}
		}
	}

	if (data_mode == "thin") {
		int drug_len;
		SDateShort2 *drug_sds2 = (SDateShort2 *)rep.get(pid, "Drug", drug_len);
		int section_id = rep.dict.section_id("Drug");
		int drug_set = rep.dict.id(section_id, "ATC_A10_____");
		for (int i=0; i<drug_len; i++) {
			if (rep.dict.is_in_set(section_id, (int)drug_sds2[i].val1, drug_set)) {
				if (first_date == 0) first_date = drug_sds2[i].date;
				sum_days += drug_sds2[i].val2;
				if (sum_days > min_days) {
					evidence += 1;
					break;
				}
			}
		}
	}

	if (sum_days > min_days) {
		events.push_back({ first_date, 2, 3, date_to_days(first_date) });
	}

	std::sort(events.begin(), events.end(), [](const vector<int> &a, const vector<int> &b) { return (a[0]<b[0]);});

	//for (int i=0; i<events.size(); i++) {
	//	MLOG("event: %d %d %d %d\n", events[i][0], events[i][1], events[i][2], events[i][3]);
	//}
	
	// detect first diabetes date
	first_diabetes_date = 0;
	int days2Y = 730;
	for (int i=0; i<events.size(); i++) {

		// drugs event mark diabetes start
		if (events[i][2] == 3) {
			first_diabetes_date = events[i][0];
			break;
		}

		if (events[i][2] == 2) {
			int found_in_2y = 0;
			for (int j=i+1; j<events.size(); j++) {
				if (events[j][3]-events[i][3] <= days2Y && events[j][2] >= 2) {
					found_in_2y = 1;
					break;
				}
			}
			if (found_in_2y) {
				first_diabetes_date = events[i][0];
				break;
			}
		}

	}

	// detect pre diabetes first date
	first_pre_diabetes_date = 0;
	int last_d = first_diabetes_date;

	for (int i=0; i<events.size(); i++) {

		if (last_d > 0 && events[i][0] >= last_d)
			break;

		if (events[i][2] >= 1) {
			int found_in_2y = 0;
			for (int j=i+1; j<events.size(); j++) {
				if (events[j][3]-events[i][3] <= days2Y && events[j][2] >= 1) {
					found_in_2y = 1;
					break;
				}
			}
			if (found_in_2y) {
				first_pre_diabetes_date = events[i][0];
				break;
			}
		}

	}

	last_healthy_date = 0;
	last_d = first_pre_diabetes_date;
	if (last_d == 0) last_d = first_diabetes_date;

	for (auto &ev : events) {
		if (last_d>0 && ev[0]>=last_d)
			break;
		if (ev[2] == 0)
			last_healthy_date = ev[0];
		else
			if (ev[2] > 1)
				break;
	}

	MLOG("Diabetes dates: pid %d healthy %d pre-diabetes %d diabetes %d evidence %d\n",pid, last_healthy_date, first_pre_diabetes_date, first_diabetes_date, evidence);
	return 0; // not censored
}

