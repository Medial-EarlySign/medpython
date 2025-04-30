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

//=================================================================================
// KFRE coefficients
//=================================================================================

double baseline_model_6 = 0.929;

double coeff_model_6[8] = {
	-0.19883,
	0.16117,
	-0.49360,
	0.35066,
	-0.33867,
	0.24197,
	-0.07429,
	-0.22129
};

double xbar_model_6[8] = {
	7.0355,
	0.5642,
	7.2216,
	5.2774,
	3.9925,
	3.9221,
	25.5441,
	9.3510
};

double baseline_model_3 = 0.924;

double coeff_model_3[4] = {
	-0.21670,
	0.26940,
	-0.55418,
	0.45608
};

double xbar_model_3[4] = {
	7.0355,
	0.5642,
	7.2216,
	5.2774
};

double baseline_model_2 = 0.916;

double coeff_model_2[3] = {
	-0.29351,
	0.37548,
	-0.61217
};

double xbar_model_2[3] = {
	7.0355,
	0.5642,
	7.2216
};

double baseline_4[8] = {
	0.9750,
	0.9751,
	0.9832,
	0.9676,
	0.9240,
	0.8996,
	0.9365,
	0.8762,
};

double baseline_6[8] = {
	0.9750,
	0.9755,
	0.9830,
	0.9707,
	0.9240,
	0.9018,
	0.9370,
	0.8839,
};

double baseline_8[8] = {
	0.9780,
	0.9757,
	0.9827,
	0.9629,
	0.9301,
	0.9096,
	0.9245,
	0.8636
};

double coeff_4[8][4] = {
	{ -0.2201, +0.2467, -0.5567, +0.4510 },
{ -0.2201, +0.2467, -0.5567, +0.4510 },
{ -0.2201, +0.2467, -0.5567, +0.4510 },
{ -0.2245, +0.3212, -0.4553, +0.4469 },
{ -0.2201, +0.2467, -0.5567, +0.4510 },
{ -0.2201, +0.2467, -0.5567, +0.4510 },
{ -0.2201, +0.2467, -0.5567, +0.4510 },
{ -0.2245, +0.3212, -0.4553, +0.4469 },
};

double coeff_6[8][6] = {
	{ -0.2218, +0.2553, -0.5541, +0.4562, -0.1475, +0.1426 },
{ -0.2218, +0.2553, -0.5541, +0.4562, -0.1475, +0.1426 },
{ -0.2218, +0.2553, -0.5541, +0.4562, -0.1475, +0.1426 },
{ -0.2401, +0.3209, -0.4650, +0.4384, +0.3018, +0.1710 },
{ -0.2218, +0.2553, -0.5541, +0.4562, -0.1475, +0.1426 },
{ -0.2218, +0.2553, -0.5541, +0.4562, -0.1475, +0.1426 },
{ -0.2218, +0.2553, -0.5541, +0.4562, -0.1475, +0.1426 },
{ -0.2401, +0.3209, -0.4650, +0.4384, +0.3018, +0.1710 },
};

double coeff_8[8][8] = {
	{ -0.1992,   +0.1602,   -0.4919,   +0.3364,   -0.3441,   +0.2604,   -0.07354,   -0.2228 },
{ -0.1992,   +0.1602,   -0.4919,   +0.3364,   -0.3441,   +0.2604,   -0.07354,   -0.2228 },
{ -0.1992,   +0.1602,   -0.4919,   +0.3364,   -0.3441,   +0.2604,   -0.07354,   -0.2228 },
{ -0.1848,   +0.2906,   -0.4156,   +0.3480,   -0.3569,   +0.1582,   -0.01199,   -0.1581 },
{ -0.1992,   +0.1602,   -0.4919,   +0.3364,   -0.3441,   +0.2604,   -0.07354,   -0.2228 },
{ -0.1992,   +0.1602,   -0.4919,   +0.3364,   -0.3441,   +0.2604,   -0.07354,   -0.2228 },
{ -0.1992,   +0.1602,   -0.4919,   +0.3364,   -0.3441,   +0.2604,   -0.07354,   -0.2228 },
{ -0.1848,   +0.2906,   -0.4156,   +0.3480,   -0.3569,   +0.1582,   -0.01199,   -0.1581 },
};

double xbar_4[4] = { 7.036, 0.5642, 7.222, 5.137 };
double xbar_6[6] = { 7.036, 0.5642, 7.222, 5.137, 0.5106, 0.8501 };
double xbar_8[8] = { 7.036, 0.5642, 7.222, 5.137, 3.997, 3.916, 25.57, 9.355 };

//---------------------------------------------------------------------------------------------------------------------------
float get_KFRE_Model_2( float age,	int gender,	float eGFR)
{
	vector <double> X(3);
	
	X[0] = gender;
	X[1] = age / 10;
	X[2] = eGFR / 5;

#ifdef KFRE_DEBUG
	for (int i = 0; i<X.size(); i++)
		cout << "X[" << i << "] = " << X[i] << endl;
#endif

	vector <double> Coeff = {
		(double)0.37548,
		(double)-0.29351,
		(double)-0.61217,
	};

	vector <double> Xbar = {
		(double)0.5642,
		(double)7.0355,
		(double)7.2216,
	};

	vector <double> betaXbar;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaXbar.push_back(Coeff[i] * Xbar[i]);
	}

	vector <double> betaX;
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

	double betaXbar_sum = std::accumulate(
		betaXbar.begin(),
		betaXbar.end(),
		decltype(betaXbar)::value_type(0)
	);

	double betaX_sum = std::accumulate(
		betaX.begin(),
		betaX.end(),
		decltype(betaX)::value_type(0)
	);

	double baseline = (double)0.916;
	float risk = 1 - (float)pow(baseline, exp(betaX_sum - betaXbar_sum));

	//printf("risk %f, age %f, gender %d, eGFR %f\n",
	//	risk,
	//	age,
	//	gender,
	//	eGFR);

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

	vector <double> X(4);

	if (gender == 1)
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

	vector <double> Coeff = {
		(double)0.26940,
		(double)-0.21670,
		(double)-0.55418,
		(double)0.45608,
	};

	vector <double> Xbar = {
		(double)0.5642,
		(double)7.0355,
		(double)7.2216,
		(double)5.2774,
	};

	vector <double> betaXbar;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaXbar.push_back(Coeff[i] * Xbar[i]);
	}

	vector <double> betaX;
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

	double betaXbar_sum = std::accumulate(
		betaXbar.begin(),
		betaXbar.end(),
		decltype(betaXbar)::value_type(0)
	);

	double betaX_sum = std::accumulate(
		betaX.begin(),
		betaX.end(),
		decltype(betaX)::value_type(0)
	);

	double baseline = (double)0.924;
	double risk = 1 - (float)pow(baseline, exp(betaX_sum - betaXbar_sum));

	//printf("risk %f, age %f, gender %d, eGFR %f, UACR %f\n",
	//	risk,
	//	age,
	//	gender,
	//	eGFR,
	//	UACR);

	return risk;

}

bool get_KFRE_Model_6(
	float &risk,
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

	vector <double> X(8);

	if (gender == 1)
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

	vector <double> Coeff = {
		(double)0.16117,
		(double)-0.19883,
		(double)-0.49360,
		(double)0.35066,
		(double)-0.22129,
		(double)0.24197,
		(double)-0.33867,
		(double)-0.07429,
	};

	vector <double> Xbar = {
		(double)0.5642,
		(double)7.0355,
		(double)7.2216,
		(double)5.2774,
		(double)9.3510,
		(double)3.9221,
		(double)3.9925,
		(double)25.5441,
	};

	vector <double> betaXbar;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaXbar.push_back(Coeff[i] * Xbar[i]);
	}

	vector <double> betaX;
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

	double betaXbar_sum = std::accumulate(
		betaXbar.begin(),
		betaXbar.end(),
		decltype(betaXbar)::value_type(0)
	);

	double betaX_sum = std::accumulate(
		betaX.begin(),
		betaX.end(),
		decltype(betaX)::value_type(0)
	);

	double baseline = (double)0.929;
	double delta = betaX_sum - betaXbar_sum;
	
	errno = 0;
	risk = 1 - (float)pow(baseline, exp(delta));

	if (errno == ERANGE) {

		//printf("exp(%f) overflows: age %f, gender %d, eGFR %f, UACR %f, Calcium %f , Phosphorus %f, Albumin %f, Bicarbonate %f\n", 
		//	delta, 
		//	age, 
		//	gender, 
		//	eGFR, 
		//	UACR, 
		//	Calcium, 
		//	Phosphorus, 
		//	Albumin, 
		//	Bicarbonate);
		
		return false;
	}
	else {
		//printf("risk %f, age %f, gender %d, eGFR %f, UACR %f, Calcium %f , Phosphorus %f, Albumin %f, Bicarbonate %f\n",  
		//	risk,
		//	age, 
		//	gender, 
		//	eGFR, 
		//	UACR, 
		//	Calcium, 
		//	Phosphorus, 
		//	Albumin, 
		//	Bicarbonate);
	}

	return true;
}

//---------------------------------------------------------------------------------------------------------------------------

bool FetchCoefficients_v1(
	int n_variables,
	double& baseline,
	vector<double>& Coeff,
	vector<double>& Xbar
)
{
	switch (n_variables) {
	case 3:
		baseline = baseline_model_2;
		for (int i = 0; i < n_variables; i++) {
			Coeff.push_back(coeff_model_2[i]);
			Xbar.push_back(xbar_model_2[i]);
		}
		break;
	case 4:
		baseline = baseline_model_3;
		for (int i = 0; i < n_variables; i++) {
			Coeff.push_back(coeff_model_3[i]);
			Xbar.push_back(xbar_model_3[i]);
		}
		break;
	case 8:
		baseline = baseline_model_6;
		for (int i = 0; i < n_variables; i++) {
			Coeff.push_back(coeff_model_6[i]);
			Xbar.push_back(xbar_model_6[i]);
		}
		break;
	default:
		break;
	}

	return true;
}

bool FetchCoefficients(
	int n_variables,
	int prediction_years,
	int region_id,
	double& baseline,
	vector<double>& Coeff,
	vector<double>& Xbar
)
{
	int offset = 0;
	if (prediction_years == 5)
		offset = 4;
	int table_row = offset + region_id;

	// This is ugly, should be replaced with setting pointers to respective table
	switch (n_variables) {
	case 4:
		baseline = baseline_4[table_row];
		for (int i = 0; i < n_variables; i++) {
			Coeff.push_back(coeff_4[table_row][i]);
			Xbar.push_back(xbar_4[i]);
		}
		break;
	case 6:
		baseline = baseline_6[table_row];
		for (int i = 0; i < n_variables; i++) {
			Coeff.push_back(coeff_6[table_row][i]);
			Xbar.push_back(xbar_6[i]);
		}
		break;
	case 8:
		baseline = baseline_8[table_row];
		for (int i = 0; i < n_variables; i++) {
			Coeff.push_back(coeff_8[table_row][i]);
			Xbar.push_back(xbar_8[i]);
		}
		break;
	default:
		cout << "ERROR: n_variables should be equal 4,6 or 8" << '\n';
		return false;
		break;
	}

	return true;
}

bool calc_KFRE_v2(float &risk,
	double baseline,
	vector <double> Coeff,
	vector <double> Xbar,
	vector <double> X
)
{
	/*
	X is a feature vector, after all scalings and log-transforms.
	We assume that lengths of Coeff, Xbar and X are identical
	*/

#ifdef KFRE_DEBUG
	for (int i = 0; i<X.size(); i++)
		cout << "X[" << i << "] = " << X[i] << endl;
#endif

	vector <double> betaXbar;
	for (unsigned int i = 0; i < Coeff.size(); i++) {
		betaXbar.push_back(Coeff[i] * Xbar[i]);
	}

	vector <double> betaX;
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

	double betaXbar_sum = std::accumulate(
		betaXbar.begin(),
		betaXbar.end(),
		decltype(betaXbar)::value_type(0)
	);

	double betaX_sum = std::accumulate(
		betaX.begin(),
		betaX.end(),
		decltype(betaX)::value_type(0)
	);

	double delta = betaX_sum - betaXbar_sum;

	errno = 0;
	risk = 1 - (float)pow(baseline, exp(delta));

	if (errno == ERANGE) {
#if 0
		printf("exp(%f) overflows: age %f, gender %d, eGFR %f, UACR %f, Calcium %f , Phosphorus %f, Albumin %f, Bicarbonate %f\n",
			delta,
			age,
			gender,
			eGFR,
			UACR,
			Calcium,
			Phosphorus,
			Albumin,
			Bicarbonate);
#endif
		return false;
	}
	else {
#if 0
		printf("risk %f, age %f, gender %d, eGFR %f, UACR %f, Calcium %f , Phosphorus %f, Albumin %f, Bicarbonate %f\n",
			risk,
			age,
			gender,
			eGFR,
			UACR,
			Calcium,
			Phosphorus,
			Albumin,
			Bicarbonate);
#endif
	}

	return true;
}

bool get_KFRE3(
	float &risk,
	double baseline,
	vector <double> Coeff,
	vector <double> Xbar,
	float age,
	int gender,
	float eGFR
)
{

	// Transform input variables [scale, log- et.c.]
	vector <double> X(3);

	// ATTENTION !!! 
	//
	// The order in v2 is DIFFERENT from what we used in v1

	X[0] = age / 10;
	if (gender == 1)
		X[1] = 1.;
	else
		X[1] = 0.;
	X[2] = eGFR / 5;

	return calc_KFRE_v2(
		risk,
		baseline,
		Coeff,
		Xbar,
		X
	);
}


bool get_KFRE4(
	float &risk,
	double baseline,
	vector <double> Coeff,
	vector <double> Xbar,
	float age,
	int gender,
	float eGFR,
	float UACR
)
{

	//	Validate ranges, e.g. UACR>0
	if (UACR <= 0)
		return false;

	// Transform input variables [scale, log- et.c.]
	vector <double> X(4);

	// ATTENTION !!! 
	//
	// The order in v2 is DIFFERENT from what we used in v1

	X[0] = age / 10;
	if (gender == 1)
		X[1] = 1.;
	else
		X[1] = 0.;
	X[2] = eGFR / 5;
	X[3] = log(UACR);

	return calc_KFRE_v2(
		risk,
		baseline,
		Coeff,
		Xbar,
		X
	);
}

bool get_KFRE8(
	float &risk,
	double baseline,
	vector <double> Coeff,
	vector <double> Xbar,
	float age,
	int gender,
	float eGFR,
	float UACR,
	float Calcium,
	float Phosphorus,
	float Albumin,
	float Bicarbonate
)
{

	//	Validate ranges, e.g. UACR>0
	if (UACR <= 0)
		return false;

	// Transform input variables [scale, log- et.c.]
	vector <double> X(8);

	// ATTENTION !!! 
	//
	// The order in v2 is DIFFERENT from what we used in v1

	X[0] = age / 10;
	if (gender == 1)
		X[1] = 1.;
	else
		X[1] = 0.;
	X[2] = eGFR / 5;
	X[3] = log(UACR);
	X[4] = Albumin;
	X[5] = Phosphorus;
	X[6] = Bicarbonate;
	X[7] = Calcium;

	return calc_KFRE_v2(
		risk,
		baseline,
		Coeff,
		Xbar,
		X
	);
}


//---------------------------------------------------------------------------------------------------------------------------
float get_eGFR_CKD_EPI(float age, float creatinine, int gender, int ethnicity)
{
	// 2021 version
	double eGFR_CKD_EPI = 142 * pow(0.9938, (double)age);

	if (gender == 1) {
		// Male
		if (creatinine <= 0.9)
			eGFR_CKD_EPI *= pow(creatinine / 0.9, -0.302);
		else
			eGFR_CKD_EPI *= pow(creatinine / 0.9, -1.200);
	}
	else {
		// Female
		eGFR_CKD_EPI *= 1.012;
		if (creatinine <= 0.7)
			eGFR_CKD_EPI *= pow(creatinine / 0.7, -0.241);
		else
			eGFR_CKD_EPI *= pow(creatinine / 0.7, -1.200);
	}


/* 2009 version, note that -0.441 was our typo, and should have been -0.411
/--------------------------------------
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

/--------------------------------------
*/
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

