//
// MedMedical.h
//
// Helper methods to generate medical info from the data - such as registries, drug groups, etc...
//

#ifndef __MED_MEDICAL_H__
#define __MED_MEDICAL_H__

#include "InfraMed/InfraMed/InfraMed.h"

//=================================================================================
// Calculated sigs
//=================================================================================
// unless otherwise stated gender is 1 for males and 2 for females
//
float get_KFRE_Model_2(float age, int gender, float eGFR);
float get_KFRE_Model_3(float age, int gender, float eGFR, float UACR);
bool get_KFRE_Model_6(float &risk, float age, int gender, float eGFR, float UACR, float Calcium, float Phosphorus, float Albumin, float Bicarbonate);

bool FetchCoefficients_v1(
	int n_variables,
	double& baseline,
	vector<double>& Coeff,
	vector<double>& Xbar
);
bool FetchCoefficients(
	int n_variables,
	int prediction_years,
	int region_id,
	double& baseline,
	vector<double>& Coeff,
	vector<double>& Xbar
);
bool get_KFRE3(
	float &risk,
	double baseline,
	vector <double> Coeff,
	vector <double> Xbar,
	float age,
	int gender,
	float eGFR
);
bool get_KFRE4(
	float &risk,
	double baseline,
	vector <double> Coeff,
	vector <double> Xbar,
	float age,
	int gender,
	float eGFR,
	float UACR
);
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
);


float get_eGFR_CKD_EPI(float age, float creatinine, int gender, int ethnicity=0);
float get_eGFR_MDRD(float age, float creatinine, int gender, int ethnicity=0);
float get_Framingham(float age, float total_cholesterol, float hdl, float bp_systolic, int smoking, int gender);




//=================================================================================
// Registries helpers
//=================================================================================

// data_mode can be mhs or thin (if left empty default is mhs)
int get_diabetes_dates(MedRepository &rep, int pid, string data_mode, int &last_healthy_date, int &first_pre_diabetes_date, int &first_diabetes_date);


#endif