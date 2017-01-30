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
float get_eGFR_CKD_EPI(float age, float creatinine, int gender, int ethnicity=0);
float get_eGFR_MDRD(float age, float creatinine, int gender, int ethnicity=0);
float get_Framingham(float age, float total_cholesterol, float hdl, float bp_systolic, int smoking, int gender);


//=================================================================================
// Registries helpers
//=================================================================================

// data_mode can be mhs or thin (if left empty default is mhs)
int get_diabetes_dates(MedRepository &rep, int pid, string data_mode, int &last_healthy_date, int &first_pre_diabetes_date, int &first_diabetes_date);


#endif