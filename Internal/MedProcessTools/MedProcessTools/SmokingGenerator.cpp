#include "SmokingGenerator.h"

void SmokingGenerator::set_names() {
	if (!names.empty())
		return;
	for (string s : {"Current_Smoker", "Ex_Smoker", "Smok_Years_Since_Quitting", "Smoking_Years", "Smok_Pack_Years"})
		names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);	
}

int SmokingGenerator::Generate(PidDynamicRec& rec, MedFeatures& features, int index, int num) {
	int missing = -1;

	for (int i = 0; i < num; i++) {
		int current_smoker, ex_smoker;
		int years_since_quitting, smoking_years;
		float pack_years;

		int len;
		bool never_smoked = true;

		string sname = "SMOKING_ENRICHED";
		SValShort4 *smx_status = (SValShort4 *)rec.get(sname, 0, len);
		if (len > 0)
			never_smoked = (smx_status[0].val1 == -1);
		assert(len <= 1);

		if (len == 0) { // No Data
			current_smoker = ex_smoker = (int)missing;
			years_since_quitting = smoking_years = (int)missing;
			pack_years = missing;
		}
		else if (never_smoked) { // Non Smoker
			current_smoker = ex_smoker = 0;
			years_since_quitting = 100;
			smoking_years = 0;
			pack_years = 0.0;
		}
		else { // (Ex)Smoker
			int start_year = smx_status[0].val1;
			int end_year = smx_status[0].val2;
			int target_year = (float)(med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time) / 10000);
			if (target_year < end_year) {
				// still in smoking period
				smoking_years = target_year - start_year;
				years_since_quitting = 0;
				current_smoker = 1;
			}
			else {
				// maybe done smoking
				current_smoker = smx_status[0].val4;
				smoking_years = end_year - start_year; // we are merciful
				if (!current_smoker)
					years_since_quitting = target_year - end_year;
				else
					years_since_quitting = 0;
			}
			pack_years = ((float)smx_status[0].val3 / 20) * smoking_years;
			ex_smoker = 1 - current_smoker;
		}

		features.data[names[0]][index + i] = (float)current_smoker; 
		features.data[names[1]][index + i] = (float)ex_smoker;
		features.data[names[2]][index + i] = (float)years_since_quitting;
		features.data[names[3]][index + i] = (float)smoking_years;
		features.data[names[4]][index + i] = (float)pack_years;

	}

	return 0;
}
