#include "AlcoholGenerator.h"
#include <boost/algorithm/string/predicate.hpp>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace boost;

void generateAlcoholRangeSignal(SDateVal2* rawSignal, SDateRangeVal *outRangeSignal) {

}

/*
void AlcoholGenerator::set_names() {
	names.clear();
	unordered_set<string> legal_features({ "Current_Drinker", "Drinking_Quantity", "PLM_Drinking_Level" });
	if (raw_feature_names.size() == 0)
		MTHROW_AND_ERR("AlcoholGenerator got no alcohol_features");
	for (string s : raw_feature_names) {
		if (legal_features.find(s) == legal_features.end())
			MTHROW_AND_ERR("AlcoholGenerator does not know how to generate [%s]", s.c_str());
		names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);
	}
}
*/

void AlcoholGenerator::set_names() {
	names.clear();
	unordered_set<string> legal_features({ "Current_Drinker", "Ex_Drinker", "Drinker_Years_Since_Quitting",
										   "Drinker_Years", "Drink_Unit_Years", "PLM_Drinking_Level","Never_Drinker",
										   "Unknown_Drinker", "Drinker_Quantity", "Current_Alcoholic","Ex_Alcoholic" });

	if (raw_feature_names.size() == 0)
		MTHROW_AND_ERR("AlcoholGenerator got no alcohol_features");
	for (string s : raw_feature_names) {
		if (legal_features.find(s) == legal_features.end())
			MTHROW_AND_ERR("AlcoholGenerator does not know how to generate [%s]", s.c_str());
		names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);
	}
}

int AlcoholGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [AlcoholGenerator::init]
		if (field == "alcohol_features")
			boost::split(raw_feature_names, entry.second, boost::is_any_of(","));
		else if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "future_ind")
			future_ind = entry.second;
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for AlcoholGenerator\n", field.c_str());
		//! [AlcoholGenerator::init]
	}
	set_names();
	return 0;
}



int AlcoholGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {


	int bdate_sid = rec.my_base_rep->sigs.sid(req_signals.back());
	int alcohol_sid = rec.my_base_rep->sigs.sid("Alcohol_quantity");
	for (int i = 0; i < num; i++) {

		float __current_drinker, __ex_drinker, __never_drinker, __unknown_drinker, __drinker_years, __drinker_quantity, __plm_drinker_level;
		int __current_alcoholist, __ex_alcoholist, __years_since_quitting;
		float __unit_years = (float)MED_MAT_MISSING_VALUE;
		int len;
		__current_drinker = __ex_drinker = __never_drinker = __unknown_drinker = __drinker_years = __drinker_quantity = __plm_drinker_level = (float)MED_MAT_MISSING_VALUE;
		__years_since_quitting = MED_MAT_MISSING_VALUE;
		__current_alcoholist = __ex_alcoholist = 0;
		int qa_print = 0;
		int plm_drinking_level = MED_MAT_MISSING_VALUE;

		SDateShort2 *drk_info = (SDateShort2 *)rec.get(alcohol_sid, i, len);

		if (len > 0) {

			if (qa_print == 1) {
				fprintf(stderr, "**************************************************\n");
				fprintf(stderr, "**************************** pid: %i  \n", rec.pid);
			}
			int bdate = medial::repository::get_value(rec, bdate_sid);
			assert(bdate != -1);
			int byear = bdate;
			if (req_signals.back() == "BDATE")
				byear = int(bdate / 10000);


			int MAX_TO_TRIM = 30;
			int MAX_TO_REMOVE = 150;

			int UNIT_SIZE = 1;
			int BACK_YEARS = 2;
			int AGE_START_DRINKING = 20;
			int DRINKING_QUANTITY_IMPUTE = 2;

			int date_age_may_start_drinking = (byear + AGE_START_DRINKING) * 10000 + 101;
			int start_age = med_time_converter.convert_times(MedTime::Date, MedTime::Days, date_age_may_start_drinking);

			int future_date = 20990101;
			if (future_ind == "0") {
				future_date = med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time);
			}


			//*********************************************** create ranges ***********************************************

			int first_drinker, last_drinker, first_ex_drinker, last_ex_drinker, first_never_drinker, last_never_drinker, first_special_group, last_special_group;
			first_drinker = last_drinker = first_ex_drinker = last_ex_drinker = first_never_drinker = last_never_drinker = first_special_group = last_special_group = -1;

			//# 1 ex # 2 current # 3 never # 4 current not ex unknown # 5 special group
			map <int, string> drink_desc;
			drink_desc[1] = "ex drinker"; drink_desc[2] = "current drinker"; drink_desc[3] = "never drinker"; drink_desc[4] = "current not ex unknown";

			vector<int> drinking_dates;
			vector<float> drinking_vals;
			vector<float> ex_drinking_vals;


			//####### read base signal
			for (int j = 0; j < len; j++) {
				int type = drk_info[j].val1;
				int quan = drk_info[j].val2;
				int date = drk_info[j].date;
				if (date == 0) continue;

				if (date > future_date) continue;

				if (quan > MAX_TO_REMOVE) quan = -9;
				if (quan > MAX_TO_TRIM) quan = MAX_TO_TRIM;

				if (qa_print == 1) fprintf(stderr, "%i %s %i \n", date, drink_desc[type].c_str(), quan);

				if (j == 0 && date > date_age_may_start_drinking) {
					int date_minus10y_n = med_time_converter.convert_times(MedTime::Date, MedTime::Days, date) - 365 * BACK_YEARS;
					int date_minus10y = med_time_converter.convert_times(MedTime::Days, MedTime::Date, date_minus10y_n);

					if (date_minus10y > date_age_may_start_drinking)
						drinking_dates.push_back(date_minus10y);
					else
						drinking_dates.push_back(date_age_may_start_drinking);
				}

				drinking_dates.push_back(date);
				if (type == 1) {
					if (first_ex_drinker == -1) first_ex_drinker = date;
					last_ex_drinker = date;
					if (quan != -9) ex_drinking_vals.push_back((float)quan);
				}
				else if (type == 2) {
					if (first_drinker == -1) first_drinker = date;
					last_drinker = date;
					if (quan != -9) drinking_vals.push_back((float)quan);
				}
				else if (type == 3 || type == 4) {
					if (first_never_drinker == -1) first_never_drinker = date;
					last_never_drinker = date;
				}
				else if (type == 5) {
					if (first_special_group == -1) first_special_group = date;
					last_special_group = date;
					//fprintf(stderr, "alcoholic : %i \n", date);
				}
			}

			//fprintf(stderr, "first_special_group %i \n", first_special_group);
			//fprintf(stderr, "last_special_group %i \n", last_special_group);

			drinking_dates.push_back(29990101);

			float drinking_quan, ex_drinking_quan;
			drinking_quan = ex_drinking_quan = missing_val;
			if (drinking_vals.size() > 0) drinking_quan = medial::stats::mean_without_cleaning(drinking_vals);
			if (ex_drinking_vals.size() > 0) ex_drinking_quan = medial::stats::mean_without_cleaning(ex_drinking_vals);


			//####### Aplly logic according to article "Development of an algorithm for determining smoking status and behavior ......
			//#date_group : 0 - no information , 1 - unknown , 2 - never drink , 3 - ex drinker , 4 - relapsed drinker , 6 - likely_drinker , 7 - drinker
			vector <pair<int, int>> drink_pre_ranges;
			for (int d = 0; d < drinking_dates.size(); d++) {
				int temp_date = drinking_dates[d];
				int date_group = 0;
				if (first_ex_drinker != -1 || first_never_drinker != -1) {
					if (first_drinker != -1 && first_drinker <= temp_date) {
						if (last_drinker != -1 && last_drinker <= temp_date) {
							if (first_never_drinker != -1 && first_never_drinker <= temp_date) {
								if (last_never_drinker <= last_drinker) {
									if (first_ex_drinker != -1 && first_ex_drinker <= temp_date) {
										if (last_ex_drinker <= last_drinker)
											date_group = 7;
										else
											date_group = 3;
									}
									else {
										date_group = 7;
									}
								}
								else {
									date_group = 3;
								}
							}
							else {
								if (first_ex_drinker != -1 && first_ex_drinker <= temp_date) {
									if (last_ex_drinker <= last_drinker)
										date_group = 7;
									else
										date_group = 3;
								}
								else {
									date_group = 7;
								}
							}
						}
						else {
							if (first_never_drinker != -1 && first_never_drinker <= temp_date) {
								if (last_never_drinker != -1 && last_never_drinker >= temp_date)
									date_group = 4;
								else {
									if (first_never_drinker != -1 && first_never_drinker <= temp_date) {
										if (first_ex_drinker != -1 && first_ex_drinker >= temp_date)
											date_group = 4;
										else
											date_group = 7;
									}
									else {
										date_group = 7;
									}
								}
							}
							else {
								if (first_ex_drinker != -1 && first_ex_drinker <= temp_date) {
									if (last_ex_drinker != -1 && last_ex_drinker >= temp_date)
										date_group = 4;
									else
										date_group = 7;
								}
								else {
									date_group = 7;
								}
							}
						}
					}
					else {
						if (first_ex_drinker != -1 && first_ex_drinker <= temp_date)
							date_group = 3;
						else {
							if (first_never_drinker != -1 && first_never_drinker <= temp_date)
								date_group = 2;
							else
								date_group = 1;
						}
					}
				}
				else {
					if (first_drinker != -1 && first_drinker <= temp_date)
						date_group = 7;
					else
						date_group = 6;
				}
				drink_pre_ranges.push_back(pair<int, int>(temp_date, date_group));
				//fprintf(stderr, "in logic : %i %i \n", temp_date, date_group);
			}

			vector<SDateRangeVal> drink_ranges;
			int start_date, end_date, pre_group;
			start_date = end_date = -1;
			pre_group = -1;
			map<int, string> final_desc;
			final_desc[0] = "no_information";
			final_desc[1] = "unknown";
			final_desc[2] = "never_drinker";
			final_desc[3] = "ex_drinker";
			final_desc[7] = "drinker";


			//date_group : 0 - no information , 1 - unknown , 2 - never smoke , 3 - ex smoker , 4 - relapsed smoker , 6 - likely_smoker , 7 - smoker
			int drink_flag = 0;
			int ex_drink_flag = 0;

			for (int kk = 0; kk < drink_pre_ranges.size(); kk++) {
				int temp_date = drink_pre_ranges[kk].first;
				int group = drink_pre_ranges[kk].second;

				if (group == 4 || group == 6) group = 7;
				if (group == 7) drink_flag = 1;
				if (group == 3) ex_drink_flag = 1;

				if (kk == 0) {
					start_date = temp_date;
					pre_group = group;
				}
				else if (kk < drink_pre_ranges.size() - 1) {
					if (group != pre_group) {
						SDateRangeVal temp;
						temp.date_start = start_date;
						temp.date_end = med_time_converter.convert_times(MedTime::Days, MedTime::Date, med_time_converter.convert_times(MedTime::Date, MedTime::Days, temp_date) - 1);
						temp.val = (float)pre_group;
						drink_ranges.push_back(temp);

						start_date = temp_date;
						pre_group = group;
					}
				}
				else {
					SDateRangeVal temp;
					temp.date_start = start_date;
					temp.date_end = temp_date;
					temp.val = (float)group;
					drink_ranges.push_back(temp);
				}
			}

			if (drinking_quan != missing_val && ex_drinking_quan == missing_val && ex_drink_flag == 1)
				ex_drinking_quan = drinking_quan;

			if (ex_drinking_quan != missing_val && drinking_quan == missing_val && drink_flag == 1)
				drinking_quan = ex_drinking_quan;

			if (qa_print == 1) {
				fprintf(stderr, "====ranges: \n");
				for (int kk = 0; kk < drink_ranges.size(); kk++) {
					if (qa_print == 1) fprintf(stderr, "date:start:%i date_end:%i val:%s \n", drink_ranges[kk].date_start, drink_ranges[kk].date_end, final_desc[(int)drink_ranges[kk].val].c_str());
				}
				fprintf(stderr, "\n");
			}

			//*********************************************** create vals for date ***********************************************

			int test_date = med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time);
			if (qa_print == 1) fprintf(stderr, "pid: %i len %i byear %i date %i  ************\n", rec.pid, len, byear, test_date);


			//***************** years drinking ***************
			int drinking_year = 0;
			int unit_years = 0;
			int any_drink_flag = 0;
			int any_info = 0;
			for (int kk = 0; kk < drink_ranges.size(); kk++) {
				int start_date = drink_ranges[kk].date_start;
				int end_date = drink_ranges[kk].date_end;

				if (test_date >= start_date) {

					any_info = 1;
					if (drink_ranges[kk].val == 7 || drink_ranges[kk].val == 3) any_drink_flag = 1;
					if (test_date < end_date)  end_date = test_date;   //for counting years
					int start_n = med_time_converter.convert_times(MedTime::Date, MedTime::Days, start_date);
					int end_n = med_time_converter.convert_times(MedTime::Date, MedTime::Days, end_date);
					if (qa_print == 1) fprintf(stderr, "years: %i %i  \n", start_date, end_date);

					if (drink_ranges[kk].val == 7) {

						int diff = end_n - start_n + 1;   //smoking
						drinking_year += diff;
						int temp_unit_years = 0;
						if (drinking_quan != missing_val) temp_unit_years += diff * ((int)drinking_quan);
						else temp_unit_years += diff * DRINKING_QUANTITY_IMPUTE;
						unit_years += temp_unit_years;
						if (qa_print == 1) fprintf(stderr, "years: smoker %f quantity  %f \n", (float)diff / 365, (float)temp_unit_years / (365 * UNIT_SIZE));
					}
					else if (drink_ranges[kk].val == 1 && drink_ranges[kk + 1].val == 7) {
						int diff = end_n - start_n;   //smoking

						drinking_year += diff;   // unknown before drinking
						if (diff >= 365 * BACK_YEARS) diff = 365 * BACK_YEARS;

						int temp_unit_years = 0;
						if (drinking_quan != missing_val) temp_unit_years += diff * ((int)drinking_quan);
						else temp_unit_years += diff * DRINKING_QUANTITY_IMPUTE;
						unit_years += temp_unit_years;
						if (qa_print == 1) fprintf(stderr, "years: unknown and drinker %f diff %f \n", (float)diff / 365, (float)temp_unit_years / (365 * UNIT_SIZE));
					}
					else if (drink_ranges[kk].val == 1 && drink_ranges[kk + 1].val == 3) {
						int diff = end_n - start_n;
						if (diff >= 365 * BACK_YEARS) diff = 365 * BACK_YEARS;

						if (diff > 0) {
							drinking_year += diff;   // unknown before smoking

							int temp_unit_years = 0;
							if (ex_drinking_quan != missing_val) temp_unit_years += diff * ((int)ex_drinking_quan);
							else temp_unit_years += diff * DRINKING_QUANTITY_IMPUTE;
							unit_years += temp_unit_years;
							if (qa_print == 1) fprintf(stderr, "years: unknown and ex drinker %f diff %f \n", (float)diff / 365, (float)temp_unit_years / (365 * UNIT_SIZE));
						}
						else {
							//fprintf(stderr, "years: unknown and ex smoker ______ \n");
						}
					}
					else if ((drink_ranges[kk].val == 3 && drink_ranges.size() == 1) || (drink_ranges[kk].val == 3 && (drink_ranges[kk - 1].val == 2 || drink_ranges[kk - 1].val == 1) && drink_ranges.size() <= 3)) {
						int diff = start_n - start_age;   //smoking
						if (diff >= 365 * BACK_YEARS) diff = 365 * BACK_YEARS;

						if (diff > 0) {
							drinking_year += diff;   // unknown before smoking

							int temp_unit_years = 0;
							if (drinking_quan != missing_val) temp_unit_years += diff * ((int)drinking_quan);
							else temp_unit_years += diff * DRINKING_QUANTITY_IMPUTE;
							unit_years += temp_unit_years;
							if (qa_print == 1) fprintf(stderr, "years: unknown and drinker %f diff %f \n", (float)diff / 365, (float)temp_unit_years / (365 * UNIT_SIZE));
						}
						else {
							drinking_year += 1;
							unit_years += 365 * DRINKING_QUANTITY_IMPUTE;
						}
					}
				}
			}


			float drinking_year_f = (float)drinking_year / 365;
			float unit_years_f = (float)unit_years / 365 / UNIT_SIZE;

			if (any_drink_flag == 1 && drinking_year_f == 0) {
				drinking_year_f = missing_val;
				unit_years_f = missing_val;
			}

			if (drinking_year_f != missing_val && drinking_year_f > 0 && drinking_year_f < 1) drinking_year_f = 1;
			float drinking_year_ff = round(drinking_year_f);
			float unit_years_ff = (float)(round(unit_years_f * 100)) / 100;



			//******************** current state and years since

			if (qa_print == 1) fprintf(stderr, "\n");
			//date_group : 0 - no information , 1 - unknown , 2 - never smoke , 3 - ex smoker , 4 - relapsed smoker , 6 - likely_smoker , 7 - smoker
			float years_since_q_f = 0;
			float never_drink_f = 0;
			float ex_drink_f = 0;
			float drink_f = 0;
			float unknown_f = 0;
			float quantity_f = 0;
			for (int kk = (int)drink_ranges.size() - 1; kk >= 0; kk--) {
				int start_date = drink_ranges[kk].date_start;
				int end_date = drink_ranges[kk].date_end;

				if (test_date >= start_date && test_date <= end_date) {
					int val = (int)drink_ranges[kk].val;
					//fprintf(stderr, "val : %i \n", val);
					if (val == 1) {
						unknown_f = 1;
						quantity_f = missing_val;
					}
					else if (val == 2) never_drink_f = 1;
					else if (val == 3) {
						ex_drink_f = 1;
						for (int jj = kk - 1; jj >= 0; jj--) {
							if (drink_ranges[jj].val == 7 || drink_ranges[jj].val == 1) {
								int drink_end_date = drink_ranges[jj].date_end;
								int diff_days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, test_date) - med_time_converter.convert_times(MedTime::Date, MedTime::Days, drink_end_date);
								years_since_q_f = (float)diff_days / 365;
								if (years_since_q_f < 1) years_since_q_f = 1;
								years_since_q_f = round(years_since_q_f);
								break;
							}
						}
						quantity_f = ex_drinking_quan;
					}
					else if (val == 7) {
						drink_f = 1;
						quantity_f = drinking_quan;
					}
					break;
				}
			}

			//********************* alcoholist

			//fprintf(stderr, "first_special_group %i \n", first_special_group);
			//fprintf(stderr, "last_special_group %i \n", last_special_group);
			//fprintf(stderr, "__current_alcoholist %i \n", __current_alcoholist);

			if (first_special_group != -1) {
				if (test_date >= first_special_group && test_date <= last_special_group) __current_alcoholist = 1;
				else if (test_date > last_special_group) __ex_alcoholist = 1;
			}

			//fprintf(stderr, "__current_alcoholist %i \n", __current_alcoholist);

			//***************************************  final ***************************************

			if (any_info == 1) {
				__current_drinker = drink_f;
				__ex_drinker = ex_drink_f;
				__never_drinker = never_drink_f;
				__unknown_drinker = unknown_f;
				__years_since_quitting = (int)years_since_q_f;
				__drinker_years = drinking_year_ff;
				__drinker_quantity = quantity_f;
				__unit_years = unit_years_ff;



				if (qa_print == 1) {
					fprintf(stderr, "====== final : \n");
					fprintf(stderr, "__years_since_quitting %i \n", __years_since_quitting);
					fprintf(stderr, "__never_drinker %f \n", __never_drinker);
					fprintf(stderr, "__ex_drinker %f \n", __ex_drinker);
					fprintf(stderr, "__current_drinker %f \n", __current_drinker);
					fprintf(stderr, "__unknown_drinker %f \n", __unknown_drinker);
					fprintf(stderr, "__drinking_quantity  %f \n", __drinker_quantity);
					fprintf(stderr, "__drinking_years: %f \n", __drinker_years);
					fprintf(stderr, "__unit_years : %f \n", __unit_years);
					fprintf(stderr, "__current_alcoholist : %i \n", __current_alcoholist);
					fprintf(stderr, "__ex_alcoholist : %i \n", __ex_alcoholist);
				}



				//************************************************** PLM ***********************************************


				int current_drinker = (int)__current_drinker;
				int drinking_level = (int)__drinker_quantity;
				if (current_drinker == 1 && drinking_level == -1)
					plm_drinking_level = 1;
				else if (drinking_level == 0)
					plm_drinking_level = 0;
				else if (drinking_level <= 7)
					plm_drinking_level = 1;
				else if (drinking_level <= 14)
					plm_drinking_level = 2;
				else if (drinking_level <= 21)
					plm_drinking_level = 3;
				else if (drinking_level <= 28)
					plm_drinking_level = 4;
				else
					plm_drinking_level = 5;
				if (__current_alcoholist)
					plm_drinking_level = 5;
			}
			//************************************************** Add to matrix  *********************************************
		}
		// Current Drinker
		if (_p_data[ALC_CURRENT_DRINKER] != NULL) _p_data[ALC_CURRENT_DRINKER][index + i] = (float)__current_drinker;
		// Ex Drinker
		if (_p_data[ALC_EX_DRINKER] != NULL) _p_data[ALC_EX_DRINKER][index + i] = (float)__ex_drinker;
		// Drinker_Years_Since_Quitting
		if (_p_data[ALC_DRINKER_YEARS_SINCE_QUITTING] != NULL) _p_data[ALC_DRINKER_YEARS_SINCE_QUITTING][index + i] = (float)__years_since_quitting;
		// Drinker_Years
		if (_p_data[ALC_DRINKING_YEARS] != NULL) _p_data[ALC_DRINKING_YEARS][index + i] = (float)__drinker_years;
		// Drink_Unit_Years
		if (_p_data[ALC_DRINKING_UNIT_YEARS] != NULL) _p_data[ALC_DRINKING_UNIT_YEARS][index + i] = (float)__unit_years;
		// PLM_Drinking_Level
		if (_p_data[ALC_PLM_DRINKING_LEVEL] != NULL) _p_data[ALC_PLM_DRINKING_LEVEL][index + i] = (float)plm_drinking_level;
		// Never_Drinker
		if (_p_data[ALC_NEVER_DRINKER] != NULL) _p_data[ALC_NEVER_DRINKER][index + i] = (float)__never_drinker;
		// Unknown_Drinker
		if (_p_data[ALC_UNKNOWN_DRINKER] != NULL) _p_data[ALC_UNKNOWN_DRINKER][index + i] = (float)__unknown_drinker;
		// Drinker_Quantity
		if (_p_data[ALC_DRINKER_QUANTITY] != NULL) _p_data[ALC_DRINKER_QUANTITY][index + i] = (float)__drinker_quantity;
		// Current_Alcoholic
		if (_p_data[ALC_CURRENT_ALCOHOLIC] != NULL) _p_data[ALC_CURRENT_ALCOHOLIC][index + i] = (float)__current_alcoholist;
		// Ex_Alcoholic
		if (_p_data[ALC_EX_ALCOHOLIC] != NULL) _p_data[ALC_EX_ALCOHOLIC][index + i] = (float)__ex_alcoholist;
	}

	return 0;
}

void AlcoholGenerator::get_p_data(MedFeatures& features, vector<float *> &_p_data) {

	p_data.resize(ALC_LAST, NULL);

	if (iGenerateWeights) {
		if (names.size() != 1)
			MTHROW_AND_ERR("Cannot generate weights using a multi-feature generator (type %d generates %d features)\n", generator_type, (int)names.size())
		else
			p_data[0] = &(features.weights[0]);
	}

	for (string &name : names) {
		if (algorithm::ends_with(name, "Current_Drinker"))
			_p_data[ALC_CURRENT_DRINKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Ex_Drinker"))
			_p_data[ALC_EX_DRINKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Drinker_Years_Since_Quitting"))
			_p_data[ALC_DRINKER_YEARS_SINCE_QUITTING] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Drinker_Years"))
			_p_data[ALC_DRINKING_YEARS] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Drink_Unit_Years"))
			_p_data[ALC_DRINKING_UNIT_YEARS] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "PLM_Drinking_Level"))
			_p_data[ALC_PLM_DRINKING_LEVEL] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Never_Drinker"))
			_p_data[ALC_NEVER_DRINKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Unknown_Drinker"))
			_p_data[ALC_UNKNOWN_DRINKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Drinker_Quantity"))
			_p_data[ALC_DRINKER_QUANTITY] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Current_Alcoholic"))
			_p_data[ALC_CURRENT_ALCOHOLIC] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Ex_Alcoholic"))
			_p_data[ALC_EX_ALCOHOLIC] = &(features.data[name][0]);
		else
			MTHROW_AND_ERR("unknown feature name [%s]", name.c_str());
	}
}