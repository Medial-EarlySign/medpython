#include "SmokingGenerator.h"
#include <boost/algorithm/string/predicate.hpp>

#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL
using namespace boost;

//iterate on id
//void generateSmokingRangeSignal(SDateVal2* rawSignal, SDateRangeVal *outRangeSignal) {
//}

void SmokingGenerator::set_names() {
	names.clear();
	unordered_set<string> legal_features({ "Current_Smoker", "Ex_Smoker", "Smok_Years_Since_Quitting", "Smoking_Years", "Smok_Pack_Years", "PLM_Smoking_Level","Never_Smoker", "Unknown_Smoker", "Smoking_Quantity" });

	if (raw_feature_names.size() == 0)
		MTHROW_AND_ERR("SmokingGenerator got no smoking_features");
	for (string s : raw_feature_names) {
		if (legal_features.find(s) == legal_features.end())
			MTHROW_AND_ERR("SmokingGenerator does not know how to generate [%s]", s.c_str());
		names.push_back("FTR_" + int_to_string_digits(serial_id, 6) + "." + s);
		tags.push_back(s);
	}
}

int SmokingGenerator::init(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [SmokingGenerator::init]
		if (field == "smoking_features")
			boost::split(raw_feature_names, entry.second, boost::is_any_of(","));
		else if (field == "smoking_method")
			smoking_method = entry.second;
		else if (field == "future_ind")
			future_ind = entry.second;
		else if (field == "tags")
			boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field == "weights_generator")
			iGenerateWeights = med_stoi(entry.second);
		else if (field != "fg_type")
			MTHROW_AND_ERR("Unknown parameter \'%s\' for SmokingGenerator\n", field.c_str());
		//! [SmokingGenerator::init]

	}
	set_names();
	req_signals.clear();
	req_signals.push_back(smoking_method);
	req_signals.push_back("BDATE");
	return 0;
}



int SmokingGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data) {

	int missing = MED_MAT_MISSING_VALUE;


	for (int i = 0; i < num; i++) {

		int __current_smoker, __ex_smoker, __never_smoker, __unknown_smoker, __years_since_quitting, __smoking_years, __smoking_quantity, __plm_smoking_level;
		float __pack_years = (float)MED_MAT_MISSING_VALUE;
		int len;
		__current_smoker = __ex_smoker = __never_smoker = __unknown_smoker = __years_since_quitting = __smoking_years = __smoking_quantity = __plm_smoking_level = MED_MAT_MISSING_VALUE;
		int qa_print = 0;

		if (smoking_method != "SMOKING_ENRICHED") {
			string sname = "Smoking_quantity";
			SDateShort2 *smx_info = (SDateShort2 *)rec.get(sname, i, len);

			if (len > 0) {

				if (qa_print == 1) fprintf(stderr, "pid: %i  \n", rec.pid);

				sname = "BDATE";
				int bdateId = rec.my_base_rep->sigs.sid(sname);
				int bdate = medial::repository::get_value(rec, bdateId);
				int byear = int(bdate / 10000);
				assert(byear != -1);

				int MAX_TO_TRIM = 60;
				int MAX_TO_REMOVE = 100;

				int AGE_START_SMOKING = 20;
				int EX_SMOKING_YEARS = 10;
				int SMOKING_QUANTITY_IMPUTE = 10;
				int PACK_SIZE = 20;

				int date_age_may_start_smoking = (byear + AGE_START_SMOKING) * 10000 + 101;
				int start_age = med_time_converter.convert_times(MedTime::Date, MedTime::Days, date_age_may_start_smoking);

				int future_date = 20990101;
				if (future_ind == "0") {
					future_date = med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time);
				}

				//*********************************************** create ranges ***********************************************

				int first_smoker, last_smoker, first_ex_smoker, last_ex_smoker, first_never_smoker, last_never_smoker, special_group;
				first_smoker = last_smoker = first_ex_smoker = last_ex_smoker = first_never_smoker = last_never_smoker = -1, special_group = -1;

				//# 1 ex
				//# 2 current
				//# 3 never
				//# 4 current not ex unknown
				//# 5 special group

				map <int, string> smoke_desc;
				smoke_desc[1] = "ex smoker";
				smoke_desc[2] = "current smoker";
				smoke_desc[3] = "never smoker";
				smoke_desc[4] = "current not ex unknown";




				vector<int>smoking_dates;
				vector<float> smoking_vals;
				vector<float> ex_smoking_vals;

				//fprintf(stderr, "====med_code: \n");
				for (int j = 0; j < len; j++) {
					int type = smx_info[j].val1;
					int quan = smx_info[j].val2;
					int date = smx_info[j].date;
					if (date == 0) continue;

					if (date > future_date) continue;

					if (quan > MAX_TO_REMOVE) quan = -9;
					if (quan > MAX_TO_TRIM) quan = MAX_TO_TRIM;

					if (j == 0 && date > date_age_may_start_smoking)
						smoking_dates.push_back(date_age_may_start_smoking);

					if (qa_print == 1) fprintf(stderr, "%i %s %i \n", date, smoke_desc[type].c_str(), quan);
					smoking_dates.push_back(date);

					if (type == 1) {
						if (first_ex_smoker == -1) first_ex_smoker = date;
						last_ex_smoker = date;
						if (quan != -9) ex_smoking_vals.push_back((float)quan);
					}
					else if (type == 2) {
						if (first_smoker == -1) first_smoker = date;
						last_smoker = date;
						if (quan != -9) smoking_vals.push_back((float)quan);
					}
					else if (type == 3 || type == 4) {
						if (first_never_smoker == -1) first_never_smoker = date;
						last_never_smoker = date;
					}
					else if (type == 5 && special_group == -1) special_group = date;
				}
				smoking_dates.push_back(29990101);

				float smoking_quan, ex_smoking_quan;
				smoking_quan = ex_smoking_quan = missing_val;
				if (smoking_vals.size() > 0) smoking_quan = medial::stats::mean_without_cleaning(smoking_vals);
				if (ex_smoking_vals.size() > 0) ex_smoking_quan = medial::stats::mean_without_cleaning(ex_smoking_vals);

				//date_group : 0 - no information , 1 - unknown , 2 - never smoke , 3 - ex smoker , 4 - relapsed smoker , 6 - likely_smoker , 7 - smoker
				vector <pair<int, int>> smoke_pre_ranges;
				for (int d = 0; d < smoking_dates.size(); d++) {
					int temp_date = smoking_dates[d];
					int date_group = 0;
					if (first_ex_smoker != -1 || first_never_smoker != -1) {
						if (first_smoker != -1 && first_smoker <= temp_date) {
							if (last_smoker != -1 && last_smoker <= temp_date) {
								if (first_never_smoker != -1 && first_never_smoker <= temp_date) {
									if (last_never_smoker <= last_smoker) {
										if (first_ex_smoker != -1 && first_ex_smoker <= temp_date) {
											if (last_ex_smoker <= last_smoker)
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
									if (first_ex_smoker != -1 && first_ex_smoker <= temp_date) {
										if (last_ex_smoker <= last_smoker)
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
								if (first_never_smoker != -1 && first_never_smoker <= temp_date) {
									if (last_never_smoker != -1 && last_never_smoker >= temp_date)
										date_group = 4;
									else {
										if (first_never_smoker != -1 && first_never_smoker <= temp_date) {
											if (first_ex_smoker != -1 && first_ex_smoker >= temp_date)
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
									if (first_ex_smoker != -1 && first_ex_smoker <= temp_date) {
										if (last_ex_smoker != -1 && last_ex_smoker >= temp_date)
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
							if (first_ex_smoker != -1 && first_ex_smoker <= temp_date)
								date_group = 3;
							else {
								if (first_never_smoker != -1 && first_never_smoker <= temp_date)
									date_group = 2;
								else
									date_group = 1;
							}
						}
					}
					else {
						if (first_smoker != -1 && first_smoker <= temp_date)
							date_group = 7;
						else
							date_group = 6;
					}
					smoke_pre_ranges.push_back(pair<int, int>(temp_date, date_group));
					//fprintf(stderr, "in logic : %i %i \n", temp_date, date_group);
				}

				vector<SDateRangeVal> smoke_ranges;
				int start_date, end_date, pre_group;
				start_date = end_date = -1;
				pre_group = -1;
				map<int, string> final_desc;
				final_desc[0] = "no_information";
				final_desc[1] = "unknown";
				final_desc[2] = "never_smoker";
				final_desc[3] = "ex_smoker";
				final_desc[7] = "smoker";


				//date_group : 0 - no information , 1 - unknown , 2 - never smoke , 3 - ex smoker , 4 - relapsed smoker , 6 - likely_smoker , 7 - smoker
				int smoke_flag = 0;
				int ex_smoke_flag = 0;

				for (int kk = 0; kk < smoke_pre_ranges.size(); kk++) {
					int temp_date = smoke_pre_ranges[kk].first;
					int group = smoke_pre_ranges[kk].second;

					if (group == 4 || group == 6) group = 7;
					if (group == 7) smoke_flag = 1;
					if (group == 3) ex_smoke_flag = 1;

					if (kk == 0) {
						start_date = temp_date;
						pre_group = group;
					}
					else if (kk < smoke_pre_ranges.size() - 1) {
						if (group != pre_group) {
							SDateRangeVal temp;
							temp.date_start = start_date;
							temp.date_end = med_time_converter.convert_times(MedTime::Days, MedTime::Date, med_time_converter.convert_times(MedTime::Date, MedTime::Days, temp_date) - 1);
							temp.val = (float)pre_group;
							smoke_ranges.push_back(temp);

							start_date = temp_date;
							pre_group = group;
						}
					}
					else {
						SDateRangeVal temp;
						temp.date_start = start_date;
						temp.date_end = temp_date;
						temp.val = (float)group;
						smoke_ranges.push_back(temp);
					}
				}

				if (smoking_quan != missing_val && ex_smoking_quan == missing_val && ex_smoke_flag == 1)
					ex_smoking_quan = smoking_quan;

				if (ex_smoking_quan != missing_val && smoking_quan == missing_val && smoke_flag == 1)
					smoking_quan = ex_smoking_quan;


				if (qa_print == 1) {
					fprintf(stderr, "====ranges: \n");
					for (int kk = 0; kk < smoke_ranges.size(); kk++) {
						if (qa_print == 1) fprintf(stderr, "date:start:%i date_end:%i val:%s \n", smoke_ranges[kk].date_start, smoke_ranges[kk].date_end, final_desc[(int)smoke_ranges[kk].val].c_str());
					}
					fprintf(stderr, "\n");
				}


				//*********************************************** create vals for date ***********************************************

				int test_date = med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time);
				if (qa_print == 1) fprintf(stderr, "********** pid: %i len %i byear %i date %i  ************\n", rec.pid, len, byear, test_date);

				int smoking_year = 0;
				int pack_years = 0;
				int any_smoke_flag = 0;
				int any_info = 0;
				for (int kk = 0; kk < smoke_ranges.size(); kk++) {
					int start_date = smoke_ranges[kk].date_start;
					int end_date = smoke_ranges[kk].date_end;

					if (start_date < date_age_may_start_smoking)
						start_date = date_age_may_start_smoking;

					if (end_date < start_date)
						end_date = start_date;



					if (test_date >= start_date) {
						any_info = 1;

						if (smoke_ranges[kk].val == 7 || smoke_ranges[kk].val == 3)
							any_smoke_flag = 1;

						if (test_date < end_date)  end_date = test_date;

						int start_n = med_time_converter.convert_times(MedTime::Date, MedTime::Days, start_date);
						int end_n = med_time_converter.convert_times(MedTime::Date, MedTime::Days, end_date);

						if (qa_print == 1) fprintf(stderr, "years: %i %i  \n", start_date, end_date);

						if (smoke_ranges[kk].val == 7) {

							int diff = end_n - start_n + 1;   //smoking
							smoking_year += diff;
							int temp_pack_years = 0;
							if (smoking_quan != missing_val) temp_pack_years += diff * ((int)smoking_quan);
							else temp_pack_years += diff * SMOKING_QUANTITY_IMPUTE;
							pack_years += temp_pack_years;
							if (qa_print == 1) fprintf(stderr, "years: smoker %f quantity  %f \n", (float)diff / 365, (float)temp_pack_years / (365 * PACK_SIZE));
						}
						else if (smoke_ranges[kk].val == 1 && smoke_ranges[kk + 1].val == 7) {
							int diff = end_n - start_n;   //smoking
							smoking_year += diff;   // unknown before smoking

							int temp_pack_years = 0;
							if (smoking_quan != missing_val) temp_pack_years += diff * ((int)smoking_quan);
							else temp_pack_years += diff * SMOKING_QUANTITY_IMPUTE;
							pack_years += temp_pack_years;
							if (qa_print == 1) fprintf(stderr, "years: unknown and smoker %f diff %f \n", (float)diff / 365, (float)temp_pack_years / (365 * PACK_SIZE));
						}
						else if (smoke_ranges[kk].val == 1 && smoke_ranges[kk + 1].val == 3) {
							int diff = end_n - start_n;
							diff -= 365 * EX_SMOKING_YEARS;

							if (diff > 0) {
								smoking_year += diff;   // unknown before smoking

								int temp_pack_years = 0;
								if (ex_smoking_quan != missing_val) temp_pack_years += diff * ((int)ex_smoking_quan);
								else temp_pack_years += diff * SMOKING_QUANTITY_IMPUTE;
								pack_years += temp_pack_years;
								if (qa_print == 1) fprintf(stderr, "years: unknown and ex smoker %f diff %f \n", (float)diff / 365, (float)temp_pack_years / (365 * PACK_SIZE));
							}
							else {
								//fprintf(stderr, "years: unknown and ex smoker ______ \n");
							}
						}
						else if ((smoke_ranges[kk].val == 3 && smoke_ranges.size() == 1) || (smoke_ranges[kk].val == 3 && (smoke_ranges[kk - 1].val == 2 || smoke_ranges[kk - 1].val == 1) && smoke_ranges.size() <= 3)) {
							int diff = start_n - 365 * EX_SMOKING_YEARS - start_age;   //smoking

							if (diff > 0) {
								smoking_year += diff;   // unknown before smoking

								int temp_pack_years = 0;
								if (smoking_quan != missing_val) temp_pack_years += diff * (int)smoking_quan;
								else temp_pack_years += diff * SMOKING_QUANTITY_IMPUTE;
								pack_years += temp_pack_years;
								if (qa_print == 1) fprintf(stderr, "years: unknown and smoker %f diff %f \n", (float)diff / 365, (float)temp_pack_years / (365 * PACK_SIZE));
							}
							else {
								smoking_year += 1;
								pack_years += 365 * SMOKING_QUANTITY_IMPUTE;
							}
						}
					}
				}


				float smoking_year_f = (float)smoking_year / 365;
				float pack_years_f = (float)pack_years / 365 / PACK_SIZE;

				if (any_smoke_flag == 1 && smoking_year_f == 0) {
					smoking_year_f = missing_val;
					pack_years_f = missing_val;
				}

				if (smoking_year_f != missing_val && smoking_year_f > 0 && smoking_year_f < 1) smoking_year_f = 1;

				float smoking_year_ff = round(smoking_year_f);
				float pack_years_ff = (float)(round(pack_years_f * 100)) / 100;


				if (qa_print == 1) fprintf(stderr, "\n");
				//date_group : 0 - no information , 1 - unknown , 2 - never smoke , 3 - ex smoker , 4 - relapsed smoker , 6 - likely_smoker , 7 - smoker
				float years_since_q_f = 0;
				float never_smoke_f = 0;
				float ex_smoke_f = 0;
				float smoke_f = 0;
				float unknown_f = 0;
				float quantity_f = 0;
				for (int kk = (int)smoke_ranges.size() - 1; kk >= 0; kk--) {
					int start_date = smoke_ranges[kk].date_start;
					int end_date = smoke_ranges[kk].date_end;

					if (test_date >= start_date && test_date <= end_date) {
						int val = (int)smoke_ranges[kk].val;
						if (val == 1) {
							unknown_f = 1;
							quantity_f = missing_val;
						}
						else if (val == 2) never_smoke_f = 1;
						else if (val == 3) {
							ex_smoke_f = 1;
							for (int jj = kk - 1; jj >= 0; jj--) {
								if (smoke_ranges[jj].val == 7 || smoke_ranges[jj].val == 1) {
									int smoke_end_date = smoke_ranges[jj].date_end;
									int diff_days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, test_date) - med_time_converter.convert_times(MedTime::Date, MedTime::Days, smoke_end_date);
									years_since_q_f = (float)diff_days / 365;
									if (years_since_q_f < 1) years_since_q_f = 1;
									years_since_q_f = round(years_since_q_f);
									break;
								}
							}
							quantity_f = ex_smoking_quan;
						}
						else if (val == 7) {
							smoke_f = 1;
							quantity_f = smoking_quan;
						}
						break;
					}
				}


				if (any_info == 1) {
					__current_smoker = (int)smoke_f;
					__ex_smoker = (int)ex_smoke_f;
					__never_smoker = (int)never_smoke_f;
					__unknown_smoker = (int)unknown_f;
					__years_since_quitting = (int)years_since_q_f;
					__smoking_years = (int)smoking_year_ff;
					__smoking_quantity = (int)quantity_f;
					__pack_years = pack_years_ff;
				}

				if (qa_print == 1) {
					fprintf(stderr, "====== final : \n");
					fprintf(stderr, "years_since_q_f %i \n", __years_since_quitting);
					fprintf(stderr, "never_smoke_f %i \n", __never_smoker);
					fprintf(stderr, "ex_smoke_f %i \n", __ex_smoker);
					fprintf(stderr, "smoke_f %i \n", __current_smoker);
					fprintf(stderr, "unknown_f %i \n", __unknown_smoker);
					fprintf(stderr, "quantity_f  %i \n", __smoking_quantity);
					fprintf(stderr, "smoking_year_f: %i \n", __smoking_years);
					fprintf(stderr, "pack_years_f : %f \n", __pack_years);
				}

			}
		}
		else {

			int len1;
			bool never_smoked = true;

			string sname = "SMOKING_ENRICHED";
			SValShort4 *smx_status = (SValShort4 *)rec.get(sname, i, len1);
			if (len1 > 0)
				never_smoked = (smx_status[0].val1 == -1); // never smoked is -1, -1, 0, 0
			assert(len1 <= 1);
			__never_smoker = (int)never_smoked;


			if (len1 == 0) { // No Data
				__current_smoker = __ex_smoker = (int)missing;
				__years_since_quitting = __smoking_years = (int)missing;
				__pack_years = (float)missing;
			}
			else if (never_smoked) { // Non Smoker
				__current_smoker = __ex_smoker = 0;
				__years_since_quitting = 100;
				__smoking_years = 0;
				__pack_years = 0.0;
			}
			else { // (Ex)Smoker
				int start_year = smx_status[0].val1;
				int end_year = smx_status[0].val2;
				int target_year = (int)(med_time_converter.convert_times(features.time_unit, MedTime::Date, features.samples[index + i].time) / 10000);
				if (target_year < end_year) {
					// still in smoking period
					__smoking_years = target_year - start_year;
					__years_since_quitting = 0;
					__current_smoker = 1;
				}
				else {
					// maybe done smoking
					__current_smoker = smx_status[0].val4;
					__smoking_years = end_year - start_year; // we are merciful
					if (!__current_smoker)
						__years_since_quitting = target_year - end_year;
					else
						__years_since_quitting = 0;
				}

				__pack_years = ((float)smx_status[0].val3 / 20) * (float)__smoking_years;
				__ex_smoker = 1 - __current_smoker;
			}

			if (qa_print == 1) {
				fprintf(stderr, "====ido \n");
				fprintf(stderr, "current_smoker : %i \n", __current_smoker);
				fprintf(stderr, "ex_smoker : %i \n", __ex_smoker);
				fprintf(stderr, "years_since_quitting : %i \n", __years_since_quitting);
				fprintf(stderr, "smoking_years : %i \n", __smoking_years);
				fprintf(stderr, "pack_years : %f \n", __pack_years);
				fprintf(stderr, "__never_smoker : %i \n", __never_smoker);
				fprintf(stderr, "\n\n");
			}
		}


		//************************************************** PLM ***********************************************
		int plm_smoking_level = missing;

		if (__current_smoker == missing)
			plm_smoking_level = missing;
		else if (__never_smoker)
			plm_smoking_level = 0;
		else if (__ex_smoker == 1) {
			if (__years_since_quitting > 5)
				plm_smoking_level = 1;
			else if (__years_since_quitting <= 5)
				plm_smoking_level = 2;
		}
		else if (__current_smoker == 1) {
			float packs_per_day = 0;
			if (__smoking_years > 0)
				packs_per_day = __pack_years / __smoking_years; //avoid dividing by zero!
			if (packs_per_day <= 0.25)
				plm_smoking_level = 3;
			else if (packs_per_day > 0.25 && packs_per_day <= 0.5)
				plm_smoking_level = 4;
			else if (packs_per_day > 0.5 && packs_per_day <= 1)
				plm_smoking_level = 5;
			else if (packs_per_day > 1)
				plm_smoking_level = 6;
		}

		//************************************************** Add to matrix  *********************************************

		// Current_Smoker
		if (_p_data[SMX_CURRENT_SMOKER] != NULL) _p_data[SMX_CURRENT_SMOKER][index + i] = (float)__current_smoker;
		// Ex_Smoker
		if (_p_data[SMX_EX_SMOKER] != NULL) _p_data[SMX_EX_SMOKER][index + i] = (float)__ex_smoker;
		// Smok_Years_Since_Quitting
		if (_p_data[SMX_YEARS_SINCE_QUITTING] != NULL) _p_data[SMX_YEARS_SINCE_QUITTING][index + i] = (float)__years_since_quitting;
		// Smoking_Years
		if (_p_data[SMX_SMOKING_YEARS] != NULL) _p_data[SMX_SMOKING_YEARS][index + i] = (float)__smoking_years;
		// Smok_Pack_Years
		if (_p_data[SMX_SMOK_PACK_YEARS] != NULL) _p_data[SMX_SMOK_PACK_YEARS][index + i] = (float)__pack_years;
		// PLM_Smoking_Level
		if (_p_data[SMX_PLM_SMOKING_LEVEL] != NULL) _p_data[SMX_PLM_SMOKING_LEVEL][index + i] = (float)plm_smoking_level;
		// Never_Smoker
		if (_p_data[SMX_NEVER_SMOKER] != NULL) _p_data[SMX_NEVER_SMOKER][index + i] = (float)__never_smoker;
		// Unknown_Smoker
		if (_p_data[SMX_UNKNOWN_SMOKER] != NULL) _p_data[SMX_UNKNOWN_SMOKER][index + i] = (float)__unknown_smoker;
		// Smoking_Quantity
		if (_p_data[SMX_SMOKING_QUANTITY] != NULL) _p_data[SMX_SMOKING_QUANTITY][index + i] = (float)__smoking_quantity;
	}
	return 0;
}

// Get pointers to data vectors
//.......................................................................................
void SmokingGenerator::get_p_data(MedFeatures& features, vector<float *> &_p_data) {

	p_data.resize(SMX_LAST, NULL);

	if (iGenerateWeights) {
		if (names.size() != 1)
			MTHROW_AND_ERR("Cannot generate weights using a multi-feature generator (type %d generates %d features)\n", generator_type, (int)names.size())
		else
			p_data[0] = &(features.weights[0]);
	}

	for (string &name : names) {
		if (algorithm::ends_with(name, "Current_Smoker"))
			_p_data[SMX_CURRENT_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Ex_Smoker"))
			_p_data[SMX_EX_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smok_Years_Since_Quitting"))
			_p_data[SMX_YEARS_SINCE_QUITTING] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smoking_Years"))
			_p_data[SMX_SMOKING_YEARS] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smok_Pack_Years"))
			_p_data[SMX_SMOK_PACK_YEARS] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "PLM_Smoking_Level"))
			_p_data[SMX_PLM_SMOKING_LEVEL] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Never_Smoker"))
			_p_data[SMX_NEVER_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Unknown_Smoker"))
			_p_data[SMX_UNKNOWN_SMOKER] = &(features.data[name][0]);
		else if (algorithm::ends_with(name, "Smoking_Quantity"))
			_p_data[SMX_SMOKING_QUANTITY] = &(features.data[name][0]);
		else
			MTHROW_AND_ERR("unknown feature name [%s]", name.c_str());
	}
}
