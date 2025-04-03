#include "RepClearSignalByDiag.h"
#include <MedUtils/MedUtils/MedUtils.h>
#include <MedUtils/MedUtils/MedMedical.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL LOG_DEF_LEVEL

int RepClearSignalByDiag::init(map<string, string> &mapper)
{

	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RepClearSignalByDiag::init]
		if (it->first == "time_window")
			time_window = med_stoi(it->second);
		else if (it->first == "max_exclusion")
			max_exclusion = med_stoi(it->second);
		else if (it->first == "diag_file")
			boost::split(diag_list, it->second, boost::is_any_of(","));
		else if (it->first == "rp_type")
		{
		}
		else
			MTHROW_AND_ERR("Unknown parameter \'%s\' for RepFilterByDiag\n", it->first.c_str())
		//! [RepClearSignalByDiag::init]
	}
	cout << "max_exclusion " << endl;
	cout << "RepClearSignalByDiag CODES:" << endl;
	cout << "###########################" << endl;
	for (int i = 0; i < diag_list.size(); ++i)
		cout << i << " " << diag_list[i] << endl;

	req_signals.clear();
	req_signals.insert(signal_name);
	req_signals.insert("BDATE");
	req_signals.insert("GENDER");
	req_signals.insert("DIAGNOSIS");
	aff_signals.clear();
	aff_signals.insert(signal_name);

	return 0;
}

void RepClearSignalByDiag::init_tables(MedDictionarySections &dict, MedSignals &sigs)
{
	int section_id = dict.section_id("DIAGNOSIS");
	dict.prep_sets_lookup_table(section_id, diag_list, lut_censor);

	sig_id = sigs.sid(signal_name);
	bdate_id = sigs.sid("BDATE");
	gender_id = sigs.sid("GENDER");
	diag_id = sigs.sid("DIAGNOSIS");
	aff_signal_ids.clear();
	aff_signal_ids.insert(sig_id);
	aff_signal_ids.insert(bdate_id);
	aff_signal_ids.insert(gender_id);
	aff_signal_ids.insert(diag_id);
	req_signal_ids.clear();
	req_signal_ids.insert(sig_id);
}

int RepClearSignalByDiag::_apply(PidDynamicRec &rec, vector<int> &time_points, vector<vector<float>> &attributes_mat)
{
	if (time_points.size() != rec.get_n_versions())
	{
		MERR("nversions mismatch\n");
		return -1;
	}
	if (bdate_id < 0)
		MTHROW_AND_ERR("Error in RepFilterByChannel::_apply - bdate_id is not initialized - bad call\n");
	if (sig_id < 0)
		MTHROW_AND_ERR("Error in RepFilterByChannel::_apply - sig_id is not initialized - bad call\n");

	// first lets fetch "static" signals without Time field:
	allVersionsIterator vit(rec, sig_id);
	UniversalSigVec diag;

	for (int iver = vit.init(); !vit.done(); iver = vit.next())
	{
		vector<int> censor_times;
		rec.uget(diag_id, iver, diag);
		for (int i = 0; i < diag.len; ++i)
		{
			// Check if should add diag time to censor times
			bool passed = lut_censor[diag.Val(i)];
			if (passed)
				censor_times.push_back(diag.Time(i, 0));
		}

		float egfr;
		int final_size = 0;
		vector<float> v_vals;
		vector<int> v_times;
		UniversalSigVec sig;
		rec.uget(sig_id, iver, sig);
		int d0 = 0;
		int num_exclusion = 0;

		for (int i = 0; i < sig.len; ++i)
		{
			if (time_points.size() != 0 && sig.Time(i) > time_points[iver])
				break;
			// if test value is OK, keep it and break
			int byear = int(medial::repository::get_value(rec, bdate_id) / 10000);
			int gender = int(medial::repository::get_value(rec, gender_id));
			int age = int(sig.Time(i) / 10000) - byear;
			
			egfr = get_eGFR_CKD_EPI(age, sig.Val(i, 0), gender);
			
			if (egfr >= 60 || num_exclusion >= max_exclusion)
			{
				v_times.push_back(sig.Time(i, 0));
				v_vals.push_back(sig.Val(i, 0));
				++final_size;
				continue;
			}
			// else ...
			int delta_days = MED_MAT_MISSING_VALUE;
			for (int d = d0; d < censor_times.size(); ++d)
			{
				delta_days = med_time_converter.convert_times(MedTime::Date, MedTime::Days, sig.Time(i, 0)) -
							 med_time_converter.convert_times(MedTime::Date, MedTime::Days, censor_times[d]);

				if (delta_days <= time_window)
				{
					d0 = d;
					break;
				}
			}

			if (delta_days >= -time_window)
			{
				num_exclusion = num_exclusion + 1;
				// cout << "exclusion " << sig.Time(i, 0) << " " << sig.Val(i, 0) << endl;
				continue;
			};

			// keep the sample
			v_times.push_back(sig.Time(i, 0));
			v_vals.push_back(sig.Val(i, 0));
			++final_size;
		}

		// pushing virtual data into rec (into orig version)
		rec.set_version_universal_data(sig_id, iver, &v_times[0], &v_vals[0], final_size);
	}
	return 0;
}