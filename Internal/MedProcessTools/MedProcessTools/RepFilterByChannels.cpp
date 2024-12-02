#include "RepFilterByChannels.h"
#include <MedUtils/MedUtils/MedUtils.h>
#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void RepFilterByChannel::register_virtual_section_name_id(MedDictionarySections& dict) {
	dict.connect_to_section(output_name, dict.section_id(signal_name));
}

int RepFilterByChannel::init(map<string, string>& mapper) {
	string prefix_str = "filter_set_by_val_channel_";
	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		//! [RepFilterByChannel::init]
		if (it->first == "signal")
			signal_name = it->second;
		else if (it->first == "output_name")
			output_name = it->second;
		else if (it->first == "signal_type")
			signal_type = it->second;
		else if (boost::starts_with(it->first, prefix_str)) {
			int val_channel_f = med_stoi(it->first.substr(prefix_str.length()));
			if (filter_set_by_val_channel.size() <= val_channel_f)
				filter_set_by_val_channel.resize(val_channel_f + 1);
			vector<string> &f_v_sets = filter_set_by_val_channel[val_channel_f];
			boost::split(f_v_sets, it->second, boost::is_any_of(","));
		}
		else if (it->first == "rp_type") {}
		else
			MTHROW_AND_ERR("Unknown parameter \'%s\' for RepFilterByChannel\n", it->first.c_str())
			//! [RepFilterByChannel::init]
	}

	if (signal_name.empty())
		MTHROW_AND_ERR("Erorr in RepFilterByChannel - Must provide signal\n");
	
	req_signals.clear();
	req_signals.insert(signal_name);
	aff_signals.clear();
	aff_signals.insert(output_name);

	virtual_signals_generic.clear();
	string const_str = signal_type;
	virtual_signals_generic.push_back(pair<string, string>(output_name, const_str));

	return 0;
}

void RepFilterByChannel::init_tables(MedDictionarySections& dict, MedSignals& sigs) {
	filter_luts.resize(filter_set_by_val_channel.size());
	filter_chan.clear();
	filter_chan.reserve(filter_set_by_val_channel.size());
	int section_id = dict.section_id(signal_name);
	for (int i = 0; i < filter_set_by_val_channel.size(); ++i)
	{
		const vector<string> &set_vals = filter_set_by_val_channel[i];
		vector<char> &lut = filter_luts[i];
		
		if (set_vals.empty())
			continue;

		dict.prep_sets_lookup_table(section_id, set_vals, lut);
		filter_chan.push_back(i);
	}

	v_out_sid = sigs.sid(output_name);
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepFilterByChannel::init_tables - virtual output signal %s not found\n",
			output_name.c_str());
	sig_id = sigs.sid(signal_name);
	aff_signal_ids.clear();
	aff_signal_ids.insert(v_out_sid);
	req_signal_ids.clear();
	req_signal_ids.insert(sig_id);

	//don't store more channels than out signals defines
	const SignalInfo &out_info = sigs.Sid2Info.at(v_out_sid);
	v_out_n_times = out_info.n_time_channels;
	v_out_n_vals = out_info.n_val_channels;
	const SignalInfo &in_info = sigs.Sid2Info.at(sig_id);
	if (in_info.n_time_channels < v_out_n_times)
		MTHROW_AND_ERR("Error RepFilterByChannel::init_tables - output signal can't have more time channels then input\n");
	if (in_info.n_val_channels < v_out_n_vals)
		MTHROW_AND_ERR("Error RepFilterByChannel::init_tables - output signal can't have more val channels then input\n");
}

void RepFilterByChannel::fit_for_repository(MedPidRepository& rep) {
	bool is_virtual = false;
	if (rep.sigs.sid(output_name) > 0) {
		const SignalInfo &si = rep.sigs.Sid2Info[rep.sigs.sid(output_name)];
		if (!si.virtual_sig)
			virtual_signals_generic.clear(); //not virtual signal
		else
			is_virtual = true;
	}
	else
		is_virtual = true;

	if (is_virtual && virtual_signals_generic.empty())
		virtual_signals_generic.push_back(pair<string, string>(output_name, signal_type));
}

void RepFilterByChannel::print() {
	MLOG("RepFilterByChannel:: output_name: %s : signal %s : req_signals %s aff_signals %s\n",
		output_name.c_str(), signal_name.c_str(), medial::io::get_list(req_signals).c_str(), 
		medial::io::get_list(aff_signals).c_str());
}

int RepFilterByChannel::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {
	if (time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}
	if (v_out_sid < 0)
		MTHROW_AND_ERR("Error in RepFilterByChannel::_apply - v_out_sid is not initialized - bad call\n");
	if (sig_id < 0)
		MTHROW_AND_ERR("Error in RepFilterByChannel::_apply - sig_id is not initialized - bad call\n");
	
	//first lets fetch "static" signals without Time field:
	allVersionsIterator vit(rec, sig_id);
	UniversalSigVec usv;
	
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		rec.uget(sig_id, iver, usv);

		int final_size = 0;
		vector<float> v_vals;
		vector<int> v_times;
		for (int i = 0; i < usv.len; ++i)
		{
			//Check filters:
			bool passed = true;
			for (int val_ch : filter_chan)
			{
				int categ_val = usv.Val<int>(i, val_ch);
				passed = filter_luts[val_ch][categ_val];
				if (!passed)
					break;
			}
			if (passed) {
				//Add val channels:
				for (int j = 0; j < usv.n_time_channels() && j < v_out_n_times; ++j)
					v_times.push_back(usv.Time(i, j));
				//Add vals:
				for (int j = 0; j < usv.n_val_channels() && j < v_out_n_vals; ++j)
					v_vals.push_back(usv.Val(i, j));
				++final_size;
			}
		}
		
		// pushing virtual data into rec (into orig version)
		rec.set_version_universal_data(v_out_sid, iver, &v_times[0], &v_vals[0], final_size);
	}
	return 0;
}