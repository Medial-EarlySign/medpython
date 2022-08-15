#include "RepProcess.h"

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void RepReoderChannels::set_signal_ids(MedSignals& sigs) {
	sid = sigs.sid(signal_name);
	if (sid < 0)
		MTHROW_AND_ERR("Error RepReoderChannels - can't find signal %s\n", signal_name.c_str());
}

int RepReoderChannels::init(map<string, string>& mapper) {
	//Assume all columns are of type float and all columns exists - can duplicate column also

	for (auto &it : mapper) {
		//! [RepReoderChannels::init]
		if (it.first == "signal")
			signal_name = it.second;
		else if (it.first == "new_order") {
			vector<string> tokens;
			boost::split(tokens, it.second, boost::is_any_of(","));
			new_order.resize(tokens.size());
			for (size_t i = 0; i < tokens.size(); ++i)
				new_order[i] = med_stoi(tokens[i]);
		}
		else if (it.first == "rp_type") {}
		else
			MTHROW_AND_ERR("Error in RepReoderChannels::init - Unsupported param \"%s\"\n",
				it.first.c_str());
		//! [RepReoderChannels::init]
	}
	req_signals = { signal_name };
	aff_signals = { signal_name };
	if (new_order.empty())
		MTHROW_AND_ERR("Error RepReoderChannels::init - please pass new_order param\n");

	return 0;
}

int RepReoderChannels::_apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat) {
	if (sid == -1) {
		MERR("RepReoderChannels::_apply - Uninitialized signalId(%s)\n", signal_name.c_str());
		return -1;
	}

	// Check that we have the correct number of dynamic-versions : one per time-point (if given)
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("RepReoderChannels::_apply - nversions mismatch\n");
		return -1;
	}
	//check that it can do it:
	GenericSigVec_mem usv_single;
	usv_single.manage = true; //clear memory in the end of this template object
	usv_single.init(rec.my_base_rep->sigs.Sid2Info[sid]);
	int struct_sz = usv_single.struct_size;

	if (new_order.size() != usv_single.n_val_channels()) {
		MERR("Error RepReoderChannels::_apply - missing columns for %s signal\n",
			signal_name.c_str());
		return -1;
	}
	for (size_t j = 0; j < new_order.size(); ++j)
		if (new_order[j] >= usv_single.n_val_channels()) {
			MERR("Error RepReoderChannels::_apply - column %d is now allowed for %s signal\n",
				new_order[j], signal_name.c_str());
			return -1;
		}
	//check all are same type
	char first_tp = 0;
	for (int j = 0; j < usv_single.n_val_channels(); ++j) {
		if (j == 0) {
			first_tp = usv_single.val_channel_types[j];
			continue;
		}
		if (usv_single.val_channel_types[j] != first_tp) {
			MERR("Error RepReoderChannels::_apply - column %d is not same type for %s signal\n",
				j, signal_name.c_str());
			return -1;
		}
	}
	vector<float> new_cols(usv_single.n_val_channels());
	vector<int> time_vec(usv_single.n_time_channels());
	usv_single.data = new char[struct_sz];
	usv_single.len = 1;

	differentVersionsIterator vit(rec, sid);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {

		// Do it 
		rec.uget(sid, iver, rec.usv); // get into the internal usv obeject - this statistically saves init time
		
		int len = rec.usv.len;
		
		vector<pair<int, vector<float>>> change(len);
		
		for (int i = 0; i < len; ++i) {
			// Construct new item with reordered fields
			for (size_t j = 0; j < new_order.size(); ++j) {
				char *src_ptr = ((char*)rec.usv.data) + i * struct_sz + rec.usv.val_channel_offsets[new_order[j]];
				char *dest_ptr = ((char *)usv_single.data + usv_single.val_channel_offsets[j]);
				memcpy(dest_ptr, src_ptr, GenericSigVec::type_enc::bytes_len(usv_single.val_channel_types[j]));
				//usv_single.setVal(0, j, rec.usv.Val(i, new_order[j])); //can't do - since it's not always float
			}
			for (int j = 0; j < time_vec.size(); ++j)
				usv_single.setTime(0, j, rec.usv.Time(i, j));

			//Change this item to the contructed object
			if (rec.change(sid, iver, i, usv_single.data) < 0) {
				MERR("Error in reordering channels for %s\n", signal_name.c_str());
				return -1;
			}
		}

	}

	return 0;
}

void RepReoderChannels::print() {
	MLOG("RepReoderChannels %s\n", signal_name.c_str());
}