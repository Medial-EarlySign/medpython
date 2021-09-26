#include "RepCategoryDescenders.h"

#define LOCAL_SECTION LOG_REPCLEANER
#define LOCAL_LEVEL	LOG_DEF_LEVEL

void RepCategoryDescenders::init_lists() {
	req_signals.insert(signalName);
	aff_signals.insert(signalName);
}

void RepCategoryDescenders::set_signal_ids(MedSignals& sigs) {
	signalId = sigs.sid(signalName);
	if (signalId < 0)
		MTHROW_AND_ERR("Error RepCategoryDescenders - can't find siganl %s\n", signalName.c_str());
}

void RepCategoryDescenders::init_tables(MedDictionarySections& dict, MedSignals& sigs) {
	int section_id = dict.section_id(signalName);
	for (const auto &it : dict.dict(section_id)->Set2Members)
		_set2Members[it.first] = it.second;

	if (val_channel >= sigs.Sid2Info[signalId].n_val_channels)
		MTHROW_AND_ERR("Error RepCategoryDescenders::init_tables - signal %s has only %d value channels, requested %d\n",
			signalName.c_str(), sigs.Sid2Info[signalId].n_val_channels, val_channel);
}

int RepCategoryDescenders::init(map<string, string>& mapper) {
	for (auto &it : mapper)
	{
		//! [RepCategoryDescenders::init]
		if (it.first == "signal")
			signalName = it.second;
		else if (it.first == "val_channel")
			val_channel = med_stoi(it.second);
		else if (it.first == "rp_type") {}
		else
			MTHROW_AND_ERR("Error in RepCategoryDescenders::init - Unsupported param \"%s\"\n",
				it.first.c_str());
		//! [RepCategoryDescenders::init]
	}

	if (signalName.empty())
		MTHROW_AND_ERR("Error RepCategoryDescenders::init - Must provide signal\n");
	init_lists();

	return 0;
}

int RepCategoryDescenders::_apply(PidDynamicRec& rec, vector<int>& time_points,
	vector<vector<float>>& attributes_mat) {
	// Sanity check
	if (signalId == -1) {
		MERR("RepCategoryDescenders::_apply - Uninitialized signalId(%s)\n", signalName.c_str());
		return -1;
	}

	// Check that we have the correct number of dynamic-versions : one per time-point (if given)
	if (time_points.size() != 0 && time_points.size() != rec.get_n_versions()) {
		MERR("nversions mismatch\n");
		return -1;
	}

	allVersionsIterator vit(rec, signalId);
	for (int iver = vit.init(); !vit.done(); iver = vit.next()) {
		int final_size = 0;
		vector<float> values;
		vector<int> times;
		UniversalSigVec usv;
		rec.uget(signalId, iver, usv);
		for (unsigned int i = 0; i < usv.len; ++i) {
			int val = usv.Val(i, val_channel);
			unordered_set<int> all_vals;
			all_vals.insert(val);
			if (_set2Members.find(val) != _set2Members.end())
				all_vals.insert(_set2Members.at(val).begin(), _set2Members.at(val).end());
			//ADD DATA for each all_vals:
			for (int v : all_vals)
			{
				//Create Time channels:
				for (int t = 0; t < usv.n_time_channels(); ++t)
					times.push_back(usv.Time(i, t));
				//Create Value channels:
				for (int v_idx = 0; v_idx < usv.n_val_channels(); ++v_idx) {
					if (v_idx != val_channel)
						values.push_back(usv.Val(i, v_idx));
					else
						values.push_back(v);
				}
				//Mark as new record
				++final_size;
			}
		}

		rec.set_version_universal_data(signalId, iver, &times[0], &values[0], final_size);
	}

	return 0;
}