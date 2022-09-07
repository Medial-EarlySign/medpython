#ifndef __REP_FILTER_BY_CHAN
#define __REP_FILTER_BY_CHAN

#include "RepProcess.h"

/**
* Filter signal by different set of values in channels
*/
class RepFilterByChannel : public RepProcessor {
private:
	vector<vector<char>> filter_luts;
	vector<int> filter_chan; ///< list of channel with filter
public:
	string output_name; ///< names of signal created by the processor or same signal name
	string signal_name; ///< names of input signal used by the processor
	string signal_type; ///< the signal type definition to create

	vector<vector<string>> filter_set_by_val_channel; ///< filter set by value channels. can be initialized by "filter_set_by_val_channel_X":"string_set_for_val_channel_X",

	RepFilterByChannel() {
		processor_type = REP_PROCESS_FILTER_BY_CHANNELS;
		output_name = "";
		signal_type = "T(i),V(f)";
		signal_name = "";
	}

	void register_virtual_section_name_id(MedDictionarySections& dict);

	/// @snippet RepFilterByChannel.cpp RepFilterByChannel::init
	int init(map<string, string>& mapper);

	void init_tables(MedDictionarySections& dict, MedSignals& sigs);
	void set_required_signal_ids(MedDictionarySections& dict) {};
	void set_affected_signal_ids(MedDictionarySections& dict) {};

	void fit_for_repository(MedPidRepository& rep);

	// Applying
	/// <summary> apply processing on a single PidDynamicRec at a set of time-points : Should be implemented for all inheriting classes </summary>
	int _apply(PidDynamicRec& rec, vector<int>& time_points, vector<vector<float>>& attributes_mat);

	void print();
	ADD_CLASS_NAME(RepFilterByChannel)
		ADD_SERIALIZATION_FUNCS(processor_type, output_name, signal_name,
			unconditional, req_signals, aff_signals, virtual_signals,
			virtual_signals_generic, signal_type, filter_set_by_val_channel)
private:
	int v_out_sid = -1;
	int sig_id = -1;
	int v_out_n_vals, v_out_n_times;
};

MEDSERIALIZE_SUPPORT(RepFilterByChannel)


#endif
