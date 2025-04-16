#ifndef __REP_FILTER_BY_OTHER_SIGNAL
#define __REP_FILTER_BY_OTHER_SIGNAL

#include "RepProcess.h"

/**
 * Filter signal by diag
 * removing signals from dates were existing diag (other, not relevant to outcome) is the cause of signal value, so signal value is not indicative
 */
class RepClearSignalByDiag : public RepProcessor
{
private:
    vector<char> lut_censor;

public:
    string signal_name = "Creatinine"; ///< Hard coded ...
    int time_window;                   ///< range in days for diag to delete signal
    int max_exclusion;                 ///< max number of tests to drop
    vector<string> diag_list;          ///< list of diags that drop tests

    RepClearSignalByDiag() 
    {
        processor_type = REP_PROCESS_FILTER_BY_DIAG;
        max_exclusion = 9999;
    }

    /// @snippet RepClearSignalByDiag.cpp RepClearSignalByDiag::init
    int init(map<string, string> &mapper);

    void init_tables(MedDictionarySections &dict, MedSignals &sigs);
    void set_required_signal_ids(MedDictionarySections &dict) {};
    void set_affected_signal_ids(MedDictionarySections &dict) {};

    // Applying
    /// <summary> apply processing on a single PidDynamicRec at a set of time-points : Should be implemented for all inheriting classes </summary>
    int _apply(PidDynamicRec &rec, vector<int> &time_points, vector<vector<float>> &attributes_mat);

    ADD_CLASS_NAME(RepClearSignalByDiag)
    ADD_SERIALIZATION_FUNCS(processor_type, 
                            unconditional, req_signals, aff_signals, virtual_signals,
                            virtual_signals_generic, time_window, diag_list)
private:
    int v_out_sid = -1;
    int sig_id = -1;
    int bdate_id = -1;
    int gender_id = -1;
    int diag_id = -1;
    int v_out_n_vals, v_out_n_times;
};

MEDSERIALIZE_SUPPORT(RepClearSignalByDiag)

#endif
