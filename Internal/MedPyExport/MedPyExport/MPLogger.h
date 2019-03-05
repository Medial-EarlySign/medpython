#ifndef __MED__MPLOGGER__H__
#define __MED__MPLOGGER__H__

#include "MedPyCommon.h"

class MPLogger {

public:
	static const int LOG_APP_;
	static const int LOG_DEF_;
	static const int LOG_INFRA_;
	static const int LOG_REP_;
	static const int LOG_INDEX_;
	static const int LOG_DICT_;
	static const int LOG_SIG_;
	static const int LOG_CONVERT_;
	static const int LOG_MED_UTILS_;
	static const int LOG_MEDMAT_;
	static const int LOG_MEDIO_;
	static const int LOG_DATA_STRUCTURES_;
	static const int LOG_MEDALGO_;
	static const int LOG_MEDFEAT_;
	static const int LOG_MEDSTAT_;
	static const int LOG_REPCLEANER_;
	static const int LOG_FTRGNRTR_;
	static const int LOG_CV_;
	static const int LOG_FEATCLEANER_;
	static const int LOG_VALCLNR_;
	static const int MED_SAMPLES_CV_;
	static const int LOG_FEAT_SELECTOR_;
	static const int LOG_SMPL_FILTER_;
	static const int LOG_SRL_;
	static const int LOG_MED_MODEL_;
	static const int LOG_REPTYPE_;
	static const int MAX_LOG_SECTION_;

	static const int NO_LOG_LEVEL_;
	static const int MIN_LOG_LEVEL_;
	static const int DEBUG_LOG_LEVEL_;
	static const int LOG_DEF_LEVEL_;
	static const int MAX_LOG_LEVEL_;
	static const int VERBOSE_LOG_LEVEL_;
	
	void init_all_levels(int level);
	int init_all_files(const string &fname);
	void init_level(int section, int level);
	int init_file(int section, const string &fname);
	int add_file(int section, const string &fname);

	void init_out(const string &fname="");

};


#endif //!__MED__MPLOGGER__H__
