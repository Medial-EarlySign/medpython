#include "MPLogger.h"
#include "Logger/Logger/Logger.h"

const int MPLogger::LOG_APP_ = LOG_APP;
const int MPLogger::LOG_DEF_ = LOG_DEF;
const int MPLogger::LOG_INFRA_ = LOG_INFRA;
const int MPLogger::LOG_REP_ = LOG_REP;
const int MPLogger::LOG_INDEX_ = LOG_INDEX;
const int MPLogger::LOG_DICT_ = LOG_DICT;
const int MPLogger::LOG_SIG_ = LOG_SIG;
const int MPLogger::LOG_CONVERT_ = LOG_CONVERT;
const int MPLogger::LOG_MED_UTILS_ = LOG_MED_UTILS;
const int MPLogger::LOG_MEDMAT_ = LOG_MEDMAT;
const int MPLogger::LOG_MEDIO_ = LOG_MEDIO;
const int MPLogger::LOG_DATA_STRUCTURES_ = LOG_DATA_STRUCTURES;
const int MPLogger::LOG_MEDALGO_ = LOG_MEDALGO;
const int MPLogger::LOG_MEDFEAT_ = LOG_MEDFEAT;
const int MPLogger::LOG_MEDSTAT_ = LOG_MEDSTAT;
const int MPLogger::LOG_REPCLEANER_ = LOG_REPCLEANER;
const int MPLogger::LOG_FTRGNRTR_ = LOG_FTRGNRTR;
const int MPLogger::LOG_CV_ = LOG_CV;
const int MPLogger::LOG_FEATCLEANER_ = LOG_FEATCLEANER;
const int MPLogger::LOG_VALCLNR_ = LOG_VALCLNR;
const int MPLogger::MED_SAMPLES_CV_ = MED_SAMPLES_CV;
const int MPLogger::LOG_FEAT_SELECTOR_ = LOG_FEAT_SELECTOR;
const int MPLogger::LOG_SMPL_FILTER_ = LOG_SMPL_FILTER;
const int MPLogger::LOG_SRL_ = LOG_SRL;
const int MPLogger::LOG_MED_MODEL_ = LOG_MED_MODEL;
const int MPLogger::LOG_REPTYPE_ = LOG_REPTYPE;
const int MPLogger::MAX_LOG_SECTION_ = MAX_LOG_SECTION;

const int MPLogger::NO_LOG_LEVEL_ = NO_LOG_LEVEL;
const int MPLogger::MIN_LOG_LEVEL_ = MIN_LOG_LEVEL;
const int MPLogger::DEBUG_LOG_LEVEL_ = DEBUG_LOG_LEVEL;
const int MPLogger::LOG_DEF_LEVEL_ = LOG_DEF_LEVEL;
const int MPLogger::MAX_LOG_LEVEL_ = MAX_LOG_LEVEL;
const int MPLogger::VERBOSE_LOG_LEVEL_ = VERBOSE_LOG_LEVEL;

void MPLogger::init_all_levels(int level) { global_logger.init_all_levels(level); };
int MPLogger::init_all_files(const string &fname) { return global_logger.init_all_files(fname); };
void MPLogger::init_level(int section, int level) { global_logger.init_level(section, level); }
int MPLogger::init_file(int section, const string &fname) { return global_logger.init_file(section, fname); };
int MPLogger::add_file(int section, const string &fname) { return global_logger.add_file(section, fname); };
void MPLogger::init_out(const string &fname) { if (fname == "") global_logger.init_out(); else global_logger.init_out(fname);};
