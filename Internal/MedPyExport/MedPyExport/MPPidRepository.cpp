#include "MPPidRepository.h"
#include "MPSigExporter.h"
#include "MPDictionary.h"

#include <time.h>
#include <string>

#include "InfraMed/InfraMed/MedConvert.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/Utils.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedModel.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"

MPPidRepository::MPPidRepository() : o(new MedPidRepository()),dict(this) {};
MPPidRepository::~MPPidRepository() { delete o; o = nullptr; };
int MPPidRepository::dict_section_id(const std::string &secName) { return o->dict.section_id(secName); };
std::string MPPidRepository::dict_name(int section_id, int id) { return o->dict.name(section_id, id); };
int MPPidRepository::read_all(const std::string &conf_fname) { return val_or_exception(o->read_all(conf_fname),string("Error: read_all() failed")); };

int MPPidRepository::read_all_i(const std::string &conf_fname, const std::vector<int> &pids_to_take, const std::vector<int> &signals_to_take)
{
	return o->read_all(conf_fname, pids_to_take, signals_to_take);
};
/*
int MPPidRepository::read_all_s(const std::string &conf_fname, const std::vector<int> &pids_to_take, const vector<std::string> &signals_to_take)
{
	return o->read_all(conf_fname, pids_to_take, signals_to_take);
};*/

int MPPidRepository::read_all(const std::string &conf_fname, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take), const std::vector<std::string> &signals_to_take)
{
	vector<int> pids_tt;
	buf_to_vector(pids_to_take, num_pids_to_take, pids_tt);
	return val_or_exception(o->read_all(conf_fname, pids_tt, signals_to_take),string("Error: read_all() failed"));
}

int MPPidRepository::loadsig(const std::string& signame) { return o->load(signame); };
int MPPidRepository::init(const std::string &conf_fname) { return o->init(conf_fname); };
const std::vector<int>& MPPidRepository::MEDPY_GET_pids() { return o->index.pids; };
int MPPidRepository::sig_id(const std::string& signame) { return o->dict.id(signame); };
int MPPidRepository::sig_type(const std::string& signame) { return o->sigs.type(signame); };
MPSigVectorAdaptor MPPidRepository::uget(int pid, int sid) { MPSigVectorAdaptor ret; o->uget(pid, sid, *((UniversalSigVec*)(ret.o))); return ret; };

std::vector<bool> MPPidRepository::dict_prep_sets_lookup_table(int section_id, const std::vector<std::string> &set_names) {
	vector<char> lut_cvec;
	o->dict.prep_sets_lookup_table(section_id, set_names, lut_cvec);
	vector<bool> lut_bvec;
	lut_bvec.reserve(lut_cvec.size());
	for (int i = 0; i < lut_cvec.size(); i++)
		lut_bvec.push_back(lut_cvec[i] != 0);
	//vector<bool> lut_bvec(lut_cvec.begin(), lut_cvec.end());
	return lut_bvec;
}

MPSigExporter MPPidRepository::export_to_numpy(string signame, MEDPY_NP_INPUT(int* pids_to_take, unsigned long long num_pids_to_take), int use_all_pids, int translate_flag) {
	return MPSigExporter(*this, signame, pids_to_take, num_pids_to_take, use_all_pids, translate_flag);
}

// ****************************** MPSig      *********************************
MPSig::MPSig(void* _o, int index) : o(_o), idx(index) {  };
MPSig::MPSig(const MPSig& other) { o = other.o; idx = other.idx; };


int MPSig::time(int chan) { return ((UniversalSigVec*)o)->Time(idx, chan); }
float MPSig::val(int chan) { return ((UniversalSigVec*)o)->Val(idx, chan); }
int MPSig::timeU(int to_time_unit) { return ((UniversalSigVec*)o)->TimeU(idx, to_time_unit); }
int MPSig::date(int chan) { return ((UniversalSigVec*)o)->Date(idx, chan); }
int MPSig::years(int chan) { return ((UniversalSigVec*)o)->Years(idx, chan); }
int MPSig::months(int chan) { return ((UniversalSigVec*)o)->Months(idx, chan); }
int MPSig::days(int chan) { return ((UniversalSigVec*)o)->Days(idx, chan); }
int MPSig::hours(int chan) { return ((UniversalSigVec*)o)->Hours(idx, chan); }
int MPSig::minutes(int chan) { return ((UniversalSigVec*)o)->Minutes(idx, chan); }

// ****************************** MPSigVectorAdaptor      *********************

MPSigVectorAdaptor::MPSigVectorAdaptor() { o = new UniversalSigVec(); };
MPSigVectorAdaptor::MPSigVectorAdaptor(const MPSigVectorAdaptor& other) { o = new UniversalSigVec(); *((UniversalSigVec*)(this->o)) = *((UniversalSigVec*)(other.o)); };
MPSigVectorAdaptor::~MPSigVectorAdaptor() { delete ((UniversalSigVec*)o); };

int MPSigVectorAdaptor::__len__() { return ((UniversalSigVec*)o)->len; };
MPSig MPSigVectorAdaptor::__getitem__(int i) { return MPSig(o, i); };

int MPSigVectorAdaptor::MEDPY_GET_type() { return ((UniversalSigVec*)o)->get_type(); }

int MPSigVectorAdaptor::MEDPY_GET_n_time_channels() { return ((UniversalSigVec*)o)->n_time_channels(); }
int MPSigVectorAdaptor::MEDPY_GET_n_val_channels() { return ((UniversalSigVec*)o)->n_val_channels(); }
int MPSigVectorAdaptor::MEDPY_GET_time_unit() { return ((UniversalSigVec*)o)->time_unit(); }
int MPSigVectorAdaptor::MEDPY_GET_size() { return (int)(((UniversalSigVec*)o)->size()); }


