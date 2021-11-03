#include "MPBootstrap.h"
#include "MedStat/MedStat/MedBootstrap.h"

/*************************************************/
MPStringBtResultMap::MPStringBtResultMap() { o = new std::map<std::string, std::map<std::string, float> >(); }
MPStringBtResultMap::MPStringBtResultMap(const MPStringBtResultMap& other)
{
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new map<string, std::map<std::string, float> >();
		*o = *other.o;
	}
};
MPStringBtResultMap::MPStringBtResultMap(std::map<std::string, std::map<std::string, float> >* ptr) { o_owned = true; o = ptr; };
MPStringBtResultMap::~MPStringBtResultMap() { if (o_owned) delete o; };
int MPStringBtResultMap::__len__() { return (int)o->size(); };
MPStringFloatMapAdaptor MPStringBtResultMap::__getitem__(const std::string& i) { return MPStringFloatMapAdaptor(&o->operator[](i)); };
void MPStringBtResultMap::__setitem__(const std::string &key, MPStringFloatMapAdaptor& val) {
	o->operator[](key) = *val.o;
}
std::vector<std::string> MPStringBtResultMap::keys() {
	std::vector<string> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
}

MPStringBtResultMap& MPStringBtResultMap::operator=(const MPStringBtResultMap& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, std::map<std::string, float> >();
		*o = *(other.o);
	}
	return *this;
}
/*************************************************/

MPBootstrap::MPBootstrap() { o = new MedBootstrap(); };
MPBootstrap::~MPBootstrap() { delete o; };

int MPBootstrap::MEDPY_GET_sample_per_pid() { return o->sample_per_pid; };
void MPBootstrap::MEDPY_SET_sample_per_pid(int _sample_per_pid) { o->sample_per_pid = _sample_per_pid; };
int MPBootstrap::MEDPY_GET_nbootstrap() { return o->loopCnt; };
void MPBootstrap::MEDPY_SET_nbootstrap(int _nbootstrap) { o->loopCnt = _nbootstrap; }
int MPBootstrap::MEDPY_GET_ROC_score_min_samples() { return o->roc_Params.score_min_samples; };
void MPBootstrap::MEDPY_SET_ROC_score_min_samples(int _ROC_score_min_samples) { o->roc_Params.score_min_samples = _ROC_score_min_samples; }
int MPBootstrap::MEDPY_GET_ROC_score_bins() { return o->roc_Params.score_bins; };
void MPBootstrap::MEDPY_SET_ROC_score_bins(int _ROC_score_bins) { o->roc_Params.score_bins = _ROC_score_bins; }

float MPBootstrap::MEDPY_GET_ROC_max_diff_working_point() { return o->roc_Params.max_diff_working_point; };
void MPBootstrap::MEDPY_SET_ROC_max_diff_working_point(float _ROC_max_diff_working_point) { o->roc_Params.max_diff_working_point = _ROC_max_diff_working_point; }
float MPBootstrap::MEDPY_GET_ROC_score_resolution() { return o->roc_Params.score_resolution; };
void MPBootstrap::MEDPY_SET_ROC_score_resolution(float _ROC_score_resolution) { o->roc_Params.score_resolution = _ROC_score_resolution; }
bool MPBootstrap::MEDPY_GET_ROC_use_score_working_points() { return o->roc_Params.use_score_working_points; };
void MPBootstrap::MEDPY_SET_ROC_use_score_working_points(bool _ROC_use_score_working_points) { o->roc_Params.use_score_working_points = _ROC_use_score_working_points; }

std::vector<float> MPBootstrap::MEDPY_GET_ROC_working_point_FPR() { return o->roc_Params.working_point_FPR; };
void MPBootstrap::MEDPY_SET_ROC_working_point_FPR(const std::vector<float> &_ROC_working_point_FPR) { o->roc_Params.working_point_FPR = _ROC_working_point_FPR; }
std::vector<float> MPBootstrap::MEDPY_GET_ROC_working_point_PR() { return o->roc_Params.working_point_PR; };
void MPBootstrap::MEDPY_SET_ROC_working_point_PR(const std::vector<float> &_ROC_working_point_PR) { o->roc_Params.working_point_PR = _ROC_working_point_PR; }
std::vector<float> MPBootstrap::MEDPY_GET_ROC_working_point_SENS() { return o->roc_Params.working_point_SENS; };
void MPBootstrap::MEDPY_SET_ROC_working_point_SENS(const std::vector<float> &_ROC_working_point_SENS) { o->roc_Params.working_point_SENS = _ROC_working_point_SENS; }
std::vector<float> MPBootstrap::MEDPY_GET_ROC_working_point_Score() { return o->roc_Params.working_point_Score; };
void MPBootstrap::MEDPY_SET_ROC_working_point_Score(const std::vector<float> &_ROC_working_point_Score) { o->roc_Params.working_point_Score = _ROC_working_point_Score; }


void MPBootstrap::parse_cohort_line(const string &line) { o->parse_cohort_line(line); }

void MPBootstrap::get_cohort_from_arg(const string &single_cohort) { o->get_cohort_from_arg(single_cohort); }

void MPBootstrap::parse_cohort_file(const string &cohorts_path) { o->parse_cohort_file(cohorts_path); }

void MPBootstrap::init_from_str(const string &str) { o->init_from_string(str); }

MPStringBtResultMap MPBootstrap::bootstrap(MPSamples *samples, const string &rep_path) {
	std::map<std::string, std::map<std::string, float> > *res = new std::map<std::string, std::map<std::string, float> >();
	*res = o->bootstrap(*samples->o, rep_path);

	return MPStringBtResultMap(res);
}

MPStringBtResultMap MPBootstrap::bootstrap(const std::vector<float> &preds, const std::vector<float> &labels) {
	if (preds.size() != labels.size()) {
		throw runtime_error(string("vectors are not in the same size: preds.size=") + 
			std::to_string(preds.size()) + string(" labels.size=") + std::to_string(labels.size()));
	}
	MedSamples samples;
	samples.idSamples.resize(preds.size());
	for (int i = 0; i < (int)preds.size(); ++i)
	{
		samples.idSamples[i].id = i;
		samples.idSamples[i].samples.resize(1);
		MedSample &s = samples.idSamples[i].samples[0];
		s.id = i;
		s.outcome = labels[i];
		s.prediction = { preds[i] };
		s.outcomeTime = 20000101;
		s.time = s.outcomeTime;
	}

	std::map<std::string, std::vector<float> > additional_data;
	std::map<std::string, std::map<std::string, float> > *res = new std::map<std::string, std::map<std::string, float> >();
	*res = o->bootstrap(samples, additional_data);

	return MPStringBtResultMap(res);
}

void MPStringBtResultMap::to_file(const std::string &file_path) {
	write_pivot_bootstrap_results(file_path, *o);
}