#include "MPSamples.h"
#include "MPFeatures.h"
#include <algorithm>

#include "InfraMed/InfraMed/MedConvert.h"
#include "InfraMed/InfraMed/InfraMed.h"
#include "InfraMed/InfraMed/Utils.h"
#include "MedUtils/MedUtils/MedUtils.h"
#include "InfraMed/InfraMed/MedPidRepository.h"
#include "MedProcessTools/MedProcessTools/MedModel.h"
#include "MedProcessTools/MedProcessTools/SampleFilter.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"


/******************** MPSample ****************/

MPSample::MPSample() { o = new MedSample(); };
MPSample::~MPSample() { if (o_owned) delete o; };
MPSample::MPSample(const MPSample& other) { o_owned = other.o_owned; o = other.o; };
MPSample::MPSample(MedSample* from_ptr) { o_owned = false; o = from_ptr; };

int MPSample::MEDPY_GET_id() { return o->id; };
void MPSample::MEDPY_SET_id(int new_id) { o->id = new_id; };

int MPSample::MEDPY_GET_split() { return o->split; };
void MPSample::MEDPY_SET_split(int new_sp) { o->split = new_sp; };

int MPSample::MEDPY_GET_time() { return o->time; };
void MPSample::MEDPY_SET_time(int new_time) { o->time = new_time; };

int MPSample::MEDPY_GET_outcome() { return o->outcome; };
void MPSample::MEDPY_SET_outcome(int new_outcome) { o->outcome = new_outcome; };

MPIdSamplesVectorAdaptor MPSamples::MEDPY_GET_idSamples() { return MPIdSamplesVectorAdaptor(&(o->idSamples)); };

int MPSample::MEDPY_GET_outcomeTime() { return o->outcomeTime; };
void MPSample::MEDPY_SET_outcomeTime(int new_outcome_time) { o->outcomeTime = new_outcome_time; };

void MPSample::MEDPY_GET_prediction(MEDPY_NP_OUTPUT(float** out_predbuf, unsigned long long* out_predbuf_len)) {
	*out_predbuf_len = o->prediction.size();
	int bufsize = sizeof(float)*(*out_predbuf_len);
	*out_predbuf = (float*)malloc(bufsize);
	memcpy((void*)(*out_predbuf), (void*)(o->prediction.data()), bufsize);
}

void MPSample::MEDPY_SET_prediction(MEDPY_NP_INPUT(float* in_predbuf, unsigned long long in_predbuf_len)) {
	if (in_predbuf_len <= 0) {
		o->prediction.clear();
		return;
	}
	o->prediction.resize(in_predbuf_len);
	memcpy((void*)(o->prediction.data()), (void*)(in_predbuf), sizeof(float)*in_predbuf_len);
}
void MPSample::print_(const string prefix) { o->print(prefix); };
void MPSample::print_() { o->print(); };

int MPSample::parse_from_string(string &s, int time_unit) { return o->parse_from_string(s, time_unit); };
void MPSample::write_to_string(string &s, int time_unit) { return o->write_to_string(s, time_unit); };

MPSample MPSample::__copy__() {
	MPSample ret;
	*(ret.o) = *(o);
	return ret;
}

MPSerializableObject MPSample::asSerializable() { return MPSerializableObject(o); }

void MPSamples::set_samples(MEDPY_NP_INPUT(int * patients, unsigned long long patient_size),
	int _time) {
	if (_time < 0) {
		time_t theTime = time(NULL);
		struct tm *now;
#if defined(__unix__) || defined(__APPLE__)
		now = localtime(&theTime);
#else
		struct tm now_m;
		now = &now_m;
		localtime_s(now, &theTime);
#endif

		_time = (now->tm_year + 1900) * 10000 + (now->tm_mon + 1) * 100 + now->tm_mday; //set to current time
	}

	o->clear();
	for (size_t i = 0; i < patient_size; ++i)
	{
		MedIdSamples smp_id(patients[i]);
		MedSample s;
		s.id = patients[i];
		s.time = _time;
		smp_id.samples.push_back(s);

		o->idSamples.push_back(smp_id);
	}
	o->sort_by_id_date();
}

/************ MPIdSamples ************/


int MPIdSamples::MEDPY_GET_id() { return o->id; };
void MPIdSamples::MEDPY_SET_id(int new_id) { o->id = new_id; };

int MPIdSamples::MEDPY_GET_split() { return o->split; };
void MPIdSamples::MEDPY_SET_split(int new_sp) { o->split = new_sp; };
MPSampleVectorAdaptor MPIdSamples::MEDPY_GET_samples() { return MPSampleVectorAdaptor(&(o->samples)); }

MPIdSamples::MPIdSamples(int _id) { o = new MedIdSamples(_id); };
MPIdSamples::MPIdSamples() { o = new MedIdSamples(); };
MPIdSamples::MPIdSamples(MedIdSamples* ptr) { o_owned = false; o = ptr; };
MPIdSamples::MPIdSamples(const MPIdSamples& other) { o_owned = false; o = other.o; };
MPIdSamples::~MPIdSamples() { if (o_owned) delete o; };

void MPIdSamples::set_split(int _split) { o->set_split(_split); };
bool MPIdSamples::same_as(MPIdSamples &other, int mode) { return o->same_as(*(other.o), mode); };
MPSerializableObject MPIdSamples::asSerializable() { return MPSerializableObject(o); }

/************ MPSampleVectorAdaptor ************/


MPSampleVectorAdaptor::MPSampleVectorAdaptor() { o = new vector<MedSample>(); };
MPSampleVectorAdaptor::MPSampleVectorAdaptor(vector<MedSample>* ptr) { o_owned = false; o = ptr; };
MPSampleVectorAdaptor::MPSampleVectorAdaptor(const MPSampleVectorAdaptor& other) {
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new vector<MedSample>();
		*o = *other.o;
	}
};
MPSampleVectorAdaptor::~MPSampleVectorAdaptor() { if (o_owned) delete o; };
int MPSampleVectorAdaptor::__len__() { return (int)o->size(); };
MPSample MPSampleVectorAdaptor::__getitem__(int i) { return MPSample(&(o->at(i))); };
void MPSampleVectorAdaptor::__setitem__(int i, MPSample& val) { o->at(i) = *(val.o); };
void MPSampleVectorAdaptor::append(MPSample& val) { o->push_back(*(val.o)); };

void MPSampleVectorAdaptor::override_splits(int nfolds) {
	map<int, int> id_folds;
	int idx = 0;
	for (auto& sample : *o) {
		int id = sample.id;
		if (id_folds.find(id) == id_folds.end()) {
			id_folds[id] = idx % nfolds;
			idx++;
		}
		sample.split = id_folds[id];
	}
}

int MPSampleVectorAdaptor::nSplits() {
	return medial::process::nSplits(*o);
}

void MPSampleVectorAdaptor::append_vec(MPSampleVectorAdaptor& other) {
	o->insert(o->end(), other.o->begin(), other.o->end());
}


/************ MPSamples ************/

MPSamples::MPSamples() { o = new MedSamples(); };
MPSamples::MPSamples(MedSamples* ptr) { o_owned = false; o = ptr; };
MPSamples::MPSamples(const MPSamples& other) {
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new MedSamples();
		*o = *other.o;
	}
};
MPSamples::~MPSamples() { if (o_owned) delete o; o = nullptr; };

void MPSamples::append(MPSamples& newSamples) { o->idSamples.insert(o->idSamples.end(), newSamples.o->idSamples.begin(), newSamples.o->idSamples.end()); }
void MPSamples::subtract(MPSamples& _dont_use) { o->subtract(*(_dont_use.o)); }
void MPSamples::split_train_test(MPSamples& train, MPSamples& test, float p_test) { o->split_train_test((*(train.o)), (*(test.o)), p_test); }

int MPSamples::read_from_bin_file(const string& file_name) { return val_or_exception(o->read_from_bin_file(file_name), "Cannot read Samples from bin file " + file_name); }
int MPSamples::write_to_bin_file(const string& file_name) { return val_or_exception(o->write_to_bin_file(file_name), "Cannot write Samples to bin file " + file_name); }

int MPSamples::read_from_file(const string& file_name) { return val_or_exception(o->read_from_file(file_name), "Cannot read Samples from file " + file_name); };
int MPSamples::write_to_file(const string& file_name) { return val_or_exception(o->write_to_file(file_name), "Cannot write Samples to file " + file_name); };

void MPSamples::get_preds(MEDPY_NP_OUTPUT(float** preds_buf, unsigned long long* preds_buf_len)) {
	vector<float> ret;
	o->get_preds(ret);
	vector_to_buf(ret, preds_buf, preds_buf_len);
}
void MPSamples::get_y(MEDPY_NP_OUTPUT(float** y_buf, unsigned long long* y_buf_len)) {
	vector<float> ret;
	o->get_y(ret);
	vector_to_buf(ret, y_buf, y_buf_len);
}
void MPSamples::get_categs(MEDPY_NP_OUTPUT(float** categs_buf, unsigned long long* categs_buf_len)) {
	vector<float> ret;
	o->get_categs(ret);
	vector_to_buf(ret, categs_buf, categs_buf_len);
}

MPSampleVectorAdaptor MPSamples::export_to_sample_vec() {
	MPSampleVectorAdaptor ret;
	o->export_to_sample_vec(*(ret.o));
	return ret;
}
void MPSamples::import_from_sample_vec(MPSampleVectorAdaptor& vec_samples, bool allow_split_inconsistency) {
	o->import_from_sample_vec(*(vec_samples.o), allow_split_inconsistency);
}

MPPandasAdaptor MPSamples::MEDPY__from_df_adaptor() {
	MPPandasAdaptor ret;
	ret.set_type_requirement("id", "int");
	ret.set_type_requirement("outcome", "float");
	ret.set_type_requirement("outcomeTime", "int");
	ret.set_type_requirement("split", "int");
	ret.set_type_requirement("time", "int");
	ret.set_type_requirement("pred_\\d+", "float");
	ret.set_type_requirement("attr_\\S+", "float");
	return ret;
}

void MPSamples::MEDPY__from_df(MPPandasAdaptor& pandas_df) {
	static const string pred_ = "pred_";
	static const string attr_ = "attr_";
	static const string str_attr_ = "str_attr_";
	vector<MedSample> vms;

	for (string col_name : pandas_df.keys()) {
		if (col_name == "id") {
			vector<int> vec;
			pandas_df.pull_col_as_vector(col_name, vec);
			vms.resize(max(vms.size(), vec.size()));
			for (int i = 0; i < vec.size(); ++i)
				vms[i].id = vec[i];
		}
		else if (col_name == "outcome") {
			vector<float> vec;
			pandas_df.pull_col_as_vector(col_name, vec);
			vms.resize(max(vms.size(), vec.size()));
			for (int i = 0; i < vec.size(); ++i)
				vms[i].outcome = vec[i];
		}
		else if (col_name == "outcomeTime") {
			vector<int> vec;
			pandas_df.pull_col_as_vector(col_name, vec);
			vms.resize(max(vms.size(), vec.size()));
			for (int i = 0; i < vec.size(); ++i)
				vms[i].outcomeTime = vec[i];
		}
		else if (col_name == "split") {
			vector<int> vec;
			pandas_df.pull_col_as_vector(col_name, vec);
			vms.resize(max(vms.size(), vec.size()));
			for (int i = 0; i < vec.size(); ++i)
				vms[i].split = vec[i];
		}
		else if (col_name == "time") {
			vector<int> vec;
			pandas_df.pull_col_as_vector(col_name, vec);
			vms.resize(max(vms.size(), vec.size()));
			for (int i = 0; i < vec.size(); ++i)
				vms[i].time = vec[i];
		}
		//import pred_[N] fields
		else if (col_name.compare(0, pred_.length(), pred_) == 0) {
			vector<float> vec;
			pandas_df.pull_col_as_vector(col_name, vec);
			vms.resize(max(vms.size(), vec.size()));
			int pred_i = stoi(col_name.substr(pred_.length()));
			int max_size = max((int)vms[0].prediction.size(), pred_i + 1);
			if (max_size > vms[0].prediction.size())
				for (int i = 0; i < vec.size(); ++i) {
					vms[i].prediction.resize(max_size);
					vms[i].prediction[pred_i] = vec[i];
				}
			else for (int i = 0; i < vec.size(); ++i)
				vms[i].prediction[pred_i] = vec[i];
		}
		//import attr_[NAME] fields
		else if (col_name.compare(0, attr_.length(), attr_) == 0) {
			vector<float> vec;
			pandas_df.pull_col_as_vector(col_name, vec);
			vms.resize(max(vms.size(), vec.size()));
			string attr_name = col_name.substr(attr_.length());
			for (int i = 0; i < vec.size(); ++i)
				vms[i].attributes[attr_name] = vec[i];
		}
	}
	o->import_from_sample_vec(vms);
}

MPPandasAdaptor MPSamples::MEDPY__to_df() {

	MPPandasAdaptor ret;
	int record_count = o->nSamples();

	int max_predictions = 0;

	unordered_set<string> attribute_names_set;
	unordered_set<string> str_attribute_names_set;

	int* id_vec = (int*)malloc(sizeof(int)*record_count);
	int* split_vec = (int*)malloc(sizeof(int)*record_count);
	int* time_vec = (int*)malloc(sizeof(int)*record_count);
	float* outcome_vec = (float*)malloc(sizeof(float)*record_count);
	int* outcome_time_vec = (int*)malloc(sizeof(int)*record_count);

	for (auto &s : o->idSamples) {
		for (auto &samp : s.samples) {
			if (samp.prediction.size() > max_predictions)
				max_predictions = (int)samp.prediction.size();
			if (samp.attributes.size() != 0)
				for (auto& entry : samp.attributes)
					attribute_names_set.insert(entry.first);
			if (samp.str_attributes.size() != 0)
				for (auto& entry : samp.str_attributes)
					str_attribute_names_set.insert(entry.first);
		}
	}

	vector<string> attribute_names(attribute_names_set.begin(), attribute_names_set.end());
	vector<string> str_attribute_names(str_attribute_names_set.begin(), str_attribute_names_set.end());

	vector<float*> pred_vecs;
	for (int i = 0; i < max_predictions; i++) {
		pred_vecs.push_back((float*)malloc(sizeof(float)*record_count));
		//data_keys.push_back((string)"pred_" + to_string(i));
	}
	map<string, float*> attr_vecs;
	for (const auto& s : attribute_names) {
		attr_vecs[s] = (float*)malloc(sizeof(float)*record_count);
		//data_keys.push_back((string)"attr_" + attribute_names[i]);
	}

	int row_i = 0;
	for (auto &s : o->idSamples) {
		for (auto &samp : s.samples) {
			id_vec[row_i] = samp.id;
			split_vec[row_i] = samp.split;
			time_vec[row_i] = samp.time;
			outcome_vec[row_i] = samp.outcome;
			outcome_time_vec[row_i] = samp.outcomeTime;
			if (max_predictions != 0)
				for (int i = 0; i < samp.prediction.size(); i++)
					pred_vecs[i][row_i] = samp.prediction[i];
			if (samp.attributes.size() != 0)
				for (const auto& s : samp.attributes)
					attr_vecs[s.first][row_i] = s.second;

			row_i++;
		}
	}

	ret.push_column("id", id_vec, record_count, "int", false);
	ret.push_column("split", split_vec, record_count, "int", false);
	ret.push_column("time", time_vec, record_count, "int", false);
	ret.push_column("outcome", outcome_vec, record_count, "float", false);
	ret.push_column("outcomeTime", outcome_time_vec, record_count, "int", false);

	for (int i = 0; i < max_predictions; i++) {
		ret.push_column(string("pred_") + to_string(i), pred_vecs[i], record_count, "float", false);
	}

	for (const auto& s : attr_vecs) {
		ret.push_column(string("attr_") + s.first, s.second, record_count, "float", false);
	}

	return ret;
}

int MPSamples::get_predictions_size() {
	int ret1, ret2;
	ret1 = o->get_predictions_size(ret2);
	if (ret1 == -1)
		return ret1;
	return ret2;
}

//int get_all_attributes(vector<string>& attributes, vector<string>& str_attributes);
std::vector<string> MPSamples::get_attributes()
{
	std::vector<string> attr, str_attr;
	o->get_all_attributes(attr, str_attr);
	return attr;
};
std::vector<string> MPSamples::get_str_attributes() {
	std::vector<string> attr, str_attr;
	o->get_all_attributes(attr, str_attr);
	return str_attr;
};


void MPSamples::dilute(float prob) { return o->dilute((float)prob); };
int MPSamples::MEDPY_GET_time_unit() { return o->time_unit; };
void MPSamples::MEDPY_SET_time_unit(int new_time_unit) { o->time_unit = new_time_unit; };

void MPSamples::get_ids(MEDPY_NP_OUTPUT(int** ids, unsigned long long* num_ids)) {
	vector<int> ids_vec;
	o->get_ids(ids_vec);
	vector_to_buf(ids_vec, ids, num_ids);
}

void MPSamples::clear() { o->clear(); };
int MPSamples::insert_preds(MPFeatures& featuresData) { return o->insert_preds(*(featuresData.o)); };

void MPSamples::sort_by_id_date() { o->sort_by_id_date(); };
void MPSamples::normalize() { o->normalize(); };
bool MPSamples::same_as(MPSamples &other, int mode) { return o->same_as(*(other.o), mode); };
int MPSamples::nSamples() { return o->nSamples(); };
int MPSamples::nSplits() { return o->nSplits(); };

void MPSamples::insertRec(int pid, int time, float outcome, int outcomeTime) { o->insertRec(pid, time, outcome, outcomeTime); };
void MPSamples::insertRec(int pid, int time, float outcome, int outcomeTime, float pred) { o->insertRec(pid, time, outcome, outcomeTime, pred); };
void MPSamples::insertRec(int pid, int time) { o->insertRec(pid, time); };
int MPSamples::version() { return o->version(); };

MPSerializableObject MPSamples::asSerializable() { return MPSerializableObject(o); }

//void MPSamples::get_ids_v(int* out_pidvec_1, int out_pidvec_n_1) {  vector<int> ids; o->get_ids(ids); memcpy(out_pidvec_1, &ids[0], out_pidvec_n_1); };
//int MPSamples::get_ids_n() { return (int)o->idSamples.size(); };

void MPSamples::override_splits(int nfolds) {
	for (int idx = 0; idx < o->idSamples.size(); idx++) {
		o->idSamples[idx].split = idx % nfolds;
		for (auto& sample : o->idSamples[idx].samples)
			sample.split = idx % nfolds;
	}
}

/************ MPIdSamplesVectorAdaptor ************/

MPIdSamplesVectorAdaptor::MPIdSamplesVectorAdaptor() { o = new vector<MedIdSamples>(); };
MPIdSamplesVectorAdaptor::MPIdSamplesVectorAdaptor(vector<MedIdSamples>* ptr) { o_owned = false; o = ptr; };
MPIdSamplesVectorAdaptor::MPIdSamplesVectorAdaptor(const MPIdSamplesVectorAdaptor& other)
{
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new vector<MedIdSamples>();
		*o = *other.o;
	}
};
MPIdSamplesVectorAdaptor::~MPIdSamplesVectorAdaptor() { if (o_owned) delete o; };
int MPIdSamplesVectorAdaptor::__len__() { return (int)o->size(); };
MPIdSamples MPIdSamplesVectorAdaptor::__getitem__(int i) { return MPIdSamples(&(o->at(i))); };
void MPIdSamplesVectorAdaptor::__setitem__(int i, MPIdSamples val) { o->at(i) = *(val.o); };
void MPIdSamplesVectorAdaptor::append(MPIdSamples val) { o->push_back(*(val.o)); };






