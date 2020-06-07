#include "MPFeatures.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"

MPFeatures::MPFeatures() { o = new MedFeatures(); }
MPFeatures::MPFeatures(int _time_unit) { o = new MedFeatures(_time_unit); }
MPFeatures::MPFeatures(MedFeatures* from_ptr) { o_owned = false; o = from_ptr; }
MPFeatures::MPFeatures(const MPFeatures& other)
{
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new MedFeatures();
		*o = *other.o;
	}
};

MPFeatures::~MPFeatures() { if (o_owned) delete o; }

int MPFeatures::MEDPY_GET_global_serial_id_cnt() { return MedFeatures::global_serial_id_cnt; };
void MPFeatures::MEDPY_SET_global_serial_id_cnt(int newval) { MedFeatures::global_serial_id_cnt = newval; };

MPStringVecFloatMapAdaptor MPFeatures::MEDPY_GET_data() { return MPStringVecFloatMapAdaptor(&(o->data)); };
void MPFeatures::MEDPY_GET_weights(MEDPY_NP_OUTPUT(float** float_out_buf, unsigned long long* float_out_buf_len)) {
	vector_to_buf(o->weights, float_out_buf, float_out_buf_len);
}
MPSampleVectorAdaptor MPFeatures::MEDPY_GET_samples() { return MPSampleVectorAdaptor(&(o->samples)); };

MPIntPairIntIntMapAdaptor MPFeatures::MEDPY_GET_pid_pos_len() { return MPIntPairIntIntMapAdaptor(&(o->pid_pos_len)); };

int MPFeatures::MEDPY_GET_time_unit() { return o->time_unit; }
void MPFeatures::MEDPY_SET_time_unit(int new_time_unit) { o->set_time_unit(new_time_unit); };

void MPFeatures::clear() { o->clear(); }
void MPFeatures::set_time_unit(int _time_unit) { o->set_time_unit(_time_unit); };
std::vector<std::string> MPFeatures::get_feature_names() {
	vector<string> names;
	o->get_feature_names(names);
	return names;
}

void MPFeatures::get_as_matrix(MPMat& mat) const { o->get_as_matrix(*(mat.o)); };
void MPFeatures::get_as_matrix(MPMat& mat, vector<string> names) const { o->get_as_matrix(*(mat.o), names); };
void MPFeatures::get_as_matrix(MPMat& mat, const vector<string> names, MEDPY_NP_INPUT(int* int_in_buf, unsigned long long int_in_buf_len)) const {
	vector<int> idx;
	buf_to_vector(int_in_buf, int_in_buf_len, idx);
	o->get_as_matrix(*(mat.o), names, idx);
}

void MPFeatures::set_as_matrix(MPMat& mat) {
	o->set_as_matrix(*(mat.o));
}

void MPFeatures::append_samples(MPIdSamples& in_samples) {
	o->append_samples(*(in_samples.o));
}

void MPFeatures::insert_samples(MPIdSamples& in_samples,int index) {
	o->insert_samples(*(in_samples.o) ,index);
}
void MPFeatures::init_all_samples(MPIdSamplesVectorAdaptor &in_samples) {
	o->init_all_samples(*(in_samples.o));
}

void MPFeatures::append_samples(MPSamples& in_samples) {
	for (auto &idSample : (*(in_samples.o)).idSamples)
	{
		o->append_samples(idSample);
	}
}


void MPFeatures::init_pid_pos_len() { o->init_pid_pos_len(); }

int MPFeatures::get_pid_pos(int pid) const { return o->get_pid_pos(pid); };
int MPFeatures::get_pid_len(int pid) const { return o->get_pid_len(pid); };
unsigned int MPFeatures::get_crc() { return o->get_crc(); };
void MPFeatures::print_csv() const { o->print_csv(); };

void MPFeatures::get_samples(MPSamples& outSamples) const {
	o->get_samples(*(outSamples.o));
}

int MPFeatures::get_max_serial_id_cnt() const { return o->get_max_serial_id_cnt();}

int MPFeatures::write_as_csv_mat(const string &csv_fname) const { return val_or_exception(o->write_as_csv_mat(csv_fname), string("Cannot write features to csv file ") + csv_fname); };
int MPFeatures::read_from_csv_mat(const string &csv_fname) { return val_or_exception(o->read_from_csv_mat(csv_fname), string("Cannot read features from csv file ") + csv_fname); };
int MPFeatures::filter(std::vector<std::string>& selectedFeatures) { 
	unordered_set<string> selectedFeatures_set;
	for (auto& sf : selectedFeatures)
		selectedFeatures_set.insert(sf);
	return o->filter(selectedFeatures_set);
}
int MPFeatures::version() { return o->version(); }
MPSerializableObject MPFeatures::asSerializable() { return MPSerializableObject(o); }

MPStringFeatureAttrMapAdaptor MPFeatures::MEDPY_GET_attributes() {
	return MPStringFeatureAttrMapAdaptor(&(o->attributes));
}

MPStringUOSetStringMapAdaptor MPFeatures::MEDPY_GET_tags() {
	return MPStringUOSetStringMapAdaptor(&(o->tags));
}

// Get learning/test matrix
void MPFeatures::split_by_fold(MPFeatures& outMatrix, int iFold, bool isLearning) {
	vector<string> feature_names;
	o->get_feature_names(feature_names);
	for (string& name : feature_names)
		outMatrix.o->data[name].clear();
	outMatrix.o->samples.clear();

	for (auto& attr : o->attributes)
		outMatrix.o->attributes[attr.first] = attr.second;

	for (unsigned int i = 0; i<o->samples.size(); i++) {
		auto& sample = o->samples[i];
		if ((isLearning && sample.split != iFold) || ((!isLearning) && sample.split == iFold)) {
			outMatrix.o->samples.push_back(sample);
			for (string& name : feature_names)
				outMatrix.o->data[name].push_back(o->data[name][i]);
		}
	}
}



/*******************  FEATURE ATTR **********************/


MPFeatureAttr::MPFeatureAttr(FeatureAttr* ptr) { o_owned = false; o = ptr; };
MPFeatureAttr::MPFeatureAttr() { o = new FeatureAttr(); };
MPFeatureAttr::MPFeatureAttr(const MPFeatureAttr& other) {
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new FeatureAttr();
		*o = *other.o;
	}
}
MPFeatureAttr::~MPFeatureAttr() { if (o_owned) delete o; };


void MPFeatureAttr::MEDPY_SET_normalized(bool _notmalized) { o->normalized = _notmalized; };
bool MPFeatureAttr::MEDPY_GET_normalized() { return o->normalized; };

void MPFeatureAttr::MEDPY_SET_imputed(bool _imputed) { o->imputed = _imputed; };
bool MPFeatureAttr::MEDPY_GET_imputed() { return o->imputed; };

void MPFeatureAttr::MEDPY_SET_denorm_mean(float _denorm_mean) { o->denorm_mean = _denorm_mean; };
float MPFeatureAttr::MEDPY_GET_denorm_mean() { return o->denorm_mean; };
void MPFeatureAttr::MEDPY_SET_denorm_sdv(float _denorm_sdv) { o->denorm_sdv = _denorm_sdv; };
float MPFeatureAttr::MEDPY_GET_denorm_sdv() { return o->denorm_sdv; };


MPSerializableObject MPFeatureAttr::asSerializable() { return MPSerializableObject(o); }

MPStringFeatureAttrMapAdaptor::MPStringFeatureAttrMapAdaptor() { o = new std::map<std::string, FeatureAttr>(); };;
MPStringFeatureAttrMapAdaptor::MPStringFeatureAttrMapAdaptor(const MPStringFeatureAttrMapAdaptor& other) {
	o_owned = other.o_owned;
	if (!other.o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, FeatureAttr>();
		*o = *other.o;
	}
};
MPStringFeatureAttrMapAdaptor::MPStringFeatureAttrMapAdaptor(std::map<std::string, FeatureAttr>* ptr) {  o_owned = false; o = ptr; };
MPStringFeatureAttrMapAdaptor::~MPStringFeatureAttrMapAdaptor() { if (o_owned) delete o; };
int MPStringFeatureAttrMapAdaptor::__len__() { return (int)o->size(); };
MPFeatureAttr MPStringFeatureAttrMapAdaptor::__getitem__(std::string key) { return MPFeatureAttr(&(o->operator[](key))); };
void MPStringFeatureAttrMapAdaptor::__setitem__(std::string key, MPFeatureAttr& val) { o->operator[](key) = *(val.o); };
std::vector<std::string> MPStringFeatureAttrMapAdaptor::keys() {
	vector<string> ret;
	ret.reserve(o->size());
	for (const auto& rec : *o) ret.push_back(rec.first);
	return ret;
};

MPStringFeatureAttrMapAdaptor& MPStringFeatureAttrMapAdaptor::operator=(const MPStringFeatureAttrMapAdaptor& other)
{
	if (&other == this)
		return *this;
	o_owned = other.o_owned;
	if (!o_owned) {
		o = other.o;
	}
	else {
		o = new std::map<std::string, FeatureAttr>();
		*o = *(other.o);
	}
	return *this;
}

