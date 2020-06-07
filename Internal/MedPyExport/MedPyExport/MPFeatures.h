#ifndef __MED__MPFEATURES__H__
#define __MED__MPFEATURES__H__

#include "MedPyCommon.h"
#include "MPMat.h"
#include "MPSamples.h"
#include "MPSerializableObject.h"

class MedFeatures;
class FeatureAttr;

class MPFeatureAttr {
	bool o_owned = true;
public:
	MEDPY_IGNORE(FeatureAttr* o);
	MPFeatureAttr();
	MEDPY_IGNORE(MPFeatureAttr(const MPFeatureAttr& other));
	~MPFeatureAttr();
	MEDPY_IGNORE(MPFeatureAttr(FeatureAttr* ptr));
	void MEDPY_SET_normalized(bool _normalized);
	bool MEDPY_GET_normalized();
	void MEDPY_SET_imputed(bool _imputed);
	bool MEDPY_GET_imputed();
	void MEDPY_SET_denorm_mean(float _denorm_mean);
	float MEDPY_GET_denorm_mean();
	void MEDPY_SET_denorm_sdv(float _denorm_sdv);
	float MEDPY_GET_denorm_sdv();


	MPSerializableObject asSerializable();
};

class MPStringFeatureAttrMapAdaptor {
	bool o_owned = true;
	std::map<std::string, FeatureAttr>* o;
public:
	MPStringFeatureAttrMapAdaptor();
	MEDPY_IGNORE(MPStringFeatureAttrMapAdaptor(const MPStringFeatureAttrMapAdaptor& other));
	MEDPY_IGNORE(MPStringFeatureAttrMapAdaptor(std::map<std::string, FeatureAttr>* ptr));
	~MPStringFeatureAttrMapAdaptor();
	int __len__();
	MPFeatureAttr __getitem__(std::string key);
	void __setitem__(std::string key, MPFeatureAttr& val);
	std::vector<std::string> keys();
	MEDPY_IGNORE(MPStringFeatureAttrMapAdaptor& operator=(const MPStringFeatureAttrMapAdaptor& other));
};


class MPFeatures {
	bool o_owned = true;
public:
	MEDPY_IGNORE(MedFeatures* o);

	MPFeatures();
	MPFeatures(int _time_unit);
	MEDPY_IGNORE(MPFeatures(MedFeatures* from_ptr));
	~MPFeatures();
	MPFeatures(const MPFeatures& other);
	
	MPStringVecFloatMapAdaptor MEDPY_GET_data();
	void MEDPY_GET_weights(MEDPY_NP_OUTPUT(float** float_out_buf, unsigned long long* float_out_buf_len));
	MPSampleVectorAdaptor MEDPY_GET_samples();
	MPIntPairIntIntMapAdaptor MEDPY_GET_pid_pos_len();

	MPStringFeatureAttrMapAdaptor MEDPY_GET_attributes();
	MPStringUOSetStringMapAdaptor MEDPY_GET_tags();
	//map<string, unordered_set<string> > tags; ///< a set of tags per feature

	int MEDPY_GET_time_unit();
	void MEDPY_SET_time_unit(int new_time_unit);

	static int MEDPY_GET_global_serial_id_cnt();
	static void MEDPY_SET_global_serial_id_cnt(int newval);

	void clear();
	void set_time_unit(int _time_unit);
	std::vector<std::string> get_feature_names();

	void get_as_matrix(MPMat& mat) const;
	void get_as_matrix(MPMat& mat, vector<string> names) const;
	void get_as_matrix(MPMat& mat, const vector<string> names, MEDPY_NP_INPUT(int* int_in_buf, unsigned long long int_in_buf_len)) const;

	void set_as_matrix(MPMat& mat);
	void append_samples(MPIdSamples& in_samples);
	void append_samples(MPSamples& in_samples);
	void insert_samples(MPIdSamples& in_samples, int index);
	void init_all_samples(MPIdSamplesVectorAdaptor &in_samples);
	void init_pid_pos_len();
	
	int get_pid_pos(int pid) const;
	int get_pid_len(int pid) const;
	unsigned int get_crc();
	void print_csv() const;
	void get_samples(MPSamples& outSamples) const;
	int get_max_serial_id_cnt() const;
	int write_as_csv_mat(const string &csv_fname) const;
	int read_from_csv_mat(const string &csv_fname);

	int filter(std::vector<std::string>& selectedFeatures);
	int version();
	MPSerializableObject asSerializable();

	void split_by_fold(MPFeatures& outMatrix, int iFold, bool isLearning);
};

#endif //!__MED__MPFEATURES__H__
