#ifndef __MP_SAMPLES_H
#define __MP_SAMPLES_H

#include "MedPyCommon.h"
#include "MPPandasAdaptor.h"
#include "MPSerializableObject.h"



class MedSample;
class MedIdSamples;
class MedSamples;
class MPFeatures;

class MPSample  {
	bool o_owned = true;
public:
	MEDPY_IGNORE(MedSample *o);
	int MEDPY_GET_id();
	void MEDPY_SET_id(int new_id) ;
	int MEDPY_GET_split();
	void MEDPY_SET_split(int new_sp);
	int MEDPY_GET_time();
	void MEDPY_SET_time(int new_time);
	int MEDPY_GET_outcome();
	void MEDPY_SET_outcome(int new_outcome);
	int MEDPY_GET_outcomeTime();
	void MEDPY_SET_outcomeTime(int new_outcome_time);

	void MEDPY_GET_prediction(MEDPY_NP_OUTPUT(float** out_predbuf, unsigned long long* out_predbuf_len));
	void MEDPY_SET_prediction(MEDPY_NP_INPUT(float* in_predbuf, unsigned long long in_predbuf_len));

/*
	map<string, float> attributes;	///< Attribute(s) - empty if non given
	map<string, string> str_attributes;	///< Attribute(s) - empty if non given

*/
	MPSample();
	~MPSample();
	
	MEDPY_IGNORE(MPSample(const MPSample& other));
	MEDPY_IGNORE(MPSample(MedSample* from_ptr));

	void print_(const string prefix);
	void print_();

	int parse_from_string(string &s, int time_unit);
	void write_to_string(string &s, int time_unit);

//	int parse_from_string(string &s, map <string, int> & pos, vector<int>& pred_pos, map<string, int>& attr_pos, int time_unit);
//	void write_to_string(string &s, const vector<string>& attr, const vector<string>& str_attr, int time_unit);

	MPSample __copy__();

	MPSerializableObject asSerializable();
};

class MPSampleVectorAdaptor {
	bool o_owned = true;
public:
	MEDPY_IGNORE(vector<MedSample>* o);
	MPSampleVectorAdaptor();
	MEDPY_IGNORE(MPSampleVectorAdaptor(const MPSampleVectorAdaptor& other));
	MEDPY_IGNORE(MPSampleVectorAdaptor(vector<MedSample>* ptr));
	~MPSampleVectorAdaptor();
	int __len__();
	MPSample __getitem__(int i);
	void __setitem__(int i, MPSample& val);
	void append(MPSample& val);
	void append_vec(MPSampleVectorAdaptor& other);
	void override_splits(int nfolds);
	int nSplits();
};


class MPIdSamples {
	bool o_owned = true;
public:
	MEDPY_IGNORE(MedIdSamples* o);
	int MEDPY_GET_id();
	void MEDPY_SET_id(int new_id);
	int MEDPY_GET_split();
	void MEDPY_SET_split(int new_id);

	MPSampleVectorAdaptor MEDPY_GET_samples();

	MPIdSamples(int _id);

	MEDPY_IGNORE(MPIdSamples(MedIdSamples* ptr));
	MEDPY_IGNORE(MPIdSamples(const MPIdSamples& other));

	MPIdSamples();
	~MPIdSamples();
	void set_split(int _split);
	bool same_as(MPIdSamples &other, int mode);
	MPSerializableObject asSerializable();
};

class MPIdSamplesVectorAdaptor {
	bool o_owned = true;
public:
	MEDPY_IGNORE(vector<MedIdSamples>* o);
	MPIdSamplesVectorAdaptor();
	MEDPY_IGNORE(MPIdSamplesVectorAdaptor(const MPIdSamplesVectorAdaptor& other));
	MPIdSamplesVectorAdaptor(vector<MedIdSamples>* ptr);
	~MPIdSamplesVectorAdaptor();
	int __len__();
	MPIdSamples __getitem__(int i);
	void __setitem__(int i, MPIdSamples val);
	void append(MPIdSamples val);
};

class MPSamples {
	bool o_owned = true;
public:
	MEDPY_IGNORE(MedSamples* o);

	int MEDPY_GET_time_unit();
	void MEDPY_SET_time_unit(int new_time_unit);

	MPIdSamplesVectorAdaptor MEDPY_GET_idSamples();

	MPSamples();
	MEDPY_IGNORE(MPSamples(const MPSamples& other));
	MEDPY_IGNORE(MPSamples(MedSamples* ptr));
	~MPSamples();

	void clear();

	int insert_preds(MPFeatures& featuresData);
	void get_ids(MEDPY_NP_OUTPUT(int** ids, unsigned long long* num_ids));

	void append(MPSamples& newSamples);
	void subtract(MPSamples& _dont_use);
	void split_train_test(MPSamples& train, MPSamples& test, float p_test);
	int read_from_bin_file(const string& file_name);
	int write_to_bin_file(const string& file_name);

	int read_from_file(const string& file_name);
	int write_to_file(const string& fname);
	
	void get_preds(MEDPY_NP_OUTPUT(float** preds_buf, unsigned long long* preds_buf_len));
	void get_y(MEDPY_NP_OUTPUT(float** y_buf, unsigned long long* y_buf_len));
	void get_categs(MEDPY_NP_OUTPUT(float** categs_buf, unsigned long long* categs_buf_len));
	
	MPSampleVectorAdaptor export_to_sample_vec();
	void import_from_sample_vec(MPSampleVectorAdaptor& vec_samples, bool allow_split_inconsistency = false);

	void MEDPY__from_df(MPPandasAdaptor& pandas_df);
	MPPandasAdaptor MEDPY__from_df_adaptor();
	MPPandasAdaptor MEDPY__to_df();


	void sort_by_id_date();
	void normalize();
	bool same_as(MPSamples &other, int mode);
	int nSamples();

	MEDPY_DOC(nSplits, "nSplits() -> int\n"
		"Return number of splits, also check mismatches between idSample and internal MedSamples and set idSamples.split if missing");
	int nSplits();
	
	int get_predictions_size();

	//int get_all_attributes(vector<string>& attributes, vector<string>& str_attributes);
	//temporary solution to inability to return two string lists:
	std::vector<string> get_attributes();
	std::vector<string> get_str_attributes();

	void dilute(float prob);
	void insertRec(int pid, int time, float outcome, int outcomeTime);
	void insertRec(int pid, int time, float outcome, int outcomeTime, float pred);
	void insertRec(int pid, int time);
	int version();
	MPSerializableObject asSerializable();

	void override_splits(int nfolds);

	/// sets sample for all the patient in specific time
	void set_samples(MEDPY_NP_INPUT(int * patients, unsigned long long patient_size), int _time = -1);
};




#endif //!__MP_SAMPLES_H
