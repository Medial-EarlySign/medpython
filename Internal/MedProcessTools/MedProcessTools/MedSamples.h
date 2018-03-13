// Classes for holding samples

#ifndef _MED_SAMPLES_H_
#define _MED_SAMPLES_H_

#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedTime/MedTime/MedTime.h"
#include "InfraMed/InfraMed/MedRepositoryType.h"
#include <unordered_set>

class MedFeatures;

//.......................................................................................
/**  MedSample represents a signle sample: id + time (date) <br>
*	 Additional (optinal) entries: outcome, outcome_date, split and prediction <br>
*/
//.......................................................................................

class MedSample : public SerializableObject {
public:
	int id = -1;			///< Patient id
	int split = -1;			///< Cross-validation split. -1 if not given. Proper use is to set the same split for all samples of a given id (MedIdSamples), but this is not enforced.
	int time = 0;			///< Time (Date)
	float outcome = 0;		///< Outcome
	int outcomeTime = 0;	///< Outcome time (date)
	vector<float> prediction;	///< Prediction(s) - empty if non given
	 
	/// <summary> Constructor </summary>
	MedSample() { prediction.clear(); }
	/// <summary> Destructor </summary>
	~MedSample() { prediction.clear(); }

	/// <summary> printing the sample (with a prefix) </summary>
	void print(const string prefix);

	/// <summary> printing the sample </summary>
	void print() { print(""); }

	/// <summary>
	/// Get sample from tab-delimited string, in old or new format (<split> and <prediction> optional, <predictions> can be several numbers (tab delimited)) <br>
	/// old format: EVENT <id> <time> <outcome> <outcomeLen(dummy here)> <outcomeTime> <split> <predictions> <br>
	/// new format: SAMPLE <id> <time> <outcome> <outcomeTime> <split> <predictions> <br>
	/// </summary>
	/// <returns> 0 upon success, -1 if string does not fit any of the formats </returns>
	int parse_from_string(string &s);
	/// <summary>
	/// Get sample from tab-delimited string, where pos indicate the position of each field (fields are id,date,outcome,outcome_date,split,pred) 
	/// if pos is empty, check old and new formats
	/// </summary>
	/// <returns> 0 upon success, -1 upon failure to parse </returns>
	int parse_from_string(string &s, map <string, int> & pos);
	/// <summary> Write to string in new format </summary>
	void write_to_string(string &s);

	// Serialization
	ADD_SERIALIZATION_FUNCS(id, split, time, outcome, outcomeTime, prediction)
};

/// <summary> Comparison functions for sorting by prediction value </summary>
inline bool comp_sample_pred(const MedSample &pr1, const MedSample &pr2) {
	return pr1.prediction[0] < pr2.prediction[0];
}

/// <summary> Comparison functions for sorting by id and date </summary>
inline bool comp_sample_id_time(const MedSample &pr1, const MedSample &pr2) {
	if (pr1.id == pr2.id)
		return pr1.time < pr2.time;
	else
		return pr1.id < pr2.id;
}

//.......................................................................................
/**  MedIdSamples represent a collection of samples of a given id <br>
*	 Additional (optinal) entries: split <br>
*/
//.......................................................................................
class MedIdSamples : public SerializableObject {
public:
	int id = -1;		///< Patient id
	/// Split for cross-validation. Note that nothing forces the id and split of each MedSample to be the same as that of MedIdSamples, though anything else is an improper use, and not guaranteed to work.
	int split = -1;		
	vector<MedSample> samples;		///< MedSamples for the given id

	/// <summary> Constructor with id </summary>
	MedIdSamples(int _id) { id = _id; split = -1; samples.clear(); }
	/// <summary> Constructor without id </summary>
	MedIdSamples() { id = -1; split = -1; samples.clear(); }
	
	/// <summary> Set split and export to all MedSample entries. </summary> 
	void set_split(int _split) { split = _split; for (auto& s : samples) s.split = _split; }

	/// <summary> Comparison function : mode 0 requires equal id/time, mode 1 requires equal outcome info, mode 2 also compares split and prediction </summary>
	/// <returns> true if equal , false otherwise </returns>
	bool same_as(MedIdSamples &other, int mode);

	// Serialization
	ADD_SERIALIZATION_FUNCS(id, split, samples)

};

/// <summary> Comparison function for sorting by id </summary>
inline bool comp_patient_id_time(const MedIdSamples &pr1, const MedIdSamples &pr2) {
	return pr1.id < pr2.id;
}

//.......................................................................................
/**  MedSamples represent a collection of samples per different id <br>
*   The data is conatined in a vector of MedIdSamples
*/
//.......................................................................................

class MedSamples : public SerializableObject {
public:
	int time_unit = MedTime::Date;	///< The time unit in which the samples are given. Default: Date
	vector<MedIdSamples> idSamples; ///< The vector of MedIdSamples

	/// <summary> Constructor. init time_unit according to global med_rep_type </summary>
	MedSamples() { time_unit = med_rep_type.basicTimeUnit; }
	~MedSamples() {}

	/// <summary> Clear data and init time_unit according to global med_rep_type </summary>
	void clear() { time_unit = MedTime::Date; idSamples.clear(); }

	/// <summary>
	/// Extract predictions from MedFeatures and insert to corresponding samples <br>
	/// Samples in MedFeatures are assumed to be of the same size and order as in MedSamples
	/// </summary>
	/// <returns> -1 if samples and features do not match in length, 0 upon success </returns>
	int insert_preds(MedFeatures& featuresData);

	/// <summary> Get all patient ids </summary>
	void get_ids(vector<int>& ids);

	/// <summary> Append new MedIdSamples at the end of current ones </summary>
	void append(MedSamples& newSamples) { idSamples.insert(idSamples.end(), newSamples.idSamples.begin(), newSamples.idSamples.end()); }

	/// <summary> Read from bin file</summary>
	/// <returns> -1 upon failure to open file, 0 upon success </returns>
	int read_from_bin_file(const string& file_name) { return SerializableObject::read_from_file(file_name); }
	/// <summary>  Write to bin file  </summary>
	/// <returns>  -1 upon failure to open file, 0 upon success  </returns>
	int write_to_bin_file(const string& file_name) { return SerializableObject::write_to_file(file_name); }

	/// <summary>
	/// Read from text file. <br>
	/// If a line starting with EVENT_FIELDS (followed by tabe-delimeted field names : id,date,outcome,outcome_date,split,pred) appears before the data lines, it is used to determine
	/// fields positions, otherwise - old or new formats are used. 
	/// </summary>
	/// <returns>  -1 upon failure to open file, 0 upon success </returns>
	int read_from_file(const string& file_name);

	/// <summary>  Write to text file in new format  </summary>
	/// <returns> -1 upon failure to open file, 0 upon success </returns>
	int write_to_file(const string &fname);

	/// <summary> Extract a single vector of concatanated predictions </summary>
	void get_preds(vector<float>& preds);
	/// <summary> Extract a vector of all outcomes  </summary>
	void get_y(vector<float>& y);
	/// <summary>  Get a list of all categories (different values) appearing in the outcome </summary>
	void get_categs(vector<float> &categs); 
	/// <summary> Get all MedSamples as a single vector </summary>
	void export_to_sample_vec(vector<MedSample> &vec_samples);

	/// <summary> Sort by id and then date </summary>
	void sort_by_id_date(); 
	/// <summary> Make sure that : (1) every pid has one idSample at most and (2) everything is sorted </summary>
	void normalize(); 

	/// <summary> Comparison function : mode 0 requires equal id/time, mode 1 requires equal outcome info, mode 2 also compares split and prediction </summary>
	/// <returns> true if equal , false otherwise </returns>
	bool same_as(MedSamples &other, int mode);

	/// <summary> Return number of samples </summary>
	int nSamples();


	// <summary> given a probability dilution prob, dilute current samples </summary>
	void dilute(float prob);

	/// <summary>  API's for online insertions : main use case is a single time point for prediction per pid </summary>
	void insertRec(int pid, int time, float outcome, int outcomeTime);
	void insertRec(int pid, int time, float outcome, int outcomeTime, float pred);
	void insertRec(int pid, int time) { insertRec(pid, time, -1, 0); }

	// Version for serialization
	int version() { return  1; };

	//Serialization, version 1: Added version, model_features, features_count to serialization
	ADD_SERIALIZATION_FUNCS(time_unit, idSamples)
};

namespace medial {
	namespace print {
		void print_samples_stats(const vector<MedSample> &samples, const string &log_file = "");
		void print_samples_stats(const MedSamples &samples, const string &log_file = "");
		void print_by_year(const vector<MedSample> &data_records, int year_bin_size, bool unique_ids = false,
			bool take_prediction_time = true, const string &log_file = "");
		void print_by_year(const MedSamples &data_records, int year_bin_size, bool unique_ids = false,
			bool take_prediction_time = true, const string &log_file = "");
	}
	namespace process {
		void down_sample(MedSamples &samples, double take_ratio);
	}
}

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(MedSample)
MEDSERIALIZE_SUPPORT(MedIdSamples)
MEDSERIALIZE_SUPPORT(MedSamples)

#endif
