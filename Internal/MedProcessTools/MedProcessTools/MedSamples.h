// Classes for holding samples

#ifndef __MED_SAMPLES_H__
#define __MED_SAMPLES_H__


#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedTime/MedTime/MedTime.h>
#include <unordered_set>
#include <json/json.hpp>

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
	map<string, float> attributes;	///< Attribute(s) - empty if non given
	map<string, string> str_attributes;	///< Attribute(s) - empty if non given
	json jrec;	///< a json record that can be built along side the sample to contain any information in a nice json format

	/// <summary> Constructor </summary>
	MedSample() { prediction.clear(); }
	MedSample(int _id, int _time, int _outcome, int _outcomeTime) { id = _id; time = _time; outcome = _outcome; outcomeTime = _outcomeTime; prediction.clear(); }
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
	int parse_from_string(string &s, int time_unit);
	/// <summary>
	/// Get sample from tab-delimited string, where pos indicate the position of each field (fields are id,date,outcome,outcome_date,split,pred or pred_,attr_NAME) 
	/// if pos is empty, check old and new formats
	/// </summary>
	/// <returns> 0 upon success, -1 upon failure to parse </returns>
	int parse_from_string(string &s, const map <string, int> & pos, const vector<int>& pred_pos, const map<string, int>& attr_pos,
		const map<string, int>& str_attr_pos, int time_unit, int raw_format, const string &delimeter = "\t");
	int parse_from_string(const vector<string> &fields, const map<string, int> & pos, const vector<int>& pred_pos, const map<string, int>& attr_pos,
		const map<string, int>& str_attr_pos, int time_unit, int raw_format, const string &delimeter);

	/// <summary> Write to string in new format </summary>
	void write_to_string(string &s, int time_unit, bool write_attrib = true, const string &delimeter = "\t", int pred_precision = -1) const;

	/// <summary> Get sample attributes </summary>
	int get_all_attributes(vector<string>& attributes, vector<string>& str_attributes) const;

	// Serialization
	ADD_CLASS_NAME(MedSample)
		ADD_SERIALIZATION_FUNCS(id, split, time, outcome, outcomeTime, prediction, attributes, str_attributes)
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

	/// <summary> get a vector of all times for the pid
	void get_times(vector<int> &times) const { times.clear(); for (auto &s : samples) times.push_back(s.time); }

	// Serialization
	ADD_CLASS_NAME(MedIdSamples)
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

class MedSamples final : public SerializableObject {
public:
	int time_unit = MedTime::Date;	///< The time unit in which the samples are given. Default: Date
	vector<MedIdSamples> idSamples; ///< The vector of MedIdSamples
	int raw_format = 0; // read times as is, no conversions

	/// <summary> Constructor. init time_unit according to default </summary>
	MedSamples() { time_unit = global_default_time_unit; }
	~MedSamples() {}

	/// <summary> Clear data and init time_unit according to default </summary>
	void clear() { time_unit = global_default_time_unit; idSamples.clear(); }

	/// <summary>
	/// Extract predictions from MedFeatures and insert to corresponding samples <br>
	/// Samples in MedFeatures are assumed to be of the same size and order as in MedSamples
	/// </summary>
	/// <returns> -1 if samples and features do not match in length, 0 upon success </returns>
	int insert_preds(MedFeatures& featuresData);

	/// <summary>
	/// Copy attributes from MedSample vector. This function is mainly used to <br>
	/// Extract post processors results from MedFeatures and insert to corresponding samples <br>
	/// Samples in MedFeatures are assumed to be of the same size and order as in MedSamples
	/// </summary>
	/// <returns> -1 if samples and features do not match in length, 0 upon success </returns>
	int copy_attributes(const vector<MedSample>& samples);

	/// <summary> Get all patient ids </summary>
	void get_ids(vector<int>& ids) const;

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
	int read_from_file(const string& file_name, bool sort_rows = true);

	/// <summary>  Write to text file in new format  </summary>
	/// <returns> -1 upon failure to open file, 0 upon success </returns>
	int write_to_file(const string &fname, int pred_precision=-1, bool print_attributes =true);
	void write_to_file(ofstream& of, int pred_precision, bool print_attributes, bool print_header);

	/// <summary> Extract a single vector of concatanated predictions </summary>
	void get_preds(vector<float>& preds) const;
	void get_preds_channel(vector<float>& preds, int channel);
	/// <summary> Extract a vector of all outcomes  </summary>
	void get_y(vector<float>& y) const;
	/// <summary>  Get a list of all categories (different values) appearing in the outcome </summary>
	void get_categs(vector<float> &categs) const;
	/// <summary> get a vector corresponding to given attr (name should include attr_) </summary>
	void get_attr_values(const string& attr_name, vector<float>& values) const;
	/// <summary> Get all MedSamples as a single vector </summary>
	void export_to_sample_vec(vector<MedSample> &vec_samples) const;
	/// <summary> Set MedSamples from a single vector </summary>
	void import_from_sample_vec(const vector<MedSample> &vec_samples, bool allow_split_inconsistency = false);

	/// <summary> Sort by id and then date </summary>
	void sort_by_id_date();
	/// <summary> Make sure that : (1) every pid has one idSample at most and (2) everything is sorted </summary>
	void normalize();

	/// <summary> Comparison function : mode 0 requires equal id/time, mode 1 requires equal outcome info, mode 2 also compares split and prediction </summary>
	/// <returns> true if equal , false otherwise </returns>
	bool same_as(MedSamples &other, int mode);

	/// <summary> Return number of samples </summary>
	int nSamples() const;

	/// <summary> Return number of splits, also check mismatches between idSample and internal MedSamples and set idSamples.split if missing  </summary>
	int nSplits();

	/// <summary> Get predictions vector size. Return -1 if not-consistent </summary>
	int get_predictions_size(int& nPreds);

	/// <summary> Get all attributes. Return -1 if not-consistent </summary>
	int get_all_attributes(vector<string>& attributes, vector<string>& str_attributes) const;

	/// <summary> given a probability dilution prob, dilute current samples </summary>
	void dilute(float prob);

	/// <summary> will dilute 0 labeled samples (traditionally controls) with p0, and all the rest with p1 </summary>
	void binary_dilute(float p0, float p1);

	/// <summary> removing all ids that appear in _dont_include from the current samples
	void subtract(MedSamples &_dont_include);

	/// <summary> gets p_test and splits by id , p_test of the ids into test, and the rest into train
	void split_train_test(MedSamples &train, MedSamples &test, float p_test);

	/// <summary> gets a split number and splits samples to lists in/off the split
	void split_by_split(MedSamples &in_split, MedSamples &off_split, int split);

	/// <summary> adding splits to the samples given in an external file
	void add_splits_from_file(string f_splits);

	/// <summary> initializing all jrecs to contain pid and time
	void init_all_jrecs();


	void flatten(vector<MedSample> &flat) const;

	/// <summary>  API's for online insertions : main use case is a single time point for prediction per pid </summary>
	void insertRec(int pid, int time, float outcome, int outcomeTime);
	void insertRec(int pid, int time, float outcome, int outcomeTime, float pred);
	void insertRec(int pid, int time) { insertRec(pid, time, -1, 0); }

	//Serialization, version 1: Added version, model_features, features_count to serialization
	//				 version 2: Added attributes
	ADD_CLASS_NAME(MedSamples)
		ADD_SERIALIZATION_FUNCS(time_unit, idSamples)
};

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	*  \brief print namespace
	*/
	namespace print {
		/// \brief print samples stats
		void print_samples_stats(const vector<MedSample> &samples, const string &log_file = "");
		/// \brief print samples stats
		void print_samples_stats(const MedSamples &samples, const string &log_file = "");
		/// \brief print samples stats by group
		void print_by(const vector<MedSample> &data_records, const vector<string> &groups,
			bool unique_ids = false, const string &log_file = "");
		/// \brief print samples stats by group
		void print_by(const MedSamples &data_records, const vector<string> &groups,
			bool unique_ids = false, const string &log_file = "");
		/// \brief print samples stats by year
		void print_by_year(const vector<MedSample> &data_records, int year_bin_size, bool unique_ids = false,
			bool take_prediction_time = true, const string &log_file = "");
		/// \brief print samples stats by year
		void print_by_year(const MedSamples &data_records, int year_bin_size, bool unique_ids = false,
			bool take_prediction_time = true, const string &log_file = "");
	}
	/*!
	*  \brief process namespace
	*/
	namespace process {
		/// \brief down sammling
		void down_sample(MedSamples &samples, double take_ratio, bool with_repeats = false);
		/// \brief down sample by selecting from pids
		void down_sample_by_pid(MedSamples &samples, double take_ratio, bool with_repeats = false);

		/// \brief down sammling
		void down_sample(MedSamples &samples, int no_more_than, bool with_repeats = false);
		/// \brief down sample by selecting from pids
		void down_sample_by_pid(MedSamples &samples, int no_more_than, bool with_repeats = false);
	}

	namespace stats {
		double kaplan_meir_on_samples(const MedSamples &incidence_samples, int time_period, const vector<pair<int, int>> *filtered_idx = NULL);
		double kaplan_meir_on_samples(const vector<MedSample> &incidence_samples, int time_unit, int time_period, const vector<int> *filtered_idx = NULL);
	}
}

//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(MedSample)
MEDSERIALIZE_SUPPORT(MedIdSamples)
MEDSERIALIZE_SUPPORT(MedSamples)

#endif
