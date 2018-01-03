// Temporary class for holding samples

#ifndef _MED_SAMPLES_H_
#define _MED_SAMPLES_H_

#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedTime/MedTime/MedTime.h"
#include "InfraMed/InfraMed/MedRepositoryType.h"
#include <unordered_set>

class MedFeatures;

//.......................................................................................
//.......................................................................................
// Samples represent id + time (date)
// Additional (optinal) entries - 
//		outcome related data: outcome + outcome_date
//		split : For cross validation. Proper use is to set the same split for all samples
//				of a given id (MedIdSamples), but this is not enforced. -1 for no-split
//		prediction : a vector of predictions (empty if non given)
//.......................................................................................
//.......................................................................................
// A single sample
class MedSample : public SerializableObject {
public:
	int id = -1;			// Patient id
	int split = -1;			// Cross-validation split. -1 if not given
	int time = 0;			// Time (Date)
	float outcome = 0;		// Outcome
	int outcomeTime = 0;	// Outcome time (date)
	vector<float> prediction;	// Predictions
	 
	// Constructor & Destructor
	MedSample() { prediction.clear(); }
	~MedSample() { prediction.clear(); }

	// Serialization
	ADD_SERIALIZATION_FUNCS(id, split, time, outcome, outcomeTime, prediction)

	// printing the sample (with an optional prefix)
	void print(const string prefix);
	void print() { print(""); }

	// Get sample from tab-delimited string, in old or new format (<split> and <prediction> optional, <predictions> can be several numbers (tab delimited))
	// old format: EVENT <id> <time> <outcome> <outcomeLen(dummy here)> <outcomeTime> <split> <predictions> 
	// new format: SAMPLE <id> <time> <outcome> <outcomeTime> <split> <predictions>
	int parse_from_string(string &s);
	// Get sample from tab-delimited string, where pos indicate the position of each field (fields are id,date,outcome,outcome_date,split,pred)
	int parse_from_string(string &s, map <string, int> & pos);
	// Write to string in new format
	int write_to_string(string &s);
};

// Comparison functions for sorting
inline bool comp_sample_pred(const MedSample &pr1, const MedSample &pr2) {
	return pr1.prediction[0] < pr2.prediction[0];
}

inline bool comp_sample_id_time(const MedSample &pr1, const MedSample &pr2) {
	if (pr1.id == pr2.id)
		return pr1.time < pr2.time;
	else
		return pr1.id < pr2.id;
}

//.......................................................................................
//.......................................................................................
// MedIdSamples represent a collection of samples of a given id
//		id = the patient id
//		split = the id's split in cross-validation
//		vector of MedSamples.
//		Note that nothing forces the id and split of each MedSample to be the same as that
//			of MedIdSamples, though anything else is an improper use, and not guaranteed
//			to work.
//.......................................................................................
//.......................................................................................
class MedIdSamples : public SerializableObject {
public:
	int id = -1;		// Patient id
	int split = -1;		// Split for cross-validation
	vector<MedSample> samples;		// List of samples for the given id

	// Constructors
	MedIdSamples(int _id) { id = _id; split = -1; samples.clear(); }
	MedIdSamples() { id = -1; split = -1; samples.clear(); }
	
	// Set split and export to all MedSample entries.
	void set_split(int _split) { split = _split; for (auto& s : samples) s.split = _split; }

	// Serialization
	ADD_SERIALIZATION_FUNCS(id, split, samples)

};

// Comparison function for sorting
inline bool comp_patient_id_time(const MedIdSamples &pr1, const MedIdSamples &pr2) {
	return pr1.id < pr2.id;
}

//.......................................................................................
//.......................................................................................
// MedSamples represent a collection of MedIdSamples
//		time_unit = the time unit in which the samples are given. Default: Date
//		vector of MedIdSamples.
//.......................................................................................
//.......................................................................................
class MedSamples : public SerializableObject {
public:
	int time_unit = MedTime::Date;	// time_unit = the time unit in which the samples are given. Default: Date
	vector<MedIdSamples> idSamples; 

	// Constructor
	MedSamples() { time_unit = med_rep_type.basicTimeUnit; }

	// Reset
	void clear() { time_unit = MedTime::Date; idSamples.clear(); }

	// Extract predictions from MedFeatures and insert to corresponding samples
	// Samples in MedFeatures are assumed to be of the same size and order as in MedSamples
	int insert_preds(MedFeatures& featuresData);

	// Get all patient ids
	void get_ids(vector<int>& ids);

	// Append new MedIdSamples at the end of current ones
	void append(MedSamples& newSamples) { idSamples.insert(idSamples.end(), newSamples.idSamples.begin(), newSamples.idSamples.end()); }

	// read/write to bin file
	int read_from_bin_file(const string& file_name) { return SerializableObject::read_from_file(file_name); }
	int write_to_bin_file(const string& file_name) { return SerializableObject::write_to_file(file_name); }

	// read from text file.
	// If a line starting with EVENT_FIELDS (followed by tabe-delimeted field names : id,date,outcome,outcome_date,split,pred) appears before the data lines, it is used to determine
	// fields positions, otherwise - old or new formats are used.
	int read_from_file(const string& file_name);
	// write to text file in new format
	int write_to_file(const string &fname);

	// Extract a single vector of concatanated predictions
	void get_preds(vector<float>& preds);
	// Extract a vector of all outcomes
	void get_y(vector<float>& y);
	// Get a list of all categories (different values) appearing in the outcome
	void get_categs(vector<float> &categs); 
	// Get all MedSamples as a single vector
	void export_to_sample_vec(vector<MedSample> &vec_samples);

	// Sort by id and then date
	void sort_by_id_date(); 
	// Make sure that : (1) every pid has one idSample at most and (2) everything is sorted
	void normalize(); 


	// Count samples
	int nSamples();


	// dilute : given a probability dilution prob, dilute current samples
	void dilute(float prob);

	// API's for online insertions : main use case is a single time point for prediction per pid
	int insertRec(int pid, int time, float outcome, int outcomeTime);
	int insertRec(int pid, int time, float outcome, int outcomeTime, float pred);
	int insertRec(int pid, int time) { return insertRec(pid, time, -1, 0); }

	// Version for serialization
	int version() { return  1; };

	//Serialization, version 1: Added version, model_features, features_count to serialization
	ADD_SERIALIZATION_FUNCS(time_unit, idSamples)
};


//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(MedSample)
MEDSERIALIZE_SUPPORT(MedIdSamples)
MEDSERIALIZE_SUPPORT(MedSamples)

#endif
