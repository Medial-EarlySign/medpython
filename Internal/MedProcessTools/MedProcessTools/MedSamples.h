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
// Samples represent id + date
// It may also include outcome related data
//.......................................................................................
//.......................................................................................
// A single sample
class MedSample : public SerializableObject {
public:
	int id = -1;		// Optional
	int split = -1;	// Optional
	int time = 0;
	float outcome = 0;
	int outcomeTime = 0;
	vector<float> prediction;

	MedSample() { prediction.clear(); }
	~MedSample() { prediction.clear(); }

	ADD_SERIALIZATION_FUNCS(id, split, time, outcome, outcomeTime, prediction)

	// print
	void print(const string prefix);
	void print() { print(""); }

	// parsing
	int parse_from_string(string &s);
	int parse_from_string(string &s, map <string, int> & pos);
	int write_to_string(string &s);
};

inline bool comp_sample_pred(const MedSample &pr1, const MedSample &pr2) {
	return pr1.prediction[0] < pr2.prediction[0];
}

inline bool comp_sample_id_time(const MedSample &pr1, const MedSample &pr2) {
	if (pr1.id == pr2.id)
		return pr1.time < pr2.time;
	else
		return pr1.id < pr2.id;
}

// A collection of samples of a given id
class MedIdSamples : public SerializableObject {
public:
	int id = -1;
	int split = -1;
	vector<MedSample> samples;

	// Constructors
	MedIdSamples(int _id) { id = _id; split = -1; samples.clear(); }
	MedIdSamples() { id = -1; split = -1; samples.clear(); }
	
	void set_split(int _split) { split = _split; for (auto& s : samples) s.split = _split; }
	ADD_SERIALIZATION_FUNCS(id, split, samples)

};

inline bool comp_patient_id_time(const MedIdSamples &pr1, const MedIdSamples &pr2) {
	return pr1.id < pr2.id;
}

// A collection of ids and relevant samples
class MedSamples : public SerializableObject {
public:
	int time_unit = MedTime::Date; // the time unit in which the samples are given. Default: Date
	vector<MedIdSamples> idSamples;

	// Constructor
	MedSamples() { time_unit = med_rep_type.basicTimeUnit; }

	void clear() { time_unit = MedTime::Date; idSamples.clear(); }
	// Functions
	int insert_preds(MedFeatures& featuresData);
	void get_ids(vector<int>& ids);
	void append(MedSamples& newSamples) { idSamples.insert(idSamples.end(), newSamples.idSamples.begin(), newSamples.idSamples.end()); }

	// bin file
	int read_from_bin_file(const string& file_name) { return SerializableObject::read_from_file(file_name); }
	int write_to_bin_file(const string& file_name) { return SerializableObject::write_to_file(file_name); }

	// text file - default option to read/write samples
	int read_from_file(const string& file_name);
	int write_to_file(const string &fname);

	void get_preds(vector<float>& preds);
	void get_y(vector<float>& y);
	void get_categs(vector<float> &categs); // gets a list of all categories appearing in the outcome

	void sort_by_id_date(); 

	ADD_SERIALIZATION_FUNCS(time_unit, idSamples)

};


//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(MedSample)
MEDSERIALIZE_SUPPORT(MedIdSamples)
MEDSERIALIZE_SUPPORT(MedSamples)

#endif
