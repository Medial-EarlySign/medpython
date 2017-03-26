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
	int id; // Optional
	int time;
	float outcome;
	vector<float> prediction;
	int outcomeTime;

	MedSample() { prediction.clear(); }
	~MedSample() { prediction.clear(); }

	// (De)Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// print
	void print(const string prefix);
	void print() { print(""); }
};

// A collection of samples of a given id
class MedIdSamples {
public:
	int id;
	int split;
	vector<MedSample> samples;

	MedIdSamples() { id = -1; split = -1; samples.clear(); }

	// De(Serialize)
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

};

// A collection of ids and relevant samples
class MedSamples {
public:
	int time_unit; // the time unit in which the samples are given. Default: Days
	vector<MedIdSamples> idSamples;


	// Constructor
	MedSamples() { time_unit = med_rep_type.basicTimeUnit; }
	// Functions
	int insert_preds(MedFeatures& featuresData);
	void get_ids(vector<int>& ids);
	void append(MedSamples& newSamples) { idSamples.insert(idSamples.end(), newSamples.idSamples.begin(), newSamples.idSamples.end()); }
	int read_from_file(const string& file_name);
	int write_to_file(const string &fname);
	void get_preds(vector<float>& preds);
	void get_y(vector<float>& y);

	// De(Serialize)
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

};


//=======================================
// Joining the MedSerialze wagon
//=======================================
MEDSERIALIZE_SUPPORT(MedSample)
MEDSERIALIZE_SUPPORT(MedIdSamples)
MEDSERIALIZE_SUPPORT(MedSamples)

#endif
