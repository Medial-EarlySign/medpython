// Temporary class for holding samples

#ifndef _MED_SAMPLES_H_
#define _MED_SAMPLES_H_

#include "MedAlgo/MedAlgo/MedAlgo.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
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
	int date;
	float outcome;
	vector<float> prediction;
	int outcomeDate;

	// (De)Serialization
	size_t get_size() { return sizeof(MedSample); }
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

};

// A collection of samples of a given id
class MedIdSamples {
public:
	int id;
	vector<MedSample> samples;

	// De(Serialize)
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

};

// A collection of ids and relevant samples
class MedSamples {
public:
	vector<MedIdSamples> idSamples;

	// Functions
	int insert_preds(MedFeatures& featuresData);
	void get_ids(vector<int>& ids);
	void append(MedSamples& newSamples) { idSamples.insert(idSamples.end(), newSamples.idSamples.begin(), newSamples.idSamples.end()); }
	int read_from_file(const string& file_name);

	// De(Serialize)
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

};

#endif
