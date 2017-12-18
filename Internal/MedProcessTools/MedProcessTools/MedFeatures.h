// MedFeatures - holding data as a map from name to vector

#ifndef _MED_FEAT_H_
#define _MED_FEAT_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"

//.......................................................................................
//.......................................................................................
// A structure holding feature attributes
//.......................................................................................
//.......................................................................................
class FeatureAttr {
public:
	bool normalized;
	bool imputed;
};


//.......................................................................................
//.......................................................................................
// A class for holding features data as a virtual matrix
//	- A vector of samples (id + date + outcome + ...)
//  - A vector of weights (one per sample)
//	- For each feature, a vector of floats (feature values - one per sample)
//
//	- Metadata per feature
//		- a FeatureAttr entry
//		- a set of tags (used by FeatureProcess to decide if to act on feature)
//
//	- General attribute:
//		- time_unit : the time-unit in which samples are given
//
// - Helpers : 
//		- pid_pos_len :  pid_pos_len[pid].first holds the first row in the matrix per id, 
//						 pid_pos_len[pid].second holds the number of relevant rows
//.......................................................................................
//.......................................................................................
class MedFeatures : public SerializableObject {
public:

	// Data
	map<string, vector<float> > data;
	vector<float> weights;
	vector<MedSample> samples;

	// feature generation assumes that all "rows" for a specific pid are adjacent.
	// pid_pos_len[pid].first holds the first position, pid_pos_len[pid].second holds its length
	map<int, pair<int, int>> pid_pos_len;

	// Attributes
	map<string, FeatureAttr> attributes;
	map<string, unordered_set<string> > tags;

	// time Unit
	int time_unit;

	// A global counter used to prevent identical names for two features by adding FTR_#_ before generated feature name.
	static int global_serial_id_cnt;

	// Functions

	// Constructor
	// Given time-unit
	MedFeatures(int _time_unit) { time_unit = _time_unit; }
	// Set time-unit according to med_rep_type
	MedFeatures() { time_unit = med_rep_type.basicTimeUnit; };
	~MedFeatures() {};

	// Initialization
	void clear() { data.clear(); samples.clear(); pid_pos_len.clear(); attributes.clear(); weights.clear(); tags.clear(); }
	void set_time_unit(int _time_unit) { time_unit = _time_unit; }

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Get a vector of feature names
	void get_feature_names(vector<string>& names) ;
	// Get data (+attributes) as matrix
	void get_as_matrix(MedMat<float>& mat);
	// Get subset of data (+attributes) as matrix : Only features in 'names'
	void get_as_matrix(MedMat<float>& mat, vector<string>& names);
	// Get subset of data (+attributes) as matrix: Only features in 'names' and rows in 'idx'
	void get_as_matrix(MedMat<float>& mat, const vector<string>& names, vector<int> &idx);


	// Append samples at end of samples vector
	void append_samples(MedIdSamples& in_samples);
	// Insert samples at position idex, assuming samples vector is properly allocated
	void insert_samples(MedIdSamples& in_samples, int index);
	// Fill samples vetor and initialize pid_pos_len according to input vector of MedIdSamples
	void init_all_samples(vector<MedIdSamples> &in_samples) { samples.clear(); for (auto& id : in_samples) append_samples(id); init_pid_pos_len(); }
	// initialize pid_pos_len vector according to samples
	void init_pid_pos_len();
	// Get first row in the virtual matrix for an id (-1 if none)
	int get_pid_pos(int pid) { if (pid_pos_len.find(pid) == pid_pos_len.end()) return -1; return (pid_pos_len[pid].first); }
	// Get the number of rows in the virtual matrix for an id (-1 if none)
	int get_pid_len(int pid) { if (pid_pos_len.find(pid) == pid_pos_len.end()) return 0; return (pid_pos_len[pid].second); }
	// Calculate a crc for the data (used for debugging mainly)
	unsigned int get_crc();
	// MLOG data in csv format
	void print_csv();

	// Write features (samples + weights + data) as csv with a header line
	int write_as_csv_mat(const string &csv_fname);
	// Read features (samples + weights + data) from a csv file with a header line
	int read_from_csv_mat(const string &csv_fname);

	// Filter data (and attributes) to include only selected features
	int filter(unordered_set<string>& selectedFeatures);
};

#endif
