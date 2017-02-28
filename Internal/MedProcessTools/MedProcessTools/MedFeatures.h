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
};


//.......................................................................................
//.......................................................................................
// A class for holding features data
//.......................................................................................
//.......................................................................................
class MedFeatures : public SerializableObject {
public:

	// Map
	map<string, vector<float> > data;
	vector<MedSample> samples;

	// feature generation assumes that all "rows" for a specific pid are adjacent.
	// pid_pos_len[pid].first holds the first position, pid_pos_len[pid].second holds its length
	map<int, pair<int, int>> pid_pos_len;

	// Attributes
	map<string, FeatureAttr> attributes;

	// time Unit
	int time_unit;

	// Constructor/Destructor
	MedFeatures(int _time_unit) { time_unit = _time_unit; }
	MedFeatures() {};
	~MedFeatures() {};

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	// Functions
	void get_feature_names(vector<string>& names) ;
	void get_as_matrix(MedMat<float>& mat);
	
	void append_samples(MedIdSamples& in_samples);
	void insert_samples(MedIdSamples& in_samples, int index);
	void init_all_samples(vector<MedIdSamples> &in_samples) { samples.clear(); for (auto& id : in_samples) append_samples(id); init_pid_pos_len(); }
	void init_pid_pos_len();
	int get_pid_pos(int pid) { if (pid_pos_len.find(pid) == pid_pos_len.end()) return -1; return (pid_pos_len[pid].first); }
	int get_pid_len(int pid) { if (pid_pos_len.find(pid) == pid_pos_len.end()) return 0; return (pid_pos_len[pid].second); }

	unsigned int get_crc(); // used for debugging mainly
	void print_csv();
};

#endif
