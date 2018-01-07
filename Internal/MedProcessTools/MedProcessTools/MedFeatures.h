// MedFeatures - holding data as a map from name to vector

#ifndef _MED_FEAT_H_
#define _MED_FEAT_H_

#include "InfraMed/InfraMed/InfraMed.h"
#include "MedProcessTools/MedProcessTools/MedSamples.h"
#include "MedProcessTools/MedProcessTools/SerializableObject.h"
#include "MedProcessTools/MedProcessTools/MedProcessUtils.h"

//.......................................................................................
/**  A structure holding feature attributes
*/
//.......................................................................................
class FeatureAttr {
public:
	bool normalized; ///< indicator that the feature has been normalized
	bool imputed; ///< indicator that the feature has been imputed and does not contain missing values
};


//.......................................................................................
/** A class for holding features data as a virtual matrix <br>
*	- A vector of samples (id + date + outcome + ...) <br>
*	- A vector of weights (one per sample) <br>
*	- For each feature, a vector of floats (feature values - one per sample) <br>
* <br>
*	- Metadata per feature <br>
*		- a FeatureAttr entry <br>
*		- a set of tags (used by FeatureProcess to decide if to act on feature) <br>
* <br>
*	- General attribute: <br>
*		- time_unit : the time-unit in which samples are given <br>
* <br>
* - Helpers :  <br>
*		- pid_pos_len :  pid_pos_len[pid].first holds the first row in the matrix per id, 
*						 pid_pos_len[pid].second holds the number of relevant rows
*/
//.......................................................................................
class MedFeatures : public SerializableObject {
public:

	// Data
	map<string, vector<float> > data; ///< the actual matrix of values per sample
	vector<float> weights; ///< a vector of weight per sample
	vector<MedSample> samples; ///< The samples representing the lines

	/// feature generation assumes that all "rows" for a specific pid are adjacent.
	/// pid_pos_len[pid].first holds the first position, pid_pos_len[pid].second holds its length
	map<int, pair<int, int>> pid_pos_len;

	// Attributes
	map<string, FeatureAttr> attributes; ///< a FeatureAttr per feature
	map<string, unordered_set<string> > tags; ///< a set of tags per feature

	// time Unit
	int time_unit; ///< the time unit of the samples 

	/// A global counter used to prevent identical names for two features by adding FTR_#_ before generated feature name.
	static int global_serial_id_cnt;

	// Functions

	/// <summary> Constructor Given time-unit </summary>
	MedFeatures(int _time_unit) { time_unit = _time_unit; }
	///<summary>  Constructor setting time-unit according to med_rep_type </summary>
	MedFeatures() { time_unit = med_rep_type.basicTimeUnit; };

	// Initialization
	/// <summary> Clear all vectors </summary>
	void clear() { data.clear(); samples.clear(); pid_pos_len.clear(); attributes.clear(); weights.clear(); tags.clear(); }
	/// <summary> set time unit </summary>
	void set_time_unit(int _time_unit) { time_unit = _time_unit; }

	/// <summary> Get a vector of feature names </summary>
	void get_feature_names(vector<string>& names) ;
	/// <summary> Get data (+attributes) as a MedMat </summary>
	void get_as_matrix(MedMat<float>& mat);
	/// <summary> Get subset of data (+attributes) as a MedMat : Only features in 'names' </summary>
	void get_as_matrix(MedMat<float>& mat, vector<string>& names);
	/// <summary> Get subset of data (+attributes) as a MetMat: Only features in 'names' and rows in 'idx' </summary>
	void get_as_matrix(MedMat<float>& mat, const vector<string>& names, vector<int> &idx);

	/// <summary> Append samples at end of samples vector (used for generating samples set before generating features) </summary>
	void append_samples(MedIdSamples& in_samples);
	/// <summary> Insert samples at position idex, assuming samples vector is properly allocated (used for generating samples set before generating features) </summary>
	void insert_samples(MedIdSamples& in_samples, int index);
	/// <summary> Fill samples vetor and initialize pid_pos_len according to input vector of MedIdSamples </summary>
	void init_all_samples(vector<MedIdSamples> &in_samples) { samples.clear(); for (auto& id : in_samples) append_samples(id); init_pid_pos_len(); }
	/// <summary> initialize pid_pos_len vector according to samples </summary>
	void init_pid_pos_len();
	/// <summary> Return the first row in the virtual matrix for an id (-1 if none) </summary>
	int get_pid_pos(int pid) { if (pid_pos_len.find(pid) == pid_pos_len.end()) return -1; return (pid_pos_len[pid].first); }
	/// <summary> Return the number of rows in the virtual matrix for an id (-1 if none) </summary>
	int get_pid_len(int pid) { if (pid_pos_len.find(pid) == pid_pos_len.end()) return 0; return (pid_pos_len[pid].second); }
	/// <summary> Calculate a crc for the data (used for debugging mainly) </summary>
	unsigned int get_crc();
	/// <summary> MLOG data in csv format </summary>
	void print_csv();
	/// <summary> Get the corresponding MedSamples object .  Assuming samples vector in features are ordered  (all id's samples are consecutive) </summary>
	void get_samples(MedSamples& outSamples);
	/// <summary> Return the max serial_id_cnt </summary>
	int get_max_serial_id_cnt();

	/// <summary> Write features (samples + weights + data) as csv with a header line  </summary>
	/// <returns> -1 upon failure to open file, 0 upon success </returns>
	int write_as_csv_mat(const string &csv_fname);
	/// <summary> Read features (samples + weights + data) from a csv file with a header line </summary>
	/// <returns> -1 upon failure to open file, 0 upon success </returns>
	int read_from_csv_mat(const string &csv_fname);

	/// <summary> Filter data (and attributes) to include only selected features </summary> 
	/// <return> -1 if any of the selected features is not present. 0 upon success  </returns>
	int filter(unordered_set<string>& selectedFeatures);

	// Serialization
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
};

#endif
