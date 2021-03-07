// MedFeatures - holding data as a map from name to vector

#ifndef __MED_FEATURES_H__
#define __MED_FEATURES_H__

#include <InfraMed/InfraMed/InfraMed.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>
#include <MedProcessTools/MedProcessTools/MedProcessUtils.h>
#include <MedProcessTools/MedProcessTools/MedSamples.h>
#include <MedMat/MedMat/MedMat.h>
#include <random>

//.......................................................................................
//  A structure holding feature attributes
//.......................................................................................
class FeatureAttr : public SerializableObject {
public:
	bool normalized = false; ///< indicator that the feature has been normalized
	bool imputed = false; ///< indicator that the feature has been imputed and does not contain missing values
	float denorm_mean = 0.0, denorm_sdv = 1.0; ///< Mean and Standard deviation for de-normalization

	unordered_map<float, string> value2Name; ///< map value to name (e.g. for naming one-hot feature processor)

	// Serialization
	ADD_CLASS_NAME(FeatureAttr)
	ADD_SERIALIZATION_FUNCS(normalized, imputed , denorm_mean, denorm_sdv, value2Name)
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
class MedFeatures final : public SerializableObject {
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

	// masks for cleaning , imputing
	const static unsigned char cleaned_mask = (unsigned char)0x01;
	const static unsigned char imputed_mask = (unsigned char)0x02;
	map<string, vector<unsigned char>> masks;
	float medf_missing_value = (float)MED_MAT_MISSING_VALUE;

	// time Unit
	int time_unit; ///< the time unit of the samples 

	/// A global counter used to prevent identical names for two features by adding FTR_#_ before generated feature name.
	static int global_serial_id_cnt;

	// Functions

	/// <summary> Constructor Given time-unit </summary>
	MedFeatures(int _time_unit) { time_unit = _time_unit; }
	///<summary>  Constructor setting time-unit to undef </summary>
	MedFeatures() { time_unit = global_default_time_unit; global_serial_id_cnt = 0; }

	// Initialization
	/// <summary> Clear all vectors </summary>
	void clear() { data.clear(); samples.clear(); pid_pos_len.clear(); attributes.clear(); weights.clear(); tags.clear(); masks.clear(); }
	/// <summary> set time unit </summary>
	void set_time_unit(int _time_unit) { time_unit = _time_unit; }

	/// <summary> Get a vector of feature names </summary>
	void get_feature_names(vector<string>& names) const;
	/// <summary> Get data (+attributes) as a MedMat </summary>
	void get_as_matrix(MedMat<float>& mat) const;
	/// <summary> Get subset of data (+attributes) as a MedMat : Only features in 'names' </summary>
	void get_as_matrix(MedMat<float>& mat, vector<string>& names) const;
	/// <summary> Get subset of data (+attributes) as a MetMat: Only features in 'names' and rows in 'idx' </summary>
	void get_as_matrix(MedMat<float>& mat, const vector<string>& names, vector<int> &idx) const;

	/// <summary> Set data (+attributes) from MedMat </summary>
	void set_as_matrix(const MedMat<float>& mat);

	/// <summary> Append samples at end of samples vector (used for generating samples set before generating features) </summary>
	void append_samples(MedIdSamples& in_samples);
	/// <summary> Insert samples at position idex, assuming samples vector is properly allocated (used for generating samples set before generating features) </summary>
	void insert_samples(MedIdSamples& in_samples, int index);
	/// <summary> Fill samples vetor and initialize pid_pos_len according to input vector of MedIdSamples </summary>
	void init_all_samples(vector<MedIdSamples> &in_samples) { samples.clear(); for (auto& id : in_samples) append_samples(id); init_pid_pos_len(); }
	/// <summary> initialize pid_pos_len vector according to samples </summary>
	void init_pid_pos_len();
	/// <summary> Return the first row in the virtual matrix for an id (-1 if none) </summary>
	int get_pid_pos(int pid) const { if (pid_pos_len.find(pid) == pid_pos_len.end()) return -1; return (pid_pos_len.at(pid).first); }
	/// <summary> Return the number of rows in the virtual matrix for an id (-1 if none) </summary>
	int get_pid_len(int pid) const { if (pid_pos_len.find(pid) == pid_pos_len.end()) return 0; return (pid_pos_len.at(pid).second); }
	/// <summary> Calculate a crc for the data (used for debugging mainly) </summary>
	unsigned int get_crc();
	/// <summary> MLOG data in csv format </summary>
	void print_csv() const;
	/// <summary> Get the corresponding MedSamples object .  Assuming samples vector in features are ordered  (all id's samples are consecutive) </summary>
	void get_samples(MedSamples& outSamples) const;
	/// <summary> Return the max serial_id_cnt </summary>
	int get_max_serial_id_cnt() const;

	/// <summary> Write features (samples + weights + data) as csv with a header line  </summary>
	/// <returns> -1 upon failure to open file or attributes inconsistency (if write_attributes is true), 0 upon success </returns>
	int write_as_csv_mat(const string &csv_fname, bool write_attributes = false) const;
	int add_to_csv_mat(const string &csv_fname, bool write_attributes, int start_idx) const;
	void write_csv_data(ofstream& out_f, bool write_attributes, vector<string>& col_names, int start_idx) const;

	/// <summary> Read features (samples + weights + data) from a csv file with a header line </summary>
	/// <returns> -1 upon failure to open file, 0 upon success </returns>
	int read_from_csv_mat(const string &csv_fname, bool read_time_raw = true);

	/// <summary> Filter data (and attributes) to include only selected features </summary> 
	/// <return> -1 if any of the selected features is not present. 0 upon success  </returns>
	int filter(unordered_set<string>& selectedFeatures);

	/// preparing a list all features that contain as a substring one of the given search strings, adds (that is not clearing selected on start)
	int prep_selected_list(vector<string>& search_str, unordered_set<string> &selected);

	// masks functions
	int init_masks();
	int get_masks_as_mat(MedMat<unsigned char> &masks_mat);
	int mark_imputed_in_masks(float _missing_val);
	int mark_imputed_in_masks() { return mark_imputed_in_masks(medf_missing_value); }

	void round_data(float r);
	void noise_data(float r);

	///\brief Sort by id and time
	void samples_sort();

	///\brief Get feature name that matches a substring
	string resolve_name(string& substr) const;

	// Serialization
	ADD_CLASS_NAME(MedFeatures)
	ADD_SERIALIZATION_FUNCS(data, weights, samples, attributes, tags, time_unit)

};

MEDSERIALIZE_SUPPORT(MedFeatures)
MEDSERIALIZE_SUPPORT(FeatureAttr)

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	*  \brief process namespace
	*/
	namespace process {
		/// \brief commit selection of indexes on vector
		template<class T> void commit_selection(vector<T> &vec, const vector<int> &idx);
		/// \brief filtering MedFeatures by selected indexes rows
		void filter_row_indexes(MedFeatures &dataMat, vector<int> &selected_indexes, bool op_flag = false);
		/// \brief filtering MedFeatures by selected indexes rows (thread safe for selected_indexes) no sort of selected_indexes
		void filter_row_indexes_safe(MedFeatures &dataMat, const vector<int> &selected_indexes, bool op_flag = false);
		/// \brief down sampling with ratio
		void down_sample(MedFeatures &dataMat, double take_ratio, bool with_repeats = false,
			vector<int> *selected_indexes = NULL);
		/// \brief reweighting method by given groups uniq values. return weights and min_factor
		double reweight_by_general(MedFeatures &data_records, const vector<string> &groups,
			vector<float> &weigths, bool print_verbose);
		/// \brief matching by given groups uniq values. returns also the row_ids filtered
		void match_by_general(MedFeatures &data_records, const vector<string> &groups,
			vector<int> &filtered_row_ids, float price_ratio, int min_grp_size, bool print_verbose);
		/// \brief matching by given groups uniq values. returns also the row_ids filtered. max_ratio is maximal allowed ratio, inf if < 0.
		void match_by_general(MedFeatures &data_records, const vector<string> &groups,
			vector<int> &filtered_row_ids, float price_ratio, float max_ratio, int min_grp_size, bool print_verbose);

		/// \brief split matrix to train test based on iFold value. folds is fold id for each sample
		void split_matrix(const MedFeatures& matrix, vector<int>& folds, int iFold,
			MedFeatures& trainMatrix, MedFeatures& testMatrix, const vector<string> *selected_features = NULL);
		/// \brief split matrix to train test based on iFold value. folds is map from patient id to fold
		void split_matrix(const MedFeatures& matrix, unordered_map<int, int>& folds, int iFold,
			MedFeatures& trainMatrix, MedFeatures& testMatrix, const vector<string> *selected_features = NULL);
		/// \brief convert feature vector to it's prctil's value in each element
		void convert_prctile(vector<float> &features_prctiles);
		/// \brief does matching to specific target_prior. 
		/// @param outcome is the outcome vector for measure prior in each group
		/// @param group_values is the groups to split the matching to. it can be year signature or age or unique combination of both
		/// @param target_prior the target prior
		/// @param the return value of selected indexes to do the matching
		void match_to_prior(const vector<float> &outcome,
			const vector<float> &group_values, float target_prior, vector<int> &sel_idx);

		/// \brief does matching to specific prior for MedSamples
		double match_to_prior(MedSamples &samples, float target_prior, vector<int> &sel_idx);

		/// \brief does matching to specific prior for MedFeatures
		double match_to_prior(MedFeatures &features, float target_prior, vector<int> &sel_idx);

		/// \brief does matching to specific target_prior. 
		/// @param features the matrix to match. will use outcome in samples
		/// @param group_values is the groups to split the matching to. it can be year signature or age or unique combination of both
		/// @param target_prior the target prior
		/// @param sel_idx the original indecies 
		void match_to_prior(MedFeatures &features,
			const vector<string> &group_values, float target_prior, vector<int> &sel_idx, bool print_verbose = true);

		/// \brief Return number of splits, also check mismatches between idSample and internal MedSamples and set idSamples.split if missing
		int nSplits(vector<MedSample>& samples);

		/// \brief multi-class matching.
		float match_multi_class(MedFeatures& data, const vector<string> &groups, vector<int> &filtered_row_ids, vector<float>& price_ratios, int nRand = 10000, int verbose = false);
		float match_multi_class(vector<MedSample>& data, const vector<string> &groups, vector<int> &filtered_row_ids, vector<float>& price_ratios, int nRand = 10000, int verbose = false);

		void match_multi_class_to_dist(MedFeatures& data, const vector<string> &groups, vector<int> &filtered_row_ids, vector<float> probs);
		void match_multi_class_to_dist(vector<MedSample>& data, const vector<string> &groups, vector<int> &filtered_row_ids, vector<float> probs);

	}

}

#endif
