#ifndef __MED__MPMODEL__H__
#define __MED__MPMODEL__H__

#include "MedPyCommon.h"
#include "MPPidRepository.h"
#include "MPFeatures.h"
#include "MPPredictor.h"

class MedModel;
class MPSamples;

class MPModelStage {

public:

	static const int LEARN_REP_PROCESSORS;
	static const int LEARN_FTR_GENERATORS;
	static const int APPLY_FTR_GENERATORS;
	static const int LEARN_FTR_PROCESSORS;
	static const int APPLY_FTR_PROCESSORS;
	static const int LEARN_PREDICTOR;
	static const int APPLY_PREDICTOR;
	static const int INSERT_PREDS;
	static const int LEARN_POST_PROCESSORS;
	static const int APPLY_POST_PROCESSORS;
	static const int END;
};


class MPModel {
public:
	MEDPY_IGNORE( MedModel* o );
	MPModel();
	~MPModel();
	void init_from_json_file(const std::string& fname);
	std::vector<std::string> init_from_json_file_with_alterations(const std::string& fname, std::vector<std::string> json_alt);
	void add_pre_processors_json_string_to_model(string in_json, string fname);
	std::vector<std::string> get_required_signal_names();
	int learn(MPPidRepository* rep, MPSamples* samples);
	int apply(MPPidRepository* rep, MPSamples* samples);
	int learn(MPPidRepository* rep, MPSamples* samples, int start_stage, int end_stage);
	int apply(MPPidRepository* rep, MPSamples* samples, int start_stage, int end_stage);

	int write_to_file(const std::string &fname);
	int read_from_file(const string &fname);
	MPFeatures MEDPY_GET_features();

	void clear();
	int MEDPY_GET_verbosity();
	void MEDPY_SET_verbosity(int new_vval);

	void add_feature_generators(string& name, vector<string>& signals);
	void add_feature_generators(string& name, vector<string>& signals, string init_string);
	void add_feature_generator(string& name, string& signal);
	void add_feature_generators(string& name, string& signal, string init_string);

	void add_age();
	void add_gender();

	void get_all_features_names(vector<string> &feat_names, int before_process_set);

	void add_normalizers();
	void add_normalizers(string init_string);
	void add_normalizers(vector<string>& features);
	void add_normalizers(vector<string>& features, string init_string);
	void add_imputers();
	void add_imputers(string init_string);
	void add_imputers(vector<string>& features);
	void add_imputers(vector<string>& features, string init_string);
	void add_rep_processor_to_set(int i_set, const string &init_string);
	void add_feature_generator_to_set(int i_set, const string &init_string);
	void add_feature_processor_to_set(int i_set, int duplicate, const string &init_string);
	void add_process_to_set(int i_set, int duplicate, const string &init_string);
	void add_process_to_set(int i_set, const string &init_string);
	void set_predictor(MPPredictor& _predictor);
	void make_predictor(string name);  //original name- set_predictor
	void set_predictor(string name, string init_string);
	int collect_and_add_virtual_signals(MPPidRepository &rep);
	int quick_learn_rep_processors(MPPidRepository& rep, MPSamples& samples);
	int learn_rep_processors(MPPidRepository& rep, MPSamples& samples);
	void filter_rep_processors();
	int learn_feature_generators(MPPidRepository &rep, MPSamples *learn_samples);
	int generate_all_features(MPPidRepository &rep, MPSamples *samples, MPFeatures &features, std::vector<std::string> req_feature_generators);
	int learn_and_apply_feature_processors(MPFeatures &features);
	int learn_feature_processors(MPFeatures &features);
	int apply_feature_processors(MPFeatures &features, bool learning);
	void fit_for_repository(MPPidRepository &rep);
	void calc_contribs(MPMat &mat, MPMat &mat_out);
	void calc_feature_contribs_conditional(MPMat &mat_x_in, const string& features_cond_string, float features_cond_float, MPMat &mat_x_out, MPMat &mat_contribs);
	/// following is for debugging, it gets a prefix, and prints it along with information on rep_processors, feature_generators, or feature_processors
	void dprint_process(const string &pref, int rp_flag, int fg_flag, int fp_flag, bool pp_flag, bool predictor_type);

	/// following is for debugging : writing the feature to a csv file as a matrix.
	int write_feature_matrix(const string mat_fname);
	MPSerializableObject asSerializable();
};



#endif // !__MED__MPMODEL__H__
