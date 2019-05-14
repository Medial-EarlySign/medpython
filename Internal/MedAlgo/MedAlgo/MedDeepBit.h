#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <algorithm>
#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//Deep Bit

struct MedDeepBitParams : public SerializableObject {
	int max_depth;
	int num_iterations;
	int num_ftrs_per_round;
	int num_vals_to_be_categorial;
	int nparts_auc;
	int niter_auc_gitter;
	int niter_coordinate_descent;
	int internal_test_ratio;
	double fraction_auc;
	double grid_fraction;
	double min_fraction_zeros_ones;
	double frac_continuous_frequent;
	double frac_categorial_frequent;
	double lambda;
	double min_cor_bin_ftr;

	MedDeepBitParams() {
		max_depth = 20; num_iterations = 100; num_ftrs_per_round = 1000; num_vals_to_be_categorial = 10; nparts_auc = 1; niter_coordinate_descent = 0;
		niter_auc_gitter = 0;  grid_fraction = 0.003; min_fraction_zeros_ones = 0.003;
		frac_continuous_frequent = 0.1; frac_categorial_frequent = 0.1; lambda = 20; fraction_auc = 1; internal_test_ratio = 10, min_cor_bin_ftr = 0.05;
	}
	ADD_CLASS_NAME(MedDeepBitParams)
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
	string to_string();
};

MEDSERIALIZE_SUPPORT(MedDeepBitParams)

class MedDeepBit : public MedPredictor {	

public:
	MedDeepBitParams params;
	void init_defaults();
	int init(void *classifier_params) { this->params = *((MedDeepBitParams*)classifier_params); return 0; };
	/// The parsed fields from init command.
	/// @snippet MedDeepBit.cpp MedDeepBit::init
	int set_params(map<string, string>& initialization_map);
	MedDeepBit() { classifier_type = MODEL_DEEP_BIT; init_defaults(); }
	~MedDeepBit() {}
	MedDeepBit(void *_params) { params = *(MedDeepBitParams *)_params; }
	MedDeepBit(MedDeepBitParams& _params) { params = _params; }
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Learn(float *x, float *y, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;
	virtual void print(FILE *fp, const string& prefix, int level=0) const;
	ADD_CLASS_NAME(MedDeepBit)
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);
	void print_model(FILE *fp, const string& prefix) const;

private:
	vector<vector<double>> x, internal_test_x, internal_test_transposed_x;
	vector<double> y, scores, r, internal_test_scores, avx, sdx;
	vector<int> label, internal_test_label;
	int nftrs, nsamples, internal_test_nsamples, num_bin_ftrs;
	vector<string> ftr_names;
	double avy;

	vector<vector<double>> ftr_grids;
	vector<int> is_categorial;
	vector<vector<double>> frequent_ftr_vals;
	// For each binary feature we document the raw feature number, a spceific value val of that feature and a boolean value bool_val. 
	// For categorial featues, we ask weather a given sample equals val or not, for bool_val in {true, false}, correpondingly.
	// For quantitative features, we ask weather a given sample is above or under val for the values of true and false for bool_val, correspondingly.
	vector<tuple<int, int, bool, bool>> bin_ftrs_map;
	vector<vector<char>> bin_ftrs;
	vector<vector<int>> bin_ftr_indexes;
	vector<vector<double>> bin_ftr_avg_sd_beta;

	void print_log();
	void train_model();
	void predict_train() {}
	double predict_sample(const vector<double>& x, int niter) const;
	double predict_sample(const vector<double>& x) const;
	void get_normalized_col(const vector<char>& col, vector<double>& normalized_col, double& av, double& std);
	double perform_lasso_iteration(const vector<double>& xk_train, const vector<double>& r, double lambda, double alpha);
	void get_col_without_na(const vector<double>& col, vector<double>& col_without_na);
	void get_avgs(const vector<vector<double>>& x, vector<double>& avx);
	void get_sds(const vector<vector<double>>& x, vector<double>& sdx);
	double avg(const vector<char>& vec);
	double avg(const vector<double>& vec);
	double sd(const vector<char>& binary_vec);
	double sd(const vector<double>& vec);
	void transpose(const vector<vector<double>>& before, vector<vector<double>>& after);
	void init(float *x1, float *y1, int nsamples1, int nftrs1);
	void predict(const vector<vector<double>>& x, vector<double>& scores) const;
	double get_normalized_val(double x_val, int j);
	void impute_x(vector<vector<double>>& x, const vector<double>& avx);
	bool is_viable_01_ratios(int count0, int count1, int count_pos0, int count_pos1);
	bool is_bin_ftr_valid(const vector<char>& bin_ftr);
	void make_bin_ftrs(int j, const vector<double>& vals, bool is_categorial);
	int step_function(int step, int i);
	void mark_grids_and_frequent_vals();
	void get_categorial_bin_ftr(const vector<double>& col, double val, bool direction, vector<char>& bin_ftr);
	void get_quant_bin_ftr(const vector<double>& col, double val, bool direction, vector<char>& bin_ftr);
	int get_categorial_bit(double x_val, double val, bool direction) const;
	int get_quant_bit(double x_val, double val, bool direction) const;
	double get_ftr_score(const vector<char>& bin_ftr);
	void mult_bin_ftrs(const vector<char>& ftr1, vector<char>& ftr2);
	void get_bin_ftr(int bin_ftr_index, vector<char>& bin_ftr);
	void print_ftr_characteristics(int index);
	void calc_bin_ftr_scores(const vector<char>& bin_ftr, double& av, double& std, double& b, vector<double>& scores1, vector<double>& r1, bool is_full_step = false);
	void gen_random_indexes(vector<int>& random_indexes);
	void score_random_ftrs(vector<double>& ftr_scores, const vector<int>& random_indexes, const vector<char>& final_bin_ftr);
	void get_bin_ftr_of_it(int it, vector<char>& bin_ftr);
	void do_auc_gittering();
	void do_coordinate_descent(int num_iterations_descent);
	double special_auc(const vector<double>& all_predictions, const vector<int>& all_label, bool is_weighted = false, int nparts = 1);
};

MEDSERIALIZE_SUPPORT(MedDeepBit)














