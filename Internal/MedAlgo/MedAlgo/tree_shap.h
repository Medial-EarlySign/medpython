/**
* Fast recursive computation of SHAP values in trees.
* See https://arxiv.org/abs/1802.03888 for details.
*
* Scott Lundberg, 2018 (independent algorithm courtesy of Hugh Chen 2018)
*/
#ifndef __TREE_SHAP_H__
#define __TREE_SHAP_H__

#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <cmath>
#include <ctime>
#if defined(_WIN32) || defined(WIN32)
#include <malloc.h>
#else
#include <alloca.h>
#endif

#include <vector>
#include <random>
#include <unordered_set>
#include <MedProcessTools/MedProcessTools/MedFeatures.h>
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include "SamplesGenerator.h"

using namespace std;

typedef double tfloat;

#define FROM_NEITHER 0
#define FROM_X_NOT_R 1
#define FROM_R_NOT_X 2

namespace FEATURE_DEPENDENCE {
	const unsigned independent = 0;
	const unsigned tree_path_dependent = 1;
	const unsigned global_path_dependent = 2;
}

namespace MODEL_TRANSFORM {
	const unsigned identity = 0;
	const unsigned logistic = 1;
	const unsigned logistic_nlogloss = 2;
	const unsigned squared_loss = 3;
}

struct ExplanationDataset {
	tfloat *X; ///< vector of all data. each row is sample of all features for that sample. cols(2nd dim) are features
	bool *X_missing; ///< bool mask to return true on missing value on matrix - same structure as X
	tfloat *y; ///< the labels
	tfloat *R;
	bool *R_missing;
	unsigned num_X; ///< number of samples
	unsigned M; ///< Features count
	unsigned num_Exp; /// number of explanation features (allowing for grouping)
	unsigned num_R;

	ExplanationDataset();
	ExplanationDataset(tfloat *X, bool *X_missing, tfloat *y, tfloat *R, bool *R_missing, unsigned num_X,
		unsigned M, unsigned num_R, unsigned num_Exp);
	ExplanationDataset(tfloat *X, bool *X_missing, tfloat *y, tfloat *R, bool *R_missing, unsigned num_X,
		unsigned M, unsigned num_R);

	void get_x_instance(ExplanationDataset &instance, const unsigned i) const;
};

struct TreeEnsemble {
	int *children_left;
	int *children_right;
	int *children_default;
	int *features;
	tfloat *thresholds;
	tfloat *values;
	tfloat *node_sample_weights;
	unsigned max_depth;
	unsigned tree_limit;
	tfloat base_offset;
	unsigned max_nodes;
	unsigned num_outputs;
	bool is_allocate;

	TreeEnsemble();
	TreeEnsemble(int *children_left, int *children_right, int *children_default, int *features,
		tfloat *thresholds, tfloat *values, tfloat *node_sample_weights,
		unsigned max_depth, unsigned tree_limit, tfloat base_offset,
		unsigned max_nodes, unsigned num_outputs);

	void get_tree(TreeEnsemble &tree, const unsigned i) const;

	void allocate(unsigned tree_limit_in, unsigned max_nodes_in, unsigned num_outputs_in);

	void free();

	void fill_adjusted_tree(int node_index, ExplanationDataset& instance, const int *mask, unsigned *feature_sets, TreeEnsemble& adjusted);
	void create_adjusted_tree(ExplanationDataset& instance, const int *mask, unsigned *feature_sets, TreeEnsemble& adjusted);
	void calc_feature_contribs_conditional(MedMat<float> &mat_x_in, unordered_map<string, float> contiditional_variables, MedMat<float> &mat_x_out, MedMat<float> &mat_contribs);
	tfloat predict(ExplanationDataset& instance, int node_index);
};

inline void tree_shap(const TreeEnsemble& tree, const ExplanationDataset &data, tfloat *out_contribs, int condition, unsigned condition_feature, unsigned *feature_sets);

// data we keep about our decision path
// note that pweight is included for convenience and is not tied with the other attributes
// the pweight of the i'th path element is the permuation weight of paths with i-1 ones in them
struct PathElement {
	int feature_index;
	tfloat zero_fraction;
	tfloat one_fraction;
	tfloat pweight;
	PathElement();
	PathElement(int i, tfloat z, tfloat o, tfloat w);
};

// Independent Tree SHAP functions below here
// ------------------------------------------
struct Node {
	short cl, cr, cd, pnode, feat, pfeat; // uint_16
	float thres, value;
	char from_flag;
};

/**
* This runs Tree SHAP with a per tree path conditional dependence assumption.
*/
void dense_tree_saabas(tfloat *out_contribs, const TreeEnsemble& trees, const ExplanationDataset &data);

/**
* Runs Tree SHAP with feature independence assumptions on dense data.
*/
void dense_independent(const TreeEnsemble& trees, const ExplanationDataset &data,
	tfloat *out_contribs, tfloat transform(const tfloat, const tfloat));


/**
* This runs Tree SHAP with a per tree path conditional dependence assumption.
*/
void dense_tree_path_dependent(const TreeEnsemble& trees, const ExplanationDataset &data,
	tfloat *out_contribs, unsigned *feature_sets, tfloat transform(const tfloat, const tfloat));

// phi = np.zeros((self._current_X.shape[1] + 1, self._current_X.shape[1] + 1, self.n_outputs))
//         phi_diag = np.zeros((self._current_X.shape[1] + 1, self.n_outputs))
//         for t in range(self.tree_limit):
//             self.tree_shap(self.trees[t], self._current_X[i,:], self._current_x_missing, phi_diag)
//             for j in self.trees[t].unique_features:
//                 phi_on = np.zeros((self._current_X.shape[1] + 1, self.n_outputs))
//                 phi_off = np.zeros((self._current_X.shape[1] + 1, self.n_outputs))
//                 self.tree_shap(self.trees[t], self._current_X[i,:], self._current_x_missing, phi_on, 1, j)
//                 self.tree_shap(self.trees[t], self._current_X[i,:], self._current_x_missing, phi_off, -1, j)
//                 phi[j] += np.true_divide(np.subtract(phi_on,phi_off),2.0)
//                 phi_diag[j] -= np.sum(np.true_divide(np.subtract(phi_on,phi_off),2.0))
//         for j in range(self._current_X.shape[1]+1):
//             phi[j][j] = phi_diag[j]
//         phi /= self.tree_limit
//         return phi

void dense_tree_interactions_path_dependent(const TreeEnsemble& trees, const ExplanationDataset &data,
	tfloat *out_contribs,
	tfloat transform(const tfloat, const tfloat));

/**
* This runs Tree SHAP with a global path conditional dependence assumption.
*
* By first merging all the trees in a tree ensemble into an equivalent single tree
* this method allows arbitrary marginal transformations and also ensures that all the
* evaluations of the model are consistent with some training data point.
*/
void dense_global_path_dependent(const TreeEnsemble& trees, const ExplanationDataset &data,
	tfloat *out_contribs, tfloat transform(const tfloat, const tfloat));


/**
* The main method for computing Tree SHAP on model using dense data.
*/
void dense_tree_shap(const TreeEnsemble& trees, const ExplanationDataset &data, tfloat *out_contribs,
	const int feature_dependence, unsigned model_transform, bool interactions);
void dense_tree_shap(const TreeEnsemble& trees, const ExplanationDataset &data, tfloat *out_contribs,
	const int feature_dependence, unsigned model_transform, bool interactions, unsigned *feature_sets);

/**
* Iterative calling to Shapley
*/
void iterative_tree_shap(const TreeEnsemble& trees, const ExplanationDataset &data, tfloat *out_contribs,
	const int feature_dependence, unsigned model_transform, bool interactions, unsigned *feature_sets, bool verbose,
	vector<string>& names, const MedMat<float>& abs_cov_mat, int iteration_cnt, bool max_in_groups);

namespace medial {
	namespace shapley {

		/// \brief nchoosek calc
		double nchoosek(long n, long k);
		/// \brief lists all binary options for mask
		void list_all_options_binary(int nfeats, vector<vector<bool>> &all_opts);
		/// \brief generate random mask
		void generate_mask(vector<bool> &mask, int nfeat, mt19937 &gen, bool uniform_rand = false, bool use_shuffle = true);
		/// \brief generate random mask using already existed mask
		void generate_mask_(vector<bool> &mask, int nfeat, mt19937 &gen, bool uniform_rand = false, float uniform_rand_p = 0.5,
			bool use_shuffle = true, int limit_zero_cnt = 0);
		/// \brief generate all masks for shapley
		void sample_options_SHAP(int nfeats, vector<vector<bool>> &all_opts, int opt_count, mt19937 &gen, bool with_repeats
			, bool uniform_rand = false, bool use_shuffle = true);
		/// \brief shapley coeff calc
		double get_c(int p1, int p2, int end_l);
		/// \brief Shapley calculation without generator
		void explain_shapley(const MedFeatures &matrix, int selected_sample, int max_tests,
			MedPredictor *predictor, float missing_value, const vector<vector<int>>& group2index, const vector<string> &groupNames,
			vector<float> &features_coeff, mt19937 &gen, bool sample_masks_with_repeats,
			float select_from_all, bool uniform_rand, bool use_shuffle, bool verbose);
		/// \brief Shapley calculation with generator
		template<typename T> void explain_shapley(const MedFeatures &matrix, int selected_sample, int max_tests,
			MedPredictor *predictor, const vector<vector<int>>& group2index, const vector<string> &groupNames,
			const SamplesGenerator<T> &sampler_gen, mt19937 &rnd_gen, int sample_per_row, void *sampling_params,
			vector<float> &features_coeff, bool use_random_sample, bool verbose = false);

		/// \brief calculates minimal set
		void explain_minimal_set(const MedFeatures &matrix, int selected_sample, int max_tests,
			MedPredictor *predictor, float missing_value, const vector<vector<int>>& group2index
			, vector<float> &features_coeff, vector<float> &scores_history, int max_set_size,
			float baseline_score, float param_all_alpha, float param_all_beta,
			float param_all_k1, float param_all_k2, bool verbose);

		/// \brief calculates minimal set using sample generator
		void explain_minimal_set(const MedFeatures &matrix, int selected_sample, int max_tests,
			MedPredictor *predictor, float missing_value, const vector<vector<int>>& group2index,
			const SamplesGenerator<float> &sampler_gen, mt19937 &rnd_gen, void *sampling_params
			, vector<float> &features_coeff, vector<float> &scores_history, int max_set_size,
			float baseline_score, float param_all_alpha, float param_all_beta,
			float param_all_k1, float param_all_k2, bool verbose);

		///< sample weights = lime (distance from orig), uniform (1), shap (shapely weights) or sum (ensuring sum of weights per # of 1's ~ 1/(k*(n-k)) 
		typedef enum {
			LimeWeightLime = 0,
			LimeWeightUniform = 1,
			LimeWeightShap = 2,
			LimeWeightSum = 3,
			LimeWeightLast
		} LimeWeightMethod;

		/// \brief Shapley Lime with generator
		void get_shapley_lime_params(const MedFeatures& data, const MedPredictor *model,
			SamplesGenerator<float> *generator, float p, int n, LimeWeightMethod weighting, float missing,
			void *params, const vector<vector<int>>& group2index, const vector<string>& group_names, vector<vector<float>>& alphas);

		/// \brief Shapley Lime with generator and forced groups
		void get_shapley_lime_params(const MedFeatures& data, const MedPredictor *model,
			SamplesGenerator<float> *generator, float p, int n, LimeWeightMethod weighting, float missing,
			void *params, const vector<vector<int>>& group2index, const vector<string>& group_names, vector<vector<int>>& forced, vector<vector<float>>& alphas);

		/// \brief Shapley Lime with generator working iteratively
		void get_iterative_shapley_lime_params(const MedFeatures& data, const MedPredictor *model,
			SamplesGenerator<float> *generator, float p, int n, LimeWeightMethod weighting, float missing,
			void *params, const vector<vector<int>>& group2index, const vector<string>& group_names, const MedMat<float>& abs_cov_mat, int iteration_cnt, vector<vector<float>>& alphas, bool max_in_groups);
	}
}

#endif