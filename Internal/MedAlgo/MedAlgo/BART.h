#ifndef BART_H
#define BART_H
#include <vector>
#include <random>
#include <unordered_map>
#include <boost/random.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;

/**
* bart tree node
*/
class bart_node {
public:
	int feature_number; ///< feature number in node for split
	int split_index; ///< the split index of the feature sorted unique value as split value in node
	float node_value; ///< the output value for input that reaches this node
	bart_node *parent; ///< the parent node
	bart_node *childs[2]; ///< the left,right childs
	vector<int> observation_indexes; ///< the indexes of observations in this node
	bool mark_change; ///< mark change for calculating node_value

	int num_feature_options; ///< number of features to select for likelihood calc
	int num_split_options; ///< number of split value to select after feature select for likelihood calc

	/// <summary>
	/// a simple default ctor
	/// </summary>
	bart_node() {
		feature_number = -1;
		split_index = -1;
		node_value = 0;
		parent = NULL; //root
					   //Childrens
		childs[0] = NULL; //smaller than or equal
		childs[1] = NULL; //biger than
		mark_change = true;
		num_feature_options = 0;
		num_split_options = 0;
	}

	/// <summary>
	/// a copy ctor shallow copy. not cloing or allocating memory for childs\parents nodes
	/// </summary>
	bart_node(const bart_node &cp) {
		feature_number = cp.feature_number;
		split_index = cp.split_index;
		node_value = cp.node_value;
		parent = cp.parent;
		childs[0] = cp.childs[0];
		childs[1] = cp.childs[1];
		observation_indexes = cp.observation_indexes;
		mark_change = cp.mark_change;
		num_feature_options = cp.num_feature_options;
		num_split_options = cp.num_split_options;
	}

	/// <summary>
	/// populating all_nodes with flatted array of all nodes from current node including this node.
	/// </summary>
	void list_all_nodes(vector<bart_node *> &all_nodes);
	/// <summary>
	/// returning node depth
	/// </summary>
	int depth();
	/// <summary>
	/// calculating node variance in the observations in the node
	/// </summary>
	float variance(const vector<float> &y);
	/// <summary>
	/// deep copying of the node and all it's acendents by allocating new copies of all nodes
	/// </summary>
	void deep_clone(bart_node *&target);
	/// <summary>
	/// printing tree from current node
	/// </summary>
	void print_tree(const vector<vector<float>> &feature_sorted_values) const; //for debug
	/// <summary>
	/// for debug - validating the tree structure is correct with childs\parents pointers and that the
	/// observations in each nodes are correct.
	/// </summary>
	void validate_tree(const vector<vector<float>> &feature_sorted_values, const vector<float> &x, int nftrs) const; //for debug
private:
};

/**
* A Class to represnet change in tree - for rollback
* or release memory in commit
*/
class tree_change_details {
public:
	vector<bart_node *>changed_nodes_before; ///< the node before changes
	vector<bart_node *>changed_nodes_after; ///<the node after changes

	int action = -1; ///< the action
	int num_node_selection = 0; ///< the number of nodes to select
};

enum bart_data_prior_type {
	regression_mean_shift = 0,
	classification = 1
};

class bart_tree;

/**
* bart tree parameters
*/
class bart_params {
public:
	int min_obs_in_node; ///< minimal allowed observations in node

	float alpha; ///< prior for tree structure: alpha * (1 + depth(node)) ^ -beta
	float beta; ///< prior for tree structure: alpha * (1 + depth(node)) ^ -beta

				//general:
	bart_data_prior_type data_prior_type;

	//params for prior:
	float k; ///< the range for bandwidth interval
	float nu; ///< the node-data dict params for sigma_i: sigma_i ~ IG(nu, mean_sigma*lambda/2)
	float lambda; ///< the node-data dict params for sigma_i: sigma_i ~ IG(mean_sigma/2, mean_sigma*lambda/2)

	/// <summary>
	/// a simple default ctor
	/// </summary>
	bart_params() {
		min_obs_in_node = 0;
		alpha = 1;
		beta = 1;
		sigsq_mu = 0;
		mean_mu = 0;
		nu = 3;
		k = 2;
		lambda = 1;
		data_prior_type = regression_mean_shift;
	}

	/// <summary>
	/// an initializer for classification problems
	/// </summary>
	void set_classification(int num_trees) {
		mean_mu = 0;
		//k = 2;
		//nu = 3.0;
		data_prior_type = bart_data_prior_type::classification;

		sigsq_mu = pow(3 / (k * sqrt(num_trees)), 2);
	}
	/// <summary>
	/// an initializer for regression problems. 
	/// @param sample_var_y_in_data it's the variance in train labels
	/// </summary>
	void set_regression(int num_trees, float sample_var_y_in_data) {
		mean_mu = 0;
		//nu = 3.0;
		//k = 2;
		data_prior_type = bart_data_prior_type::regression_mean_shift;

		sigsq_mu = pow(1 / (2 * k * sqrt(num_trees)), 2);

		boost::math::chi_squared_distribution<> chi_dist(nu);
		double ten_pctile_chisq_df_hyper_nu = boost::math::cdf(chi_dist, 1 - 0.9);
		lambda = sample_var_y_in_data * ten_pctile_chisq_df_hyper_nu / nu;
		//lambda = 2.7e-3;
		//printf("ten_pctile_chisq_df_hyper_nu=%2.5f, var=%2.5f\n", ten_pctile_chisq_df_hyper_nu, sample_var_y_in_data);
	}
private:
	float mean_mu; ///< the node-data dist params for Mi: (Mi,sigma_i) ~ N(mean_mu, sigma_i/a); mean_mu set to zero always
	float sigsq_mu; ///< the node-data dist params for Mi: (Mi,sigma_i) ~ N(mean_mu, sigma_i/a);

	friend class bart_tree;
};

class BART;

/**
* bart tree
*/
class bart_tree {
public:
	bart_node *root; ///< the tree root
	double tree_loglikelihood; ///< the tree likelihood based on tree prior and tree match to data
	bart_params params; ///< the barat params

	/// <summary>
	/// creating next move in MCMC from current tree using metropolis hasting algorithm
	/// @param x the matrix observation by observation in one vector
	/// @param y the labels for each observation
	/// </summary>
	void next_gen_tree(const vector<float> &x, const vector<float> &y);

	/// <summary>
	/// a function to clone bart_tree from root node deep clone - creating
	/// a copy of each of the nodes data
	/// </summary>
	void clone_tree(bart_tree &tree);

	/// <summary>
	/// setting the parameter sigma - should happen before generation new tree using
	/// metropolis hasting
	/// </summary>
	void set_sigma(double sig) {
		sigma = sig;
	}

	/// <summary>
	/// a copy ctor preforming a shallow copy: not allocating memory for all tree nodes
	/// again. usibng the same pointers
	/// </summary>
	bart_tree(const bart_tree &other) {
		root = other.root;
		tree_loglikelihood = other.tree_loglikelihood;
		_rnd_gen = other._rnd_gen;
		action_priors = other.action_priors;
		feature_to_sorted_vals = other.feature_to_sorted_vals;
		feature_to_val_index = other.feature_to_val_index;
		params = other.params;
	}

	/// <summary>
	/// a simple default ctor
	/// </summary>
	bart_tree() {
		root = NULL;
		tree_loglikelihood = -1;
		_rnd_gen = mt19937(_rd());

		//Default Values:
		action_priors = { (float)0.25, (float)0.25, (float)0.4, (float)0.1 };
		params.min_obs_in_node = 5;
		params.alpha = (float)0.95;
		params.beta = 1.0;

		params.data_prior_type = bart_data_prior_type::regression_mean_shift;
		params.mean_mu = 0;
		params.nu = 3;
		params.sigsq_mu = 0;
		params.lambda = 1.0;
	}
private:
	//general variables
	random_device _rd;
	mt19937 _rnd_gen;
	vector<vector<float>> feature_to_sorted_vals;
	vector<unordered_map<float, int>> feature_to_val_index;
	double sigma;

	vector<float> action_priors;
	//actions:
	tree_change_details do_grow(const vector<float> &x, const vector<float> &y);
	tree_change_details do_prune(const vector<float> &x, const vector<float> &y);
	tree_change_details do_change(const vector<float> &x, const vector<float> &y);
	tree_change_details do_swap(const vector<float> &x, const vector<float> &y);

	float score_leaf(const vector<float> &y, const vector<int> &obs_indexes);

	double node_data_likelihood(const vector<bart_node *> &leaf_node, const vector<float> &x, int nftrs, const vector<float> &y);
	void calc_likelihood(const vector<float> &x, int nftrs, const vector<float> &y); //will calc before mean in all leaves

																					 //helper for change:
	void get_avalible_change(const vector<float> &x, int nftrs, vector<bart_node *> &good_idx);
	void get_avalible_feats_change(const bart_node *selected_node, const vector<float> &x, int nftrs, vector<int> &good_idx);
	void propogate_change_down(bart_node *current_node, const vector<float> &x, int nftrs, vector<bart_node *> &list_nodes_after);

	//helper for grow:
	void get_avalible_grow(const vector<float> &x, int nftrs, vector<bart_node *> &good_idx) const;
	void get_avalible_feats_grow(const bart_node *selected_node, const vector<float> &x, int nftrs, vector<int> &good_idx);

	//helper for prune;
	void get_avalible_prune(const vector<float> &x, int nftrs, vector<bart_node *> &good_idx);

	//helper for swap:
	void get_avalible_swap(const vector<float> &x, int nftrs, vector<bart_node *> &good_idx);

	bool has_split(const bart_node *current_node, const vector<float> &x, int nftrs, int feature_num) const;
	bool select_split(bart_node *current_node, const vector<float> &x, int nftrs,
		vector<vector<int>> &split_obx_indexes);

	void commit_change(const tree_change_details &change);
	void rollback_change(const tree_change_details &change);

	int clear_tree_mem(bart_node *node); //erase count
protected:
	void predict(const vector<float> &x, int nSamples, vector<float> &scores) const;
	void predict_on_train(const vector<float> &x, int nSamples, vector<float> &scores) const; //faster

	friend class BART;
};

/**
* Bayesian Additive Regression Trees.\n
* A Monte-Carlo-Markov-Chain process to create trees based on bayesian assumtions\n
* to reach maximum likelihood of tree based on the observations with Metropolis_Hastings algorithm.\n
* It's very usefull for casual inference:\n
* In 2016 at Atlantic Casual Infernce Confernce - 1st place in data competition
*/
class BART {
public:
	int ntrees; ///< The nubmer of trees/restarts
	int iter_count; ///< the number of steps to call next_gen_tree for each tree
	int burn_count; ///< the burn count
	int restart_count; ///< number of restarts
	bart_params tree_params; ///< additional tree parameters

	/// <summary>
	/// learning on x vector which represents matrix. y is the labels
	/// @param x a vector which represnts matrix. the data is ordered by observations. 
	/// we first see first observation all features and than second obseravtion all features...
	/// @param y labels vector for each observation in x
	/// </summary>
	void learn(const vector<float> &x, const vector<float> &y);

	/// <summary>
	/// prediction on x vector which represents matrix
	/// @param x a vector which represnts matrix. the data is ordered by observations. 
	/// we first see first observation all features and than second obseravtion all features...
	/// @param nSamples the number of samples in x
	/// </summary>
	/// <returns>
	/// @param scores the result scores for each observation
	/// </returns>
	void predict(const vector<float> &x, int nSamples, vector<float> &scores) const;

	/// <summary>
	/// a simple default ctor
	/// </summary>
	BART(int ntrees, int iterations, int burn_cnt, int restart_cnt, bart_params &tree_pr) {
		nftrs = 0;
		//default:
		this->ntrees = ntrees;
		this->iter_count = iterations;
		this->burn_count = burn_cnt;
		this->tree_params = tree_pr;
		this->restart_count = restart_cnt;

		_trees.resize(ntrees);
		for (size_t i = 0; i < ntrees; ++i)
			_trees[i].params = this->tree_params;
		trans_y_b = 0;
		trans_y_max = 0;
	}

	/// <summary>
	/// a simple assignment operator to shallow copy all BART model with all trees.
	/// not allocating new memory for trees. pointing to same objects
	/// </summary>
	void operator=(const BART &other) {
		ntrees = other.ntrees; iter_count = other.iter_count;
		burn_count = other.burn_count; restart_count = other.restart_count;
		tree_params = other.tree_params;
		nftrs = other.nftrs;
		_trees.resize(other._trees.size());
		for (size_t i = 0; i < other._trees.size(); ++i) {
			_trees[i].action_priors = other._trees[i].action_priors;
			_trees[i].root = other._trees[i].root;
			_trees[i].tree_loglikelihood = other._trees[i].tree_loglikelihood;
			_trees[i]._rnd_gen = other._trees[i]._rnd_gen;
			_trees[i].feature_to_sorted_vals = other._trees[i].feature_to_sorted_vals;
			_trees[i].feature_to_val_index = other._trees[i].feature_to_val_index;
			_trees[i].params = other._trees[i].params;
		}
	}

	/// <summary>
	/// a dctor to free all tree memory
	/// </summary>
	~BART() {
		for (size_t i = 0; i < _trees.size(); ++i)
			_trees[i].clear_tree_mem(_trees[i].root);
	}

private:
	int nftrs;
	vector<bart_tree> _trees;
	float trans_y_b; //movement of y values
	float trans_y_max;

	void transform_y(vector<float> &y);
	void untransform_y(vector<float> &y) const;

	void init_hyper_parameters(const vector<float> &residuals);
	void update_sigma_param(boost::mt19937 &rng, const vector<float> &residuals, double &sigma); //for regression
	void update_latent_z_params(boost::random::random_number_generator<boost::mt19937> &rng_gen,
		const vector<float> &x, const vector<float> &y, const vector<bart_tree> &forest_trees,
		vector<float> &residuals); //for classification
};

#endif // !BART_H