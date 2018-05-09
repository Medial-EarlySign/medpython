#ifndef _MED_BART_H_
#define _MED_BART_H_
#include "MedAlgo.h"
#include "BART.h"

/**
* a wrapper for BART class model.
* BART is forest of binary trees - sum of trees.
* the model created MCMC and tries to maximize likelihood using stochstic process
*/
class MedBART : public MedPredictor {
public:

	/// <summary>
	/// an initialization for model
	/// @snippet MedBART.cpp MedBART::init
	/// </summary>
	void init_defaults() {};
	int set_params(map<string, string>& mapper);

	/// <summary>
	/// learning on x vector which represents matrix. y is the labels
	/// @param x a vector which represnts matrix. the data is ordered by observations. 
	/// we first see first observation all features and than second obseravtion all features...
	/// @param y labels vector for each observation in x
	/// </summary>
	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	/// <summary>
	/// prediction on x vector which represents matrix
	/// @param x a vector which represnts matrix. the data is ordered by observations. 
	/// we first see first observation all features and than second obseravtion all features...
	/// @param nsamples the number of samples in x
	/// @param nftrs the number of features in x
	/// </summary>
	/// <returns>
	/// @param preds the result scores for each observation
	/// </returns>
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	/// <summary>
	/// a simple default ctor
	/// </summary>
	MedBART() : _model(0, 0, 0, 0, tree_params) {
		normalize_for_learn = false;
		transpose_for_learn = false;
		normalize_y_for_learn = false;
		transpose_for_predict = false;
		normalize_for_predict = false;
		classifier_type = MODEL_BART;

		//default params:
		ntrees = 50;
		iter_count = 1000;
		burn_count = 250;
		restart_count = 5;
		tree_params.alpha = (float)0.95;
		tree_params.beta = 1;
		tree_params.min_obs_in_node = 5;

		tree_params.k = 2;
		tree_params.nu = 3;
		tree_params.lambda = 1;
	}

	ADD_CLASS_NAME(MedBART)
	ADD_SERIALIZATION_FUNCS(classifier_type, ntrees, iter_count, burn_count, restart_count)
private:
	int ntrees; ///< The nubmer of trees/restarts
	int iter_count; ///< the number of steps to call next_gen_tree for each tree
	int burn_count; ///< the burn count
	int restart_count; ///< number of restarts
	bart_params tree_params; ///< additional tree parameters

	BART _model;
};

MEDSERIALIZE_SUPPORT(MedBART)

#endif // !_MED_BART_H_

