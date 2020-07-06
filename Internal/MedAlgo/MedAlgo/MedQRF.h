#pragma once
#include <MedAlgo/MedAlgo/MedAlgo.h>


//======================================================================================
// QRF: Quantized Regression/Classification random forest
//======================================================================================
#define MED_QRF_DEF_NTREES 100
#define MED_QRF_DEF_MAXQ 200
#define MED_QRF_DEF_MIN_NODE 50
#define MED_QRF_DEF_LEARN_NTHREADS 8
#define MED_QRF_DEF_PREDICT_NTHREADS 8
#define MED_QRF_DEF_SPREAD	0.1

struct MedQRFParams : public SerializableObject {

	// Required
	int ntrees;
	int maxq;
	int learn_nthreads, predict_nthreads;
	QRF_TreeType type;

	// Optional
	int max_samp; ///<M if > 0 & sampsize is NULL : the maximal sampsize we will take from each category
	float samp_factor; ///< if > 0 & sampsize if NULL : the maximal factor of samples between the 2 largest categories
	vector<int> samp_vec; ///< to be used when sampsize is NULL and max_samp,samp_vector > 0
	int *sampsize;
	int ntry; ///< if ntry <= 0: ntry = (int)(sqrt((double)nfeat) + 1.0);
	int get_only_this_categ;
	int max_depth; ///<maximial depth of tree branches - if 0 no limit
	bool take_all_samples; ///<use all samples - no sampling in building tree

						   // Regression
	float spread;
	bool keep_all_values; ///< For quantile regression
	bool sparse_values; ///< For keeping all values as a value-index(int):count(char) vector

						// categorical
	int min_node;
	int n_categ;

	int collect_oob;

	// For Prediction
	int get_count;
	vector<float> quantiles; ///< For quantile regression

	ADD_CLASS_NAME(MedQRFParams)
		ADD_SERIALIZATION_FUNCS(ntrees, maxq, learn_nthreads, predict_nthreads, type, max_samp, samp_factor, samp_vec,
			ntry, get_only_this_categ, max_depth, take_all_samples, spread, keep_all_values, sparse_values, min_node, n_categ, collect_oob, get_count, quantiles)
		void post_deserialization() { if (samp_vec.size() == 0) sampsize = NULL;  else sampsize = &samp_vec[0]; }

};

class MedQRF : public MedPredictor {
public:
	/// Model 
	QRF_Forest qf;

	/// Parameters
	MedQRFParams params;

	// Function
	MedQRF();
	~MedQRF() {};
	MedQRF(void *params);
	MedQRF(MedQRFParams& params);
	int init(void *params);
	/// The parsed fields from init command.
	/// @snippet MedQRF.cpp MedQRF::init
	virtual int set_params(map<string, string>& mapper);
	//	int init(const string &init_str); // allows init of parameters from a string. Format is: param=val,... , for sampsize: 0 is NULL, a list of values is separated by ; (and not ,)
	void init_defaults();

	/// @snippet MedQRF.cpp MedQRF_get_types
	QRF_TreeType get_tree_type(string name);

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	//int denormalize_model(float *f_avg, float *f_std, float lavel_avg, float label_std) {return 0;};

	// (De)Desrialize - virtual class methods that do the actuale (De)Serializing. Should be created for each predictor
	ADD_CLASS_NAME(MedQRF)
		ADD_SERIALIZATION_FUNCS(classifier_type, qf, params, model_features, features_count)

		// Print
		void print(FILE *fp, const string& prefix, int level = 0) const;
	void printTrees(const vector<string> &modelSignalNames, const string &outputPath) const;
	void calc_feature_importance(vector<float> &features_importance_scores, const string &general_params, const MedFeatures *features);

	// Predictions per sample
	int n_preds_per_sample() const;

	void prepare_predict_single();
	void predict_single(const vector<float> &x, vector<float> &preds) const;

private:
	void set_sampsize(float *y, int nsamples); // checking if there's a need to prep sampsize based on max_samp and samp_factor
	int Predict(float *x, float *&preds, int nsamples, int nftrs, int get_count) const;

	vector<pair<float, int>> _indexd_quantiles;
	vector<float> _sorted_quantiles;
	qrf_scoring_thread_params _single_pred_args;
	bool prepared_single;
};

MEDSERIALIZE_SUPPORT(MedQRFParams)
MEDSERIALIZE_SUPPORT(MedQRF)