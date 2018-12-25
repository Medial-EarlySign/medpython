#include "MedBART.h"

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL LOG_DEF_LEVEL

int MedBART::Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
	_model = BART(ntrees, iter_count, burn_count, restart_count, tree_params);

	vector<float> x_wrapper(x, x + nsamples*nftrs);
	vector<float> y_wrapper(y, y + nsamples);
	_model.learn(x_wrapper, y_wrapper);

	return 0;
}

int MedBART::Predict(float *x, float *&preds, int nsamples, int nftrs) const {
	vector<float> x_wrapper(x, x + nsamples*nftrs);
	vector<float> scores; //already allocated
	_model.predict(x_wrapper, nsamples, scores);
	//save values in preds:
	for (size_t i = 0; i < nsamples; ++i)
		preds[i] = scores[i];

	return 0;
}

int MedBART::set_params(map<string, string>& mapper) {
	//! [MedBART::init]
	unordered_map<string, int> map_types = {
		{ "classification", bart_data_prior_type::classification },
		{ "regression_mean_shift", bart_data_prior_type::regression_mean_shift }
	};
	string all_options = medial::io::get_list(map_types);

	for (auto it = mapper.begin(); it != mapper.end(); ++it)
	{
		if (it->first == "ntrees")
			ntrees = stoi(it->second);
		else if (it->first == "iter_count")
			iter_count = stoi(it->second);
		else if (it->first == "burn_count")
			burn_count = stoi(it->second);
		else if (it->first == "restart_count")
			restart_count = stoi(it->second);
		else if (it->first == "tree_params") {
			map<string, string> tree_args;
			MedSerialize::initialization_text_to_map(it->second, tree_args);
			for (auto ti = tree_args.begin(); ti != tree_args.end(); ++ti)
			{
				if (ti->first == "k")
					tree_params.k = stof(ti->second);
				else if (ti->first == "alpha")
					tree_params.alpha = stof(ti->second);
				else if (ti->first == "beta")
					tree_params.beta = stof(ti->second);
				else if (ti->first == "lambda")
					tree_params.lambda = stof(ti->second);
				else if (ti->first == "nu")
					tree_params.nu = stof(ti->second);
				else if (ti->first == "min_obs_in_node")
					tree_params.min_obs_in_node = stoi(ti->second);
				else if (ti->first == "data_prior_type") {
					if (map_types.find(ti->second) == map_types.end()) {
						MTHROW_AND_ERR("unsupported data_prior_type \"%s\". options are: %s\n",
							ti->second.c_str(), all_options.c_str());
					}
					else
						tree_params.data_prior_type = bart_data_prior_type(map_types.at(ti->second));
				}
				else
					MTHROW_AND_ERR("Unsupported argument \"%s\" for tree_params\n", ti->first.c_str());
			}
		}
		else
			MTHROW_AND_ERR("Unsupported argument \"%s\" for MedBart\n", it->first.c_str());
		if (tree_params.data_prior_type == classification)
			tree_params.set_classification(ntrees);
		//else if (tree_params.data_prior_type == regression_mean_shift)
		//	tree_params.set_regression(ntrees, );
	}
	//! [MedBART::init]
	return 0;
}