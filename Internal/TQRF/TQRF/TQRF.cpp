#include "TQRF.h"

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//================================================================================================
// TQRF_Params
//================================================================================================
//------------------------------------------------------------------------------------------------
int TQRF_Params::init(map<string, string>& map)
{
	for (auto &f : map) {

		// general
		if (f.first == "samples_time_unit") samples_time_unit = f.second;

		// time slices
		else if (f.first == "time_slice_unit") time_slice_unit = f.second;
		else if (f.first == "time_slice_size") time_slice_size = stoi(f.second);
		else if (f.first == "time_slices") {
			vector<string> fields;
			boost::split(fields, f.second, boost::is_any_of(","));
			for (auto &ts : fields)
				time_slices.push_back(stoi(ts));
		}

		// quantization
		else if (f.first == "max_q") max_q = stoi(f.second);
		//else if (f.first == "max_q_sample") max_q_sample = stoi(f.second);

		// trees and stopping criteria
		if (f.first == "tree_type") tree_type = f.second;
		else if (f.first == "ntrees") ntrees = stoi(f.second);
		else if (f.first == "max_depth") max_depth = stoi(f.second);
		else if (f.first == "min_node_last_slice") min_node_last_slice = stoi(f.second);
		else if (f.first == "min_node") min_node = stoi(f.second);

		// speedup by subsample control
		else if (f.first == "max_node_test_samples") max_node_test_samples = stoi(f.second);

		// bagging control
		else if (f.first == "single_sample_per_pid") single_sample_per_pid = stoi(f.second);
		else if (f.first == "bag_with_repeats") bag_with_repeats = stoi(f.second);
		else if (f.first == "bag_prob") bag_prob = stof(f.second);
		else if (f.first == "bag_ratio") bag_ratio = stof(f.second);
		else if (f.first == "bag_feat") bag_feat = stof(f.second);

		// missing value
		else if (f.first == "missing_val") missing_val = stof(f.second);
		else if (f.first == "missing_method_str") {
			missing_method_str = f.second;
			missing_method = TQRF_Forest::get_missing_value_method(missing_method_str);
		}

		// categorial features
		else if (f.first == "nvals_for_categorial") nvals_for_categorial = stoi(f.second);
		else if (f.first == "categorial_str") { boost::split(categorial_str, f.second, boost::is_any_of(",")); }
		else if (f.first == "categorial_tags") { boost::split(categorial_tags, f.second, boost::is_any_of(",")); }

		// sanities
		else if (f.first == "test_for_inf") test_for_inf = stoi(f.second);
		else if (f.first == "test_for_missing") test_for_missing = stoi(f.second);

		// sanities
		else if (f.first == "verbosity") verbosity = stoi(f.second);

		else {
			MERR("TQRF_Params::init(): No such parameter \'%s\'\n", f.first.c_str());
			return -1;
		}
	}

	tree_type_i = TQRF_Forest::get_tree_type(tree_type);

	return 0;
}

//================================================================================================
// Quantized_Feat
//================================================================================================
//------------------------------------------------------------------------------------------------
// quantize_feat(int i_feat, TQRF_Params &params) :
// assumes: nfeat, orig_data, ncateg, y are initialized and qx, q_to_val initialized to right size.
// actually doing the quantization:
// (1) prepare a value list and sort it.
// (2) count the number of different values.
// (3??) in case of categorial and adjacent values with a single same category - unite them. (Not relevant at the moment due to time slicing)
// (4) if left with less than max_q values : we finished
// (5) Otherwise: we split given the split strategy (currently - even splits. Future Idea: better subsampling of tails/information rich areas)
//------------------------------------------------------------------------------------------------
int Quantized_Feat::quantize_feat(int i_feat, TQRF_Params &params)
{
	vector<pair<float,int>> curr_data;
	vector<float> *data;
	data = orig_data[i_feat];

	if (params.verbosity > 0) MLOG("Quantized_Feat::quantize_feat :: feat %d %s :: nvals %d\n", i_feat, feat_names[i_feat].c_str(), data->size());

	curr_data.resize(data->size());
	for (int i=0; i<(int)data->size(); i++)
		curr_data[i] = { (*data)[i], i };

	// sorting for easier quantization later
	sort(curr_data.begin(), curr_data.end());

	// count diff values
	vector<pair<float, int>> v_counts;
	float prev = curr_data[0].first;
	int n = 0;
	for (auto v : curr_data) {
		if (v.first == prev)
			n++;
		else {
			v_counts.push_back({ prev,n });
			prev = v.first;
			n = 1;
		}
	}
	v_counts.push_back({ prev,n }); // last count

	if (params.verbosity > 0) MLOG("Quantized_Feat::quantize_feat :: feat %d %s :: ndiff %d (max_q %d)\n", i_feat, feat_names[i_feat].c_str(), v_counts.size(), params.max_q);

	// Step2 of calculating which features are categorial
	if (v_counts.size() < params.nvals_for_categorial) is_categorial_feat[i_feat] = 1;

	// prepare q_to_val
	// the convention is that if a value is in the range (q_to_val[i-1],q_to_val[i]] then its quantized value will be i
	if (v_counts.size() <= params.max_q) {
		q_to_val[i_feat].clear();
		for (auto &v : v_counts)
			q_to_val[i_feat].push_back(v.first);
	}
	else {
		// need more work to be done, we have too many different values... 
		// our algorithm is simply trying to make the i-th qval such that there are (i+1)/max_q * len , of the elements
		// up to some fixes that are caused due to integer numbers
		int delta = (int)((float)curr_data.size()/(float)params.max_q);
		int j = 0;
		int len = (int)curr_data.size();
		while (j < len-1) {
			j = j + delta;
			if (j >= len) j = len - 1;
			float q = curr_data[j].first;
			while ((j < len-1) && (curr_data[j+1].first == q)) j++;
			q_to_val[i_feat].push_back(q);
		}
	}

	if (params.verbosity > 1) {
		MLOG("Quantized_Feat::quantize_feat :: %d %s :: q_to_val size %d :: ", i_feat, feat_names[i_feat].c_str(), q_to_val.size());
		for (auto q : q_to_val[i_feat])
			MLOG(" %f :", q);
		MLOG("\n");
	}

	// now q_to_val is ready and all that is left is to actually use it to create qx
	int q_i = 0;
	int q_size = (int)q_to_val[i_feat].size();
	qx[i_feat].resize(curr_data.size());
	for (auto &v : curr_data) {
		if (q_i < q_size && v.first >= q_to_val[i_feat][q_i])
			q_i++;
		qx[i_feat][v.second] = (short)q_i;
	}

	return 0;
}

//------------------------------------------------------------------------------------------------
// Quantized_Feat::init(const MedFeatures &medf, TQRF_Params &params) :
// initializations needed for all trees together, mainly data quantization , time slices, and
// filling in all the needed variables in Quantized_Feat
//------------------------------------------------------------------------------------------------
int Quantized_Feat::init(MedFeatures &medf, TQRF_Params &params)
{
	// filling in needed variables
	orig_medf = &medf;
	orig_data.clear();
	for (auto &df : medf.data) {
		orig_data.push_back(&df.second);
	}
	medf.get_feature_names(feat_names);
	nfeat = (int)(medf.data.size());

	// init y and get ncateg if needed
	y.resize(medf.samples.size());
#pragma omp parallel
	for (int i=0; i<medf.samples.size(); i++)
		y[i] = medf.samples[i].outcome;

	if (params.tree_type_i != TQRF_TREE_REGRESSION) {
		float maxy = *max_element(y.begin(), y.end());
		ncateg = (int)maxy + 1;
		// TBD :: sanity test that the y's in this case are all in the 0,ncateg-1 integers range ... currently assuming this
	}

	// time slices
	if (params.time_slices.size() == 0) {
		if (params.time_slice_size > 0)
			for (int i=0; i<params.n_time_slices; i++)
				params.time_slices.push_back(i*params.time_slice_size);
		else
			params.time_slices.push_back(TQRF_MAX_TIME_SLICE); // single infinite cell
	}
	// slices for each sample
	// we assume that '0' outcome is the control, for binary and multicategory cases
	// The regression case for time slices is much more complex, and needs more input from the user: currently a project to the future
	last_time_slice.resize(medf.samples.size());
#pragma omp parallel
	for (int i=0; i<medf.samples.size(); i++) {
		for (int j=0; j<params.time_slices.size(); j++)
			if (medf.samples[i].outcomeTime < params.time_slices[j]) {
				last_time_slice[i] = j;
				break;
			}
	}

	// categorial features Step1 (Step2 will be done in quantize_feat)
	is_categorial_feat.clear();
	is_categorial_feat.resize(nfeat, 0);
	for (int i=0; i<nfeat; i++) {
		for (int j=0; j<params.categorial_str.size(); j++) {
			if (feat_names[i].find(params.categorial_str[j]) != string::npos) {
				is_categorial_feat[i] = 1;
				break;
			}
		}

		for (int j=0; j<params.categorial_tags.size(); j++) {
			if (medf.tags[feat_names[i]].find(params.categorial_tags[j]) != medf.tags[feat_names[i]].end()) {
				is_categorial_feat[i] = 1;
				break;
			}
		}
	}	

	int rc = 0;
	// quantization of features
#pragma omp parallel
	for (int i=0; i<nfeat; i++)
		if (quantize_feat(i, params) < 0)
			rc = -1;

	// rc
	return rc;
}


//================================================================================================
// TQRF_Forest
//================================================================================================
//------------------------------------------------------------------------------------------------
int TQRF_Forest::get_tree_type(const string &str)
{
	if (boost::iequals(str, "entropy")) return TQRF_TREE_ENTROPY;
	if (boost::iequals(str, "logrank")) return TQRF_TREE_LOGRANK;
	if (boost::iequals(str, "regression")) return TQRF_TREE_REGRESSION;

	return TQRF_TREE_UNDEFINED;
}

//------------------------------------------------------------------------------------------------
int TQRF_Forest::get_missing_value_method(const string &str)
{
	if (boost::iequals(str, "mean")) return TQRF_MISSING_VALUE_MEAN;
	if (boost::iequals(str, "median")) return TQRF_MISSING_VALUE_MEDIAN;
	if (boost::iequals(str, "larger")) return TQRF_MISSING_VALUE_LARGER_NODE;
	if (boost::iequals(str, "left")) return TQRF_MISSING_VALUE_LEFT;
	if (boost::iequals(str, "right")) return TQRF_MISSING_VALUE_RIGHT;
	if (boost::iequals(str, "rand")) return TQRF_MISSING_VALUE_RAND_ALL;
	if (boost::iequals(str, "rand_each_sample")) return TQRF_MISSING_VALUE_RAND_EACH_SAMPLE;

	return TQRF_MISSING_VALUE_NOTHING;
}


//================================================================================================
// TQRF_Tree
//================================================================================================
//--------------------------------------------------------------------------------------------------------------------
int TQRF_Tree::get_bagged_indexes()
{
	unordered_set<int> s_pids[2] , s_indexes[2];
	vector<int> v_pids[2], v_indexes[2];
	unordered_map<int, vector<int>> pid2indexes;

	bool is_regression = (tree_type != TQRF_TREE_REGRESSION);

	for (int i=0; i<_qfeat.orig_medf->samples.size(); i++) {
		MedSample &s = _qfeat.orig_medf->samples[i];
		int j = 0;
		if (is_regression || (s.outcome != 0)) j=1;

		if (_params.single_sample_per_pid) {
			if (s_pids[j].find(s.id) == s_pids[j].end()) { v_pids[j].push_back(s.id); s_pids[j].insert(s.id); }
			if (pid2indexes.find(s.id) == pid2indexes.end()) pid2indexes[s.id] = vector<int>();
			pid2indexes[s.id].push_back(i);
		}
		else {
			if (s_indexes[j].find(s.id) == s_indexes[j].end()) { v_indexes[j].push_back(i); s_indexes[j].insert(i); }
		}
	}

	if (tree_type != TQRF_TREE_REGRESSION) {

		// calculate the bagging probabilities for 0 and 1
		int n0, n1;
		float p0, p1;

		if (_params.bag_ratio < 0) {
			p0 = _params.bag_prob;
			p1 = _params.bag_prob;
		}
		else {
			if (_params.single_sample_per_pid) {
				n0 = (int)s_pids[0].size();
				n1 = (int)s_pids[1].size();
			}
			else {
				n0 = (int)s_indexes[0].size();
				n1 = (int)s_indexes[1].size();
			}

			//
			// following calculations use the fact that we want : bag_ratio = (n0*p0)/(n1*p1)
			//
			if (n0 > n1) {
				p1 = _params.bag_prob;
				p0 = _params.bag_ratio * (float)(n1+1) * p1 / (float)(n0+1);
			}
			else {
				p0 = _params.bag_prob;
				p1 = p0 * (float)(n0+1) / (_params.bag_ratio * (float)(n1+1));
			}

			p0 = min(p0, (float)1);
			p1 = min(p1, (float)1);
		}

		indexes.clear();

		bag_chooser(_params.bag_with_repeats, _params.single_sample_per_pid, p0, v_pids[0], v_indexes[0], pid2indexes, indexes);
		bag_chooser(_params.bag_with_repeats, _params.single_sample_per_pid, p1, v_pids[1], v_indexes[1], pid2indexes, indexes);

	}
	else {
		// in regression all samples are with j==1
		bag_chooser(_params.bag_with_repeats, _params.single_sample_per_pid, _params.bag_prob, v_pids[1], v_indexes[1], pid2indexes, indexes);
	}

}

//--------------------------------------------------------------------------------------------------------------------
int TQRF_Tree::bag_chooser(int choose_with_repeats, int single_sample_per_id, float p, vector<int> &pids, vector<int> &idx, unordered_map<int, vector<int>> &pid2idx, /* OUT APPEND */ vector<int> &_indexes)
{
	int n_pids = (int)pids.size();
	int n_idx = (int)idx.size();

	if (single_sample_per_id) {

		if (choose_with_repeats) {
			unordered_map<int, int> pid_idx;
			int n_choose = (int)(p*(float)n_pids);
			for (int i=0; i<n_choose; i++) {
				int rand_pid = pids[rand_N(n_pids)];
				if (pid_idx.find(rand_pid) == pid_idx.end()) {
					int len_pid = (int)pid2idx[rand_pid].size();
					int rand_idx = rand_N(len_pid);
					pid_idx[rand_pid] = pid2idx[rand_pid][rand_idx];
				}
				_indexes.push_back(pid_idx[rand_pid]); // making sure to choose the SAME sample per id
			}
		}
		else {
			for (int i=0; i<n_pids; i++)
				if (rand_1() < p)
					_indexes.push_back(pid2idx[pids[i]][rand_N((int)pid2idx[pids[i]].size())]);
		}

	}
	else {
		if (choose_with_repeats) {

			int n_choose = (int)(p*(float)n_idx);
			for (int i=0; i<n_choose; i++)
				_indexes.push_back(idx[rand_N(n_idx)]);

		}
		else {
			for (int i=0; i<n_idx; i++)
				if (rand_1() < p)
					_indexes.push_back(idx[i]);
		}
	}
}
