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

		MLOG("%s => %s %d\n", f.first.c_str(), f.second.c_str(), (f.first == "time_slice_size"));

		// general
		if (f.first == "samples_time_unit") samples_time_unit = f.second;

		// time slices
		else if (f.first == "time_slice_unit") time_slice_unit = f.second;
		else if (f.first == "time_slice_size") {
			time_slice_size = stoi(f.second); MLOG("#####\n");		}
		else if (f.first == "n_time_slices") n_time_slices = stoi(f.second);
		else if (f.first == "time_slices") {
			vector<string> fields;
			boost::split(fields, f.second, boost::is_any_of(","));
			for (auto &ts : fields)
				time_slices.push_back(stoi(ts));
			n_time_slices = (int)time_slices.size();
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
		else if (f.first == "ids_to_print") ids_to_print = stoi(f.second);

		else {
			MERR("TQRF_Params::init(): No such parameter \'%s\'\n", f.first.c_str());
			//return -1;
		}
	}

	tree_type_i = TQRF_Forest::get_tree_type(tree_type);
	samples_time_unit_i = med_time_converter.string_to_type(samples_time_unit);
	time_slice_unit_i = med_time_converter.string_to_type(time_slice_unit);

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

	//if (params.verbosity > 0) MLOG("Quantized_Feat::quantize_feat :: feat %d %s :: nvals %d\n", i_feat, feat_names[i_feat].c_str(), data->size());

	curr_data.resize(data->size());
	for (int i=0; i<(int)data->size(); i++)
		curr_data[i] = { (*data)[i], i };


	if (params.verbosity > 2) {
		for (int i=0; i<params.ids_to_print; i++)
			MLOG("Quantuized Feat (unsorted) :: feat %d :: i %d :: curr_data %f %d\n", i_feat, i, curr_data[i].first, curr_data[i].second);
	}

	// sorting for easier quantization later
	sort(curr_data.begin(), curr_data.end());

	if (params.verbosity > 2) {
		for (int i=0; i<params.ids_to_print; i++)
			MLOG("Quantuized Feat (sorted) :: feat %d :: i %d :: curr_data %f %d\n", i_feat, i, curr_data[i].first, curr_data[i].second);
	}


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
			int prev_j = j;
			j = j + delta;
			if (j >= len) j = len - 1;
			float q = curr_data[j].first;
			while ((j < len-1) && (curr_data[j+1].first == q)) j++;
			if (params.verbosity > 1) MLOG("%s : %f : j %d size %d\n", feat_names[i_feat].c_str(), q, j, j-prev_j);
			q_to_val[i_feat].push_back(q);
		}
	}

	if (params.verbosity > 1) {
		MLOG("Quantized_Feat::quantize_feat :: %d %s :: q_to_val size %d :: ", i_feat, feat_names[i_feat].c_str(), q_to_val[i_feat].size());
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
	feat_names.clear();
	for (auto &df : medf.data) {
		orig_data.push_back(&df.second);
		feat_names.push_back(df.first);
		if (params.verbosity > 1) MLOG("Quantized_Feat:: %s :: %d elements :: %f %f %f %f ....\n", feat_names.back().c_str(), orig_data.back()->size(),
			(*orig_data.back())[0], (*orig_data.back())[1], (*orig_data.back())[2], (*orig_data.back())[3]);
	}


	nfeat = (int)(medf.data.size());

	// init y and get ncateg if needed
	y.resize(medf.samples.size());
#pragma omp parallel for
	for (int i=0; i<medf.samples.size(); i++)
		y[i] = medf.samples[i].outcome;

	if (params.tree_type_i != TQRF_TREE_REGRESSION) {
		float maxy = *max_element(y.begin(), y.end());
		ncateg = (int)maxy + 1;
		// TBD :: sanity test that the y's in this case are all in the 0,ncateg-1 integers range ... currently assuming this
	}
	else
		ncateg = 0;

	// time slices
	if (init_time_slices(medf, params) < 0)
		return -1;

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


	// init q arrays to right size
	qx.resize(nfeat);
	q_to_val.resize(nfeat);

	int rc = 0;
	// quantization of features
#pragma omp parallel for
	for (int i=0; i<nfeat; i++)
		if (quantize_feat(i, params) < 0)
			rc = -1;

	// rc
	return rc;
}


//------------------------------------------------------------------------------------------------
int Quantized_Feat::init_time_slices(MedFeatures &medf, TQRF_Params &params)
{
	if (params.tree_type_i != TQRF_TREE_REGRESSION) {
		if (params.time_slices.size() == 0) {
			if (params.time_slice_size > 0)
				for (int i=0; i<params.n_time_slices; i++)
					params.time_slices.push_back((i+1)*params.time_slice_size);
			else
				params.time_slices.push_back(TQRF_MAX_TIME_SLICE); // single infinite cell
		}
	}
	else {
		// regression tree case - currently using a single infinity time slice
		params.time_slices.push_back(TQRF_MAX_TIME_SLICE); // single infinite cell
	}

	n_time_slices = (int)params.time_slices.size();

	if (params.verbosity > 0) {
		MLOG("time_slices (%d) : ", params.time_slices.size());
		for (auto t : params.time_slices)
			MLOG("%d ", t);
		MLOG("\n");
	}

	// slices for each sample
	// we assume that '0' outcome is the control, for binary and multicategory cases
	// The regression case for time slices is much more complex, and needs more input from the user: currently a project to the future
	last_time_slice.resize(medf.samples.size());
#pragma omp parallel for
	for (int i=0; i<medf.samples.size(); i++) {
		MedSample &s = medf.samples[i];
		int d1 = med_time_converter.convert_times(params.samples_time_unit_i, params.time_slice_unit_i, (int)s.outcomeTime);
		int d2 = med_time_converter.convert_times(params.samples_time_unit_i, params.time_slice_unit_i, (int)s.time);
		int time_diff = d1-d2;
		int j = 0;
		for (j=0; j<params.time_slices.size(); j++)
			if (time_diff < params.time_slices[j])
				break;
		if (y[i] == 0 || j<params.time_slices.size()) {
			last_time_slice[i] = min(j, (int)params.time_slices.size()-1);
		}
		else {
			// interesting case in which we have y[i] != 0 that ended AFTER the LAST time slice.
			// in this case we have to treat the sample as a 0 outcome (!)
			last_time_slice[i] = (int)params.time_slices.size()-1;
			y[i] = 0;
		}

		if (params.verbosity > 0 && i < params.ids_to_print)
			MLOG("sample %d : pid %d otime %d outcome %f time %d : d1 %d d2 %d diff %d : last_time_slice %d y %f (%d)\n", i, s.id, s.outcomeTime, s.outcome, s.time, d1, d2, time_diff, last_time_slice[i], y[i], y[i] != s.outcome);
	}

	// getting stats for each time slice
	vector<int> cnt[2];
	cnt[0].resize(n_time_slices, 0);
	cnt[1].resize(n_time_slices, 0);
	for (int i=0; i<last_time_slice.size(); i++) {
		int j;
		for (j=0; j<last_time_slice[i]; j++)
			cnt[0][j]++;
		if (y[i] == 0)
			cnt[0][j]++;
		else
			cnt[1][j]++;
	}
	for (int j=0; j<n_time_slices; j++) {
		if (params.verbosity > 0)
			MLOG("Slice %d : %d : 0: %d 1: %d p: %f\n", j, params.time_slices[j], cnt[0][j], cnt[1][j], (float)cnt[1][j]/(1+cnt[0][j]));
		if (cnt[0][j] < MIN_ELEMENTS_IN_TIME_SLICE || cnt[1][j] < MIN_ELEMENTS_IN_TIME_SLICE) {
			MWARN("TQRF: WARNING : time slice %d (%d) too small or non variable :0: %d 1: %d p: %f\n", j, params.time_slices[j], cnt[0][j], cnt[1][j], (float)cnt[1][j]/(1+cnt[0][j]));
		}
	}

	return 0;
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


//------------------------------------------------------------------------------------------------
int TQRF_Forest::Train(MedFeatures &medf, const MedMat<float> &Y)
{
	// first - quantifying data
	Quantized_Feat qfeat;
	qfeat.init(medf, params);

	// creating the trees and parallelize train on each tree
	trees.resize(params.ntrees);
	for (int i=0; i<params.ntrees; i++) {
		trees[i].id = i;
		trees[i].tree_type = params.tree_type_i;
		trees[i].init(qfeat, params);
		trees[i].Train();
	}

	return 0;
}

//================================================================================================
// TQRF_Tree
//================================================================================================
//--------------------------------------------------------------------------------------------------------------------
int TQRF_Tree::get_bagged_indexes()
{
	unordered_set<int> s_pids[2];
	vector<int> v_pids[2], v_indexes[2];
	unordered_map<int, vector<int>> pid2indexes;

	bool is_regression = (tree_type == TQRF_TREE_REGRESSION);

	if (_params->verbosity > 0) MLOG("Tree %d %s : bagging : params: bag_prob %f bag_ratio %f single_per_pid %d bag_with_repeats %d\n", id, _params->tree_type.c_str(), _params->bag_prob, _params->bag_ratio, _params->single_sample_per_pid, _params->bag_with_repeats);

	int cnt[2] = { 0,0 };
	for (int i=0; i<_qfeat->orig_medf->samples.size(); i++) {
		MedSample &s = _qfeat->orig_medf->samples[i];

		//if (_params->verbosity > 0 && i<_params->ids_to_print) MLOG("i %d : s : pid %d outcome_time %d outcome %f time %d\n", i, s.id, s.outcomeTime, s.outcome, s.time);

		int j = 0;
		if (is_regression || (s.outcome != 0)) j=1;
		cnt[j]++;
		if (_params->single_sample_per_pid) {
			if (s_pids[j].find(s.id) == s_pids[j].end()) { v_pids[j].push_back(s.id); s_pids[j].insert(s.id); }
			if (pid2indexes.find(s.id) == pid2indexes.end()) pid2indexes[s.id] = vector<int>();
			pid2indexes[s.id].push_back(i);
		}
		else {
			//if (s_indexes[j].find(s.id) == s_indexes[j].end()) { v_indexes[j].push_back(i); s_indexes[j].insert(i); }
			v_indexes[j].push_back(i);
		}
	}
	if (_params->verbosity > 0) MLOG("Tree %d : bagging : cnt %d %d : s_pids %d %d : v_pids %d %d : v_indexes %d %d\n", id, cnt[0], cnt[1], s_pids[0].size(), s_pids[1].size(), v_pids[0].size(), v_pids[1].size(), v_indexes[0].size(), v_indexes[1].size());

	if (tree_type != TQRF_TREE_REGRESSION) {

		// calculate the bagging probabilities for 0 and 1
		int n0 = 0, n1 = 0;
		float p0 = 0, p1 = 0;

		if (_params->bag_ratio < 0) {
			p0 = _params->bag_prob;
			p1 = _params->bag_prob;
		}
		else {
			if (_params->single_sample_per_pid) {
				n0 = (int)s_pids[0].size();
				n1 = (int)s_pids[1].size();
			}
			else {
				n0 = (int)v_indexes[0].size();
				n1 = (int)v_indexes[1].size();
			}

			//
			// following calculations use the fact that we want : bag_ratio = (n0*p0)/(n1*p1)
			//
			if (n0 > n1) {
				p1 = _params->bag_prob;
				p0 = _params->bag_ratio * (float)(n1+1) * p1 / (float)(n0+1);
			}
			else {
				p0 = _params->bag_prob;
				p1 = p0 * (float)(n0+1) / (_params->bag_ratio * (float)(n1+1));
			}

			p0 = min(p0, (float)1);
			p1 = min(p1, (float)1);
		}

		indexes.clear();

		bag_chooser(_params->bag_with_repeats, _params->single_sample_per_pid, p0, v_pids[0], v_indexes[0], pid2indexes, indexes);
		int nc0 = (int)indexes.size();
		bag_chooser(_params->bag_with_repeats, _params->single_sample_per_pid, p1, v_pids[1], v_indexes[1], pid2indexes, indexes);
		int nc1 = (int)indexes.size() - nc0;
		float actual_ratio = (float)nc0/(nc1+1);
		if (_params->verbosity > 0) MLOG("Tree %d :: bagging : categorial : n0 %d n1 %d p0 %f p1 %f : actual 0: %d 1: %d ratio %f\n", id, n0, n1, p0, p1, nc0, nc1, actual_ratio);

	}
	else {
		// in regression all samples are with j==1
		bag_chooser(_params->bag_with_repeats, _params->single_sample_per_pid, _params->bag_prob, v_pids[1], v_indexes[1], pid2indexes, indexes);
		if (_params->verbosity > 0) MLOG("Tree %d :: bagging : regression : indexes %d\n", id, indexes.size());
	}

	return 0;
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

	return 0;
}

//--------------------------------------------------------------------------------------------------------------------
// simple initiations that needs to be done for the first node into the tree, and in general just right before
// we start the learning process.
//--------------------------------------------------------------------------------------------------------------------
int TQRF_Tree::init_root_node()
{
	nodes.clear();
	TQRF_Node root;

	root.from_idx = 0;
	root.to_idx = (int)indexes.size()-1;

	root.state = TQRF_Node_State_Initiated;

	nodes.push_back(root);

}

//--------------------------------------------------------------------------------------------------------------------
int TQRF_Tree::Train()
{

	// creating the bag
	get_bagged_indexes();


	return 0;
}