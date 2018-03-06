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
		else if (f.first == "qpoints_per_split") qpoints_per_split = stoi(f.second);

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
//
// Our basic convention is that missing value for each feature is mapped to q=0.
//------------------------------------------------------------------------------------------------
int Quantized_Feat::quantize_feat(int i_feat, TQRF_Params &params)
{
	vector<pair<float,int>> curr_data;
	vector<float> *data;
	data = orig_data[i_feat];

	//if (params.verbosity > 0) MLOG("Quantized_Feat::quantize_feat :: feat %d %s :: nvals %d\n", i_feat, feat_names[i_feat].c_str(), data->size());

	curr_data.resize(data->size());
	for (int i=0; i<(int)data->size(); i++) {
		if ((*data)[i] != params.missing_val)
		curr_data[i] ={ (*data)[i], i };
	}


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
	int n = 1; // 0 is kept for missing values
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
	if (v_counts.size() < params.max_q) {
		q_to_val[i_feat] = { params.missing_val }; // 0 place is special : always there to signal us the missing val cases
		for (auto &v : v_counts)
			q_to_val[i_feat].push_back(v.first);
	}
	else {
		// need more work to be done, we have too many different values... 
		// our algorithm is simply trying to make the i-th qval such that there are (i+1)/max_q * len , of the elements
		// up to some fixes that are caused due to integer numbers
		int delta = (int)((float)curr_data.size()/((float)params.max_q-1));
		int j = 0;
		int len = (int)curr_data.size();
		q_to_val[i_feat] ={ params.missing_val };
		while (j < len-1) {
			int prev_j = j;
			j = j + delta;
			if (j >= len) j = len - 1;
			float q = curr_data[j].first;
			while ((j < len-1) && (curr_data[j+1].first == q)) j++;
			if (params.verbosity > 2) MLOG("%s : %f : j %d size %d\n", feat_names[i_feat].c_str(), q, j, j-prev_j);
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
	int q_i = 1; // missing values will be filled later
	int q_size = (int)q_to_val[i_feat].size();
	qx[i_feat].resize(data->size());
	for (auto &v : curr_data) {
		while (q_i < q_size && v.first > q_to_val[i_feat][q_i])
			q_i++;
		qx[i_feat][v.second] = (short)q_i;
	}

	// fill in missing vals
	for (int i=0; i<(int)data->size(); i++)
		if ((*data)[i] == params.missing_val)
			qx[i_feat][i] = 0;

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
		if (params.verbosity > 2) MLOG("Quantized_Feat:: %s :: %d elements :: %f %f %f %f ....\n", feat_names.back().c_str(), orig_data.back()->size(),
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

	medf.get_feature_names(feature_names);

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
			MLOG("Slice %d : %d : 0: %d 1: %d p: %f\n", j, params.time_slices[j], cnt[0][j], cnt[1][j], (float)cnt[1][j]/(1+cnt[1][j]+cnt[0][j]));
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
void TQRF_Forest::init_tables(Quantized_Feat &qfeat)
{
	if (params.tree_type_i == TQRF_TREE_ENTROPY) {

		// initializing a precomputed n * log(n) table for the sake of computing entropy scores
		int max_samples = (int)qfeat.y.size()+1;
		params.log_table.resize(max_samples);
		params.log_table[0] = 0; // we will always use it in an nlogn manner hence the 0 rather than -inf
		for (int i = 1; i < max_samples; i++)
			params.log_table[i] = (double)i * log((double)i)/log(2.0);

	}
}

//------------------------------------------------------------------------------------------------
int TQRF_Forest::Train(MedFeatures &medf, const MedMat<float> &Y)
{
	MLOG("====================================TQRF==================================\n");
	MLOG("TQRF_Forest: Running with params: %s\n", params.init_string.c_str());
	MLOG("TQRF_Forest: Train: medf : %d x %d \n", medf.data.size(), medf.samples.size());

	MedTimer timer;

	timer.start();
	// first - quantifying data
	Quantized_Feat qfeat;
	qfeat.init(medf, params);

	// additional initializations of needed lookup tables (which will be kept in params)
	// for example: log tables for entropy scores etc...
	init_tables(qfeat);

	timer.take_curr_time();
	MLOG("TQRF_Forest: Init qfeat and tables time: %f sec\n", timer.diff_sec());

	MLOG("TQRF_Forest: Starting run on %d trees\n", params.ntrees);

	timer.start();
	// creating the trees and parallelize train on each tree
	trees.resize(params.ntrees);
#pragma omp parallel for
	for (int i=0; i<params.ntrees; i++) {
		trees[i].id = i;
		trees[i].tree_type = params.tree_type_i;
		trees[i].init(qfeat, params);
		trees[i].Train();
		if (params.verbosity > 0) MLOG("TQRF: Trained Tree %d : type %d : indexes %d : feats %d : n_nodes %d\n", trees[i].id, trees[i].tree_type, trees[i].indexes.size(), trees[i].i_feats.size(), trees[i].nodes.size());
	}

	timer.take_curr_time();
	MLOG("TQRF_Forest: Trees training time : %f sec \n", timer.diff_sec());

	return 0;
}

//================================================================================================
// TQRF_Split_Stat family of classes
//================================================================================================
TQRF_Split_Stat *TQRF_Split_Stat::make_tqrf_split_stat(int tree_type)
{
	if (tree_type == TQRF_TREE_LOGRANK)
		return new TQRF_Split_LogRank;
	if (tree_type == TQRF_TREE_ENTROPY)
		return new TQRF_Split_Entropy;
	if (tree_type == TQRF_TREE_REGRESSION)
		return new TQRF_Split_Regression;
	 
	return NULL;
}

//--------------------------------------------------------------------------------------------------------------------
// another source for randomization and speedup is calculating the split scores not at all quantized points
// but only on a subset of them.
// This method selects a random set of qpoints to work on, and returns them:
// an empty qpoints, means : test all points
// otherwise it will contain the q values sorted.
//--------------------------------------------------------------------------------------------------------------------
int TQRF_Split_Stat::get_q_test_points(int feat_max_q, TQRF_Params &params, vector<int> &_qpoints)
{
	int n = feat_max_q - 2; // possible split points are 1,2,...,feat_max_q-2 : 
							// 0 is not a split point: it is missing values , and feat_max_q-1 means we take it all to the left side, so there's no split


	if (params.qpoints_per_split == 0 || params.qpoints_per_split >= n) {
		_qpoints.clear();
		return 0;
	}

	_qpoints.resize(n);
	
	for (int i=0; i<n; i++)
		_qpoints[i] = i+1;

	for (int i=0; i<params.qpoints_per_split; i++)
		swap(_qpoints[i], _qpoints[i + rand_N(n-i)]);

	_qpoints.resize(params.qpoints_per_split);
	sort(_qpoints.begin(), _qpoints.end());
	return 0;
}

//--------------------------------------------------------------------------------------------------------------------
int TQRF_Split_Categorial::init(Quantized_Feat &qf, TQRF_Params &params)
{
	// initializing counts[t][q][c] , sums[t][c]
	ncateg = qf.ncateg;
	maxq = params.max_q;
	nslices = qf.n_time_slices;

	counts.resize(nslices, vector<vector<int>>(maxq , vector<int>(ncateg)));
	sums.resize(nslices, vector<int>(ncateg));

	return 0;
}

//--------------------------------------------------------------------------------------------------------------------
// preparing counts for the current i_feat and node
int TQRF_Split_Categorial::prep_histograms(int i_feat, TQRF_Node &node, vector<int> &indexes, Quantized_Feat &qf, TQRF_Params &params)
{
	int feat_maxq = (int)qf.q_to_val[i_feat].size();

	// zero counts (is it needed at all??)
	for (int t=0; t<nslices; t++)
		for (int q=0; q<feat_maxq; q++)
			for (int c=0; c<ncateg; c++)
				counts[t][q][c] = 0;
	for (int t=0; t<nslices; t++)
		for (int c=0; c<ncateg; c++)
			sums[t][c] = 0;

	// now going over each sample in the node (sometimes we have to sample)
	// and adding counts for each category.
	
	float p_choice = 1;
	
	if (params.max_node_test_samples>0 && params.max_node_test_samples < node.size())
		p_choice = (float)params.max_node_test_samples/(float)node.size();

	vector<short> &feat_qvals = qf.qx[i_feat];
	for (int i=node.from_idx; i<=node.to_idx; i++) {
		int idx = indexes[i];
		int q = feat_qvals[idx];
		int c = (int)(qf.y[idx]);
		int t = qf.last_time_slice[idx];

		counts[t][q][c]++;
		// in case c is 0 we need now to also add a ++ for all cells [0...(t-1)][q][0] . However, we will do those in one swip
		// at the end which should be more efficient.
		if (q > 0) sums[t][c]++; // not summing q=0 cases: missing values are separate
	}

	// now we consider the cases of fewer qpoints 

	// get qpoints
	get_q_test_points(feat_maxq, params, qpoints);

	counts_q = feat_maxq;
	if (qpoints.size() != 0) {

		// checking fewer q points - given orderes in qpoints
		// we squeeze all q values to qpoints.size() values
		int prev = 0;
		for (int i=0; i<qpoints.size(); i++) {
			if (prev > i) {
				for (int t=0; t<nslices; t++)
					for (int c=0; c<ncateg; c++)
						counts[t][i][c] = counts[t][prev][c];
			}

			for (int j=prev+1; j<qpoints[i]; j++) {
				for (int t=0; t<nslices; t++)
					for (int c=0; c<ncateg; c++)
						counts[t][i][c] += counts[t][j][c];
				prev = qpoints[i];
			}
		}
		counts_q = (int)qpoints.size() + 1; // +1 for q=0
	}


	// summing all 0 counts in one swip
	for (int t=nslices-2; t>=0; t--) {
		sums[t][0] += sums[t+1][0];
		for (int q=1; q<counts_q; q++)
			counts[t][q][0] += counts[t+1][q][0];
	}

	return 0;
}

//--------------------------------------------------------------------------------------------------------------------
void TQRF_Split_Categorial::print_histograms()
{
	MLOG("nslices %d ncateg %d maxq %d counts_q %d\n", nslices, ncateg, maxq, counts_q);

	for (int q=0; q<counts_q; q++) {
		MLOG("counts q[%d] :", q);
		for (int c=0; c<ncateg; c++)
			for (int t=0; t<nslices; t++)
				MLOG(" c[%1d] t[%1d] %d :", c, t, counts[t][q][c]);
		MLOG("\n");
	}


	MLOG("sums :");
	for (int c=0; c<ncateg; c++)
		for (int t=0; t<nslices; t++)
			MLOG(" c[%1d] t[%1d] %d :", c, t, sums[t][c]);
	MLOG("\n");
}


//--------------------------------------------------------------------------------------------------------------------
int TQRF_Split_Entropy::get_best_split(TQRF_Params &params, int &best_q, double &best_score)
{
	// the scenario is that we have counts and sums ready and squeezed to the qpoints we want to test
	// we need to go over each of the possible split points, and if it is valid in the sense of 
	// number of samples left on each side (we always test this WITHOUT the missing vals) we get the score

	best_q = -1;			// signing no split point found
	best_score = -1;

	vector<int> left_sums(ncateg, 0), right_sums(ncateg, 0);
	
	// below is a version for a single time slice only... to begin with debugging and build code
	// will be changed to sum over slices later
	//if (nslices != 1)
	//	MTHROW_AND_ERR("ERROR: Running on debug version that supports only a single time slice !!!!!!\n");

	if (params.verbosity > 2) MLOG("TQRF_Split_Entropy::get_best_split counts_q=%d ncateg=%d\n", counts_q, ncateg);

	for (int q=1; q<counts_q-1; q++) {

		int lsum = 0, rsum = 0;
		for (int c=0; c<ncateg; c++) {
			left_sums[c] += counts[0][q][c];
			right_sums[c] = sums[0][c] - left_sums[c];
			lsum += left_sums[c];
			rsum += right_sums[c];
		}

		if (params.verbosity > 2) {
			MLOG("Left  q %d : lsum %d :", q, lsum);
			for (int c=0; c<ncateg; c++)
				MLOG("c[%d] %d :", c, left_sums[c]);
			MLOG("\n");
			MLOG("Right q %d : rsum %d :", q, rsum);
			for (int c=0; c<ncateg; c++)
				MLOG("c[%d] %d :", c, right_sums[c]);
			MLOG("\n");
		}

		if (lsum >= params.min_node && rsum >= params.min_node) {
			double H = 0;

			// add left and right side entropy
			for (int c=0; c<ncateg; c++) {
				H += params.log_table[left_sums[c]];
				H += params.log_table[right_sums[c]];
				H -= params.log_table[left_sums[c]+right_sums[c]];
			}

			// subtract overall sum entropy
			H -= params.log_table[lsum];
			H -= params.log_table[rsum];

			H += params.log_table[lsum+rsum];

			H /= (double)(lsum+rsum)/1000.0; // actual information gain (x1000 for better resolution)

			if (params.verbosity > 2) {
				MLOG("      q %d : H %f best_score %f best_q %d\n", q, H, best_score, best_q);
			}

			//if (best_score < 0 || (H > 0 && H < best_score)) {
			if (H > best_score) {
				best_score = H;
				best_q = q;
			}
		}

	}

	if (qpoints.size() > 0 && best_q>0) best_q = qpoints[best_q-1];

	return 0;
}

//================================================================================================
// TQRF_Tree
//================================================================================================
//--------------------------------------------------------------------------------------------------------------------
// major stage in algorithm:
// we finished the work on deciding if and how to split our node, and need to actually do it.
// list of issues handled in this stage:
// (1) close work on our node, wheather needed a split or not.
// (2) add info from the splitting tqs into our node (distributions etc)
// (3) split node if needed, create the new nodes.
// (4) update indexes as needed
// (5) decide what to do with missing values !!
int TQRF_Tree::node_splitter(int i_curr_node, int i_best, int q_best)
{
	// 
	// finish the work on current node
	//

	TQRF_Node *cnode = &nodes[i_curr_node];

	if (_params->verbosity > 1) MLOG("TQRF: node_splitter : Tree %d : node %d / %d : %d - %d : i_best %d : q_best %d : feat %s\n", 
		id, i_curr_node, nodes.size(), cnode->from_idx, cnode->to_idx, i_best, q_best, (i_best<0)? "NONE" :_qfeat->feature_names[i_best].c_str());
	if (i_best >= 0) {

		// We found a point to split the node
		TQRF_Node Left, Right;

		cnode->i_feat = i_best;
		cnode->bound = _qfeat->q_to_val[i_best][q_best];

		// need to calc node sizes, and general average in order to decide for missing value strategy
		int n_missing = 0, n_left = 0, n_right = 0;
		float sum_vals = 0;
		for (int i=cnode->from_idx; i<=cnode->to_idx; i++) {
			int idx = indexes[i];
			int q = _qfeat->qx[i_best][idx];
			if (q > 0) {
				if (q <= q_best)
					n_left++;
				else
					n_right++;
				sum_vals += _qfeat->q_to_val[i_best][q];
			}
			else
				n_missing++;
		}
		if (_params->verbosity > 1) MLOG("TQRF: node_splitter : Tree %d : node %d : n_missing %d n_left %d n_right %d\n",
											id, i_curr_node, n_missing, n_left, n_right);

		// decide missing val strategy
		if (_params->missing_method == TQRF_MISSING_VALUE_LEFT) cnode->missing_direction = TQRF_MISSING_DIRECTION_LEFT;
		else if (_params->missing_method == TQRF_MISSING_VALUE_RAND_ALL) {
			if (rand_1() < 0.5)
				cnode->missing_direction = TQRF_MISSING_DIRECTION_LEFT;
			else
				cnode->missing_direction = TQRF_MISSING_DIRECTION_RIGHT;
		}
		else if (_params->missing_method == TQRF_MISSING_VALUE_LARGER_NODE || _params->missing_method == TQRF_MISSING_VALUE_MEDIAN) {
			if (n_left >= n_right)
				cnode->missing_direction = TQRF_MISSING_DIRECTION_LEFT;
			else
				cnode->missing_direction = TQRF_MISSING_DIRECTION_RIGHT;
		}
		else if (_params->missing_method == TQRF_MISSING_VALUE_MEAN) {
			float node_avg = sum_vals / ((float)(n_right + n_left) + (float)1e-3);
			if (node_avg <= cnode->bound)
				cnode->missing_direction = TQRF_MISSING_DIRECTION_LEFT;
			else
				cnode->missing_direction = TQRF_MISSING_DIRECTION_RIGHT;
		}
		else if (_params->missing_method == TQRF_MISSING_VALUE_RAND_EACH_SAMPLE)
			cnode->missing_direction = TQRF_MISSING_DIRECTION_RAND_EACH_SAMPLE;


		if (_params->verbosity > 1) MLOG("TQRF: node_splitter : Tree %d : node %d : missing direction %d\n", id, cnode->node_idx, cnode->missing_direction);

		// making the split , first we rearange indexes
		int n_in_left = 0;
		vector<int> left_inds, right_inds;
		//left_inds.reserve(cnode.size());
		//right_inds.reserve(cnode.size());
		for (int i=cnode->from_idx; i<=cnode->to_idx; i++) {
			int idx = indexes[i];
			int q = _qfeat->qx[i_best][idx];
			if (q > 0) {
				if (q <= q_best)
					left_inds.push_back(idx);
				else
					right_inds.push_back(idx);
			}
			else {
				if (cnode->missing_direction == TQRF_MISSING_DIRECTION_LEFT) left_inds.push_back(idx);
				else if (cnode->missing_direction == TQRF_MISSING_DIRECTION_RIGHT) right_inds.push_back(idx);
				else { // if (cnode.missing_direction == TQRF_MISSING_DIRECTION_RAND_EACH_SAMPLE) case ...
					if (rand_1() < (float)0.5)
						left_inds.push_back(idx);
					else
						right_inds.push_back(idx);
				}
			}
		}

		int curr_i = cnode->from_idx;
		for (auto idx : left_inds) indexes[curr_i++] = idx;
		for (auto idx : right_inds) indexes[curr_i++] = idx;
		Left.from_idx = cnode->from_idx;
		Left.to_idx = cnode->from_idx + (int)left_inds.size() - 1;
		Right.from_idx = Left.to_idx + 1;
		Right.to_idx = cnode->to_idx;

		Left.depth = cnode->depth + 1;
		Right.depth = cnode->depth + 1;

		// we may be running in threads over nodes... hence we make sure the following part is protected
#pragma omp critical
		{
			int n_nodes = (int)nodes.size();

			Left.node_idx = n_nodes;
			Right.node_idx = n_nodes+1;

			cnode->left_node = n_nodes;
			cnode->right_node = n_nodes+1;
			cnode->is_terminal = 0;

			nodes.push_back(Left);
			nodes.push_back(Right);

			cnode = &nodes[i_curr_node]; // reassigning as it may have changed in the push !!!
			if (_params->verbosity > 1) MLOG("TQRF: Tree %d Node %d ( s %d d %d ) : split: feat %d : q %d : qval %f : left %d (%d) , right %d (%d)\n", 
				id, cnode->node_idx, cnode->size(), cnode->depth, i_best, q_best, cnode->bound, Left.node_idx, Left.size(), Right.node_idx, Right.size());
		}
	}

	// we now have to finalize the work on cnode (current node) no matter if it was split or not.
	// we need to make sure it has the needed counts etc in order to be able to give predictions
	// plus we change its state

	cnode = &nodes[i_curr_node]; // reassigning as it may have changed in the push !!!
	if (tree_type == TQRF_TREE_ENTROPY || tree_type == TQRF_TREE_LOGRANK) {

		assert(cnode->node_idx == nodes[i_curr_node].node_idx);
		// we go over the y values for the indexes in our node and collect them for each time slice
		cnode->time_categ_count.resize(_qfeat->n_time_slices, vector<int>(_qfeat->ncateg, 0));
		for (int i=cnode->from_idx; i<cnode->to_idx; i++) {
			int idx = indexes[i];
			int c = (int)_qfeat->y[idx];
			int t = _qfeat->last_time_slice[idx];
			cnode->time_categ_count[t][c]++;
		}
		// reverse adding 0 (=control) categs
		for (int t=_qfeat->n_time_slices-2; t>=0; t--)
			cnode->time_categ_count[t][0] += cnode->time_categ_count[t+1][0];
	}
	else {
		MTHROW_AND_ERR("TQRF::node_splitter(): tree_type %d doesn't know how to finalize nodes yet !!... sayonara...\n", tree_type);
	}
	if (cnode->depth >= _params->max_depth || cnode->size() <= _params->min_node) {
		nodes[i_curr_node].is_terminal = 1;
		if (_params->verbosity > 1) MLOG("TQRF: node_splitter : Tree %d : node %d : depth %d size %d : terminal %d\n", id, i_curr_node, nodes[i_curr_node].depth, nodes[i_curr_node].size(), nodes[i_curr_node].is_terminal);
	}

	cnode->state = TQRF_Node_State_Done;

	if (_params->verbosity > 1) MLOG("TQRF: node_splitter : Tree %d : node %d : Done\n", id, i_curr_node);

	return 0;
}
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


	// now choosing the features to be used in this tree (in case feature bagging is requested)
	i_feats.clear();
	for (int i=0; i<_qfeat->nfeat; i++)
		if (rand_1() <= _params->bag_feat)
			i_feats.push_back(i);

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
	root.depth = 0;

	root.node_idx = 0;
	nodes.push_back(root);

	return 0;
}

//--------------------------------------------------------------------------------------------------------------------
// returns the next node with state TQRF_Node_State_Initiated
int TQRF_Tree::get_next_node(int curr_node)
{
	for (int i=max(curr_node+1, 0); i<nodes.size(); i++)
		if (nodes[i].state == TQRF_Node_State_Initiated) {
			nodes[i].state = TQRF_Node_State_In_Progress;
			return i;
		}
	return -1; // none found
}

//--------------------------------------------------------------------------------------------------------------------
// getting a list of features to test, based on tree parameters
int TQRF_Tree::get_feats_to_test(vector<int> &feats_to_test)
{
	// first we need to know how many features we need to choose
	int n_to_choose = max(_params->ntry, (int)(_params->ntry_prob * (float)i_feats.size()));

	if (n_to_choose > (int)i_feats.size())
		n_to_choose = (int)i_feats.size();

	// we now go through n_to_choose steps of random swapping on i_feats
	feats_to_test.clear();

	if (_params->verbosity > 1) MLOG("TQRF: get_feats_to_test: n_to_choose %d/%d\n", n_to_choose, i_feats.size());
	for (int i=0; i<n_to_choose; i++) {
		int n = (int)i_feats.size() - i - 1;
		int j = rand_N(n) + i;
		int f = i_feats[i];
		//MLOG("i %d n %d j %d f %d\n", i, n, j, f);
		i_feats[i] = i_feats[j];
		i_feats[j] = f;
		feats_to_test.push_back(i_feats[i]);
	}
	if (_params->verbosity > 1) {
		MLOG("TQRF: Tree %d : chose %d features to split by : ", id, feats_to_test.size());
		for (auto f : feats_to_test) MLOG(" %d", f);
		MLOG("\n");
	}

	return 0;
}

//--------------------------------------------------------------------------------------------------------------------
void TQRF_Tree::free_split_stats(vector<TQRF_Split_Stat *> &tqs)
{
	for (auto &t : tqs) {
		delete t;
	}

	tqs.clear();
}

//--------------------------------------------------------------------------------------------------------------------
int TQRF_Tree::init_split_stats(vector<TQRF_Split_Stat *> &tqs)
{
	free_split_stats(tqs);
	tqs.resize(i_feats.size());

	for (int i=0; i<tqs.size(); i++) {
		tqs[i] = TQRF_Split_Stat::make_tqrf_split_stat(tree_type);
		if (tqs[i] == NULL)
			return -1;
		tqs[i]->init((*_qfeat), (*_params));
	}


	return 0;
}




//--------------------------------------------------------------------------------------------------------------------
// non threaded version
// threading a specific tree is much more complex.... we have to do it over nodes.
//--------------------------------------------------------------------------------------------------------------------
int TQRF_Tree::Train()
{
	if (_params->verbosity > 1) MLOG("Tree %d  Train() before bagging\n", id);

	// creating the bag
	get_bagged_indexes();

	if (_params->verbosity > 1) MLOG("After get_bagged_indexes()\n");
	// initializing tree and root
	init_root_node();
	if (_params->verbosity > 1) MLOG("After init_root_node()\n");

	//bool go_on = true;
	int i_curr_node = -1;

	vector<int> feats_to_test;
	vector<TQRF_Split_Stat *> tqs;

	init_split_stats(tqs);
	if (_params->verbosity > 1) MLOG("After init_split_stats()\n");

	vector<pair<int, double>> best_q;
	while (1) {

		if ((i_curr_node = get_next_node(i_curr_node)) < 0)
			break; // finished work on this tree - no more nodes to work on.

		if (_params->verbosity > 1) MLOG("TQRF: Tree %d Working on =================>>>>> node %d\n", id, i_curr_node);

		// getting a list of features to test, based on tree parameters
		get_feats_to_test(feats_to_test);

		if (feats_to_test.size() > 0)
			best_q.resize(feats_to_test.size()); //, { -1,-1 });

		// optional "easy" threading right here !
		for (int i=0; i<feats_to_test.size(); i++) {
			best_q[i] ={ -1, -1.0 };
			int i_f = feats_to_test[i];
			if (_params->verbosity > 2) MLOG("TQRF: Tree %d node %d feat[%d] = %d : %s before histogram\n", id, i_curr_node, i, i_f, _qfeat->feature_names[i_f].c_str());
			tqs[i]->prep_histograms(i_f, nodes[i_curr_node], indexes, (*_qfeat), (*_params));
			if (_params->verbosity > 2) MLOG("TQRF: Tree %d node %d feat[%d] = %d : %s after histogram\n", id, i_curr_node, i, i_f, _qfeat->feature_names[i_f].c_str());
			if (_params->verbosity > 2) tqs[i]->print_histograms();
			if (nodes[i_curr_node].depth <= _params->max_depth)
				tqs[i]->get_best_split((*_params), best_q[i].first, best_q[i].second);
			if (_params->verbosity > 2) MLOG("TQRF: Tree %d node %d feat[%d] = %d : after get_best_split %s : %d %f : cut off val %f\n", id, i_curr_node, i, i_f, _qfeat->feature_names[i_f].c_str(), best_q[i].first, best_q[i].second, (best_q[i].first < 0) ? -1 : _qfeat->q_to_val[i_f][best_q[i].first]);
		}

		// choose best choice : scores are ALWAYS for maximum
		int i_best = -1;
		int q_best = -1;
		double q_best_score = -1e10;

		for (int i=0; i<feats_to_test.size(); i++) {
			if (_params->verbosity > 2) MLOG("TQRF: after features scan : %d : feat %d %s : q %d score %f\n", i, feats_to_test[i], _qfeat->feature_names[feats_to_test[i]].c_str(), best_q[i].first, best_q[i].second);
			if (best_q[i].first > 0 && best_q[i].second > q_best_score) {
				q_best_score = best_q[i].second;
				q_best = best_q[i].first;
				i_best = feats_to_test[i];
			}
		}

		if (_params->verbosity > 1) MLOG("TQRF: Tree %d Node %d : best feature %d %s , q %d qval %f score %f\n", id, nodes[i_curr_node].node_idx, i_best, (i_best < 0) ? "" : _qfeat->feature_names[i_best].c_str(), q_best, (i_best < 0) ? 0 : _qfeat->q_to_val[i_best][q_best], q_best_score);

		node_splitter(i_curr_node, i_best, q_best);


	}

	// need to free split stats
	free_split_stats(tqs);

	return 0;
}
