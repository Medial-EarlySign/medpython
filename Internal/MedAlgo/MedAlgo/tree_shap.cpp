#include "tree_shap.h"
#include <omp.h>

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL LOG_DEF_LEVEL

bool comp_score_flt_str(const pair<float, string> &pr1, const pair<float, string> &pr2) {
	return abs(pr1.first) > abs(pr2.first); //bigger is better
}

TreeEnsemble::TreeEnsemble() {
	is_allocate = false;
}

TreeEnsemble::TreeEnsemble(int *children_left, int *children_right, int *children_default, int *features,
	tfloat *thresholds, tfloat *values, tfloat *node_sample_weights,
	unsigned max_depth, unsigned tree_limit, tfloat base_offset,
	unsigned max_nodes, unsigned num_outputs) :
	children_left(children_left), children_right(children_right),
	children_default(children_default), features(features), thresholds(thresholds),
	values(values), node_sample_weights(node_sample_weights),
	max_depth(max_depth), tree_limit(tree_limit),
	base_offset(base_offset), max_nodes(max_nodes), num_outputs(num_outputs), is_allocate(true) {}

void TreeEnsemble::get_tree(TreeEnsemble &tree, const unsigned i) const {
	const unsigned d = i * max_nodes;

	tree.children_left = children_left + d;
	tree.children_right = children_right + d;
	tree.children_default = children_default + d;
	tree.features = features + d;
	tree.thresholds = thresholds + d;
	tree.values = values + d * num_outputs;
	tree.node_sample_weights = node_sample_weights + d;
	tree.max_depth = max_depth;
	tree.tree_limit = 1;
	tree.base_offset = base_offset;
	tree.max_nodes = max_nodes;
	tree.num_outputs = num_outputs;
	tree.is_allocate = true;
}

void TreeEnsemble::allocate(unsigned tree_limit_in, unsigned max_nodes_in, unsigned num_outputs_in) {
	tree_limit = tree_limit_in;
	max_nodes = max_nodes_in;
	num_outputs = num_outputs_in;
	children_left = new int[tree_limit * max_nodes];
	children_right = new int[tree_limit * max_nodes];
	children_default = new int[tree_limit * max_nodes];
	features = new int[tree_limit * max_nodes];
	thresholds = new tfloat[tree_limit * max_nodes];
	values = new tfloat[tree_limit * max_nodes * num_outputs];
	node_sample_weights = new tfloat[tree_limit * max_nodes];
	is_allocate = true;
}

void TreeEnsemble::free() {
	if (is_allocate) {
		delete[] children_left;
		delete[] children_right;
		delete[] children_default;
		delete[] features;
		delete[] thresholds;
		delete[] values;
		delete[] node_sample_weights;
		children_left = NULL;
		children_right = NULL;
		children_default = NULL;
		features = NULL;
		thresholds = NULL;
		values = NULL;
		node_sample_weights = NULL;
	}
	is_allocate = false;
}

ExplanationDataset::ExplanationDataset() {}

ExplanationDataset::ExplanationDataset(tfloat *X, bool *X_missing, tfloat *y, tfloat *R, bool *R_missing, unsigned num_X,
	unsigned M, unsigned num_R) :
	X(X), X_missing(X_missing), y(y), R(R), R_missing(R_missing), num_X(num_X), M(M), num_R(num_R) {}

void ExplanationDataset::get_x_instance(ExplanationDataset &instance, const unsigned i) const {
	instance.M = M;
	instance.X = X + i * M;
	instance.X_missing = X_missing + i * M;
	instance.num_X = 1;
}

PathElement::PathElement() {}
PathElement::PathElement(int i, tfloat z, tfloat o, tfloat w) :
	feature_index(i), zero_fraction(z), one_fraction(o), pweight(w) {}


inline tfloat logistic_transform(const tfloat margin, const tfloat y) {
	return 1 / (1 + exp(-margin));
}

inline tfloat logistic_nlogloss_transform(const tfloat margin, const tfloat y) {
	return log(1 + exp(margin)) - y * margin; // y is in {0, 1}
}

inline tfloat squared_loss_transform(const tfloat margin, const tfloat y) {
	return (margin - y) * (margin - y);
}


inline tfloat *tree_predict(unsigned i, const TreeEnsemble &trees, const tfloat *x, const bool *x_missing) {
	const unsigned offset = i * trees.max_nodes;
	unsigned node = 0;
	while (true) {
		const unsigned pos = offset + node;
		const unsigned feature = trees.features[pos];

		// we hit a leaf so return a pointer to the values
		if (trees.children_left[pos] < 0) {
			return trees.values + pos * trees.num_outputs;
		}

		// otherwise we are at an internal node and need to recurse
		if (x_missing[feature]) {
			node = trees.children_default[pos];
		}
		else if (x[feature] <= trees.thresholds[pos]) {
			node = trees.children_left[pos];
		}
		else {
			node = trees.children_right[pos];
		}
	}
}

inline void dense_tree_predict(tfloat *out, const TreeEnsemble &trees, const ExplanationDataset &data, unsigned model_transform) {
	tfloat *row_out = out;
	const tfloat *x = data.X;
	const bool *x_missing = data.X_missing;

	// see what transform (if any) we have
	tfloat(*transform)(const tfloat margin, const tfloat y) = NULL;
	switch (model_transform) {
	case MODEL_TRANSFORM::logistic:
		transform = logistic_transform;
		break;

	case MODEL_TRANSFORM::logistic_nlogloss:
		transform = logistic_nlogloss_transform;
		break;

	case MODEL_TRANSFORM::squared_loss:
		transform = squared_loss_transform;
		break;
	}

	for (unsigned i = 0; i < data.num_X; ++i) {

		// add the base offset
		for (unsigned k = 0; k < trees.num_outputs; ++k) {
			row_out[k] += trees.base_offset;
		}

		// add the leaf values from each tree
		for (unsigned j = 0; j < trees.tree_limit; ++j) {
			const tfloat *leaf_value = tree_predict(j, trees, x, x_missing);

			for (unsigned k = 0; k < trees.num_outputs; ++k) {
				row_out[k] += leaf_value[k];
			}
		}

		// apply any needed transform
		if (transform != NULL) {
			const tfloat y_i = data.y == NULL ? 0 : data.y[i];
			for (unsigned k = 0; k < trees.num_outputs; ++k) {
				row_out[k] = transform(row_out[k], y_i);
			}
		}

		x += data.M;
		x_missing += data.M;
		row_out += trees.num_outputs;
	}
}

inline void tree_update_weights(unsigned i, TreeEnsemble &trees, const tfloat *x, const bool *x_missing) {
	const unsigned offset = i * trees.max_nodes;
	unsigned node = 0;
	while (true) {
		const unsigned pos = offset + node;
		const unsigned feature = trees.features[pos];

		// Record that a sample passed through this node
		trees.node_sample_weights[pos] += 1.0;

		// we hit a leaf so return a pointer to the values
		if (trees.children_left[pos] < 0) break;

		// otherwise we are at an internal node and need to recurse
		if (x_missing[feature]) {
			node = trees.children_default[pos];
		}
		else if (x[feature] <= trees.thresholds[pos]) {
			node = trees.children_left[pos];
		}
		else {
			node = trees.children_right[pos];
		}
	}
}

inline void dense_tree_update_weights(TreeEnsemble &trees, const ExplanationDataset &data) {
	const tfloat *x = data.X;
	const bool *x_missing = data.X_missing;

	for (unsigned i = 0; i < data.num_X; ++i) {

		// add the leaf values from each tree
		for (unsigned j = 0; j < trees.tree_limit; ++j) {
			tree_update_weights(j, trees, x, x_missing);
		}

		x += data.M;
		x_missing += data.M;
	}
}

inline void tree_saabas(tfloat *out, const TreeEnsemble &tree, const ExplanationDataset &data) {
	unsigned curr_node = 0;
	unsigned next_node = 0;
	while (true) {

		// we hit a leaf and are done
		if (tree.children_left[curr_node] < 0) return;

		// otherwise we are at an internal node and need to recurse
		const unsigned feature = tree.features[curr_node];
		if (data.X_missing[feature]) {
			next_node = tree.children_default[curr_node];
		}
		else if (data.X[feature] <= tree.thresholds[curr_node]) {
			next_node = tree.children_left[curr_node];
		}
		else {
			next_node = tree.children_right[curr_node];
		}

		// assign credit to this feature as the difference in values at the current node vs. the next node
		for (unsigned i = 0; i < tree.num_outputs; ++i) {
			out[feature * tree.num_outputs + i] += tree.values[next_node * tree.num_outputs + i] - tree.values[curr_node * tree.num_outputs + i];
		}

		curr_node = next_node;
	}
}

/**
* This runs Tree SHAP with a per tree path conditional dependence assumption.
*/
void dense_tree_saabas(tfloat *out_contribs, const TreeEnsemble& trees, const ExplanationDataset &data) {
	MedTimer tm;
	tm.start();
	chrono::high_resolution_clock::time_point tm_prog = chrono::high_resolution_clock::now();
	int progress = 0;
	int max_loop = data.num_X;

	// build explanation for each sample
	for (int i = 0; i < data.num_X; ++i) {
		TreeEnsemble tree;
		ExplanationDataset instance;
		tfloat *instance_out_contribs = out_contribs + i * (data.M + 1) * trees.num_outputs;
		data.get_x_instance(instance, i);

		// aggregate the effect of explaining each tree
		// (this works because of the linearity property of Shapley values)
		for (unsigned j = 0; j < trees.tree_limit; ++j) {
			trees.get_tree(tree, j);
			tree_saabas(instance_out_contribs, tree, instance);
		}

		// apply the base offset to the bias term
		for (unsigned j = 0; j < trees.num_outputs; ++j) {
			instance_out_contribs[data.M * trees.num_outputs + j] += trees.base_offset;
		}

#pragma omp atomic
		++progress;
		double duration = (unsigned long long)(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
			- tm_prog).count()) / 1000000.0;
		if (duration > 15 && progress % 50 == 0) {
#pragma omp critical
			tm_prog = chrono::high_resolution_clock::now();
			double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
				- tm.t[0]).count()) / 1000000.0;
			double estimate_time = int(double(max_loop - progress) / double(progress) * double(time_elapsed));
			MLOG("SHAPLEY Processed %d out of %d(%2.2f%%) time elapsed: %2.1f Minutes, "
				"estimate time to finish %2.1f Minutes\n", progress, max_loop, 100.0*(progress / float(max_loop)), time_elapsed / 60,
				estimate_time / 60.0);
		}
	}
}


// extend our decision path with a fraction of one and zero extensions
inline void extend_path(PathElement *unique_path, unsigned unique_depth,
	tfloat zero_fraction, tfloat one_fraction, int feature_index) {
	unique_path[unique_depth].feature_index = feature_index;
	unique_path[unique_depth].zero_fraction = zero_fraction;
	unique_path[unique_depth].one_fraction = one_fraction;
	unique_path[unique_depth].pweight = (unique_depth == 0 ? 1.0f : 0.0f);
	for (int i = unique_depth - 1; i >= 0; i--) {
		unique_path[i + 1].pweight += one_fraction * unique_path[i].pweight * (i + 1)
			/ static_cast<tfloat>(unique_depth + 1);
		unique_path[i].pweight = zero_fraction * unique_path[i].pweight * (unique_depth - i)
			/ static_cast<tfloat>(unique_depth + 1);
	}
}

// undo a previous extension of the decision path
inline void unwind_path(PathElement *unique_path, unsigned unique_depth, unsigned path_index) {
	const tfloat one_fraction = unique_path[path_index].one_fraction;
	const tfloat zero_fraction = unique_path[path_index].zero_fraction;
	tfloat next_one_portion = unique_path[unique_depth].pweight;

	for (int i = unique_depth - 1; i >= 0; --i) {
		if (one_fraction != 0) {
			const tfloat tmp = unique_path[i].pweight;
			unique_path[i].pweight = next_one_portion * (unique_depth + 1)
				/ static_cast<tfloat>((i + 1) * one_fraction);
			next_one_portion = tmp - unique_path[i].pweight * zero_fraction * (unique_depth - i)
				/ static_cast<tfloat>(unique_depth + 1);
		}
		else {
			unique_path[i].pweight = (unique_path[i].pweight * (unique_depth + 1))
				/ static_cast<tfloat>(zero_fraction * (unique_depth - i));
		}
	}

	for (unsigned i = path_index; i < unique_depth; ++i) {
		unique_path[i].feature_index = unique_path[i + 1].feature_index;
		unique_path[i].zero_fraction = unique_path[i + 1].zero_fraction;
		unique_path[i].one_fraction = unique_path[i + 1].one_fraction;
	}
}

// determine what the total permuation weight would be if
// we unwound a previous extension in the decision path
inline tfloat unwound_path_sum(const PathElement *unique_path, unsigned unique_depth,
	unsigned path_index) {
	const tfloat one_fraction = unique_path[path_index].one_fraction;
	const tfloat zero_fraction = unique_path[path_index].zero_fraction;
	tfloat next_one_portion = unique_path[unique_depth].pweight;
	tfloat total = 0;

	if (one_fraction != 0) {
		for (int i = unique_depth - 1; i >= 0; --i) {
			const tfloat tmp = next_one_portion / static_cast<tfloat>((i + 1) * one_fraction);
			total += tmp;
			next_one_portion = unique_path[i].pweight - tmp * zero_fraction * (unique_depth - i);
		}
	}
	else {
		for (int i = unique_depth - 1; i >= 0; --i) {
			total += unique_path[i].pweight / (zero_fraction * (unique_depth - i));
		}
	}
	return total * (unique_depth + 1);
}

// recursive computation of SHAP values for a decision tree
inline void tree_shap_recursive(const unsigned num_outputs, const int *children_left,
	const int *children_right,
	const int *children_default, const int *features,
	const tfloat *thresholds, const tfloat *values,
	const tfloat *node_sample_weight,
	const tfloat *x, const bool *x_missing, tfloat *phi,
	unsigned node_index, unsigned unique_depth,
	PathElement *parent_unique_path, tfloat parent_zero_fraction,
	tfloat parent_one_fraction, int parent_feature_index,
	int condition, unsigned condition_feature,
	tfloat condition_fraction) {

	// stop if we have no weight coming down to us
	if (condition_fraction == 0) return;

	// extend the unique path
	PathElement *unique_path = parent_unique_path + unique_depth + 1;
	std::copy(parent_unique_path, parent_unique_path + unique_depth + 1, unique_path);

	if (condition == 0 || condition_feature != static_cast<unsigned>(parent_feature_index)) {
		extend_path(unique_path, unique_depth, parent_zero_fraction,
			parent_one_fraction, parent_feature_index);
	}
	const unsigned split_index = features[node_index];

	// leaf node
	if (children_right[node_index] < 0) {
		for (unsigned i = 1; i <= unique_depth; ++i) {
			const tfloat w = unwound_path_sum(unique_path, unique_depth, i);
			const PathElement &el = unique_path[i];
			const unsigned phi_offset = el.feature_index * num_outputs;
			const unsigned values_offset = node_index * num_outputs;
			const tfloat scale = w * (el.one_fraction - el.zero_fraction) * condition_fraction;
			for (unsigned j = 0; j < num_outputs; ++j) {
				phi[phi_offset + j] += scale * values[values_offset + j];
			}
		}

		// internal node
	}
	else {
		// find which branch is "hot" (meaning x would follow it)
		unsigned hot_index = 0;
		if (x_missing[split_index]) {
			hot_index = children_default[node_index];
		}
		else if (x[split_index] <= thresholds[node_index]) {
			hot_index = children_left[node_index];
		}
		else {
			hot_index = children_right[node_index];
		}
		const unsigned cold_index = (static_cast<int>(hot_index) == children_left[node_index] ?
			children_right[node_index] : children_left[node_index]);
		const tfloat w = node_sample_weight[node_index];
		const tfloat hot_zero_fraction = node_sample_weight[hot_index] / w;
		const tfloat cold_zero_fraction = node_sample_weight[cold_index] / w;
		tfloat incoming_zero_fraction = 1;
		tfloat incoming_one_fraction = 1;

		// see if we have already split on this feature,
		// if so we undo that split so we can redo it for this node
		unsigned path_index = 0;
		for (; path_index <= unique_depth; ++path_index) {
			if (static_cast<unsigned>(unique_path[path_index].feature_index) == split_index) break;
		}
		if (path_index != unique_depth + 1) {
			incoming_zero_fraction = unique_path[path_index].zero_fraction;
			incoming_one_fraction = unique_path[path_index].one_fraction;
			unwind_path(unique_path, unique_depth, path_index);
			unique_depth -= 1;
		}

		// divide up the condition_fraction among the recursive calls
		tfloat hot_condition_fraction = condition_fraction;
		tfloat cold_condition_fraction = condition_fraction;
		if (condition > 0 && split_index == condition_feature) {
			cold_condition_fraction = 0;
			unique_depth -= 1;
		}
		else if (condition < 0 && split_index == condition_feature) {
			hot_condition_fraction *= hot_zero_fraction;
			cold_condition_fraction *= cold_zero_fraction;
			unique_depth -= 1;
		}

		tree_shap_recursive(
			num_outputs, children_left, children_right, children_default, features, thresholds, values,
			node_sample_weight, x, x_missing, phi, hot_index, unique_depth + 1, unique_path,
			hot_zero_fraction * incoming_zero_fraction, incoming_one_fraction,
			split_index, condition, condition_feature, hot_condition_fraction
		);

		tree_shap_recursive(
			num_outputs, children_left, children_right, children_default, features, thresholds, values,
			node_sample_weight, x, x_missing, phi, cold_index, unique_depth + 1, unique_path,
			cold_zero_fraction * incoming_zero_fraction, 0,
			split_index, condition, condition_feature, cold_condition_fraction
		);
	}
}

inline int compute_expectations(TreeEnsemble &tree, int i = 0, int depth = 0) {
	unsigned max_depth = 0;

	if (tree.children_right[i] >= 0) {
		const unsigned li = tree.children_left[i];
		const unsigned ri = tree.children_right[i];
		const unsigned depth_left = compute_expectations(tree, li, depth + 1);
		const unsigned depth_right = compute_expectations(tree, ri, depth + 1);
		const tfloat left_weight = tree.node_sample_weights[li];
		const tfloat right_weight = tree.node_sample_weights[ri];
		const unsigned li_offset = li * tree.num_outputs;
		const unsigned ri_offset = ri * tree.num_outputs;
		const unsigned i_offset = i * tree.num_outputs;
		for (unsigned j = 0; j < tree.num_outputs; ++j) {
			const tfloat v = (left_weight * tree.values[li_offset + j] + right_weight * tree.values[ri_offset + j]) / (left_weight + right_weight);
			tree.values[i_offset + j] = v;
		}
		max_depth = std::max(depth_left, depth_right) + 1;
	}

	if (depth == 0) tree.max_depth = max_depth;

	return max_depth;
}

inline void tree_shap(const TreeEnsemble& tree, const ExplanationDataset &data,
	tfloat *out_contribs, int condition, unsigned condition_feature) {

	// update the reference value with the expected value of the tree's predictions
	if (condition == 0) {
		for (unsigned j = 0; j < tree.num_outputs; ++j) {
			out_contribs[data.M * tree.num_outputs + j] += tree.values[j];
		}
	}

	// Pre-allocate space for the unique path data
	const unsigned maxd = tree.max_depth + 2; // need a bit more space than the max depth
	PathElement *unique_path_data = new PathElement[(maxd * (maxd + 1)) / 2];

	tree_shap_recursive(
		tree.num_outputs, tree.children_left, tree.children_right, tree.children_default,
		tree.features, tree.thresholds, tree.values, tree.node_sample_weights, data.X,
		data.X_missing, out_contribs, 0, 0, unique_path_data, 1, 1, -1, condition,
		condition_feature, 1
	);

	delete[] unique_path_data;
}


unsigned build_merged_tree_recursive(TreeEnsemble &out_tree, const TreeEnsemble &trees,
	const tfloat *data, const bool *data_missing, int *data_inds,
	const unsigned num_background_data_inds, unsigned num_data_inds,
	unsigned M, unsigned row = 0, unsigned i = 0, unsigned pos = 0,
	tfloat *leaf_value = NULL) {
	//tfloat new_leaf_value[trees.num_outputs];
	tfloat *new_leaf_value = (tfloat *)alloca(sizeof(tfloat) * trees.num_outputs); // allocate on the stack
	unsigned row_offset = row * trees.max_nodes;

	// we have hit a terminal leaf!!!
	if (trees.children_left[row_offset + i] < 0 && row + 1 == trees.tree_limit) {

		// create the leaf node
		const tfloat *vals = trees.values + (row * trees.max_nodes + i) * trees.num_outputs;
		if (leaf_value == NULL) {
			for (unsigned j = 0; j < trees.num_outputs; ++j) {
				out_tree.values[pos * trees.num_outputs + j] = vals[j];
			}
		}
		else {
			for (unsigned j = 0; j < trees.num_outputs; ++j) {
				out_tree.values[pos * trees.num_outputs + j] = leaf_value[j] + vals[j];
			}
		}
		out_tree.children_left[pos] = -1;
		out_tree.children_right[pos] = -1;
		out_tree.children_default[pos] = -1;
		out_tree.features[pos] = -1;
		out_tree.thresholds[pos] = 0;
		out_tree.node_sample_weights[pos] = num_background_data_inds;

		return pos;
	}

	// we hit an intermediate leaf (so just add the value to our accumulator and move to the next tree)
	if (trees.children_left[row_offset + i] < 0) {

		// accumulate the value of this original leaf so it will land on all eventual terminal leaves
		const tfloat *vals = trees.values + (row * trees.max_nodes + i) * trees.num_outputs;
		if (leaf_value == NULL) {
			for (unsigned j = 0; j < trees.num_outputs; ++j) {
				new_leaf_value[j] = vals[j];
			}
		}
		else {
			for (unsigned j = 0; j < trees.num_outputs; ++j) {
				new_leaf_value[j] = leaf_value[j] + vals[j];
			}
		}
		leaf_value = new_leaf_value;

		// move forward to the next tree
		row += 1;
		row_offset += trees.max_nodes;
		i = 0;
	}

	// split the data inds by this node's threshold
	const tfloat t = trees.thresholds[row_offset + i];
	const int f = trees.features[row_offset + i];
	const bool right_default = trees.children_default[row_offset + i] == trees.children_right[row_offset + i];
	int low_ptr = 0;
	int high_ptr = num_data_inds - 1;
	unsigned num_left_background_data_inds = 0;
	int low_data_ind;
	while (low_ptr <= high_ptr) {
		low_data_ind = data_inds[low_ptr];
		const int data_ind = std::abs(low_data_ind) * M + f;
		const bool is_missing = data_missing[data_ind];
		if ((!is_missing && data[data_ind] > t) || (right_default && is_missing)) {
			data_inds[low_ptr] = data_inds[high_ptr];
			data_inds[high_ptr] = low_data_ind;
			high_ptr -= 1;
		}
		else {
			if (low_data_ind >= 0) ++num_left_background_data_inds; // negative data_inds are not background samples
			low_ptr += 1;
		}
	}
	int *left_data_inds = data_inds;
	const unsigned num_left_data_inds = low_ptr;
	int *right_data_inds = data_inds + low_ptr;
	const unsigned num_right_data_inds = num_data_inds - num_left_data_inds;
	const unsigned num_right_background_data_inds = num_background_data_inds - num_left_background_data_inds;

	// all the data went right, so we skip creating this node and just recurse right
	if (num_left_data_inds == 0) {
		return build_merged_tree_recursive(
			out_tree, trees, data, data_missing, data_inds,
			num_background_data_inds, num_data_inds, M, row,
			trees.children_right[row_offset + i], pos, leaf_value
		);

		// all the data went left, so we skip creating this node and just recurse left
	}
	else if (num_right_data_inds == 0) {
		return build_merged_tree_recursive(
			out_tree, trees, data, data_missing, data_inds,
			num_background_data_inds, num_data_inds, M, row,
			trees.children_left[row_offset + i], pos, leaf_value
		);

		// data went both ways so we create this node and recurse down both paths
	}
	else {

		// build the left subtree
		const unsigned new_pos = build_merged_tree_recursive(
			out_tree, trees, data, data_missing, left_data_inds,
			num_left_background_data_inds, num_left_data_inds, M, row,
			trees.children_left[row_offset + i], pos + 1, leaf_value
		);

		// fill in the data for this node
		out_tree.children_left[pos] = pos + 1;
		out_tree.children_right[pos] = new_pos + 1;
		if (trees.children_left[row_offset + i] == trees.children_default[row_offset + i]) {
			out_tree.children_default[pos] = pos + 1;
		}
		else {
			out_tree.children_default[pos] = new_pos + 1;
		}

		out_tree.features[pos] = trees.features[row_offset + i];
		out_tree.thresholds[pos] = trees.thresholds[row_offset + i];
		out_tree.node_sample_weights[pos] = num_background_data_inds;

		// build the right subtree
		return build_merged_tree_recursive(
			out_tree, trees, data, data_missing, right_data_inds,
			num_right_background_data_inds, num_right_data_inds, M, row,
			trees.children_right[row_offset + i], new_pos + 1, leaf_value
		);
	}
}


void build_merged_tree(TreeEnsemble &out_tree, const ExplanationDataset &data, const TreeEnsemble &trees) {

	// create a joint data matrix from both X and R matrices
	tfloat *joined_data = new tfloat[(data.num_X + data.num_R) * data.M];
	std::copy(data.X, data.X + data.num_X * data.M, joined_data);
	std::copy(data.R, data.R + data.num_R * data.M, joined_data + data.num_X * data.M);
	bool *joined_data_missing = new bool[(data.num_X + data.num_R) * data.M];
	std::copy(data.X_missing, data.X_missing + data.num_X * data.M, joined_data_missing);
	std::copy(data.R_missing, data.R_missing + data.num_R * data.M, joined_data_missing + data.num_X * data.M);

	// create an starting array of data indexes we will recursively sort
	int *data_inds = new int[data.num_X + data.num_R];
	for (unsigned i = 0; i < data.num_X; ++i) data_inds[i] = i;
	for (unsigned int i = data.num_X; i < data.num_X + data.num_R; ++i) {
		data_inds[i] = -(int)i; // a negative index means it won't be recorded as a background sample
	}

	build_merged_tree_recursive(
		out_tree, trees, joined_data, joined_data_missing, data_inds, data.num_R,
		data.num_X + data.num_R, data.M
	);

	delete[] joined_data;
	delete[] joined_data_missing;
	delete[] data_inds;
}




// https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
inline int bin_coeff(int n, int k) {
	int res = 1;
	if (k > n - k)
		k = n - k;
	for (int i = 0; i < k; ++i) {
		res *= (n - i);
		res /= (i + 1);
	}
	return res;
}

// note this only handles single output models, so multi-output models get explained using multiple passes
inline void tree_shap_indep(const unsigned max_depth, const unsigned num_feats,
	const unsigned num_nodes, const tfloat *x,
	const bool *x_missing, const tfloat *r,
	const bool *r_missing, tfloat *out_contribs,
	float *pos_lst, float *neg_lst, signed short *feat_hist,
	float *memoized_weights, int *node_stack, Node *mytree) {

	//     const bool DEBUG = true;
	//     ofstream myfile;
	//     if (DEBUG) {
	//       myfile.open ("/homes/gws/hughchen/shap/out.txt",fstream::app);
	//       myfile << "Entering tree_shap_indep\n";
	//     }
	int ns_ctr = 0;
	std::fill_n(feat_hist, num_feats, 0);
	short node = 0, feat, cl, cr, cd, pnode, pfeat = -1;
	short next_xnode = -1, next_rnode = -1;
	short next_node = -1, from_child = -1;
	float thres, pos_x = 0, neg_x = 0, pos_r = 0, neg_r = 0;
	char from_flag;
	unsigned M = 0, N = 0;

	Node curr_node = mytree[node];
	feat = curr_node.feat;
	thres = curr_node.thres;
	cl = curr_node.cl;
	cr = curr_node.cr;
	cd = curr_node.cd;

	// short circut when this is a stump tree (with no splits)
	if (cl < 0) {
		out_contribs[num_feats] += curr_node.value;
		return;
	}

	//     if (DEBUG) {
	//       myfile << "\nNode: " << node << "\n";
	//       myfile << "x[feat]: " << x[feat] << ", r[feat]: " << r[feat] << "\n";
	//       myfile << "thres: " << thres << "\n";
	//     }

	if (x_missing[feat]) {
		next_xnode = cd;
	}
	else if (x[feat] > thres) {
		next_xnode = cr;
	}
	else if (x[feat] <= thres) {
		next_xnode = cl;
	}

	if (r_missing[feat]) {
		next_rnode = cd;
	}
	else if (r[feat] > thres) {
		next_rnode = cr;
	}
	else if (r[feat] <= thres) {
		next_rnode = cl;
	}

	if (next_xnode != next_rnode) {
		mytree[next_xnode].from_flag = FROM_X_NOT_R;
		mytree[next_rnode].from_flag = FROM_R_NOT_X;
	}
	else {
		mytree[next_xnode].from_flag = FROM_NEITHER;
	}

	// Check if x and r go the same way
	if (next_xnode == next_rnode) {
		next_node = next_xnode;
	}

	// If not, go left
	if (next_node < 0) {
		next_node = cl;
		if (next_rnode == next_node) { // rpath
			N = N + 1;
			feat_hist[feat] -= 1;
		}
		else if (next_xnode == next_node) { // xpath
			M = M + 1;
			N = N + 1;
			feat_hist[feat] += 1;
		}
	}
	node_stack[ns_ctr] = node;
	ns_ctr += 1;
	while (true) {
		node = next_node;
		curr_node = mytree[node];
		feat = curr_node.feat;
		thres = curr_node.thres;
		cl = curr_node.cl;
		cr = curr_node.cr;
		cd = curr_node.cd;
		pnode = curr_node.pnode;
		pfeat = curr_node.pfeat;
		from_flag = curr_node.from_flag;



		//         if (DEBUG) {
		//           myfile << "\nNode: " << node << "\n";
		//           myfile << "N: " << N << ", M: " << M << "\n";
		//           myfile << "from_flag==FROM_X_NOT_R: " << (from_flag==FROM_X_NOT_R) << "\n";
		//           myfile << "from_flag==FROM_R_NOT_X: " << (from_flag==FROM_R_NOT_X) << "\n";
		//           myfile << "from_flag==FROM_NEITHER: " << (from_flag==FROM_NEITHER) << "\n";
		//           myfile << "feat_hist[feat]: " << feat_hist[feat] << "\n";
		//         }

		// At a leaf
		if (cl < 0) {
			//      if (DEBUG) {
			//        myfile << "At a leaf\n";
			//      }

			if (M == 0) {
				out_contribs[num_feats] += mytree[node].value;
			}

			// Currently assuming a single output
			if (N != 0) {
				if (M != 0) {
					pos_lst[node] = mytree[node].value * memoized_weights[N + max_depth * (M - 1)];
				}
				if (M != N) {
					neg_lst[node] = -mytree[node].value * memoized_weights[N + max_depth * M];
				}
			}
			//             if (DEBUG) {
			//               myfile << "pos_lst[node]: " << pos_lst[node] << "\n";
			//               myfile << "neg_lst[node]: " << neg_lst[node] << "\n";
			//             }
			// Pop from node_stack
			ns_ctr -= 1;
			next_node = node_stack[ns_ctr];
			from_child = node;
			// Unwind
			if (feat_hist[pfeat] > 0) {
				feat_hist[pfeat] -= 1;
			}
			else if (feat_hist[pfeat] < 0) {
				feat_hist[pfeat] += 1;
			}
			if (feat_hist[pfeat] == 0) {
				if (from_flag == FROM_X_NOT_R) {
					N = N - 1;
					M = M - 1;
				}
				else if (from_flag == FROM_R_NOT_X) {
					N = N - 1;
				}
			}
			continue;
		}

		const bool x_right = x[feat] > thres;
		const bool r_right = r[feat] > thres;

		if (x_missing[feat]) {
			next_xnode = cd;
		}
		else if (x_right) {
			next_xnode = cr;
		}
		else if (!x_right) {
			next_xnode = cl;
		}

		if (r_missing[feat]) {
			next_rnode = cd;
		}
		else if (r_right) {
			next_rnode = cr;
		}
		else if (!r_right) {
			next_rnode = cl;
		}

		if (next_xnode >= 0) {
			if (next_xnode != next_rnode) {
				mytree[next_xnode].from_flag = FROM_X_NOT_R;
				mytree[next_rnode].from_flag = FROM_R_NOT_X;
			}
			else {
				mytree[next_xnode].from_flag = FROM_NEITHER;
			}
		}

		// Arriving at node from parent
		if (from_child == -1) {
			//      if (DEBUG) {
			//        myfile << "Arriving at node from parent\n";
			//      }
			node_stack[ns_ctr] = node;
			ns_ctr += 1;
			next_node = -1;

			//      if (DEBUG) {
			//        myfile << "feat_hist[feat]" << feat_hist[feat] << "\n";
			//      }
			// Feature is set upstream
			if (feat_hist[feat] > 0) {
				next_node = next_xnode;
				feat_hist[feat] += 1;
			}
			else if (feat_hist[feat] < 0) {
				next_node = next_rnode;
				feat_hist[feat] -= 1;
			}

			// x and r go the same way
			if (next_node < 0) {
				if (next_xnode == next_rnode) {
					next_node = next_xnode;
				}
			}

			// Go down one path
			if (next_node >= 0) {
				continue;
			}

			// Go down both paths, but go left first
			next_node = cl;
			if (next_rnode == next_node) {
				N = N + 1;
				feat_hist[feat] -= 1;
			}
			else if (next_xnode == next_node) {
				M = M + 1;
				N = N + 1;
				feat_hist[feat] += 1;
			}
			from_child = -1;
			continue;
		}

		// Arriving at node from child
		if (from_child != -1) {
			//             if (DEBUG) {
			//               myfile << "Arriving at node from child\n";
			//             }
			next_node = -1;
			// Check if we should unroll immediately
			if ((next_rnode == next_xnode) || (feat_hist[feat] != 0)) {
				next_node = pnode;
			}

			// Came from a single path, so unroll
			if (next_node >= 0) {
				//                 if (DEBUG) {
				//                   myfile << "Came from a single path, so unroll\n";
				//                 }
				// At the root node
				if (node == 0) {
					break;
				}
				// Update and unroll
				pos_lst[node] = pos_lst[from_child];
				neg_lst[node] = neg_lst[from_child];

				//                 if (DEBUG) {
				//                   myfile << "pos_lst[node]: " << pos_lst[node] << "\n";
				//                   myfile << "neg_lst[node]: " << neg_lst[node] << "\n";
				//                 }
				from_child = node;
				ns_ctr -= 1;

				// Unwind
				if (feat_hist[pfeat] > 0) {
					feat_hist[pfeat] -= 1;
				}
				else if (feat_hist[pfeat] < 0) {
					feat_hist[pfeat] += 1;
				}
				if (feat_hist[pfeat] == 0) {
					if (from_flag == FROM_X_NOT_R) {
						N = N - 1;
						M = M - 1;
					}
					else if (from_flag == FROM_R_NOT_X) {
						N = N - 1;
					}
				}
				continue;
				// Go right - Arriving from the left child
			}
			else if (from_child == cl) {
				//                 if (DEBUG) {
				//                   myfile << "Go right - Arriving from the left child\n";
				//                 }
				node_stack[ns_ctr] = node;
				ns_ctr += 1;
				next_node = cr;
				if (next_xnode == next_node) {
					M = M + 1;
					N = N + 1;
					feat_hist[feat] += 1;
				}
				else if (next_rnode == next_node) {
					N = N + 1;
					feat_hist[feat] -= 1;
				}
				from_child = -1;
				continue;
				// Compute stuff and unroll - Arriving from the right child
			}
			else if (from_child == cr) {
				//                 if (DEBUG) {
				//                   myfile << "Compute stuff and unroll - Arriving from the right child\n";
				//                 }
				pos_x = 0;
				neg_x = 0;
				pos_r = 0;
				neg_r = 0;
				if ((next_xnode == cr) && (next_rnode == cl)) {
					pos_x = pos_lst[cr];
					neg_x = neg_lst[cr];
					pos_r = pos_lst[cl];
					neg_r = neg_lst[cl];
				}
				else if ((next_xnode == cl) && (next_rnode == cr)) {
					pos_x = pos_lst[cl];
					neg_x = neg_lst[cl];
					pos_r = pos_lst[cr];
					neg_r = neg_lst[cr];
				}
				// out_contribs needs to have been initialized as all zeros
				// if (pos_x + neg_r != 0) {
				//   std::cout << "val " << pos_x + neg_r << "\n";
				// }
				out_contribs[feat] += pos_x + neg_r;
				pos_lst[node] = pos_x + pos_r;
				neg_lst[node] = neg_x + neg_r;

				//                 if (DEBUG) {
				//                   myfile << "out_contribs[feat]: " << out_contribs[feat] << "\n";
				//                   myfile << "pos_lst[node]: " << pos_lst[node] << "\n";
				//                   myfile << "neg_lst[node]: " << neg_lst[node] << "\n";
				//                 }

				// Check if at root
				if (node == 0) {
					break;
				}

				// Pop
				ns_ctr -= 1;
				next_node = node_stack[ns_ctr];
				from_child = node;

				// Unwind
				if (feat_hist[pfeat] > 0) {
					feat_hist[pfeat] -= 1;
				}
				else if (feat_hist[pfeat] < 0) {
					feat_hist[pfeat] += 1;
				}
				if (feat_hist[pfeat] == 0) {
					if (from_flag == FROM_X_NOT_R) {
						N = N - 1;
						M = M - 1;
					}
					else if (from_flag == FROM_R_NOT_X) {
						N = N - 1;
					}
				}
				continue;
			}
		}
	}
	//  if (DEBUG) {
	//    myfile.close();
	//  }
}


/**
* Runs Tree SHAP with feature independence assumptions on dense data.
*/
void dense_independent(const TreeEnsemble& trees, const ExplanationDataset &data,
	tfloat *out_contribs, tfloat transform(const tfloat, const tfloat)) {

	// reformat the trees for faster access
	Node *node_trees = new Node[trees.tree_limit * trees.max_nodes];
	for (unsigned i = 0; i < trees.tree_limit; ++i) {
		Node *node_tree = node_trees + i * trees.max_nodes;
		for (unsigned j = 0; j < trees.max_nodes; ++j) {
			const unsigned en_ind = i * trees.max_nodes + j;
			node_tree[j].cl = trees.children_left[en_ind];
			node_tree[j].cr = trees.children_right[en_ind];
			node_tree[j].cd = trees.children_default[en_ind];
			if (j == 0) {
				node_tree[j].pnode = 0;
			}
			if (trees.children_left[en_ind] >= 0) { // relies on all unused entries having negative values in them
				node_tree[trees.children_left[en_ind]].pnode = j;
				node_tree[trees.children_left[en_ind]].pfeat = trees.features[en_ind];
			}
			if (trees.children_right[en_ind] >= 0) { // relies on all unused entries having negative values in them
				node_tree[trees.children_right[en_ind]].pnode = j;
				node_tree[trees.children_right[en_ind]].pfeat = trees.features[en_ind];
			}

			node_tree[j].thres = (float)trees.thresholds[en_ind];
			node_tree[j].feat = trees.features[en_ind];
		}
	}

	// preallocate arrays needed by the algorithm
	float *pos_lst = new float[trees.max_nodes];
	float *neg_lst = new float[trees.max_nodes];
	int *node_stack = new int[(unsigned)trees.max_depth];
	signed short *feat_hist = new signed short[data.M];
	tfloat *tmp_out_contribs = new tfloat[(data.M + 1)];

	// precompute all the weight coefficients
	float *memoized_weights = new float[(trees.max_depth + 1) * (trees.max_depth + 1)];
	for (unsigned n = 0; n <= trees.max_depth; ++n) {
		for (unsigned m = 0; m <= trees.max_depth; ++m) {
			memoized_weights[n + trees.max_depth * m] = float(1.0 / (n * bin_coeff(n - 1, m)));
		}
	}

	// compute the explanations for each sample
	tfloat *instance_out_contribs;
	tfloat rescale_factor = 1.0;
	tfloat margin_x = 0;
	tfloat margin_r = 0;
	for (unsigned oind = 0; oind < trees.num_outputs; ++oind) {
		// set the values int he reformated tree to the current output index
		for (unsigned i = 0; i < trees.tree_limit; ++i) {
			Node *node_tree = node_trees + i * trees.max_nodes;
			for (unsigned j = 0; j < trees.max_nodes; ++j) {
				const unsigned en_ind = i * trees.max_nodes + j;
				node_tree[j].value = (float)trees.values[en_ind * trees.num_outputs + oind];
			}
		}

		// loop over all the samples
		for (unsigned i = 0; i < data.num_X; ++i) {
			const tfloat *x = data.X + i * data.M;
			const bool *x_missing = data.X_missing + i * data.M;
			instance_out_contribs = out_contribs + i * (data.M + 1) * trees.num_outputs;
			const tfloat y_i = data.y == NULL ? 0 : data.y[i];

			//print_progress_bar(last_print, start_time, oind * data.num_X + i, data.num_X * trees.num_outputs);

			// compute the model's margin output for x
			if (transform != NULL) {
				margin_x = trees.base_offset;
				for (unsigned k = 0; k < trees.tree_limit; ++k) {
					margin_x += tree_predict(k, trees, x, x_missing)[oind];
				}
			}

			for (unsigned j = 0; j < data.num_R; ++j) {
				const tfloat *r = data.R + j * data.M;
				const bool *r_missing = data.R_missing + j * data.M;
				std::fill_n(tmp_out_contribs, (data.M + 1), 0);

				// compute the model's margin output for r
				if (transform != NULL) {
					margin_r = trees.base_offset;
					for (unsigned k = 0; k < trees.tree_limit; ++k) {
						margin_r += tree_predict(k, trees, r, r_missing)[oind];
					}
				}

				for (unsigned k = 0; k < trees.tree_limit; ++k) {
					tree_shap_indep(
						trees.max_depth, data.M, trees.max_nodes, x, x_missing, r, r_missing,
						tmp_out_contribs, pos_lst, neg_lst, feat_hist, memoized_weights,
						node_stack, node_trees + k * trees.max_nodes
					);
				}

				// compute the rescale factor
				if (transform != NULL) {
					if (margin_x == margin_r) {
						rescale_factor = 1.0;
					}
					else {
						rescale_factor = (*transform)(margin_x, y_i) - (*transform)(margin_r, y_i);
						rescale_factor /= margin_x - margin_r;
					}
				}

				// add the effect of the current reference to our running total
				// this is where we can do per reference scaling for non-linear transformations
				for (unsigned k = 0; k < data.M; ++k) {
					instance_out_contribs[k * trees.num_outputs + oind] += tmp_out_contribs[k] * rescale_factor;
				}

				// Add the base offset
				if (transform != NULL) {
					instance_out_contribs[data.M * trees.num_outputs + oind] += (*transform)(trees.base_offset + tmp_out_contribs[data.M], 0);
				}
				else {
					instance_out_contribs[data.M * trees.num_outputs + oind] += trees.base_offset + tmp_out_contribs[data.M];
				}
			}

			// average the results over all the references.
			for (unsigned j = 0; j < (data.M + 1); ++j) {
				instance_out_contribs[j * trees.num_outputs + oind] /= data.num_R;
			}

			// apply the base offset to the bias term
			// for (unsigned j = 0; j < trees.num_outputs; ++j) {
			//     instance_out_contribs[data.M * trees.num_outputs + j] += (*transform)(trees.base_offset, 0);
			// }
		}
	}

	delete[] tmp_out_contribs;
	delete[] node_trees;
	delete[] pos_lst;
	delete[] neg_lst;
	delete[] node_stack;
	delete[] feat_hist;
	delete[] memoized_weights;
}


/**
* This runs Tree SHAP with a per tree path conditional dependence assumption.
*/
void dense_tree_path_dependent(const TreeEnsemble& trees, const ExplanationDataset &data,
	tfloat *out_contribs, tfloat transform(const tfloat, const tfloat)) {

	MedTimer tm;
	tm.start();
	chrono::high_resolution_clock::time_point tm_prog = chrono::high_resolution_clock::now();
	int progress = 0;
	int max_loop = data.num_X;
	// build explanation for each sample
#pragma omp parallel for
	for (int i = 0; i < data.num_X; ++i) {
		tfloat *instance_out_contribs = out_contribs + i * (data.M + 1) * trees.num_outputs;
		ExplanationDataset instance;
		TreeEnsemble tree;
		data.get_x_instance(instance, i);

		// aggregate the effect of explaining each tree
		// (this works because of the linearity property of Shapley values)
		for (unsigned j = 0; j < trees.tree_limit; ++j) {
			trees.get_tree(tree, j);
			tree_shap(tree, instance, instance_out_contribs, 0, 0);
		}

		// apply the base offset to the bias term
		for (unsigned j = 0; j < trees.num_outputs; ++j) {
			instance_out_contribs[data.M * trees.num_outputs + j] += trees.base_offset;
		}

#pragma omp atomic
		++progress;
		double duration = (unsigned long long)(chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
			- tm_prog).count()) / 1000000.0;
		if (duration > 15 && progress % 50 == 0) {
#pragma omp critical
			tm_prog = chrono::high_resolution_clock::now();
			double time_elapsed = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now()
				- tm.t[0]).count()) / 1000000.0;
			double estimate_time = int(double(max_loop - progress) / double(progress) * double(time_elapsed));
			MLOG("SHAPLEY Processed %d out of %d(%2.2f%%) time elapsed: %2.1f Minutes, "
				"estimate time to finish %2.1f Minutes\n", progress, max_loop, 100.0*(progress / float(max_loop)), time_elapsed / 60,
				estimate_time / 60.0);
		}
	}
}

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
	tfloat transform(const tfloat, const tfloat)) {

	// build a list of all the unique features in each tree
	int *unique_features = new int[trees.tree_limit * trees.max_nodes];
	std::fill(unique_features, unique_features + trees.tree_limit * trees.max_nodes, -1);
	for (unsigned j = 0; j < trees.tree_limit; ++j) {
		const int *features_row = trees.features + j * trees.max_nodes;
		int *unique_features_row = unique_features + j * trees.max_nodes;
		for (unsigned k = 0; k < trees.max_nodes; ++k) {
			for (unsigned l = 0; l < trees.max_nodes; ++l) {
				if (features_row[k] == unique_features_row[l]) break;
				if (unique_features_row[l] < 0) {
					unique_features_row[l] = features_row[k];
					break;
				}
			}
		}
	}

	// build an interaction explanation for each sample
	tfloat *instance_out_contribs;
	TreeEnsemble tree;
	ExplanationDataset instance;
	const unsigned contrib_row_size = (data.M + 1) * trees.num_outputs;
	tfloat *diag_contribs = new tfloat[contrib_row_size];
	tfloat *on_contribs = new tfloat[contrib_row_size];
	tfloat *off_contribs = new tfloat[contrib_row_size];
	for (unsigned i = 0; i < data.num_X; ++i) {
		instance_out_contribs = out_contribs + i * (data.M + 1) * contrib_row_size;
		data.get_x_instance(instance, i);

		// aggregate the effect of explaining each tree
		// (this works because of the linearity property of Shapley values)
		std::fill(diag_contribs, diag_contribs + contrib_row_size, 0);
		for (unsigned j = 0; j < trees.tree_limit; ++j) {
			trees.get_tree(tree, j);
			tree_shap(tree, instance, diag_contribs, 0, 0);

			const int *unique_features_row = unique_features + j * trees.max_nodes;
			for (unsigned k = 0; k < trees.max_nodes; ++k) {
				const int ind = unique_features_row[k];
				if (ind < 0) break; // < 0 means we have seen all the features for this tree

									// compute the shap value with this feature held on and off
				std::fill(on_contribs, on_contribs + contrib_row_size, 0);
				std::fill(off_contribs, off_contribs + contrib_row_size, 0);
				tree_shap(tree, instance, on_contribs, 1, ind);
				tree_shap(tree, instance, off_contribs, -1, ind);

				// save the difference between on and off as the interaction value
				for (unsigned l = 0; l < contrib_row_size; ++l) {
					const tfloat val = (on_contribs[l] - off_contribs[l]) / 2;
					instance_out_contribs[ind * contrib_row_size + l] += val;
					diag_contribs[l] -= val;
				}
			}
		}

		// set the diagonal
		for (unsigned j = 0; j < data.M + 1; ++j) {
			const unsigned offset = j * contrib_row_size + j * trees.num_outputs;
			for (unsigned k = 0; k < trees.num_outputs; ++k) {
				instance_out_contribs[offset + k] = diag_contribs[j * trees.num_outputs + k];
			}
		}

		// apply the base offset to the bias term
		const unsigned last_ind = (data.M * (data.M + 1) + data.M) * trees.num_outputs;
		for (unsigned j = 0; j < trees.num_outputs; ++j) {
			instance_out_contribs[last_ind + j] += trees.base_offset;
		}
	}

	delete[] diag_contribs;
	delete[] on_contribs;
	delete[] off_contribs;
	delete[] unique_features;
}

/**
* This runs Tree SHAP with a global path conditional dependence assumption.
*
* By first merging all the trees in a tree ensemble into an equivalent single tree
* this method allows arbitrary marginal transformations and also ensures that all the
* evaluations of the model are consistent with some training data point.
*/
void dense_global_path_dependent(const TreeEnsemble& trees, const ExplanationDataset &data,
	tfloat *out_contribs, tfloat transform(const tfloat, const tfloat)) {

	// allocate space for our new merged tree (we save enough room to totally split all samples if need be)
	TreeEnsemble merged_tree;
	merged_tree.allocate(1, (data.num_X + data.num_R) * 2, trees.num_outputs);

	// collapse the ensemble of trees into a single tree that has the same behavior
	// for all the X and R samples in the dataset
	build_merged_tree(merged_tree, data, trees);

	// compute the expected value and depth of the new merged tree
	compute_expectations(merged_tree);

	// explain each sample using our new merged tree
	ExplanationDataset instance;
	tfloat *instance_out_contribs;
	for (unsigned i = 0; i < data.num_X; ++i) {
		instance_out_contribs = out_contribs + i * (data.M + 1) * trees.num_outputs;
		data.get_x_instance(instance, i);

		// since we now just have a single merged tree we can just use the tree_path_dependent algorithm
		tree_shap(merged_tree, instance, instance_out_contribs, 0, 0);

		// apply the base offset to the bias term
		for (unsigned j = 0; j < trees.num_outputs; ++j) {
			instance_out_contribs[data.M * trees.num_outputs + j] += trees.base_offset;
		}
	}

	merged_tree.free();
}


/**
* The main method for computing Tree SHAP on model using dense data.
*/
void dense_tree_shap(const TreeEnsemble& trees, const ExplanationDataset &data, tfloat *out_contribs,
	const int feature_dependence, unsigned model_transform, bool interactions) {

	// see what transform (if any) we have
	tfloat(*transform)(const tfloat margin, const tfloat y) = NULL;
	switch (model_transform) {
	case MODEL_TRANSFORM::logistic:
		transform = logistic_transform;
		break;

	case MODEL_TRANSFORM::logistic_nlogloss:
		transform = logistic_nlogloss_transform;
		break;

	case MODEL_TRANSFORM::squared_loss:
		transform = squared_loss_transform;
		break;
	}

	// dispatch to the correct algorithm handler
	switch (feature_dependence) {
	case FEATURE_DEPENDENCE::independent:
		if (interactions) {
			std::cerr << "FEATURE_DEPENDENCE::independent does not support interactions!\n";
		}
		else dense_independent(trees, data, out_contribs, transform);
		return;

	case FEATURE_DEPENDENCE::tree_path_dependent:
		if (interactions) dense_tree_interactions_path_dependent(trees, data, out_contribs, transform);
		else dense_tree_path_dependent(trees, data, out_contribs, transform);
		return;

	case FEATURE_DEPENDENCE::global_path_dependent:
		if (interactions) {
			std::cerr << "FEATURE_DEPENDENCE::global_path_dependent does not support interactions!\n";
		}
		else dense_global_path_dependent(trees, data, out_contribs, transform);
		return;
	}
}

long medial::shapley::nchoosek(long n, long k) {
	if (k == 0)
		return 1;
	return (n * nchoosek(n - 1, k - 1)) / k;
}

void medial::shapley::list_all_options_binary(int nfeats, vector<vector<bool>> &all_opts) {
	//int all_opts_count = nchoosek(nfeats, k);
	//all_opts.reserve(all_opts_count);
	vector<vector<int>> all_opts_num;
	for (size_t k = 1; k < nfeats; ++k)
	{
		vector<vector<int>> all_opts_num_batch; //option_num, list of selected idx
		all_opts_num_batch.push_back({}); //builds all options for k by bottom to up
		for (int i = 1; i <= k; ++i)
		{
			//advance to i step when passing on all current options:
			int curr_opts_count = nchoosek(nfeats, i);
			vector<vector<int>> curr_options(curr_opts_count);
			int opt_num = 0;
			for (size_t j = 0; j < all_opts_num_batch.size(); ++j)
			{
				const vector<int> &curr_opt = all_opts_num_batch[j];
				int st = 0;
				if (!curr_opt.empty())
					st = curr_opt.back() + 1;
				for (int k = st; k < nfeats; ++k) {
					//add new option:
					curr_options[opt_num] = curr_opt;
					curr_options[opt_num].push_back(k);
					++opt_num;
				}
			}
			//update all_opts:
			all_opts_num_batch.swap(curr_options);
		}
		all_opts_num.insert(all_opts_num.end(), all_opts_num_batch.begin(), all_opts_num_batch.end());
	}
	//converts all_opts_num to all_opts
	all_opts.resize(all_opts_num.size());
	for (size_t i = 0; i < all_opts.size(); ++i)
	{
		all_opts[i].resize(nfeats, false);
		for (int sel : all_opts_num[i])
			all_opts[i][sel] = true;
	}
}

void medial::shapley::generate_mask(vector<bool> &mask, int nfeat, mt19937 &gen, bool uniform_rand, bool use_shuffle) {
	mask.clear(); //If has values

	if (uniform_rand) {
		mask.reserve(nfeat);
		uniform_int_distribution<> rnd_coin(0, 1);
		for (size_t i = 0; i < nfeat; ++i)
			mask.push_back(rnd_coin(gen) > 0);
		return;
	}

	uniform_int_distribution<> rnd_dist(0, nfeat);
	int zero_count = rnd_dist(gen);
	if (zero_count == nfeat) {
		mask.resize(nfeat);
		return;
	}
	else if (zero_count == 0) {
		mask.resize(nfeat, true);
		return;
	}
	mask.resize(nfeat, true);


	if (use_shuffle) {
		//calc shuffle of indexes and takes top - may be faster for small arrays
		vector<int> sel_idx; sel_idx.reserve(nfeat);
		for (int i = 0; i < nfeat; ++i)
			sel_idx.push_back(i);
		shuffle(sel_idx.begin(), sel_idx.end(), gen);
		for (size_t i = 0; i < zero_count; ++i)
			mask[sel_idx[i]] = false;
	}
	else {
		//calc with no reapets till converges:
		vector<bool> seen_idx(nfeat);
		uniform_int_distribution<> rnd_feat(0, nfeat - 1);
		for (size_t i = 0; i < zero_count; ++i)
		{
			int sel = rnd_feat(gen);
			while (seen_idx[sel])
				sel = rnd_feat(gen);
			seen_idx[sel] = true;
			mask[sel] = false;
		}

	}
}

void medial::shapley::generate_mask_(vector<bool> &mask, int nfeat, mt19937 &gen, bool uniform_rand, bool use_shuffle) {
	//mask is not empty - sample from non-empty
	if (uniform_rand) {
		uniform_int_distribution<> rnd_coin(0, 1);
		for (size_t i = 0; i < nfeat; ++i)
			if (mask[i])
				mask[i] = (rnd_coin(gen) > 0);
		return;
	}
	int curr_cnt = 0;
	for (size_t i = 0; i < mask.size(); ++i)
		curr_cnt += int(!mask[i]); //how many zeros now

	uniform_int_distribution<> rnd_dist(0, nfeat);
	int zero_count = rnd_dist(gen);
	if (zero_count == nfeat) {
		mask.clear(); //all zeros
		mask.resize(nfeat);
		return;
	}
	else if (zero_count <= curr_cnt) {
		return; //no change in mask is needed, has more zeros than requested
	}

	if (use_shuffle) {
		//calc shuffle of indexes and takes top - may be faster for small arrays
		vector<int> sel_idx; sel_idx.reserve(nfeat);
		for (int i = 0; i < nfeat; ++i)
			if (mask[i]) //list only potentials
				sel_idx.push_back(i);
		shuffle(sel_idx.begin(), sel_idx.end(), gen);
		for (size_t i = 0; i < zero_count - curr_cnt; ++i)
			mask[sel_idx[i]] = false;
	}
	else {
		//calc with no reapets till converges:
		vector<bool> seen_idx(nfeat);
		for (size_t i = 0; i < mask.size(); ++i)
			if (!mask[i])
				seen_idx[i] = true;

		uniform_int_distribution<> rnd_feat(0, nfeat - 1);
		for (size_t i = 0; i < zero_count - curr_cnt; ++i)
		{
			int sel = rnd_feat(gen);
			while (seen_idx[sel])
				sel = rnd_feat(gen);
			seen_idx[sel] = true;
			mask[sel] = false;
		}

	}
}

/**
* sample uniformally masks - first by choosing number of 0's in the mask(set missings) and than select them randomally
*/
void medial::shapley::sample_options_SHAP(int nfeats, vector<vector<bool>> &all_opts, int opt_count, mt19937 &gen, bool with_repeats
	, bool uniform_rand, bool use_shuffle) {
	all_opts.resize(opt_count);
	unordered_set<vector<bool>> seen_mask;
	for (size_t i = 0; i < all_opts.size(); ++i)
	{
		vector<bool> &cr_opt = all_opts[i];
		//generate mask:
		generate_mask(cr_opt, nfeats, gen);
		while (!with_repeats && seen_mask.find(cr_opt) != seen_mask.end())
			generate_mask(cr_opt, nfeats, gen);
		//finished - already manipulated memory
		if (!with_repeats)
			seen_mask.insert(cr_opt);
	}
}

double medial::shapley::get_c(int p1, int p2, int end_l) {
	//c := (end_l)! / ( p1! * p2! )
	//returns 1/c 
	double c = 1, d = 1;
	if (p1 > p2) {
		for (int i = p1 + 1; i <= end_l; ++i)
			c *= i;
		for (int i = 2; i <= p2; ++i)
			d *= i;
		c /= d;
	}
	else {
		for (int i = p2 + 1; i <= end_l; ++i)
			c *= i;
		for (int i = 2; i <= p1; ++i)
			d *= i;
		c /= d;
	}

	return 1 / c;
}

void medial::shapley::explain_shapley(const MedFeatures &matrix, int selected_sample, int max_tests,
	MedPredictor *predictor, float missing_value, const vector<vector<int>>& group2index, const vector<string> &groupNames
	, vector<float> &features_coeff,
	bool sample_masks_with_repeats, float select_from_all, bool uniform_rand, bool use_shuffle,
	bool verbose) {
	random_device rd;
	mt19937 gen(rd());

	//int ngrps = (int)group2index.size();

	int tot_feat_cnt = (int)matrix.data.size();

	vector<string> full_feat_ls;
	matrix.get_feature_names(full_feat_ls);
	vector<float> fast_access(full_feat_ls.size());
	for (size_t i = 0; i < full_feat_ls.size(); ++i)
		fast_access[i] = matrix.data.at(full_feat_ls[i])[selected_sample];

	features_coeff.resize(tot_feat_cnt);

	//calc shapley for each variable
	if (verbose)
		MLOG("Start explain_shapely\n");
	MedTimer tm_taker;
	tm_taker.start();
	bool warn_shown = false;
	MedProgress progress_full("shapley", tot_feat_cnt, 15);
	//bool sample_masks_with_repeats = false;
	//float select_from_all = (float)0.8; //If 80% of all options than sample from all
	//bool uniform_rand = false;
	//bool use_shuffle = true;

	for (size_t param_i = 0; param_i < tot_feat_cnt; ++param_i)
	{
		double phi_i = 0;
		if (fast_access[param_i] == missing_value)
			continue;
		//iterate on all other features  execpt param_i, and other features that are already missing in the given example
		vector<string> candidates; //for sampling options from - skip missings
		for (size_t i = 0; i < full_feat_ls.size(); ++i)
		{
			if (i == param_i || fast_access[i] == missing_value)
				continue;
			candidates.push_back(full_feat_ls[i]);
		}

		vector<vector<bool>> all_opts;
		bool iter_all = true;
		int max_loop = max_tests;
		double nchoose = pow(2, (int)candidates.size());
		if (nchoose < max_loop)
			max_loop = nchoose;
		else {
			iter_all = false;
			if (!warn_shown && verbose)
				MLOG("Warning have %d options, and max_test is %d\n", (int)nchoose, max_loop);
		}

		if (iter_all || float(nchoose) / max_tests >= select_from_all) {
			list_all_options_binary((int)candidates.size(), all_opts);
			vector<bool> empty_vec(candidates.size()), full_vec(candidates.size(), true);
			all_opts.push_back(empty_vec);
			all_opts.push_back(full_vec);

			//select random masks from all options when not iterating all options:
			if (!iter_all) {
				if (!warn_shown && verbose) {
					MLOG("Warning: not iterating all in feature %zu has %zu candidates, has %d options, max_test=%d\n",
						param_i, candidates.size(), (int)nchoose, max_loop);
#pragma omp critical
					warn_shown = true;
				}
				shuffle(all_opts.begin(), all_opts.end(), gen);
				all_opts.resize(max_loop);
			}
		}
		else
			sample_options_SHAP((int)candidates.size(), all_opts, max_loop, gen, sample_masks_with_repeats, uniform_rand, use_shuffle);

		//complete all_opts to nfeats size:
		bool deafult_not_selected = true; //mark all the rest(missing values that aren't tested) as fixed to missing value
		for (size_t i = 0; i < all_opts.size(); ++i)
		{
			vector<bool> mask(full_feat_ls.size(), deafult_not_selected);
			vector<bool> &truncated_mask = all_opts[i];
			int cand_idx = 0;
			for (size_t j = 0; j < full_feat_ls.size(); ++j)
			{
				if (cand_idx < candidates.size() && full_feat_ls[j] == candidates[cand_idx]) {
					mask[j] = truncated_mask[cand_idx];
					++cand_idx;
				}
				//else keeps default (was missing value) - which means keep as missing value or iterate through
			}
			mask[param_i] = false; //mark always as false the selected feature to test
			truncated_mask = move(mask);
		}

		//build test matrix from all_opts:
		MedTimer tm;
		tm.start();

		//collect score for each permutition of missing values:
		int end_l = (int)candidates.size();
		vector<float> full_pred_all_masks_without(max_loop * tot_feat_cnt), full_pred_all_masks_with(max_loop* tot_feat_cnt);
		for (int i = 0; i < max_loop; ++i) {
			float *mat_without = &full_pred_all_masks_without[i * tot_feat_cnt];
			float *mat_with = &full_pred_all_masks_with[i * tot_feat_cnt];
			for (size_t j = 0; j < tot_feat_cnt; ++j)
				if (!all_opts[i][j])
					mat_without[j] = missing_value;
				else
					mat_without[j] = fast_access[j];
			for (size_t j = 0; j < tot_feat_cnt; ++j)
				if (!all_opts[i][j] && j != param_i)
					mat_with[j] = missing_value;
				else
					mat_with[j] = fast_access[j];
		}

		vector<float> preds_with, preds_without;
		predictor->predict(full_pred_all_masks_without, preds_without, max_loop, (int)full_feat_ls.size());
		predictor->predict(full_pred_all_masks_with, preds_with, max_loop, (int)full_feat_ls.size());

		for (int i = 0; i < max_loop; ++i) {
			int f_cnt = 0;
			for (size_t j = 0; j < all_opts[i].size(); ++j)
				f_cnt += int(all_opts[i][j]);

			float f_diff = preds_with[i] - preds_without[i];
			int p1 = f_cnt;
			int p2 = end_l - p1;
			double c = get_c(p1, p2, end_l + 1);
			phi_i += c * f_diff;
		}

		if (verbose)
			progress_full.update();
		features_coeff[param_i] = phi_i;
	}

	tm_taker.take_curr_time();
	if (verbose)
		MLOG("Done explain_shapely. took %2.1f seconds\n", tm_taker.diff_sec());
}

void collect_score_mask(const vector<float> &x, const vector<bool> &mask, SamplesGenerator<float> &sampler_gen
	, int sample_per_row, void *sampling_params, const vector<string> &feat_names, const map<string, FeatureAttr> &attr,
	const MedPredictor *predictor, vector<float> &preds) {

	if (sampler_gen.use_vector_api) {
		vector<vector<float>> mat_inp = { x };
		vector<vector<bool>> masks = { mask };
		vector<vector<float>> res; //the result matrix
		sampler_gen.get_samples(res, sample_per_row, sampling_params, masks, mat_inp);

		MedFeatures gen_matrix;
		gen_matrix.samples.resize(res.size());
		preds.resize(res.size());
		for (size_t i = 0; i < feat_names.size(); ++i)
		{
			gen_matrix.data[feat_names[i]].resize(res.size());
			for (size_t k = 0; k < res.size(); ++k)
				gen_matrix.data[feat_names[i]][k] = res[k][i];
		}
		gen_matrix.attributes = attr;
		gen_matrix.init_pid_pos_len();
		predictor->predict(gen_matrix); //todo: use vector api

		for (size_t i = 0; i < preds.size(); ++i)
			preds[i] = gen_matrix.samples[i].prediction[0];
	}
	else {
		MedFeatures gen_matrix;

		sampler_gen.get_samples(gen_matrix.data, sampling_params, mask, x);
		gen_matrix.samples.resize(gen_matrix.data.begin()->second.size());
		preds.resize(gen_matrix.samples.size());

		gen_matrix.attributes = attr;
		gen_matrix.init_pid_pos_len();
		predictor->predict(gen_matrix); //todo: use vector api

		for (size_t i = 0; i < preds.size(); ++i)
			preds[i] = gen_matrix.samples[i].prediction[0];
	}

}

int collect_mask(const vector<float> &x, const vector<bool> &mask, SamplesGenerator<float> &sampler_gen
	, int sample_per_row, void *sampling_params, const vector<string> &feat_names, map<string, vector<float>> &gen_matrix) {

	if (sampler_gen.use_vector_api) {
		int size_before = (int)gen_matrix[feat_names.front()].size();
		vector<vector<float>> mat_inp = { x };
		vector<vector<bool>> masks = { mask };
		vector<vector<float>> res; //the result matrix
		sampler_gen.get_samples(res, sample_per_row, sampling_params, masks, mat_inp);


#pragma omp critical
		for (size_t i = 0; i < feat_names.size(); ++i) {
			gen_matrix[feat_names[i]].resize(size_before + res.size());
			for (size_t k = 0; k < res.size(); ++k)
				gen_matrix[feat_names[i]][size_before + k] = res[k][i];
		}

		return (int)res.size();
	}
	else {
		//no parallel:
		int size_before = (int)gen_matrix[feat_names.front()].size();
		sampler_gen.get_samples(gen_matrix, sampling_params, mask, x);
		return (int)gen_matrix[feat_names.front()].size() - size_before;
	}
}

template<typename T> double mean_vec(const T *v, int len) {

	double s = 0;
	for (size_t i = 0; i < len; ++i)
		s += v[i];

	if (len == 0)
		MTHROW_AND_ERR("No values given for mean_vec. Cannot return anything valid\n");

	return s / len;
}

template<typename T> void medial::shapley::explain_shapley(const MedFeatures &matrix, int selected_sample, int max_tests,
	MedPredictor *predictor, const vector<vector<int>>& group2index, const vector<string> &groupNames,
	SamplesGenerator<T> &sampler_gen, int sample_per_row, void *sampling_params,
	vector<float> &features_coeff, bool verbose) {
	random_device rd;
	mt19937 gen(rd());

	int tot_feat_cnt = (int)matrix.data.size();
	vector<string> full_feat_ls;
	matrix.get_feature_names(full_feat_ls);
	vector<float> fast_access(tot_feat_cnt);
	for (size_t i = 0; i < full_feat_ls.size(); ++i)
		fast_access[i] = matrix.data.at(full_feat_ls[i])[selected_sample];

	int ngrps = (int)group2index.size();
	features_coeff.resize(ngrps);

	//calc shapley for each variable
	if (verbose)
		MLOG("Start explain_shapely\n");
	MedTimer tm_taker;
	tm_taker.start();
	bool warn_shown = false;
	float select_from_all = (float)0.8;
	MedProgress tm_full("", ngrps, 15);

	for (size_t param_i = 0; param_i < ngrps; ++param_i)
	{
		double phi_i = 0;
		int grps_opts = ngrps - 1;
		//iterate on all other features  execpt param_i, and other features that are already missing in the given example
		bool iter_all = true;
		int max_loop = max_tests;
		double nchoose = grps_opts <= 20 ? pow(2, grps_opts) : -1;
		if (grps_opts <= 20 && nchoose < max_loop)
			max_loop = nchoose;
		else {
			iter_all = false;
			if (!warn_shown && verbose && grps_opts <= 20)
				MLOG("Warning have %d options, and max_test is %d\n", (int)nchoose, max_loop);
		}

		vector<vector<bool>> all_opts;
		if (grps_opts <= 20 && (iter_all || float(nchoose) / max_tests >= select_from_all)) {
			list_all_options_binary(grps_opts, all_opts);
			vector<bool> empty_vec(grps_opts), full_vec(grps_opts, true);
			all_opts.push_back(empty_vec);
			all_opts.push_back(full_vec);

			//select random masks from all options when not iterating all options:
			if (!iter_all) {
				if (!warn_shown && verbose) {
					MLOG("Warning: not iterating all in feature %zu has %zu candidates, has %d options, max_test=%d\n",
						param_i, grps_opts, (int)nchoose, max_loop);
#pragma omp critical
					warn_shown = true;
				}
				shuffle(all_opts.begin(), all_opts.end(), gen);
				all_opts.resize(max_loop);
			}
		}
		else
			sample_options_SHAP(grps_opts, all_opts, max_loop, gen, false, true, false);

		//complete all_opts to nfeats size using groups:
		bool deafult_not_selected = true; //mark all the rest(missing values that aren't tested) as fixed to missing value
		for (size_t i = 0; i < all_opts.size(); ++i)
		{
			vector<bool> mask(full_feat_ls.size(), deafult_not_selected);
			vector<bool> &truncated_mask = all_opts[i];

			for (int j = 0; j < grps_opts; ++j)
			{
				bool mask_val = truncated_mask[j];
				int cand_idx = j + int(j >= param_i);
				for (int ind : group2index[cand_idx]) //set all group indexes as in mask
					mask[ind] = mask_val;
				//else keeps default (was missing value) - which means keep as missing value or iterate through
			}
			for (int ind : group2index[param_i])
				mask[ind] = false; //mark always as false the selected feature to test

			truncated_mask = move(mask);
		}

		//build test matrix from all_opts:
		MedProgress progress("VV_Shapley(Feat " + to_string(param_i + 1) +
			" out of " + to_string(ngrps) + ")", max_loop, 15, 1);
		//collect score for each permutition of missing values:
		int end_l = grps_opts;

		MedFeatures full_gen_samples;
		full_gen_samples.attributes = matrix.attributes;
		vector<int> splits_without(max_loop), splits_with(max_loop);

		for (int i = 0; i < max_loop; ++i) {

			int add_cnt = collect_mask(fast_access, all_opts[i], sampler_gen, sample_per_row, sampling_params,
				full_feat_ls, full_gen_samples.data);

			splits_without[i] = add_cnt;

			vector<bool> with_mask = all_opts[i];
			for (int ind : group2index[param_i])
				with_mask[ind] = true; //mark always as false the selected feature to test
			int add_cnt_with = collect_mask(fast_access, with_mask, sampler_gen, sample_per_row, sampling_params,
				full_feat_ls, full_gen_samples.data);
			splits_with[i] = add_cnt_with;

			if (verbose)
				progress.update();
		}
		full_gen_samples.samples.resize(full_gen_samples.data.begin()->second.size());
		full_gen_samples.init_pid_pos_len();

		predictor->predict(full_gen_samples);
		vector<float> preds_with, preds_without;
		preds_with.reserve(full_gen_samples.samples.size() / 2);
		preds_without.reserve(full_gen_samples.samples.size() / 2);
		vector<int> cumsum_without(max_loop), cumsum_with(max_loop);
		int curr_smp_pos = 0;
		for (int i = 0; i < max_loop; ++i) {
			int cnt_1 = splits_without[i];
			int cnt_2 = splits_with[i];

			int without_st = (int)preds_without.size();
			int with_st = (int)preds_with.size();
			cumsum_without[i] = without_st;
			cumsum_with[i] = with_st;

			for (size_t j = 0; j < cnt_1; ++j)
				preds_without.push_back(full_gen_samples.samples[curr_smp_pos + j].prediction[0]);
			for (size_t j = 0; j < cnt_2; ++j)
				preds_with.push_back(full_gen_samples.samples[curr_smp_pos + cnt_1 + j].prediction[0]);
			curr_smp_pos += cnt_1 + cnt_2;
		}

		for (int i = 0; i < max_loop; ++i) {
			float score_without, score_with;
			int f_cnt = 0;
			for (size_t j = 0; j < all_opts[i].size(); ++j)
				f_cnt += int(all_opts[i][j]);

			score_without = mean_vec(preds_without.data() + cumsum_without[i], splits_without[i]);
			score_with = mean_vec(preds_with.data() + cumsum_with[i], splits_with[i]);

			float f_diff = score_with - score_without;

			int p1 = f_cnt;
			int p2 = end_l - p1;
			double c = get_c(p1, p2, end_l + 1);
			phi_i += c * f_diff;
		}

		if (verbose)
			tm_full.update();
		features_coeff[param_i] = phi_i;
	}

	tm_taker.take_curr_time();
	if (verbose)
		MLOG("Done explain_shapely. took %2.1f seconds\n", tm_taker.diff_sec());

	if (verbose) {
		vector<pair<float, string>> feat_rank(ngrps);
		for (size_t i = 0; i < feat_rank.size(); ++i)
		{
			feat_rank[i].first = features_coeff[i];
			feat_rank[i].second = groupNames[i];
		}
		sort(feat_rank.begin(), feat_rank.end(), comp_score_flt_str);

		for (int i = 0; i < ngrps; ++i)
			MLOG("EXPLAIN #%d by %s : feat_score=%f\n",
				i + 1, feat_rank[i].second.c_str(), feat_rank[i].first);

		double sm = 0;
		for (int i = 0; i < ngrps - 1; ++i) {
			MLOG("%f +", feat_rank[i].first);
			sm += feat_rank[i].first;
		}
		sm += feat_rank[ngrps - 1].first;
		MLOG("%f=%f\n", feat_rank[ngrps - 1].first, sm);
	}
}

template void medial::shapley::explain_shapley<float>(const MedFeatures &matrix, int selected_sample, int max_tests,
	MedPredictor *predictor, const vector<vector<int>>& group2index, const vector<string> &groupNames,
	SamplesGenerator<float> &sampler_gen, int sample_per_row, void *sampling_params,
	vector<float> &features_coeff, bool verbose);

// Generate sampled matrix
void generate_samples(const MedFeatures& data, int isample, const vector<vector<bool>>& masks, SamplesGenerator<float> *generator,
	void *params, MedFeatures *out_data) {

	if (generator->use_vector_api) {
		int ncols = (int)data.data.size();
		vector<vector<float>> in(masks.size(), vector<float>(ncols)), out(masks.size(), vector<float>(ncols));

		int icol = 0;
		for (auto& rec : data.data) {
#pragma omp parallel for
			for (int irow = 0; irow < masks.size(); irow++)
				in[irow][icol] = rec.second[isample];
			icol++;
		}

		MLOG("Generate\n");
		generator->get_samples(out, 1, params, masks, in);

		icol = 0;
		for (auto& rec : data.data) {
			out_data->attributes[rec.first] = data.attributes.at(rec.first);
			out_data->data[rec.first].resize(masks.size());
#pragma omp parallel for
			for (int irow = 0; irow < masks.size(); irow++)
				(out_data->data)[rec.first][irow] = out[irow][icol];
			icol++;
		}
	}
	else {
		vector<string> all_n;
		data.get_feature_names(all_n);
		vector<const vector<float>*> data_p(data.data.size());
		for (size_t i = 0; i < data_p.size(); ++i)
			data_p[i] = &data.data.at(all_n[i]);
		vector<float> init_data(out_data->data.size());
		for (unsigned int icol = 0; icol < init_data.size(); ++icol)
			init_data[icol] = data_p[icol]->at(isample);

		for (int i = 0; i < masks.size(); ++i)
			generator->get_samples(out_data->data, params, masks[i], init_data);
	}
}

// Learn a Shapely-Lime model
void medial::shapley::get_shapley_lime_params(const MedFeatures& data, const MedPredictor *model,
	SamplesGenerator<float> *generator, float p, int n, float missing,
	void *params, const vector<vector<int>>& group2index, const vector<string>& group_names, vector<vector<float>>& alphas) {
	random_device rd;
	int N_TH = omp_get_max_threads();
	vector<mt19937> gen(N_TH);
	for (size_t i = 0; i < N_TH; ++i)
		gen[i] = mt19937(rd());
	uniform_real_distribution<> coin_dist(0, 1);
	uniform_int_distribution<> smp_choose(0, (int)data.samples.size() - 1);

	int nsamples = (int)data.samples.size();
	if (nsamples == 0) return;

	vector<string> features;
	data.get_feature_names(features);
	int nftrs = (int)features.size();
	int ngrps = (int)group_names.size();

	alphas.resize(nsamples);

	MedFeatures p_features;
	p_features.attributes = data.attributes;
	p_features.samples = data.samples;



	// for predictions - dummy Rep + samples + features
	for (int i = 0; i < n; i++) {
		//random select:
		int sel = smp_choose(gen[0]);
		p_features.samples.push_back(data.samples[sel]);
		//Add DATA:
		for (auto it = data.data.begin(); it != data.data.end(); ++it)
			p_features.data[it->first].push_back(it->second[sel]);
	}
	p_features.init_pid_pos_len();

	// Generate samples and predict
	MedTimer tm;
	for (int isample = 0; isample < nsamples; isample++) {
		tm.start();
		MLOG("Working on sample %d\n", isample);

		// Generate random masks
		MedMat<float> train(ngrps, n);
		vector<float> wgts(n);
		vector < vector<bool>> masks(n, vector<bool>(nftrs));

		vector<bool> missing_v(nftrs, false);
		int nMissing = 0;
		for (int i = 0; i < nftrs; i++) {
			if (data.data.at(features[i])[isample] == missing) {
				missing_v[i] = true;
				nMissing++;
			}
		}

		if (nMissing == nftrs)
			MTHROW_AND_ERR("All values are missing for entry %d , Cannot explain\n", isample);

#pragma omp parallel for
		for (int irow = 0; irow < n; irow++) {
			int n_th = omp_get_thread_num();
			// Mask
			int S = 0;
			while (S == 0 || S == ngrps) {
				S = 0;
				for (int igrp = 0; igrp < ngrps; igrp++) {
					bool grp_mask = coin_dist(gen[n_th]) > p;
					if (grp_mask) { // Keep, unless all are missing
						size_t nMissing = 0;
						for (int iftr : group2index[igrp]) {
							if (missing_v[iftr]) {
								nMissing++;
								masks[irow][iftr] = false;
							}
							else
								masks[irow][iftr] = true;
						}

						if (nMissing == group2index[igrp].size()) // All are missing ...
							train(igrp, irow) = 0.0;
						else {
							train(igrp, irow) = 1.0;
							S++;
						}
					}
					else { // Mask
						for (int iftr : group2index[igrp])
							masks[irow][iftr] = false;
						train(igrp, irow) = 0.0;
					}
				}
			}

			// Weights
			wgts[irow] = (ngrps - 1.0) / (medial::shapley::nchoosek(ngrps, S)*S*(ngrps - S));
		}
		train.transposed_flag = 1;

		// Generate sampled data
		generate_samples(data, isample, masks, generator, params, &p_features);

		// Get predictions for samples
		model->predict(p_features);

		vector<float> preds(n);
		double sum = 0;
		for (int irow = 0; irow < n; ++irow) {
			preds[irow] = p_features.samples[irow].prediction[0];
			sum += preds[irow];
		}

		MLOG("sample=%d, mean_pred=%f\n", isample, sum / n);

		// Learn linear model
		MedLM lm;
		lm.params.rfactor = (float)0.98;
		lm.learn(train, preds, wgts);

		// Extract alphas
		alphas[isample].resize(ngrps);
		for (int igrp = 0; igrp < ngrps; igrp++)
			alphas[isample][igrp] = lm.b[igrp];

		tm.take_curr_time();
		MLOG("Explaining sample took %2.1f sec\n", tm.diff_sec());

	}
}