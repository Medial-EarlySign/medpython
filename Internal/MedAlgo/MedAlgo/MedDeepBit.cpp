#include "MedDeepBit.h"
#include <MedUtils/MedUtils/MedGlobalRNG.h>
#include <SerializableObject/SerializableObject/SerializableObject.h>

//Initialization of the model
void MedDeepBit::init(float *x1, float *y1, int nsamples1, int nftrs1) {
	nftrs = nftrs1;
	vector<int> indexes(nsamples1);
	for (int i = 0; i < nsamples1; i++)
		indexes[i] = i;
	shuffle(indexes.begin(), indexes.end(), globalRNG::get_engine());
	vector<bool> is_train(nsamples1, true);
	internal_test_nsamples = (int) ((double) nsamples1 / params.internal_test_ratio);
	cout << "internal test size = " << internal_test_nsamples << endl;
	nsamples = nsamples1 - internal_test_nsamples;
	cout << "Train size = " << nsamples << endl;
	for (int i = 0; i < internal_test_nsamples; i++)
		is_train[indexes[i]] = false;
	x.resize(nftrs); 
	internal_test_x.resize(nftrs);
	for (int j = 0; j < nftrs; j++) {
		for (int i = 0; i < nsamples1; i++) {
			if (is_train[i])
				x[j].push_back(x1[i * nftrs + j]);
			else
				internal_test_x[j].push_back(x1[i * nftrs + j]);
		}
	}
	transpose(internal_test_x, internal_test_transposed_x);
	for (int i = 0; i < nsamples1; i++) {
		if (is_train[i]) {
			y.push_back(y1[i]);
			label.push_back((int)y1[i]);
		}
		else
			internal_test_label.push_back((int)y1[i]);
	}
	get_avgs(x, avx);
	get_sds(x, sdx);
	avy = avg(y);
	impute_x(x, avx);
	impute_x(internal_test_x, avx);
	for (double& val : y)
		val -= avy;
	r = vector<double>(y);
	scores.resize(nsamples);
	internal_test_scores.resize(internal_test_nsamples);
	ftr_grids.resize(nftrs);
	is_categorial.resize(nftrs);
	frequent_ftr_vals.resize(nftrs);
	bin_ftr_indexes.resize(params.num_iterations);
	bin_ftr_avg_sd_beta.resize(params.num_iterations);
	mark_grids_and_frequent_vals();
	cout << endl << "Deep Bit:" << endl << "The number of binary featrues is:" << num_bin_ftrs << endl;
	//print_log();
}

void MedDeepBit::impute_x(vector<vector<double>>& x, const vector<double>& avx) {
	int nftrs = (int)x.size();
	int nsamples = (int)x[0].size();
	for (int j = 0; j < nftrs; j++)
		for (int i = 0; i < nsamples; i++)
			if (x[j][i] == MED_MAT_MISSING_VALUE)
				x[j][i] = avx[j];
}

int MedDeepBit::Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
	return Learn(x, y, nsamples, nftrs);
}

int MedDeepBit::Learn(float *x, float *y, int nsamples, int nftrs) {
	init(x, y, nsamples, nftrs);
	train_model();
	return 1;
}

void MedDeepBit::init_defaults() {
	classifier_type = MODEL_DEEP_BIT;
	transpose_for_learn = false;
	transpose_for_predict = false;
	normalize_for_learn = false;
	normalize_for_predict = false;
	normalize_y_for_learn = false;
}

int MedDeepBit::set_params(map<string, string>& mapper) {
	for (auto& entry : mapper) {
		string field = entry.first;
		//! [MedDeepBit::init]
		if (field == "max_depth") params.max_depth = stoi(entry.second);
		else if (field == "num_iterations") params.num_iterations = stoi(entry.second);
		else if (field == "num_ftrs_per_round") params.num_ftrs_per_round = stoi(entry.second);
		else if (field == "num_vals_to_be_categorial") params.num_vals_to_be_categorial = stoi(entry.second);
		else if (field == "nparts_auc") params.nparts_auc = stoi(entry.second);
		else if (field == "niter_coordinate_descent") params.niter_coordinate_descent = stoi(entry.second);
		else if (field == "niter_auc_gitter") params.niter_auc_gitter = stoi(entry.second);
		else if (field == "internal_test_ratio") params.internal_test_ratio = stoi(entry.second);
		else if (field == "fraction_auc") params.fraction_auc = stof(entry.second);
		else if (field == "grid_fraction") params.grid_fraction = stof(entry.second);
		else if (field == "min_fraction_zeros_ones") params.grid_fraction = stof(entry.second);
		else if (field == "frac_continuous_frequent") params.frac_continuous_frequent = stof(entry.second);
		else if (field == "frac_categorial_frequent") params.frac_categorial_frequent = stof(entry.second);
		else if (field == "lambda") params.lambda = stof(entry.second);
		else if (field == "min_cor_bin_ftr") params.min_cor_bin_ftr = stof(entry.second);

		else MLOG("Unknonw parameter \'%s\' for XGB\n", field.c_str());
		//! [MedDeepBit::init]
	}
	return 0;
}

bool MedDeepBit::is_viable_01_ratios(int count0, int count1, int count_pos0, int count_pos1) {
	int min_num_zeros_ones = (int)(params.min_fraction_zeros_ones * nsamples);
	if (min(count1 + 9 * count_pos1, count0 + 9 * count_pos0) < min_num_zeros_ones)
		return false;
	return true;
}

bool MedDeepBit::is_bin_ftr_valid(const vector<char>& bin_ftr) {
	vector<double> ftr(bin_ftr.size());
	for (int i = 0; i < ftr.size(); i++)
		ftr[i] = bin_ftr[i];
	double ftr_cor = medial::performance::pearson_corr_without_cleaning(ftr, y);
	return(abs(ftr_cor) >= params.min_cor_bin_ftr);
}

void MedDeepBit::make_bin_ftrs(int j, const vector<double>& vals, bool is_categorial) {
	int nvals = (int)vals.size();
	vector<vector<char>> temp_bin_ftrs(2 * nvals);
	vector<tuple<int, int, bool, bool>> temp_ftrs_map(2 * nvals);
#pragma omp parallel for schedule(dynamic)
	for (int k = 0; k < nvals; k++) {
		vector<char> bin_ftr(nsamples);
		if (is_categorial)
			get_categorial_bin_ftr(x[j], vals[k], true, bin_ftr);
		else
			get_quant_bin_ftr(x[j], vals[k], true, bin_ftr);
		if (is_bin_ftr_valid(bin_ftr)) {
			temp_bin_ftrs[2 * k] = bin_ftr;
			bin_ftr = vector<char>(nsamples, 0);
			if (is_categorial)
				get_categorial_bin_ftr(x[j], vals[k], false, bin_ftr);
			else
				get_quant_bin_ftr(x[j], vals[k], false, bin_ftr);
			temp_bin_ftrs[2 * k + 1] = bin_ftr;
			temp_ftrs_map[2 * k] = tuple<int, int, bool, bool>(j, k, true, is_categorial);
			temp_ftrs_map[2 * k + 1] = tuple<int, int, bool, bool>(j, k, false, is_categorial);
			num_bin_ftrs += 2;
		}
	}
	for (int k = 0; k < 2 * nvals; k++) {
		if (temp_bin_ftrs[k].size() == nsamples) {
			bin_ftrs.push_back(temp_bin_ftrs[k]);
			bin_ftrs_map.push_back(temp_ftrs_map[k]);
		}
	}
}

int MedDeepBit::step_function(int step, int i) {
	int n = nsamples, s = step;
	return -16 * s / (n * n) * (i * i) + 16 * s / n * i + s;
}

void MedDeepBit::mark_grids_and_frequent_vals() {
	num_bin_ftrs = 0;
	int num_vals_to_be_frequent;
	cout << "Constructing binary features." << endl;
	for (int j = 0; j < nftrs; j++) {
		cout << j << " ";
		frequent_ftr_vals[j] = vector<double>();
		vector<double> col{ x[j] };
		map<double, int> count;
		for (double val : col) {
			if (count.count(val) == 0)
				count[val] = 1;
			else
				count[val]++;
		}
		if (count.size() <= params.num_vals_to_be_categorial) {
			is_categorial[j] = true;
			num_vals_to_be_frequent = (int)round(params.frac_categorial_frequent * nsamples);
		}
		else {
			is_categorial[j] = false;
			num_vals_to_be_frequent = (int)round(params.frac_continuous_frequent * nsamples);
		}
		for (auto const& e : count)
			if (e.second > num_vals_to_be_frequent)
				frequent_ftr_vals[j].push_back(e.first);
		frequent_ftr_vals[j].push_back(avx[j]);
		if (frequent_ftr_vals[j].size() > 0)
			sort(frequent_ftr_vals[j].begin(), frequent_ftr_vals[j].end());
		make_bin_ftrs(j, frequent_ftr_vals[j], true);
		if (!is_categorial[j]) {
			vector<double> grid;
			vector<double> sorted_col(col);
			sort(sorted_col.begin(), sorted_col.end());
			int step = (int)round(params.grid_fraction * nsamples);
			step += (int)round(rand() % (int)(params.grid_fraction * nsamples / 2)) - (int)round(params.grid_fraction * nsamples) / 4; //Different jumps for each feature.
			//cout << "jump = " << step << endl;
			int start = (int)round((rand() % nsamples) * params.grid_fraction);
			grid.push_back(sorted_col[start]);
			for (int i = start + step; i < nsamples; i += step_function(step, i)) {
				int k = i;
				while (k < nsamples && sorted_col[k] == grid.back()) k++;
				if (k < nsamples)
					grid.push_back(sorted_col[k]);
				i = k;
			}
			make_bin_ftrs(j, grid, false);
			ftr_grids[j] = grid;
		}
	}
}

void MedDeepBit::print_log() {
	for (int j = 0; j < nftrs; j++) {
		if (is_categorial[j]) {
			cout << "The feature: " << j << " is used as a categorial feature with " << frequent_ftr_vals[j].size() << " values:" << endl;
			for (double val : frequent_ftr_vals[j])
				cout << val << ", ";
			cout << endl << endl;
		}
		else {
			cout << "The feature: " << j << " is used as a continuous feature with the following grid:" << endl;
			for (double val : ftr_grids[j])
				cout << val << ", ";
			cout << endl;
			cout << "The common vals we use as categories are:" << endl;
			for (double val : frequent_ftr_vals[j])
				cout << val << ", ";
			cout << endl << endl;
		}
	}
}


void MedDeepBit::get_categorial_bin_ftr(const vector<double>& col, double val, bool direction, vector<char>& bin_ftr) {
	if (direction == true) {
		for (int i = 0; i < nsamples; i++)
			if (col[i] == val)
				bin_ftr[i] = 1;
	}
	else {
		for (int i = 0; i < nsamples; i++)
			if (col[i] != val)
				bin_ftr[i] = 1;
	}
}

void MedDeepBit::get_quant_bin_ftr(const vector<double>& col, double val, bool direction, vector<char>& bin_ftr) {
	if (direction == true) {
		for (int i = 0; i < nsamples; i++) {
			if (col[i] >= val || (col[i] == MED_MAT_MISSING_VALUE && 0 >= val))
				bin_ftr[i] = 1;
		}
	}
	else {
		for (int i = 0; i < nsamples; i++)
			if (col[i] < val || (col[i] == MED_MAT_MISSING_VALUE && 0 < val))
				bin_ftr[i] = 1;
	}
}

int MedDeepBit::get_categorial_bit(double x_val, double val, bool direction) const {
	if ((direction == true && x_val == val) || (direction == false && x_val != val)) return 1;
	return 0;
}

int MedDeepBit::get_quant_bit(double x_val, double val, bool direction) const {
	if ((direction == true && x_val >= val) || (direction == false && x_val < val)) return 1;
	return 0;
}

double MedDeepBit::get_ftr_score(const vector<char>& bin_ftr) {
	double s0 = 0, s1 = 0;
	int count0 = 0, count1, count_pos0 = 0, count_pos1 = 0;
	for (int i = 0; i < nsamples; i++) {
		if (bin_ftr[i] == 0) {
			count0++;
			s0 += r[i];
			count_pos0 += label[i];
		}
		else {
			s1 += r[i];
			count_pos1 += label[i];
		}
	}
	count1 = nsamples - count0;
	if (min(count0, count1) < params.min_fraction_zeros_ones * nsamples)
		//if (!is_viable_01_ratios(count0, count1, count_pos0, count_pos1))
		return 0;
	double p = ((double)count1) / nsamples;
	double v1 = ((double)1) / nsamples * ((1 - p) / p * s1 - s0);
	double v0 = ((double)1) / nsamples * (p / (1 - p) * s0 - s1);
	double ans = 0;
	for (int i = 0; i < nsamples; i++) {
		double v = v1;
		double w = 1;
		if (label[i] == 1) w = 1;
		if (bin_ftr[i] == 0) v = v0;
		ans += w * (r[i] - v) * (r[i] - v);
	}
	ans /= nsamples;
	//return abs(av0 - av1);
	//cout << 1e10 - ans;
	return 1e10 - ans;
}

//double MedDeepBit::get_ftr_score(const vector<char>& bin_ftr) {
//	double s0 = 0, s1 = 0;
//	int count0 = 0, count1, count_pos0 = 0, count_pos1 = 0;
//	for (int i = 0; i < nsamples; i++) {
//		if (bin_ftr[i] == 0) {
//			count0++;
//			s0 += r[i];
//			count_pos0 += label[i];
//		}
//		else {
//			s1 += r[i];
//			count_pos1 += label[i];
//		}
//	}
//	count1 = nsamples - count0;
//	//if (min(count0, count1) < min_fraction_zeros_ones * nsamples)
//	if (!is_viable_01_ratios(count0, count1, count_pos0, count_pos1))
//		return 0;
//	double av, std, b;
//	vector<double> scores1(nsamples), r1(nsamples);
//	calc_bin_ftr_scores(bin_ftr, av, std, b, scores1, r1, true);
//	double ans = 0;
//	for (int i = 0; i < nsamples; i++)	ans += r1[i] * r1[i];
//	ans /= nsamples;
//	//cout << ans << endl;
//	//return abs(av0 - av1);
//	return 1e10 - ans;
//}



void MedDeepBit::mult_bin_ftrs(const vector<char>& ftr1, vector<char>& ftr2) {
	for (int i = 0; i < nsamples; i++)
		ftr2[i] *= ftr1[i];
}

void MedDeepBit::get_bin_ftr(int bin_ftr_index, vector<char>& bin_ftr) {
	int j, k;
	bool direction, is_frequent;
	tie(j, k, direction, is_frequent) = bin_ftrs_map[bin_ftr_index];
	if (is_categorial[j] || is_frequent) {
		double val = frequent_ftr_vals[j][k];
		get_categorial_bin_ftr(x[j], val, direction, bin_ftr);
	}
	else {
		double val = ftr_grids[j][k];
		get_quant_bin_ftr(x[j], val, direction, bin_ftr);
	}
}


void MedDeepBit::print_ftr_characteristics(int index) {
	int j, k;
	bool direction, is_frequent;
	tie(j, k, direction, is_frequent) = bin_ftrs_map[index];

	if (is_categorial[j] || is_frequent)
		cout << j << " " << direction << " " << frequent_ftr_vals[j][k] << " " << is_frequent << endl;
	else
		cout << j << " " << direction << " " << ftr_grids[j][k] << " " << is_frequent << endl;
}

double MedDeepBit::avg(const vector<char>& vec) {
	double ans = 0;
	for (char val : vec) ans += val;
	ans /= vec.size();
	return(ans);
}

double MedDeepBit::sd(const vector<char>& binary_vec) {
	double av = avg(binary_vec);
	return sqrt(av * (1 - av));
	//return (av * (1 - av));
}

void MedDeepBit::get_normalized_col(const vector<char>& col, vector<double>& normalized_col, double& av, double& std) {
	int nrow = (int)col.size();
	av = avg(col);
	std = sd(col);
	if (std > 0)
		for (int i = 0; i < nrow; i++)
			normalized_col[i] = (col[i] - av) / std;
}

double MedDeepBit::avg(const vector<double>& vec) {
	double ans = 0;
	for (double val : vec) ans += val;
	ans /= vec.size();
	return(ans);
}

double MedDeepBit::sd(const vector<double>& vec) {
	double sum = 0;
	double av = avg(vec);
	for (double val : vec)
		sum += pow((val - av), 2);
	return sqrt(sum / vec.size());
}

void MedDeepBit::get_col_without_na(const vector<double>& col, vector<double>& col_without_na) {
	for (int i = 0; i < col.size(); i++)
		if (col[i] != MED_MAT_MISSING_VALUE)
			col_without_na.push_back(col[i]);
}

void MedDeepBit::get_avgs(const vector<vector<double>>& x, vector<double>& avx) {
	for (const vector<double>& col : x) {
		vector<double> col_without_na;
		get_col_without_na(col, col_without_na);
		double av = avg(col_without_na);
		avx.push_back(av);
		cout << av << " ";
	}
	cout << endl;
}
void MedDeepBit::get_sds(const vector<vector<double>>& x, vector<double>& sdx) {
	for (const vector<double>& col : x) {
		vector<double> col_without_na;
		get_col_without_na(col, col_without_na);
		sdx.push_back(sd(col_without_na));
	}
}

void MedDeepBit::transpose(const vector<vector<double>>& before, vector<vector<double>>& after) {
	int nrow_old = (int)before.size();
	int ncol_old = (int)before[0].size();
	after.resize(ncol_old);
	cout << "ftr_num = " << ncol_old << " sample_num = " << nrow_old << endl;
	for (int j = 0; j < ncol_old; j++)
		after[j].resize(nrow_old);
	for (int i = 0; i < nrow_old; i++) {
		for (int j = 0; j < ncol_old; j++)
			after[j][i] = before[i][j];
	}
}


double MedDeepBit::perform_lasso_iteration(const vector<double>& xk_train, const vector<double>& r, double lambda, double alpha) {
	double bk_hat = 0;
	if (alpha < 1) {
		int nsamples = (int)xk_train.size();
		for (int i = 0; i < nsamples; i++)
			bk_hat += r[i] * xk_train[i];
		bk_hat /= nsamples;
		if (bk_hat - lambda * (1 - alpha) > 0)
			bk_hat -= lambda * (1 - alpha);
		else if (bk_hat + lambda * (1 - alpha) < 0)
			bk_hat += lambda * (1 - alpha);
		else
			bk_hat = 0;
		if (bk_hat != 0)
			bk_hat /= (1 + alpha * lambda);
	}
	else {
		int nsamples = (int)xk_train.size();
		for (int i = 0; i < nsamples; i++)
			bk_hat += r[i] * xk_train[i];
		bk_hat /= nsamples;
		bk_hat /= (1 + lambda);
	}
	return bk_hat;

}


void MedDeepBit::calc_bin_ftr_scores(const vector<char>& bin_ftr, double& av, double& std, double& b, vector<double>& scores1, vector<double>& r1, bool is_full_step) {
	vector<double> normalized_ftr(nsamples);
	get_normalized_col(bin_ftr, normalized_ftr, av, std);
	if (is_full_step)
		b = perform_lasso_iteration(normalized_ftr, r, 0, 1);
	else
		b = perform_lasso_iteration(normalized_ftr, r, params.lambda, 1);
	for (int i = 0; i < nsamples; i++) {
		scores1[i] = scores[i] + b * normalized_ftr[i];
		r1[i] = r[i] - b * normalized_ftr[i];
	}
}

void MedDeepBit::gen_random_indexes(vector<int>& random_indexes) {
	for (int bin_ftr_num = 0; bin_ftr_num < params.num_ftrs_per_round; bin_ftr_num++)
		random_indexes[bin_ftr_num] = rand() % num_bin_ftrs;
	//random_indexes[bin_ftr_num] = count_bin_ftrs++ % num_bin_ftrs;
}

void MedDeepBit::score_random_ftrs(vector<double>& ftr_scores, const vector<int>& random_indexes, const vector<char>& final_bin_ftr) {
#pragma omp parallel for schedule(dynamic)
	for (int bin_ftr_num = 0; bin_ftr_num < params.num_ftrs_per_round; bin_ftr_num++) {
		int bin_ftr_index = random_indexes[bin_ftr_num];
		vector<char> bin_ftr(bin_ftrs[bin_ftr_index]);
		mult_bin_ftrs(final_bin_ftr, bin_ftr);
		ftr_scores[bin_ftr_num] = get_ftr_score(bin_ftr);
	}
}

void MedDeepBit::train_model() {
	double cur_auc = 0;
	for (int it = 0; it < params.num_iterations; it++) {
		vector<char> final_bin_ftr(nsamples, 1);
		bin_ftr_indexes[it] = vector<int>(params.max_depth);
		double best_score = 0;
		for (int depth = 0; depth < params.max_depth; depth++) {
			vector<double> ftr_scores(params.num_ftrs_per_round);
			vector<int> random_indexes(params.num_ftrs_per_round);
			gen_random_indexes(random_indexes);
			score_random_ftrs(ftr_scores, random_indexes, final_bin_ftr);
			int best_ftr_index = -1;
			bool improved_score = false;
			for (int i = 0; i < ftr_scores.size(); i++) {
				if (ftr_scores[i] > best_score) {
					improved_score = true;
					best_score = ftr_scores[i];
					best_ftr_index = random_indexes[i];
				}
			}
			if (improved_score) {
				mult_bin_ftrs(bin_ftrs[best_ftr_index], final_bin_ftr);
				bin_ftr_indexes[it][depth] = best_ftr_index;
				//print_ftr_characteristics(best_ftr_index);
			}
			else {
				//cout << "Feature proportion too small or can't improve bin feature, finishing the iteration. depth = " << depth << endl;
				for (; depth < params.max_depth; depth++)
					bin_ftr_indexes[it][depth] = -1;
			}
		}
		double av, std, b;
		vector<double> scores1(nsamples), r1(nsamples);
		calc_bin_ftr_scores(final_bin_ftr, av, std, b, scores1, r1);
		double new_auc = 0;
		if (abs(b) > 0)
			new_auc = special_auc(scores1, label, false, params.nparts_auc);
		if (new_auc > cur_auc * params.fraction_auc) {
			cur_auc = new_auc;
			scores = scores1;
			r = r1;
			bin_ftr_avg_sd_beta[it] = vector<double>{ av, std, b };
			if (it % 10 == 0) {
				cout << "Iteration #" << it << endl;
				cout << "p = " << av << endl;
				cout << "Train AUC = " << cur_auc << endl;
				if (internal_test_nsamples > 0) {
#pragma omp parallel for schedule(dynamic)
					for (int i = 0; i < internal_test_nsamples; i++) 
						internal_test_scores[i] = predict_sample(internal_test_transposed_x[i], it);
					cout << "Test AUC = " << special_auc(internal_test_scores, internal_test_label) << endl << endl;
				}
			}
		}
		else {
			it--;
			//cout << "Iteration #" << it << ", failed to improve the auc. " << endl;
		}
	}
	do_coordinate_descent(params.num_iterations - 1);
	do_auc_gittering();
}

void MedDeepBit::get_bin_ftr_of_it(int it, vector<char>& bin_ftr) {
	bin_ftr.resize(nsamples, 1);
	for (int depth = 0; depth < params.max_depth; depth++) {
		int bin_ftr_index = bin_ftr_indexes[it][depth];
		if (bin_ftr_index != -1)
			mult_bin_ftrs(bin_ftrs[bin_ftr_index], bin_ftr);
		else
			break;
	}
}

void MedDeepBit::do_coordinate_descent(int num_iterations_descent) {
	int niter_descent = 0;
	double cur_auc = special_auc(scores, label, true);
	for (int descent_i = 0; descent_i < niter_descent; descent_i++) {
		for (int it = 0; it <= num_iterations_descent; it++) {
			vector<char> bin_ftr;
			get_bin_ftr_of_it(it, bin_ftr);
			double av, std, b;
			vector<double> r1(nsamples), scores1(nsamples);
			calc_bin_ftr_scores(bin_ftr, av, std, b, scores1, r1);
			double new_auc = special_auc(scores1, label, true);
			if (new_auc > cur_auc) {
				cur_auc = new_auc;
				scores = scores1;
				r = r1;
				bin_ftr_avg_sd_beta[it][2] = b;
			}
			if (it % 10 == 0)
				cout << "Coordinate descent iteration #" << it << ", the special auc is: " << cur_auc << endl << "The auc is: " << special_auc(scores, label, false) << endl;
		}
	}
}

void MedDeepBit::do_auc_gittering() {
	double cur_auc = special_auc(scores, label);
	int niter_gitter = 0;
	double step = 0.1, from = 0.6;
	cout << "Starting AUC gittering." << endl;
	for (int auc_it = 0; auc_it < niter_gitter; auc_it++) {
		for (int it = 0; it < params.num_iterations; it++) {
			vector<char> bin_ftr;
			get_bin_ftr_of_it(it, bin_ftr);
			vector<double> normalized_ftr(nsamples);
			double av, std;
			get_normalized_col(bin_ftr, normalized_ftr, av, std);
			double b = bin_ftr_avg_sd_beta[it][2];
			vector<double> grid;
			for (double pt = b * from; abs(pt) <= abs((2 - from) * b); pt += b * step)
				//for (double pt = b * from; abs(pt) <= abs(b)*1.01; pt += b * step)
				grid.push_back(pt);
			vector<double> grid_aucs(grid.size());

#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < grid.size(); j++) {
				double b0 = grid[j];
				vector<double> scores1(scores);
				for (int i = 0; i < nsamples; i++)
					scores1[i] += (b0 - b) * normalized_ftr[i];
				grid_aucs[j] = special_auc(scores1, label);
			}

			double best_auc = cur_auc;
			int best_j = -1;
			for (int j = 0; j < grid.size(); j++) {
				if (grid_aucs[j] > best_auc) {
					best_auc = grid_aucs[j];
					best_j = j;
				}
			}
			double best_b;
			if (best_j >= 0) {
				cur_auc = grid_aucs[best_j];
				best_b = grid[best_j];
				bin_ftr_avg_sd_beta[it][2] = best_b;
				for (int i = 0; i < nsamples; i++) {
					scores[i] += (best_b - b) * normalized_ftr[i];
					r[i] -= (best_b - b) * normalized_ftr[i];
				}
			}
			if (it % 10 == 0) {
				cout << "Auc iteration#" << auc_it << ", deep bit iteration# " << it << endl;
				cout << "auc gitter it# " << it << ", auc = " << cur_auc << endl;
			}
		}
	}
}

double MedDeepBit::predict_sample(const vector<double>& x, int niter) const {
	double ans = 0;
	for (int it = 0; it < niter; it++) {
		int final_bit = 1;
		for (int depth = 0; depth < params.max_depth; depth++) {
			int bit;
			int bin_ftr_index = bin_ftr_indexes[it][depth];
			if (bin_ftr_index == -1) break;
			int j, k;
			bool direction, is_frequent;
			tie(j, k, direction, is_frequent) = bin_ftrs_map[bin_ftr_index];
			double x_val = x[j];
			if (x_val == MED_MAT_MISSING_VALUE)
				x_val = avx[j];
			if (is_categorial[j] || is_frequent) {
				double val = frequent_ftr_vals[j][k];
				bit = get_categorial_bit(x_val, val, direction);
			}
			else {
				double val = ftr_grids[j][k];
				bit = get_quant_bit(x_val, val, direction);
			}
			final_bit *= bit;
			if (final_bit == 0) break;
		}
		double av = bin_ftr_avg_sd_beta[it][0], std = bin_ftr_avg_sd_beta[it][1], b = bin_ftr_avg_sd_beta[it][2];
		if (std > 0)
			ans += b * (final_bit - av) / std;
	}
	return ans;
}

double MedDeepBit::predict_sample(const vector<double>& x) const {
	return predict_sample(x, params.num_iterations);
}


double MedDeepBit::get_normalized_val(double x_val, int j) {
	double std = sdx[j];
	if (std == 0 || x_val == MED_MAT_MISSING_VALUE) return 0;
	return (x_val - avx[j]) / std;
}

void MedDeepBit::predict(const vector<vector<double>>& x, vector<double>& scores) const {
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < x.size(); i++) {
		scores[i] = predict_sample(x[i]);
	}
}

int MedDeepBit::Predict(float *x, float *&preds, int nsamples, int nftrs) const {
	if (preds == NULL)
		preds = new float[nsamples];
	vector<vector<double>> x1(nsamples);
	for (int i = 0; i < nsamples; i++) {
		x1[i].resize(nftrs);
		for (int j = 0; j < nftrs; j++) 
			x1[i][j] = x[i * nftrs + j];
	}
	vector<double> scores(nsamples);
	predict(x1, scores);
	for (int i = 0; i < nsamples; i++)
		preds[i] = (float)scores[i];
	return 1;
}


struct sort_special {
	bool operator()(const pair<double, int>& left, const pair<double, int>& right) {
		return left.first > right.first;
	}
};

//Returns the minimal result by parts.
double MedDeepBit::special_auc(const vector<double>& all_predictions, const vector<int>& all_label, bool is_weighted, int nparts) {
	vector<double> results(nparts);
	for (int part = 0; part < nparts; part++) {
		int first = (int)all_predictions.size() / nparts * part, last = (int)all_predictions.size() / nparts * (part + 1) - 1;
		vector<double> predictions(all_predictions.begin() + first, all_predictions.begin() + last);
		vector<int> label(all_label.begin() + first, all_label.begin() + last);
		int len = (int)label.size();
		if (predictions.size() != len) {
			cout << "Can't calculate auc for vectors of different sizes!";
			exit(0);
		}
		vector<pair<double, int>> pred_and_label(len);
		double sum[2] = { 0, 0 };
		for (int i = 0; i < len; i++) {
			pred_and_label[i] = pair<double, int>{ predictions[i], label[i] };
			sum[label[i]] += 1;
		}
		stable_sort(pred_and_label.begin(), pred_and_label.end(), sort_special());
		double ones = 0, auc_cnt = 0;
		for (int i = 0; i < len; ) {
			double val = pred_and_label[i].first;
			double delta_ones = 0;
			int j = i;
			for (; j < len && pred_and_label[j].first == val; j++)
				if (pred_and_label[j].second == 1)
					delta_ones++;
			/*double dist = (double)(len - ((i + j) / 2));
			double delta = (dist / len) * (j - i);
			auc_cnt += delta * (ones + delta_ones / 2);*/
			double dist = (double)(len - ((i + j) / 2)) / len;
			double dx;
			if (is_weighted)
				dx = dist * ((j - i) - delta_ones);
			else
				dx = ((j - i) - delta_ones);
			double dy = ones + delta_ones / 2;
			auc_cnt += dx * dy;
			i = j;
			ones += delta_ones;
		}
		if (sum[1] == 0 || sum[0] == 0)
			results[part] = 0;
		results[part] = auc_cnt / (sum[0] * sum[1]);
	}
	double result = 1;
	for (int part = 0; part < nparts; part++) {
		if (results[part] < result)	result = results[part];
	}
	return result;
}

size_t MedDeepBit::get_size() {
	size_t ptr = 0;
	ptr += params.get_size();
	ptr += sizeof(int);
	ptr += sizeof(int);
	ptr += sizeof(int);
	ptr += sizeof(double);
	ptr += nftrs * sizeof(double);
	//memcpy(blob + ptr, &(ftr_names[0]), nftrs * sizeof(string)); ptr += nftrs * sizeof(string);
	//ptr += MedSerialize::serialize(blob + ptr, ftr_names);
	ptr += nftrs * sizeof(int);
	for (int j = 0; j < nftrs; j++) {
		int grid_size = (int)ftr_grids[j].size();
		ptr += sizeof(int);
		ptr += grid_size * sizeof(double);
		int nfreq_vals = (int)frequent_ftr_vals[j].size();
		ptr += sizeof(int);
		ptr += nfreq_vals * sizeof(double);
	}

	for (int i = 0; i < num_bin_ftrs; i++) {
		int j, k;
		bool direction, is_frequent;
		tie(j, k, direction, is_frequent) = bin_ftrs_map[i];
		ptr += sizeof(int);
		ptr += sizeof(int);
		ptr += sizeof(bool);
		ptr += sizeof(bool);
	}

	for (int it = 0; it < params.num_iterations; it++) {
		ptr += params.max_depth * sizeof(int);
	}

	for (int i = 0; i < params.num_iterations; i++)
		ptr += 3 * sizeof(double);
	return ptr;
}

void MedDeepBit::print_model(FILE *fp, const string& prefix) const {
	fprintf(fp, "%s %f %d %d %d\n", prefix.c_str(), params.lambda, params.num_ftrs_per_round, params.num_iterations, params.max_depth);
	fprintf(fp, "%s nftrs = %d num-bin-ftrs = %d avy = %lf\n", prefix.c_str(), nftrs, num_bin_ftrs, avy);
	
	fprintf(fp, "%s ", prefix.c_str());
	for (double val : avx) 
		fprintf(fp, "%lf ", val);
	fprintf(fp, "\n");

	fprintf(fp, "%s ", prefix.c_str());
	for (int val : is_categorial)
		fprintf(fp, "%d ", ((bool)(val != 0)));
	fprintf(fp, "\n");
}

size_t MedDeepBit::serialize(unsigned char *blob) {
	size_t ptr = 0;
	ptr += params.serialize(blob);
	memcpy(blob + ptr, &classifier_type, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &nftrs, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &num_bin_ftrs, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(avy), sizeof(double)); ptr += sizeof(double);
	memcpy(blob + ptr, &(avx[0]), nftrs * sizeof(double)); ptr += nftrs * sizeof(double);
	//memcpy(blob + ptr, &(ftr_names[0]), nftrs * sizeof(string)); ptr += nftrs * sizeof(string);
	//ptr += MedSerialize::serialize(blob + ptr, ftr_names);
	memcpy(blob + ptr, &(is_categorial[0]), nftrs * sizeof(int)); ptr += nftrs * sizeof(int);
	for (int j = 0; j < nftrs; j++) {
		int grid_size = (int)ftr_grids[j].size();
		memcpy(blob + ptr, &(grid_size), sizeof(int)); ptr += sizeof(int);
		memcpy(blob + ptr, &(ftr_grids[j][0]), grid_size * sizeof(double)); ptr += grid_size * sizeof(double);
		int nfreq_vals = (int)frequent_ftr_vals[j].size();
		memcpy(blob + ptr, &(nfreq_vals), sizeof(int)); ptr += sizeof(int);
		memcpy(blob + ptr, &(frequent_ftr_vals[j][0]), nfreq_vals * sizeof(double)); ptr += nfreq_vals * sizeof(double);
	}

	for (int i = 0; i < num_bin_ftrs; i++) {
		int j, k;
		bool direction, is_frequent;
		tie(j, k, direction, is_frequent) = bin_ftrs_map[i];
		memcpy(blob + ptr, &(j), sizeof(int)); ptr += sizeof(int);
		memcpy(blob + ptr, &(k), sizeof(int)); ptr += sizeof(int);
		memcpy(blob + ptr, &(direction), sizeof(bool)); ptr += sizeof(bool);
		memcpy(blob + ptr, &(is_frequent), sizeof(bool)); ptr += sizeof(bool);
	}

	for (int it = 0; it < params.num_iterations; it++) {
		memcpy(blob + ptr, &(bin_ftr_indexes[it][0]), params.max_depth * sizeof(int)); ptr += params.max_depth * sizeof(int);
	}

	for (int i = 0; i < params.num_iterations; i++) {
		memcpy(blob + ptr, &(bin_ftr_avg_sd_beta[i][0]), 3 * sizeof(double)); ptr += 3 * sizeof(double);
	}
	//print_model();
	return ptr;
}

size_t MedDeepBit::deserialize(unsigned char *blob) {
	size_t ptr = 0;
	ptr += params.deserialize(blob);
	memcpy(&classifier_type, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&nftrs, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&num_bin_ftrs, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(avy), blob + ptr, sizeof(double)); ptr += sizeof(double);
	avx.resize(nftrs); memcpy(&(avx[0]), blob + ptr, nftrs * sizeof(double)); ptr += nftrs * sizeof(double);
	//ptr += MedSerialize::deserialize(blob + ptr, ftr_names);
	is_categorial.resize(nftrs); memcpy(&(is_categorial[0]), blob + ptr, nftrs * sizeof(int)); ptr += nftrs * sizeof(int);
	ftr_grids.resize(nftrs);
	frequent_ftr_vals.resize(nftrs);
	for (int j = 0; j < nftrs; j++) {
		int grid_size; memcpy(&(grid_size), blob + ptr, sizeof(int)); ptr += sizeof(int);
		ftr_grids[j].resize(grid_size); memcpy(&(ftr_grids[j][0]), blob + ptr, grid_size * sizeof(double)); ptr += grid_size * sizeof(double);
		int nfreq_vals;	memcpy(&(nfreq_vals), blob + ptr, sizeof(int)); ptr += sizeof(int);
		frequent_ftr_vals[j].resize(nfreq_vals); 
		memcpy(&(frequent_ftr_vals[j][0]), blob + ptr, nfreq_vals * sizeof(double));
		ptr += nfreq_vals * sizeof(double);
	}
	bin_ftrs_map.resize(num_bin_ftrs);
	for (int i = 0; i < num_bin_ftrs; i++) {
		int j, k;
		bool direction, is_frequent;
		memcpy(&(j), blob + ptr, sizeof(int)); ptr += sizeof(int);
		memcpy(&(k), blob + ptr, sizeof(int)); ptr += sizeof(int);
		memcpy(&(direction), blob + ptr, sizeof(bool)); ptr += sizeof(bool);
		memcpy(&(is_frequent), blob + ptr, sizeof(bool)); ptr += sizeof(bool);
		bin_ftrs_map[i] = tuple<int, int, bool, bool>(j, k, direction, is_frequent);
	}
	bin_ftr_indexes.resize(params.num_iterations);
	for (int it = 0; it < params.num_iterations; it++) {
		bin_ftr_indexes[it].resize(params.max_depth); memcpy(&(bin_ftr_indexes[it][0]), blob + ptr, params.max_depth * sizeof(int)); ptr += params.max_depth * sizeof(int);
	}
	bin_ftr_avg_sd_beta.resize(params.num_iterations);
	for (int i = 0; i < params.num_iterations; i++) {
		bin_ftr_avg_sd_beta[i].resize(3);
		memcpy(&(bin_ftr_avg_sd_beta[i][0]), blob + ptr, 3 * sizeof(double)); ptr += 3 * sizeof(double);
	}
	//print_model();
	return ptr;
}

//----------------------------------------------------------------------------------------------------------------------------------
//MedDeepBitParams:

size_t MedDeepBitParams::get_size() {
	return 8 * sizeof(int) + 7 * sizeof(double);
}

size_t MedDeepBitParams::serialize(unsigned char *blob) {
	size_t ptr = 0;
	memcpy(blob + ptr, &(max_depth), sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(num_iterations), sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(num_ftrs_per_round), sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(num_vals_to_be_categorial), sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(nparts_auc), sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(niter_auc_gitter), sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(niter_coordinate_descent), sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(internal_test_ratio), sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &(fraction_auc), sizeof(double)); ptr += sizeof(double);
	memcpy(blob + ptr, &(grid_fraction), sizeof(double)); ptr += sizeof(double);
	memcpy(blob + ptr, &(min_fraction_zeros_ones), sizeof(double)); ptr += sizeof(double);
	memcpy(blob + ptr, &(frac_continuous_frequent), sizeof(double)); ptr += sizeof(double);
	memcpy(blob + ptr, &(frac_categorial_frequent), sizeof(double)); ptr += sizeof(double);
	memcpy(blob + ptr, &(lambda), sizeof(double)); ptr += sizeof(double);
	memcpy(blob + ptr, &(min_cor_bin_ftr), sizeof(double)); ptr += sizeof(double);
	return ptr;
}

size_t MedDeepBitParams::deserialize(unsigned char *blob) {
	size_t ptr = 0;
	memcpy(&(max_depth), blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(num_iterations), blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(num_ftrs_per_round), blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(num_vals_to_be_categorial), blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(nparts_auc), blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(niter_auc_gitter), blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(niter_coordinate_descent), blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(internal_test_ratio), blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&(fraction_auc), blob + ptr, sizeof(double)); ptr += sizeof(double);
	memcpy(&(grid_fraction), blob + ptr, sizeof(double)); ptr += sizeof(double);
	memcpy(&(min_fraction_zeros_ones), blob + ptr, sizeof(double)); ptr += sizeof(double);
	memcpy(&(frac_continuous_frequent), blob + ptr, sizeof(double)); ptr += sizeof(double);
	memcpy(&(frac_categorial_frequent), blob + ptr, sizeof(double)); ptr += sizeof(double);
	memcpy(&(lambda), blob + ptr, sizeof(double)); ptr += sizeof(double);
	memcpy(&(min_cor_bin_ftr), blob + ptr, sizeof(double)); ptr += sizeof(double);
	return ptr;
}

void MedDeepBit::print(FILE *fp, const string& prefix, int level) const {
	if (level == 0)
		fprintf(fp, "%s: MedDeepBit ()\n", prefix.c_str());
	else
		print_model(fp, prefix);
}

string MedDeepBitParams::to_string() {
	stringstream strm;
	strm << "max_depth:" << max_depth << ",niter:" << num_iterations << "nftrs_round:" << num_ftrs_per_round << ",nvals_categorial:" << num_vals_to_be_categorial << "nparts_auc:" << nparts_auc <<
		"niter_descent:" << niter_coordinate_descent << "niter_auc_gitter:" << niter_auc_gitter << "grid_fraction:" << grid_fraction << "min_fraction_zeros_ones:" << min_fraction_zeros_ones <<
		"frac_continuous:" << frac_continuous_frequent << "frac_categorial" << frac_categorial_frequent << "lambda:" << lambda << "fraction_auc:" << fraction_auc << "internal_test_ratio:" << internal_test_ratio
		<< "min_cor_bin_ftr=" << min_cor_bin_ftr;
	return strm.str();
}
