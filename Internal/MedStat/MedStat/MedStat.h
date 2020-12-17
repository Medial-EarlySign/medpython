//
// MedStat - Statistics utilities
//

#ifndef _MED_STAT_H_
#define _MED_STAT_H_

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include "assert.h"
#include "math.h"

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cstring>

#include "string.h"
#include <MedMat/MedMat/MedMat.h>
#include <bitset>

#define MED_DEFAULT_MISSING_VALUE		-1
#define MED_DEFAULT_MIN_TRIM			-1e9
#define MED_DEFAULT_MAX_TRIM			 1e9

using namespace std;

// general useful routines for various statistics given two vectors (usually score + labels)

namespace medial {
	namespace performance {
		/// <summary> Pearson correlation, options: clean missing values, add weights</summary>
		/// <returns> Pearson's R, -2.0 if cannot calculate. outputs n = number of cleaned values </returns>
		template <typename T> float pearson_corr_without_cleaning(const vector<T> &v1, const vector<T> &v2, const vector<float> *weights = NULL);
		template <typename T> float pearson_corr(const vector<T> &v1, const vector<T> &v2, T missing_value, int& n, const vector<float> *weights = NULL);

		/// <summary> Spearman correlation, options: clean missing values, add weights </summary>
		/// <returns> Spearman's R, -2.0 if cannot calculate. outputs n = number of cleaned values </returns>
		template <typename T, typename S> float spearman_corr_without_cleaning(const vector<T> &v1, const vector<S> &v2, const vector<float> *weights = NULL);
		template <typename T> float spearman_corr(const vector<T> &v1, const vector<T> &v2, T missing_val, int &n, const vector<float> *weights = NULL);

		/// <sumary> RMSE of two vectors. options: clean missing values, add weights </summary>
		/// <returns> RMSE. Outputs n = number of cleaned values </returns>
		template <typename T>float rmse_without_cleaning(const vector<T> &preds, const vector<T> &y, const vector<float> *weights = NULL);
		template <typename T>float rmse(const vector<T> &preds, const vector<T> &y, T missing_val, int &n, const vector<float> *weights = NULL);

		/// <sumary> L1 distance of two vectors. options: clean missing values, add weights </summary>
		/// <returns> L1 distance (average of absolute distance). Outputs n = number of cleaned values </returns>
		template <typename T>float L1_dist_without_cleaning(const vector<T> &preds, const vector<T> &y, const vector<float> *weights = NULL);
		template <typename T>float L1_dist(const vector<T> &preds, const vector<T> &y, T missing_val, int &n, const vector<float> *weights = NULL);

		/// <sumary> Relative L1 distance of two vectors. options: clean missing values, add weights </summary>
		/// <returns> Relative L1 distance (average of (abs(preds[i]-y[i])/abs(preds[i])). Outputs n = number of cleaned values </returns>
		template <typename T>float relative_L1_dist_without_cleaning(const vector<T> &preds, const vector<T> &y, const vector<float> *weights = NULL);
		template <typename T>float relative_L1_dist(const vector<T> &preds, const vector<T> &y, T missing_val, int &n, const vector<float> *weights = NULL);

		/// <summary> Kendall rank correlation, options: clean missing values</summary>
		/// <summary> Use is01Vec1,2 to indicated that v2/v2 are 0/1 for faster implementation </summary>
		/// <returns> Kendall's Tau, -2.0 if cannot calculate, output n = number of cleaned values  </returns>
		template <typename T, typename S> double kendall_tau_without_cleaning(const vector<T> &v1, const vector<S> &v2, bool is01Vec1 = false, bool is01Vec2 = false);
		template <typename T> double kendall_tau(const vector<T> &v1, const vector<T> &v2, T missing_value, int &n, bool is01Vec1 = false, bool is01Vec2 = false);

		/// <summary> _q version for a more efficient version if there are only few possibilities fo v1/v2 </summary>
		/// <summary> This version also has optional weights </summary>
		template <typename T, typename S> double kendall_tau_without_cleaning_q(const vector<T> &v1, const vector<S> &v2, const vector<float> *weights = NULL);
		template <typename T, typename S> double kendall_tau_q(const vector<T> &v1, const vector<S> &v2, T missing_val1, S missing_val2, int& n,
			const vector<float> *weights = NULL);

		/// <summary> calculate mutual information between quantized vectors. n = number of non-empty bins. x,y are binnned features </summary>
		/// <returns> mutual information, n= of non-empty bins </returns>
		float mutual_information(const vector<float>& x, const vector<float>& y, int &n);
		/// <summary> calculate mutual information given vectors of counts and co-counts. n = number of non-empty bins </summary>
		/// <returns> mutual information, -1.0 if cannot calculate </returns>
		float mutual_information(vector<int>& xCounts, vector<int>& yCounts, vector<int> coCounts, int n);
		/// <summary> calculate mutual information between quantized vectors. n = number of non-empty bins. x,y are feature bin indexes </summary>
		/// <returns> mutual information, n= of non-empty bins </returns>
		float mutual_information(vector<int>& x, vector<int>& y, int &n);

		/// <summary> calculate AUC <</summary>
		/// <returns> AUC </returns>
		template<typename T> float auc(vector<T> &preds, vector<float> &y);
		/// <summary> calculate AUC on quantized predictions , optinally with weights </summary>
		/// <returns> AUC </returns>
		template<typename T> float auc_q(const vector<T> &preds, const vector<float> &y, const vector<float>* weights = NULL);

		/// <sumary> Collect cnts : TP,FP,FN,TN per positive rate </summary>
		/// <summary> input : vector 'size' of positive rates ; direction>0 indicate that higher score are positive, otherwise, lower scores are positive </summary>
		/// <returns> vector cnts of 4-tuples : {true-positive, false-positive, false-negative, true-nagative} per positive rate = size[i] </returns>
		template<typename T> void get_preds_perf_cnts(vector<T> &preds, vector<float> &y, vector<float> &size, int direction, vector<vector<int>> &cnts);
		/// <summary> translate a 4-tuple {true-positive, false-positive, false-negative, true-nagative}  into performance measures </summary>
		/// <summary> output = sensitivity, specificity, ppv and relative-risk
		void cnts_to_perf(vector<int> &cnt, float &sens, float &spec, float &ppv, float &rr);

		/// <sumary> prediction accuracy. optinally weighted </summary>
		/// <returns> percentage of predictions which are identical to the labels </returns>
		template <typename T> float accuracy(const vector<T> &preds, const vector<float> &y, const vector<float> *weights = NULL);
		template <typename T> float approx_accuracy(const vector<T> &preds, const vector<float> &y, T epsilon, const vector<float> *weights = NULL);

		// Functions for distance correlation
		/// <summary> get the normalized distance matrix of a vector </summary>
		template <typename T> void get_dMatrix(vector<T>& values, MedMat<T>& dMatrix, T missing_value);
		/// <summary> get the variance of distance matrix </summary>
		template <typename T> float get_dVar(MedMat<T>& dMatrix);
		/// <summary> get the covariance of two distance matrices </summary>
		template <typename T> float get_dCov(MedMat<T>& xDistMat, MedMat<T>& yDistMat);

		/// <summary> multi category helpers </summary>
		/// <summary>given multicategory probs (or probs-like) predictions generates a single prediction of the categ with max prob for each sample </summary>
		template <typename T> void multicateg_get_max_pred(vector<T> &probs, int nsamples, int ncateg, vector<float> &max_pred);
		/// <summary>given multicategory probs (or probs-like) predictions generates a single prediction of the categ with average prob for each sample </summary>
		template <typename T> void multicateg_get_avg_pred(vector<T> &probs, int nsamples, int ncateg, vector<T> &avg_pred);
		/// <summary> given multicategory probs (or probs-like) predictions gets the classification error rate and the rms (also for the avg preds) </summary>
		template <typename T> void multicateg_get_error_rate(vector<T> &probs, vector<float> &y, int nsamples, int ncateg, float &err_rate, T &rms, T &avg_rms);

		/// <sumary>given two vectors and a vector of quantization bounds, create the 'confusion' matrix counts </summary>
		template <typename T> void get_quantized_breakdown(vector<T> &preds, vector<T> &y, vector<T> &bounds, MedMat<int> &counts);
		/// <summary>print the 'confusion' matrix counts </summary>
		template <typename T> void print_quantized_breakdown(MedMat<int> &cnt, vector<T> &bounds);
		/// <summary>calculate ICI - calibration index </summary>
		template <typename T> double integrated_calibration_index(const vector<T> &predicted_prob, const vector<float> &y_label, const vector<float>* weights = NULL);
		/// <summary>calculate jaccard similarity   </summary>
		float jaccard_similarity(int val1, int val2);
		/// <summary>calculate jaccard distance   </summary>
		float jaccard_distance(int val1, int val2);
		/// <summary>calculate jaccard distance matrix
		void get_jaccard_matrix(int n, vector<vector<float>>& jaccard_dist);
	}

	namespace stats {
		/// <summary> calculate n x m contigency table chi2 score </summary>
		/// <returns> X^2 </returns>
		double chi2_n_x_m(vector<int> &cnts, int n, int m);
		/// <summary> calculate n x m contigency table chi2 score </summary>
		/// <returns> X^2, exp = expected table </returns>
		double chi2_n_x_m(vector<int> &cnts, int n, int m, vector<double> &exp);

		// Moments
		/// <summary> Mean, options: clean missing values, add weights</summary>
		/// <returns> Mean. outputs n = number of cleaned values. if no values given(left) return missing-val or throw an error if not given </returns>
		template <typename T> double mean_without_cleaning(const vector<T> &v1, const vector<float> *weights = NULL);
		template <typename T> double mean(const vector<T> &v1, T missing_value, int& n, const vector<float> *weights = NULL);

		/// <summary> Standard deviation, options: clean missing values, add weights</summary>
		/// <returns> Standard deviation. outputs n = number of cleaned values.  </returns>
		/// <returns> if no values given(left) return missing-val or throw an error if not given. If a single value is left, return 1.0 </returns>
		template <typename T> double std_without_cleaning(const vector<T> &v, T mean, const vector<float> *weights = NULL);
		template <typename T> double std(const vector<T> &v, T mean, T missing_value, int& n, const vector<float> *weights = NULL);

		/// <summary> Envelopes for mean and std</summary>
		template<typename T> void get_mean_and_std_without_cleaning(const vector<T> &v, T& mean, T& std, const vector<float> *weights = NULL) {
			mean = medial::stats::mean_without_cleaning(v, weights); std = medial::stats::std_without_cleaning(v, mean, weights);
		}
		template<typename T> void get_mean_and_std(const vector<T> &v, T missing_value, int& n, T& mean, T& std, const vector<float> *weights = NULL) {
			mean = medial::stats::mean(v, missing_value, n, weights); std = medial::stats::std(v, mean, missing_value, n, weights);
		}

		template<typename T> void get_z_transform(const vector<T> &v, T missing_value_in_v, T missing_value_z, vector<T> &z);

		/// <summary> Backward compatible version for mean and std</summary>
		void get_mean_and_std(float *values, const float* wgts, int size, float missing_value, float& mean, float&sd, int& n, bool do_missing);

		/// <summary> Median, options : cleaning missing values </summary>
		/// <returns> Median. if no values given(left) return missing-val or throw an error if not given. </returns>
		template<typename T> T median_without_cleaning(vector<T>& v, bool in_place = false);
		template<typename T> T median(vector<T>& v, T missing_value, int& n);

		/// <summary> Most-Common value, options : cleaning missing values </summary>
		/// <returns> Most-Common value. if no values given(left) return missing-val or throw an error if not given. </returns>
		template<typename T> T most_common_without_cleaning(vector<T> &v);
		template<typename T> T most_common(vector<T>& v, T missing_value, int& n);

		/// <summary> Build histogram of values probabilities. option : cleaning missing values </summary>
		template<typename T> void get_histogram_without_cleaning(vector<T>& v, vector<pair<T, float> >& hist, bool in_place = false);
		template<typename T> void get_histogram(vector<T>& v, T missing_value, int& n, vector<pair<T, float> >& hist);
		template<typename T> T sample_from_histogram(vector<pair<T, float> >& hist);


		/// <summary>get a vector of values, a vector of probabilities, and returning a matching vector of values such that Prob(x<=pvals[i])=p[i] </summary>
		/// <summary>currently, there is no in-place version </summary>
		template<class T> void get_percentiles(vector<T> &vals, vector<float> &p, vector<T> &out_pvals, int only_positive_flag = false);
		/// <summary>get a vector of values and weights, a float 0<=q<=1, and returning the   element of the vector that the weights to its left(after sorting) sum to q of the total weights. 
		/// <summary>All this done in average linear time
		template<typename T> T get_quantile(vector<T> vals, vector<float> w, float q);
		/// <summary> gets a vector of values, and checks the best way to round it at the 10^-3 to 10^3 (at 10x steps) range. </summary>
		float get_best_rounding(vector<float>& vals, vector<float>& res, vector<int>& counts, float missing_value = -1);

		double chi_square_table(double grp1_cntrl, double grp1_cases, double grp2_cntrl, double grp2_cases,
			int smooth_balls = 0, float allowed_error = 0);

		/// <summary> standart t_test - assume equal sample size form both groups </summary>
		template<typename T> void t_test(const vector<T> &grp1, const vector<T> &grp2, double &t_value, double &degree_of_freedom, double &p_value);

		/// <summary> standart t_test - assume unequal sample size form both groups, but similar variance - no more than factor 2 </summary>
		template<typename T> void t_test_unequal_sample_size(const vector<T> &grp1, const vector<T> &grp2, double &t_value, double &degree_of_freedom, double &p_value);

		/// <summary> welch t_test when the samples sizes are unequal or the variance of the 2 samples is different </summary>
		template<typename T> void welch_t_test(const vector<T> &grp1, const vector<T> &grp2, double &t_value, double &degree_of_freedom, double &p_value);

		/// <summary> Kullback-Leibler divergence. Epsilon used to correct for zeros
		template<typename T> double KL_divergence(const vector<T>& p, const vector<T>& q, T epsilon = 1e-8);
	}
}

#endif
