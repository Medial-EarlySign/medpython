#include "SimpleButWhy.h"
#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_MEDFEAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//-------------------------------------------------------------------------------------------------
// this method can be used whenever we get contributions summarizing the additive add of each feature to the score.
// assumes the contributions are positive when contributing positively to case prediction
void SimpleButWhyRes::get_from_contribs(SimpleButWhyParams &bw_params, vector<float> &contribs)
{
	float pos_sum = 0, neg_sum = 0;

	for (auto c : contribs)
		if (c >= 0)
			pos_sum += c;
		else
			neg_sum += c;

	float sum_all = pos_sum + neg_sum;
	float epsilon = 1e-20;
	if (abs(sum_all) < epsilon) sum_all = 1.0f;
	if (abs(pos_sum) < epsilon) pos_sum = 1.0f;
	if (abs(neg_sum) < epsilon) neg_sum = -1.0f;
	vector<pair<int, float>> relative_contrib;

	float s_all = pos_sum - neg_sum;
	for (int j = 0; j<contribs.size(); j++) {
		pair<int, float> p;
		p.first = j;
		if (contribs[j] >= 0)
			p.second = contribs[j] / s_all;
		else
			p.second = contribs[j] / s_all;
		//p.second = contribs[j] / sum_all;
		relative_contrib.push_back(p);
	}

	vector<pair<int, float>> sorted_relative_contrib = relative_contrib;

	auto sort_pairs_lambda = [](const pair<int, float> &a, const pair<int, float> &b)->bool { return abs(a.second) > abs(b.second); };
	sort(sorted_relative_contrib.begin(), sorted_relative_contrib.end(), sort_pairs_lambda);

	MLOG("Top Features: pos_sum %f neg_sum %f\n", pos_sum, neg_sum);
	float acc = 0;
	for (int i = 0; i < 2*bw_params.top_features; i++) {
		int j = sorted_relative_contrib[i].first;
		float c = relative_contrib[j].second;
//		if (c > 0) {
			acc += abs(relative_contrib[j].second);
			MLOG("[%d] rel %6.3f : acc %5.3f : contrib %6.3f : feat %3d : %s\n", i, c, acc, contribs[j], j, bw_params.feature_names[j].c_str());
//		}
	}
/*
	acc = 0;
	for (int i = 1; i <= bw_params.top_features;  i++) {
		int j = sorted_relative_contrib[sorted_relative_contrib.size() - i].first;
		float c = relative_contrib[j].second;
//		if (c > 0) {
			acc += relative_contrib[j].second;
			MLOG("[-%d] rel %f : acc %f : contrib %f : feat %d : %s\n", i, c, acc, contribs[j], j, bw_params.feature_names[j].c_str());
//		}
	}
*/
}
