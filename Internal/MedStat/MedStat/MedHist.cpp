//
// MedStat
//

#include "MedStat.h"
#include "MedHist.h"

#include <MedUtils/MedUtils/MedGlobalRNG.h>
#include <MedMat/MedMat/MedMat.h>
#include <Logger/Logger/Logger.h>
#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL

extern MedLogger global_logger;

//...............................................................................................................................
// MedHist
//...............................................................................................................................

// Helper: rounding
#define ROUNDING_EPSILON 0.00001
inline float roundf(float val, float rounder) {
	return ((float)((int)((val + (float)ROUNDING_EPSILON) / rounder)) * rounder);
}

void MedHist::clear() {
	hist.clear(); from.clear(); to.clear(); missing_value = (float)MED_DEFAULT_MISSING_VALUE; from_val = 0; to_val = -1; ncells = 200; rounder = -1; positive_only_flag = 1;
	min_percentile = (float)0.001; max_percentile = (float)0.999;
}

int MedHist::get_hist(vector<float> &vals)
{
	if (rounder < 0) rounder = (float)0.01;

	// get boundaries
	if (to_val < from_val) {
		vector<float> p(2), pvals(2);
		p[0] = min_percentile;
		p[1] = max_percentile;
		medial::stats::get_percentiles(vals, p, pvals, positive_only_flag);
		from_val = roundf(pvals[0], rounder);
		to_val = roundf(pvals[1], rounder);
	}

	// reevaluate ncells
	int ncells_in_bound = (int)((to_val - from_val) / rounder);
	if (ncells >= ncells_in_bound)
		ncells = ncells_in_bound + 2;
	else {
		int n_per_cell = 1 + ncells_in_bound / ncells;
		ncells = 2 + ncells_in_bound / n_per_cell;
	}

	if (ncells_in_bound > ncells && ((float)(ncells_in_bound - ncells) / (float)ncells) < 1) ncells = ncells_in_bound + 2;

	// set cells boundaries
	from.resize(ncells);
	to.resize(ncells);

	from[0] = MED_DEFAULT_MIN_TRIM;
	to[0] = from_val + (float)0.5*rounder;
	from[ncells - 1] = to_val + (float)0.5*rounder;
	to[ncells - 1] = MED_DEFAULT_MAX_TRIM;

	float dc = roundf((to_val - from_val) / ((float)(ncells - 2)), rounder);
	//	float dc = (to_val - from_val)/((float)(ncells-2));
	for (int i = 1; i<ncells - 1; i++) {
		from[i] = to[i - 1];
		//		to[i] = from[1] + roundf((float)i*dc,rounder);
		to[i] = from[1] + (float)i*dc;
	}
	to[ncells - 2] = from[ncells - 1];

	// do actual hist
	hist.resize(ncells);
	fill(hist.begin(), hist.end(), 0);
	float epsilon = (float)0.00001; // fighting numerical issues....
	for (int i = 0; i<vals.size(); i++) {
		if (vals[i] != missing_value && (positive_only_flag != 1 || vals[i]>0)) {
			if (vals[i] < (to[0] - epsilon)) hist[0]++;
			else if (vals[i] >= (from[ncells - 1] - epsilon)) hist[ncells - 1]++;
			else {
				int j = (int)(roundf((vals[i] - from[1]), rounder) / dc);
				//				int j = (int)(roundf(((vals[i]-from[1])/dc),rounder));
				if (j<ncells - 2)
					hist[1 + j]++;
				else {
					j = ncells - 3;
					hist[1 + j]++;
				}
			}
		}
	}

	MLOG("DC is %f\n", dc);

	return ncells;
}

//...............................................................................................................................
void MedHist::print(string &prefix)
{
	for (int i = 0; i<ncells; i++) {
		MOUT("%s ", prefix.c_str());
		MOUT("hist %3d : %6.3f - %6.3f : %5d\n", i, from[i], to[i], hist[i]);
	}
}

//...............................................................................................................................
int MedHist::get_cdf(vector<float> &cdf)
{
	int len = (int)hist.size();

	if (len > 0) {
		cdf.resize(len);
		cdf[0] = (float)hist[0];
		for (int i = 1; i<len; i++) {
			cdf[i] = cdf[i - 1] + (float)hist[i];
		}
		if (cdf[len - 1] > 0) {
			for (int i = 0; i<len; i++)
				cdf[i] = cdf[i] / cdf[len - 1];
		}
	}
	return 0;
}

//...............................................................................................................................
int MedHist::sum()
{
	int s = 0;
	for (int i = 0; i<hist.size(); i++)
		s += hist[i];
	return s;
}