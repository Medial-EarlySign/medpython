//
// MedStat
//

#include "MedStat.h"

#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_MEDSTAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//...............................................................................................................................
// get a vector of values, a vector of probabilities, and returning a matching vector of values such that Prob(x<=pvals[i])=p[i]
void get_percentiles(vector<float> &vals, vector<float> &p, vector<float> &pvals)
{
	//sort(vals.begin(),vals.end());
	//
	//int n = (int) vals.size();

	//pvals.resize(p.size());
	//for (int i=0; i<p.size(); i++) {
	//	int k = (int)((float)n*p[i]);
	//	if (k < 0) k = 0;
	//	if (k >= n) k = n-1;
	//	pvals[i] = vals[k];
	//}
	get_percentiles(vals,p,pvals,0);
}

//...............................................................................................................................
// get a vector of values, a vector of probabilities, and returning a matching vector of values such that Prob(x<=pvals[i])=p[i]
void get_percentiles(vector<float> &vals, vector<float> &p, vector<float> &pvals, int only_positive_flag)
{
	sort(vals.begin(),vals.end());
	
	int n = (int) vals.size();
	int n0 = 0;
	if (only_positive_flag) {
		while (n0<n && vals[n0]<=0) n0++;
	}
	pvals.resize(p.size());
	pvals.clear();
	for (int i=0; i<p.size(); i++) {
		int k = n0 + (int)((float)(n-n0)*p[i]);
		if (k < 0) k = 0;
		if (k >= n) k = n-1;
		pvals[i] = vals[k];
	}
}
//...............................................................................................................................
// gets a vector of values, and checks the best way to round it at the 10^-3 to 10^3 (at 10x steps) range.
float get_best_rounding(vector<float> &vals, vector<float>& res, vector<int> &counts, float missing_value)
{
	int n_r = 7;
	res = {(float)0.001,(float)0.01,(float)0.1,1,10,100,1000};
	vector<int> cnts(n_r,0);

	// find best rounder for each value and count them
	for (int i=0; i<vals.size(); i++) {
		for (int j=n_r-1; j>=0; j--) {
			if (vals[i] != 0 && vals[i] != missing_value) {
				float vr = ((float)((int)((vals[i]/res[j]))));
				if (abs((float)(vals[i]/res[j])-vr) == 0) {
					cnts[j]++;
					break;
				}
			}
		}

	}


	// return most common rounder
	int best = 0, best_cnt = 0;
	for (int j=0; j<n_r; j++) {		
		if (cnts[j] > best_cnt) {
			best = j;
			best_cnt = cnts[j];
		}
	}
	counts = cnts;
	return res[best];
}

//...............................................................................................................................
void get_mean(vector<float> &vals, float &mean)
{
	double sx = 0;

	for (int i = 0; i<vals.size(); i++)
		sx += vals[i];

	if (vals.size() > 0)
		mean = (float)sx / (float)vals.size();
	else
		mean = 0;
}

//...............................................................................................................................
void get_common(vector<float>& vals, float &common) {

	map<float, int> counter;
	for (float val : vals)
		counter[val]++;

	int count = 0;
	for (auto& rec : counter) {
		if (rec.second > count) {
			count = rec.second;
			common = rec.first;
		}
	}
}

//...............................................................................................................................
void get_mean_and_std(vector<float> &vals, float &mean, float &std)
{
	if (vals.size() == 0) {
		mean = 0;
		std = 1;
	}
	else {
		double sum = 0;
		int n = (int)vals.size();
		for (int i = 0; i < vals.size(); i++)
			sum += vals[i];
		mean = (float) (sum / n);

		if (n == 1)
			std = 1.0;
		else {
			sum = 0;
			for (int i = 0; i < vals.size(); i++)
				sum += (vals[i] - mean)*(vals[i] - mean);
			std = sqrt((float)(sum / n));
		}
	}

}

//...............................................................................................................................
void get_mean_and_std(vector<float> &vals, float &mean, float &std, float missing_val)
{
	double sum = 0;
	int n = 0;
	for (int i = 0; i < vals.size(); i++) {
		if (vals[i] != missing_val) {
			n++;
			sum += vals[i];
		}
	}

	if (n == 0) {
		mean = 0;
		std = 1.0;
	}
	else {
		mean = (float) (sum / n);

		if (n == 1)
			std = 1.0;
		else {
			sum = 0;
			for (int i = 0; i < vals.size(); i++) {
				if (vals[i] != missing_val)
					sum += (vals[i] - mean)*(vals[i] - mean);
			}
			std = sqrt((float)(sum / n));
		}
	}
}


//...............................................................................................................................
// MedHist
//...............................................................................................................................
int MedHist::get_hist(vector<float> &vals)
{
	if (rounder < 0) rounder = (float)0.01;

	// get boundaries
	if (to_val < from_val) {
		vector<float> p(2),pvals(2);
		p[0] = min_percentile;
		p[1] = max_percentile;
		get_percentiles(vals,p,pvals,positive_only_flag);
		from_val = roundf(pvals[0],rounder);
		to_val = roundf(pvals[1],rounder);
	}

	// reevaluate ncells
	int ncells_in_bound = (int)((to_val-from_val)/rounder);
	if (ncells >= ncells_in_bound) 
		ncells = ncells_in_bound + 2;
	else {
		int n_per_cell = 1 + ncells_in_bound/ncells;
		ncells = 2+ncells_in_bound/n_per_cell;
	}

	if (ncells_in_bound > ncells && ((float)(ncells_in_bound-ncells)/(float)ncells) < 1) ncells = ncells_in_bound+2;

	// set cells boundaries
	from.resize(ncells);
	to.resize(ncells);

	from[0] = MED_DEFAULT_MIN_TRIM;
	to[0] = from_val + (float)0.5*rounder;
	from[ncells-1] = to_val + (float)0.5*rounder;
	to[ncells-1] = MED_DEFAULT_MAX_TRIM;
	
	float dc = roundf((to_val - from_val)/((float)(ncells-2)),rounder);
//	float dc = (to_val - from_val)/((float)(ncells-2));
	for (int i=1; i<ncells-1; i++) {
		from[i] = to[i-1];
//		to[i] = from[1] + roundf((float)i*dc,rounder);
		to[i] = from[1] + (float)i*dc;
	}
	to[ncells-2] = from[ncells-1];

	// do actual hist
	hist.resize(ncells);
	fill(hist.begin(), hist.end(), 0);
	float epsilon = (float)0.00001; // fighting numerical issues....
	for (int i=0; i<vals.size(); i++) {
		if (vals[i] != missing_value && (positive_only_flag != 1 || vals[i]>0)) {
			if (vals[i] < (to[0]-epsilon)) hist[0]++;
			else if (vals[i] >= (from[ncells-1]-epsilon)) hist[ncells-1]++;
			else {
				int j = (int)(roundf((vals[i]-from[1]),rounder)/dc);
//				int j = (int)(roundf(((vals[i]-from[1])/dc),rounder));
				if (j<ncells-2)
					hist[1+j]++;
				else {
					j = ncells-3;
					hist[1+j]++;
				}
			}
		}
	}

	MLOG("DC is %f\n",dc);

	return ncells;
}

//...............................................................................................................................
void MedHist::print(string &prefix)
{
	for (int i=0; i<ncells; i++) {
		MOUT("%s ",prefix.c_str());
		MOUT("hist %3d : %6.3f - %6.3f : %5d\n",i,from[i],to[i],hist[i]);
	}
}

//...............................................................................................................................
int MedHist::get_cdf(vector<float> &cdf)
{
	int len = (int)hist.size();

	if (len > 0) {
		cdf.resize(len);
		cdf[0] = (float)hist[0];
		for (int i=1; i<len; i++) {
			cdf[i] = cdf[i-1] + (float)hist[i];
		}
		if (cdf[len-1] > 0) {
			for (int i=0; i<len; i++)
				cdf[i] = cdf[i]/cdf[len-1];
		}
	}
	return 0;
}

//...............................................................................................................................
int MedHist::sum()
{
	int s = 0;
	for (int i=0; i<hist.size(); i++)
		s += hist[i];
	return s;
}