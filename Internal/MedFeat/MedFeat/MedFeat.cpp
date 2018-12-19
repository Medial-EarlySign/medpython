//
// MedFeat - Tools to help build features, and for feature selection
//

#if 0

#include "MedFeat.h"
#include "boost/foreach.hpp"
#include <vector>
#include <MedIO/MedIO/MedIO.h>

#define LOCAL_SECTION LOG_MEDFEAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;
//==========================================================================
// loadFromRepository
//==========================================================================
int MedFeaturesData::loadFromRepository(MedRepository &rep, string readlabel_name, vector<string> sigs_name)
{


	clear();
	label_name = readlabel_name;

	
	
	// put the name sin
	for (auto name:sigs_name)
		signals.push_back(name);

	MLOG("going over %d patient records\n",rep.index.pids.size());
	//put the labels in
	int n_in = 0;
	for (auto pid:rep.index.pids) {
		
		
		int len;
		SDateVal *sdv = (SDateVal *)rep.get(pid,label_name,len);

		for(int r =0;r<len;r++)
			label.push_back(sdv[r].val);
		

//		
		for (auto name:signals) {
		
			sdv = (SDateVal *)rep.get(pid,name,len);
			for(int r =0;r<len;r++)
				data[name].push_back(sdv[r].val);
		}
	}		

	return 0;
}
//==========================================================================
// MedFeaturesData
//==========================================================================
int MedFeaturesData::get_learning_nrows(int isplit) {
	int nrows = 0;
	for (unsigned int i = 0; i<splits.size(); i++) {
		if (splits[i] != isplit)
			nrows++;
	}
	return nrows;
}

//..........................................................................
int MedFeaturesData::get_testing_nrows(int isplit) {

	int nrows = 0 ;
	for (unsigned int i=0; i<splits.size(); i++) {
		if (splits[i] == isplit )
			nrows ++ ;
	}

	return nrows ;
}


void MedFeaturesData::setup_cleaners_for_all_folds(const string& signal, bool trim_flag, bool remove_flag, bool replace_missing_to_mean_flag, bool normalize_flag, float missing_v) {

	cleaners[signal].resize(nsplits + 1);
	for (int isplit = 0; isplit < nsplits; isplit++) {
		vector<float> learn_data;

		//MLOG("Learning cleaner for sig %s , split %d\n", signal.c_str(), isplit);
		for (unsigned int i = 0; i<data[signal].size(); i++) {
			if ((splits[i] != isplit) || nsplits == 1) {
				learn_data.push_back(data[signal][i]);
			}
		}
		cleaners[signal][isplit].missing_value = missing_v;
		cleaners[signal][isplit].trim_flag = trim_flag;
		cleaners[signal][isplit].remove_flag = remove_flag;
		cleaners[signal][isplit].replace_missing_to_mean_flag = replace_missing_to_mean_flag;
		cleaners[signal][isplit].normalize_flag = normalize_flag;
		cleaners[signal][isplit].get_limits_iteratively(learn_data); // this one sets also mean and sdv, so no need to call it after.

		//MLOG("cleaners: signal %s , isplit %d\n", signal.c_str(), isplit);
		//cleaners[signal][isplit].print(signal);
	}
	if (nsplits > 1) {
		//MLOG("Learning cleaner for sig %s , split %d\n", signal.c_str(), nsplits);
		vector<float> learn_data = data[signal];
		cleaners[signal][nsplits].missing_value = missing_v;
		cleaners[signal][nsplits].trim_flag = trim_flag;
		cleaners[signal][nsplits].remove_flag = remove_flag;
		cleaners[signal][nsplits].replace_missing_to_mean_flag = replace_missing_to_mean_flag;
		cleaners[signal][nsplits].normalize_flag = normalize_flag;
		cleaners[signal][nsplits].get_limits_iteratively(learn_data);
	}
	else
		cleaners[signal][1] = cleaners[signal][0];
	//MLOG("cleaners: signal %s , isplit %d\n", signal.c_str(), nsplits);
	//cleaners[signal][nsplits].print(signal);
}


void MedFeaturesData::get_normalization(const string& signal) {
	setup_cleaners_for_all_folds(signal, false, false, true, true);
}

void MedFeaturesData::get_label_normalization() {

	label_cleaners.resize(nsplits+1) ;
	for (int isplit=0; isplit < nsplits; isplit ++) {
		vector<float> learn_data ;

		for (unsigned int i=0; i<label.size(); i++) {
			if (splits[i] != isplit)
				learn_data.push_back(label[i]) ;
		}

		label_cleaners[isplit].trim_flag = false ;
		label_cleaners[isplit].remove_flag = false ;
		label_cleaners[isplit].replace_missing_to_mean_flag = true;
		label_cleaners[isplit].normalize_flag = true;
		label_cleaners[isplit].get_mean_and_sdv(learn_data) ;
	}

	label_cleaners[nsplits].trim_flag = false ;
	label_cleaners[nsplits].remove_flag = false ;
	label_cleaners[nsplits].replace_missing_to_mean_flag = true;
	label_cleaners[nsplits].normalize_flag = true;
	label_cleaners[nsplits].get_mean_and_sdv(label) ;
}

void MedFeaturesData::get_normalization_and_cleaning(const string& signal, float missing_v) {
	setup_cleaners_for_all_folds(signal, true, true, true, true, missing_v);
}

void MedFeaturesData::get_all_normalization_and_cleaning(float missing_v)
{
	// signals
	for (auto &it : data)  {
		//MLOG("getting normalization and cleaning for feature %s\n", it.first.c_str());
		get_normalization_and_cleaning(it.first,missing_v);
	}
	// label
	//MLOG("getting label normalization");
	get_label_normalization();
}

void MedFeaturesData::get_features_as_mat(vector<string> &f_names, MedMat<float> &x, int _isplit, int flags)
{
	int nsamples,nftrs;
	int n_all = 0;
	x.clear();

	if (_isplit == nsplits)
		nsamples = (int)splits.size();
	else {
		if (flags & Split_Equal)
			nsamples = get_testing_nrows(_isplit);
		else
			nsamples = get_learning_nrows(_isplit);
	}

	nftrs = (int)f_names.size();

	if (flags & Transpose) {
		x.resize(nftrs,nsamples);
		x.transposed_flag = 1;
	} else
		x.resize(nsamples,nftrs);

	if (nftrs > 0)
		n_all = (int)data[f_names[0]].size();

	int k;
	MLOG_D("get_mat:: nsamples %d nftrs %d ::  x: %d x %d :: n_all %d :: split %d nsplits %d\n",nsamples,nftrs,x.nrows,x.ncols,n_all,_isplit,nsplits); 
	for (int j=0; j<nftrs; j++) {
		k = 0;
		MedCleaner *cln = &cleaners[f_names[j]][_isplit];
		float *d = VEC_DATA(data[f_names[j]]);
		float val;
		for (int i=0; i<n_all; i++) {
			if ((_isplit == nsplits) || ((flags & Split_Equal) && (splits[i] == _isplit)) || (!(flags & Split_Equal) && (splits[i] != _isplit))) {

//				float val = data[f_names[j]][i];
				val = d[i];

				if (flags & Clean)
//					cleaners[f_names[j]][_isplit].single_remove_trim_replace(val);
					cln->single_remove_trim_replace(val);

				if (flags & Normalize)
//					cleaners[f_names[j]][_isplit].single_normalize(val);
					cln->single_normalize(val);

				if (flags & Transpose)
					x(j,k) = val;
				else
					x(k,j) = val;

				k++;

			}
		}
	}

	if (flags & Normalize) {
		x.normalized_flag = 1;
		x.avg.clear();
		x.std.clear();
		for (int j=0; j<nftrs; j++) {
			x.avg.push_back(cleaners[f_names[j]][_isplit].mean);
			x.std.push_back(cleaners[f_names[j]][_isplit].sdv);
		}
	}


}

void MedFeaturesData::get_label_as_mat(MedMat<float> &y, int _isplit, int flags)
{
	int nsamples;
	int n_all = 0;
	y.clear();

	if (_isplit == nsplits)
		nsamples = (int)splits.size();
	else {
		if (flags & Split_Equal)
			nsamples = get_testing_nrows(_isplit);
		else
			nsamples = get_learning_nrows(_isplit);
	}
	y.resize(nsamples,1);

	n_all = (int)label.size();

	int k = 0;
	for (int i=0; i<n_all; i++) {
		if ((_isplit == nsplits) || ((flags & Split_Equal) && (splits[i] == _isplit)) || (!(flags & Split_Equal) && (splits[i] != _isplit))) {

			float val = label[i];

			if (flags & Clean)
				label_cleaners[_isplit].single_remove_trim_replace(val);

			if (flags & Normalize)
				label_cleaners[_isplit].single_normalize(val);

			y(k,0) = val;
			k++;

		}
	}
	y.resize(k,1);

	if (flags & Normalize) {
		y.normalized_flag = 1;
		y.avg.clear();
		y.std.clear();
		y.avg.push_back(label_cleaners[_isplit].mean);
		y.std.push_back(label_cleaners[_isplit].sdv);
	}


}

// note that when write/read cleaners one has to make sure to call the matching read with exactly the same nsplits and exactly the same order of signals in the vector sigs.
int MedFeaturesData::write_cleaners_to_file(const string &fname, vector<string> &sigs)
{
	unsigned long long tot_size = cleaners[sigs[0]][0].get_size()*(nsplits+1)*sigs.size();

	vector<unsigned char> serialized_cleaners(tot_size);

	unsigned char *buf = &serialized_cleaners[0];

	size_t size;
	unsigned char *curr_buf = buf;
	for (int i=0; i<sigs.size(); i++)
		for (int j=0; j<=nsplits; j++) {
			size = cleaners[sigs[i]][j].serialize(curr_buf);
			curr_buf += size;
		}
	
	if (write_binary_data(fname,buf,tot_size) < 0) {
		MERR("Error: Failed writing file %s\n",fname.c_str());
		return -1;
	}

	return 0;
}

int MedFeaturesData::read_cleaners_from_file(const string &fname, vector<string> &sigs)
{
	unsigned char *buf;
	unsigned long long size;

	if (read_binary_data_alloc(fname,buf,size) < 0) {
		MERR("Error: failed reading from file %s\n",fname.c_str());
		return -1;
	}

	vector<MedCleaner> empty;
	size_t s;
	unsigned char *curr_buf = buf;
	for (int i=0; i<sigs.size(); i++) {
		if (cleaners.find(sigs[i]) == cleaners.end())
			cleaners[sigs[i]] = empty;
		else
			cleaners[sigs[i]].clear();

		for (int j=0; j<=nsplits; j++) {
			MedCleaner c;
			s = c.deserialize(curr_buf);
			cleaners[sigs[i]].push_back(c);
			curr_buf += s;
		}
	}

	delete [] buf;
	return 0;
}

int MedFeaturesData::write_label_cleaners_to_file(const string &fname)
{
	unsigned long long tot_size = label_cleaners[0].get_size()*(nsplits+1);

	vector<unsigned char> serialized_cleaners(tot_size);

	unsigned char *buf = &serialized_cleaners[0];

	size_t size;
	unsigned char *curr_buf = buf;
	for (int j=0; j<=nsplits; j++) {
		size = label_cleaners[j].serialize(curr_buf);
		curr_buf += size;
	}
	
	if (write_binary_data(fname,buf,tot_size) < 0) {
		MERR("Error: Failed writing file %s\n",fname.c_str());
		return -1;
	}

	return 0;
}

int MedFeaturesData::read_label_cleaners_to_file(const string &fname)
{
	unsigned char *buf;
	unsigned long long size;

	if (read_binary_data_alloc(fname,buf,size) < 0) {
		MERR("Error: failed reading from file %s\n",fname.c_str());
		return -1;
	}

	label_cleaners.clear();
	size_t s;
	unsigned char *curr_buf = buf;
	for (int j=0; j<=nsplits; j++) {
		MedCleaner c;
		s = c.deserialize(curr_buf);
		label_cleaners.push_back(c);
		curr_buf += s;
	}

	delete [] buf;
	return 0;
}

//............................................................................................................................
int MedFeaturesData::apply_clean_and_normalize(int i_split, float missing_v, string &sig, int normalize_only_flag) // for a single split - if == nsplits then on all
{
	for (int i=0; i<data[sig].size(); i++) {
		
		if (i_split == nsplits || splits[i] == i_split) {

			float val = data[sig][i];

			if (normalize_only_flag == 0) {
				if (val != missing_v) {
					// remove
					if (val < cleaners[sig][i_split].remove_min)
						val = missing_v;
					if (val > cleaners[sig][i_split].remove_max)
						val = missing_v;
				}

				if (val != missing_v) {
					// trim
					if (val < cleaners[sig][i_split].trim_min)
						val = cleaners[sig][i_split].trim_min;
					if (val > cleaners[sig][i_split].trim_max)
						val = cleaners[sig][i_split].trim_max;

				}
			}

			if (val != missing_v) {
				// normalize
				val = (val - cleaners[sig][i_split].mean)/cleaners[sig][i_split].sdv;
			} else
				val = 0;

			data[sig][i] = val;
		}

	}

	return 0;
}


//............................................................................................................................
int MedFeaturesData::apply_clean_and_normalize(float missing_v, string &sig, int i_cleaner, int normalize_only_flag) // for all splits
{
	int i_split;

	if (data.find(sig) == data.end()) {MERR("MedFeaturesData: apply_clean_and_normalize() ERROR !!!!!!! sig %s not in mf.\n",sig.c_str());}

	float *vdata;
	vector<MedCleaner> *vcln;
	MedCleaner *cln;

	vdata = VEC_DATA(data[sig]);
	vcln = &cleaners[sig];
	int nsamples = (int)data[sig].size();
	int printed = 0;
	register float val;
	for (int i=0; i<nsamples; i++) {
		
		if (i_cleaner < 0) {
			i_split = splits[i];
			cln = &((*vcln)[i_split]);
		} else {
			cln = &((*vcln)[i_cleaner]);
		}	
		val = vdata[i];

		if (normalize_only_flag == 0) {
			if (val != missing_v) {
				// remove
				if (val < cln->remove_min)
					val = missing_v;
				else if (val > cln->remove_max)
					val = missing_v;
			}

			if (val != missing_v) {
				// trim
				if (val < cln->trim_min)
					val = cln->trim_min;
				else if (val > cln->trim_max)
					val = cln->trim_max;

			}
		}

		if (val != missing_v) {
			// normalize
			val = (val - (cln->mean))/(cln->sdv);
		} else
			val = 0;


		vdata[i] = val;

	}

	return 0;

}



//==========================================================================
// SDateVal Features
//==========================================================================

//............................................................................................................................
// returns -1 if len is 0
float sdv_get_min(SDateVal *sdv, int len)
{
	if (sdv == NULL || len == 0)
		return -1;

	float minval = sdv[0].val;
	for (int i=1; i<len; i++)
		if (sdv[i].val < minval)
			minval = sdv[i].val;

	return minval;
}

//............................................................................................................................
// returns -1 if len is 0
float sdv_get_max(SDateVal *sdv, int len)
{
	if (sdv == NULL || len == 0)
		return -1;

	float maxval = sdv[0].val;
	for (int i=1; i<len; i++)
		if (sdv[i].val > maxval)
			maxval = sdv[i].val;

	return maxval;
}

//............................................................................................................................
// returns -1 if len is 0
float sdv_get_avg(SDateVal *sdv, int len)
{
	if (sdv == NULL || len == 0)
		return -1;

	float sum = 0;
	for (int i=0; i<len; i++)
		sum += sdv[i].val;
	return(sum/(float)len);
}

//............................................................................................................................
float sdv_get_time_adjusted_avg(SDateVal *sdv, int len)
{
	if (sdv == NULL || len == 0)
		return -1;

	if (len == 1)
		return sdv[0].val;

	float sum = 0;
	float tot_time_diff = (float)(date_to_days(sdv[len-1].date) - date_to_days(sdv[0].date));
	if (tot_time_diff == 0) tot_time_diff++;
	for (int i=1; i<len; i++) {
		float dt = (float)(date_to_days(sdv[i].date) - date_to_days(sdv[i-1].date));
		dt = dt/tot_time_diff;
		sum += dt*(float)0.5*(sdv[i].val+sdv[i-1].val);
	}

	return sum;

}

//............................................................................................................................
// returns -1 if len is 0
float sdv_get_std(SDateVal *sdv, int len)
{
	if (sdv == NULL || len == 0)
		return -1;

	float sx = 0;
	float sxx = 0;
	for (int i=0; i<len; i++) {
		sx += sdv[i].val;
		sxx += sdv[i].val*sdv[i].val;
	}

	sxx = sxx/(float)len;
	sx = sx/(float)len;

	return(sqrt(abs(sxx-sx*sx)));
}

//............................................................................................................................
// slope is 0 (= no change) if there are less than 2 samples
float sdv_get_slope(SDateVal *sdv, int len)
{
	double sx,sy,sxx,sxy,n;

	if (len < 2)
		return 0;

	sx = 0;
	sy = 0;
	sxx = 0;
	sxy = 0;
	for (int i=0; i<len; i++) {
		double x = (double)(date_to_days(sdv[i].date) - date_to_days(sdv[0].date))/(double)365;
		sx += x;
		sy += sdv[i].val;
		sxx += x*x;
		sxy += x*sdv[i].val;
	}
	n = (double)len;

	double cov = sxy - sx*sy/n;
	double var = sxx - sx*sx/n;

	if (var < 0.1)
		return 0;

	return ((float)(cov/var));
}

//............................................................................................................................
float sdv_get_linear_delta(SDateVal *sdv, int len, int date)
{
	if (len < 1)
		return 0;

	float val = sdv_get_linear_val(sdv,len,date);
	return (val - sdv[len-1].val);
}

//............................................................................................................................
float sdv_get_linear_val(SDateVal *sdv, int len, int date)
{
	if (len < 1)
		return 0;

	float beta = sdv_get_slope(sdv,len);

	float dt = (float)(date_to_days(date) - date_to_days(sdv[len-1].date))/(float)365;

	float val = sdv[len-1].val + dt*beta;

	return val;
}

//............................................................................................................................
float sdv_get_fraction_below(SDateVal *sdv, int len, float bound)
{
	if (sdv == NULL || len == 0)
		return 0;

	int nbelow = 0;
	for (int i=0; i<len; i++)
		if (sdv[i].val <= bound)
			nbelow++;
	return((float)nbelow/(float)len);
}

//............................................................................................................................
float sdv_get_fraction_above(SDateVal *sdv, int len, float bound)
{
	if (sdv == NULL || len == 0)
		return 0;

	int nabove = 0;
	for (int i=0; i<len; i++)
		if (sdv[i].val >= bound)
			nabove++;
	return((float)nabove/(float)len);
}

#endif
