#define _CRT_SECURE_NO_WARNINGS

#include <thread>
#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedAlgo/MedAlgo/MedLM.h>
#include <External/Eigen/Core>

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

using namespace Eigen;

//============================================================================
// Linear Model
//============================================================================
//..............................................................................
void init_default_lm_params(MedLMParams& _params) {

	_params.eiter = (float) LM_EITER;
	_params.niter = LM_NITER;

	_params.rfactor = -1 ;
	_params.rfactors = NULL ;
	_params.corrs = NULL ;
	_params.sumxx = NULL ;
}

//..............................................................................
void MedLM::init_defaults()
{
	classifier_type = MODEL_LINEAR_MODEL ;
	transpose_for_learn = true ; 
	transpose_for_predict = false;
	normalize_for_learn = false;
	normalize_y_for_learn = false;
	normalize_for_predict = false;
	
	init_default_lm_params(params) ;

	b0 = 0 ;
	b.clear() ;
}

//..............................................................................
int MedLM::init(void *_in_params) 
{
	init_defaults();

	MedLMParams *in_params = (MedLMParams *) _in_params ;
	 
	params.eiter = in_params->eiter ;
	params.niter = in_params->niter ;

	params.rfactor = in_params->rfactor ;
	params.rfactors = in_params->rfactors ;
	params.corrs = in_params->corrs ;
	params.sumxx = in_params->sumxx ;

	return 0 ;
}

//..............................................................................
int MedLM::set_params(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [MedLM::init]
		if (field == "eiter") params.eiter = stof(entry.second);
		else if (field == "niter") params.niter = stoi(entry.second);
		else if (field == "rfactor") params.rfactor = stof(entry.second);
		else if (field == "get_col") params.get_col = stoi(entry.second);
		else if (field == "rfactors") {
			if (init_farray(entry.second, &(params.rfactors)) == -1) {
				fprintf(stderr, "Cannot initialize rfactors for LM\n");
				return -1;
			}
		} 
		else if (field == "corrs") {
			if (init_farray(entry.second, &(params.corrs)) == -1) {
				fprintf(stderr, "Cannot initialize corrs for LM\n");
				return -1;
			}
		}
		else if (field == "sumxx") {
			if (init_farray(entry.second, &(params.sumxx)) == -1) {
				fprintf(stderr, "Cannot initialize sumxx for LM\n");
				return -1;
			}
		}
		else MLOG("Unknonw parameter \'%s\' for LM\n", field.c_str());
		//! [MedLM::init]
	}
	return 0;
}

//..............................................................................
int init_farray(string& in, float **out) {

	vector<string> vals;
	split(vals, in, boost::is_any_of(","));

	*out = (float *)malloc(sizeof(float)*vals.size());
	if ((*out) == NULL)
		return -1;

	for (unsigned int i = 0; i < vals.size(); i++)
		(*out)[i] = stof(vals[i]);

	return 0;
}

int init_darray(string& in, int **out) {

	vector<string> vals;
	split(vals, in, boost::is_any_of(","));

	*out = (int *)malloc(sizeof(float)*vals.size());
	if ((*out) == NULL)
		return -1;

	for (unsigned int i = 0; i < vals.size(); i++)
		(*out)[i] = stoi(vals[i]);

	return 0;
}


//..............................................................................
MedLM::MedLM () 
{
	init_defaults();
}


//..............................................................................
MedLM::MedLM(MedLMParams& _in_params) 
{
	init((void *) &_in_params);
}

//..............................................................................
MedLM::MedLM(void *_in_params) 
{
	init(_in_params);
}

//..............................................................................
void MedLM::calc_feature_contribs(MedMat<float> &x, MedMat<float> &contribs) {
	contribs.resize(x.nrows, x.ncols + 1); // last is bo
	contribs.signals = x.signals;
	contribs.signals.push_back("b0");
	for (size_t i = 0; i < x.nrows; ++i) {
		for (size_t j = 0; j < x.ncols; ++j)
			contribs(i, j) = x(i, j)*b[j];
		contribs(i, x.ncols) = b0;
	}
	
}

//..............................................................................
int MedLM::Learn (float *x, float *y, int nsamples, int nftrs) {

	vector<float> weights(nsamples,1.0) ;
	return Learn(x,y,&(weights[0]),nsamples,nftrs) ;
	
}
		
//..............................................................................
int MedLM::Learn (float *x, float *y, const float *w, int nsamples, int nftrs) {

	if (params.get_col >= 0) {
		b0 = 0;
		b.resize(nftrs, 0);
		b[params.get_col] = 1;

		return 0;
	}


	if (w == NULL) 
		return (Learn(x,y,nsamples,nftrs));

	// Normalization
	vector<float> x_avg(nftrs), x_std(nftrs);
	float y_avg, y_std;
	normalize_x_and_y(x, y, w, nsamples, nftrs, x_avg, x_std, y_avg, y_std);

	vector<float> _rfactors(nftrs) ;
	vector<float> _corrs(nftrs) ;

	float _rfactor = (params.rfactor >= 0) ? params.rfactor : (float) 1.0 ;

	if (params.rfactors == NULL) {
		for (int i=0; i<nftrs; i++)
			_rfactors[i] = _rfactor ;
	} else {
		_rfactors.assign(params.rfactors,params.rfactors+nftrs) ;
	}

	if (params.corrs == NULL) {
		for (int i=0; i<nftrs; i++)
			_corrs[i] = 0 ;
	} else {
		_corrs.assign(params.corrs,params.corrs+nftrs) ;
	}

	// Learn
	n_ftrs = nftrs ;
	b.resize(nftrs) ;

	int rc;
	if (params.sumxx == NULL)
		rc = learn_lm(x,y,w,nsamples,nftrs,params.niter,params.eiter,&(_rfactors[0]),&(b[0]),&err,&(_corrs[0])) ;
	else
		rc = learn_lm(x,y,w,nsamples,nftrs,params.niter,params.eiter,&(_rfactors[0]),&(b[0]),&err,&(_corrs[0]),params.sumxx) ;

	if (rc == 0)
		denormalize_model(&(x_avg[0]), &(x_std[0]), y_avg, y_std);

	return rc;
}

//..............................................................................
int MedLM::Predict(float *x, float *&preds, int nsamples, int nftrs) const {

	return Predict(x,preds,nsamples,nftrs,0);
}

//..............................................................................
int MedLM::Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) const {

	if (preds == NULL)
		preds = new float[nsamples];
/*
	if (params.get_col >= 0) {

		//MLOG("nsamples %d get_col %d nftrs %d transposed %d\n", nsamples, params.get_col, nftrs, transposed_flag);
		if (transposed_flag) {
			for (int i=0; i<nsamples; i++)
				preds[i] = x[XIDX(params.get_col, i, nsamples)];
		}
		else {
			for (int i=0; i<nsamples; i++)
				preds[i] = x[XIDX(i, params.get_col, nftrs)];
		}
		return 0;
	}
*/
	if (preds == NULL)
		preds = new float[nsamples];

	memset(preds,0,nsamples*sizeof(float)) ;

	int ncores = std::thread::hardware_concurrency();
	Eigen::setNbThreads(3*ncores/4);

	Map<const MatrixXf> bf(&b[0],nftrs,1);
	Map<MatrixXf> pf(preds,nsamples,1);
	if (transposed_flag) {
		Map<MatrixXf> xf(x,nsamples,nftrs);
		//for (int j=0; j<nftrs; j++) {
		//	for (int i=0; i<nsamples; i++)
		//		preds[i] += b[j] * x[XIDX(j,i,nsamples)] ;
		//}
		pf = xf*bf;
	} else {
		Map<MatrixXf> xf(x,nftrs,nsamples);
		//for (int i=0; i<nsamples; i++) {
		//	for (int j=0; j<nftrs; j++)
		//		preds[i] += b[j] * x[XIDX(i,j,nftrs)] ;
		//}
		pf = xf.transpose() * bf;
	}

	pf.array() += b0;


	//for (int i=0; i<nsamples; i++)
	//	preds[i] += b0 ;

	return 0;
}

//..............................................................................
void MedLM::normalize_x_and_y(float *x, float *y, const float *w, int nsamples, int nftrs, vector<float>& x_avg, vector<float>& x_std, float& y_avg, float& y_std) {

	// Get moments
	int n_clean;
	for (int i = 0; i < nftrs; i++)
		medial::stats::get_mean_and_std(x + nsamples * i, w, nsamples, (float) -1.0, x_avg[i], x_std[i], n_clean, false);
	medial::stats::get_mean_and_std(y, w, nsamples, (float)-1.0, y_avg, y_std, n_clean, false);

	// Normalize
	for (int i = 0; i < nftrs; i++) {
		float *xi = x + i*nsamples;
		for (int j = 0; j < nsamples; j++)
			xi[j] = (xi[j] - x_avg[i]) / x_std[i];
	}

	for (int j = 0; j < nsamples; j++)
		y[j] = (y[j] - y_avg) / y_std;

}

//..............................................................................
int MedLM::denormalize_model(float *f_avg, float *f_std, float label_avg, float label_std)
{
	float new_b0;
	vector<float> new_b(n_ftrs);

	new_b0 = b0*label_std + label_avg;
	fill(new_b.begin(),new_b.end(),(float)0);
	for (int j=0; j<n_ftrs; j++) {
		new_b[j] = label_std*b[j]/f_std[j];
		new_b0 -= label_std*f_avg[j]*b[j]/f_std[j];
	}

	b0 = new_b0;
	for (int j=0; j<n_ftrs; j++)
		b[j] = new_b[j];

	transpose_for_predict = false;
	normalize_for_predict = false;
	return 0;
}


//..............................................................................
int learn_lm (float *x, float *_y, const float *w, int nsamples, int nftrs, int niter, float eiter , float *rfactors, float *b, float *err, float *corrs) {

	// Prepare
	float *sumxx ;
	if ((sumxx = (float *) malloc (nftrs*sizeof(float)))==NULL) {
		MERR("Sumxx allocation failed\n") ;
		return -1 ;
	}
	
	for (int j=0; j<nftrs; j++) {
		sumxx[j] = 0 ;
		for (int i=0; i<nsamples; i++)
			sumxx[j] += w[i]*x[XIDX(j,i,nsamples)]*x[XIDX(j,i,nsamples)];
	}
	int rc = learn_lm(x,_y,w,nsamples,nftrs,niter,eiter,rfactors,b,err,corrs,sumxx) ;

	free(sumxx) ;

	return rc ;

}

//..............................................................................
int learn_lm (float *x, float *_y, const float *w, int nsamples, int nftrs, int niter, float eiter , float *rfactors, float *b, float *err, float *corrs, float *sumxx)  {

	// Prepare
	float *y ;
	if ((y = (float *) malloc (nsamples*sizeof (float))) == NULL)
		return -1 ;

	memcpy(y,_y,nsamples*sizeof(float)) ;	
	memset(b,0,nftrs*sizeof(float)) ;

	// Do the iterations
	int it = 0 ;
	float  prev_err,alpha,oldb ;
	float sumxy ;
	while (it < niter && (it <= 1 || abs((*err)-prev_err)/(*err)>eiter)) {
	
		it++ ;
		prev_err = *err ;

		for(int j=0; j<nftrs; j++) {
			sumxy = 0 ;

			for (int i = 0; i < nsamples; i++) {
				sumxy += w[i] * x[XIDX(j, i, nsamples)] * y[i];
				//fprintf(stderr, "%f %f %f \n", w[i], x[XIDX(j, i, nsamples)], y[i]);
			}

			alpha = (sumxx[j]>0) ? (sumxy/sumxx[j]) : 0 ;
			if (alpha*corrs[j] < 0)
				alpha = 0 ;

			oldb=b[j];

			b[j] = rfactors[j] * (b[j] + alpha) ;
			alpha = b[j] - oldb;

			for(int i=0; i<nsamples; i++)
				y[i] -= alpha * x[XIDX(j,i,nsamples)] ;
		}

		*err = 0 ;
		for(int i=0; i<nsamples; i++)
			(*err)+=w[i]*y[i]*y[i];
	}	

	free(y) ;
	return 0 ;
}




//..............................................................................
void MedLM::print(FILE *fp, const string& prefix, int level) const {

	if (level == 0) 
		fprintf(fp, "%s: MedLM ()\n", prefix.c_str());
	else {
		fprintf(fp, "%s : Linear Model : Nftrs = %d\n", prefix.c_str(), n_ftrs);
		fprintf(fp, "%s : Linear Model b0 = %f\n", prefix.c_str(), b0);

		for (int i = 0; i < n_ftrs; i++)
			fprintf(fp, "%s : Linear Model b[%d] = %f\n", prefix.c_str(), i, b[i]);
	}
}

//..............................................................................
