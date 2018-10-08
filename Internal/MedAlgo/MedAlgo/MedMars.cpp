//
// MedMars - wrapper for the Mars Library
//

#define _CRT_SECURE_NO_WARNINGS

#include "MedAlgo.h"
#include "Mars/Mars/earth.hpp"


#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//..............................................................................
void init_default_mars_params(MedMarsParams& _params) {

	_params.MaxTerms = MED_MARS_DEF_MAXTERMS;
	_params.MaxDegree = MED_MARS_DEF_DEGREE;
	_params.Penalty = (double)MED_MARS_DEF_PENALTY;
	_params.Thresh = (double)MED_MARS_DEF_THRESH;
	_params.MinSpan = MED_MARS_DEF_MINSPAN;
	_params.Prune = (bool)MED_MARS_DEF_PRUNE;
	_params.FastK = MED_MARS_DEF_FASTK;
	_params.FastBeta = (double)MED_MARS_DEF_FASTBETA;
	_params.NewVarPenalty = (double)MED_MARS_DEF_NEWVAR_PENALTY;
	_params.UseBetaCache = (bool)MED_MARS_DEF_USEBETACACHE;
	_params.Trace = (double)MED_MARS_DEF_TRACE;
}

//..............................................................................
void MedMars::init_defaults()
{
	classifier_type = MODEL_MARS;
	transpose_for_learn = false;
	normalize_for_learn = false;
	normalize_y_for_learn = false;
	transpose_for_predict = false;
	normalize_for_predict = false;

	nMaxTerms = 0;
	nPreds = 0;
//	BestSet.clear();
	BestSet = NULL;
	Dirs.clear();
	Cuts.clear();
	Betas.clear();
	
	// Model Inner quality measures
	BestGcv = 0;
	bx.clear();
	Residuals.clear();

	init_default_mars_params(params) ;

}

//..............................................................................
int MedMars::set_params(map<string, string>& mapper) {

	for (auto entry : mapper) {
		string field = entry.first;
		//! [MedMars::init]
		if (field == "MaxTerms") params.MaxTerms = stoi(entry.second);
		else if (field == "MaxDegree") params.MaxDegree = stoi(entry.second);
		else if (field == "Penalty") params.Penalty = stod(entry.second);
		else if (field == "Thresh") params.Thresh = stod(entry.second);
		else if (field == "MinSpan") params.MinSpan = stoi(entry.second);
		else if (field == "Prune") params.Prune = (bool)(stoi(entry.second) != 0);
		else if (field == "FastK") params.FastK = stoi(entry.second);
		else if (field == "FastBeta") params.FastBeta = stod(entry.second);
		else if (field == "NewVarPenalty") params.NewVarPenalty = stod(entry.second);
		else if (field == "UseBetaCache") params.UseBetaCache = (bool)(stoi(entry.second) != 0);
		else if (field == "Trace") params.Trace = stod(entry.second);
		else MLOG("Unknonw parameter \'%s\' for Mars\n", field.c_str());
		//! [MedMars::init]
	}
	return 0;
}

//..............................................................................
int MedMars::init(void *_in_params) 
{
	init_defaults();

	MedMarsParams *in_params = (MedMarsParams *) _in_params ;
	 
	params = (*in_params);

	return 0 ;
}

//..............................................................................
MedMars::MedMars() 
{
	init_defaults();
}

//..............................................................................
MedMars::MedMars(MedMarsParams& _in_params) 
{
	init((void *) &_in_params);
}

//..............................................................................
MedMars::MedMars(void *_in_params) 
{
	init(_in_params);
}

//..............................................................................
int MedMars::Learn (float *x, float *y, const float *w, int nsamples, int nftrs) {

	if (w != NULL)
		MWARN("Weights are not implemented for Mars. Ignoring\n") ;
	
	// allocate space for model
	nMaxTerms = params.MaxTerms;
	nPreds = nftrs;
	if (BestSet != NULL) free(BestSet);
	BestSet = (bool *)calloc(nMaxTerms,sizeof(bool));
//	BestSet.resize(nMaxTerms); fill(BestSet.begin(),BestSet.end(),false);
	Dirs.resize(nMaxTerms*nPreds); fill(Dirs.begin(), Dirs.end(), 0);
	Cuts.resize(nMaxTerms*nPreds); fill(Cuts.begin(), Cuts.end(), (double)0);
	Betas.resize(nMaxTerms); fill(Betas.begin(), Betas.end(), (double)0);

	// transpose x and make it a double to match inner package req
	MedMat<double> xd;
	xd.load_transposed(x,nsamples,nftrs);

	// transfer y to double
	MedMat<double> yd;
	yd.load(y,nsamples,1);

	// a few more necessary inputs for inner learn
	bx.resize(nsamples*nMaxTerms); fill(bx.begin(),bx.end(),(double)0);
	Residuals.resize(nsamples); fill(Residuals.begin(), Residuals.end(), (double)0);

	vector<int> LinPreds(nftrs,0);
	const char **pred_names = NULL;
	double *weights = NULL; // unsupported yet in internal imp

	
	// Actual Learn !
    Earth(&BestGcv, &nTerms, BestSet, VEC_DATA(bx), VEC_DATA(Dirs), VEC_DATA(Cuts), VEC_DATA(Residuals), VEC_DATA(Betas),
        xd.data_ptr(), yd.data_ptr(), weights, nsamples, 1, nftrs,
        params.MaxDegree, nMaxTerms, params.Penalty, params.Thresh, params.MinSpan, params.Prune,
        params.FastK, params.FastBeta, params.NewVarPenalty, VEC_DATA(LinPreds), params.UseBetaCache, params.Trace, pred_names);
	
	fflush(stdout);
/*
    printf("Expression:\n");
    FormatEarth(BestSet, VEC_DATA(Dirs), VEC_DATA(Cuts), VEC_DATA(Betas), nPreds, 1, nTerms, nMaxTerms, 3, 0);
	fflush(stdout);
*/
	// currently releasing unused arrays
	bx.clear();
	Residuals.clear();
	return 0;
}

//..............................................................................
int MedMars::Predict(float *x, float *&preds, int nsamples, int nftrs) const {

	// transfer x to double
	MedMat<double> xd;
	xd.load(x,nsamples,nftrs);
	double *xp = (double *)xd.data_ptr();
	double y;
	for (int i=0; i<nsamples; i++) {
		PredictEarth(&y, &xp[i*nftrs], BestSet, VEC_DATA(Dirs), VEC_DATA(Cuts), VEC_DATA(Betas), nPreds, 1, nTerms, nMaxTerms);
		preds[i] = (float)y;
	}
	return 0;
}

//..............................................................................
size_t MedMars::get_size() {

	size_t size = 0;

	size += sizeof(int);							// nMaxTerms;
	size += sizeof(int);							// nTerms;
	size += sizeof(int);							// nPreds;
	size += sizeof(bool) * nMaxTerms;				// BestSet
	size += sizeof(int) * nMaxTerms * nPreds;		// Dirs
	size += sizeof(double) * nMaxTerms * nPreds;	// Cuts
	size += sizeof(double) * nMaxTerms;				// Betas

	return size;
}

//..............................................................................
size_t MedMars::serialize(unsigned char *blob) {

	size_t ptr = 0;

	memcpy(blob+ptr,&nMaxTerms,sizeof(int)); ptr += sizeof(int);
	memcpy(blob+ptr,&nTerms,sizeof(int)); ptr += sizeof(int);
	memcpy(blob+ptr,&nPreds,sizeof(int)); ptr += sizeof(int);
	memcpy(blob+ptr,BestSet,sizeof(bool) * nMaxTerms); ptr += sizeof(bool) * nMaxTerms;
	memcpy(blob+ptr,VEC_DATA(Dirs),sizeof(int) * nMaxTerms * nPreds); ptr += sizeof(int) * nMaxTerms * nPreds;
	memcpy(blob+ptr,VEC_DATA(Cuts),sizeof(double) * nMaxTerms * nPreds); ptr += sizeof(double) * nMaxTerms * nPreds;
	memcpy(blob+ptr,VEC_DATA(Betas),sizeof(double) * nMaxTerms); ptr += sizeof(double) * nMaxTerms;

	return ptr;
}

//..............................................................................
size_t MedMars::deserialize(unsigned char *blob) {
	
	size_t ptr = 0;

	memcpy(&nMaxTerms,blob+ptr,sizeof(int)); ptr += sizeof(int);
	memcpy(&nTerms,blob+ptr,sizeof(int)); ptr += sizeof(int);
	memcpy(&nPreds,blob+ptr,sizeof(int)); ptr += sizeof(int);

	// allocating
	if (BestSet != NULL) free(BestSet);
	BestSet = (bool *)malloc(sizeof(bool)*nMaxTerms);
	Dirs.resize(nMaxTerms*nPreds);
	Cuts.resize(nMaxTerms*nPreds);
	Betas.resize(nMaxTerms);

	// reading arrays
	memcpy(BestSet, blob+ptr, sizeof(bool) * nMaxTerms); ptr += sizeof(bool) * nMaxTerms;
	memcpy(VEC_DATA(Dirs), blob+ptr, sizeof(int) * nMaxTerms * nPreds); ptr += sizeof(int) * nMaxTerms * nPreds;
	memcpy(VEC_DATA(Cuts), blob+ptr, sizeof(double) * nMaxTerms * nPreds); ptr += sizeof(double) * nMaxTerms * nPreds;
	memcpy(VEC_DATA(Betas), blob+ptr, sizeof(double) * nMaxTerms); ptr += sizeof(double) * nMaxTerms;

	return ptr;
}

// Printing
void MedMars::print(FILE *fp, const string& prefix) const {
	fprintf(fp,"%s: MedMars : nMaxTerms %d nTerms %d nPreds %d\n",prefix.c_str(), nMaxTerms, nTerms, nPreds) ;
}

// Prdictions per sample
int MedMars::n_preds_per_sample() const
{
	return 1;
}