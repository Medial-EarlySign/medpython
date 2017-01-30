//
// MedGBM - wrapper for the GBM Library
//

#define _CRT_SECURE_NO_WARNINGS

#include "MedAlgo.h"
#include "gbm/gbm/gbm_utils.h"


#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//..............................................................................
void init_default_gbm_params(MedGBMParams& _params) {

	_params.bag_p = MED_GBM_DEF_BAG_P ;
	_params.depth = MED_GBM_DEF_DEPTH ;
	_params.min_obs_in_node = MED_GBM_DEF_MIN_OBS_IN_NODE ;
	_params.ntrees = MED_GBM_DEF_NTREES ;
	_params.shrinkage = MED_GBM_DEF_SHRINKAGE ;
	_params.take_all_pos = MED_GBM_DEF_TAKE_ALL_POS ;
}

//..............................................................................
GBM_LossFunctions MedGBM::get_loss_function(string name) {

	boost::algorithm::to_lower(name);
	if (name == "adaBoost")
		return GBM_Loss_AdaBoost;
	else if (name == "bernulli")
		return GBM_Loss_Bernulli;
	else if (name == "gaussian")
		return GBM_Loss_Gaussian;
	else if (name == "laplace")
		return GBM_Loss_Laplace;
	else if (name == "poisson")
		return GBM_Loss_Poisson;
	else if (name == "quantile")
		return GBM_Loss_Quantile;
	else
		return GBM_Loss_Last;
}
	
//..............................................................................
int MedGBM::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "bag_p") params.bag_p = stod(entry.second);
		else if (field == "depth") params.depth = stoi(entry.second);
		else if (field == "min_obs_in_node") params.min_obs_in_node = stoi(entry.second);
		else if (field == "ntrees") params.ntrees = stoi(entry.second);
		else if (field == "shrinkage") params.shrinkage = stod(entry.second);
		else if (field == "take_all_pos") params.take_all_pos = (bool)(stoi(entry.second) != 0);
		else if (field == "loss_function") loss_function = get_loss_function(entry.second);
		else if (field == "alpha_quantile") alpha_quantile = stod(entry.second);
		else MLOG("Unknonw parameter \'%s\' for GBM\n", field.c_str());
	}

	return 0;
}

//..............................................................................
void MedGBM::init_defaults()
{
	classifier_type = MODEL_GBM ;
	transpose_for_learn = false;
	normalize_for_learn = false;
	normalize_y_for_learn = false;
	transpose_for_predict = false;
	normalize_for_predict = false;

	init_default_gbm_params(params) ;

	predict_ntrees = params.ntrees ;
	loss_function = GBM_Loss_AdaBoost ;
	alpha_quantile = 0.25;

	gbm_model.ntrees = 0 ;
}

//..............................................................................
int MedGBM::init(void *_in_params) 
{
	init_defaults();

	MedGBMParams *in_params = (MedGBMParams *) _in_params ;

	params = (*in_params);
	predict_ntrees = params.ntrees ;

	return 0 ;
}

//..............................................................................
MedGBM::MedGBM() 
{
	init_defaults();
}

//..............................................................................
MedGBM::MedGBM(MedGBMParams& _in_params) 
{
	init((void *) &_in_params);
}

//..............................................................................
MedGBM::MedGBM(void *_in_params) 
{
	init(_in_params);
}

//..............................................................................
MedGBM::~MedGBM() 
{
	if (gbm_model.ntrees)
		clear_gbm_info(&gbm_model) ;
}

//..............................................................................
int MedGBM::Learn (float *x, float *y, int nsamples, int nftrs) {

	vector<float> weights(nsamples,1.0) ;
	return Learn(x,y,&(weights[0]),nsamples,nftrs) ;
}

//..............................................................................
int MedGBM::Learn (float *x, float *y, float *w, int nsamples, int nftrs) {

	if (w == NULL) 
		return (Learn(x,y,nsamples,nftrs));

	// transpose x and make it a double to match inner package req
	MedMat<double> xd;
	xd.load_transposed(x,nsamples,nftrs);

	// transfer y to double
	MedMat<double> yd(nsamples,1) ;
	for (int i=0; i<nsamples; i++)
		yd(i,0) = (y[i] > 0) ? 1.0 : -1.0 ;

	// transfer w to double
	MedMat<double> wd;
	wd.load(w,nsamples,1);

	// Learn
	if (get_gbm_predictor(xd.data_ptr(),yd.data_ptr(),wd.data_ptr(),nsamples,nftrs,params.shrinkage,params.bag_p,params.take_all_pos,params.ntrees,params.depth,params.min_obs_in_node,&gbm_model, loss_function, alpha_quantile) < 0) {
		fprintf(stderr,"GBM learning failed\n") ;
		return -1 ;
	}

	return 0;
}

//..............................................................................
int MedGBM::Predict(float *x, float *&preds, int nsamples, int nftrs) {

	// transpose x and make it a double to match inner package req
	MedMat<double> xd;
	xd.load_transposed(x,nsamples,nftrs);

	// Predict
	MedMat<double> predsd(nsamples,1) ;
	if (gbm_predict(xd.data_ptr(),nsamples,nftrs,predict_ntrees,&gbm_model,predsd.data_ptr())==-1) {
		fprintf(stderr,"GBM prediction failed\n") ;
		return -1 ;
	}

	for (int i=0; i<nsamples; i++)
		preds[i] = (float) predsd(i,0) ;

	return 0;
}

//..............................................................................
size_t MedGBM::get_size() {

	return get_gbm_info_size(&gbm_model) ;

}

//..............................................................................
size_t MedGBM::serialize(unsigned char *blob) {

	return (size_t) gbm_serialize(&gbm_model,blob) ;

}

//..............................................................................
size_t MedGBM::deserialize(unsigned char *blob) {

	return (size_t) gbm_deserialize(blob, &gbm_model) ;

}

// Printing
void MedGBM::print(FILE *fp, const string& prefix) {
	fprintf(fp, "%s: MedGBM (ntrees = %d)\n", prefix.c_str(), gbm_model.ntrees);
	//write_full_gbm_info(prefix,&gbm_model,stdout) ;
		
}
