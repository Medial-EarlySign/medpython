//
// MedAlgo - unified wrappers for prediction algorithms
//

#include "MedAlgo.h"
#include "MedXGB.h"
#include "MedDeepBit.h"
#include "MedUtils/MedUtils/MedIO.h"
#include "MedProcessTools/MedProcessTools/MedFeatures.h"
#include "MedLightGBM.h"

#include <thread>

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//=======================================================================================
// MedPredictor
//=======================================================================================
// Model types
MedPredictorTypes predictor_name_to_type(const string& model_name) {

	if (model_name == "linear_model")
		return  MODEL_LINEAR_MODEL;
	else if (model_name == "gdlm")
		return MODEL_GD_LINEAR;
	else if (model_name == "qrf")
		return MODEL_QRF ;
	else if (model_name == "knn")
		return MODEL_KNN;
	else if (model_name == "mars")
		return MODEL_MARS;
	else if (model_name == "gbm")
		return MODEL_GBM ;
	else if (model_name == "BP")
		return MODEL_BP;
	else if (model_name == "multi_class")
		return MODEL_MULTI_CLASS ;
	else if (model_name == "xgb")
		return MODEL_XGB;
	else if (model_name == "lasso")
		return MODEL_LASSO;
	else if (model_name == "micNet")
		return MODEL_MIC_NET;
	else if (model_name == "booster")
		return MODEL_BOOSTER;
	else if (model_name == "deep_bit")
		return MODEL_DEEP_BIT;
	else if (model_name == "lightgbm")
		return MODEL_LIGHTGBM;
	else
		return MODEL_LAST ;
}

// Initialization
//.......................................................................................
MedPredictor * MedPredictor::make_predictor(string model_type) {

	return make_predictor(predictor_name_to_type(model_type)) ;
}

//.......................................................................................
MedPredictor * MedPredictor::make_predictor(string model_type, string init_string) {

	return make_predictor(predictor_name_to_type(model_type),init_string);
}

//.......................................................................................
MedPredictor * MedPredictor::make_predictor(MedPredictorTypes model_type) {

	if (model_type == MODEL_LINEAR_MODEL)
		return new MedLM;
	else if (model_type == MODEL_GD_LINEAR)
		return new MedGDLM;
	else if (model_type == MODEL_QRF)
		return new MedQRF ;
	else if (model_type== MODEL_KNN)
		return new MedKNN;
	else if (model_type == MODEL_MARS)
		return new MedMars;
	else if (model_type == MODEL_GBM)
		return new MedGBM;
	else if (model_type == MODEL_BP)
		return new MedGBM;
	else if (model_type == MODEL_MULTI_CLASS)
		return new MedMultiClass ;
	else if (model_type == MODEL_XGB)
		return new MedXGB;
	else if (model_type == MODEL_LASSO)
		return new MedLasso;
	else if (model_type == MODEL_MIC_NET)
		return new MedMicNet;
	else if (model_type == MODEL_BOOSTER) 
		return new MedBooster;
	else if (model_type == MODEL_DEEP_BIT)
		return new MedDeepBit;
	else if (model_type == MODEL_LIGHTGBM)
		return new MedLightGBM;
	else
		return NULL ;

}
MedPredictor * MedPredictor::make_predictor(MedPredictorTypes model_type, string init_string) {

	MedPredictor *newPred = make_predictor(model_type);
	newPred->init_from_string(init_string);
	return newPred;
}
//.......................................................................................
int MedPredictor::init_from_string(string text) {
	
	cerr << "MedPredictor init from string (classifier type is " << classifier_type << " )\n";
	// parse text of the format "Name = Value ; Name = Value ; ..."

	if (classifier_type == MODEL_MIC_NET) {
		cerr << "But we are going to call mic net version directly\n";
		MedMicNet *mic = (MedMicNet *)this;
		cerr << "before\n";
		int rc =  mic->init_from_string(text);
		cerr << "after " << rc << "\n";
		return rc;
	}

	if (classifier_type == MODEL_BOOSTER) {
		cerr << "But we are going to call booster version directly\n";
		MedBooster *med_b = (MedBooster *)this;
		cerr << "before\n";
		int rc =  med_b->init_from_string(text);
		cerr << "after " << rc << "\n";
		return rc;
	}

	if (classifier_type == MODEL_LIGHTGBM) {
		MedLightGBM *med_light = (MedLightGBM *)this;
		return med_light->init_from_string(text);
	}

	// remove white spaces
	text.erase(remove_if(text.begin(), text.end(), ::isspace), text.end());

	map<string, string> init_map;
	if (initialization_text_to_map(text, init_map) == -1)
		return -1;

	for (auto rec : init_map)
		MLOG("Initializing predictor with \'%s\' = \'%s\'\n", rec.first.c_str(), rec.second.c_str());

	init(init_map);

	return 0;
}

//.......................................................................................
void MedPredictor::prepare_x_mat(MedMat<float> &x, vector<float> &wgts, int &nsamples, int &nftrs, bool transpose_needed)
{

	if ((transpose_needed && !x.transposed_flag) || (!transpose_needed && x.transposed_flag)) {
//		MLOG("transposing matrix\n");
		x.transpose();
	}

	if (x.transposed_flag) {
		nsamples = x.ncols;
		nftrs = x.nrows;
	} else {
		nsamples = x.nrows;
		nftrs = x.ncols;
	}
}

//.......................................................................................
int MedPredictor::learn(MedMat<float> &x, MedMat<float> &y, vector<float> &wgts)
{
	int nsamples, nftrs;

	// patch for micNet
	if (classifier_type == MODEL_MIC_NET) {
		MedMicNet *mic = (MedMicNet *)this;
		cerr << "running micNet learn()\n";
		return mic->learn(x, y);
	}
	// patch for booster
	if (classifier_type == MODEL_BOOSTER) {
		MedBooster *med_b = (MedBooster *)this;
		cerr << "running MedBooster learn()\n";
		return med_b->learn(x, y);
	}

	// ToDo : sanity check of sizes (nonzero, matching x,y dimensions)

	if (normalize_for_learn && !x.normalized_flag) {
		MERR("Learner Requires normalized matrix. Quitting\n");
		return -1;
	}

	prepare_x_mat(x,wgts,nsamples,nftrs,transpose_for_learn);
	if (normalize_y_for_learn && !y.normalized_flag)
		y.normalize();

	return Learn(x.data_ptr(),y.data_ptr(),VEC_DATA(wgts),nsamples,nftrs) ;
}
	
//.......................................................................................
int MedPredictor::learn(MedMat<float> &x, vector<float> &y, vector<float> &wgts) {

	int nsamples,nftrs ;

	if (normalize_for_learn && !x.normalized_flag) {
		MERR("Learner Requires normalized matrix. Quitting\n");
		return -1;
	}

	prepare_x_mat(x,wgts,nsamples,nftrs,transpose_for_learn);

	return Learn(x.data_ptr(),y.data(),VEC_DATA(wgts),nsamples,nftrs) ;
}
	
//.......................................................................................
int MedPredictor::learn(vector<float> &x, vector<float> &y, vector<float> &w, int n_samples, int n_ftrs) 
{
	return(Learn(VEC_DATA(x),VEC_DATA(y),VEC_DATA(w),n_samples,n_ftrs)) ;
}

//.......................................................................................
int MedPredictor::predict(MedMat<float> &x, vector<float> &preds) {

	int nsamples,nftrs ;
	vector<float> w ;

	// patch for micNet
	if (classifier_type == MODEL_MIC_NET) {
		MedMicNet *mic = (MedMicNet *)this;
		cerr << "running micNet predict()\n";
		return mic->predict(x, preds);
	}
	// patch for booster
	if (classifier_type == MODEL_BOOSTER) {
		MedBooster *med_b = (MedBooster *)this;
		cerr << "running MedBooster predict()\n";
		return med_b->predict(x, preds);
	}


	if (normalize_for_predict && !x.normalized_flag) {
		MERR("Predictor Requires normalized matrix. Quitting\n");
		return -1;
	}


	prepare_x_mat(x,w,nsamples,nftrs,transpose_for_predict);

	preds.resize(nsamples*n_preds_per_sample()) ;
	float *_preds = &(preds[0]);
//	MLOG("MedMat,vector call: preds size is %d n_preds_per_sample %d nsamples %d\n",preds.size(),n_preds_per_sample(),nsamples);
	return Predict(x.data_ptr(),_preds,nsamples,nftrs) ;
}

struct pred_thread_info {
	int id;
	int from_sample;
	int to_sample;
	float *preds;
	MedMat<float> *x;
	int nftrs;
	int n_preds_per_sample;
	int rc;
	int state;
};

void MedPredictor::predict_thread(void *p)
//void MedPredictor::predict_thread()
{
	
	pred_thread_info *tp = (pred_thread_info *)p;

	//MLOG("Start thread %d : from: %d to: %d\n",tp->id,tp->from_sample,tp->to_sample);
	float *x = &(tp->x->m[tp->from_sample * tp->nftrs]);
	float *preds = &(tp->preds[tp->from_sample * n_preds_per_sample()]);
	int nsamples = tp->to_sample - tp->from_sample + 1;
	int nftrs = tp->nftrs;

	tp->rc = Predict(x,preds,nsamples,nftrs);
	//MLOG("End thread %d : from: %d to: %d\n",tp->id,tp->from_sample,tp->to_sample);
	
	// signing job ended
	tp->state = 1;
}

//.......................................................................................
int MedPredictor::threaded_predict(MedMat<float> &x, vector<float> &preds, int nthreads) {

	int nsamples,nftrs ;
	vector<float> w ;

	if (transpose_for_predict) {
		MERR("!!!!!! UNSUPORTED !!!! --> Currently threaded_predict does not support transposed matrices for predictions");
		return -1;
	}

	if (normalize_for_predict && !x.normalized_flag) {
		MERR("Predictor Requires normalized matrix. Quitting\n");
		return -1;
	}

	prepare_x_mat(x,w,nsamples,nftrs,transpose_for_predict);
	preds.resize(nsamples*n_preds_per_sample()) ;
	
	int th_nsamples = nsamples/nthreads;
	vector<pred_thread_info> tp(nthreads);
	for (int i=0; i<nthreads; i++) {
		tp[i].id = i;
		tp[i].from_sample = i*th_nsamples;
		tp[i].to_sample = min((i+1)*th_nsamples - 1,nsamples-1);
		tp[i].preds = VEC_DATA(preds);
		tp[i].x = &x;
		tp[i].nftrs = nftrs;
		tp[i].rc = 0;
		tp[i].state = 0;
	}


	// sending threads
	vector<thread> th_handle(nthreads);
	for (int i=0; i<nthreads; i++) {
//		MLOG("Sending Thread %d\n",i);
//		th_handle[i] = std::thread(&MedPredictor::predict_thread, (void *)&tp[i]);
		th_handle[i] = std::thread(&MedPredictor::predict_thread,this, (void *)&tp[i]);
	}

	int n_state = 0;
	while (n_state < nthreads) {
		this_thread::sleep_for(chrono::milliseconds(10));
		n_state = 0;
		for (int i=0; i<nthreads; i++)
			n_state += tp[i].state;
	}
	for (int i=0; i<nthreads; i++)
		th_handle[i].join();

	for (int i=0; i<nthreads; i++)
		if (tp[i].rc != 0)
			return -1;
	return 0;

}


//.......................................................................................
int MedPredictor::predict(vector<float> &x, vector<float> &preds, int n_samples, int n_ftrs) {
	preds.resize(n_samples*n_preds_per_sample()) ;
	float *_preds = &(preds[0]);
	return Predict(VEC_DATA(x),_preds,n_samples,n_ftrs) ;
}


// (De)Serialize
//.......................................................................................
size_t MedPredictor::get_predictor_size()  {
	return sizeof(classifier_type) + get_size() ;
}
	
//.......................................................................................
size_t MedPredictor::predictor_serialize(unsigned char *blob) {

	size_t ptr = 0 ;
 	memcpy(blob+ptr,&classifier_type,sizeof(MedPredictorTypes)) ; ptr += sizeof(MedPredictorTypes) ;
	ptr += serialize(blob+ptr) ;

	return ptr ;
}

//.......................................................................................

void MedPredictor::print(FILE *fp, const string& prefix) { fprintf(fp,"%s : No Printing method defined\n",prefix.c_str());}

//===========================================================================================
// MedPredictor MedFeaturesData API
//===========================================================================================
// Cross Validation
int MedPredictor::cross_validate_splits(MedFeaturesData& data) {

	data.split_preds.resize(data.nsplits) ;
	data.split_preds_on_train.resize(data.nsplits);
	
	for (int isplit=0; isplit<data.nsplits; isplit++) {
		if (learn(data,isplit) == -1) {
			fprintf(stderr,"Learning failed\n") ;
			return -1;
		}

		if (data.predict_on_train) {
			if (predict_on_train(data, isplit) == -1) {
				fprintf(stderr, "Predict on train failed\n");
				return -1;
			}
		}

		if (predict(data,isplit) == -1) {
			fprintf(stderr,"Predict failed\n") ;
			return -1 ;
		}
	}

	return 0 ;
}

void MedPredictor::build_learning_x_mat_for_split(MedFeaturesData& ftrs_data, vector<float>& signal, int isplit, MedMat<float>& x) {
	x.normalized_flag = 1;

	for (int icol = 0; icol<x.ncols; icol++) {
		string& name = ftrs_data.signals[icol];

		int split_irow = 0;
		for (int irow = 0; irow<ftrs_data.splits.size(); irow++) {
			if (ftrs_data.splits[irow] != isplit) 
				signal[split_irow++] = ftrs_data.data[name][irow];
		}

		ftrs_data.cleaners[name][isplit].clear(signal);
		ftrs_data.cleaners[name][isplit].normalize(signal);

		for (unsigned irow = 0; irow<signal.size(); irow++)
			x(irow, icol) = signal[irow];
	}
}
// Learning ....
int MedPredictor::learn(MedFeaturesData& ftrs_data, int isplit) {
	
	MedMat<float> x((int) ftrs_data.get_learning_nrows(isplit),(int) ftrs_data.signals.size()) ;
	MedMat<float> y(x.nrows,1) ;
	vector<float> signal(x.nrows);

	// Build X
	build_learning_x_mat_for_split(ftrs_data, signal, isplit, x);

	// Build Y
	int split_irow = 0 ;
	for (int irow=0; irow<ftrs_data.splits.size(); irow++) {
		if (ftrs_data.splits[irow] != isplit)
			signal[split_irow++] = ftrs_data.label[irow] ;
	}
	if (ftrs_data.label_cleaners.size() > 0) {
		ftrs_data.label_cleaners[isplit].clear(signal);
		ftrs_data.label_cleaners[isplit].normalize(signal);
	}

	for (unsigned irow=0; irow<signal.size(); irow++)
		y(irow,0) = signal[irow] ;
	y.normalized_flag = 1 ;

	// Learn
	if (learn(x,y) == -1)
		return -1;

	ftrs_data.n_preds_per_sample = n_preds_per_sample();

	return 0;
}

int MedPredictor::learn(MedFeaturesData& ftrs_data) {
	
	// Build X
	MedMat<float> x((int) ftrs_data.label.size(), (int) ftrs_data.signals.size()) ;
	MedMat<float> y(x.nrows,1) ;
	x.normalized_flag = 1 ;

	vector<float> signal(x.nrows) ;
	for (int icol=0; icol<x.ncols; icol++) {
		string& name = ftrs_data.signals[icol] ;
		
		signal = ftrs_data.data[name] ;
		ftrs_data.cleaners[name][ftrs_data.nsplits].clear(signal) ;
		ftrs_data.cleaners[name][ftrs_data.nsplits].normalize(signal) ;

		for (unsigned irow=0; irow<signal.size(); irow++)
				x(irow,icol) = signal[irow] ;
	}

	// Build Y
	signal = ftrs_data.label ;
	ftrs_data.label_cleaners[ftrs_data.nsplits].clear(signal) ;
	ftrs_data.label_cleaners[ftrs_data.nsplits].normalize(signal) ;

	for (unsigned irow=0; irow<signal.size(); irow++)
		y(irow,0) = signal[irow] ;
	y.normalized_flag = 1;

	// Learn
	if (learn(x,y) == -1)
		return -1;

	return 0;
}

int MedPredictor::learn(MedFeatures& ftrs_data) {

	// Build X
	MedMat<float> x;
	ftrs_data.get_as_matrix(x);

	MLOG("MedPredictor::learn() from MedFeatures, got train matrix of %d x %d\n", x.nrows, x.ncols);

	// Labels
	MedMat<float> y(x.nrows, 1);
	for (int i = 0; i < y.nrows; i++)
		y(i,0) = ftrs_data.samples[i].outcome ;
	return learn(x, y);
}

// Predicting
int MedPredictor::predict_on_train(MedFeaturesData& ftrs_data, int isplit) {
	MedMat<float> x((int)ftrs_data.get_learning_nrows(isplit), (int)ftrs_data.signals.size());
	vector<float> signal(x.nrows);

	// Build X
	build_learning_x_mat_for_split(ftrs_data, signal, isplit, x);

	// Predict
	if (predict(x, ftrs_data.split_preds_on_train[isplit]) == -1)
		return -1;

	if (ftrs_data.label_cleaners.size() > 0) {
		for (unsigned int i = 0; i < ftrs_data.split_preds_on_train[isplit].size(); i++) {
			ftrs_data.split_preds_on_train[isplit][i] += ftrs_data.label_cleaners[isplit].mean;
		}
	}

	ftrs_data.n_preds_per_sample = n_preds_per_sample();

	return 0;
}


// Predicting
int MedPredictor::predict(MedFeaturesData& ftrs_data, int isplit) {

	// Build X
	MedMat<float> x((int) ftrs_data.get_testing_nrows(isplit), (int) ftrs_data.signals.size()) ;

	vector<float> signal(x.nrows) ;
	for (int icol=0; icol<ftrs_data.signals.size(); icol++) {
		string& name = ftrs_data.signals[icol] ;

		int split_irow = 0 ;
		for (int irow=0; irow<ftrs_data.splits.size(); irow++) {
			if (ftrs_data.splits[irow] == isplit)
				signal[split_irow++] = ftrs_data.data[name][irow] ;
		}

		ftrs_data.cleaners[name][isplit].clear(signal) ;
		ftrs_data.cleaners[name][isplit].normalize(signal) ;

		for (unsigned irow=0; irow<signal.size(); irow++) 
			x(irow,icol) = signal[irow] ;

	}

	x.normalized_flag = 1 ;

	// Predict
	if (predict(x,ftrs_data.split_preds[isplit])  == -1)
		return -1 ;

	if (ftrs_data.label_cleaners.size() > 0) {
		for (unsigned int i = 0; i < ftrs_data.split_preds[isplit].size(); i++) {
			ftrs_data.split_preds[isplit][i] += ftrs_data.label_cleaners[isplit].mean;
		}
	}

	ftrs_data.n_preds_per_sample = n_preds_per_sample();

	return 0 ;

}

int MedPredictor::predict(MedFeaturesData& ftrs_data) {
	
	// Build X
	MedMat<float> x((int) ftrs_data.label.size(),(int) ftrs_data.signals.size()) ;
	MedMat<float> y(x.nrows,1) ;

	vector<float> signal(x.nrows) ;
	for (int icol=0; icol<x.ncols; icol++) {
		string& name = ftrs_data.signals[icol] ;
		
		signal = ftrs_data.data[name] ;
		ftrs_data.cleaners[name][ftrs_data.nsplits].clear(signal) ;

		for (unsigned irow=0; irow<signal.size(); irow++)
				x(icol,irow) = signal[irow] ;

	}
	x.normalized_flag = 1 ;

	// Predict
	ftrs_data.n_preds_per_sample = n_preds_per_sample();
	return predict(x,ftrs_data.preds)  ;
}

int MedPredictor::predict(MedFeatures& ftrs_data) {

	// Build X
	MedMat<float> x;
	ftrs_data.get_as_matrix(x);

	// Predict
	vector<float> preds;
	if (predict(x, preds) < 0)
		return -1;

	int n = n_preds_per_sample();
	ftrs_data.samples.resize(preds.size()/n);
	for (int i = 0; i < x.nrows; i++) {
		ftrs_data.samples[i].prediction.resize(n);
		for (int j = 0; j < n; j++)
			ftrs_data.samples[i].prediction[j] = preds[i*n + j];
	}

	return 0;
}

int MedPredictor::read_from_file(const string &fname) // read and deserialize model
{
	unsigned char *blob;
	unsigned long long size;

	if (read_binary_data_alloc(fname,blob,size) < 0) {
		MERR("Error reading model from file %s\n",fname.c_str());
		return -1;
	}

	deserialize(blob);
	if (size > 0) delete[] blob;
	return 0;
}

int MedPredictor::write_to_file(const string &fname)  // serialize model and write to file
{
	unsigned char *blob;
	unsigned long long size;

	size = get_size();
	vector<unsigned char> local_buf(size+1);
	//blob = new unsigned char[size];
	blob = &local_buf[0];

	//MLOG("write_to_file 1: size %lld blob %llx\n", size, blob);
	serialize(blob);

	if (write_binary_data(fname,blob,size) < 0) {
		MERR("Error writing model to file %s\n",fname.c_str());
		return -1;
	}

	//MLOG("write_to_file 1: size %lld blob %llx\n", size, blob);
	//MLOG("Before Release blob size = %lld\n",size);
	//if (size > 0) delete[] blob;
	return 0;
}
