//
// MedMultiClass - multi-class prediction
//

#include <MedAlgo/MedAlgo/MedAlgo.h>
#include <MedAlgo/MedAlgo/MedMultiClass.h>


#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//..............................................................................
void  MedMultiClass::set_internal_method(MedPredictorTypes method) {
	
	params.method = method ;

	MedPredictor *temp_predictor = MedPredictor::make_predictor(method) ;
	normalize_for_learn = temp_predictor->normalize_for_learn ;
	transpose_for_learn = temp_predictor->transpose_for_learn ;
	normalize_y_for_learn = temp_predictor->normalize_y_for_learn ;
	transpose_for_predict = temp_predictor->transpose_for_predict ;
	normalize_for_predict = temp_predictor->normalize_for_predict ;
}

//..............................................................................
void MedMultiClass::init_defaults()
{
	classifier_type = MODEL_MULTI_CLASS ;
	params.method = MODEL_LAST ;

	params.internal_params = NULL ;
	params.class_values.clear() ;

	internal_predictors.clear() ;
}

//..............................................................................
int MedMultiClass::init(void *_in_params) 
{
	init_defaults();

	MedMultiClassParams *in_params = (MedMultiClassParams *) _in_params ;
	 
	params = (*in_params);

	return 0 ;
}

//..............................................................................
MedMultiClass::MedMultiClass() 
{
	init_defaults();
}

//..............................................................................
MedMultiClass::MedMultiClass(MedMultiClassParams& _in_params) 
{
	init((void *) &_in_params);
}

//..............................................................................
MedMultiClass::MedMultiClass(void *_in_params) 
{
	init(_in_params);
}

//..............................................................................
MedMultiClass::~MedMultiClass() 
{

	if (params.internal_params != NULL)
		free(params.internal_params) ;

}

//..............................................................................
int MedMultiClass::init_classifiers() {

	if (params.multi_class_type == MULTI_CLASS_ONE_VS_ALL) {
		internal_predictors.resize(params.class_values.size()) ;
	} else {
		MERR("Unknown multi-class-type %d\n",params.multi_class_type) ;
		return -1 ;
	}

	for (unsigned int i=0; i<internal_predictors.size(); i++) {
		if (init_classifier(i) < 0) {
			MERR("Initialization of classifier %d failed\n",i) ;
			return -1 ; ;
		}
	}

	return 0 ;

}
	
int MedMultiClass::init_classifier(int index) {

	internal_predictors[index] = MedPredictor::make_predictor(params.method) ;
	if (internal_predictors[index] == NULL) {
		MERR("Initiailization of predictor of type %d failed\n",params.method) ;
		return -1 ;
	}

	if (internal_predictors[index]->init(params.internal_params))
		return -1 ;

	return 0 ;
}


//..............................................................................
int MedMultiClass::Learn(float *x, float *y, int nsamples, int nftrs) {

	vector<float> weights(nsamples,1.0) ;
	return Learn(x,y,&(weights[0]),nsamples,nftrs) ;
}

//..............................................................................
int MedMultiClass::Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {

	if (params.method == MODEL_LAST || params.multi_class_type == MULTI_CLASS_LAST) {
		MERR("MedMultiClass:Learn not properly initialized\n") ;
		return -1 ;
	}

	// Get values of Y per class + initialize models
	if (params.class_values.empty()) {
		map<float,float> y_values ;
		for (int i=0; i<nsamples; i++)
			y_values[y[i]] = 1 ;

		for (auto it = y_values.begin(); it != y_values.end(); it++)
			params.class_values.push_back(it->first) ;
		
		init_classifiers() ;

	} else if (internal_predictors.empty()) {
		init_classifiers() ;
	}

	// Learn 
	if (params.multi_class_type == MULTI_CLASS_ONE_VS_ALL) {

		for (unsigned int iclass=0; iclass<params.class_values.size(); iclass++) {

			// Create local Y
			vector<float> class_y(nsamples) ;

			for (int i=0; i<nsamples; i++) {
				if (y[i] == params.class_values[iclass])
					class_y[i] = 1.0 ;
				else
					class_y[i] = -1.0 ;
			}

			// Learn Model
			if (internal_predictors[iclass]->Learn(x,&(class_y[0]),w,nsamples,nftrs) < 0) {
				MERR("Learning failed for class %d\n",iclass) ;
				return -1 ;
			}
		}
	}

	return 0 ;

}

//..............................................................................
int MedMultiClass::n_preds_per_sample() const {
	
	if (params.multi_class_type == MULTI_CLASS_ONE_VS_ALL) {
		return (int) internal_predictors.size() ;
	} else
		return -1 ;
}

//..............................................................................
int MedMultiClass::Predict(float *x, float *&preds, int nsamples, int nftrs) const {

	vector<float> class_pred(nsamples) ;
	float *_preds = &(class_pred[0]);
	int npreds = n_preds_per_sample() ;

	if (params.multi_class_type == MULTI_CLASS_ONE_VS_ALL) {

		for (unsigned int i=0; i<internal_predictors.size(); i++) {
			// Local Predictions
			internal_predictors[i]->Predict(x,_preds,nsamples,nftrs) ;
			for (int j=0; j<nsamples; j++)
				preds[j*npreds + i] = class_pred[j] ;
		}
	}

	return 0 ;

}

// Print
//..............................................................................
void MedMultiClass::print(FILE *fp, const string& prefix, int level) const {

	fprintf(fp,"%s: MedMultiClass - Number of predictors  = %d\n",prefix.c_str(),(int)internal_predictors.size()) ;
	if (level > 0) {
		for (unsigned int i = 0; i < internal_predictors.size(); i++) {
			fprintf(fp, "%s: MedMultiClass Predictor %d\n", prefix.c_str(), i);
			internal_predictors[i]->print(fp, prefix);
		}
	}
}