#include "classifiers.h"
#include "QRF/QRF/QRF.h"

#define QRF_PRED_NTHREADS 8
int qrf_predict(double *tx, double *preds, int nrows, int ncols, unsigned char *model) 
{

	return qrf_predict(tx,preds,nrows,ncols,model,QRF_PRED_NTHREADS) ;
}


int qrf_predict(double *tx, double *preds, int nrows, int ncols, unsigned char *model, int nthreads) 
{
	QRF_Forest qf;
	int size;

	fprintf(stderr,"deserializing qrf_forest\n"); fflush(stderr);
	size = (int) qf.deserialize(model);

	qf.nthreads = nthreads;
	if (qf.score_samples_t(tx,ncols,nrows,preds) < 0) {
		fprintf(stderr,"failed predicting qrf forest\n"); fflush(stderr);
		return -1;
	}

	fprintf(stderr,"scored %d samples using qrf model, size is %d\n",nrows,size); fflush(stderr);
	return size;
}