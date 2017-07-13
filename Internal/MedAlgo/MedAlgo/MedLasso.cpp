#define _CRT_SECURE_NO_WARNINGS

#include "MedAlgo.h"

#define LOCAL_SECTION LOG_MEDALGO
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//============================================================================
// LASSO Model
//============================================================================
//..............................................................................

void init_default_lasso_params(MedLassoParams& _params) {
	_params.lambda = LASSO_LAMBDA;
	_params.num_iterations = LASSO_NITER;
}

int MedLasso::init(map<string, string>& mapper) {

	init_defaults();

	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "lambda") params.lambda = stod(entry.second);
		else if (field == "num_iterations") params.num_iterations = stoi(entry.second);
		else MLOG("Unknonw parameter \'%s\' for Lasso\n", field.c_str());
	}

	return 0;
}

void MedLasso::init_defaults()
{
	classifier_type = MODEL_LASSO;
	transpose_for_learn = false;
	transpose_for_predict = false;
	normalize_for_learn = true;
	normalize_y_for_learn = true;
	normalize_for_predict = true;
	init_default_lasso_params(params);
}


void MedLasso::Initb() {
	b = vector<float>();
	b0 = 0;
}

//..............................................................................
int MedLasso::init(void *_in_params)
{
	init_defaults();
	MedLassoParams *in_params = (MedLassoParams *)_in_params;
	params.lambda = in_params->lambda;
	params.num_iterations = in_params->num_iterations;
	Initb();
	return 0;
}
//..............................................................................
MedLasso::MedLasso()
{
	init_defaults();
	params.lambda = LASSO_LAMBDA;
	params.num_iterations = LASSO_NITER;
	Initb();
}


//..............................................................................
MedLasso::MedLasso(MedLassoParams& _in_params)
{
	init((void *)&_in_params);
}

//..............................................................................
MedLasso::MedLasso(void *_in_params)
{
	init(_in_params);
}


//..............................................................................
int MedLasso::Learn(float *x, float *y, int nsamples, int nftrs) {
	vector<float> weights(nsamples, 1.0);
	return Learn(x, y, &(weights[0]), nsamples, nftrs);
}

double perform_lasso_iteration(double* xk_train, vector<double>& r, float bk, int nrow_train, double lambda) {
	double bk_hat = 0;
	for (int i = 0; i < nrow_train; i++)
		bk_hat += r[i] * xk_train[i];
	bk_hat /= nrow_train;
	bk_hat += bk;
	if (bk_hat - lambda > 0)
		bk_hat -= lambda;
	else if (bk_hat + lambda < 0)
		bk_hat += lambda;
	else
		bk_hat = 0;
	if (bk_hat != bk) {
		for (int i = 0; i < nrow_train; i++)
			r[i] += xk_train[i] * (bk - bk_hat);
	}
	return bk_hat;
}

void lasso_regression(double **trainx, double *y, vector<float>& b, int nrow_train, int n_ftrs, double lambda, int num_iterations) {
	vector<double> r(nrow_train);
	for (int i = 0; i < nrow_train; i++)
		r[i] = y[i];
	for (int it = 0; it < num_iterations; it++)
		for (int k = 0; k < n_ftrs; k++)
			b[k] = (float)perform_lasso_iteration(trainx[k], r, b[k], nrow_train, lambda);
}


void initialize_vars(double **&trainx, double *&y, float *x_in, float *y_in, float *w, vector<float>& b, int nrow_train, int n_ftrs) {
	trainx = (double **)malloc(n_ftrs * sizeof(double*));
	for (int j = 0; j < n_ftrs; j++)
		trainx[j] = (double *)malloc(nrow_train * sizeof(double));
	y = (double *)malloc(nrow_train * sizeof(double));
	
	b.resize(n_ftrs);
	for (int i = 0; i < n_ftrs; i++) b[i] = 0;

	for (int i = 0; i < nrow_train; i++)
		for (int j = 0; j < n_ftrs; j++)
			trainx[j][i] = sqrt(w[i]) * x_in[i * n_ftrs + j];
	for (int i = 0; i < nrow_train; i++)
		y[i] = sqrt(w[i]) * y_in[i];
}

//void initialize_test(double **&testx, float *test_in, int nrow_test, int n_ftrs) {
//	testx = (double**)malloc(n_ftrs * sizeof(double *));
//	for()
//	for (int i = 0; i < nrow_test; i++) {
//		testx[j] = (double *)malloc(nrow_test * sizeof(double));
//		for (int j = 0; j < n_ftrs; j++)
//			testx[i][j] = test_in[i * n_ftrs + j];
//	}
//}


//..............................................................................
int MedLasso::Learn(float *x, float *y, float *w, int nsamples, int nftrs) {
	MedLasso::n_ftrs = nftrs;
	if (w == NULL)
		return (Learn(x, y, nsamples, nftrs));
	double **trainx;
	double *y1;
	initialize_vars(trainx, y1, x, y, w, b, nsamples, n_ftrs);
	lasso_regression(trainx, y1, b, nsamples, n_ftrs, params.lambda, params.num_iterations);
	return 0;
}

//..............................................................................
int MedLasso::Predict(float *x, float *&preds, int nsamples, int nftrs) {
	return Predict(x, preds, nsamples, nftrs, 0);
}

//..............................................................................
int MedLasso::Predict(float *x, float *&preds, int nsamples, int nftrs, int transposed_flag) {

	if (preds == NULL)
		preds = new float[nsamples];

	memset(preds, 0, nsamples*sizeof(float));
	/*double **testx;
	initialize_test(testx, x, nsamples, nftrs);*/

	if (transposed_flag) {
		for (int j = 0; j < nftrs; j++) {
			for (int i = 0; i < nsamples; i++)
				preds[i] += b[j] * x[j * nsamples + i];
		}
	}
	else {
		for (int i = 0; i < nsamples; i++) {
			for (int j = 0; j < nftrs; j++)
				preds[i] += b[j] * x[i * nftrs + j];
		}
	}

	for (int i = 0; i < nsamples; i++)
		preds[i] += b0;

	return 0;
}

//..............................................................................
size_t MedLasso::get_size() {
	return sizeof(int) + (n_ftrs + 1) * sizeof(float);
}

//..............................................................................
size_t MedLasso::serialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(blob + ptr, &n_ftrs, sizeof(int)); ptr += sizeof(int);
	memcpy(blob + ptr, &b0, sizeof(float)); ptr += sizeof(float);
	memcpy(blob + ptr, &(b[0]), n_ftrs*sizeof(float)); ptr += n_ftrs*sizeof(float);

	return ptr;
}

size_t MedLasso::deserialize(unsigned char *blob) {

	size_t ptr = 0;
	memcpy(&n_ftrs, blob + ptr, sizeof(int)); ptr += sizeof(int);
	memcpy(&b0, blob + ptr, sizeof(float)); ptr += sizeof(float);

	b.resize(n_ftrs);
	memcpy(&(b[0]), blob + ptr, n_ftrs*sizeof(float)); ptr += n_ftrs*sizeof(float);

	return ptr;
}

//..............................................................................
int MedLasso::denormalize_model(float *f_avg, float *f_std, float label_avg, float label_std)
{
	float new_b0;
	vector<float> new_b(n_ftrs);

	new_b0 = b0*label_std + label_avg;
	fill(new_b.begin(), new_b.end(), (float)0);
	for (int j = 0; j < n_ftrs; j++) {
		new_b[j] = label_std*b[j] / f_std[j];
		new_b0 -= label_std*f_avg[j] * b[j] / f_std[j];
	}

	b0 = new_b0;
	for (int j = 0; j < n_ftrs; j++)
		b[j] = new_b[j];

	transpose_for_predict = true;
	normalize_for_predict = false;
	return 0;
}







//..............................................................................
void MedLasso::print(FILE *fp, const string& prefix) {

	fprintf(fp, "%s : Linear Model : Nftrs = %d\n", prefix.c_str(), n_ftrs);
	fprintf(fp, "%s : Linear Model b0 = %f\n", prefix.c_str(), b0);

	for (int i = 0; i < n_ftrs; i++)
		fprintf(fp, "%s : Linear Model b[%d] = %f\n", prefix.c_str(), i, b[i]);
}