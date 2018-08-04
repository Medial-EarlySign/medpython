#include "MedVW.h"
#include "Logger/Logger/Logger.h"

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL	LOG_DEF_LEVEL

#if NEW_COMPLIER
#pragma message ( "You Are Compiling with the new compiler" )
#include "External/vowpal_wabbit/vowpalwabbit/ezexample.h"
#include "External/vowpal_wabbit/vowpalwabbit/parse_regressor.h"

void MedVW::init_defaults() {
	transpose_for_learn = false;transpose_for_predict = false;
	normalize_for_learn = false;
	normalize_y_for_learn = false;
	normalize_for_predict = false;

	classifier_type = MODEL_VW;
	_v = NULL;
}

MedVW::MedVW() {
	init_defaults();
}

int MedVW::init_from_string(string text) {
	_v = VW::initialize(text);
	_v->vw_is_main = true;
	_v->sd->min_label = (float)INT_MIN;
	_v->sd->max_label = (float)INT_MAX;

	return 0;
}

int MedVW::Learn(float *x, float *y, const float *w, int nsamples, int nftrs) {
	char buff[500];
	for (int i = 0; i < nsamples; ++i)
	{
		string example_string = to_string(int(y[i])) + " |";
		for (size_t k = 0; k < nftrs; ++k)
		{
			snprintf(buff, 1000, " FTR_%d:%f", (int)k, x[i*nftrs + k]);
			example_string += buff;
		}

		example *e = VW::read_example(*_v, example_string);
		e->example_counter = (int)i;
		_v->learn(e); //v->l->learn(*e);
		_v->l->finish_example(*_v, *e);
	}

	//_v->l->end_examples();
	return 0;
}

int MedVW::Predict(float *x, float *&preds, int nsamples, int nftrs) const {
	char buff[500];
	for (size_t i = 0; i < nsamples; ++i)
	{
		string example_string = "|";
		for (size_t k = 0; k < nftrs; ++k)
		{
			snprintf(buff, 1000, " FTR_%d:%f", (int)k, x[i*nftrs + k]);
			example_string += buff;
		}

		//ezexample ex(_v);
		//ex(vw_namespace('s'))(example_string);
		//preds[i] = ex.predict();

		example *e = VW::read_example(*_v, example_string);
		e->in_use = false;
		e->test_only = true;

		_v->l->predict(*e);
		preds[i] = e->pred.scalar;
		//preds[i] = e->partial_prediction;
	}
	return 0;
}

int MedVW::write_to_file(const string &path) {
	_v->save_per_pass = false;
	save_predictor(*_v, path, 0);
	return 0;
}

int MedVW::read_from_file(const string &path) {
	_v->l = NULL;
	io_buf buf;
	buf.open_file(path.c_str(), _v->stdin_off, io_buf::READ);
	save_load_header(*_v, buf, true, false);
	return 0;
}

size_t MedVW::get_size() {
	return 0;
}
size_t MedVW::serialize(unsigned char *blob) {

	MTHROW_AND_ERR("ERROR: Use save_model directly\n");
}
size_t MedVW::deserialize(unsigned char *blob) {
	MTHROW_AND_ERR("ERROR: Use load_model directly\n");
}
#endif