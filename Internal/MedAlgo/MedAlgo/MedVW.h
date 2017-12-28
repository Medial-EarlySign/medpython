#ifndef __MED_VOWPAL_WABBIT__H_
#define __MED_VOWPAL_WABBIT__H_
#include "MedAlgo.h"
#if NEW_COMPLIER
#include "vowpal_wabbit/vowpalwabbit/vw.h"
#include "vowpal_wabbit/vowpalwabbit/vwdll.h"

class MedVW : public MedPredictor {
public:
	void init_defaults();

	// Function
	MedVW();

	int init_from_string(string text);

	int Learn(float *x, float *y, float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs);

	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	void save_model(const string &path);
	void load_model(const string &path);

private:
	vw* _v;

};

#endif
#endif
