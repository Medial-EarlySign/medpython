#ifndef __MED_VOWPAL_WABBIT__H_
#define __MED_VOWPAL_WABBIT__H_
#include <MedAlgo/MedAlgo/MedAlgo.h>

#if NEW_COMPLIER
#include <vowpal_wabbit/vowpalwabbit/vw.h>
#include <vowpal_wabbit/vowpalwabbit/vwdll.h>

/**
* Vowpal Wabbit library who have been developed in Yahoo Research, It has
* various machine learning models specially for text.
* It may be initialize by string
*/
class MedVW : public MedPredictor {
public:
	void init_defaults();

	// Function
	MedVW();

	///The initialization parameters - please look at vowpal wabbit documentation
	int init_from_string(string text);

	int Learn(float *x, float *y, const float *w, int nsamples, int nftrs);
	int Predict(float *x, float *&preds, int nsamples, int nftrs) const;

	ADD_CLASS_NAME(MedVW)
	size_t get_size();
	size_t serialize(unsigned char *blob);
	size_t deserialize(unsigned char *blob);

	///The Serialization function
	int write_to_file(const string &path);
	///The Serialization function
	int read_from_file(const string &path);

private:
	vw* _v;

};

MEDSERIALIZE_SUPPORT(MedVW)

#endif
#endif
